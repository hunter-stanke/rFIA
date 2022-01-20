growMortStarter <- function(x,
                            db,
                            grpBy_quo = NULL,
                            polys = NULL,
                            returnSpatial = FALSE,
                            bySpecies = FALSE,
                            bySizeClass = FALSE,
                            landType = 'forest',
                            treeType = 'all',
                            method = 'TI',
                            lambda = .5,
                            stateVar = 'TPA',
                            treeDomain = NULL,
                            areaDomain = NULL,
                            totals = FALSE,
                            byPlot = FALSE,
                            treeList = FALSE,
                            nCores = 1,
                            remote,
                            mr){



  ## Read required data, prep the database -------------------------------------
  if ('SUBP_COND_CHNG_MTRX' %in% names(db)) {
    reqTables <- c('PLOT', 'COND', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT',
                   'SUBP_COND_CHNG_MTRX',
                   'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                   'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  } else {
    reqTables <- c('PLOT', 'COND', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT',
                   'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                   'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  }


  ## If remote, read in state by state. Otherwise, drop all unnecessary tables
  db <- readRemoteHelper(x, db, remote, reqTables, nCores)

  ## Handle TX issues - we only keep inventory years that are present in BOTH
  ## EAST AND WEST TX
  db <- handleTX(db)




  ## Some warnings if inputs are bogus -----------------------------------------
  if (!is.null(polys) &
      dplyr::first(class(polys)) %in%
      c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (treeType %in% c('gs', 'all') == FALSE){
    stop('treeType must be one of: "all" or "gs".')
  }
  if (any(reqTables[!c(reqTables %in% 'SUBP_COND_CHNG_MTRX')] %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '),
               'not found in object db.'))
  }
  if (stringr::str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method,
                  'unknown. Defaulting to Temporally Indifferent (TI).'))
  }
  ## No EXP_GROW available for Western States, make sure we warn that values will be returned as 0
  # These states do not allow temporal queries. Things are extremely weird with their eval groups
  noGrow <- c(02,03,04,07,08,11,14,15,16, 30, 32, 35,43,49, 78)
  if(any(unique(db$PLOT$STATECD) %in% noGrow)){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    warning(paste('Recruitment, growth, and net change data unavailable for: ', toString(fancyName) , '. Returning 0 for all such estimates which include these states.', sep = ''))
  }
  # These states do not allow change estimates.
  if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    stop(paste('Growth & Mortality Estimates unavailable for: ', paste(as.character(fancyName), collapse = ', '), sep = ''))
  }




  ## Prep other variables ------------------------------------------------------
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>%
    dplyr::mutate(PLT_CN = CN,
           pltID = stringr::str_c(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))
  db$TREE <- db$TREE %>%
    dplyr::mutate(TRE_CN = CN)

  ## Convert grpBy to character
  grpBy <- grpByToChar(db, grpBy_quo)


  # I like a unique ID for a plot through time
  if (byPlot | byPlot) {grpBy <- c('pltID', grpBy)}


  ## Intersect plots with polygons if polygons are given
  if (!is.null(polys)){

    ## Add shapefile names to grpBy
    grpBy = c(grpBy, names(polys)[names(polys) != 'geometry'])
    ## Do the intersection
    db <- arealSumPrep2(db, grpBy, polys, nCores, remote)

    ## If there's nothing there, skip the state
    if (is.null(db)) return('no plots in polys')
  }

  ## If we want to return spatial plots
  if (byPlot & returnSpatial){
    grpBy <- c(grpBy, 'LON', 'LAT')
  }

  ### HANDLE THE STATE VARIABLE, only applying to the midpoint table for consistency
  if (stringr::str_to_upper(stateVar) == 'TPA'){
    db$TREE_GRM_MIDPT$state <- 1
    db$TREE$state_recr <- 1
  } else if (stringr::str_to_upper(stateVar) == 'BAA'){
    db$TREE_GRM_MIDPT$state <- basalArea(db$TREE_GRM_MIDPT$DIA)
    db$TREE$state_recr <- basalArea(db$TREE$DIA)
  } else if (stringr::str_to_upper(stateVar) == 'SAWVOL'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLCSNET
    db$TREE$state_recr <- db$TREE$VOLCSNET
  } else if (stringr::str_to_upper(stateVar) == 'SAWVOL_BF'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLBFNET
    db$TREE$state_recr <- db$TREE$VOLBFNET
  } else if (stringr::str_to_upper(stateVar) == 'NETVOL'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLCFNET
    db$TREE$state_recr <- db$TREE$VOLCFNET
  } else if (stringr::str_to_upper(stateVar) == 'SNDVOL'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLCFSND
    db$TREE$state_recr <- db$TREE$VOLCFSND
  } else if (stringr::str_to_upper(stateVar) == 'BIO_AG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_AG
    db$TREE$state_recr <- db$TREE$DRYBIO_AG
  } else if (stringr::str_to_upper(stateVar) == 'BIO_BG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_BG
    db$TREE$state_recr <- db$TREE$DRYBIO_BG
  } else if (stringr::str_to_upper(stateVar) == 'BIO'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_BG + db$TREE_GRM_MIDPT$DRYBIO_AG
    db$TREE$state_recr <- db$TREE$DRYBIO_BG + db$TREE$DRYBIO_AG
  } else if (stringr::str_to_upper(stateVar) == 'CARB_AG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_AG * .5
    db$TREE$state_recr <- db$TREE$DRYBIO_AG * .5
  } else if (stringr::str_to_upper(stateVar) == 'CARB_BG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_BG * .5
    db$TREE$state_recr <- db$TREE$DRYBIO_BG * .5
  } else if (stringr::str_to_upper(stateVar) == 'CARB'){
    db$TREE_GRM_MIDPT$state <- (db$TREE_GRM_MIDPT$DRYBIO_AG + db$TREE_GRM_MIDPT$DRYBIO_BG) * .5
    db$TREE$state_recr <- (db$TREE$DRYBIO_AG + db$TREE$DRYBIO_BG) * .5
  } else {
    stop(paste0('Method not known for stateVar: ', stateVar, '. Please choose one of: TPA, BAA, SAWVOL, SAWVOL_BF, NETVOL, BIO_AG, BIO_BG, BIO, CARB_AG, CARB_BG, or CARB.' ))
  }






  ## Build a domain indicator for each observation (1 or 0) --------------------
  ## Land type and tree type combined
  db <- typeDomain_grow(db, treeType, landType, type = 'gm', stateVar)

  ## Spatial boundary
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(!is.na(db$PLOT$polyID), 1, 0)
  } else {
    db$PLOT$sp <- 1
  }

  # User defined domain indicator for area (ex. specific forest type)
  db <- udAreaDomain(db, areaDomain)

  # User defined domain indicator for tree (ex. trees > 20 ft tall)
  db <- udTreeDomain(db, treeDomain)




  ## Handle population tables --------------------------------------------------
  ## Filtering out all inventories that are not relevant to the current estimation
  ## type. If using estimator other than TI, handle the differences in P2POINTCNT
  ## and in assigning YEAR column (YEAR = END_INVYR if method = 'TI')
  pops <- handlePops(db, evalType = c('GROW', 'MORT', 'REMV'), method, mr, ga = TRUE)

  ## Generate a list that tells us if a plot is ever associated with a growth
  ## accounting inventory. If so, we will be more strict with domain definitions
  plt.ga <- pops %>%
    dtplyr::lazy_dt() %>%
    dplyr::mutate(ga = case_when(GROWTH_ACCT == 'Y' ~ 1,
                                 TRUE ~ 0)) %>%
    dplyr::group_by(PLT_CN) %>%
    dplyr::summarise(ga = sum(ga, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ga = case_when(ga > 0 ~ 1, TRUE ~ 0)) %>%
    as.data.frame()

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (stringr::str_to_upper(method) != 'TI') {
    pops <- mergeSmallStrata(db, pops)
  }




  ## Canned groups -------------------------------------------------------------
  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      dplyr::left_join(dplyr::select(intData$REF_SPECIES_2018,
                       c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      dplyr::mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' ')) %>%
      dplyr::mutate_if(is.factor,
                as.character)
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }

  ## Break into size classes
  if (bySizeClass){
    grpBy <- c(grpBy, 'sizeClass')
    db$TREE$sizeClass <- makeClasses(db$TREE$DIA, interval = 2, numLabs = TRUE)
    db$TREE <- db$TREE[!is.na(db$TREE$sizeClass),]
  }




  ## Prep the tree list --------------------------------------------------------
  ## Narrow up the tables to the necessary variables
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy &
                           !c(names(db$COND) %in% grpP)]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy &
                           !c(names(db$TREE) %in% c(grpP, grpC))]

  ## Dropping irrelevant rows and columns
  db$PLOT <- db$PLOT %>%
    dplyr::select(c(PLT_CN, STATECD, MACRO_BREAKPOINT_DIA,
                    INVYR, MEASYEAR, PLOT_STATUS_CD,
                    dplyr::all_of(grpP), sp, COUNTYCD,
                    PREV_PLT_CN, REMPER)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    dplyr::filter(PLOT_STATUS_CD == 1 & sp == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)

  db$COND <- db$COND %>%
    dplyr::select(c(PLT_CN, CONDPROP_UNADJ, PROP_BASIS,
                    COND_STATUS_CD, CONDID,
                    dplyr::all_of(grpC), aD, landD)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    #dplyr::filter(aD == 1 & landD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))

  db$TREE_GRM_COMPONENT <- db$TREE_GRM_COMPONENT %>%
    dplyr::select(c(PLT_CN, TRE_CN, SUBPTYP_GRM, TPAGROW_UNADJ, TPARECR_UNADJ,
                    TPAREMV_UNADJ, TPAMORT_UNADJ, COMPONENT)) %>%
    dplyr::filter(TPAGROW_UNADJ > 0 | TPARECR_UNADJ > 0 | TPAREMV_UNADJ > 0 | TPAMORT_UNADJ > 0) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% db$PLOT$PLT_CN) %>%
    dplyr::select(-c(PLT_CN))

  db$TREE <- db$TREE %>%
    dplyr::select(c(PLT_CN, CONDID, PREVCOND, TRE_CN,
                    PREV_TRE_CN, SUBP, TREE, dplyr::all_of(grpT), tD,
                    typeD, state_recr, TPA_UNADJ,
                    STATUSCD, DIA)) %>%
    ## Drop plots outside our domain of interest
    dplyr::filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
    #dplyr::filter(TRE_CN %in% c(db$TREE_GRM_COMPONENT$TRE_CN, db$TREE_GRM_COMPONENT$PREV_TRE_CN))

  db$TREE_GRM_MIDPT <- db$TREE_GRM_MIDPT %>%
    select(c(TRE_CN, DIA, state)) %>%
    filter(TRE_CN %in% db$TREE_GRM_COMPONENT$TRE_CN)

  if ('SUBP_COND_CHNG_MTRX' %in% names(db)) {
    db$SUBP_COND_CHNG_MTRX <- dplyr::select(db$SUBP_COND_CHNG_MTRX, PLT_CN, PREV_PLT_CN,
                                     SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID) %>%
      dplyr::filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  }

  # Separate area grouping names from tree grouping names
  if (!is.null(polys)){
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND), names(polys))]
  } else {
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND))]
  }


  ## Full tree list
  data <- db$PLOT %>%
    dtplyr::lazy_dt() %>%
    dplyr::select(c(PLT_CN, STATECD, MACRO_BREAKPOINT_DIA, INVYR,
                                   MEASYEAR, PLOT_STATUS_CD, PREV_PLT_CN, REMPER,
                                   dplyr::all_of(grpP), sp)) %>%
    dplyr::left_join(dplyr::select(db$COND, c(PLT_CN, CONDPROP_UNADJ, PROP_BASIS,
                                              COND_STATUS_CD, CONDID, dplyr::all_of(grpC),
                                              aD, landD)),
                     by = c('PLT_CN')) %>%
    dplyr::left_join(dplyr::select(db$TREE, c(PLT_CN, CONDID, PREVCOND, TRE_CN,
                                              PREV_TRE_CN, SUBP, TREE, dplyr::all_of(grpT),
                                              tD, typeD, state_recr, TPA_UNADJ, STATUSCD, DIA)),
                     by = c('PLT_CN', 'CONDID')) %>%
    dplyr::left_join(dplyr::select(db$TREE_GRM_COMPONENT, c(TRE_CN, SUBPTYP_GRM,
                                                            TPAGROW_UNADJ, TPARECR_UNADJ,
                                                            TPAREMV_UNADJ, TPAMORT_UNADJ,
                                                            COMPONENT)),
                     by = c('TRE_CN')) %>%
    dplyr::left_join(dplyr::select(db$TREE_GRM_MIDPT, c(TRE_CN, DIA, state)),
                     by = c('TRE_CN'), suffix = c('', '.mid')) %>%
    dplyr::left_join(dplyr::select(db$PLOT, c(PLT_CN, dplyr::all_of(grpP), sp)),
                     by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    dplyr::left_join(dplyr::select(db$COND, c(PLT_CN, CONDID, landD, aD, dplyr::all_of(grpC),
                                              COND_STATUS_CD)),
                     by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'),
                     suffix = c('', '.prev')) %>%
    dplyr::left_join(dplyr::select(db$TREE, c(TRE_CN, dplyr::all_of(grpT),
                                              typeD, tD, DIA, STATUSCD, TPA_UNADJ, state.prev = state_recr)),
                     by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
    # dplyr::mutate_if(is.factor,
    #           as.character) %>%
    dplyr::mutate(TPAREMV_UNADJ = TPAREMV_UNADJ * state,
                  TPAMORT_UNADJ = TPAMORT_UNADJ * state,
                  TPARECR_UNADJ = TPARECR_UNADJ * state_recr / REMPER,
                  ## State recruit is the state variable adjustment for ALL TREES at T2,
                  ## So we can estimate live TPA at t2 (t1 unavailable w/out growth accounting) with:
                  TPA_UNADJ = TPAGROW_UNADJ * state_recr * ifelse(COMPONENT %in% c('SURVIVOR', 'INGROWTH'), 1, 0),
                  #TPA_UNADJ = TPA_UNADJ * state_recr * ifelse(STATUSCD == 1 & DIA >= 5, 1, 0),
                  TPA_UNADJ.prev = TPAGROW_UNADJ * state.prev * ifelse(COMPONENT %in% c('SURVIVOR'), 1, 0),
                  ## Add our indicator of whether or not a plot is ever associated with a,
                  #TPA_UNADJ.prev = TPA_UNADJ.prev * state.prev * ifelse(STATUSCD.prev == 1 & STATUSCD != 0 & DIA.prev >= 5, 1, 0)
                  ) %>%
    ## Add our indicator of whether or not a plot is ever associated with a
    ## growth accounting inventory
    dplyr::left_join(plt.ga, by = 'PLT_CN') %>%
    ## Set up our domain indicators accordingly
    dplyr::mutate(aChng = dplyr::case_when(ga == 1 & COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1 & !is.null(CONDPROP_UNADJ) ~ 1,
                                           ga == 0 ~ 1,
                                           TRUE ~ 0),
                  tChng = dplyr::case_when(ga == 1 & COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1 ~ 1,
                                           ga == 0 ~ 1,
                                           TRUE ~ 0),
                  landD.prev = dplyr::case_when(ga == 1 & landD == 1 & landD.prev == 1 ~ 1,
                                                ga == 0 & is.na(landD) & is.na(landD.prev) ~ 0,
                                                ga == 0 & is.na(landD.prev) ~ landD,
                                                ga == 0 ~ landD.prev,
                                                TRUE ~ 0)) %>%
    # If previous attributes are unavailable for trees, default to current
    dplyr::mutate(tD.prev = dplyr::case_when(is.na(tD.prev) ~ tD, TRUE ~ tD.prev),
                  typeD.prev = dplyr::case_when(is.na(typeD.prev) ~ typeD, TRUE ~ typeD.prev),
                  aD.prev = dplyr::case_when(is.na(aD.prev) ~ aD, TRUE ~ aD.prev),
                  sp.prev = dplyr::case_when(is.na(sp.prev) ~ sp, TRUE ~ sp.prev)) %>%
    # Comprehensive domain indicators
    dplyr::mutate(tDI = landD.prev * aD.prev * tD.prev * typeD.prev * sp.prev * tChng,
                  tDI_r = landD * aD * tD * typeD * sp * tChng, # All previous attributes NA for recruitment
                  aDI = landD * aD * sp * aChng) %>%
    as.data.frame()


  if ('SUBP_COND_CHNG_MTRX' %in% names(db)) {
    ### DOING AREA SEPARATELY NOW FOR GROWTH ACCOUNTING PLOTS
    aData <- dplyr::select(db$PLOT, c(PLT_CN, STATECD, MACRO_BREAKPOINT_DIA,
                                      INVYR, MEASYEAR, PLOT_STATUS_CD, PREV_PLT_CN,
                                      REMPER, dplyr::all_of(grpP), sp)) %>%
      dplyr::left_join(dplyr::select(db$SUBP_COND_CHNG_MTRX, PLT_CN, PREV_PLT_CN,
                                     SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID),
                       by = c('PLT_CN', 'PREV_PLT_CN')) %>%
      dplyr::left_join(dplyr::select(db$COND, c(PLT_CN, CONDPROP_UNADJ, PROP_BASIS,
                                                COND_STATUS_CD, CONDID, dplyr::all_of(grpC),
                                                aD, landD)),
                       by = c('PLT_CN', 'CONDID')) %>%
      dplyr::left_join(dplyr::select(db$COND, c(PLT_CN, PROP_BASIS, COND_STATUS_CD,
                                                CONDID, dplyr::all_of(grpC), aD, landD)),
                       by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
      dplyr::mutate(aChng = dplyr::if_else(COND_STATUS_CD == 1 &
                                             COND_STATUS_CD.prev == 1 &
                                             !is.null(CONDPROP_UNADJ) &
                                             SUBPTYP == 1,
                                           1, 0),
                    SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25)

    aData$landD <- ifelse(aData$landD == 1 & aData$landD.prev == 1, 1, 0)
    aData$aDI_ga <- aData$landD * aData$aD * aData$sp * aData$aChng

  }



  ## Plot-level summaries ------------------------------------------------------
  if (byPlot & !treeList){

    grpBy <- c('YEAR', grpBy)
    grpSyms <- syms(grpBy)
    aGrpSyms <- syms(aGrpBy)


    if ('SUBP_COND_CHNG_MTRX' %in% names(db)) {
      ## Forested land area w growth accounting
      a.ga <- aData %>%
        ## Will be lots of trees here, so CONDPROP listed multiple times
        ## Adding PROP_BASIS so we can handle adjustment factors at strata level
        #distinct(PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
        dtplyr::lazy_dt() %>%
        dplyr::group_by(PLT_CN, !!!aGrpSyms) %>%
        dplyr::summarize(fa_ga = sum(SUBPTYP_PROP_CHNG * aDI_ga, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        as.data.frame()

      a <- data %>%
        ## Will be lots of trees here, so CONDPROP listed multiple times
        ## Adding PROP_BASIS so we can handle adjustment factors at strata level
        dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
        dtplyr::lazy_dt() %>%
        dplyr::group_by(PLT_CN, !!!aGrpSyms) %>%
        dplyr::summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(dplyr::select(a.ga, PLT_CN, !!!aGrpSyms, fa_ga),
                         by = c('PLT_CN', aGrpBy)) %>%
        dplyr::left_join(plt.ga, by = 'PLT_CN') %>%
        dplyr::mutate(PROP_FOREST = case_when(ga == 1 ~ fa_ga,
                                              TRUE ~ fa)) %>%
        dplyr::select(PLT_CN, !!!aGrpSyms, PROP_FOREST) %>%
        as.data.frame()

    } else {
      ## Forested land area w/out growth accounting
      a <- data %>%
        ## Will be lots of trees here, so CONDPROP listed multiple times
        ## Adding PROP_BASIS so we can handle adjustment factors at strata level
        dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
        dtplyr::lazy_dt() %>%
        dplyr::group_by(PLT_CN, !!!aGrpSyms) %>%
        dplyr::summarize(PROP_FOREST = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }

    t <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      dplyr::distinct(PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      # Compute estimates at plot level
      dplyr::group_by(!!!grpSyms, PLT_CN, REMPER) %>%
      dplyr::summarise(RECR_TPA = sum(TPARECR_UNADJ * tDI_r, na.rm = TRUE),
                       MORT_TPA = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
                       REMV_TPA = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
                       CURR_TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
                       PREV_TPA = sum(TPA_UNADJ.prev * tDI, na.rm = TRUE)) %>%
      dplyr::mutate(PREV_TPA = PREV_TPA + (MORT_TPA + REMV_TPA)*REMPER) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(CHNG_TPA = (CURR_TPA - PREV_TPA) / REMPER,
                    GROW_TPA = CHNG_TPA - RECR_TPA + MORT_TPA + REMV_TPA,
                    RECR_PERC = RECR_TPA / PREV_TPA * 100,
                    MORT_PERC = MORT_TPA / PREV_TPA * 100,
                    REMV_PERC = REMV_TPA / PREV_TPA * 100,
                    GROW_PERC = GROW_TPA / PREV_TPA * 100,
                    CHNG_PERC = CHNG_TPA / PREV_TPA * 100) %>%
      dplyr::select(PLT_CN, !!!grpSyms, REMPER, RECR_TPA:REMV_TPA, GROW_TPA, CHNG_TPA, RECR_PERC:CHNG_PERC, PREV_TPA, CURR_TPA) %>%
      as.data.frame() %>%
      ## Rounding errors will generate tiny values, make them zero
      dplyr::mutate(GROW_TPA = dplyr::case_when(abs(GROW_TPA) < 1e-5 ~ 0,
                                                TRUE ~ GROW_TPA)) %>%
      dplyr::left_join(a, by = c('PLT_CN', aGrpBy)) %>%
      dplyr::distinct()


    ## Make it spatial
    if (returnSpatial){
      t <- t %>%
        dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
        sf::st_as_sf(coords = c('LON', 'LAT'),
                     crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }

    out <- list(tEst = t, grpBy = grpBy, aGrpBy = aGrpBy)




    ## Population estimation -----------------------------------------------------
  } else {

    aGrpSyms <- syms(aGrpBy)
    ## Forested area
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dplyr::mutate(fa = CONDPROP_UNADJ * aDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!aGrpSyms, fa)


    if ('SUBP_COND_CHNG_MTRX' %in% names(db)) {
      ### Plot-level estimates -- growth accounting
      a_ga <- aData %>%
        ## Will be lots of trees here, so CONDPROP listed multiple times
        ## Adding PROP_BASIS so we can handle adjustment factors at strata level
        #distinct(PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
        dtplyr::lazy_dt() %>%
        dplyr::filter(!is.na(PROP_BASIS)) %>%
        dplyr::group_by(PLT_CN, PROP_BASIS, CONDID, !!!aGrpSyms) %>%
        dplyr::summarize(fa_ga = sum(SUBPTYP_PROP_CHNG * aDI_ga, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        as.data.frame()

      a <- a %>%
        dplyr::left_join(dplyr::select(a_ga, PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!aGrpSyms, fa_ga),
                  by = c('PLT_CN', 'AREA_BASIS', 'CONDID', aGrpBy)) %>%
        dplyr::left_join(plt.ga, by = 'PLT_CN') %>%
        dplyr::mutate(fa = case_when(ga == 1 ~ fa_ga,
                                              TRUE ~ fa)) %>%
        dplyr::select(PLT_CN, AREA_BASIS, CONDID, !!!aGrpSyms, fa) %>%
        dplyr::filter(fa > 0)
    }

    ## Tree list
    grpSyms <- syms(grpBy)
    ## Tree list
    t <- data %>%
      dplyr::distinct(PLT_CN, TRE_CN, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::filter(!is.na(SUBPTYP_GRM)) %>%
      dplyr::filter(tDI > 0 | tDI_r > 0) %>%
      # Compute estimates at plot level
      dplyr::mutate(rPlot = TPARECR_UNADJ * tDI_r,
                    mPlot = TPAMORT_UNADJ * tDI,
                    hPlot = TPAREMV_UNADJ * tDI,
                    tPlot = TPA_UNADJ * tDI,
                    pPlot = (TPA_UNADJ.prev * tDI) + ((mPlot + hPlot)*REMPER),
                    cPlot = (tPlot - pPlot) / REMPER,
                    gPlot = cPlot - rPlot + mPlot + hPlot) %>%
      dplyr::mutate(TREE_BASIS = case_when(SUBPTYP_GRM == 0 ~ NA_character_,
                                           SUBPTYP_GRM == 1 ~ 'SUBP',
                                           SUBPTYP_GRM == 2 ~ 'MICR',
                                           SUBPTYP_GRM == 3 ~ 'MACR')) %>%
      as.data.frame() %>%
      dplyr::select(PLT_CN, TREE_BASIS, SUBP, TREE, !!!grpSyms, rPlot:gPlot)


    if (treeList) {

      tEst <- a %>%
        dplyr::left_join(t, by = c('PLT_CN', aGrpBy)) %>%
        dplyr::mutate(EVAL_TYP = list(c('GROW', 'MORT', 'REMV'))) %>%
        dplyr::select(PLT_CN, EVAL_TYP, TREE_BASIS, AREA_BASIS,
                      !!!grpSyms, CONDID, SUBP, TREE,
                      RECR_TPA = rPlot,
                      MORT_TPA = mPlot,
                      REMV_TPA = hPlot,
                      GROW_TPA = gPlot,
                      CHNG_TPA = cPlot,
                      CURR_TPA = tPlot,
                      PREV_TPA = pPlot,
                      PROP_FOREST = fa)
      out <- list(tEst = tEst, aEst = NULL, grpBy = grpBy, aGrpBy = aGrpBy)

    } else {

      ## Sum variable(s) up to plot-level and adjust for non-response
      tPlt <- sumToPlot(t, pops, grpBy)
      aPlt <- sumToPlot(a, pops, aGrpBy)

      ## Adding YEAR to groups
      grpBy <- c('YEAR', grpBy)
      aGrpBy <- c('YEAR', aGrpBy)


      ## Sum variable(s) up to strata then estimation unit level
      eu.sums <- sumToEU(db, tPlt, aPlt, pops, grpBy, aGrpBy, method)
      tEst <- eu.sums$x
      aEst <- eu.sums$y

      ## Have to repeat this with tree totals as the denominator
      eu.sums <- sumToEU(db, dplyr::select(tPlt, -c(tPlot, pPlot)), dplyr::select(tPlt, -c(rPlot, mPlot, hPlot, gPlot, cPlot, tPlot)), pops, grpBy, grpBy, method)
      ttEst <- eu.sums$x %>%
        dplyr::select(ESTN_UNIT_CN, all_of(grpBy),
                      rPlot_cv_t = rPlot_cv, mPlot_cv_t = mPlot_cv, hPlot_cv_t = hPlot_cv,
                      gPlot_cv_t = gPlot_cv, cPlot_cv_t = cPlot_cv)
      tEst <- dplyr::left_join(tEst, ttEst, by = c('ESTN_UNIT_CN', grpBy))

      out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy)
    }


  }


  return(out)

}


#' @export
growMort <- function(db,
                     grpBy = NULL,
                     polys = NULL,
                     returnSpatial = FALSE,
                     bySpecies = FALSE,
                     bySizeClass = FALSE,
                     landType = 'forest',
                     treeType = 'all',
                     method = 'TI',
                     lambda = .5,
                     stateVar = 'TPA',
                     treeDomain = NULL,
                     areaDomain = NULL,
                     totals = FALSE,
                     variance = FALSE,
                     byPlot = FALSE,
                     treeList = FALSE,
                     nCores = 1) {

  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)
  areaDomain <- rlang::enquo(areaDomain)
  treeDomain <- rlang::enquo(treeDomain)


  ## Handle iterator if db is remote
  remote <- ifelse(class(db) == 'Remote.FIA.Database', 1, 0)
  iter <- remoteIter(db, remote)

  ## Check for a most recent subset
  mr <- checkMR(db, remote)

  ## prep for areal summary
  polys <- arealSumPrep1(polys)



  ## Run the main portion
  out <- lapply(X = iter, FUN = growMortStarter, db,
                grpBy_quo, polys, returnSpatial,
                bySpecies, bySizeClass,
                landType, treeType, method,
                lambda, stateVar, treeDomain, areaDomain,
                totals, byPlot, treeList,
                nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
  aEst <- dplyr::bind_rows(out[names(out) == 'aEst'])
  tEst <- dplyr::bind_rows(out[names(out) == 'tEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  aGrpBy <- out[names(out) == 'aGrpBy'][[1]]
  grpSyms <- dplyr::syms(grpBy)
  aGrpSyms <- dplyr::syms(aGrpBy)


  ## Summarize population estimates across estimation units
  if (!byPlot & !treeList){

    ## Combine most-recent population estimates across states with potentially
    ## different reporting schedules, i.e., if 2016 is most recent in MI and 2017 is
    ## most recent in WI, combine them and label as 2017
    if (mr) {
      tEst <- combineMR(tEst, grpBy)
      aEst <- combineMR(aEst, aGrpBy)
    }



    ## Totals and ratios -------------------------------------------------------
    aEst <- aEst %>%
      dplyr::group_by( !!!aGrpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::select(!!!aGrpSyms, fa_mean, fa_var, nPlots.y)


    tEst <- tEst %>%
      dplyr::group_by(!!!grpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::left_join(aEst, by = aGrpBy) %>%
      dplyr::mutate(CURR_TOTAL = tPlot_mean,
                    PREV_TOTAL = pPlot_mean,
                    RECR_TOTAL = rPlot_mean,
                    MORT_TOTAL = mPlot_mean,
                    REMV_TOTAL = hPlot_mean,
                    GROW_TOTAL = gPlot_mean,
                    CHNG_TOTAL = cPlot_mean,
                    AREA_TOTAL = fa_mean,
                    # Ratios
                    RECR_TPA = rPlot_mean / fa_mean,
                    MORT_TPA = mPlot_mean / fa_mean,
                    REMV_TPA = hPlot_mean / fa_mean,
                    GROW_TPA = gPlot_mean / fa_mean,
                    CHNG_TPA = cPlot_mean / fa_mean,
                    RECR_PERC = rPlot_mean / pPlot_mean,
                    MORT_PERC = mPlot_mean / pPlot_mean,
                    REMV_PERC = hPlot_mean / pPlot_mean,
                    GROW_PERC = gPlot_mean / pPlot_mean,
                    CHNG_PERC = cPlot_mean / pPlot_mean,
                    # Variances
                    CURR_TOTAL_VAR = tPlot_var,
                    PREV_TOTAL_VAR = pPlot_var,
                    RECR_TOTAL_VAR = rPlot_var,
                    MORT_TOTAL_VAR = mPlot_var,
                    REMV_TOTAL_VAR = hPlot_var,
                    GROW_TOTAL_VAR = gPlot_var,
                    CHNG_TOTAL_VAR = cPlot_var,
                    AREA_TOTAL_VAR = fa_var,
                    RECR_TPA_VAR = ratioVar(rPlot_mean, fa_mean, rPlot_var, fa_var, rPlot_cv),
                    MORT_TPA_VAR = ratioVar(mPlot_mean, fa_mean, mPlot_var, fa_var, mPlot_cv),
                    REMV_TPA_VAR = ratioVar(hPlot_mean, fa_mean, hPlot_var, fa_var, hPlot_cv),
                    GROW_TPA_VAR = ratioVar(gPlot_mean, fa_mean, gPlot_var, fa_var, gPlot_cv),
                    CHNG_TPA_VAR = ratioVar(cPlot_mean, fa_mean, cPlot_var, fa_var, cPlot_cv),
                    RECR_PERC_VAR = ratioVar(rPlot_mean, pPlot_mean, rPlot_var, pPlot_var, rPlot_cv_t),
                    MORT_PERC_VAR = ratioVar(mPlot_mean, pPlot_mean, mPlot_var, pPlot_var, mPlot_cv_t),
                    REMV_PERC_VAR = ratioVar(hPlot_mean, pPlot_mean, hPlot_var, pPlot_var, hPlot_cv_t),
                    GROW_PERC_VAR = ratioVar(gPlot_mean, pPlot_mean, gPlot_var, pPlot_var, gPlot_cv_t),
                    CHNG_PERC_VAR = ratioVar(cPlot_mean, pPlot_mean, cPlot_var, pPlot_var, cPlot_cv_t),

                    # Convert to percentages
                    RECR_PERC = RECR_PERC * 100,
                    MORT_PERC = MORT_PERC * 100,
                    REMV_PERC = REMV_PERC * 100,
                    GROW_PERC = GROW_PERC * 100,
                    CHNG_PERC = CHNG_PERC * 100,
                    RECR_PERC_VAR = RECR_PERC_VAR * 100^2,
                    MORT_PERC_VAR = MORT_PERC_VAR * 100^2,
                    REMV_PERC_VAR = REMV_PERC_VAR * 100^2,
                    GROW_PERC_VAR = GROW_PERC_VAR * 100^2,
                    CHNG_PERC_VAR = CHNG_PERC_VAR * 100^2,

                    # Sampling Errors
                    CURR_TOTAL_SE = sqrt(tPlot_var) / tPlot_mean * 100,
                    PREV_TOTAL_SE = sqrt(pPlot_var) / pPlot_mean * 100,
                    RECR_TOTAL_SE = sqrt(rPlot_var) / rPlot_mean * 100,
                    MORT_TOTAL_SE = sqrt(mPlot_var) / mPlot_mean * 100,
                    REMV_TOTAL_SE = sqrt(hPlot_var) / rPlot_mean * 100,
                    GROW_TOTAL_SE = sqrt(gPlot_var) / gPlot_mean * 100,
                    CHNG_TOTAL_SE = sqrt(cPlot_var) / cPlot_mean * 100,
                    AREA_TOTAL_SE = sqrt(fa_var) / fa_mean * 100,
                    RECR_TPA_SE = sqrt(RECR_TPA_VAR) / RECR_TPA * 100,
                    MORT_TPA_SE = sqrt(MORT_TPA_VAR) / MORT_TPA * 100,
                    REMV_TPA_SE = sqrt(REMV_TPA_VAR) / REMV_TPA * 100,
                    RECR_PERC_SE = sqrt(RECR_PERC_VAR) / RECR_PERC * 100,
                    MORT_PERC_SE = sqrt(MORT_PERC_VAR) / MORT_PERC * 100,
                    REMV_PERC_SE = sqrt(REMV_PERC_VAR) / REMV_PERC * 100,
                    GROW_PERC_SE = sqrt(GROW_PERC_VAR) / GROW_PERC * 100,
                    CHNG_PERC_SE = sqrt(CHNG_PERC_VAR) / CHNG_PERC * 100,

                    # Plot counts
                    nPlots_TREE = nPlots.x,
                    nPlots_AREA = nPlots.y,
                    N = P2PNTCNT_EU) %>%
      dplyr::select(!!!grpSyms, RECR_TPA:CHNG_PERC,
                    RECR_TOTAL:CHNG_TOTAL, PREV_TOTAL, CURR_TOTAL, AREA_TOTAL,
                    RECR_TPA_VAR:CHNG_PERC_VAR,
                    RECR_TOTAL_VAR:CHNG_TOTAL_VAR, PREV_TOTAL_VAR, CURR_TOTAL_VAR, AREA_TOTAL_VAR,
                    RECR_TPA_SE:CHNG_PERC_SE,
                    RECR_TOTAL_SE:CHNG_TOTAL_SE, PREV_TOTAL_SE, CURR_TOTAL_SE, AREA_TOTAL_SE,
                    nPlots_TREE, nPlots_AREA, N) %>%
      ## Rounding errors can cause GROW_TPA to take an extremely small value instead of zero
      ## Make it zero when this happens
      dplyr::mutate(dplyr::across(c(GROW_TPA, GROW_PERC, GROW_TOTAL,
                                    GROW_TPA_VAR, GROW_PERC_VAR, GROW_TOTAL_VAR),
                                  .fns = ~case_when(abs(.x) < 1e-5 ~ 0,
                                                   TRUE ~ .x)))

    ## Drop totals unless told not to
    if (!totals) {
      tEst <- tEst[,!stringr::str_detect(names(tEst), '_TOTAL')]
    }

    ## Select either variance or SE, depending on input
    if (variance) {
      tEst <- tEst[,!stringr::str_detect(names(tEst), '_SE')]
    } else {
      tEst <- tEst[,!stringr::str_detect(names(tEst), '_VAR')]
    }

  }

  ## Modify some names if a different state variable was given
  if (stateVar != 'TPA') {
    names(tEst) <- str_replace(names(tEst), 'TPA', paste(stateVar, 'ACRE', sep = '_'))
    #names(tEst) <- str_replace(names(tEst), 'TREE', ifelse(stateVar == 'BAA', 'BA', stateVar))
  }
  names(tEst) <- str_replace(names(tEst), 'BAA_ACRE', 'BAA')

  ## Pretty output
  tEst <- tEst %>%
    dplyr::ungroup() %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    as_tibble()

  # We don't include YEAR in treeList output, and NA groups will be important
  # for retaining non-treed forestland
  if (!treeList) {
    tEst <- tEst %>%
      tidyr::drop_na(grpBy[!c(grpBy %in% names(polys))]) %>%
      dplyr::arrange(YEAR)
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## For spatial polygons
  if (returnSpatial & !byPlot) {
    tEst <- dplyr::left_join(tEst,
                             as.data.frame(dplyr::select(polys, polyID, geometry)),
                             by = 'polyID')
  }

  ## Above converts to tibble
  if (returnSpatial) tEst <- sf::st_sf(tEst)

  return(tEst)

}
