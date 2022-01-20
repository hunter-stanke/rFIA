vrStarter <- function(x,
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
                      treeDomain = NULL,
                      areaDomain = NULL,
                      totals = FALSE,
                      byPlot = FALSE,
                      treeList = FALSE,
                      nCores = 1,
                      remote,
                      mr){





  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'TREE', 'COND',
                 'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT', 'TREE_GRM_BEGIN',
                 'SUBP_COND_CHNG_MTRX',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')


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
  if (treeType %in% c('live', 'gs', 'all') == FALSE){
    stop('treeType must be one of: "live", "gs", or "all".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '),
               'not found in object db.'))
  }
  if (stringr::str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method,
                  'unknown. Defaulting to Temporally Indifferent (TI).'))
  }
  ## No EXP_GROW available for most Western States, make sure we warn that values will be returned as 0
  # These states do not allow temporal queries. Things are extremely weird with their eval groups
  noGrow <- c(02,03,04,07,08,11,14,15,16, 30, 32, 35,43,49, 78)
  if(any(unique(db$PLOT$STATECD) %in% noGrow)){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    warning(paste('Growth data unavailable for: ', toString(fancyName) , '. Returning 0 for all growth estimates which include these states.', sep = ''))
  }

  # These states do not allow change estimates.
  if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    stop(paste('Growth & Mortality Estimates unavailable for: ', as.character(fancyName), sep = ''))
  }


  ## Prep other variables ------------------------------------------------------
  ## Need a plotCN, and a new ID
  db$TREE <- db[['TREE']] %>% dplyr::mutate(TRE_CN = CN)
  db$PLOT <- db$PLOT %>% dplyr::mutate(PLT_CN = CN,
                                pltID = stringr::str_c(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ## Convert grpBy to character
  grpBy <- grpByToChar(db, grpBy_quo)


  # I like a unique ID for a plot through time
  if (byPlot | treeList) {grpBy <- c('pltID', grpBy)}


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






  ## Build a domain indicator for each observation (1 or 0) --------------------

  ## Land type and tree type combined
  db <- typeDomain_grow(db, treeType, landType, type = 'vr')

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
  pops <- handlePops(db, evalType = c('GROW'), method, mr)

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
  db$TREE <- db$TREE %>%
    dplyr::select(c(PLT_CN, CONDID, PREVCOND, TRE_CN,
                               PREV_TRE_CN, SUBP, TREE, dplyr::all_of(grpT), tD,
                               typeD, DIA, DRYBIO_AG, VOLCFNET,
                               VOLBFNET, STATUSCD)) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  db$TREE_GRM_COMPONENT <- db$TREE_GRM_COMPONENT %>%
    dplyr::select(c(TRE_CN, SUBPTYP_GRM, TPAGROW_UNADJ,
             TPAREMV_UNADJ, TPAMORT_UNADJ, COMPONENT)) %>%
    dplyr::filter(TRE_CN %in% db$TREE$TRE_CN)
  db$TREE_GRM_MIDPT <- db$TREE_GRM_MIDPT %>%
    dplyr::select(c(TRE_CN, DIA, VOLCFNET, VOLBFNET, DRYBIO_AG)) %>%
    dplyr::filter(TRE_CN %in% db$TREE$TRE_CN)
  db$TREE_GRM_BEGIN <- db$TREE_GRM_BEGIN %>%
    dplyr::select(c(TRE_CN, DIA, VOLCFNET, VOLBFNET, DRYBIO_AG)) %>%
    dplyr::filter(TRE_CN %in% db$TREE$TRE_CN)
  db$SUBP_COND_CHNG_MTRX <- db$SUBP_COND_CHNG_MTRX %>%
    dplyr::select(PLT_CN, PREV_PLT_CN, SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID) %>%
    dplyr::filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))


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
                                              COND_STATUS_CD, CONDID,
                                              dplyr::all_of(grpC), aD, landD)),
                     by = c('PLT_CN')) %>%
    dplyr::left_join(dplyr::select(db$TREE, c(PLT_CN, CONDID, PREVCOND, TRE_CN,
                                              PREV_TRE_CN, SUBP, TREE,
                                              dplyr::all_of(grpT), tD, typeD,
                                              DIA, DRYBIO_AG, VOLCFNET,
                                              VOLBFNET, STATUSCD)),
                     by = c('PLT_CN', 'CONDID')) %>%
    dplyr::left_join(dplyr::select(db$TREE_GRM_COMPONENT, c(TRE_CN, SUBPTYP_GRM,
                                                            TPAGROW_UNADJ, TPAREMV_UNADJ,
                                                            TPAMORT_UNADJ, COMPONENT)),
                     by = c('TRE_CN')) %>%
    dplyr::left_join(dplyr::select(db$TREE_GRM_MIDPT, c(TRE_CN, DIA, VOLCFNET,
                                                        VOLBFNET, DRYBIO_AG)),
                     by = c('TRE_CN'), suffix = c('', '.mid')) %>%
    dplyr::left_join(dplyr::select(db$TREE_GRM_BEGIN, c(TRE_CN, DIA, VOLCFNET,
                                                        VOLBFNET, DRYBIO_AG)),
                     by = c('TRE_CN'), suffix = c('', '.beg')) %>%
    dplyr::left_join(dplyr::select(db$PLOT, c(PLT_CN, dplyr::all_of(grpP), sp)),
                     by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    dplyr::left_join(dplyr::select(db$COND, c(PLT_CN, CONDID, landD, aD,
                                              dplyr::all_of(grpC), COND_STATUS_CD)),
                     by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
    dplyr::left_join(dplyr::select(db$TREE, c(TRE_CN, dplyr::all_of(grpT), typeD,
                                              tD, DIA,  DRYBIO_AG, VOLCFNET,
                                              VOLBFNET, STATUSCD)),
                     by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
    ## Some domain indicators
    dplyr::mutate(aChng = dplyr::case_when(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1 & !is.null(CONDPROP_UNADJ) ~ 1,
                                           TRUE ~ 0),
                  tChng = dplyr::case_when(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1 ~ 1,
                                           TRUE ~ 0),
                  landD.prev = dplyr::case_when(landD == 1 & landD.prev == 1 ~ 1,
                                                TRUE ~ 0),
                  status = case_when(COMPONENT == 'SURVIVOR' ~ 1,
                                     TRUE ~ 0)) %>%
    # If previous attributes are unavailable for trees, default to current
    dplyr::mutate(tD.prev = dplyr::case_when(is.na(tD.prev) ~ tD, TRUE ~ tD.prev),
                  typeD.prev = dplyr::case_when(is.na(typeD.prev) ~ typeD, TRUE ~ typeD.prev),
                  aD.prev = dplyr::case_when(is.na(aD.prev) ~ aD, TRUE ~ aD.prev),
                  sp.prev = dplyr::case_when(is.na(sp.prev) ~ sp, TRUE ~ sp.prev)) %>%
    # Comprehensive domain indicators
    dplyr::mutate(tDI = landD.prev * aD.prev * tD.prev * typeD.prev * sp.prev * tChng,
                  aDI = landD.prev * aD * sp * aChng) %>%
    as.data.frame() %>%
    distinct()


  ## Only if we're only considering live trees that stayed live
  if (tolower(treeType) == 'live') {data$tDI <- data$tDI * data$status}


  ## Modify  attributes depending on component (mortality uses midpoint)
  data <- data %>%
    dplyr::mutate(DIA2 = vrAttHelper(DIA, DIA.prev, DIA.mid, DIA.beg, COMPONENT, REMPER, 2) * tDI,
                  DIA1 = vrAttHelper(DIA, DIA.prev, DIA.mid, DIA.beg, COMPONENT, REMPER, 1) * tDI,
                  BA2 = vrAttHelper(basalArea(DIA), basalArea(DIA.prev), basalArea(DIA.mid), basalArea(DIA.beg), COMPONENT, REMPER, 2) * tDI,
                  BA1 = vrAttHelper(basalArea(DIA), basalArea(DIA.prev), basalArea(DIA.mid), basalArea(DIA.beg), COMPONENT, REMPER, 1) * tDI,
                  VOLCFNET2 = vrAttHelper(VOLCFNET, VOLCFNET.prev, VOLCFNET.mid, VOLCFNET.beg, COMPONENT, REMPER, 2) * tDI,
                  VOLCFNET1 = vrAttHelper(VOLCFNET, VOLCFNET.prev, VOLCFNET.mid, VOLCFNET.beg, COMPONENT, REMPER, 1) * tDI,
                  VOLBFNET2 = vrAttHelper(VOLBFNET, VOLBFNET.prev, VOLBFNET.mid, VOLBFNET.beg, COMPONENT, REMPER, 2) * tDI,
                  VOLBFNET1 = vrAttHelper(VOLBFNET, VOLBFNET.prev, VOLBFNET.mid, VOLBFNET.beg, COMPONENT, REMPER, 1) * tDI,
                  DRYBIO_AG2 = vrAttHelper(DRYBIO_AG, DRYBIO_AG.prev, DRYBIO_AG.mid, DRYBIO_AG.beg, COMPONENT, REMPER, 2) * tDI,
                  DRYBIO_AG1 = vrAttHelper(DRYBIO_AG, DRYBIO_AG.prev, DRYBIO_AG.mid, DRYBIO_AG.beg, COMPONENT, REMPER, 1) * tDI) %>%
    dplyr::select(-c(DIA.mid, VOLCFNET.mid, VOLBFNET.mid, DRYBIO_AG.mid,
                     DIA.prev, VOLCFNET.prev, VOLBFNET.prev, DRYBIO_AG.prev,
                     DIA.beg, VOLCFNET.beg, VOLBFNET.beg, DRYBIO_AG.beg,
                     DIA, VOLCFNET, VOLBFNET, DRYBIO_AG))

  ## Just what we need
  data <- data %>%
    dplyr::select(PLT_CN, TRE_CN, SUBP, CONDID, TREE, tDI,
                  grpP, grpC, grpT, TPAGROW_UNADJ, PROP_BASIS, SUBPTYP_GRM, PLOT_STATUS_CD,
                  DIA2, DIA1, BA2, BA1, DRYBIO_AG2, DRYBIO_AG1, VOLCFNET2, VOLCFNET1, VOLBFNET2, VOLBFNET1, MEASYEAR) %>%
    ## Rearrange previous values as observations
    tidyr::pivot_longer(cols = DIA2:VOLBFNET1,
                        names_to = c(".value", 'ONEORTWO'),
                        names_sep = -1)

  ### DOING AREA SEPARATELY NOW FOR GROWTH ACCOUNTING PLOTS
  aData <- db$PLOT %>%
    dplyr::select(c(PLT_CN, STATECD, MACRO_BREAKPOINT_DIA, INVYR, MEASYEAR,
                    PLOT_STATUS_CD, PREV_PLT_CN, REMPER, dplyr::all_of(grpP), sp)) %>%
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
    dplyr::mutate(aChng = dplyr::case_when(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1 & !is.null(CONDPROP_UNADJ) & SUBPTYP == 1 ~ 1 ,
                                           TRUE ~ 0),
                  SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25) %>%
    ## Comprehensive domain indicator
    dplyr::mutate(aDI = landD * landD.prev * aD * sp * aChng)







  ## Plot-level summaries ------------------------------------------------------
  if (byPlot & !treeList){

    grpBy <- c('YEAR', grpBy)
    grpSyms <- syms(grpBy)
    aGrpSyms <- syms(aGrpBy)

    ### Plot-level estimates
    a <- aData %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(PLT_CN, !!!aGrpSyms) %>%
      dplyr::summarize(PROP_FOREST = sum(SUBPTYP_PROP_CHNG * aDI, na.rm = TRUE)) %>%
      as.data.frame()

    t <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!grpSyms, PLT_CN) %>%
      dplyr::summarize(t = sum(TPAGROW_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), ## Previous only
                       d = sum(DIA * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                       ba = sum(BA * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                       vol = sum(VOLCFNET * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                       svol = sum(VOLBFNET * TPAGROW_UNADJ * tDI, na.rm = TRUE) / 1000,
                       bio = sum(DRYBIO_AG * TPAGROW_UNADJ * tDI, na.rm = TRUE)) %>%
      dplyr::mutate(DIA_GROW = d / t,
                    BA_GROW = ba / t,
                    NETVOL_GROW = vol / t,
                    SAWVOL_GROW = svol / t,
                    BIO_GROW = bio / t,
                    BAA_GROW = ba,
                    NETVOL_GROW_AC = vol,
                    SAWVOL_GROW_AC = svol,
                    BIO_GROW_AC = bio,
                    PREV_TPA = t) %>%
      dplyr::select(-c(t:bio)) %>%
      as.data.frame() %>%
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
    ### Condition list
    a <- aData %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      #dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dplyr::mutate(fa = SUBPTYP_PROP_CHNG * aDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!aGrpSyms, PROP_BASIS, fa)


    grpSyms <- syms(grpBy)
    ## Tree list
    t <- data %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(tPlot = dplyr::case_when(ONEORTWO == 1 ~ TPAGROW_UNADJ * tDI,
                                             TRUE ~ 0),## Previous only
                    dPlot = DIA * TPAGROW_UNADJ * tDI,
                    bPlot = BA * TPAGROW_UNADJ * tDI,
                    gPlot = VOLCFNET * TPAGROW_UNADJ * tDI,
                    sPlot = VOLBFNET * TPAGROW_UNADJ * tDI / 1000,
                    bioPlot = DRYBIO_AG * TPAGROW_UNADJ * tDI / 2000) %>%
      ## Need a code that tells us where the tree was measured
      ## macroplot, microplot, subplot
      dplyr::mutate(TREE_BASIS = case_when(SUBPTYP_GRM == 0 ~ NA_character_,
                                           SUBPTYP_GRM == 1 ~ 'SUBP',
                                           SUBPTYP_GRM == 2 ~ 'MICR',
                                           SUBPTYP_GRM == 3 ~ 'MACR')) %>%
      dplyr::filter(!is.na(TREE_BASIS)) %>%
      dplyr::select(PLT_CN, TREE_BASIS, SUBP, TREE, ONEORTWO, !!!grpSyms, tPlot:bioPlot) %>%
      as.data.frame()

    ## Return a tree/condition list ready to be handed to `customPSE`
    if (treeList) {

      tEst <- a %>%
        dtplyr::lazy_dt() %>%
        dplyr::group_by(PLT_CN, CONDID, !!!aGrpSyms, AREA_BASIS) %>%
        dplyr::summarize(fa = sum(fa, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        as.data.frame() %>%
        dplyr::left_join(t, by = c('PLT_CN', aGrpBy)) %>%
        dplyr::mutate(EVAL_TYP = 'GROW') %>%
        dplyr::select(PLT_CN, EVAL_TYP, TREE_BASIS, AREA_BASIS,
                      !!!grpSyms, CONDID, SUBP, TREE, ONEORTWO,
                      DIA_GROW = dPlot,
                      BAA_GROW = bPlot,
                      NETVOL_GROW = gPlot,
                      SAWVOL_GROW = sPlot,
                      BIO_GROW = bioPlot,
                      PREV_TPA = tPlot,
                      PROP_FOREST = fa)
      out <- list(tEst = tEst, aEst = NULL, grpBy = grpBy, aGrpBy = aGrpBy)

    ## Otherwise, proceed to population estimation
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
      eu.sums <- sumToEU(db, dplyr::select(tPlt, -c(tPlot)), dplyr::select(tPlt, -c(dPlot:bioPlot)), pops, grpBy, grpBy, method)
      ttEst <- eu.sums$x %>%
        dplyr::select(ESTN_UNIT_CN, all_of(grpBy), dPlot_cv_t = dPlot_cv,
                      gPlot_cv_t = gPlot_cv, bPlot_cv_t = bPlot_cv,
                      sPlot_cv_t = sPlot_cv, bioPlot_cv_t = bioPlot_cv,)
      tEst <- dplyr::left_join(tEst, ttEst, by = c('ESTN_UNIT_CN', grpBy))

      out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy)
    }
  }


  return(out)

}


#' @export
vitalRates <- function(db,
                       grpBy = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       bySpecies = FALSE,
                       bySizeClass = FALSE,
                       landType = 'forest',
                       treeType = 'all',
                       method = 'TI',
                       lambda = .5,
                       treeDomain = NULL,
                       areaDomain = NULL,
                       totals = FALSE,
                       variance = FALSE,
                       byPlot = FALSE,
                       treeList = FALSE,
                       nCores = 1) {

  ##  don't have to change original code
  grpBy_quo <- rlang::enquo(grpBy)
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
  out <- lapply(X = iter, FUN = vrStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                bySpecies, bySizeClass,
                landType, treeType, method,
                lambda, treeDomain, areaDomain,
                totals, byPlot, treeList,
                nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
  aEst <- bind_rows(out[names(out) == 'aEst'])
  tEst <- bind_rows(out[names(out) == 'tEst'])
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
      dplyr::mutate(TREE_TOTAL = tPlot_mean,
                    DIA_TOTAL = dPlot_mean,
                    BA_TOTAL = bPlot_mean,
                    NETVOL_TOTAL = gPlot_mean,
                    SAWVOL_TOTAL = sPlot_mean,
                    BIO_TOTAL = bioPlot_mean,
                    AREA_TOTAL = fa_mean,

                    # Ratios
                    DIA_GROW = dPlot_mean / tPlot_mean,
                    BA_GROW = bPlot_mean / tPlot_mean,
                    NETVOL_GROW = gPlot_mean / tPlot_mean,
                    SAWVOL_GROW = sPlot_mean / tPlot_mean,
                    BIO_GROW = bioPlot_mean / tPlot_mean,
                    BA_GROW_AC = bPlot_mean / fa_mean,
                    NETVOL_GROW_AC = gPlot_mean / fa_mean,
                    SAWVOL_GROW_AC = sPlot_mean / fa_mean,
                    BIO_GROW_AC = bioPlot_mean / fa_mean,

                    # Variances
                    TREE_TOTAL_VAR = tPlot_var,
                    DIA_TOTAL_VAR = dPlot_var,
                    BA_TOTAL_VAR = bPlot_var,
                    NETVOL_TOTAL_VAR = gPlot_var,
                    SAWVOL_TOTAL_VAR = sPlot_var,
                    BIO_TOTAL_VAR = bioPlot_var,
                    AREA_TOTAL_VAR = fa_var,

                    DIA_GROW_VAR = ratioVar(dPlot_mean, tPlot_mean, dPlot_var, tPlot_var, dPlot_cv_t),
                    BA_GROW_VAR = ratioVar(bPlot_mean, tPlot_mean, bPlot_var, tPlot_var, bPlot_cv_t),
                    NETVOL_GROW_VAR = ratioVar(gPlot_mean, tPlot_mean, gPlot_var, tPlot_var, gPlot_cv_t),
                    SAWVOL_GROW_VAR = ratioVar(sPlot_mean, tPlot_mean, sPlot_var, tPlot_var, sPlot_cv_t),
                    BIO_GROW_VAR = ratioVar(bioPlot_mean, tPlot_mean, bioPlot_var, tPlot_var, bioPlot_cv_t),
                    BA_GROW_AC_VAR = ratioVar(bPlot_mean, fa_mean, bPlot_var, fa_var, bPlot_cv),
                    NETVOL_GROW_AC_VAR = ratioVar(gPlot_mean, fa_mean, gPlot_var, fa_var, gPlot_cv),
                    SAWVOL_GROW_AC_VAR = ratioVar(sPlot_mean, fa_mean, sPlot_var, fa_var, sPlot_cv),
                    BIO_GROW_AC_VAR = ratioVar(bioPlot_mean, fa_mean, bioPlot_var, fa_var, bioPlot_cv),

                    # Sampling Errors
                    TREE_TOTAL_SE = sqrt(tPlot_var) / tPlot_mean * 100,
                    DIA_TOTAL_SE = sqrt(dPlot_var) / abs(dPlot_mean) * 100,
                    BA_TOTAL_SE = sqrt(bPlot_var) / abs(bPlot_mean) * 100,
                    NETVOL_TOTAL_SE = sqrt(gPlot_var) / abs(gPlot_mean) * 100,
                    SAWVOL_TOTAL_SE = sqrt(sPlot_var) / abs(sPlot_mean) * 100,
                    BIO_TOTAL_SE = sqrt(bioPlot_var) / abs(bioPlot_mean) * 100,
                    AREA_TOTAL_SE = sqrt(fa_var) / fa_mean * 100,

                    DIA_GROW_SE = sqrt(DIA_GROW_VAR) / abs(DIA_GROW) * 100,
                    BA_GROW_SE = sqrt(BA_GROW_VAR) / abs(BA_GROW) * 100,
                    NETVOL_GROW_SE = sqrt(NETVOL_GROW_VAR) / abs(NETVOL_GROW) * 100,
                    SAWVOL_GROW_SE = sqrt(SAWVOL_GROW_VAR) / abs(SAWVOL_GROW) * 100,
                    BIO_GROW_SE = sqrt(BIO_GROW_VAR) / abs(BIO_GROW) * 100,
                    BA_GROW_AC_SE = sqrt(BA_GROW_AC_VAR) / abs(BA_GROW_AC) * 100,
                    NETVOL_GROW_AC_SE = sqrt(NETVOL_GROW_AC_VAR) / abs(NETVOL_GROW_AC) * 100,
                    SAWVOL_GROW_AC_SE = sqrt(SAWVOL_GROW_AC_VAR) / abs(SAWVOL_GROW_AC) * 100,
                    BIO_GROW_AC_SE = sqrt(BIO_GROW_AC_VAR) / abs(BIO_GROW_AC) * 100,


                    # Plot counts
                    nPlots_TREE = nPlots.x,
                    nPlots_AREA = nPlots.y,
                    N = P2PNTCNT_EU) %>%
      dplyr::select(!!!grpSyms,
                    DIA_GROW:BIO_GROW_AC,
                    DIA_TOTAL:BIO_TOTAL, TREE_TOTAL, AREA_TOTAL,
                    DIA_GROW_VAR:BIO_GROW_AC_VAR,
                    DIA_TOTAL_VAR:BIO_TOTAL_VAR, TREE_TOTAL_VAR, AREA_TOTAL_VAR,
                    DIA_GROW_SE:BIO_GROW_AC_SE,
                    DIA_TOTAL_SE:BIO_TOTAL_SE, TREE_TOTAL_SE, AREA_TOTAL_SE,
                    nPlots_TREE, nPlots_AREA, N)

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


