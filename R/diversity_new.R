diversityStarter <- function(x,
                             db,
                             grpBy_quo = NULL,
                             polys = NULL,
                             returnSpatial = FALSE,
                             bySizeClass = FALSE,
                             landType = 'forest',
                             treeType = 'live',
                             method = 'TI',
                             lambda = .5,
                             stateVar = TPA_UNADJ,
                             grpVar = SPCD,
                             treeDomain = NULL,
                             areaDomain = NULL,
                             byPlot = FALSE,
                             condList = FALSE,
                             totals = FALSE,
                             nCores = 1,
                             remote,
                             mr){


  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN',
                 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')

  ## If remote, read in state by state. Otherwise, drop all unnecessary tables
  db <- readRemoteHelper(x, db, remote, reqTables, nCores)

  ## IF the object was clipped
  if ('prev' %in% names(db$PLOT)){
    ## Only want the current plots, no grm
    db$PLOT <- dplyr::filter(db$PLOT, prev == 0)
  }

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
  if (treeType %in% c('live', 'dead', 'gs', 'all') == FALSE){
    stop('treeType must be one of: "live", "dead", "gs", or "all".')
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


  ## Prep other variables ------------------------------------------------------
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>%
    dplyr::mutate(PLT_CN = CN,
                  pltID = stringr::str_c(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ## Convert grpBy to character
  grpBy <- grpByToChar(db, grpBy_quo)

  ## Handle the modifier if it was given
  db$TREE <- db$TREE %>%
    dplyr::mutate(TRE_CN = CN,
                  state = !!stateVar,
                  grp = !!grpVar)

  # I like a unique ID for a plot through time
  if (byPlot | condList) {grpBy <- c('pltID', grpBy)}


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
  ## Land type
  db$COND$landD <- landTypeDomain(landType,
                                  db$COND$COND_STATUS_CD,
                                  db$COND$SITECLCD,
                                  db$COND$RESERVCD)
  ## Tree type
  db$TREE$typeD <- treeTypeDomain(treeType,
                                  db$TREE$STATUSCD,
                                  db$TREE$DIA,
                                  db$TREE$TREECLCD)

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
  pops <- handlePops(db, evalType = c('VOL'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (stringr::str_to_upper(method) != 'TI') {
    pops <- mergeSmallStrata(db, pops)
  }



  ## Canned groups -------------------------------------------------------------

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
                    dplyr::all_of(grpP), sp, COUNTYCD)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    dplyr::filter(PLOT_STATUS_CD == 1 & sp == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)

  db$COND <- db$COND %>%
    dplyr::select(c(PLT_CN, CONDPROP_UNADJ, PROP_BASIS,
                    COND_STATUS_CD, CONDID,
                    dplyr::all_of(grpC), aD, landD)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    dplyr::filter(aD == 1 & landD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)

  db$TREE <- db$TREE %>%
    dplyr::select(c(PLT_CN, CONDID, DIA, SPCD, TPA_UNADJ,
                    SUBP, TREE, dplyr::all_of(grpT), tD, typeD,
                    grp, state)) %>%
    ## Drop plots outside our domain of interest
    dplyr::filter(!is.na(DIA) & TPA_UNADJ > 0 & tD == 1 & typeD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% db$PLOT$PLT_CN)


  ### Full tree list
  data <- db$PLOT %>%
    dplyr::left_join(db$COND, by = c('PLT_CN')) %>%
    dplyr::left_join(db$TREE, by = c('PLT_CN', 'CONDID'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp
  data$tDI <- data$landD * data$aD * data$tD * data$typeD * data$sp




  ## Plot-level summaries ------------------------------------------------------
  if (byPlot & !condList){

    grpBy <- c('YEAR', grpBy)
    grpSyms <- syms(grpBy)

    ### Plot-level estimates
    a <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(PLT_CN, !!!grpSyms) %>%
      dplyr::summarize(PROP_FOREST = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
      as.data.frame()

    t <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!grpSyms, PLT_CN) %>%
      dplyr::summarize(H = divIndex(grp, state  * tDI, index = 'H'),
                       S = divIndex(grp, state * tDI, index = 'S'),
                       Eh = divIndex(grp, state * tDI, index = 'Eh')) %>%
      as.data.frame() %>%
      dplyr::left_join(a, by = c('PLT_CN', grpBy)) %>%
      dplyr::distinct()

    ## Make it spatial
    if (returnSpatial){
      t <- t %>%
        dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
        sf::st_as_sf(coords = c('LON', 'LAT'),
                     crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }

    out <- list(tEst = t, grpBy = grpBy)




    ## Population estimation -----------------------------------------------------
  } else {

    grpSyms <- syms(grpBy)

    ## Condition list
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dplyr::mutate(fa = CONDPROP_UNADJ * aDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!grpSyms, PROP_BASIS, fa)


    ## Tree list --> condition list w/ diversity indices
    t <- data %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      dplyr::rename(AREA_BASIS = PROP_BASIS) %>%
      dplyr::group_by(PLT_CN, CONDID, !!!grpSyms, CONDPROP_UNADJ, aDI, AREA_BASIS) %>%
      dplyr::summarize(H = divIndex(grp, state  * tDI, index = 'H'),
                       S = divIndex(grp, state * tDI, index = 'S'),
                       Eh = divIndex(grp, state * tDI, index = 'Eh')) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(across(H:Eh, .fns = ~ .x * CONDPROP_UNADJ)) %>%
      as.data.frame()


    ## Return a tree/condition list ready to be handed to `customPSE`
    if (condList) {

      tEst <- a %>%
        dplyr::left_join(t, by = c('PLT_CN', 'CONDID', 'AREA_BASIS', grpBy)) %>%
        dplyr::mutate(EVAL_TYP = 'VOL') %>%
        dplyr::select(PLT_CN, EVAL_TYP, AREA_BASIS,
                      !!!grpSyms, CONDID,
                      H:Eh,
                      PROP_FOREST = fa)

      out <- list(tEst = tEst, aEst = NULL, grpBy = grpBy, full = NULL)

    ## Otherwise, proceed to population estimation
    } else {

      ## Sum variable(s) up to plot-level and adjust for non-response
      tPlt <- sumToPlot(t, pops, grpBy)
      aPlt <- sumToPlot(a, pops, grpBy)

      ## Adding YEAR to groups
      grpBy <- c('YEAR', grpBy)


      ## Sum variable(s) up to strata then estimation unit level
      eu.sums <- sumToEU(db, tPlt, aPlt, pops, grpBy, grpBy, method)
      tEst <- eu.sums$x
      aEst <- eu.sums$y

      ## Using this to return a tree list for gamma and beta
      full <- data %>%
        dplyr::mutate(state = state * tDI) %>%
        dplyr::distinct(PLT_CN, !!!grpSyms, SUBP, TREE, grp, state) %>%
        dplyr::inner_join(dplyr::select(pops, c(YEAR, PLT_CN)), by = 'PLT_CN') %>%
        dplyr::filter(!is.na(YEAR) & !is.na(state) & !is.na(grp))

      out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, full = full)

    }
  }





  return(out)

}


#' @export
diversity <- function(db,
                      grpBy = NULL,
                      polys = NULL,
                      returnSpatial = FALSE,
                      bySizeClass = FALSE,
                      landType = 'forest',
                      treeType = 'live',
                      method = 'TI',
                      lambda = .5,
                      stateVar = TPA_UNADJ,
                      grpVar = SPCD,
                      treeDomain = NULL,
                      areaDomain = NULL,
                      byPlot = FALSE,
                      condList = FALSE,
                      totals = FALSE,
                      variance = FALSE,
                      nCores = 1) {

  ##  don't have to change original code
  grpBy_quo <- rlang::enquo(grpBy)
  areaDomain <- rlang::enquo(areaDomain)
  treeDomain <- rlang::enquo(treeDomain)
  stateVar <- rlang::enquo(stateVar)
  grpVar <- rlang::enquo(grpVar)


  ## Handle iterator if db is remote
  remote <- ifelse(class(db) == 'Remote.FIA.Database', 1, 0)
  iter <- remoteIter(db, remote)

  ## Check for a most recent subset
  mr <- checkMR(db, remote)

  ## prep for areal summary
  polys <- arealSumPrep1(polys)



  ## Run the main portion
  out <- lapply(X = iter, FUN = diversityStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                bySizeClass,
                landType, treeType, method,
                lambda, stateVar, grpVar,
                treeDomain, areaDomain,
                byPlot, condList,
                totals, nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
  tEst <- dplyr::bind_rows(out[names(out) == 'tEst'])
  aEst <- dplyr::bind_rows(out[names(out) == 'aEst'])
  full <- dplyr::bind_rows(out[names(out) == 'full'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  grpSyms <- dplyr::syms(grpBy)


  ## Summarize population estimates across estimation units
  if (!byPlot & !condList){

    ## Combine most-recent population estimates across states with potentially
    ## different reporting schedules, i.e., if 2016 is most recent in MI and 2017 is
    ## most recent in WI, combine them and label as 2017
    if (mr) {
      tEst <- combineMR(tEst, grpBy)
      aEst <- combineMR(aEst, aGrpBy)
    }



    ## Totals and ratios -------------------------------------------------------
    aEst <- aEst %>%
      dplyr::group_by( !!!grpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::select(!!!grpSyms, fa_mean, fa_var, nPlots.y)


    tEst <- tEst %>%
      dplyr::group_by(!!!grpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::left_join(aEst, by = grpBy) %>%
      dplyr::mutate(H_a = H_mean / fa_mean,
                    Eh_a = Eh_mean / fa_mean,
                    S_a = S_mean / fa_mean,

                    AREA_TOTAL = fa_mean,

                    # Variances
                    AREA_TOTAL_VAR = fa_var,
                    H_a_VAR = ratioVar(H_mean, fa_mean, H_var, fa_var, H_cv),
                    Eh_a_VAR = ratioVar(Eh_mean, fa_mean, Eh_var, fa_var, Eh_cv),
                    S_a_VAR = ratioVar(S_mean, fa_mean, S_var, fa_var, S_cv),

                    # Sampling Errors
                    H_a_SE = sqrt(H_a_VAR) / H_a * 100,
                    Eh_a_SE = sqrt(Eh_a_VAR) / Eh_a * 100,
                    S_a_SE = sqrt(S_a_VAR) / S_a * 100,
                    AREA_TOTAL_SE = sqrt(fa_var) / fa_mean * 100,

                    # Plot counts
                    nPlots_TREE = nPlots.x,
                    nPlots_AREA = nPlots.y,
                    N = P2PNTCNT_EU)

    ## Up a few spatial scales
    fullGrps <- dplyr::syms(grpBy[!c(grpBy %in% 'lambda')])
    full <- full %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!fullGrps) %>%
      dplyr::summarize(H_g = divIndex(grp, state, index = 'H'),
                       Eh_g = divIndex(grp, state, index = 'Eh'),
                       S_g = divIndex(grp, state, index = 'S')) %>%
      as.data.frame()
    tEst <- tEst %>%
      dplyr::left_join(full, by = grpBy[!c(grpBy %in% 'lambda')]) %>%
      dplyr::mutate(H_b = H_g - H_a,
                    Eh_b = Eh_g - Eh_a,
                    S_b = S_g - S_a) %>%
      dplyr::select(!!!grpSyms,
                    H_a, Eh_a, S_a,
                    H_b, Eh_b, S_b,
                    H_g, Eh_g, S_g,
                    AREA_TOTAL,
                    H_a_VAR, Eh_a_VAR, S_a_VAR, AREA_TOTAL_VAR,
                    H_a_SE, Eh_a_SE, S_a_SE, AREA_TOTAL_SE,
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
  if (!condList) {
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



