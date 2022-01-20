standStructStarter <- function(x,
                               db,
                               grpBy_quo = NULL,
                               polys = NULL,
                               returnSpatial = FALSE,
                               landType = 'forest',
                               method = 'TI',
                               lambda = .5,
                               areaDomain = NULL,
                               totals = FALSE,
                               byPlot = FALSE,
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
  grpBy <- grpByToChar(db[!c(names(db) %in% 'TREE')], grpBy_quo)


  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}


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

  ## Spatial boundary
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(!is.na(db$PLOT$polyID), 1, 0)
  } else {
    db$PLOT$sp <- 1
  }

  # User defined domain indicator for area (ex. specific forest type)
  db <- udAreaDomain(db, areaDomain)




  ## Handle population tables --------------------------------------------------
  ## Filtering out all inventories that are not relevant to the current estimation
  ## type. If using estimator other than TI, handle the differences in P2POINTCNT
  ## and in assigning YEAR column (YEAR = END_INVYR if method = 'TI')
  pops <- handlePops(db, evalType = c('CURR'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (str_to_upper(method) != 'TI') {
    pops <- mergeSmallStrata(db, pops)
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
                    SUBP, TREE, dplyr::all_of(grpT), CCLCD)) %>%
    ## Drop plots outside our domain of interest
    dplyr::filter(!is.na(DIA) & TPA_UNADJ > 0) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% db$PLOT$PLT_CN)


  ## Full tree list
  data <- db$PLOT %>%
    dplyr::left_join(db$COND, by = c('PLT_CN')) %>%
    dplyr::left_join(db$TREE, by = c('PLT_CN', 'CONDID'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp



  ## Plot-level summaries ------------------------------------------------------
  if (byPlot){

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
      dplyr::group_by(!!!grpSyms, PLT_CN, CONDID, CONDPROP_UNADJ, aDI) %>%
      dplyr::summarize(STAGE = structHelper(DIA, CCLCD)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(STAGE = stringr::str_to_upper(STAGE)) %>%
      dplyr::group_by(!!!grpSyms, PLT_CN, STAGE) %>%
      dplyr::summarize(PROP_STAGE = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
      as.data.frame() %>%
      dplyr::distinct() %>%
      dplyr::left_join(a, by = c('PLT_CN', grpBy))


    ## Make it spatial
    if (returnSpatial){
      t <- t %>%
        dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
        sf::st_as_sf(coords = c('LON', 'LAT'),
                     crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }

    out <- list(tEst = t, aEst = NULL, grpBy = grpBy, aGrpBy = NULL)




    ## Population estimation -----------------------------------------------------
  } else {

    grpSyms <- syms(grpBy)
    ### Condition list
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dplyr::mutate(fa = CONDPROP_UNADJ * aDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!grpSyms, PROP_BASIS, fa)


    ## condition stage list
    t <- data %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!grpSyms, PLT_CN, PROP_BASIS, CONDID, CONDPROP_UNADJ, aDI) %>%
      dplyr::summarize(STAGE = structHelper(DIA, CCLCD)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(STAGE = stringr::str_to_upper(STAGE),
                    cover = CONDPROP_UNADJ * aDI) %>%
      dplyr::rename(AREA_BASIS = PROP_BASIS) %>%
      as.data.frame()

    ## Adding stage to grpBy
    aGrpBy <- grpBy
    grpBy <- c(grpBy, 'STAGE')
    grpSyms <- dplyr::syms(grpBy)
    aGrpSyms <- dplyr::syms(aGrpBy)


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

    out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy)
  }


  return(out)

}


#' @export
standStruct <- function(db,
                        grpBy = NULL,
                        polys = NULL,
                        returnSpatial = FALSE,
                        landType = 'forest',
                        method = 'TI',
                        lambda = .5,
                        areaDomain = NULL,
                        totals = FALSE,
                        variance = FALSE,
                        byPlot = FALSE,
                        nCores = 1) {

  ##  don't have to change original code
  grpBy_quo <- rlang::enquo(grpBy)
  areaDomain <- rlang::enquo(areaDomain)

  ## Handle iterator if db is remote
  remote <- ifelse(class(db) == 'Remote.FIA.Database', 1, 0)
  iter <- remoteIter(db, remote)

  ## Check for a most recent subset
  mr <- checkMR(db, remote)

  ## prep for areal summary
  polys <- arealSumPrep1(polys)


  ## Run the main portion
  out <- lapply(X = iter, FUN = standStructStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                landType, method,
                lambda, areaDomain,
                totals, byPlot, nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
  tEst <- bind_rows(out[names(out) == 'tEst'])
  aEst <- bind_rows(out[names(out) == 'aEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  aGrpBy <- out[names(out) == 'aGrpBy'][[1]]
  grpSyms <- dplyr::syms(grpBy)
  aGrpSyms <- dplyr::syms(aGrpBy)


  ## Summarize population estimates across estimation units
  if (!byPlot){

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
      dplyr::mutate(STAGE_AREA_TOTAL = cover_mean,
                    AREA_TOTAL = fa_mean,
                    # Ratios
                    COVER_PCT = cover_mean / fa_mean,
                    # Variances
                    STAGE_AREA_TOTAL_VAR = cover_var,
                    AREA_TOTAL_VAR = fa_var,
                    COVER_PCT_VAR = ratioVar(cover_mean, fa_mean, cover_var, fa_var, cover_cv),
                    COVER_PCT = COVER_PCT * 100,
                    COVER_PCT_VAR = COVER_PCT_VAR * (100^2),
                    # Sampling Errors
                    STAGE_AREA_TOTAL_SE = sqrt(cover_var) / cover_mean * 100,
                    AREA_TOTAL_SE = sqrt(fa_var) / fa_mean * 100,
                    COVER_PCT_SE = sqrt(COVER_PCT_VAR) / COVER_PCT * 100,
                    # Plot counts
                    nPlots_AREA = nPlots.x,
                    N = P2PNTCNT_EU) %>%
      dplyr::select(!!!grpSyms, COVER_PCT, STAGE_AREA_TOTAL, AREA_TOTAL,
                    COVER_PCT_VAR, STAGE_AREA_TOTAL_VAR, AREA_TOTAL_VAR,
                    COVER_PCT_SE, STAGE_AREA_TOTAL_SE, AREA_TOTAL_SE,
                    nPlots_AREA, N)

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
    tidyr::drop_na(grpBy[!c(grpBy %in% names(polys))]) %>%
    dplyr::arrange(YEAR) %>%
    as_tibble()



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

