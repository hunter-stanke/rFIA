areaStarter <- function(x,
                        db,
                        grpBy_quo = NULL,
                        polys = NULL,
                        returnSpatial = FALSE,
                        byLandType = FALSE,
                        landType = 'forest',
                        method = 'TI',
                        lambda = .5,
                        treeDomain = NULL,
                        areaDomain = NULL,
                        totals = FALSE,
                        byPlot = FALSE,
                        condList = FALSE,
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
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest', 'water', 'non-forest', 'census water', 'non-census water', 'all') == FALSE){
    stop('landType must be one of: "forest", "timber", "non-forest", "water", or "all".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method, 'unknown. Defaulting to Temporally Indifferent (TI).'))
  }





  ## Prep other variables ------------------------------------------------------
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>%
    dplyr::mutate(PLT_CN = CN,
           pltID = stringr::str_c(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ## Convert grpBy to character
  grpBy <- grpByToChar(db, grpBy_quo)


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
  pops <- handlePops(db, evalType = c('CURR'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (str_to_upper(method) != 'TI') {
    pops <- mergeSmallStrata(db, pops)
  }



  ## Canned groups -------------------------------------------------------------
  # Make a new column that describes the land type and hold in COND
  if (byLandType){
    grpBy <- c(grpBy, 'landType')
    db$COND <- db$COND %>%
      dplyr::mutate(landType = dplyr::case_when(
        COND_STATUS_CD == 1 & SITECLCD %in% c(1:6) & RESERVCD == 0 ~ 'Timber',
        COND_STATUS_CD == 1 ~ 'Non-Timber Forest',
        COND_STATUS_CD == 2 ~ 'Non-Forest',
        COND_STATUS_CD == 3 | COND_STATUS_CD == 4 ~ 'Water'),
        landD = 1) # Reset the land basis to all
    db$COND <- db$COND[!is.na(db$COND$landType),]
  }







  ## Slim down the database before we hand it off to the estimators ------------
  ## Reduces memory requirements and speeds up processing ----------------------

  ## Only the necessary plots for EVAL of interest
  db$PLOT <- dplyr::filter(db$PLOT, PLT_CN %in% pops$PLT_CN)

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
                    SUBP, TREE, dplyr::all_of(grpT), tD)) %>%
    ## Drop plots outside our domain of interest
    dplyr::filter(!is.na(DIA) & TPA_UNADJ > 0 & tD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% db$PLOT$PLT_CN)

  ## was treeDomain NULL? If so, replace NAs w/ 1
  treeD <- ifelse(mean(db$TREE$tD, na.rm = TRUE) == 1, 1, 0)



  ## Summarize variables of interest to plot level ------------------------------
  ## Char to syms
  grpSyms <- dplyr::syms(grpBy)

  ### Only joining tables necessary to produce plot level estimates,
  ## adjusted for non-response
  data <- db$PLOT %>%
    dplyr::left_join(db$COND, by = c('PLT_CN')) %>%
    dplyr::left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
    dtplyr::lazy_dt() %>%
    dplyr::mutate(tD = tidyr::replace_na(tD, treeD)) %>%
    dplyr::group_by(PLT_CN, PROP_BASIS, CONDID, !!!grpSyms) %>%
    dplyr::mutate(tD = ifelse(sum(tD, na.rm = TRUE) > 0, 1, 0)) %>%
    dplyr::ungroup() %>%
    as.data.frame()


  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp * data$tD
  data$pDI <- data$landD * data$aD

  if (byPlot & !condList){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- dplyr::syms(grpBy)

    tEst <- data %>%
      dplyr::mutate(YEAR = INVYR) %>%
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!grpSyms, PLT_CN, CONDID) %>%
      dplyr::summarize(CONDPROP_UNADJ = data.table::first(CONDPROP_UNADJ),
                aDI = data.table::first(aDI)) %>%
      dplyr::group_by(!!!grpSyms, PLT_CN) %>%
      dplyr::summarize(PROP_FOREST = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
      as.data.frame()

    ## Make it spatial
    if (returnSpatial){
      tEst <- tEst %>%
        dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
        sf::st_as_sf(coords = c('LON', 'LAT'),
                     crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }

    out <- list(tEst = tEst, aEst = NULL, grpBy = grpBy)

  } else {

    grpSyms <- syms(grpBy)

    ### Plot-level estimates
    t <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at stratum level
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(fa = CONDPROP_UNADJ * aDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!grpSyms, fa) %>%
      as.data.frame()

    ## Total land area in areaDomain and landType, for proportions
    a <- data %>%
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(fad = CONDPROP_UNADJ * pDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, fad) %>%
      as.data.frame()

    ## If any grpBy are NA in t, then those plots need to be NA in a as well
    ## to simplify interpretation of the proportions
    naPlts <- t %>%
      dplyr::distinct(PLT_CN, !!!grpSyms) %>%
      tidyr::drop_na() %>%
      dplyr::mutate(good = 1) %>%
      dplyr::distinct(PLT_CN, good)
    a <- a %>%
      dplyr::left_join(naPlts, by = 'PLT_CN') %>%
      ## Make fa NA when grps are NA
      dplyr::mutate(fad = dplyr::case_when(good == 1 ~ fad,
                                   TRUE ~ NA_real_)) %>%
      dplyr::select(-c(good))


    ## Return a tree/condition list ready to be handed to `customPSE`
    if (condList) {

      tEst <- t %>%
        dplyr::mutate(EVAL_TYP = 'CURR') %>%
        dplyr::select(PLT_CN, EVAL_TYP, AREA_BASIS, !!!grpSyms, CONDID,
                      PROP_FOREST = fa)
      out <- list(tEst = tEst, aEst = NULL, grpBy = grpBy)


    ## Otherwise, proceed to population estimation
    } else {

      ## Sum variable(s) up to plot-level and adjust for non-response
      tPlt <- sumToPlot(t, pops, grpBy)
      aPlt <- sumToPlot(a, pops, NULL)

      ## Adding YEAR to groups
      grpBy <- c('YEAR', grpBy)
      aGrpBy <- c('YEAR')


      ## Sum variable(s) up to strata then estimation unit level
      eu.sums <- sumToEU(db, tPlt, aPlt, pops, grpBy, aGrpBy, method)
      tEst <- eu.sums$x
      aEst <- eu.sums$y

      out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy)
    }
  }

  return(out)

}




#'@export
area <- function(db,
                 grpBy = NULL,
                 polys = NULL,
                 returnSpatial = FALSE,
                 byLandType = FALSE,
                 landType = 'forest',
                 method = 'TI',
                 lambda = .5,
                 treeDomain = NULL,
                 areaDomain = NULL,
                 totals = FALSE,
                 variance = FALSE,
                 byPlot = FALSE,
                 condList = FALSE,
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
  out <- lapply(X = iter, FUN = areaStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                byLandType, landType, method,
                lambda, treeDomain, areaDomain,
                totals, byPlot, condList,
                nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
  tEst <- bind_rows(out[names(out) == 'tEst'])
  aEst <- bind_rows(out[names(out) == 'aEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  grpSyms <- dplyr::syms(grpBy)




  ## Summarize population estimates across estimation units
  if (!byPlot & !condList){

    ## Combine most-recent population estimates across states with potentially
    ## different reporting schedules, i.e., if 2016 is most recent in MI and 2017 is
    ## most recent in WI, combine them and label as 2017
    if (mr) {
      tEst <- combineMR(tEst, grpBy)
    }


    ## Totals and ratios -------------------------------------------------------
    tEst <- tEst %>%
      dplyr::group_by(!!!grpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE))
    aEst <- aEst %>%
      dplyr::group_by(YEAR) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::select(YEAR, fad_mean, fad_var)


    suppressWarnings({
      ## Bring them together
      tEst <- tEst %>%
        dplyr::left_join(aEst, by = 'YEAR') %>%
        # Renaming, computing ratios, and SE
        dplyr::mutate(PERC_AREA = fa_mean / fad_mean,
                      AREA_TOTAL = fa_mean,
                      AREA_TOTAL_SE = sqrt(fa_var) / AREA_TOTAL *100,

                      ## ratio variance
                      rVar = ratioVar(fa_mean, fad_mean, fa_var, fad_var, fa_cv),

                      ## Convert to percentage
                      PERC_AREA = PERC_AREA * 100,
                      rVar = rVar * (100^2),

                      ## Ratio variances
                      # These aren't truly negative values, but come from rounding errors
                      # when PERC_AREA = 100, i.e., estimated variance is 0
                      PERC_AREA_SE = dplyr::case_when(rVar < 0 ~ 0,
                                                      TRUE ~ sqrt(rVar) / PERC_AREA * 100),
                      PERC_AREA_VAR = dplyr::case_when(rVar < 0 ~ 0,
                                                       TRUE ~ rVar),

                      AREA_TOTAL_VAR = fa_var,
                      nPlots_AREA = nPlots.x,
                      N = P2PNTCNT_EU) %>%
        dplyr::select(grpBy, PERC_AREA, AREA_TOTAL, PERC_AREA_SE, AREA_TOTAL_SE,
                      PERC_AREA_VAR, AREA_TOTAL_VAR, nPlots_AREA, N)
    })



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

