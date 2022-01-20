carbonStarter <- function(x,
                          db,
                          grpBy_quo = NULL,
                          polys = NULL,
                          returnSpatial = FALSE,
                          byPool = TRUE,
                          byComponent = FALSE,
                          modelSnag = TRUE,
                          landType = 'forest',
                          method = 'TI',
                          lambda = .5,
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
  if (!is.null(polys) &
      dplyr::first(class(polys)) %in%
      c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest', 'all') == FALSE){
    stop('landType must be one of: "forest", "timber", or "all".')
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
  grpBy <- grpByToChar(db[names(db) != 'TREE'], grpBy_quo)


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





  ## Prep the tree list --------------------------------------------------------
  ## Narrow up the tables to the necessary variables
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy &
                           !c(names(db$COND) %in% grpP)]


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
                    dplyr::all_of(grpC), aD, landD,
                    CARBON_DOWN_DEAD, CARBON_LITTER,
                    CARBON_SOIL_ORG, CARBON_STANDING_DEAD,
                    CARBON_UNDERSTORY_AG, CARBON_UNDERSTORY_BG)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    dplyr::filter(aD == 1 & landD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)

  db$TREE <- db$TREE %>%
    dplyr::select(c(PLT_CN, CONDID, DIA, SPCD, TPA_UNADJ,
                    SUBP, TREE, STATUSCD,
                    CARBON_AG, CARBON_BG)) %>%
    ## Drop plots outside our domain of interest
    dplyr::filter(!is.na(DIA) & TPA_UNADJ > 0) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% db$PLOT$PLT_CN)



  ## Full tree list
  data <- db$PLOT %>%
    dplyr::left_join(db$COND, by = c('PLT_CN')) %>%
    dplyr::left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
    dplyr::mutate(live = dplyr::case_when(STATUSCD == 1 ~ 1,
                                   is.na(DIA) ~ NA_real_,
                                   TRUE ~ 0),
           dead = dplyr::case_when(STATUSCD == 2 ~ 1,
                                   is.na(DIA) ~ NA_real_,
                                   TRUE ~ 0))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp


  ## Plot-level summaries ------------------------------------------------------
  if (byPlot & !condList){

    grpBy <- c('YEAR', grpBy)
    grpSyms <- dplyr::syms(grpBy)

    ### Condition-level estimates
    a <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(PLT_CN, !!!grpSyms) %>%
      dplyr::summarize(AG_UNDER_LIVE = sum(CONDPROP_UNADJ * CARBON_UNDERSTORY_AG * aDI, na.rm = TRUE),
                       BG_UNDER_LIVE = sum(CONDPROP_UNADJ * CARBON_UNDERSTORY_BG * aDI, na.rm = TRUE),
                       DOWN_DEAD = sum(CONDPROP_UNADJ * CARBON_DOWN_DEAD * aDI, na.rm = TRUE),
                       STAND_DEAD_MOD = sum(CONDPROP_UNADJ * CARBON_STANDING_DEAD * aDI, na.rm = TRUE),
                       LITTER = sum(CONDPROP_UNADJ * CARBON_LITTER * aDI, na.rm = TRUE),
                       SOIL_ORG = sum(CONDPROP_UNADJ * CARBON_SOIL_ORG * aDI, na.rm = TRUE),
                       PROP_FOREST = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
      as.data.frame()

    t <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!grpSyms, PLT_CN) %>%
      dplyr::summarize(AG_OVER_LIVE = sum(CARBON_AG * TPA_UNADJ * live * aDI, na.rm = TRUE) / 2000,
                       BG_OVER_LIVE = sum(CARBON_BG * TPA_UNADJ * live * aDI, na.rm = TRUE) / 2000,
                       AG_OVER_DEAD = sum(CARBON_AG * TPA_UNADJ * dead * aDI, na.rm = TRUE) / 2000,
                       BG_OVER_DEAD = sum(CARBON_BG * TPA_UNADJ * dead * aDI, na.rm = TRUE) / 2000) %>%
      as.data.frame()

    ## Join these back together
    t <- a %>%
      dplyr::left_join(t, by = c('PLT_CN', grpBy))

    ## Decide which estimate to use for snags
    if (modelSnag){
      t <- t %>%
        dplyr::mutate(STAND_DEAD = STAND_DEAD_MOD) %>%
        dplyr::select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))

    } else {
      t <- t %>%
        dplyr::mutate(STAND_DEAD = AG_OVER_DEAD + BG_OVER_DEAD) %>%
        dplyr::select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))
    }


    ## Convert to long format, where rows are ecosystem components
    t <- t %>%
      tidyr::pivot_longer(cols = -c(PLT_CN, !!!grpSyms, PROP_FOREST),
                          names_to = 'COMPONENT',
                          values_to = 'CARB_ACRE') %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(POOL = dplyr::case_when(COMPONENT %in% c('AG_UNDER_LIVE', 'AG_OVER_LIVE') ~ 'AG_LIVE',
                                            COMPONENT %in% c('BG_UNDER_LIVE', 'BG_OVER_LIVE') ~ 'BG_LIVE',
                                            COMPONENT %in% c('STAND_DEAD', 'DOWN_DEAD') ~ 'DEAD_WOOD',
                                            TRUE ~ COMPONENT)) %>%
      dplyr::mutate(CARB_ACRE = CARB_ACRE * 0.90718474) %>%

      as.data.frame()


    ## Add either component or pool to grpBy, depending on user input
    if (byComponent) grpBy <- c(grpBy, 'COMPONENT')
    if (byPool) grpBy <- c(grpBy, 'POOL')
    grpSyms <- dplyr::syms(grpBy)

    ## Summarize across COMPONENTS, if necessary
    if (!byComponent) {
      t <- t %>%
        dtplyr::lazy_dt() %>%
        dplyr::group_by(PLT_CN, !!!grpSyms, PROP_FOREST) %>%
        dplyr::summarise(CARB_ACRE = sum(CARB_ACRE, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::relocate(PROP_FOREST, .after = CARB_ACRE) %>%
        as.data.frame()
    } else {
      t <- t %>%
        dplyr::select(PLT_CN, !!!grpSyms, CARB_ACRE, PROP_FOREST)
    }



    ## Make it spatial
    if (returnSpatial){
      t <- t %>%
        dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
        sf::st_as_sf(coords = c('LON', 'LAT'),
                     crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }

    out <- list(tEst = t, grpBy = grpBy, grpBy = grpBy, aGrpBy = NULL)




  ## Population estimation ---------------------------------------------------
  } else {

    grpSyms <- dplyr::syms(grpBy)

    ### Condition-level estimates
    a <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(AG_UNDER_LIVE = CONDPROP_UNADJ * CARBON_UNDERSTORY_AG * aDI,
                       BG_UNDER_LIVE = CONDPROP_UNADJ * CARBON_UNDERSTORY_BG * aDI,
                       DOWN_DEAD = CONDPROP_UNADJ * CARBON_DOWN_DEAD * aDI,
                       STAND_DEAD_MOD = CONDPROP_UNADJ * CARBON_STANDING_DEAD * aDI,
                       LITTER = CONDPROP_UNADJ * CARBON_LITTER * aDI,
                    SOIL_ORG = CONDPROP_UNADJ * CARBON_SOIL_ORG * aDI,
                    PROP_FOREST = CONDPROP_UNADJ * aDI) %>%
      dplyr::select(PLT_CN, CONDID, !!!grpSyms, AREA_BASIS = PROP_BASIS, AG_UNDER_LIVE:PROP_FOREST) %>%
      as.data.frame()

    t <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(AG_OVER_LIVE = CARBON_AG * TPA_UNADJ * live * aDI / 2000,
                    BG_OVER_LIVE = CARBON_BG * TPA_UNADJ * live * aDI / 2000,
                    AG_OVER_DEAD = CARBON_AG * TPA_UNADJ * dead * aDI / 2000,
                    BG_OVER_DEAD = CARBON_BG * TPA_UNADJ * dead * aDI / 2000) %>%
      ## Need a code that tells us where the tree was measured
      ## macroplot, microplot, subplot
      dplyr::mutate(
        TREE_BASIS = dplyr::case_when(
          ## When DIA is na, adjustment is NA
          is.na(DIA) ~ NA_character_,
          ## When DIA is less than 5", use microplot value
          DIA < 5 ~ 'MICR',
          ## When DIA is greater than 5", use subplot value
          DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
          DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
          DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) %>%
      dplyr::filter(!is.na(TREE_BASIS)) %>%
      dplyr::select(!!!grpSyms, PLT_CN, CONDID, TREE_BASIS,
                    AG_OVER_LIVE, BG_OVER_LIVE, AG_OVER_DEAD, BG_OVER_DEAD) %>%
      as.data.frame()


    ## Return a tree/condition list ready to be handed to `customPSE`
    if (condList) {

      ## Sum tree biomass up to conditions
      t <- t %>%
        dtplyr::lazy_dt() %>%
        dplyr::group_by(PLT_CN, CONDID, !!!grpSyms) %>%
        dplyr::summarise(dplyr::across(AG_OVER_LIVE:BG_OVER_DEAD, sum, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        as.data.frame()

      t <- a %>%
        dplyr::left_join(t, by = c('PLT_CN', grpBy, 'CONDID'))


      ## Decide which estimate to use for snags
      if (modelSnag){
        t <- t %>%
          dplyr::mutate(STAND_DEAD = STAND_DEAD_MOD) %>%
          dplyr::select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))

      } else {
        t <- t %>%
          dplyr::mutate(STAND_DEAD = AG_OVER_DEAD + BG_OVER_DEAD) %>%
          dplyr::select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))
      }

      ## Convert to long format, where rows are ecosystem components
      t <- t %>%
        tidyr::pivot_longer(cols = -c(PLT_CN, AREA_BASIS, CONDID, PROP_FOREST, !!!grpSyms),
                            names_to = 'COMPONENT',
                            values_to = 'CARB_ACRE') %>%
        dtplyr::lazy_dt() %>%
        dplyr::mutate(POOL = dplyr::case_when(COMPONENT %in% c('AG_UNDER_LIVE', 'AG_OVER_LIVE') ~ 'AG_LIVE',
                                              COMPONENT %in% c('BG_UNDER_LIVE', 'BG_OVER_LIVE') ~ 'BG_LIVE',
                                              COMPONENT %in% c('STAND_DEAD', 'DOWN_DEAD') ~ 'DEAD_WOOD',
                                              TRUE ~ COMPONENT)) %>%
        dplyr::mutate(CARB_ACRE = CARB_ACRE * 0.90718474) %>%
        as.data.frame()

      ## Add either component or pool to grpBy, depending on user input
      if (byComponent) grpBy <- c(grpBy, 'POOL', 'COMPONENT')
      if (byPool & !byComponent) grpBy <- c(grpBy, 'POOL')
      grpSyms <- dplyr::syms(grpBy)
      aGrpBy <- grpBy[!c(grpBy %in% c('POOL', 'COMPONENT'))]
      aGrpSyms <- dplyr::syms(aGrpBy)

      ## Summarize across COMPONENTS, if necessary
      if (!byComponent) {
        t <- t %>%
          dtplyr::lazy_dt() %>%
          dplyr::group_by(PLT_CN, CONDID, !!!grpSyms, PROP_FOREST, AREA_BASIS) %>%
          dplyr::summarise(CARB_ACRE = sum(CARB_ACRE, na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          as.data.frame()
      }

      ## Reorder variable names
      t <- t %>%
        dplyr::mutate(EVAL_TYP = 'VOL') %>%
        dplyr::select(PLT_CN, EVAL_TYP, AREA_BASIS,
                      CONDID, !!!grpSyms, CARB_ACRE,
                      PROP_FOREST)


      out <- list(tEst = t, aEst = NULL, grpBy = grpBy, aGrpBy = aGrpBy)


    ## Otherwise, proceed to population estimation
    } else {

      ## Sum variable(s) up to plot-level and adjust for non-response
      tPlt <- sumToPlot(t, pops, grpBy)
      aPlt <- sumToPlot(a, pops, grpBy)
      tPlt <- aPlt %>%
        dplyr::select(-c(PROP_FOREST)) %>%
        dplyr::left_join(tPlt, by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'PLT_CN', grpBy))


      ## Decide which estimate to use for snags
      if (modelSnag){
        tPlt <- tPlt %>%
          dplyr::mutate(STAND_DEAD = STAND_DEAD_MOD) %>%
          dplyr::select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))

      } else {
        tPlt <- tPlt %>%
          dplyr::mutate(STAND_DEAD = AG_OVER_DEAD + BG_OVER_DEAD) %>%
          dplyr::select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))
      }

      ## Convert to long format, where rows are ecosystem components
      tPlt <- tPlt %>%
        tidyr::pivot_longer(cols = -c(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, !!!grpSyms),
                            names_to = 'COMPONENT',
                            values_to = 'cPlot') %>%
        dtplyr::lazy_dt() %>%
        dplyr::mutate(POOL = dplyr::case_when(COMPONENT %in% c('AG_UNDER_LIVE', 'AG_OVER_LIVE') ~ 'AG_LIVE',
                                              COMPONENT %in% c('BG_UNDER_LIVE', 'BG_OVER_LIVE') ~ 'BG_LIVE',
                                              COMPONENT %in% c('STAND_DEAD', 'DOWN_DEAD') ~ 'DEAD_WOOD',
                                              TRUE ~ COMPONENT)) %>%
        dplyr::mutate(cPlot = cPlot * 0.90718474) %>%
        as.data.frame()

      ## Add either component or pool to grpBy, depending on user input
      if (byComponent) grpBy <- c(grpBy, 'POOL', 'COMPONENT')
      if (byPool & !byComponent) grpBy <- c(grpBy, 'POOL')
      grpSyms <- dplyr::syms(grpBy)
      aGrpBy <- grpBy[!c(grpBy %in% c('POOL', 'COMPONENT'))]
      aGrpSyms <- dplyr::syms(aGrpBy)

      ## Summarize across COMPONENTS, if necessary
      if (!byComponent) {
        tPlt <- tPlt %>%
          dtplyr::lazy_dt() %>%
          dplyr::group_by(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, !!!grpSyms) %>%
          dplyr::summarise(cPlot = sum(cPlot, na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          as.data.frame()
      }

      ## Drop all variables except forested area from a
      aPlt <- aPlt %>%
        dplyr::select(ESTN_UNIT_CN, STRATUM_CN, PLT_CN,
                      !!!aGrpSyms, fa = PROP_FOREST)

      ## Adding YEAR to groups
      grpBy <- c('YEAR', grpBy)
      aGrpBy <- c('YEAR', aGrpBy)


      ## Sum variable(s) up to strata then estimation unit level
      eu.sums <- sumToEU(db, tPlt, aPlt, pops, grpBy, aGrpBy, method)
      tEst <- eu.sums$x
      aEst <- eu.sums$y

      out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy)

    }
  }

  return(out)

}


#' @export
carbon <- function(db,
                   grpBy = NULL,
                   polys = NULL,
                   returnSpatial = FALSE,
                   byPool = TRUE,
                   byComponent = FALSE,
                   modelSnag = TRUE,
                   landType = 'forest',
                   method = 'TI',
                   lambda = .5,
                   areaDomain = NULL,
                   totals = FALSE,
                   variance = FALSE,
                   byPlot = FALSE,
                   condList = FALSE,
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
  out <- lapply(X = iter, FUN = carbonStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                byPool, byComponent, modelSnag,
                landType, method,
                lambda, areaDomain,
                totals, byPlot, condList,
                nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
  tEst <- dplyr::bind_rows(out[names(out) == 'tEst'])
  aEst <- dplyr::bind_rows(out[names(out) == 'aEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  aGrpBy <- out[names(out) == 'aGrpBy'][[1]]
  grpSyms <- dplyr::syms(grpBy)
  aGrpSyms <- dplyr::syms(aGrpBy)



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
      dplyr::group_by( !!!aGrpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::select(!!!aGrpSyms, fa_mean, fa_var, nPlots.y)


    tEst <- tEst %>%
      dplyr::group_by(!!!grpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::left_join(aEst, by = aGrpBy) %>%
      dplyr::mutate(CARB_TOTAL = cPlot_mean,
                    AREA_TOTAL = fa_mean,
                    # Ratios
                    CARB_ACRE = CARB_TOTAL / AREA_TOTAL,
                    # Variances
                    CARB_TOTAL_VAR = cPlot_var,
                    AREA_TOTAL_VAR = fa_var,
                    CARB_ACRE_VAR = ratioVar(cPlot_mean, fa_mean, cPlot_var, fa_var, cPlot_cv),
                    # Sampling Errors
                    CARB_TOTAL_SE = sqrt(cPlot_var) / cPlot_mean * 100,
                    AREA_TOTAL_SE = sqrt(fa_var) / fa_mean * 100,
                    CARB_ACRE_SE = sqrt(CARB_ACRE_VAR) / CARB_ACRE * 100,
                    # Plot counts
                    nPlots_AREA = nPlots.y,
                    N = P2PNTCNT_EU) %>%
      dplyr::select(!!!grpSyms, CARB_ACRE, CARB_TOTAL, AREA_TOTAL,
                    CARB_ACRE_VAR, CARB_TOTAL_VAR, AREA_TOTAL_VAR,
                    CARB_ACRE_SE, CARB_TOTAL_SE, AREA_TOTAL_SE,
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




