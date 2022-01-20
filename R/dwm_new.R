dwmStarter <- function(x,
                       db,
                       grpBy_quo = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       landType = 'forest',
                       method = 'TI',
                       lambda = .5,
                       areaDomain = NULL,
                       byPlot = FALSE,
                       condList = FALSE,
                       totals = FALSE,
                       byFuelType = TRUE,
                       nCores = 1,
                       remote,
                       mr){



  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'COND_DWM_CALC', 'COND', 'POP_PLOT_STRATUM_ASSGN',
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
  db$COND_DWM_CALC <- db[['COND_DWM_CALC']] %>% dplyr::mutate(DWM_CN = CN)
  db$COND <- db[['COND']] %>% dplyr::mutate(CND_CN = CN)

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





  ## Handle population tables --------------------------------------------------
  ## Filtering out all inventories that are not relevant to the current estimation
  ## type. If using estimator other than TI, handle the differences in P2POINTCNT
  ## and in assigning YEAR column (YEAR = END_INVYR if method = 'TI')
  pops <- handlePops(db, evalType = c('DWM'), method, mr)

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
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy &
                           !c(names(db$TREE) %in% c(grpP, grpC))]

  ## Dropping irrelevant rows and columns
  db$PLOT <- db$PLOT %>%
    dplyr::select(c(PLT_CN, pltID, STATECD, MACRO_BREAKPOINT_DIA,
                    INVYR, MEASYEAR, PLOT_STATUS_CD,
                    dplyr::all_of(grpP), sp, COUNTYCD)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    dplyr::filter(PLOT_STATUS_CD == 1 & sp == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)
  db$COND <- dplyr::select(db$COND, c(PLT_CN, CONDPROP_UNADJ, PROP_BASIS,
                                      COND_STATUS_CD, CONDID,
                                      dplyr::all_of(grpC), aD, landD)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    dplyr::filter(aD == 1 & landD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)
  db$COND_DWM_CALC <- dplyr::select(db$COND_DWM_CALC, -c( STATECD, COUNTYCD,
                                                          UNITCD, INVYR,
                                                          MEASYEAR, PLOT,
                                                          EVALID)) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)



  ## Full condition list
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$COND_DWM_CALC, by = c('PLT_CN', 'CONDID'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp


  ## Plot-level summaries ------------------------------------------------------
  if (byPlot & !condList){

    grpBy <- c('YEAR', grpBy)
    grpSyms <- syms(grpBy)

    t <- data %>%
      dplyr::mutate(YEAR = INVYR) %>%
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!grpSyms, PLT_CN) %>%
      dplyr::summarize(VOL_1HR = sum(FWD_SM_VOLCF_ADJ * aDI, na.rm = TRUE),
                       VOL_10HR = sum(FWD_MD_VOLCF_ADJ * aDI, na.rm = TRUE),
                       VOL_100HR = sum(FWD_LG_VOLCF_ADJ * aDI, na.rm = TRUE),
                       VOL_1000HR = sum(CWD_VOLCF_ADJ * aDI, na.rm = TRUE),
                       VOL_PILE = sum(PILE_VOLCF_ADJ * aDI, na.rm = TRUE),
                       BIO_DUFF = sum(DUFF_BIOMASS* aDI / 2000, na.rm = TRUE),
                       BIO_LITTER = sum(LITTER_BIOMASS * aDI / 2000, na.rm = TRUE),
                       BIO_1HR = sum(FWD_SM_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                       BIO_10HR = sum(FWD_MD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                       BIO_100HR = sum(FWD_LG_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                       BIO_1000HR = sum(CWD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                       BIO_PILE = sum(PILE_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                       CARB_DUFF = sum(DUFF_CARBON* aDI / 2000, na.rm = TRUE),
                       CARB_LITTER = sum(LITTER_CARBON * aDI / 2000, na.rm = TRUE),
                       CARB_1HR = sum(FWD_SM_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                       CARB_10HR = sum(FWD_MD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                       CARB_100HR = sum(FWD_LG_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                       CARB_1000HR = sum(CWD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                       CARB_PILE = sum(PILE_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                       PROP_FOREST = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      as.data.frame() %>%
      tidyr::pivot_longer(cols = -c(PLT_CN, !!!grpSyms, PROP_FOREST),
                          names_to = c('.value', 'FUEL_TYPE'),
                          names_sep = '_') %>%
      dplyr::rename(VOL_ACRE = VOL,
                    BIO_ACRE = BIO,
                    CARB_ACRE = CARB) %>%
      dplyr::relocate(PROP_FOREST, .after = CARB_ACRE) %>%
      dplyr::mutate(FUEL_TYPE = factor(FUEL_TYPE, levels = c('DUFF', 'LITTER',
                                                             '1HR', '10HR', '100HR',
                                                             '1000HR', 'PILE')))

    ## If by fuel type, add to grpBy
    if (byFuelType) {
      grpBy <- c(grpBy, 'FUEL_TYPE')
      t <- t %>%
        dplyr::arrange(PLT_CN, !!!grpSyms, FUEL_TYPE)
    } else {
      ## Otherwise summarize over fuel types for totals
      t <- t %>%
        dplyr::ungroup() %>%
        dtplyr::lazy_dt() %>%
        dplyr::select(-c(FUEL_TYPE)) %>%
        dplyr::group_by(PLT_CN, !!!grpSyms) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }

    ## Make it spatial
    if (returnSpatial){
      t <- t %>%
        dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
        sf::st_as_sf(coords = c('LON', 'LAT'),
                     crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }

    out <- list(tEst = t, grpBy = grpBy, aGrpBy = NULL)

  } else {

    grpSyms <- dplyr::syms(grpBy)

    ### Condition list for forested area
    a <- data %>%
      dplyr::select(PLT_CN, PROP_BASIS, CONDID, CONDPROP_UNADJ, aDI, !!!grpSyms) %>%
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      dplyr::distinct() %>%
      dplyr::mutate(fa = CONDPROP_UNADJ * aDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!grpSyms, PROP_BASIS, fa)


    ## Return a tree/condition list ready to be handed to `customPSE`
    if (condList) {

      ## All DWM variables have already been adjusted for non-response, so we can
      ## just sum them up here.
      tPlt <- data %>%
        dplyr::distinct(PLT_CN, CONDID, COND_STATUS_CD, .keep_all = TRUE) %>%
        dtplyr::lazy_dt() %>%
        dplyr::mutate(VOL_1HR = FWD_SM_VOLCF_ADJ * aDI,
                         VOL_10HR = FWD_MD_VOLCF_ADJ * aDI,
                         VOL_100HR = FWD_LG_VOLCF_ADJ * aDI,
                         VOL_1000HR = CWD_VOLCF_ADJ * aDI,
                         VOL_PILE = PILE_VOLCF_ADJ * aDI,
                         BIO_DUFF = DUFF_BIOMASS* aDI / 2000,
                         BIO_LITTER = LITTER_BIOMASS * aDI / 2000,
                         BIO_1HR = FWD_SM_DRYBIO_ADJ * aDI / 2000,
                         BIO_10HR = FWD_MD_DRYBIO_ADJ * aDI / 2000,
                         BIO_100HR = FWD_LG_DRYBIO_ADJ * aDI / 2000,
                         BIO_1000HR = CWD_DRYBIO_ADJ * aDI / 2000,
                         BIO_PILE = PILE_DRYBIO_ADJ * aDI / 2000,
                         CARB_DUFF = DUFF_CARBON* aDI / 2000,
                         CARB_LITTER = LITTER_CARBON * aDI / 2000,
                         CARB_1HR = FWD_SM_CARBON_ADJ * aDI / 2000,
                         CARB_10HR = FWD_MD_CARBON_ADJ * aDI / 2000,
                         CARB_100HR = FWD_LG_CARBON_ADJ * aDI / 2000,
                         CARB_1000HR = CWD_CARBON_ADJ * aDI / 2000,
                         CARB_PILE = PILE_CARBON_ADJ * aDI / 2000) %>%
        dplyr::select(PLT_CN, CONDID, !!!grpSyms, VOL_1HR:CARB_PILE) %>%
        as.data.frame() %>%
        tidyr::pivot_longer(cols = -c(PLT_CN, CONDID, !!!grpSyms),
                            names_to = c('.value', 'FUEL_TYPE'),
                            names_sep = '_') %>%
        dplyr::mutate(FUEL_TYPE = factor(FUEL_TYPE, levels = c('DUFF', 'LITTER',
                                                               '1HR', '10HR', '100HR',
                                                               '1000HR', 'PILE'))) %>%
        dplyr::left_join(a, by = c('PLT_CN', 'CONDID', grpBy)) %>%
        dplyr::rename(PROP_FOREST = fa)

      aGrpBy <- grpBy
      ## If by fuel type, add to grpBy
      if (byFuelType) {
        grpBy <- c(grpBy, 'FUEL_TYPE')
        grpSyms <- syms(grpBy)
      } else {
        ## Otherwise summarize over fuel types for totals
        tPlt <- tPlt %>%
          dtplyr::lazy_dt() %>%
          dplyr::select(-c(FUEL_TYPE)) %>%
          dplyr::group_by(PLT_CN, CONDID, !!!grpSyms, AREA_BASIS, PROP_FOREST) %>%
          dplyr::summarise(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          as.data.frame()
      }


      ## Re-order some columns
      tPlt <- tPlt %>%
        dplyr::mutate(EVAL_TYP = 'DWM') %>%
        dplyr::select(PLT_CN, EVAL_TYP, AREA_BASIS,
                      !!!grpSyms, CONDID,
                      VOL_ACRE = VOL,
                      BIO_ACRE = BIO,
                      CARB_ACRE = CARB,
                      PROP_FOREST)

      out <- list(tEst = tPlt, aEst = NULL, grpBy = grpBy, aGrpBy = aGrpBy)


    ## Otherwise, proceed to population estimation
    } else {

      ## Sum variable(s) up to plot-level and adjust for non-response
      aPlt <- sumToPlot(a, pops, grpBy)

      ## All DWM variables have already been adjusted for non-response, so we can
      ## just sum them up here.
      tPlt <- data %>%
        dplyr::distinct(STRATUM_CN, PLT_CN, CONDID, COND_STATUS_CD, .keep_all = TRUE) %>%
        dtplyr::lazy_dt() %>%
        dplyr::group_by(STRATUM_CN, PLT_CN, !!!grpSyms) %>%
        dplyr::summarize(VOL_1HR = sum(FWD_SM_VOLCF_ADJ * aDI, na.rm = TRUE),
                         VOL_10HR = sum(FWD_MD_VOLCF_ADJ * aDI, na.rm = TRUE),
                         VOL_100HR = sum(FWD_LG_VOLCF_ADJ * aDI, na.rm = TRUE),
                         VOL_1000HR = sum(CWD_VOLCF_ADJ * aDI, na.rm = TRUE),
                         VOL_PILE = sum(PILE_VOLCF_ADJ * aDI, na.rm = TRUE),
                         BIO_DUFF = sum(DUFF_BIOMASS* aDI / 2000, na.rm = TRUE),
                         BIO_LITTER = sum(LITTER_BIOMASS * aDI / 2000, na.rm = TRUE),
                         BIO_1HR = sum(FWD_SM_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                         BIO_10HR = sum(FWD_MD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                         BIO_100HR = sum(FWD_LG_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                         BIO_1000HR = sum(CWD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                         BIO_PILE = sum(PILE_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                         CARB_DUFF = sum(DUFF_CARBON* aDI / 2000, na.rm = TRUE),
                         CARB_LITTER = sum(LITTER_CARBON * aDI / 2000, na.rm = TRUE),
                         CARB_1HR = sum(FWD_SM_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                         CARB_10HR = sum(FWD_MD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                         CARB_100HR = sum(FWD_LG_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                         CARB_1000HR = sum(CWD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                         CARB_PILE = sum(PILE_CARBON_ADJ * aDI / 2000, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(dplyr::distinct(dplyr::select(pops, STRATUM_CN, ESTN_UNIT_CN)), by = 'STRATUM_CN') %>%
        as.data.frame() %>%
        tidyr::pivot_longer(cols = -c(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, !!!grpSyms),
                            names_to = c('.value', 'FUEL_TYPE'),
                            names_sep = '_') %>%
        dplyr::rename(VOL = VOL,
                      BIO = BIO,
                      CARB = CARB) %>%
        dplyr::mutate(FUEL_TYPE = factor(FUEL_TYPE, levels = c('DUFF', 'LITTER',
                                                               '1HR', '10HR', '100HR',
                                                               '1000HR', 'PILE')))

      aGrpBy <- grpBy
      ## If by fuel type, add to grpBy
      if (byFuelType) {
        grpBy <- c(grpBy, 'FUEL_TYPE')
      } else {
        ## Otherwise summarize over fuel types for totals
        tPlt <- tPlt %>%
          dtplyr::lazy_dt() %>%
          dplyr::select(-c(FUEL_TYPE)) %>%
          dplyr::group_by(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, !!!grpSyms) %>%
          dplyr::summarise(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          as.data.frame()
      }

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
dwm <- function(db,
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
                condList = FALSE,
                byFuelType = TRUE,
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
  out <- lapply(X = iter, FUN = dwmStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                landType, method,
                lambda, areaDomain,
                byPlot, condList,
                totals, byFuelType,
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
      dplyr::mutate(VOL_TOTAL = VOL_mean,
                    BIO_TOTAL = BIO_mean,
                    CARB_TOTAL = CARB_mean,
                    AREA_TOTAL = fa_mean,
                    # Ratios
                    VOL_ACRE = VOL_TOTAL / AREA_TOTAL,
                    BIO_ACRE = BIO_TOTAL / AREA_TOTAL,
                    CARB_ACRE = CARB_TOTAL / AREA_TOTAL,
                    # Variances
                    VOL_TOTAL_VAR = VOL_var,
                    BIO_TOTAL_VAR = BIO_var,
                    CARB_TOTAL_VAR = CARB_var,
                    AREA_TOTAL_VAR = fa_var,
                    VOL_ACRE_VAR = ratioVar(VOL_mean, fa_mean, VOL_var, fa_var, VOL_cv),
                    BIO_ACRE_VAR = ratioVar(BIO_mean, fa_mean, BIO_var, fa_var, BIO_cv),
                    CARB_ACRE_VAR = ratioVar(CARB_mean, fa_mean, CARB_var, fa_var, CARB_cv),
                    # Sampling Errors
                    VOL_TOTAL_SE = sqrt(VOL_var) / VOL_mean * 100,
                    BIO_TOTAL_SE = sqrt(BIO_var) / BIO_mean * 100,
                    CARB_TOTAL_SE = sqrt(CARB_var) / CARB_mean * 100,
                    AREA_TOTAL_SE = sqrt(fa_var) / fa_mean * 100,
                    VOL_ACRE_SE = sqrt(VOL_ACRE_VAR) / VOL_ACRE * 100,
                    BIO_ACRE_SE = sqrt(BIO_ACRE_VAR) / BIO_ACRE * 100,
                    CARB_ACRE_SE = sqrt(CARB_ACRE_VAR) / CARB_ACRE * 100,
                    # Plot counts
                    nPlots_DWM = nPlots.x,
                    nPlots_AREA = nPlots.y,
                    N = P2PNTCNT_EU) %>%
      dplyr::select(!!!grpSyms, VOL_ACRE, BIO_ACRE, CARB_ACRE,
                    VOL_TOTAL, BIO_TOTAL, CARB_TOTAL, AREA_TOTAL,
                    VOL_ACRE_VAR, BIO_ACRE_VAR, CARB_ACRE_VAR,
                    VOL_TOTAL_VAR, BIO_TOTAL_VAR, CARB_TOTAL_VAR, AREA_TOTAL_VAR,
                    VOL_ACRE_SE, BIO_ACRE_SE, CARB_ACRE_SE,
                    VOL_TOTAL_SE, BIO_TOTAL_SE, CARB_TOTAL_SE, AREA_TOTAL_SE,
                    nPlots_DWM, nPlots_AREA, N)

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

