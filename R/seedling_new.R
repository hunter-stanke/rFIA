seedStarter <- function(x,
                        db,
                        grpBy_quo = NULL,
                        polys = NULL,
                        returnSpatial = FALSE,
                        bySpecies = FALSE,
                        landType = 'forest',
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
  reqTables <- c('PLOT', 'SEEDLING', 'COND', 'POP_PLOT_STRATUM_ASSGN',
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
      first(class(polys)) %in%
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
  db <- udSeedDomain(db, treeDomain)




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
  ## Add species to groups
  if (bySpecies) {
    db$SEEDLING <- db$SEEDLING %>%
      dplyr::left_join(dplyr::select(intData$REF_SPECIES_2018,
                                     c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      dplyr::mutate(SCIENTIFIC_NAME = stringr::str_c(GENUS, SPECIES, sep = ' ')) %>%
      dplyr::mutate_if(is.factor,
                       as.character)
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }




  ## Prep the tree list --------------------------------------------------------
  ## Narrow up the tables to the necessary variables
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy &
                           !c(names(db$COND) %in% grpP)]
  grpT <- names(db$SEEDLING)[names(db$SEEDLING) %in% grpBy &
                           !c(names(db$SEEDLING) %in% c(grpP, grpC))]

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

  db$SEEDLING <- db$SEEDLING %>%
    dplyr::select(c(PLT_CN, CONDID, SPCD, TPA_UNADJ,
                    SUBP, TREECOUNT_CALC, dplyr::all_of(grpT), tD)) %>%
    ## Drop plots outside our domain of interest
    dplyr::filter(TPA_UNADJ > 0 & tD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% db$PLOT$PLT_CN)


  # Separate area grouping names from tree grouping names
  if (!is.null(polys)){
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND), names(polys))]
  } else {
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND))]
  }


  ## Full tree list
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$SEEDLING, by = c('PLT_CN', 'CONDID'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp
  data$tDI <- data$landD * data$aD * data$tD * data$sp


  ## Plot-level summaries ------------------------------------------------------
  if (byPlot & !treeList){

    grpBy <- c('YEAR', grpBy)
    grpSyms <- syms(grpBy)
    aGrpSyms <- syms(aGrpBy)

    ### Plot-level estimates
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(PLT_CN, !!!aGrpSyms) %>%
      dplyr::summarize(PROP_FOREST = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      as.data.frame()

    t <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      dplyr::distinct(PLT_CN, SUBP, SPCD, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!grpSyms, PLT_CN) %>%
      dplyr::summarize(TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
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
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dplyr::mutate(fa = CONDPROP_UNADJ * aDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!aGrpSyms, PROP_BASIS, fa)


    grpSyms <- syms(grpBy)
    ## Tree list
    t <- data %>%
      dplyr::distinct(PLT_CN, SUBP, SPCD, .keep_all = TRUE) %>%
      dplyr::mutate(tPlot = TPA_UNADJ * tDI,
                    TREE_BASIS = 'MICR') %>%
      dplyr::select(PLT_CN, TREE_BASIS, SPCD, !!!grpSyms, tPlot)


    ## Return a tree/condition list ready to be handed to `customPSE`
    if (treeList) {
      tEst <- a %>%
        dplyr::left_join(t, by = c('PLT_CN', aGrpBy)) %>%
        ## customPSE requires these be present, so simply set to NA for seedling
        dplyr::mutate(SUBP = NA,
                      TREE = NA,
                      PROP_FOREST = fa,
                      EVAL_TYP = 'VOL') %>%
        dplyr::group_by(PLT_CN, EVAL_TYP, TREE_BASIS, AREA_BASIS,
                        !!!grpSyms, CONDID, SUBP, TREE,
                        PROP_FOREST) %>%
        dplyr::summarise(TPA = sum(tPlot, na.rm = TRUE)) %>%
        dplyr::ungroup()
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

      out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy)
    }
  }


  return(out)

}

#' @export
seedling <- function(db,
                     grpBy = NULL,
                     polys = NULL,
                     returnSpatial = FALSE,
                     bySpecies = FALSE,
                     landType = 'forest',
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
  out <- lapply(X = iter, FUN = seedStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                bySpecies,
                landType, method,
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
                    AREA_TOTAL = fa_mean,
                    # Ratios
                    TPA = TREE_TOTAL / AREA_TOTAL,
                    # Variances
                    TREE_TOTAL_VAR = tPlot_var,
                    AREA_TOTAL_VAR = fa_var,
                    TPA_VAR = ratioVar(tPlot_mean, fa_mean, tPlot_var, fa_var, tPlot_cv),
                    # Sampling Errors
                    TREE_TOTAL_SE = sqrt(tPlot_var) / tPlot_mean * 100,
                    AREA_TOTAL_SE = sqrt(fa_var) / fa_mean * 100,
                    TPA_SE = sqrt(TPA_VAR) / TPA * 100,
                    # Plot counts
                    nPlots_TREE = nPlots.x,
                    nPlots_AREA = nPlots.y,
                    N = P2PNTCNT_EU) %>%
      dplyr::select(!!!grpSyms, TPA, TREE_TOTAL, AREA_TOTAL,
                    TPA_VAR, TREE_TOTAL_VAR, AREA_TOTAL_VAR,
                    TPA_SE, TREE_TOTAL_SE, AREA_TOTAL_SE,
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




