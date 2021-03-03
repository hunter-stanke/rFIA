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
                            nCores = 1,
                            remote,
                            mr){



  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'COND', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT',
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
      first(class(polys)) %in%
      c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (treeType %in% c('live', 'dead', 'gs', 'all') == FALSE){
    stop('treeType must be one of: "live", "dead", "gs", or "all".')
  }
  if (any(reqTables[!c(reqTables %in% 'SUBP_COND_CHNG_MTRX')] %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '),
               'not found in object db.'))
  }
  if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method,
                  'unknown. Defaulting to Temporally Indifferent (TI).'))
  }
  ## No EXP_GROW available for Western States, make sure we warn that values will be returned as 0
  # These states do not allow temporal queries. Things are extremely weird with their eval groups
  noGrow <- c(02,03,04,07,08,11,14,15,16, 30, 32, 35,43,49, 78)
  if(any(unique(db$PLOT$STATECD) %in% noGrow)){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    warning(paste('Recruitment data unavailable for: ', toString(fancyName) , '. Returning 0 for all recruitment estimates which include these states.', sep = ''))
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
    mutate(PLT_CN = CN,
           pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))
  db$TREE <- db$TREE %>%
    mutate(TRE_CN = CN)

  ## Convert grpBy to character
  grpBy <- grpByToChar(db, grpBy_quo)

  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy


  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}


  ## Intersect plots with polygons if polygons are given
  if (!is.null(polys)){

    ## Add shapefile names to grpBy
    grpBy = c(grpBy, 'polyID')
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
  if (str_to_upper(stateVar) == 'TPA'){
    db$TREE_GRM_MIDPT$state <- 1
    db$TREE$state_recr <- 1
  } else if (str_to_upper(stateVar) == 'BAA'){
    db$TREE_GRM_MIDPT$state <- basalArea(db$TREE_GRM_MIDPT$DIA)
    db$TREE$state_recr <- basalArea(db$TREE$DIA)
  } else if (str_to_upper(stateVar) == 'SAWVOL'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLCSNET
    db$TREE$state_recr <- db$TREE$VOLCSNET
  } else if (str_to_upper(stateVar) == 'NETVOL'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLCFNET
    db$TREE$state_recr <- db$TREE$VOLCFNET
  } else if (str_to_upper(stateVar) == 'SNDVOL'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLCFSND
    db$TREE$state_recr <- db$TREE$VOLCFSND
  } else if (str_to_upper(stateVar) == 'BIO_AG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_AG
    db$TREE$state_recr <- db$TREE$DRYBIO_AG
  } else if (str_to_upper(stateVar) == 'BIO_BG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_BG
    db$TREE$state_recr <- db$TREE$DRYBIO_BG
  } else if (str_to_upper(stateVar) == 'BIO'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_BG + db$TREE_GRM_MIDPT$DRYBIO_AG
    db$TREE$state_recr <- db$TREE$DRYBIO_BG + db$TREE$DRYBIO_AG
  } else if (str_to_upper(stateVar) == 'CARB_AG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_AG * .5
    db$TREE$state_recr <- db$TREE$DRYBIO_AG * .5
  } else if (str_to_upper(stateVar) == 'CARB_BG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_BG * .5
    db$TREE$state_recr <- db$TREE$DRYBIO_BG * .5
  } else if (str_to_upper(stateVar) == 'CARB'){
    db$TREE_GRM_MIDPT$state <- (db$TREE_GRM_MIDPT$DRYBIO_AG + db$TREE_GRM_MIDPT$DRYBIO_BG) * .5
    db$TREE$state_recr <- (db$TREE$DRYBIO_AG + db$TREE$DRYBIO_BG) * .5
  } else {
    stop(paste0('Method not known for stateVar: ', stateVar, '. Please choose one of: TPA, BAA, SAWVOL, NETVOL, BIO_AG, BIO_BG, BIO, CARB_AG, CARB_BG, or CARB.' ))
  }






  ## Build a domain indicator for each observation (1 or 0) --------------------
  ## Land type and tree type combined
  db <- typeDomain_grow(db, treeType, landType, type = 'gm')

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
  pops <- handlePops(db, evalType = c('EXPGROW', 'EXPMORT', 'EXPREMV'), method, mr, ga = TRUE)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (str_to_upper(method) != 'TI') {
    pops <- mergeSmallStrata(db, pops)
  }




  ## Canned groups -------------------------------------------------------------
  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      left_join(select(intData$REF_SPECIES_2018,
                       c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' ')) %>%
      mutate_if(is.factor,
                as.character)
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }

  ## Break into size classes
  if (bySizeClass){
    grpBy <- c(grpBy, 'sizeClass')
    grpByOrig <- c(grpByOrig, 'sizeClass')
    db$TREE$sizeClass <- makeClasses(db$TREE$DIA, interval = 2, numLabs = TRUE)
    db$TREE <- db$TREE[!is.na(db$TREE$sizeClass),]
  }





  ## Slim down the database for we hand it off to the estimators ---------------
  ## Reduces memory requirements and speeds up processing ----------------------

  ## Only the necessary plots for EVAL of interest
  db$PLOT <- filter(db$PLOT, PLT_CN %in% pops$PLT_CN)

  ## Narrow up the tables to the necessary variables
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy &
                           !c(names(db$COND) %in% grpP)]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy &
                           !c(names(db$TREE) %in% c(grpP, grpC))]

  ### Only joining tables necessary to produce plot level estimates
  db$PLOT <- select(db$PLOT, c('PLT_CN', 'STATECD', 'COUNTYCD',
                               'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR',
                               'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER',
                               all_of(grpP), 'aD_p', 'sp'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS',
                               'COND_STATUS_CD', 'CONDID', all_of(grpC),
                               'aD_c', 'landD')) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  db$TREE <- select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN',
                               'PREV_TRE_CN', 'SUBP', 'TREE', all_of(grpT), 'tD',
                               'typeD', 'state_recr', TPA_UNADJ,
                               STATUSCD, DIA)) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  db$TREE_GRM_COMPONENT <- db$TREE_GRM_COMPONENT %>%
    select(c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPARECR_UNADJ',
             'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')) %>%
    filter(TRE_CN %in% db$TREE$TRE_CN)
  db$TREE_GRM_MIDPT <- db$TREE_GRM_MIDPT %>%
    select(c('TRE_CN', 'DIA', 'state')) %>%
    filter(TRE_CN %in% db$TREE$TRE_CN)

  if ('SUBP_COND_CHNG_MTRX' %in% names(db)) {
    db$SUBP_COND_CHNG_MTRX <- select(db$SUBP_COND_CHNG_MTRX, PLT_CN, PREV_PLT_CN,
                                     SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID) %>%
      filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  }



  # Separate area grouping names from tree grouping names
  if (!is.null(polys)){
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND), names(polys))]
  } else {
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND))]
  }


  ## Compute plot-level summaries ----------------------------------------------
  ## An iterator for plot-level summaries
  plts <- split(db$PLOT, as.factor(paste(db$PLOT$COUNTYCD, db$PLOT$STATECD, sep = '_')))
  suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      clusterEvalQ(cl, {
        library(dplyr)
        library(stringr)
        library(rFIA)
      })
      out <- parLapply(cl, X = names(plts), fun = gmHelper1, plts,
                       db[names(db) %in% c('COND', 'TREE', 'TREE_GRM_COMPONENT',
                                           'TREE_GRM_MIDPT', 'SUBP_COND_CHNG_MTRX')],
                       grpBy, aGrpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = gmHelper1, plts,
                      db[names(db) %in% c('COND', 'TREE', 'TREE_GRM_COMPONENT',
                                          'TREE_GRM_MIDPT', 'SUBP_COND_CHNG_MTRX')],
                      grpBy, aGrpBy, byPlot, mc.cores = nCores)
    }
  })






  ## If byPlot, return plot-level estimates ------------------------------------
  ## Otherwise continue to population estimation -------------------------------
  if (byPlot){
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tOut <- bind_rows(out[names(out) == 't'])
    ## Make it spatial
    if (returnSpatial){
      tOut <- tOut %>%
        filter(!is.na(LAT) & !is.na(LON)) %>%
        st_as_sf(coords = c('LON', 'LAT'),
                 crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

      ## Modify some names if a different state variable was given
      names(tOut) <- str_replace(names(tOut), 'TPA', paste(stateVar, 'ACRE', sep = '_'))

    }
    out <- list(tEst = tOut, grpBy = grpBy, aGrpBy = aGrpBy, grpByOrig = grpByOrig)



  ## Population estimation
  } else {
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    a <- bind_rows(out[names(out) == 'a'])
    t <- bind_rows(out[names(out) == 't'])


    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)
    aGrpBy <- c('YEAR', aGrpBy)

    ## Splitting up by ESTN_UNIT
    popState <- split(pops, as.factor(pops$STATECD))

    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        out <- parLapply(cl, X = names(popState), fun = gmHelper2, popState, a, t, grpBy, aGrpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = gmHelper2, popState, a, t, grpBy, aGrpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    aEst <- bind_rows(out[names(out) == 'aEst'])
    tEst <- bind_rows(out[names(out) == 'tEst'])





    ## Compute moving average weights if not TI ----------------------------------
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){

      ## Compute the weights
      wgts <- maWeights(pops, method, lambda)

      ## If moving average ribbons, add lambda to grpBy for easier summary
      if (str_to_upper(method) == 'EMA' & length(lambda) > 1){
        grpBy <- c('lambda', grpBy)
        aGrpBy <- c('lambda', aGrpBy)
      }


      ## Apply the weights
      if (str_to_upper(method) %in% c('LMA', 'EMA')){
        joinCols <- c('YEAR', 'STATECD', 'INVYR')
      } else {
        joinCols <- c('YEAR', 'STATECD')
      }
      aEst <- aEst %>%
        left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
        left_join(wgts, by = joinCols) %>%
        mutate(across(aEst, ~(.*wgt))) %>%
        mutate(across(aVar, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
        summarize(across(aEst:plotIn_AREA, sum, na.rm = TRUE))
      tEst <- tEst %>%
        left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
        left_join(wgts, by = joinCols) %>%
        mutate(across(tEst:hEst, ~(.*wgt))) %>%
        mutate(across(tVar:cvEst_hp, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize(across(tEst:plotIn_h, sum, na.rm = TRUE))



      ## If using an ANNUAL estimator --------------------------------------------
    } else if (str_to_upper(method) == 'ANNUAL') {

      # If INVYR is in YEAR, choose the estimates when INVYR == YEAR
      # Otherwise, choose the estimates produced with the most plots
      aEst <- filterAnnual(aEst, aGrpBy, plotIn_AREA, db$POP_ESTN_UNIT)
      tEst <- filterAnnual(tEst, grpBy, plotIn_t, db$POP_ESTN_UNIT)
    }

    out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy, grpByOrig = grpByOrig)

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
                totals, byPlot, nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
  aEst <- bind_rows(out[names(out) == 'aEst'])
  tEst <- bind_rows(out[names(out) == 'tEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  aGrpBy <- out[names(out) == 'aGrpBy'][[1]]
  grpByOrig <- out[names(out) == 'grpByOrig'][[1]]



  ## Plot-level estimates
  if (byPlot){

    ## Name change for consistency below - awkward, I know. I don't care.
    tOut <- tEst

    ## Population estimates
  } else {

    ## Combine most-recent population estimates across states with potentially
    ## different reporting schedules, i.e., if 2016 is most recent in MI and 2017 is
    ## most recent in WI, combine them and label as 2017
    if (mr) {
      tEst <- combineMR(tEst, grpBy)
      aEst <- combineMR(aEst, aGrpBy)
    }



    ## Totals and ratios -------------------------------------------------------
    aTotal <- aEst %>%
      group_by(.dots = aGrpBy) %>%
      summarize_all(sum,na.rm = TRUE)
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    suppressWarnings({
      ## Bring them together
      tOut <- tTotal %>%
        left_join(aTotal, by = aGrpBy) %>%
        # Renaming, computing ratios, and SE
        mutate(TREE_TOTAL = tEst,
               RECR_TREE_TOTAL = rEst,
               MORT_TREE_TOTAL = mEst,
               REMV_TREE_TOTAL = hEst,
               AREA_TOTAL = aEst,
               RECR_TPA = RECR_TREE_TOTAL / AREA_TOTAL,
               MORT_TPA = MORT_TREE_TOTAL / AREA_TOTAL,
               REMV_TPA = REMV_TREE_TOTAL / AREA_TOTAL,
               RECR_PERC = RECR_TREE_TOTAL / TREE_TOTAL * 100,
               MORT_PERC = MORT_TREE_TOTAL / TREE_TOTAL * 100,
               REMV_PERC = REMV_TREE_TOTAL / TREE_TOTAL * 100,
               ## Ratio Var
               raVar = (1/AREA_TOTAL^2) * (rVar + (RECR_TPA^2 * aVar) - 2 * RECR_TPA * cvEst_r),
               maVar = (1/AREA_TOTAL^2) * (mVar + (MORT_TPA^2 * aVar) - 2 * MORT_TPA * cvEst_m),
               haVar = (1/AREA_TOTAL^2) * (hVar + (REMV_TPA^2 * aVar) - 2 * REMV_TPA * cvEst_h),
               rpVar = (1/TREE_TOTAL^2) * (rVar + (RECR_PERC^2 * tVar) - 2 * RECR_PERC * cvEst_rp),
               mpVar = (1/TREE_TOTAL^2) * (mVar + (MORT_PERC^2 * tVar) - 2 * MORT_PERC * cvEst_mp),
               hpVar = (1/TREE_TOTAL^2) * (hVar + (REMV_PERC^2 * tVar) - 2 * REMV_PERC * cvEst_hp),
               ## SE RATIO
               RECR_TPA_SE = sqrt(raVar) / RECR_TPA * 100,
               MORT_TPA_SE = sqrt(maVar) / MORT_TPA * 100,
               REMV_TPA_SE = sqrt(haVar) / REMV_TPA * 100,
               RECR_PERC_SE = sqrt(rpVar) / RECR_PERC * 100,
               MORT_PERC_SE = sqrt(mpVar) / MORT_PERC * 100,
               REMV_PERC_SE = sqrt(hpVar) / REMV_PERC * 100,
               ## Var ratio
               RECR_TPA_VAR = raVar,
               MORT_TPA_VAR = maVar,
               REMV_TPA_VAR = haVar,
               RECR_PERC_VAR = rpVar,
               MORT_PERC_VAR = mpVar,
               REMV_PERC_VAR = hpVar,
               ## SE TOTAL
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL *100,
               TREE_TOTAL_SE = sqrt(tVar) / TREE_TOTAL *100,
               RECR_TREE_TOTAL_SE = sqrt(rVar) / RECR_TREE_TOTAL *100,
               MORT_TREE_TOTAL_SE = sqrt(mVar) / MORT_TREE_TOTAL *100,
               REMV_TREE_TOTAL_SE = sqrt(hVar) / REMV_TREE_TOTAL *100,
               ## VAR TOTAL
               AREA_TOTAL_VAR = aVar,
               TREE_TOTAL_VAR = tVar,
               RECR_TREE_TOTAL_VAR = rVar,
               MORT_TREE_TOTAL_VAR = mVar,
               REMV_TREE_TOTAL_VAR = hVar,
               ## nPlots
               # Non-zero plots
               nPlots_TREE = plotIn_t,
               nPlots_RECR = plotIn_r,
               nPlots_MORT = plotIn_m,
               nPlots_REMV = plotIn_h,
               nPlots_AREA = plotIn_AREA)
    })

    # Make some columns go away
    if (totals) {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
                 TREE_TOTAL, RECR_TREE_TOTAL, MORT_TREE_TOTAL, REMV_TREE_TOTAL, AREA_TOTAL,
                 RECR_TPA_VAR, MORT_TPA_VAR, REMV_TPA_VAR, RECR_PERC_VAR, MORT_PERC_VAR, REMV_PERC_VAR,
                 TREE_TOTAL_VAR, RECR_TREE_TOTAL_VAR, MORT_TREE_TOTAL_VAR, REMV_TREE_TOTAL_VAR, AREA_TOTAL_VAR,
                 nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV, nPlots_AREA, N)
      } else {
        tOut <- tOut %>%
          select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
                 TREE_TOTAL, RECR_TREE_TOTAL, MORT_TREE_TOTAL, REMV_TREE_TOTAL, AREA_TOTAL,
                 RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE, RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE,
                 TREE_TOTAL_SE, RECR_TREE_TOTAL_SE, MORT_TREE_TOTAL_SE, REMV_TREE_TOTAL_SE, AREA_TOTAL_SE,
                 nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV, nPlots_AREA)
      }

    } else {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
                 RECR_TPA_VAR, MORT_TPA_VAR, REMV_TPA_VAR, RECR_PERC_VAR, MORT_PERC_VAR, REMV_PERC_VAR,
                 nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV,nPlots_AREA, N)
      } else {
        tOut <- tOut %>%
          select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
                 RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE, RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE,
                 nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV,nPlots_AREA)
      }

    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]
    ## Modify some names if a different state variable was given
    #names(tOut) <- str_replace(names(tOut), 'TPA', paste(stateVar, 'ACRE', sep = '_'))

  }

  ## Modify some names if a different state variable was given
  if (stateVar != 'TPA') {
    names(tOut) <- str_replace(names(tOut), 'TPA', paste(stateVar, 'ACRE', sep = '_'))
    names(tOut) <- str_replace(names(tOut), 'TREE', ifelse(stateVar == 'BAA', 'BA', stateVar))
  }
  names(tOut) <- str_replace(names(tOut), 'BAA_ACRE', 'BAA')

  # Snag the names
  tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

  ## Pretty output
  tOut <- tOut %>%
    ungroup() %>%
    mutate_if(is.factor, as.character) %>%
    drop_na(grpBy) %>%
    arrange(YEAR) %>%
    as_tibble()


  ## Make implicit NA explicit for spatial summaries
  ## Not sure if I like this or not, but I'm going with it for now
  tOut <- prettyNamesSF(tOut, polys, byPlot, grpBy, grpByOrig, tNames, returnSpatial)


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)

  ## remove any duplicates in byPlot
  ## Also make PLOT_STATUS_CD more informative
  if (byPlot) {
    tOut <- unique(tOut)
    tOut <- tOut %>%
      mutate(PLOT_STATUS = case_when(is.na(PLOT_STATUS_CD) ~ NA_character_,
                                     PLOT_STATUS_CD == 1 ~ 'Forest',
                                     PLOT_STATUS_CD == 2 ~ 'Non-forest',
                                     PLOT_STATUS_CD == 3 ~ 'Non-sampled')) %>%
      relocate(PLOT_STATUS, .after = PLOT_STATUS_CD)
  }

  return(tOut)
}
