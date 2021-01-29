acStarter <- function(x,
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
                      chngType = 'net',
                      nCores = 1,
                      remote,
                      mr){

  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'TREE', 'COND', 'SUBP_COND_CHNG_MTRX',
                 'POP_PLOT_STRATUM_ASSGN',
                 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')

  ## If remote, read in state by state. Otherwise, drop all unneccesary tables
  db <- readRemoteHelper(x, db, remote, reqTables, nCores)


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
    mutate(PLT_CN = CN,
           pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

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
    db <- arealSumPrep2(db, grpBy, polys, nCores)
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

  ## Make land status a grouping variable if estimating change components
  if (chngType == 'component') {
    grpBy <- c('STATUS', grpBy)

    if (landType == 'forest') {
      db$COND <- db$COND %>%
        mutate(STATUS = case_when(landD == 1 ~ 'Forest',
                                  TRUE ~ 'Non-forest'))
    } else {
      db$COND <- db$COND %>%
        mutate(STATUS = case_when(landD == 1 ~ 'Timber',
                                  TRUE ~ 'Non-timber'))
    }
  }

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
  pops <- handlePops(db, evalType = c('EXPCHNG'), method, mr)

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
    grpByOrig <- c(grpByOrig, 'landType')
    db$COND <- db$COND %>%
      mutate(landType = case_when(
        COND_STATUS_CD == 1 & SITECLCD %in% c(1:6) & RESERVCD ==0 ~ 'Timber',
        COND_STATUS_CD == 1 ~ 'Non-Timber Forest',
        COND_STATUS_CD == 2 ~ 'Non-Forest',
        COND_STATUS_CD == 3 | COND_STATUS_CD == 4 ~ 'Water'))
    db$COND <- db$COND[!is.na(db$COND$landType),]
  }







  ## Slim down the database before we hand it off to the estimators ------------
  ## Reduces memory requirements and speeds up processing ----------------------

  ## Only the necessary plots for EVAL of interest
  ## We want the plots in the inventory and their previous measurements
  keepThese <- distinct(pops, PLT_CN) %>% left_join(select(db$PLOT, PLT_CN, PREV_PLT_CN), by = 'PLT_CN')
  db$PLOT <- filter(db$PLOT, PLT_CN %in% c(keepThese$PLT_CN, keepThese$PREV_PLT_CN))

  ## Narrow up the tables to the necessary variables
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy &
                           !c(names(db$COND) %in% grpP)]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy &
                           !c(names(db$TREE) %in% c(grpP, grpC))]

  ### Only joining tables necessary to produce plot level estimates
  db$PLOT <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA',
                               'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', COUNTYCD,
                               'PREV_PLT_CN', 'REMPER',
                               all_of(grpP), 'aD_p', 'sp'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS',
                               'CONDID',
                               all_of(grpC), 'aD_c', 'landD')) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  db$TREE <- select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ',
                               'SUBP', 'TREE', all_of(grpT), 'tD')) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  db$SUBP_COND_CHNG_MTRX <- select(db$SUBP_COND_CHNG_MTRX,
                                   c('PLT_CN', 'SUBP', 'CONDID', 'PREVCOND',
                                     'SUBPTYP', 'SUBPTYP_PROP_CHNG')) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))





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
      out <- parLapply(cl, X = names(plts), fun = acHelper1, plts,
                       db[names(db) %in% c('COND', 'TREE', 'SUBP_COND_CHNG_MTRX')],
                       grpBy, byPlot, keepThese, chngType, areaDomain, treeDomain)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = acHelper1, plts,
                      db[names(db) %in% c('COND', 'TREE', 'SUBP_COND_CHNG_MTRX')],
                      grpBy, byPlot, keepThese, chngType, areaDomain, treeDomain,
                      mc.cores = nCores)
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

    }

    out <- list(tEst = tOut, grpBy = grpBy, grpByOrig = grpByOrig)
    ## Population estimation
  } else {
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    t <- bind_rows(out[names(out) == 't'])
    grpBy <- out[names(out) == 'grpBy'][[1]]

    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)

    ## An iterator for population estimation
    popState <- split(pops, as.factor(pops$STATECD))

    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        out <- parLapply(cl, X = names(popState), fun = acHelper2, popState, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = acHelper2, popState, t, grpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
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
      tEst <- tEst %>%
        left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
        left_join(wgts, by = joinCols) %>%
        mutate(across(c(pEst, cEst), ~(.*wgt))) %>%
        mutate(across(c(pVar, cVar), ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize(across(pEst:plotIn_AREA, sum, na.rm = TRUE))




      ## If using an ANNUAL estimator --------------------------------------------
    } else if (str_to_upper(method) == 'ANNUAL') {

      # If INVYR is in YEAR, choose the estimates when INVYR == YEAR
      # Otherwise, choose the estimates produced with the most plots
      tEst <- filterAnnual(tEst, grpBy, plotIn_AREA, db$POP_ESTN_UNIT)

    }


    out <- list(tEst = tEst, grpBy = grpBy, grpByOrig = grpByOrig)
  }

  return(out)

}


#' @export
areaChange <- function (db,
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
                        chngType = 'net',
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
  out <- lapply(X = iter, FUN = acStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                byLandType, landType, method,
                lambda, treeDomain, areaDomain,
                totals, byPlot, chngType, nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  tEst <- bind_rows(out[names(out) == 'tEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
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
    }




    ## Totals and ratios -------------------------------------------------------
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    if (chngType == 'net') {
      suppressWarnings({
        ## Bring them together
        tOut <- tTotal %>%
          #left_join(aTotal, by = grpBy) %>%
          # Renaming, computing ratios, and SE
          mutate(PERC_CHNG = cEst / pEst * 100,
                 AREA_CHNG = cEst,
                 PREV_AREA = pEst,

                 AREA_CHNG_SE = sqrt(cVar) / abs(AREA_CHNG) *100,
                 PREV_AREA_SE = sqrt(pVar) / PREV_AREA *100,

                 ## ratio variance
                 rVar = (1/pEst^2) * (cVar + (PERC_CHNG^2 * pVar) - 2 * PERC_CHNG * cCV),
                 PERC_CHNG_SE = sqrt(rVar) / abs(PERC_CHNG) * 100,
                 PERC_CHNG_VAR = rVar,
                 #N = sum(N),
                 AREA_CHNG_VAR = cVar,
                 PREV_AREA_VAR = pVar,
                 nPlots_AREA = plotIn_AREA) %>%
          select(grpBy, PERC_CHNG, AREA_CHNG, PREV_AREA,
                 PERC_CHNG_SE, AREA_CHNG_SE, PREV_AREA_SE,
                 PERC_CHNG_VAR, AREA_CHNG_VAR, PREV_AREA_VAR,
                 nPlots_AREA, N)
      })



      if (variance) {
        tOut <- tOut %>%
          select(-c(PERC_CHNG_SE, AREA_CHNG_SE, PREV_AREA_SE))
      } else {
        tOut <- tOut %>%
          select(-c(PERC_CHNG_VAR, AREA_CHNG_VAR, PREV_AREA_VAR, N))
      }

      # Snag the names
      tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

      ## Pretty output
      tOut <- tOut %>%
        ungroup() %>%
        mutate_if(is.factor, as.character) %>%
        drop_na(grpBy) %>%
        arrange(YEAR) %>%
        as_tibble()


    ## Change by component
    } else {
      suppressWarnings({
        ## Bring them together
        tOut <- tTotal %>%
          filter(cEst > 0) %>%
          mutate(AREA_CHNG = cEst,
                 AREA_CHNG_SE = sqrt(cVar) / abs(AREA_CHNG) *100,
                 AREA_CHNG_VAR = cVar,
                 nPlots_AREA = plotIn_AREA) %>%
          select(grpBy, AREA_CHNG,
                 AREA_CHNG_SE,
                 AREA_CHNG_VAR,
                 nPlots_AREA, N)
      })

      # Snag the names
      tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

      if (variance) {
        tOut <- tOut %>%
          select(-c(AREA_CHNG_SE))
      } else {
        tOut <- tOut %>%
          select(-c(AREA_CHNG_VAR, N))
      }
    }

    ## Pretty output
    tOut <- tOut %>%
      ungroup() %>%
      mutate_if(is.factor, as.character) %>%
      arrange(YEAR) %>%
      as_tibble()


  }





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
