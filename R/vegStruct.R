vegStructStarter <- function(x,
                            db,
                            grpBy_quo = NULL,
                            polys = NULL,
                            returnSpatial = FALSE,
                            landType = 'forest',
                            method = 'TI',
                            lambda = .5,
                            areaDomain = NULL,
                            byPlot = FALSE,
                            totals = FALSE,
                            nCores = 1,
                            remote = NULL,
                            mr){


  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'P2VEG_SUBP_STRUCTURE', 'SUBP_COND', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')

  ## If remote, read in state by state. Otherwise, drop all unnecessary tables
  db <- readRemoteHelper(x, db, remote, reqTables, nCores)

  ## IF the object was clipped
  if ('prev' %in% names(db$PLOT)){
    ## Only want the current plots, no grm
    db$PLOT <- filter(db$PLOT, prev == 0)
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
  if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method,
                  'unknown. Defaulting to Temporally Indifferent (TI).'))
  }


  ## Prep other variables ------------------------------------------------------
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>%
    mutate(PLT_CN = CN,
           pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ## Reduce our sample right off the bat, only plots sampled for veg
  db$PLOT <- filter(db$PLOT, P2VEG_SAMPLING_STATUS_CD %in% 1:3)

  ## Convert grpBy to character
  grpBy <- grpByToChar(db[!c(names(db) %in% c('TREE'))], grpBy_quo)

  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy


  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}

  ## Add variables to grouping
  grpBy <- c(grpBy, 'LAYER', 'GROWTH_HABIT')

  ## ADDING names of id columns for layer and growth habit
  db$P2VEG_SUBP_STRUCTURE <- db$P2VEG_SUBP_STRUCTURE %>%
    mutate(LAYER = case_when(is.na(LAYER) ~ NA_character_,
                             LAYER == 1 ~ '0 to 2.0 feet',
                             LAYER == 2 ~ '2.1 to 6.0 feet',
                             LAYER == 3 ~ '6.1 to 16.0 feet',
                             LAYER == 4 ~ 'Greater than 16 feet',
                             LAYER == 5 ~ 'Areal: all layers'),
           GROWTH_HABIT = case_when(is.na(GROWTH_HABIT_CD) ~ NA_character_,
                                    GROWTH_HABIT_CD == 'TT' ~ 'Tally tree',
                                    GROWTH_HABIT_CD == 'NT' ~ 'Non-tally tree',
                                    GROWTH_HABIT_CD == 'SH' ~ 'Shrubs/ vines',
                                    GROWTH_HABIT_CD == 'FB' ~ 'Forbs',
                                    GROWTH_HABIT_CD == 'GR' ~ 'Graminoids'))


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
  pops <- handlePops(db, evalType = c('EXPCURR'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (str_to_upper(method) != 'TI') {
    pops <- mergeSmallStrata(db, pops)
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

  ### Only joining tables necessary to produce plot level estimates
  db$PLOT <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA',
                               'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD',
                               all_of(grpP), 'aD_p', 'sp', 'COUNTYCD'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS',
                               'COND_STATUS_CD', 'CONDID',
                               all_of(grpC), 'aD_c', 'landD')) %>%
    filter(PLT_CN %in% db$PLOT$PLT_CN)
  db$SUBP_COND <- select(db$SUBP_COND, c(PLT_CN, SUBP, CONDID, SUBPCOND_PROP))%>%
    filter(PLT_CN %in% db$PLOT$PLT_CN)
  db$P2VEG_SUBP_STRUCTURE <- select(db$P2VEG_SUBP_STRUCTURE,
                                    c('PLT_CN', 'COVER_PCT', 'SUBP', 'CONDID',
                                      'LAYER', 'GROWTH_HABIT')) %>%
    filter(PLT_CN %in% db$PLOT$PLT_CN)

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
      out <- parLapply(cl, X = names(plts), fun = vegStructHelper1, plts,
                       db[names(db) %in% c('COND', 'SUBP_COND', 'P2VEG_SUBP_STRUCTURE')],
                       grpBy, aGrpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = vegStructHelper1, plts,
                      db[names(db) %in% c('COND', 'SUBP_COND', 'P2VEG_SUBP_STRUCTURE')],
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

    ## An iterator for population estimation
    popState <- split(pops, as.factor(pops$STATECD))

    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        out <- parLapply(cl, X = names(popState), fun = vegStructHelper2, popState, a, t, grpBy, aGrpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = vegStructHelper2, popState, a, t, grpBy, aGrpBy, method, mc.cores = nCores)
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
        mutate(across(cEst, ~(.*wgt))) %>%
        mutate(across(cVar:cvEst_c, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, A, .dots = grpBy) %>%
        summarize(across(cEst:cvEst_c, sum, na.rm = TRUE))



      ## If using an ANNUAL estimator --------------------------------------------
    } else if (str_to_upper(method) == 'ANNUAL') {

      # If INVYR is in YEAR, choose the estimates when INVYR == YEAR
      # Otherwise, choose the estimates produced with the most plots
      aEst <- filterAnnual(aEst, aGrpBy, plotIn_AREA, db$POP_ESTN_UNIT)
      tEst <- filterAnnual(tEst, grpBy, plotIn_VEG, db$POP_ESTN_UNIT)
    }

    out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy, grpByOrig = grpByOrig)
  }

  return(out)

}





#' @export
vegStruct <- function(db,
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
  out <- lapply(X = iter, FUN = vegStructStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                landType, method,
                lambda, areaDomain,
                byPlot, totals, nCores, remote, mr)
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
    suppressWarnings({
      aTotal <- aEst %>%
        group_by(.dots = aGrpBy) %>%
        summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                  aVar = sum(aVar, na.rm = TRUE),
                  AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
                  AREA_TOTAL_VAR = aVar,
                  nPlots_AREA = sum(plotIn_AREA, na.rm = TRUE))
      tOut <- tEst %>%
        group_by(.dots = grpBy) %>%
        #left_join(aTotal, by = c(aGrpBy)) %>%
        summarize(TOTAL_COVER_AREA = sum(cEst, na.rm = TRUE),
                  cVar = sum(cVar, na.rm = TRUE),
                  cvEst_c = sum(cvEst_c, na.rm = TRUE),
                  ## Sampling Errors
                  TOTAL_COVER_AREA_SE = sqrt(cVar) / TOTAL_COVER_AREA * 100,
                  TOTAL_COVER_AREA_VAR = cVar,
                  N = sum(N, na.rm = TRUE),
                  nPlots_VEG = sum(plotIn_VEG, na.rm = TRUE)) %>%
        left_join(aTotal, by = aGrpBy) %>%
        mutate(COVER_PCT = TOTAL_COVER_AREA / AREA_TOTAL * 100,
               cpVar = (1/AREA_TOTAL^2) * (cVar + (COVER_PCT^2 * aVar) - 2 * COVER_PCT * cvEst_c),
               COVER_PCT_SE = sqrt(cpVar) / COVER_PCT * 100,
               COVER_PCT_VAR = cpVar)

    })

    if (totals) {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, "COVER_PCT","TOTAL_COVER_AREA", "AREA_TOTAL",
                 "COVER_PCT_VAR","TOTAL_COVER_AREA_VAR", "AREA_TOTAL_VAR",
                 "nPlots_VEG", "nPlots_AREA", N)
      } else {
        tOut <- tOut %>%
          select(grpBy, "COVER_PCT","TOTAL_COVER_AREA", "AREA_TOTAL",
                 "COVER_PCT_SE","TOTAL_COVER_AREA_SE", "AREA_TOTAL_SE",
                 "nPlots_VEG", "nPlots_AREA", N)
      }

    } else {
      if (variance){
        tOut <- tOut %>%
          select(grpBy,"COVER_PCT","COVER_PCT_VAR","nPlots_VEG", "nPlots_AREA", N)
      } else {
        tOut <- tOut %>%
          select(grpBy,"COVER_PCT","COVER_PCT_SE","nPlots_VEG", "nPlots_AREA", N)
      }

    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

  } # End byPlot

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


