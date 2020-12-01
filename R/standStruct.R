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
                               tidy = TRUE,
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

  ## Convert grpBy to character
  grpBy <- grpByToChar(db[!c(names(db) %in% 'TREE')], grpBy_quo)

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
  if (method != 'TI') {
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
  db$TREE <- select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'TPA_UNADJ',
                               'CCLCD', 'SUBP', 'TREE')) %>%
    filter(PLT_CN %in% db$PLOT$PLT_CN)




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
      out <- parLapply(cl, X = names(plts), fun = ssHelper1, plts,
                       db[names(db) %in% c('COND', 'TREE')],
                       grpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = ssHelper1, plts,
                      db[names(db) %in% c('COND', 'TREE')],
                      grpBy, byPlot, mc.cores = nCores)
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


    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)

    ## An iterator for population estimation
    popState <- split(pops, as.factor(pops$STATECD))
    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        out <- parLapply(cl, X = names(popState), fun = ssHelper2, popState, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = ssHelper2, popState, t, grpBy, method, mc.cores = nCores)
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
        mutate(across(aEst:moEst, ~(.*wgt))) %>%
        mutate(across(aVar:cvEst_mo, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize(across(aEst:plotIn_AREA, sum, na.rm = TRUE))






      ## If using an ANNUAL estimator --------------------------------------------
    } else if (str_to_upper(method) == 'ANNUAL') {

      # If INVYR is in YEAR, choose the estimates when INVYR == YEAR
      # Otherwise, choose the estimates produced with the most plots
      tEst <- filterAnnual(tEst, grpBy, plotIn_AREA)
    }

    out <- list(tEst = tEst, grpBy = grpBy, grpByOrig = grpByOrig)
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
                          tidy = TRUE,
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
                totals, byPlot, tidy, nCores, remote, mr)
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
    suppressWarnings({
      ## Bring them together
      tOut <- ungroup(tEst) %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        # Renaming, computing ratios, and SE
        mutate(POLE_AREA = pEst,
               MATURE_AREA = maEst,
               LATE_AREA = lEst,
               MOSAIC_AREA = moEst,
               AREA_TOTAL = aEst,
               ## Ratios
               POLE_PERC = POLE_AREA / AREA_TOTAL * 100,
               MATURE_PERC = MATURE_AREA / AREA_TOTAL * 100,
               LATE_PERC = LATE_AREA / AREA_TOTAL * 100,
               MOSAIC_PERC = MOSAIC_AREA / AREA_TOTAL * 100,
               ## Ratio Var
               ppVar = (1/AREA_TOTAL^2) * (pVar + (POLE_PERC^2 * aVar) - 2 * POLE_PERC * cvEst_p),
               mapVar = (1/AREA_TOTAL^2) * (maVar + (MATURE_PERC^2 * aVar) - 2 * MATURE_PERC * cvEst_ma),
               lpVar = (1/AREA_TOTAL^2) * (lVar + (LATE_PERC^2 * aVar) - 2 * LATE_PERC * cvEst_l),
               mopVar = (1/AREA_TOTAL^2) * (moVar + (MOSAIC_PERC^2 * aVar) - 2 * MOSAIC_PERC * cvEst_mo),
               ## SE RATIO
               POLE_PERC_SE = sqrt(ppVar) / POLE_PERC * 100,
               MATURE_PERC_SE = sqrt(mapVar) / MATURE_PERC * 100,
               LATE_PERC_SE = sqrt(lpVar) / LATE_PERC * 100,
               MOSAIC_PERC_SE = sqrt(mopVar) / MOSAIC_PERC * 100,
               ## VAR RATIO
               POLE_PERC_VAR = ppVar,
               MATURE_PERC_VAR = mapVar,
               LATE_PERC_VAR = lpVar,
               MOSAIC_PERC_VAR = mopVar,
               ## SE TOTAL
               POLE_AREA_SE = sqrt(pVar) / POLE_AREA * 100,
               MATURE_AREA_SE = sqrt(maVar) / MATURE_AREA * 100,
               LATE_AREA_SE = sqrt(lVar) / LATE_AREA * 100,
               MOSAIC_AREA_SE = sqrt(moVar) / MOSAIC_AREA * 100,
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL *100,
               ## VAR TOTAL
               POLE_AREA_VAR = pVar,
               MATURE_AREA_VAR = maVar,
               LATE_AREA_VAR = lVar,
               MOSAIC_AREA_VAR = moVar,
               AREA_TOTAL_VAR = aVar,
               ## nPlots
               nPlots_AREA = plotIn_AREA)
    })

    if (totals) {
      if (variance) {
        tOut <- tOut %>%
          select(grpBy,"POLE_PERC","MATURE_PERC", "LATE_PERC", "MOSAIC_PERC",
                 "POLE_AREA","MATURE_AREA","LATE_AREA","MOSAIC_AREA","AREA_TOTAL",
                 "POLE_PERC_VAR","MATURE_PERC_VAR","LATE_PERC_VAR",   "MOSAIC_PERC_VAR",
                 "POLE_AREA_VAR",   "MATURE_AREA_VAR", "LATE_AREA_VAR",   "MOSAIC_AREA_VAR",
                 "AREA_TOTAL_VAR","nPlots_AREA", 'N')
      } else {
        tOut <- tOut %>%
          select(grpBy,"POLE_PERC","MATURE_PERC", "LATE_PERC", "MOSAIC_PERC",
                 "POLE_AREA","MATURE_AREA","LATE_AREA","MOSAIC_AREA","AREA_TOTAL",
                 "POLE_PERC_SE","MATURE_PERC_SE","LATE_PERC_SE",   "MOSAIC_PERC_SE",
                 "POLE_AREA_SE",   "MATURE_AREA_SE", "LATE_AREA_SE",   "MOSAIC_AREA_SE",
                 "AREA_TOTAL_SE","nPlots_AREA")
      }

    } else {
      if (variance){
        tOut <- tOut %>%
          select(grpBy,"POLE_PERC","MATURE_PERC", "LATE_PERC", "MOSAIC_PERC",
                 "POLE_PERC_VAR","MATURE_PERC_VAR","LATE_PERC_VAR",   "MOSAIC_PERC_VAR",
                 "nPlots_AREA", 'N')
      } else {
        tOut <- tOut %>%
          select(grpBy,"POLE_PERC","MATURE_PERC", "LATE_PERC", "MOSAIC_PERC",
                 "POLE_PERC_SE","MATURE_PERC_SE","LATE_PERC_SE",   "MOSAIC_PERC_SE",
                 "nPlots_AREA")
      }

    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]


    ## Tidy things up if they didn't specify polys, returnSpatial
    if (tidy & is.null(polys) & returnSpatial == FALSE){
      ## Writing the variance options after the below, and don't want to change below
      ## So instead, do a temporary name swap
      if (variance) names(tOut) <- str_replace(names(tOut), '_VAR', '_SE')

      ## pivot longer
      stagePERC <- pivot_longer(select(tOut, grpBy, POLE_PERC:MOSAIC_PERC, nPlots_AREA), names_to = 'STAGE', values_to = 'PERC', cols = POLE_PERC:MOSAIC_PERC) %>%
        mutate(STAGE = str_split(STAGE,pattern= '_', simplify = TRUE,)[,1])
      stagePERC_SE <- pivot_longer(select(tOut, grpBy, POLE_PERC_SE:MOSAIC_PERC_SE), names_to = 'STAGE', values_to = 'PERC_SE', cols = POLE_PERC_SE:MOSAIC_PERC_SE) %>%
        mutate(STAGE = str_split(STAGE,pattern= '_', simplify = TRUE,)[,1])
      ## rejoin
      stage <- left_join(stagePERC, stagePERC_SE, by = c(grpBy, 'STAGE')) %>%
        select(grpBy, STAGE, PERC, PERC_SE, nPlots_AREA)

      if (totals){
        stageAREA <- pivot_longer(select(tOut, grpBy, POLE_AREA:MOSAIC_AREA), names_to = 'STAGE', values_to = 'AREA', cols = POLE_AREA:MOSAIC_AREA) %>%
          mutate(STAGE = str_split(STAGE,pattern= '_', simplify = TRUE,)[,1])
        stageAREA_SE <- pivot_longer(select(tOut, grpBy, POLE_AREA_SE:MOSAIC_AREA_SE), names_to = 'STAGE', values_to = 'AREA_SE', cols = POLE_AREA_SE:MOSAIC_AREA_SE) %>%
          mutate(STAGE = str_split(STAGE,pattern= '_', simplify = TRUE,)[,1])
        ## Rejoin
        stage <- stage %>%
          left_join(stageAREA, by = c(grpBy, 'STAGE')) %>%
          left_join(stageAREA_SE, by = c(grpBy, 'STAGE'))%>%
          select(grpBy, STAGE, PERC, AREA, PERC_SE, AREA_SE, nPlots_AREA)
      }
      tOut <- stage

      ## Writing the variance options after the below, and don't want to change below
      ## So instead, do a temporary name swap
      if (variance) names(tOut) <- str_replace(names(tOut), '_SE', '_VAR')

      }
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

