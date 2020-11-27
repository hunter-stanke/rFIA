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
                             totals = FALSE,
                             nCores = 1,
                             remote,
                             mr){


  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN',
                 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')

  ## If remote, read in state by state. Otherwise, drop all unnecessary tables
  db <- readRemoteHelper(db, remote, reqTables, nCores)

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
  if (treeType %in% c('live', 'dead', 'gs', 'all') == FALSE){
    stop('treeType must be one of: "live", "dead", "gs", or "all".')
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
  grpBy <- grpByToChar(db, grpBy_quo)

  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy

  ## Handle the modifier if it was given
  db$TREE <- db$TREE %>%
    mutate(TRE_CN = CN,
           state = !!stateVar,
           grp = !!grpVar)

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
  pops <- handlePops(db, evalType = c('EXPVOL', 'EXPCURR'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (method != 'TI') {
    pops <- mergeSmallStrata(db, pops)
  }




  ## Canned groups -------------------------------------------------------------

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
  db$PLOT <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA',
                               'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD',
                               all_of(grpP), 'aD_p', 'sp', 'COUNTYCD'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS',
                               'COND_STATUS_CD', 'CONDID',
                               all_of(grpC), 'aD_c', 'landD')) %>%
    filter(PLT_CN %in% db$PLOT$PLT_CN)
  db$TREE <- select(db$TREE, c('PLT_CN', 'TRE_CN', 'CONDID', 'DIA', 'grp',
                                'state', 'SUBP', 'TREE',
                                all_of(grpT), 'tD', 'typeD')) %>%
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
      out <- parLapply(cl, X = names(plts), fun = divHelper1, plts,
                       db[names(db) %in% c('COND', 'TREE')],
                       grpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = divHelper1, plts,
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
    full <- bind_rows(out[names(out) == 'full'])


    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)


    ## An iterator for population estimation
    popState <- split(pops, as.factor(pops$STATECD))
    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        out <- parLapply(cl, X = names(popState), fun = divHelper2, popState, t, full, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = divHelper2, popState, t, full, grpBy, method,  mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])
    full <- bind_rows(out[names(out) == 'full'])




    ## Compute moving average weights if not TI ----------------------------------
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){

      ## Compute the weights
      wgts <- maWeights(pops, method, lambda)

      ## If moving average ribbons, add lambda to grpBy for easier summary
      if (str_to_upper(method) == 'EMA' & length(lambda) > 1){
        grpBy <- c('lambda', grpBy)
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
        mutate(across(aEst:sEst, ~(.*wgt))) %>%
        mutate(across(aVar:cvEst_s, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize(across(aEst:plotIn_AREA, sum, na.rm = TRUE))


      ## If using an ANNUAL estimator --------------------------------------------
    } else if (str_to_upper(method) == 'ANNUAL') {

      ## ANNUAL ESTIMATOR is when END_INVYR = INVYR
      tEst <- tEst %>%
        group_by(INVYR, .dots = grpBy) %>%
        summarize(across(.cols = everything(),  sum, na.rm = TRUE)) %>%
        filter(YEAR == INVYR)%>%
        mutate(YEAR = INVYR)


      ## Rather than choose the annual panel estimate when INVYR = END_INVYR,
      ## choose the END_INVYR that has the highest N for each INVYR. Doing this
      ## because maybe not all 2018 data had been entered by the time the 2018
      ## END_INVYR cycle was produced. Maybe 2019 has more info on 2018 plots.
      ## So, ideally we would choose the cycle with the most plots for a given
      ## panel. Doing that here, important distinction from previous.
      ## NOT USED CURRENTLY --------------------------------------------------

      # tEst <- tEst %>%
      #   left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      #   group_by(STATECD, INVYR, .dots = aGrpBy[aGrpBy != 'STATECD']) %>%
      #   summarize(across(.cols = everything(),  sum, na.rm = TRUE)) %>%
      #   group_by(STATECD, INVYR, .dots = aGrpBy[aGrpBy %in% c('STATECD', 'YEAR') == FALSE]) %>%
      #   filter(plotIn_TREE == max(plotIn_TREE, na.rm = TRUE)) %>%
      #   filter(YEAR == max(YEAR, na.rm = TRUE)) %>%
      #   select(-c(YEAR)) %>%
      #   mutate(YEAR = INVYR)
    }


    out <- list(tEst = tEst, full = full, grpBy = grpBy, grpByOrig = grpByOrig)
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
                byPlot, totals, nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  tEst <- bind_rows(out[names(out) == 'tEst'])
  full <- bind_rows(out[names(out) == 'full'])
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
      tOut <- tEst %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        mutate(H_a = hEst / aEst,
               Eh_a = ehEst / aEst,
               S_a = sEst / aEst,
               AREA_TOTAL = aEst,
               ## Ratio variance
               haVar = (1/AREA_TOTAL^2) * (hVar + (H_a^2 * aVar) - 2 * H_a * cvEst_h),
               ehaVar = (1/AREA_TOTAL^2) * (ehVar + (Eh_a^2 * aVar) - (2 * Eh_a * cvEst_eh)),
               saVar = (1/AREA_TOTAL^2) * (sVar + (S_a^2 * aVar) - 2 * S_a * cvEst_s),
               H_a_VAR = haVar,
               Eh_a_VAR = ehaVar,
               S_a_VAR = saVar,
               ## RATIO SE
               H_a_SE = sqrt(haVar) / H_a * 100,
               Eh_a_SE = sqrt(ehaVar) / Eh_a * 100,
               S_a_SE = sqrt(saVar) / S_a * 100,
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
               AREA_TOTAL_VAR = aVar,
               nPlots = plotIn_AREA)
    })

    fullGrps <- syms(grpBy[!c(grpBy %in% 'lambda')])

    ### Up a few spatial scales
    full <- full %>%
      lazy_dt() %>%
      group_by(!!!fullGrps) %>%
      summarize(H_g = divIndex(grp, state, index = 'H'),
                Eh_g = divIndex(grp, state, index = 'Eh'),
                S_g = divIndex(grp, state, index = 'S')) %>%
    as.data.frame()

    tOut <- tOut %>%
      left_join(full, by = grpBy[!c(grpBy %in% 'lambda')]) %>%
      mutate(H_b = H_g - H_a,
             Eh_b = Eh_g - Eh_a,
             S_b = S_g - S_a)


    if (totals) {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, H_a, H_b, H_g, Eh_a, Eh_b, Eh_g, S_a, S_b, S_g, AREA_TOTAL, H_a_VAR, Eh_a_VAR, S_a_VAR, AREA_TOTAL_VAR, nPlots, N)

      } else {
        tOut <- tOut %>%
          select(grpBy, H_a, H_b, H_g, Eh_a, Eh_b, Eh_g, S_a, S_b, S_g, AREA_TOTAL, H_a_SE, Eh_a_SE, S_a_SE, AREA_TOTAL_SE, nPlots)
      }

    } else {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, H_a, H_b, H_g, Eh_a, Eh_b, Eh_g, S_a, S_b, S_g,  H_a_VAR, Eh_a_VAR, S_a_VAR,  nPlots, N)

      } else {
        tOut <- tOut %>%
          select(grpBy, H_a, H_b, H_g, Eh_a, Eh_b, Eh_g, S_a, S_b, S_g,  H_a_SE, Eh_a_SE, S_a_SE,  nPlots)
      }
    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

  }

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



