standStructNew <- function(db,
                          grpBy = NULL,
                          polys = NULL,
                          returnSpatial = FALSE,
                          landType = 'forest',
                          areaDomain = NULL,
                          byPlot = FALSE,
                          totals = FALSE,
                          tidy = TRUE,
                          SE = TRUE,
                          nCores = 1) {

  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        0
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        select(!!grpBy_quo)
    )
    # If column doesnt exist, just returns 0, not a dataframe
    if (is.null(nrow(d_quo))){
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT, TREE, or COND tables. Did you accidentally quote the variables names? e.g. use grpBy = ECOSUBCD (correct) instead of grpBy = "ECOSUBCD". ', collapse = ', '))
    } else {
      # Convert to character
      grpBy <- names(d_quo)
    }
  }

  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(first(class(polys))) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (tidy & returnSpatial & !is.null(polys)){
    warning('Returning multiple observations for each areal unit. If returnSpatial = TRUE, tidy = FALSE is recommended.')
  }
  if (landType %in% c('timber', 'forest', 'all') == FALSE){
    stop('landType must be one of: "forest", "timber", or "all".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }

  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}
  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy

  ## IF the object was clipped
  if ('prev' %in% names(db$PLOT)){
    ## Only want the current plots, no grm
    db$PLOT <- filter(db$PLOT, prev == 0)
  }

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    ## A unique ID
    polys$polyID <- 1:nrow(polys)

    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)

    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('LON', 'LAT', pltID)) %>%
      filter(!is.na(LAT) & !is.na(LON)) %>%
      distinct(pltID, .keep_all = TRUE)
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)

    ## Split up polys
    polyList <- split(polys, as.factor(polys$polyID))
    suppressWarnings({suppressMessages({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        cl <- makeCluster(nCores)
        clusterEvalQ(cl, {
          library(dplyr)
          library(stringr)
          library(rFIA)
        })
        out <- parLapply(cl, X = names(polyList), fun = areal_par, pltSF, polyList)
        #stopCluster(cl) # Keep the cluster active for the next run
      } else { # Unix systems
        out <- mclapply(names(polyList), FUN = areal_par, pltSF, polyList, mc.cores = nCores)
      }
    })})
    pltSF <- bind_rows(out)

    # A warning
    if (length(unique(pltSF$pltID)) < 1){
      stop('No plots in db overlap with polys.')
    }
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'pltID')

    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        right_join(select(db$PLOT, PLT_CN, pltID), by = pltID) %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    grpBy <- c(grpBy, 'LON', 'LAT')
  } # END AREAL

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
  } else if (tolower(landType) == 'all') {
    db$COND$landD <- 1
  }

  # update spatial domain indicator
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(db$PLOT$pltID %in% pltSF$pltID, 1, 0)
  } else {
    db$PLOT$sp <- 1
  }

  # User defined domain indicator for area (ex. specific forest type)
  pcEval <- left_join(db$PLOT, select(db$COND, -c('STATECD', 'UNITCD', 'COUNTYCD', 'INVYR', 'PLOT')), by = 'PLT_CN')
  areaDomain <- substitute(areaDomain)
  pcEval$aD <- eval(areaDomain, pcEval) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(pcEval$aD)) pcEval$aD[is.na(pcEval$aD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(pcEval$aD)) pcEval$aD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  pcEval$aD <- as.numeric(pcEval$aD)
  db$COND <- left_join(db$COND, select(pcEval, c('PLT_CN', 'CONDID', 'aD')), by = c('PLT_CN', 'CONDID')) %>%
    mutate(aD_c = aD)
  aD_p <- pcEval %>%
    group_by(PLT_CN) %>%
    summarize(aD_p = as.numeric(any(aD > 0)))
  db$PLOT <- left_join(db$PLOT, aD_p, by = 'PLT_CN')
  rm(pcEval)

  ### Snag the EVALIDs that are needed
  db$POP_EVAL<- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPCURR') %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)# %>%
  #group_by(END_INVYR) %>%
  #summarise(id = list(EVALID)

  ## Make an annual panel ID, associated with an INVYR

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'END_INVYR')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
    mutate_if(is.factor,
              as.character)

  ### Which estimator to use?
  if (str_to_upper(method) %in% c('ANNUAL', "SMA", 'EMA')){
    ## Keep an original
    popOrig <- pops
    ## Want to use the year where plots are measured, no repeats
    ## Breaking this up into pre and post reporting becuase
    ## Estimation units get weird on us otherwise
    #pops <- distinct(pops, INVYR, PLT_CN, .keep_all = TRUE)
    pops <- pops %>%
      group_by(STATECD) %>%
      filter(END_INVYR == INVYR)

    prePops <- popOrig %>%
      group_by(STATECD) %>%
      filter(INVYR < min(END_INVYR, na.rm = TRUE)) %>%
      distinct(PLT_CN, .keep_all = TRUE)

    pops <- bind_rows(pops, prePops)

    ## P2POINTCNT column is NOT consistent for annnual estimates, plots
    ## within individual strata and est units are related to different INVYRs
    p2 <- pops %>%
      group_by(STRATUM_CN, INVYR) %>%
      summarize(P2POINTCNT = length(unique(PLT_CN)))
    ## Rejoin
    pops <- pops %>%
      mutate(YEAR = INVYR) %>%
      select(-c(P2POINTCNT)) %>%
      left_join(p2, by = c('STRATUM_CN', 'YEAR' = 'INVYR'))

  } else {     # Otherwise temporally indifferent
    pops <- rename(pops, YEAR = END_INVYR)
  }

  ## Want a count of p2 points / eu, gets screwed up with grouping below
  p2eu <- pops %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, P2POINTCNT) %>%
    group_by(ESTN_UNIT_CN) %>%
    summarize(p2eu = sum(P2POINTCNT, na.rm = TRUE))

  ## Rejoin
  pops <- left_join(pops, p2eu, by = 'ESTN_UNIT_CN')


  ## Recode a few of the estimation methods to make things easier below
  pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')



  ## Only the necessary plots for EVAL of interest
  db$PLOT <- filter(db$PLOT, PLT_CN %in% pops$PLT_CN)

  ## Merging state and county codes
  plts <- split(db$PLOT, as.factor(paste(db$PLOT$COUNTYCD, db$PLOT$STATECD, sep = '_')))
  #plts <- split(db$PLOT, as.factor(db$PLOT$STATECD))

  suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      clusterEvalQ(cl, {
        library(dplyr)
        library(stringr)
        library(rFIA)
      })
      out <- parLapply(cl, X = names(plts), fun = ssHelper1, plts, db, grpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = ssHelper1, plts, db, grpBy, byPlot, mc.cores = nCores)
    }
  })

  ## back to dataframes
  out <- unlist(out, recursive = FALSE)
  t <- bind_rows(out[names(out) == 't'])


  ## Adding YEAR to groups
  grpBy <- c('YEAR', grpBy)


  ## Splitting up by ESTN_UNIT
  popState <- split(pops, as.factor(pops$STATECD))

  suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      ## Use the same cluster as above
      # cl <- makeCluster(nCores)
      # clusterEvalQ(cl, {
      #   library(dplyr)
      #   library(stringr)
      #   library(rFIA)
      # })
      out <- parLapply(cl, X = names(popState), fun = ssHelper2, popState, t, grpBy)
      stopCluster(cl)
    } else { # Unix systems
      out <- mclapply(names(popState), FUN = ssHelper2, popState, t, grpBy, mc.cores = nCores)
    }
  })
  ## back to dataframes
  out <- unlist(out, recursive = FALSE)
  tEst <- bind_rows(out[names(out) == 'tEst'])



}















#
#
#   ## Which grpByNames are in which table? Helps us subset below
#   grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
#   grpC <- names(db$COND)[names(db$COND) %in% grpBy]
#
#   ### Snag the EVALIDs that are needed
#   ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
#   ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
#   if (!is.null(polys)){
#     ids <- pltSF %>%
#       left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
#       left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
#       inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
#       filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
#       filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
#       distinct(polyID, END_INVYR, EVALID)
#     ## Must be a most recent subset for mergeYears to exists
#     ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
#     ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
#     if (exists('mergeYears')){
#       ids <- ids %>%
#         left_join(mergeYears, by = 'polyID') %>%
#         mutate(END_INVYR = maxYear)
#     }
#     ## Snag all the EVALIDs for each poly
#     ids <- ids %>%
#       group_by(polyID, END_INVYR) %>%
#       summarise(id = list(EVALID))
#   } else {
#     ids <- db$POP_EVAL %>%
#       select('CN', 'END_INVYR', 'EVALID') %>%
#       inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
#       filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
#       filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
#       distinct(END_INVYR, EVALID) %>%
#       group_by(END_INVYR) %>%
#       summarise(id = list(EVALID))
#   }
#
#   ## Add a progress bar with ETA
#   pb <- progress_bar$new(
#     format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
#     total = nrow(ids), clear = FALSE, width= 100)
#
#   # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
#   out <- list()
#   for (y in 1:nrow(ids)){
#     ## Clip out the necessary data
#     if (!is.null(polys)){
#       ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
#       ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
#       ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
#       if (y == 1){
#         ## First iteration, do both spatial and temporal
#         db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
#         db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
#       } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
#         ## Just rerun temporal clip with the same spatially clipped object
#         db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
#       } else {
#         ## Hit a new poly, make a new spatial clip
#         db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
#         db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
#       }
#       # update spatial domain indicator
#       db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
#
#       ## Polys not specified, just temporal
#     } else {
#       db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
#       # update spatial domain indicator
#       db_clip$PLOT$sp <- 1
#     }
#
#     ## Prep joins and filters
#     data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', grpP, 'aD_p', 'sp')) %>%
#       left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', 'landD', 'aD_c', grpC)), by = c('PLT_CN')) %>%
#       left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
#       left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
#       left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
#       right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
#       #left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
#       #left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
#       left_join(select(db_clip$TREE, c('PLT_CN', 'CONDID', 'DIA', 'STATUSCD', 'CCLCD', 'TREECLCD', 'STANDING_DEAD_CD', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE')), by = c('PLT_CN', 'CONDID')) %>%
#       mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
#       mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
#       rename(YEAR = END_INVYR,
#              YEAR_RANGE = REPORT_YEAR_NM) %>%
#       mutate_if(is.factor,
#                 as.character)%>%
#       filter(!is.na(YEAR)) %>%
#       distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE)
#
#     ## Recode a few of the estimation methods to make things easier below
#     data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
#                               `Post-Stratification` = 'strat',
#                               `Stratified random sampling` = 'strat',
#                               `Double sampling for stratification` = 'double',
#                               `Simple random sampling` = 'simple',
#                               `Subsampling units of unequal size` = 'simple')
#
#     # Test if any polygons cross state boundaries w/ different recent inventory years
#     if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
#       # Replace YEAR from above w/ max year so that data is pooled across states
#       data <- left_join(data, mergeYears, by = 'polyID') %>%
#         select(-c(YEAR)) %>%
#         mutate(YEAR = maxYear)
#     }
#
#     ## Comprehensive indicator function
#     data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
#     data$tDI <- data$landD * data$aD_p * data$aD_c * data$sp
#
#
#
#     ####################  COMPUTE ESTIMATES  ###########################
#     ### -- BYPLOT -- TPA Estimates at each plot location
#     if (byPlot) {
#       sOut <- data %>%
#         mutate(YEAR = INVYR) %>%
#         distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
#         group_by(.dots = grpBy, PLT_CN) %>%
#         summarize(stage = structHelper(DIA, CCLCD),
#                   nStems = length(which(tDI == 1)))
#
#       if (returnSpatial){
#         sOut <- sOut %>%
#           filter(!is.na(LAT) & !is.na(LON)) %>%
#           st_as_sf(coords = c('LON', 'LAT'),
#                    crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
#
#       }
#
#       ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
#     } else {
#       ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
#       if(is.null(polys)){
#         combos <- select(data, c(grpBy)) %>%
#           as.data.frame() %>%
#           group_by(.dots = grpBy) %>%
#           summarize() %>%
#           filter(!is.na(YEAR))
#       } else {
#         ## Non spatial combos
#         combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
#           as.data.frame() %>%
#           group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
#           summarize() %>%
#           filter(!is.na(YEAR))
#         combosNSpoly <- combosNS %>%
#           mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])
#
#         # New combos for spatial objects
#         combos <- pltSF %>%
#           select(-c(PLT_CN)) %>%
#           distinct(polyID, .keep_all = TRUE) %>%
#           right_join(combosNSpoly, by = 'polyID')
#       }
#       # List of rows for lapply
#       combos <- split(combos, seq(nrow(combos)))
#       suppressWarnings({
#         ## Compute estimates in parallel -- Clusters in windows, forking otherwise
#         if (Sys.info()['sysname'] == 'Windows'){
#           cl <- makeCluster(nCores)
#           clusterEvalQ(cl, {
#             library(dplyr)
#             library(stringr)
#             library(tidyr)
#           })
#           sOut <- parLapply(cl, X = names(combos), fun = standStructHelper, combos, data, grpBy, tidy, totals, SE)
#           stopCluster(cl)
#         } else { # Unix systems
#           sOut <- mclapply(X = names(combos), FUN = standStructHelper, combos, data, grpBy, totals, tidy, SE, mc.cores = nCores)
#         }
#       })
#
#       if (SE){
#         # Convert from list to dataframe
#         sOut <- do.call(rbind,sOut) %>% #bind_rows(sOut, .id = NULL) %>%
#           as.data.frame()
#
#         ## IF the user wants a tidy dataframe at the end, handle it for them
#         if (tidy){
#           # Gather up all those rando columns
#           stage <- gather(sOut, key = 'STAGE', value = 'PERC_AREA', POLE_PERC:MOSAIC_PERC)
#           stageSE <- gather(sOut, key = 'STAGE', value = 'PERC_AREA_SE', POLE_PERC_SE:MOSAIC_PERC_SE)
#           # Join them back up all nice like
#           sTidy <- bind_cols(select(stage, c(names(combos[[1]]), 'STAGE', 'PERC_AREA'), nPlots),
#                              select(stageSE, PERC_AREA_SE))
#           if(totals){
#             stageT <- gather(sOut, key = 'STAGE', value = 'AREA', POLE_AREA:MOSAIC_AREA)
#             stageTSE <- gather(sOut, key = 'STAGE', value = 'AREA_SE', POLE_AREA_SE:MOSAIC_AREA_SE)
#             # Join them back up all nice like
#             sTidy <- bind_cols(sTidy,
#                                select(stageT, AREA),
#                                select(stageTSE, AREA_SE))
#           }
#           sOut <- sTidy %>%
#             select(-nPlots, nPlots) %>%
#             mutate(STAGE = str_split(STAGE, "_", simplify = TRUE)[,1]) %>%
#             arrange(YEAR)
#         }
#       } else {
#         # Pull out dataframe
#         sOut <- sOut[[1]] %>%
#           ungroup() %>%
#           as.data.frame()
#
#         if (tidy){
#           # Gather up all those rando columns
#           stage <- gather(sOut, key = 'STAGE', value = 'PERC_AREA', POLE_PERC:MOSAIC_PERC)
#           # Join them back up all nice like
#           sTidy <- bind_cols(select(stage, c(names(combos[[1]]), 'STAGE', 'PERC_AREA'), nPlots))
#           if(totals){
#             stageT <- gather(sOut, key = 'STAGE', value = 'AREA', POLE_AREA:MOSAIC_AREA)
#             # Join them back up all nice like
#             sTidy <- bind_cols(sTidy,
#                                select(stageT, AREA))
#           }
#           sOut <- sTidy %>%
#             select(-nPlots, nPlots) %>%
#             mutate(STAGE = str_split(STAGE, "_", simplify = TRUE)[,1]) %>%
#             arrange(YEAR)
#         }
#       }
#
#       # Names for below
#       sNames <- names(sOut)[names(sOut) %in% grpBy == FALSE]
#
#       # Return a spatial object
#       if ('YEAR' %in% names(sOut)){
#         # Return a spatial object
#         if (!is.null(polys) & returnSpatial) {
#           suppressMessages({suppressWarnings({sOut <- left_join(polys, sOut) %>%
#             select(c(grpByOrig, sNames, names(polys))) %>%
#             filter(!is.na(polyID))})})
#         } else if (!is.null(polys) & returnSpatial == FALSE){
#           sOut <- select(sOut, c(grpByOrig, sNames, everything())) %>%
#             filter(!is.na(polyID))
#         }
#       } else { ## Function found no plots within the polygon, so it panics
#         combos <- data %>%
#           as.data.frame() %>%
#           group_by(.dots = grpBy) %>%
#           summarize()
#         sOut <- data.frame("YEAR" = combos$YEAR,
#                            'POLE_PERC' = rep(NA, nrow(combos)), 'MATURE_PERC' = rep(NA, nrow(combos)),
#                            'LATE_PERC' = rep(NA, nrow(combos)),
#                            'MOSAIC_PERC' = rep(NA, nrow(combos)), 'POLE_PERC_SE' = rep(NA, nrow(combos)),
#                            'MATURE_PERC_SE' = rep(NA, nrow(combos)),
#                            'LATE_PERC_SE' = rep(NA, nrow(combos)),
#                            'MOSAIC_PERC_SE' = rep(NA, nrow(combos)),
#                            "nPlots" = rep(NA, nrow(combos)))
#         if (!is.null(polys) & returnSpatial) {
#           suppressMessages({suppressWarnings({
#             polys = left_join(polys, combos)
#             sOut <- left_join(polys, sOut)})})
#         } else if (!is.null(polys) & returnSpatial == FALSE){
#           sOut <- data.frame(select(sOut, -c('YEAR')), combos)
#         }
#       }
#
#
#     } # End byPlot = FALSE
#     out[[y]] <- sOut
#     pb$tick()
#   }
#   sOut <- do.call(rbind, out)
#   ## For spatial plots
#   if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
#   sOut <- drop_na(sOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
#     arrange(YEAR) %>%
#     as_tibble()
#
#   ## Above converts to tibble
#   if (returnSpatial) sOut <- st_sf(sOut)
#
#   # ## remove any duplicates in byPlot (artifact of END_INYR loop)
#   if (byPlot) sOut <- unique(sOut)
#
#   return(sOut)
# }
