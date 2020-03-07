sustIndex <- function(db,
                      grpBy = NULL,
                      polys = NULL,
                      returnSpatial = FALSE,
                      bySpecies = FALSE,
                      bySizeClass = FALSE,
                      landType = 'forest',
                      treeType = 'live',
                      minLive = 0,
                      method = 'annual',
                      lambda = .5,
                      treeDomain = NULL,
                      areaDomain = NULL,
                      totals = TRUE,
                      byPlot = FALSE,
                      nCores = 1) {

  ## Need a plotCN
  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[10,] %>% # Just the first row
        left_join(select(db$COND, PLT_CN, names(db$COND)[names(db$COND) %in% names(db$PLOT) == FALSE]), by = 'PLT_CN') %>%
        inner_join(select(db$TREE, PLT_CN, names(db$TREE)[names(db$TREE) %in% c(names(db$PLOT), names(db$COND)) == FALSE]), by = 'PLT_CN') %>%
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

  reqTables <- c('PLOT', 'TREE', 'TREE_GRM_COMPONENT', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method, 'unknown. Defaulting to Temporally Indifferent (TI).'))
  }

  # ## No EXP_GROW available for Western States, make sure we warn that values will be returned as 0
  # # These states do not allow temporal queries. Things are extremely weird with their eval groups
  # noGrow <- c(02,03,04,07,08,11,14,15,16, 23, 30, 32, 35,43,49, 78)
  # if(any(unique(db$PLOT$STATECD) %in% noGrow)){
  #   vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
  #   fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
  #   warning(paste('Recruitment data unavailable for: ', toString(fancyName) , '. Returning 0 for all recruitment estimates which include these states.', sep = ''))
  # }
  # # These states do not allow change estimates.
  # if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
  #   vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
  #   fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
  #   stop(paste('Growth & Mortality Estimates unavailable for: ', paste(as.character(fancyName), collapse = ', '), sep = ''))
  # }



  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}
  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy

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
    #grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
    grpBy = c(grpBy, 'polyID')

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
        right_join(select(db$PLOT, PLT_CN, pltID), by = 'pltID') %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  }
  if (byPlot & returnSpatial){
    grpBy <- c(grpBy, 'LON', 'LAT')
  } # END AREAL

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
  }

  # ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # # Land type domain indicator
  # if (tolower(landType) == 'forest'){
  #   db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
  #   # Tree Type domain indicator
  #   if (tolower(treeType) == 'live'){
  #     db$TREE$typeD <- 1
  #     ## Rename some variables in grm
  #     db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
  #                                     TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_FOREST,
  #                                     TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_FOREST,
  #                                     TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_FOREST,
  #                                     SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_FOREST,
  #                                     COMPONENT = SUBP_COMPONENT_AL_FOREST) %>%
  #       mutate(TPARECR_UNADJ = case_when(
  #         is.na(COMPONENT) ~ NA_real_,
  #         COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2') ~ TPAGROW_UNADJ,
  #         TRUE ~ 0))
  #
  #   } else if (tolower(treeType) == 'gs'){
  #     db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
  #     db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
  #                                     TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_FOREST,
  #                                     TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_FOREST,
  #                                     TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_FOREST,
  #                                     SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_FOREST,
  #                                     COMPONENT = SUBP_COMPONENT_GS_FOREST)%>%
  #       mutate(TPARECR_UNADJ = case_when(
  #         is.na(COMPONENT) ~ NA_real_,
  #         COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2') ~ TPAGROW_UNADJ,
  #         TRUE ~ 0))
  #   }
  # } else if (tolower(landType) == 'timber'){
  #   db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
  #   # Tree Type domain indicator
  #   if (tolower(treeType) == 'live'){
  #     db$TREE$typeD <- 1
  #     ## Rename some variables in grm
  #     db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
  #                                     TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_TIMBER,
  #                                     TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_TIMBER,
  #                                     TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_TIMBER,
  #                                     SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_TIMBER,
  #                                     COMPONENT = SUBP_COMPONENT_AL_TIMBER)%>%
  #       mutate(TPARECR_UNADJ = case_when(
  #         is.na(COMPONENT) ~ NA_real_,
  #         COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2') ~ TPAGROW_UNADJ,
  #         TRUE ~ 0))
  #
  #   } else if (tolower(treeType) == 'gs'){
  #     db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
  #     db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
  #                                     TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_TIMBER,
  #                                     TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_TIMBER,
  #                                     TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_TIMBER,
  #                                     SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_TIMBER,
  #                                     COMPONENT = SUBP_COMPONENT_GS_TIMBER)%>%
  #       mutate(TPARECR_UNADJ = case_when(
  #         is.na(COMPONENT) ~ NA_real_,
  #         COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2') ~ TPAGROW_UNADJ,
  #         TRUE ~ 0))
  #   }
  # }


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

  # Same as above for tree (ex. trees > 20 ft tall)
  treeDomain <- substitute(treeDomain)
  tD <- eval(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$TREE$tD <- as.numeric(tD)


  ### Snag the EVALIDs that are needed
  db$POP_EVAL  <- db$POP_EVAL %>%
    #left_join(ga, by = 'END_INVYR') %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD', 'GROWTH_ACCT') %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP %in% c('EXPVOL')) %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)# %>%
  #group_by(END_INVYR) %>%
  #summarise(id = list(EVALID)

  ## Make an annual panel ID, associated with an INVYR

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'GROWTH_ACCT', 'END_INVYR', 'EVAL_TYP')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
    ## Join on REMPER PLOTS
    left_join(select(db$PLOT, PLT_CN, REMPER, PREV_PLT_CN, DESIGNCD, PLOT_STATUS_CD), by = 'PLT_CN') %>%
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & DESIGNCD == 1 & PLOT_STATUS_CD != 3) %>%
    mutate_if(is.factor,
              as.character)

  ### Which estimator to use?
  if (str_to_upper(method) %in% c('ANNUAL')){
    ## Want to use the year where plots are measured, no repeats
    ## Breaking this up into pre and post reporting becuase
    ## Estimation units get weird on us otherwise
    popOrig <- pops
    pops <- pops %>%
      group_by(STATECD) %>%
      filter(END_INVYR == INVYR) %>%
      ungroup()

    prePops <- popOrig %>%
      group_by(STATECD) %>%
      filter(INVYR < min(END_INVYR, na.rm = TRUE)) %>%
      distinct(PLT_CN, .keep_all = TRUE) %>%
      ungroup()

    pops <- bind_rows(pops, prePops) %>%
      mutate(YEAR = INVYR)

  } else {     # Otherwise temporally indifferent
    pops <- rename(pops, YEAR = END_INVYR)
  }

  ## P2POINTCNT column is NOT consistent for annnual estimates, plots
  ## within individual strata and est units are related to different INVYRs
  p2_INVYR <- pops %>%
    group_by(ESTN_UNIT_CN, STRATUM_CN, INVYR) %>%
    summarize(P2POINTCNT_INVYR = length(unique(PLT_CN)))
  ## Want a count of p2 points / eu, gets screwed up with grouping below
  p2eu_INVYR <- p2_INVYR %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, INVYR, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN, INVYR) %>%
    summarize(p2eu_INVYR = sum(P2POINTCNT_INVYR, na.rm = TRUE))
  p2eu <- pops %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN) %>%
    summarize(p2eu = sum(P2POINTCNT, na.rm = TRUE))

  ## Rejoin
  pops <- pops %>%
    left_join(p2_INVYR, by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'INVYR')) %>%
    left_join(p2eu_INVYR, by = c('ESTN_UNIT_CN', 'INVYR')) %>%
    left_join(p2eu, by = 'ESTN_UNIT_CN')


  ## Recode a few of the estimation methods to make things easier below
  pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')


  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
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


  # # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
  # if (!is.null(polys)){
  #   aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
  # } else {
  #   aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
  # }

  ## Only the necessary plots for EVAL of interest
  db$PLOT <- filter(db$PLOT, PLT_CN %in% pops$PLT_CN | PLT_CN %in% pops$PREV_PLT_CN)

  ## Reduce the memory load for others
  db <- clipFIA(db, mostRecent = FALSE)

  ## Merging state and county codes
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
      out <- parLapply(cl, X = names(plts), fun = sustIndexHelper1, plts, db, grpBy, byPlot, minLive)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = sustIndexHelper1, plts, db, grpBy, byPlot, minLive, mc.cores = nCores)
    }
  })


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
    ## Population estimation
  } else {
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    #a <- bind_rows(out[names(out) == 'a'])
    t <- bind_rows(out[names(out) == 't'])


    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)
    #aGrpBy <- c('YEAR', aGrpBy)


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
        out <- parLapply(cl, X = names(popState), fun = sustIndexHelper2, popState, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = sustIndexHelper2, popState, t, grpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])


    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
      ### ---- SIMPLE MOVING AVERAGE
      if (str_to_upper(method) == 'SMA'){
        ## Assuming a uniform weighting scheme
        wgts <- pops %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(wgt = 1 / length(unique(INVYR)))

        #aEst <- left_join(aEst, wgts, by = 'ESTN_UNIT_CN')
        tEst <- left_join(tEst, wgts, by = 'ESTN_UNIT_CN')

        #### ----- Linear MOVING AVERAGE
      } else if (str_to_upper(method) == 'LMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))

        ## Want a number of INVYRs per EU
        neu <- wgts %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(n = sum(rank, na.rm = TRUE))

        ## Rejoining and computing wgts
        wgts <- wgts %>%
          left_join(neu, by = 'ESTN_UNIT_CN') %>%
          mutate(wgt = rank / n) %>%
          ungroup() %>%
          select(ESTN_UNIT_CN, INVYR, wgt)

        #aEst <- left_join(aEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

        #### ----- EXPONENTIAL MOVING AVERAGE
      } else if (str_to_upper(method) == 'EMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))


        if (length(lambda) < 2){
          ## Want sum of weighitng functions
          neu <- wgts %>%
            mutate(l = lambda) %>%
            group_by(ESTN_UNIT_CN) %>%
            summarize(l = first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = 'ESTN_UNIT_CN') %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(ESTN_UNIT_CN, INVYR, wgt)
        } else {
          grpBy <- c('lambda', grpBy)
          #aGrpBy <- c('lambda', aGrpBy)
          ## Duplicate weights for each level of lambda
          yrWgts <- list()
          for (i in 1:length(unique(lambda))) {
            yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
          }
          wgts <- bind_rows(yrWgts)
          ## Want sum of weighitng functions
          neu <- wgts %>%
            group_by(lambda, ESTN_UNIT_CN) %>%
            summarize(l = first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = c('lambda', 'ESTN_UNIT_CN')) %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(lambda, ESTN_UNIT_CN, INVYR, wgt)
        }

        #aEst <- left_join(aEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

      }

      ### Applying the weights
      # Area
      # aEst <- aEst %>%
      #   mutate_at(vars(aEst), ~(.*wgt)) %>%
      #   mutate_at(vars(aVar), ~(.*(wgt^2))) %>%
      #   group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
      #   summarize_at(vars(aEst:plotIn_AREA), sum, na.rm = TRUE)


      tEst <- tEst %>%
        mutate_at(vars(ctEst:silvEst), ~(.*wgt)) %>%
        mutate_at(vars(ctVar:cvEst_silv), ~(.*(wgt^2))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize_at(vars(ctEst:plotIn_t), sum, na.rm = TRUE)

    }

    ##---------------------  TOTALS and RATIOS
    # Tree
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    ##---------------------  TOTALS and RATIOS
    suppressWarnings({
      tOut <- tTotal %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        mutate(CHNG_TPA = ctEst,
               CHNG_BAA = cbEst,
               PREV_TPA = ptEst,
               PREV_BAA = pbEst,
               TPA_RATE = ctEst / ptEst,
               BAA_RATE = cbEst / pbEst,
               INSECT_RATE = bugEst / pbEst,
               DISEASE_RATE = diseaseEst / pbEst,
               FIRE_RATE = fireEst / pbEst,
               ANIMAL_RATE = animalEst / pbEst,
               WEATHER_RATE = weatherEst / pbEst,
               VEG_RATE = vegEst / pbEst,
               UNKNOWN_RATE = unEst / pbEst,
               SILV_RATE = silvEst / pbEst,
               x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
               y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
               M = sqrt(x^2 + y^2),
               SUST_INDEX = if_else(x < 0, -M, M),
               ## TOTAL SE
               CHNG_TPA_SE = sqrt(ctVar) / abs(ctEst) * 100,
               CHNG_BAA_SE = sqrt(cbVar) / abs(cbEst) * 100,
               PREV_TPA_SE = sqrt(ptVar) / abs(ptEst) * 100,
               PREV_BAA_SE = sqrt(pbVar) / abs(pbEst) * 100,
               ## Ratio variance
               ctVar = (1/PREV_TPA^2) * (ctVar + (TPA_RATE^2 * ptVar) - 2 * TPA_RATE * cvEst_ct),
               cbVar = (1/PREV_BAA^2) * (cbVar + (BAA_RATE^2 * pbVar) - 2 * BAA_RATE * cvEst_cb),
               bugVar = (1/PREV_BAA^2) * (bugVar + (INSECT_RATE^2 * pbVar) - 2 * INSECT_RATE * cvEst_bug),
               diseaseVar = (1/PREV_BAA^2) * (diseaseVar + (DISEASE_RATE^2 * pbVar) - 2 * DISEASE_RATE * cvEst_disease),
               fireVar = (1/PREV_BAA^2) * (fireVar + (FIRE_RATE^2 * pbVar) - 2 * FIRE_RATE * cvEst_fire),
               animalVar = (1/PREV_BAA^2) * (animalVar + (ANIMAL_RATE^2 * pbVar) - 2 * ANIMAL_RATE * cvEst_animal),
               weatherVar = (1/PREV_BAA^2) * (weatherVar + (WEATHER_RATE^2 * pbVar) - 2 * WEATHER_RATE * cvEst_weather),
               vegVar = (1/PREV_BAA^2) * (vegVar + (VEG_RATE^2 * pbVar) - 2 * VEG_RATE * cvEst_veg),
               unVar = (1/PREV_BAA^2) * (unVar + (UNKNOWN_RATE^2 * pbVar) - 2 * UNKNOWN_RATE * cvEst_un),
               silvVar = (1/PREV_BAA^2) * (silvVar + (SILV_RATE^2 * pbVar) - 2 * SILV_RATE * cvEst_silv),
               ## RATIO SE
               TPA_RATE_SE = sqrt(ctVar) / abs(TPA_RATE) * 100,
               BAA_RATE_SE = sqrt(cbVar) / abs(BAA_RATE) * 100,
               INSECT_RATE_SE = sqrt(bugVar) / abs(INSECT_RATE) * 100,
               DISEASE_RATE_SE = sqrt(diseaseVar) / abs(DISEASE_RATE) * 100,
               FIRE_RATE_SE = sqrt(fireVar) / abs(FIRE_RATE) * 100,
               ANIMAL_RATE_SE = sqrt(animalVar) / abs(ANIMAL_RATE) * 100,
               WEATHER_RATE_SE = sqrt(weatherVar) / abs(WEATHER_RATE) * 100,
               VEG_RATE_SE = sqrt(vegVar) / abs(VEG_RATE) * 100,
               UNKNOWN_RATE_SE = sqrt(unVar) / abs(UNKNOWN_RATE) * 100,
               SILV_RATE_SE = sqrt(silvVar) / abs(SILV_RATE) * 100,
               #SUST_INDEX_SE = sqrt(Mvar) / abs(SUST_INDEX) * 100,
               nPlots = plotIn_t,
               nTotal = nh,
               TPA_RATE_INT = abs(TPA_RATE) * 1.96 * TPA_RATE_SE / 100,
               BAA_RATE_INT = abs(BAA_RATE) * 1.96 * BAA_RATE_SE / 100,
               TPA_STATUS = case_when(
                 TPA_RATE < 0 & TPA_RATE + TPA_RATE_INT < 0 ~ 'Decline',
                 TPA_RATE < 0 & TPA_RATE + TPA_RATE_INT > 0 ~ 'Stable',
                 TPA_RATE > 0 & TPA_RATE - TPA_RATE_INT > 0 ~ 'Expand',
                 TRUE ~ 'Stable'
               ),
               BAA_STATUS = case_when(
                 BAA_RATE < 0 & BAA_RATE + BAA_RATE_INT < 0 ~ 'Decline',
                 BAA_RATE < 0 & BAA_RATE + BAA_RATE_INT > 0 ~ 'Stable',
                 BAA_RATE > 0 & BAA_RATE - BAA_RATE_INT > 0 ~ 'Expand',
                 TRUE ~ 'Stable'
               )
               ) %>%
        mutate(SI_STATUS = case_when(
          TPA_STATUS == 'Expand' & BAA_STATUS == 'Expand' ~ 'Expand',
          TPA_STATUS == 'Expand' & BAA_STATUS == 'Stable' ~ 'Marginal Expand',
          TPA_STATUS == 'Stable' & BAA_STATUS == 'Expand' ~ 'Marginal Expand',
          TPA_STATUS == 'Stable' & BAA_STATUS == 'Stable' ~ 'Stable',
          TPA_STATUS == 'Stable' & BAA_STATUS == 'Decline' ~ 'Marginal Decline',
          TPA_STATUS == 'Decline' & BAA_STATUS == 'Stable' ~ 'Marginal Decline',
          TPA_STATUS == 'Decline' & BAA_STATUS == 'Decline' ~ 'Decline',
          TPA_STATUS != BAA_STATUS & SUST_INDEX < 0  ~ 'Marginal Decline',
          TRUE ~ 'Marginal Expand'
        ))
    })



    if (totals) {
      tOut <- tOut %>%
        select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE, TPA_RATE_INT, BAA_RATE_INT,
               TPA_STATUS, SI_STATUS, BAA_STATUS, CHNG_TPA, CHNG_BAA, PREV_TPA, PREV_BAA,
               INSECT_RATE, DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE, UNKNOWN_RATE, SILV_RATE,
               TPA_RATE_SE, BAA_RATE_SE, CHNG_TPA_SE, CHNG_BAA_SE, PREV_TPA_SE, PREV_BAA_SE,
               INSECT_RATE_SE, DISEASE_RATE_SE, FIRE_RATE_SE, ANIMAL_RATE_SE, WEATHER_RATE_SE, VEG_RATE_SE, UNKNOWN_RATE_SE, SILV_RATE_SE,
               nPlots)

    } else {
      tOut <- tOut %>%
        select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE,TPA_RATE_INT, BAA_RATE_INT,
               SI_STATUS, TPA_STATUS, BAA_STATUS,
               INSECT_RATE, DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE, UNKNOWN_RATE, SILV_RATE,
               TPA_RATE_SE, BAA_RATE_SE,
               INSECT_RATE_SE, DISEASE_RATE_SE, FIRE_RATE_SE, ANIMAL_RATE_SE, WEATHER_RATE_SE, VEG_RATE_SE, UNKNOWN_RATE_SE, SILV_RATE_SE,
               nPlots)
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

  # Return a spatial object
  if (!is.null(polys) & byPlot == FALSE) {
    ## NO IMPLICIT NA
    nospGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE])
    nospSym <- syms(nospGrp)
    tOut <- complete(tOut, !!!nospSym)
    ## If species, we don't want unique combos of variables related to same species
    ## but we do want NAs in polys where species are present
    if (length(nospGrp) < length(grpBy)){
      spGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME')])
      spSym <- syms(spGrp)
      tOut <- complete(tOut, nesting(!!!nospSym))
    }

    suppressMessages({suppressWarnings({
      tOut <- left_join(tOut, polys, by = 'polyID') %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})

    ## Makes it horrible to work with as a dataframe
    if (returnSpatial == FALSE) tOut <- select(tOut, -c(geometry))
  } else if (!is.null(polys) & byPlot){
    polys <- as.data.frame(polys)
    tOut <- left_join(tOut, select(polys, -c(geometry)), by = 'polyID')
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
  return(tOut)
}


### RUNS CLIMATE DATA AS WELL -- > symmetric percent change
si <- function(db,
               grpBy = NULL,
               polys = NULL,
               returnSpatial = FALSE,
               bySpecies = FALSE,
               bySizeClass = FALSE,
               landType = 'forest',
               treeType = 'live',
               minLive = 0,
               method = 'annual',
               lambda = .5,
               treeDomain = NULL,
               areaDomain = NULL,
               totals = TRUE,
               byPlot = FALSE,
               nCores = 1) {

  if (!('grow_drought_sev' %in% names(db$PLOT))){
    stop("Need climate data?")
  }

  ## Need a plotCN
  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[10,] %>% # Just the first row
        left_join(select(db$COND, PLT_CN, names(db$COND)[names(db$COND) %in% names(db$PLOT) == FALSE]), by = 'PLT_CN') %>%
        inner_join(select(db$TREE, PLT_CN, names(db$TREE)[names(db$TREE) %in% c(names(db$PLOT), names(db$COND)) == FALSE]), by = 'PLT_CN') %>%
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

  reqTables <- c('PLOT', 'TREE', 'TREE_GRM_COMPONENT', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method, 'unknown. Defaulting to Temporally Indifferent (TI).'))
  }

  # ## No EXP_GROW available for Western States, make sure we warn that values will be returned as 0
  # # These states do not allow temporal queries. Things are extremely weird with their eval groups
  # noGrow <- c(02,03,04,07,08,11,14,15,16, 23, 30, 32, 35,43,49, 78)
  # if(any(unique(db$PLOT$STATECD) %in% noGrow)){
  #   vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
  #   fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
  #   warning(paste('Recruitment data unavailable for: ', toString(fancyName) , '. Returning 0 for all recruitment estimates which include these states.', sep = ''))
  # }
  # # These states do not allow change estimates.
  # if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
  #   vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
  #   fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
  #   stop(paste('Growth & Mortality Estimates unavailable for: ', paste(as.character(fancyName), collapse = ', '), sep = ''))
  # }



  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', 'PLOT_STATUS_CD', grpBy)}
  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy

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
    #grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
    grpBy = c(grpBy, 'polyID')

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
        right_join(select(db$PLOT, PLT_CN, pltID), by = 'pltID') %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  }
  if (byPlot & returnSpatial){
    grpBy <- c(grpBy, 'LON', 'LAT')
  } # END AREAL

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
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

  # Same as above for tree (ex. trees > 20 ft tall)
  treeDomain <- substitute(treeDomain)
  tD <- eval(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$TREE$tD <- as.numeric(tD)


  ### Snag the EVALIDs that are needed
  db$POP_EVAL  <- db$POP_EVAL %>%
    #left_join(ga, by = 'END_INVYR') %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD', 'GROWTH_ACCT') %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP %in% c('EXPVOL')) %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)# %>%
  #group_by(END_INVYR) %>%
  #summarise(id = list(EVALID)

  ## Make an annual panel ID, associated with an INVYR

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'GROWTH_ACCT', 'END_INVYR', 'EVAL_TYP')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
    ## Join on REMPER PLOTS
    left_join(select(db$PLOT, PLT_CN, REMPER, PREV_PLT_CN, DESIGNCD, PLOT_STATUS_CD), by = 'PLT_CN') %>%
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & DESIGNCD == 1 & PLOT_STATUS_CD != 3) %>%
    mutate_if(is.factor,
              as.character)

  ### Which estimator to use?
  if (str_to_upper(method) %in% c('ANNUAL')){
    ## Want to use the year where plots are measured, no repeats
    ## Breaking this up into pre and post reporting becuase
    ## Estimation units get weird on us otherwise
    popOrig <- pops
    pops <- pops %>%
      group_by(STATECD) %>%
      filter(END_INVYR == INVYR) %>%
      ungroup()

    prePops <- popOrig %>%
      group_by(STATECD) %>%
      filter(INVYR < min(END_INVYR, na.rm = TRUE)) %>%
      distinct(PLT_CN, .keep_all = TRUE) %>%
      ungroup()

    pops <- bind_rows(pops, prePops) %>%
      mutate(YEAR = INVYR)

  } else {     # Otherwise temporally indifferent
    pops <- rename(pops, YEAR = END_INVYR)
  }

  ## P2POINTCNT column is NOT consistent for annnual estimates, plots
  ## within individual strata and est units are related to different INVYRs
  p2_INVYR <- pops %>%
    group_by(ESTN_UNIT_CN, STRATUM_CN, INVYR) %>%
    summarize(P2POINTCNT_INVYR = length(unique(PLT_CN)))
  ## Want a count of p2 points / eu, gets screwed up with grouping below
  p2eu_INVYR <- p2_INVYR %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, INVYR, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN, INVYR) %>%
    summarize(p2eu_INVYR = sum(P2POINTCNT_INVYR, na.rm = TRUE))
  p2eu <- pops %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN) %>%
    summarize(p2eu = sum(P2POINTCNT, na.rm = TRUE))

  ## Rejoin
  pops <- pops %>%
    left_join(p2_INVYR, by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'INVYR')) %>%
    left_join(p2eu_INVYR, by = c('ESTN_UNIT_CN', 'INVYR')) %>%
    left_join(p2eu, by = 'ESTN_UNIT_CN')


  ## Recode a few of the estimation methods to make things easier below
  pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')


  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
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
  db$TREE$htClass <- makeClasses(db$TREE$HT, interval = 5)


  # # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
  # if (!is.null(polys)){
  #   aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
  # } else {
  #   aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
  # }

  ## Only the necessary plots for EVAL of interest
  db$PLOT <- filter(db$PLOT, PLT_CN %in% pops$PLT_CN | PLT_CN %in% pops$PREV_PLT_CN)

  ## Reduce the memory load for others
  db <- clipFIA(db, mostRecent = FALSE)

  #### Need to scale our variables globally
  db$TREE <- db$TREE %>%
    mutate(BAA = basalArea(DIA) * TPA_UNADJ)
  # ## Scaling factors
  # tpaMean <- mean(db$TREE$TPA_UNADJ, na.rm = TRUE)
  # tpaSD <- sd(db$TREE$TPA_UNADJ, na.rm = TRUE)
  # baaMean <- mean(db$TREE$BAA, na.rm = TRUE)
  # baaSD <- sd(db$TREE$BAA, na.rm = TRUE)
  # ## Apply them
  # db$TREE <- db$TREE %>%
  #   mutate(TPA_UNADJ = scale(TPA_UNADJ),
  #          BAA = scale(BAA))


  ## Merging state and county codes
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
      out <- parLapply(cl, X = names(plts), fun = siHelper1, plts, db, grpBy, byPlot, minLive)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = siHelper1, plts, db, grpBy, byPlot, minLive, mc.cores = nCores)
    }
  })


  if (byPlot){
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tOut <- bind_rows(out[names(out) == 't'])

    ## Standardize the changes in each state variable
    tOut <- tOut %>%
      mutate(#TPA_RATE = CHNG_TPA / REMPER / abs(mean(CHNG_TPA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE)),
             #BAA_RATE = CHNG_BAA / REMPER / abs(mean(CHNG_BAA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE)),
             TPA_RATE = (scale(CHNG_TPA) +
               (mean(CHNG_TPA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE) / sd(CHNG_TPA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE))) /
               REMPER,
             BAA_RATE = (scale(CHNG_BAA) +
                           (mean(CHNG_BAA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE) / sd(CHNG_BAA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE))) /
               REMPER,

             x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
             y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
             M = sqrt(x^2 + y^2),
             SI = if_else(x < 0, -M, M)) %>%
      select(-c(x,y,M)) %>%
      select(YEAR, PLT_CN, PREV_PLT_CN, REMPER, grpBy[grpBy != 'YEAR'], SI, TPA_RATE, BAA_RATE,
             everything())

    ## Save these
    tpaRateMean <- mean(tOut$CHNG_TPA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE)
    baaRateMean <- mean(tOut$CHNG_BAA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE)
    tpaRateSD <- sd(tOut$CHNG_BAA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE)
    baaRateSD <- sd(tOut$CHNG_BAA[tOut$PLOT_STATUS_CD == 1], na.rm = TRUE)


    ## Make it spatial
    if (returnSpatial){
      tOut <- tOut %>%
        filter(!is.na(LAT) & !is.na(LON)) %>%
        st_as_sf(coords = c('LON', 'LAT'),
                 crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }
    ## Population estimation
  } else {
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    a <- bind_rows(out[names(out) == 'a'])
    t <- bind_rows(out[names(out) == 't'])
    full <- bind_rows(out[names(out) == 'full'])

    ## Standardize the changes in each state variable
    t <- t %>%
      mutate(#TPA_RATE = CHNG_TPA / REMPER / abs(mean(CHNG_TPA[t$plotIn == 1], na.rm = TRUE)),
             #BAA_RATE = CHNG_BAA / REMPER / abs(mean(CHNG_BAA[t$plotIn == 1], na.rm = TRUE)),
             TPA_RATE = (scale(CHNG_TPA) +
                           (mean(CHNG_TPA[t$plotIn == 1], na.rm = TRUE) / sd(CHNG_TPA[t$plotIn == 1], na.rm = TRUE))) /
               REMPER,
             BAA_RATE = (scale(CHNG_BAA) +
                           (mean(CHNG_BAA[t$plotIn == 1], na.rm = TRUE) / sd(CHNG_BAA[t$plotIn == 1], na.rm = TRUE))) /
               REMPER)
    ## Save these
    tpaRateMean <- mean(t$CHNG_TPA[t$plotIn == 1], na.rm = TRUE)
    baaRateMean <- mean(t$CHNG_BAA[t$plotIn == 1], na.rm = TRUE)
    tpaRateSD <- sd(t$CHNG_TPA[t$plotIn == 1], na.rm = TRUE)
    baaRateSD <- sd(t$CHNG_BAA[t$plotIn == 1], na.rm = TRUE)
    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)
    #aGrpBy <- c('YEAR', aGrpBy)


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
        out <- parLapply(cl, X = names(popState), fun = siHelper2, popState, a, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = siHelper2, popState, a, t, grpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])

    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
      ### ---- SIMPLE MOVING AVERAGE
      if (str_to_upper(method) == 'SMA'){
        ## Assuming a uniform weighting scheme
        wgts <- pops %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(wgt = 1 / length(unique(INVYR)))

        #aEst <- left_join(aEst, wgts, by = 'ESTN_UNIT_CN')
        tEst <- left_join(tEst, wgts, by = 'ESTN_UNIT_CN')

        #### ----- Linear MOVING AVERAGE
      } else if (str_to_upper(method) == 'LMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))

        ## Want a number of INVYRs per EU
        neu <- wgts %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(n = sum(rank, na.rm = TRUE))

        ## Rejoining and computing wgts
        wgts <- wgts %>%
          left_join(neu, by = 'ESTN_UNIT_CN') %>%
          mutate(wgt = rank / n) %>%
          ungroup() %>%
          select(ESTN_UNIT_CN, INVYR, wgt)

        #aEst <- left_join(aEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

        #### ----- EXPONENTIAL MOVING AVERAGE
      } else if (str_to_upper(method) == 'EMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))


        if (length(lambda) < 2){
          ## Want sum of weighitng functions
          neu <- wgts %>%
            mutate(l = lambda) %>%
            group_by(ESTN_UNIT_CN) %>%
            summarize(l = first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = 'ESTN_UNIT_CN') %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(ESTN_UNIT_CN, INVYR, wgt)
        } else {
          grpBy <- c('lambda', grpBy)
          #aGrpBy <- c('lambda', aGrpBy)
          ## Duplicate weights for each level of lambda
          yrWgts <- list()
          for (i in 1:length(unique(lambda))) {
            yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
          }
          wgts <- bind_rows(yrWgts)
          ## Want sum of weighitng functions
          neu <- wgts %>%
            group_by(lambda, ESTN_UNIT_CN) %>%
            summarize(l = first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = c('lambda', 'ESTN_UNIT_CN')) %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(lambda, ESTN_UNIT_CN, INVYR, wgt)
        }

        #aEst <- left_join(aEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

      }

      ### Applying the weights
      # Area
      # aEst <- aEst %>%
      #   mutate_at(vars(aEst), ~(.*wgt)) %>%
      #   mutate_at(vars(aVar), ~(.*(wgt^2))) %>%
      #   group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
      #   summarize_at(vars(aEst:plotIn_AREA), sum, na.rm = TRUE)


      tEst <- tEst %>%
        mutate_at(vars(ctEst:sspEst), ~(.*wgt)) %>%
        mutate_at(vars(ctVar:cvEst_ssp), ~(.*(wgt^2))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize_at(vars(ctEst:plotIn_t), sum, na.rm = TRUE)

    }

    ##---------------------  TOTALS and RATIOS
    # Tree
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    ##---------------------  TOTALS and RATIOS
    suppressWarnings({
      tOut <- tTotal %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        mutate(CHNG_TPA = ctEst,
               CHNG_BAA = cbEst,
               PREV_TPA = ptEst,
               PREV_BAA = pbEst,
               TPA_RATE = ctEst / ptEst,
               BAA_RATE = cbEst / pbEst,
               SI = siEst / faEst,
               ## Components of SI
               TPA_MORT = tmortEst / ptEst,
               TPA_RECR = trecrEst / ptEst,
               BAA_MORT = bmortEst / pbEst,
               BAA_RECR = brecrEst / pbEst,
               BAA_GROW = bgrowEst / pbEst,
               ELEV = elevEst / faEst,
               INSECT_RATE = bugEst / pbEst,
               DISEASE_RATE = diseaseEst / pbEst,
               FIRE_RATE = fireEst / pbEst,
               ANIMAL_RATE = animalEst / pbEst,
               WEATHER_RATE = weatherEst / pbEst,
               VEG_RATE = vegEst / pbEst,
               UNKNOWN_RATE = unEst / pbEst,
               SILV_RATE = silvEst / pbEst,
               MORT_RATE = mortEst / ptEst,
               DROUGHT_SEV = dEst / faEst,
               WET_SEV = wEst / faEst,
               ALL_SEV = aEst / faEst,
               GROW_DROUGHT_SEV = gdEst / faEst,
               GROW_WET_SEV = gwEst / faEst,
               GROW_ALL_SEV = gaEst / faEst,
               VPD_ANOM = vpdEst / faEst,
               TMAX_ANOM = tmaxEst / faEst,
               TMEAN_ANOM = tmeanEst / faEst,
               GROW_VPD_ANOM = gvpdEst / faEst,
               GROW_TMAX_ANOM = gtmaxEst / faEst,
               GROW_TMEAN_ANOM = gtmeanEst / faEst,
               H_a_species = hspEst / faEst,
               Eh_a_species = espEst / faEst,
               S_a_species = sspEst / faEst,
               H_a_struct = hstEst / faEst,
               Eh_a_struct = estEst / faEst,
               S_a_struct = sstEst / faEst,

               hstVar = (1/faEst^2) * (hstVar + (H_a_struct^2 * faVar) - 2 * H_a_struct * cvEst_hst),
               estVar = (1/faEst^2) * (estVar + (Eh_a_struct^2 * faVar) - 2 * Eh_a_struct * cvEst_est),
               sstVar = (1/faEst^2) * (sstVar + (S_a_struct^2 * faVar) - 2 * S_a_struct * cvEst_sst),
               hspVar = (1/faEst^2) * (hspVar + (H_a_species^2 * faVar) - 2 * H_a_species * cvEst_hsp),
               espVar = (1/faEst^2) * (espVar + (Eh_a_species^2 * faVar) - 2 * Eh_a_species * cvEst_esp),
               sspVar = (1/faEst^2) * (sspVar + (S_a_species^2 * faVar) - 2 * S_a_species * cvEst_ssp),


               ## TOTAL SE
               CHNG_TPA_SE = sqrt(ctVar) / abs(ctEst) * 100,
               CHNG_BAA_SE = sqrt(cbVar) / abs(cbEst) * 100,
               PREV_TPA_SE = sqrt(ptVar) / abs(ptEst) * 100,
               PREV_BAA_SE = sqrt(pbVar) / abs(pbEst) * 100,
               ## Ratio variance
               ctVar = (1/PREV_TPA^2) * (ctVar + (TPA_RATE^2 * ptVar) - 2 * TPA_RATE * cvEst_ct),
               cbVar = (1/PREV_BAA^2) * (cbVar + (BAA_RATE^2 * pbVar) - 2 * BAA_RATE * cvEst_cb),
               siVar = (1/faEst^2) * (siVar + (SI^2 * faVar) - 2 * SI * cvEst_si),
               tmortVar = (1/PREV_TPA^2) * (tmortVar + (TPA_MORT^2 * ptVar) - 2 * TPA_MORT * cvEst_tmort),
               trecrVar = (1/PREV_TPA^2) * (trecrVar + (TPA_RECR^2 * ptVar) - 2 * TPA_RECR * cvEst_trecr),
               bmortVar = (1/PREV_BAA^2) * (bmortVar + (BAA_MORT^2 * pbVar) - 2 * BAA_MORT * cvEst_bmort),
               brecrVar = (1/PREV_BAA^2) * (brecrVar + (BAA_RECR^2 * pbVar) - 2 * BAA_RECR * cvEst_brecr),
               bgrowVar = (1/PREV_BAA^2) * (bgrowVar + (BAA_GROW^2 * pbVar) - 2 * BAA_GROW * cvEst_bgrow),
               elevVar = (1/faEst^2) * (elevVar + (ELEV^2 * faVar) - 2 * ELEV * cvEst_elev),

               bugVar = (1/PREV_BAA^2) * (bugVar + (INSECT_RATE^2 * pbVar) - 2 * INSECT_RATE * cvEst_bug),
               diseaseVar = (1/PREV_BAA^2) * (diseaseVar + (DISEASE_RATE^2 * pbVar) - 2 * DISEASE_RATE * cvEst_disease),
               fireVar = (1/PREV_BAA^2) * (fireVar + (FIRE_RATE^2 * pbVar) - 2 * FIRE_RATE * cvEst_fire),
               animalVar = (1/PREV_BAA^2) * (animalVar + (ANIMAL_RATE^2 * pbVar) - 2 * ANIMAL_RATE * cvEst_animal),
               weatherVar = (1/PREV_BAA^2) * (weatherVar + (WEATHER_RATE^2 * pbVar) - 2 * WEATHER_RATE * cvEst_weather),
               vegVar = (1/PREV_BAA^2) * (vegVar + (VEG_RATE^2 * pbVar) - 2 * VEG_RATE * cvEst_veg),
               unVar = (1/PREV_BAA^2) * (unVar + (UNKNOWN_RATE^2 * pbVar) - 2 * UNKNOWN_RATE * cvEst_un),
               silvVar = (1/PREV_BAA^2) * (silvVar + (SILV_RATE^2 * pbVar) - 2 * SILV_RATE * cvEst_silv),
               mortVar = (1/PREV_TPA^2) * (mortVar + (MORT_RATE^2 * ptVar) - 2 * MORT_RATE * cvEst_mort),
               dVar = (1/faEst^2) * (dVar + (DROUGHT_SEV^2 * faVar) - 2 * DROUGHT_SEV * cvEst_d),
               wVar = (1/faEst^2) * (wVar + (WET_SEV^2 * faVar) - 2 * WET_SEV * cvEst_w),
               aVar = (1/faEst^2) * (aVar + (ALL_SEV^2 * faVar) - 2 * ALL_SEV * cvEst_a),
               gdVar = (1/faEst^2) * (gdVar + (GROW_DROUGHT_SEV^2 * faVar) - 2 * GROW_DROUGHT_SEV * cvEst_gd),
               gwVar = (1/faEst^2) * (gwVar + (GROW_WET_SEV^2 * faVar) - 2 * GROW_WET_SEV * cvEst_gw),
               gaVar = (1/faEst^2) * (gaVar + (GROW_ALL_SEV^2 * faVar) - 2 * GROW_ALL_SEV * cvEst_ga),
               vpdVar = (1/faEst^2) * (vpdVar + (VPD_ANOM^2 * faVar) - 2 * VPD_ANOM * cvEst_vpd),
               tmaxVar = (1/faEst^2) * (tmaxVar + (TMAX_ANOM^2 * faVar) - 2 * TMAX_ANOM * cvEst_tmax),
               tmeanVar = (1/faEst^2) * (tmeanVar + (TMEAN_ANOM^2 * faVar) - 2 * TMEAN_ANOM * cvEst_tmean),
               gvpdVar = (1/faEst^2) * (gvpdVar + (GROW_VPD_ANOM^2 * faVar) - 2 * GROW_VPD_ANOM * cvEst_gvpd),
               gtmaxVar = (1/faEst^2) * (gtmaxVar + (GROW_TMAX_ANOM^2 * faVar) - 2 * GROW_TMAX_ANOM * cvEst_gtmax),
               gtmeanVar = (1/faEst^2) * (gtmeanVar + (GROW_TMEAN_ANOM^2 * faVar) - 2 * GROW_TMEAN_ANOM * cvEst_gtmean),

               H_a_species_SE = sqrt(hspVar) / abs(H_a_species) * 100,
               Eh_a_species_SE = sqrt(espVar) / abs(Eh_a_species) * 100,
               S_a_species_SE = sqrt(sspVar) / abs(S_a_species) * 100,
               H_a_struct_SE = sqrt(hstVar) / abs(H_a_struct) * 100,
               Eh_a_struct_SE = sqrt(estVar) / abs(Eh_a_struct) * 100,
               S_a_struct_SE = sqrt(sstVar) / abs(S_a_struct) * 100,
               ## RATIO SE
               TPA_RATE_SE = sqrt(ctVar) / abs(TPA_RATE) * 100,
               BAA_RATE_SE = sqrt(cbVar) / abs(BAA_RATE) * 100,
               SI_SE = sqrt(siVar) / abs(SI) * 100,
               TPA_MORT_SE = sqrt(tmortVar) / abs(TPA_MORT) * 100,
               TPA_RECR_SE = sqrt(trecrVar) / abs(TPA_RECR) * 100,
               BAA_MORT_SE = sqrt(bmortVar) / abs(BAA_MORT) * 100,
               BAA_RECR_SE = sqrt(brecrVar) / abs(BAA_RECR) * 100,
               BAA_GROW_SE = sqrt(bgrowVar) / abs(BAA_GROW) * 100,
               ELEV_SE = sqrt(elevVar) / abs(ELEV) * 100,
               INSECT_RATE_SE = sqrt(bugVar) / abs(INSECT_RATE) * 100,
               DISEASE_RATE_SE = sqrt(diseaseVar) / abs(DISEASE_RATE) * 100,
               FIRE_RATE_SE = sqrt(fireVar) / abs(FIRE_RATE) * 100,
               ANIMAL_RATE_SE = sqrt(animalVar) / abs(ANIMAL_RATE) * 100,
               WEATHER_RATE_SE = sqrt(weatherVar) / abs(WEATHER_RATE) * 100,
               VEG_RATE_SE = sqrt(vegVar) / abs(VEG_RATE) * 100,
               UNKNOWN_RATE_SE = sqrt(unVar) / abs(UNKNOWN_RATE) * 100,
               SILV_RATE_SE = sqrt(silvVar) / abs(SILV_RATE) * 100,
               MORT_RATE_SE = sqrt(mortVar) / abs(MORT_RATE) * 100,
               DROUGHT_SEV_SE = sqrt(dVar) / abs(DROUGHT_SEV) * 100,
               WET_SEV_SE = sqrt(wVar) / abs(WET_SEV) * 100,
               ALL_SEV_SE = sqrt(aVar) / abs(ALL_SEV) * 100,
               GROW_DROUGHT_SEV_SE = sqrt(gdVar) / abs(GROW_DROUGHT_SEV) * 100,
               GROW_WET_SEV_SE = sqrt(gwVar) / abs(GROW_WET_SEV) * 100,
               GROW_ALL_SEV_SE = sqrt(gaVar) / abs(GROW_WET_SEV) * 100,
               VPD_ANOM_SE = sqrt(vpdVar) / abs(VPD_ANOM) * 100,
               TMAX_ANOM_SE = sqrt(tmaxVar) / abs(TMAX_ANOM) * 100,
               TMEAN_ANOM_SE = sqrt(tmeanVar) / abs(TMEAN_ANOM) * 100,
               GROW_VPD_ANOM_SE = sqrt(gvpdVar) / abs(GROW_VPD_ANOM) * 100,
               GROW_TMAX_ANOM_SE = sqrt(gtmaxVar) / abs(GROW_TMAX_ANOM) * 100,
               GROW_TMEAN_ANOM_SE = sqrt(gtmeanVar) / abs(GROW_TMEAN_ANOM) * 100,
               #SI_SE = sqrt(Mvar) / abs(SI) * 100,
               nPlots = plotIn_t,
               nTotal = nh,
               TPA_RATE_INT = abs(TPA_RATE) * 1.96 * TPA_RATE_SE / 100,
               BAA_RATE_INT = abs(BAA_RATE) * 1.96 * BAA_RATE_SE / 100,
               SI_INT = abs(SI) * 1.96 * SI_SE / 100,
               TPA_STATUS = case_when(
                 TPA_RATE < 0 & TPA_RATE + TPA_RATE_INT < 0 ~ 'Decline',
                 TPA_RATE < 0 & TPA_RATE + TPA_RATE_INT > 0 ~ 'Stable',
                 TPA_RATE > 0 & TPA_RATE - TPA_RATE_INT > 0 ~ 'Expand',
                 TRUE ~ 'Stable'
               ),
               BAA_STATUS = case_when(
                 BAA_RATE < 0 & BAA_RATE + BAA_RATE_INT < 0 ~ 'Decline',
                 BAA_RATE < 0 & BAA_RATE + BAA_RATE_INT > 0 ~ 'Stable',
                 BAA_RATE > 0 & BAA_RATE - BAA_RATE_INT > 0 ~ 'Expand',
                 TRUE ~ 'Stable'
               )
        ) %>%
        mutate(SI_STATUS = case_when(
          SI < 0 & SI + SI_INT < 0 ~ 'Decline',
          SI < 0 & SI + SI_INT > 0 ~ 'Stable',
          SI > 0 & SI - SI_INT > 0  ~ 'Expand',
          TRUE ~ 'Stable'
        ))
    })

    ### REALLY DON'T LIKE WORKING WITH FULL TABLES, WORK ON A FIX
    ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
    fullArea <- full %>%
      left_join(select(pops, PLT_CN, YEAR, INVYR), by = 'PLT_CN') %>%
      group_by(.dots = grpBy) %>%
      summarize(H_g_struct = divIndex(htClass, -BAA * tDI, index = 'H'),
                Eh_g_struct = divIndex(htClass, -BAA * tDI, index = 'Eh'),
                S_g_struct = divIndex(htClass, -BAA * tDI, index = 'S'),
                H_g_species = divIndex(SPCD, -TPA_UNADJ * tDI, index = 'H'),
                Eh_g_species = divIndex(SPCD, -TPA_UNADJ * tDI, index = 'Eh'),
                S_g_species = divIndex(SPCD, -TPA_UNADJ * tDI, index = 'S'))

    tOut <- left_join(tOut, fullArea, by = grpBy) %>%
      mutate(H_b_struct = H_g_struct - H_a_struct,
             Eh_b_struct = Eh_g_struct - Eh_a_struct,
             S_b_struct = S_g_struct - S_a_struct,
             H_b_species = H_g_species - H_a_species,
             Eh_b_species = Eh_g_species - Eh_a_species,
             S_b_species = S_g_species - S_a_species)



    if (totals) {
      tOut <- tOut %>%
        select(grpBy, SI, TPA_RATE, BAA_RATE, SI_STATUS, TPA_STATUS, BAA_STATUS,
               TPA_MORT, TPA_RECR, BAA_MORT, BAA_RECR, BAA_GROW, ELEV,
               SI_INT, TPA_RATE_INT, BAA_RATE_INT, CHNG_TPA, CHNG_BAA, PREV_TPA, PREV_BAA,
               INSECT_RATE, DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE, UNKNOWN_RATE, SILV_RATE,
               DROUGHT_SEV, WET_SEV, ALL_SEV, GROW_DROUGHT_SEV, GROW_WET_SEV, GROW_ALL_SEV,
               VPD_ANOM, TMAX_ANOM, TMEAN_ANOM, GROW_VPD_ANOM, GROW_TMAX_ANOM, GROW_TMEAN_ANOM,
               H_a_struct, H_b_struct, H_g_struct, S_a_struct, S_b_struct, S_g_struct, Eh_a_struct, Eh_b_struct, Eh_g_struct,
               H_a_species, H_b_species, H_g_species, S_a_species, S_b_species, S_g_species, Eh_a_species, Eh_b_species, Eh_g_species,
               SI_SE, TPA_RATE_SE, BAA_RATE_SE,
               TPA_MORT_SE, TPA_RECR_SE, BAA_MORT_SE, BAA_RECR_SE, BAA_GROW_SE, ELEV_SE,
               CHNG_TPA_SE, CHNG_BAA_SE, PREV_TPA_SE, PREV_BAA_SE,
               INSECT_RATE_SE, DISEASE_RATE_SE, FIRE_RATE_SE, ANIMAL_RATE_SE, WEATHER_RATE_SE, VEG_RATE_SE, UNKNOWN_RATE_SE, SILV_RATE_SE,
               DROUGHT_SEV_SE, WET_SEV_SE, ALL_SEV_SE, GROW_DROUGHT_SEV_SE, GROW_WET_SEV_SE, GROW_ALL_SEV_SE,
               VPD_ANOM_SE, TMAX_ANOM_SE, TMEAN_ANOM_SE, GROW_VPD_ANOM_SE, GROW_TMAX_ANOM_SE, GROW_TMEAN_ANOM_SE,
               H_a_struct_SE, Eh_a_struct_SE, S_a_struct_SE, H_a_species_SE, Eh_a_species_SE, S_a_species_SE,

               nPlots)

    } else {
      tOut <- tOut %>%
        select(grpBy, SI, TPA_RATE, BAA_RATE, SI_STATUS, TPA_STATUS, BAA_STATUS,
               TPA_MORT, TPA_RECR, BAA_MORT, BAA_RECR, BAA_GROW, ELEV,
               INSECT_RATE, DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE, UNKNOWN_RATE, SILV_RATE,
               DROUGHT_SEV, WET_SEV, ALL_SEV, GROW_DROUGHT_SEV, GROW_WET_SEV, GROW_ALL_SEV,
               VPD_ANOM, TMAX_ANOM, TMEAN_ANOM, GROW_VPD_ANOM, GROW_TMAX_ANOM, GROW_TMEAN_ANOM,
               H_a_struct, H_b_struct, H_g_struct, S_a_struct, S_b_struct, S_g_struct, Eh_a_struct, Eh_b_struct, Eh_g_struct,
               H_a_species, H_b_species, H_g_species, S_a_species, S_b_species, S_g_species, Eh_a_species, Eh_b_species, Eh_g_species,
               SI_SE, TPA_RATE_SE, BAA_RATE_SE,
               TPA_MORT_SE, TPA_RECR_SE, BAA_MORT_SE, BAA_RECR_SE, BAA_GROW_SE, ELEV_SE,
               INSECT_RATE_SE, DISEASE_RATE_SE, FIRE_RATE_SE, ANIMAL_RATE_SE, WEATHER_RATE_SE, VEG_RATE_SE, UNKNOWN_RATE_SE, SILV_RATE_SE,
               DROUGHT_SEV_SE, WET_SEV_SE, ALL_SEV_SE, GROW_DROUGHT_SEV_SE, GROW_WET_SEV_SE, GROW_ALL_SEV_SE,
               VPD_ANOM_SE, TMAX_ANOM_SE, TMEAN_ANOM_SE, GROW_VPD_ANOM_SE, GROW_TMAX_ANOM_SE, GROW_TMEAN_ANOM_SE,
               H_a_struct_SE, Eh_a_struct_SE, S_a_struct_SE, H_a_species_SE, Eh_a_species_SE, S_a_species_SE,
               nPlots)
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

  # Return a spatial object
  if (!is.null(polys) & byPlot == FALSE) {
    ## NO IMPLICIT NA
    nospGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE])
    nospSym <- syms(nospGrp)
    tOut <- complete(tOut, !!!nospSym)
    ## If species, we don't want unique combos of variables related to same species
    ## but we do want NAs in polys where species are present
    if (length(nospGrp) < length(grpBy)){
      spGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME')])
      spSym <- syms(spGrp)
      tOut <- complete(tOut, nesting(!!!nospSym))
    }

    suppressMessages({suppressWarnings({
      tOut <- left_join(tOut, polys, by = 'polyID') %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})

    ## Makes it horrible to work with as a dataframe
    if (returnSpatial == FALSE) tOut <- select(tOut, -c(geometry))
  } else if (!is.null(polys) & byPlot){
    polys <- as.data.frame(polys)
    tOut <- left_join(tOut, select(polys, -c(geometry)), by = 'polyID')
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)

  ## Standardization factors
  tOut$tpaRateMean <- tpaRateMean
  tOut$baaRateMean <- baaRateMean
  tOut$tpaRateSD <- tpaRateSD
  tOut$baaRateSD <- baaRateSD

  return(tOut)
}





### RUNS CLIMATE DATA AS WELL -- > symmetric percent change
si_backup <- function(db,
               grpBy = NULL,
               polys = NULL,
               returnSpatial = FALSE,
               bySpecies = FALSE,
               bySizeClass = FALSE,
               landType = 'forest',
               treeType = 'live',
               minLive = 0,
               method = 'annual',
               lambda = .5,
               treeDomain = NULL,
               areaDomain = NULL,
               totals = TRUE,
               byPlot = FALSE,
               nCores = 1) {

  if (!('grow_drought_sev' %in% names(db$PLOT))){
    stop("Need climate data?")
  }

  ## Need a plotCN
  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[10,] %>% # Just the first row
        left_join(select(db$COND, PLT_CN, names(db$COND)[names(db$COND) %in% names(db$PLOT) == FALSE]), by = 'PLT_CN') %>%
        inner_join(select(db$TREE, PLT_CN, names(db$TREE)[names(db$TREE) %in% c(names(db$PLOT), names(db$COND)) == FALSE]), by = 'PLT_CN') %>%
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

  reqTables <- c('PLOT', 'TREE', 'TREE_GRM_COMPONENT', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method, 'unknown. Defaulting to Temporally Indifferent (TI).'))
  }

  # ## No EXP_GROW available for Western States, make sure we warn that values will be returned as 0
  # # These states do not allow temporal queries. Things are extremely weird with their eval groups
  # noGrow <- c(02,03,04,07,08,11,14,15,16, 23, 30, 32, 35,43,49, 78)
  # if(any(unique(db$PLOT$STATECD) %in% noGrow)){
  #   vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
  #   fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
  #   warning(paste('Recruitment data unavailable for: ', toString(fancyName) , '. Returning 0 for all recruitment estimates which include these states.', sep = ''))
  # }
  # # These states do not allow change estimates.
  # if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
  #   vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
  #   fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
  #   stop(paste('Growth & Mortality Estimates unavailable for: ', paste(as.character(fancyName), collapse = ', '), sep = ''))
  # }



  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}
  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy

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
    #grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
    grpBy = c(grpBy, 'polyID')

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
        right_join(select(db$PLOT, PLT_CN, pltID), by = 'pltID') %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  }
  if (byPlot & returnSpatial){
    grpBy <- c(grpBy, 'LON', 'LAT')
  } # END AREAL

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
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

  # Same as above for tree (ex. trees > 20 ft tall)
  treeDomain <- substitute(treeDomain)
  tD <- eval(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$TREE$tD <- as.numeric(tD)


  ### Snag the EVALIDs that are needed
  db$POP_EVAL  <- db$POP_EVAL %>%
    #left_join(ga, by = 'END_INVYR') %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD', 'GROWTH_ACCT') %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP %in% c('EXPVOL')) %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)# %>%
  #group_by(END_INVYR) %>%
  #summarise(id = list(EVALID)

  ## Make an annual panel ID, associated with an INVYR

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'GROWTH_ACCT', 'END_INVYR', 'EVAL_TYP')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
    ## Join on REMPER PLOTS
    left_join(select(db$PLOT, PLT_CN, REMPER, PREV_PLT_CN, DESIGNCD, PLOT_STATUS_CD), by = 'PLT_CN') %>%
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & DESIGNCD == 1 & PLOT_STATUS_CD != 3) %>%
    mutate_if(is.factor,
              as.character)

  ### Which estimator to use?
  if (str_to_upper(method) %in% c('ANNUAL')){
    ## Want to use the year where plots are measured, no repeats
    ## Breaking this up into pre and post reporting becuase
    ## Estimation units get weird on us otherwise
    popOrig <- pops
    pops <- pops %>%
      group_by(STATECD) %>%
      filter(END_INVYR == INVYR) %>%
      ungroup()

    prePops <- popOrig %>%
      group_by(STATECD) %>%
      filter(INVYR < min(END_INVYR, na.rm = TRUE)) %>%
      distinct(PLT_CN, .keep_all = TRUE) %>%
      ungroup()

    pops <- bind_rows(pops, prePops) %>%
      mutate(YEAR = INVYR)

  } else {     # Otherwise temporally indifferent
    pops <- rename(pops, YEAR = END_INVYR)
  }

  ## P2POINTCNT column is NOT consistent for annnual estimates, plots
  ## within individual strata and est units are related to different INVYRs
  p2_INVYR <- pops %>%
    group_by(ESTN_UNIT_CN, STRATUM_CN, INVYR) %>%
    summarize(P2POINTCNT_INVYR = length(unique(PLT_CN)))
  ## Want a count of p2 points / eu, gets screwed up with grouping below
  p2eu_INVYR <- p2_INVYR %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, INVYR, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN, INVYR) %>%
    summarize(p2eu_INVYR = sum(P2POINTCNT_INVYR, na.rm = TRUE))
  p2eu <- pops %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN) %>%
    summarize(p2eu = sum(P2POINTCNT, na.rm = TRUE))

  ## Rejoin
  pops <- pops %>%
    left_join(p2_INVYR, by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'INVYR')) %>%
    left_join(p2eu_INVYR, by = c('ESTN_UNIT_CN', 'INVYR')) %>%
    left_join(p2eu, by = 'ESTN_UNIT_CN')


  ## Recode a few of the estimation methods to make things easier below
  pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')


  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
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

  ## No tree groups go into are computations
  aGrpBy = grpBy[grpBy %in% c('SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE]

  ## Make some height classes
  db$TREE$htClass <- makeClasses(db$TREE$HT, interval = 5, numLabs = TRUE)


  # # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
  # if (!is.null(polys)){
  #   aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
  # } else {
  #   aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
  # }

  ## Only the necessary plots for EVAL of interest
  db$PLOT <- filter(db$PLOT, PLT_CN %in% pops$PLT_CN | PLT_CN %in% pops$PREV_PLT_CN)

  ## Reduce the memory load for others
  db <- clipFIA(db, mostRecent = FALSE)

  ## Merging state and county codes
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
      out <- parLapply(cl, X = names(plts), fun = siHelper1, plts, db, grpBy, byPlot, minLive, aGrpBy)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = siHelper1, plts, db, grpBy, byPlot, minLive, aGrpBy, mc.cores = nCores)
    }
  })


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
    ## Population estimation
  } else {
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    a <- bind_rows(out[names(out) == 'a'])
    t <- bind_rows(out[names(out) == 't'])
    full <- bind_rows(out[names(out) == 'full'])


    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)
    #aGrpBy <- c('YEAR', aGrpBy)


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
        out <- parLapply(cl, X = names(popState), fun = siHelper2, popState, a, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = siHelper2, popState, a, t, grpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])

    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
      ### ---- SIMPLE MOVING AVERAGE
      if (str_to_upper(method) == 'SMA'){
        ## Assuming a uniform weighting scheme
        wgts <- pops %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(wgt = 1 / length(unique(INVYR)))

        #aEst <- left_join(aEst, wgts, by = 'ESTN_UNIT_CN')
        tEst <- left_join(tEst, wgts, by = 'ESTN_UNIT_CN')

        #### ----- Linear MOVING AVERAGE
      } else if (str_to_upper(method) == 'LMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))

        ## Want a number of INVYRs per EU
        neu <- wgts %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(n = sum(rank, na.rm = TRUE))

        ## Rejoining and computing wgts
        wgts <- wgts %>%
          left_join(neu, by = 'ESTN_UNIT_CN') %>%
          mutate(wgt = rank / n) %>%
          ungroup() %>%
          select(ESTN_UNIT_CN, INVYR, wgt)

        #aEst <- left_join(aEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

        #### ----- EXPONENTIAL MOVING AVERAGE
      } else if (str_to_upper(method) == 'EMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))


        if (length(lambda) < 2){
          ## Want sum of weighitng functions
          neu <- wgts %>%
            mutate(l = lambda) %>%
            group_by(ESTN_UNIT_CN) %>%
            summarize(l = first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = 'ESTN_UNIT_CN') %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(ESTN_UNIT_CN, INVYR, wgt)
        } else {
          grpBy <- c('lambda', grpBy)
          #aGrpBy <- c('lambda', aGrpBy)
          ## Duplicate weights for each level of lambda
          yrWgts <- list()
          for (i in 1:length(unique(lambda))) {
            yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
          }
          wgts <- bind_rows(yrWgts)
          ## Want sum of weighitng functions
          neu <- wgts %>%
            group_by(lambda, ESTN_UNIT_CN) %>%
            summarize(l = first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = c('lambda', 'ESTN_UNIT_CN')) %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(lambda, ESTN_UNIT_CN, INVYR, wgt)
        }

        #aEst <- left_join(aEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

      }

      ### Applying the weights
      # Area
      # aEst <- aEst %>%
      #   mutate_at(vars(aEst), ~(.*wgt)) %>%
      #   mutate_at(vars(aVar), ~(.*(wgt^2))) %>%
      #   group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
      #   summarize_at(vars(aEst:plotIn_AREA), sum, na.rm = TRUE)


      tEst <- tEst %>%
        mutate_at(vars(ctEst:sspEst), ~(.*wgt)) %>%
        mutate_at(vars(ctVar:cvEst_ssp), ~(.*(wgt^2))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize_at(vars(ctEst:plotIn_t), sum, na.rm = TRUE)

    }

    ##---------------------  TOTALS and RATIOS
    # Tree
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    ##---------------------  TOTALS and RATIOS
    suppressWarnings({
      tOut <- tTotal %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        mutate(CHNG_TPA = ctEst,
               CHNG_BAA = cbEst,
               PREV_TPA = ptEst,
               PREV_BAA = pbEst,
               TPA_RATE = ctEst / ptEst,
               BAA_RATE = cbEst / pbEst,
               SI = siEst / faEst,
               ## Components of SI
               TPA_MORT = tmortEst / ptEst,
               TPA_RECR = trecrEst / ptEst,
               BAA_MORT = bmortEst / pbEst,
               BAA_RECR = brecrEst / pbEst,
               BAA_GROW = bgrowEst / pbEst,
               ELEV = elevEst / faEst,
               INSECT_RATE = bugEst / pbEst,
               DISEASE_RATE = diseaseEst / pbEst,
               FIRE_RATE = fireEst / pbEst,
               ANIMAL_RATE = animalEst / pbEst,
               WEATHER_RATE = weatherEst / pbEst,
               VEG_RATE = vegEst / pbEst,
               UNKNOWN_RATE = unEst / pbEst,
               SILV_RATE = silvEst / pbEst,
               MORT_RATE = mortEst / ptEst,
               DROUGHT_SEV = dEst / faEst,
               WET_SEV = wEst / faEst,
               ALL_SEV = aEst / faEst,
               GROW_DROUGHT_SEV = gdEst / faEst,
               GROW_WET_SEV = gwEst / faEst,
               GROW_ALL_SEV = gaEst / faEst,
               VPD_ANOM = vpdEst / faEst,
               TMAX_ANOM = tmaxEst / faEst,
               TMEAN_ANOM = tmeanEst / faEst,
               GROW_VPD_ANOM = gvpdEst / faEst,
               GROW_TMAX_ANOM = gtmaxEst / faEst,
               GROW_TMEAN_ANOM = gtmeanEst / faEst,
               H_a_species = hspEst / faEst,
               Eh_a_species = espEst / faEst,
               S_a_species = sspEst / faEst,
               H_a_struct = hstEst / faEst,
               Eh_a_struct = estEst / faEst,
               S_a_struct = sstEst / faEst,
               STOCKING = stkEst / faEst,

               ## TOTAL SE
               CHNG_TPA_SE = sqrt(ctVar) / abs(ctEst) * 100,
               CHNG_BAA_SE = sqrt(cbVar) / abs(cbEst) * 100,
               PREV_TPA_SE = sqrt(ptVar) / abs(ptEst) * 100,
               PREV_BAA_SE = sqrt(pbVar) / abs(pbEst) * 100,
               ## Ratio variance
               ctVar = (1/PREV_TPA^2) * (ctVar + (TPA_RATE^2 * ptVar) - 2 * TPA_RATE * cvEst_ct),
               cbVar = (1/PREV_BAA^2) * (cbVar + (BAA_RATE^2 * pbVar) - 2 * BAA_RATE * cvEst_cb),
               siVar = (1/faEst^2) * (siVar + (SI^2 * faVar) - 2 * SI * cvEst_si),
               tmortVar = (1/PREV_TPA^2) * (tmortVar + (TPA_MORT^2 * ptVar) - 2 * TPA_MORT * cvEst_tmort),
               trecrVar = (1/PREV_TPA^2) * (trecrVar + (TPA_RECR^2 * ptVar) - 2 * TPA_RECR * cvEst_trecr),
               bmortVar = (1/PREV_BAA^2) * (bmortVar + (BAA_MORT^2 * pbVar) - 2 * BAA_MORT * cvEst_bmort),
               brecrVar = (1/PREV_BAA^2) * (brecrVar + (BAA_RECR^2 * pbVar) - 2 * BAA_RECR * cvEst_brecr),
               bgrowVar = (1/PREV_BAA^2) * (bgrowVar + (BAA_GROW^2 * pbVar) - 2 * BAA_GROW * cvEst_bgrow),
               elevVar = (1/faEst^2) * (elevVar + (ELEV^2 * faVar) - 2 * ELEV * cvEst_elev),

               bugVar = (1/PREV_BAA^2) * (bugVar + (INSECT_RATE^2 * pbVar) - 2 * INSECT_RATE * cvEst_bug),
               diseaseVar = (1/PREV_BAA^2) * (diseaseVar + (DISEASE_RATE^2 * pbVar) - 2 * DISEASE_RATE * cvEst_disease),
               fireVar = (1/PREV_BAA^2) * (fireVar + (FIRE_RATE^2 * pbVar) - 2 * FIRE_RATE * cvEst_fire),
               animalVar = (1/PREV_BAA^2) * (animalVar + (ANIMAL_RATE^2 * pbVar) - 2 * ANIMAL_RATE * cvEst_animal),
               weatherVar = (1/PREV_BAA^2) * (weatherVar + (WEATHER_RATE^2 * pbVar) - 2 * WEATHER_RATE * cvEst_weather),
               vegVar = (1/PREV_BAA^2) * (vegVar + (VEG_RATE^2 * pbVar) - 2 * VEG_RATE * cvEst_veg),
               unVar = (1/PREV_BAA^2) * (unVar + (UNKNOWN_RATE^2 * pbVar) - 2 * UNKNOWN_RATE * cvEst_un),
               silvVar = (1/PREV_BAA^2) * (silvVar + (SILV_RATE^2 * pbVar) - 2 * SILV_RATE * cvEst_silv),
               mortVar = (1/PREV_TPA^2) * (mortVar + (MORT_RATE^2 * ptVar) - 2 * MORT_RATE * cvEst_mort),
               dVar = (1/faEst^2) * (dVar + (DROUGHT_SEV^2 * faVar) - 2 * DROUGHT_SEV * cvEst_d),
               wVar = (1/faEst^2) * (wVar + (WET_SEV^2 * faVar) - 2 * WET_SEV * cvEst_w),
               aVar = (1/faEst^2) * (aVar + (ALL_SEV^2 * faVar) - 2 * ALL_SEV * cvEst_a),
               gdVar = (1/faEst^2) * (gdVar + (GROW_DROUGHT_SEV^2 * faVar) - 2 * GROW_DROUGHT_SEV * cvEst_gd),
               gwVar = (1/faEst^2) * (gwVar + (GROW_WET_SEV^2 * faVar) - 2 * GROW_WET_SEV * cvEst_gw),
               gaVar = (1/faEst^2) * (gaVar + (GROW_ALL_SEV^2 * faVar) - 2 * GROW_ALL_SEV * cvEst_ga),
               vpdVar = (1/faEst^2) * (vpdVar + (VPD_ANOM^2 * faVar) - 2 * VPD_ANOM * cvEst_vpd),
               tmaxVar = (1/faEst^2) * (tmaxVar + (TMAX_ANOM^2 * faVar) - 2 * TMAX_ANOM * cvEst_tmax),
               tmeanVar = (1/faEst^2) * (tmeanVar + (TMEAN_ANOM^2 * faVar) - 2 * TMEAN_ANOM * cvEst_tmean),
               gvpdVar = (1/faEst^2) * (gvpdVar + (GROW_VPD_ANOM^2 * faVar) - 2 * GROW_VPD_ANOM * cvEst_gvpd),
               gtmaxVar = (1/faEst^2) * (gtmaxVar + (GROW_TMAX_ANOM^2 * faVar) - 2 * GROW_TMAX_ANOM * cvEst_gtmax),
               gtmeanVar = (1/faEst^2) * (gtmeanVar + (GROW_TMEAN_ANOM^2 * faVar) - 2 * GROW_TMEAN_ANOM * cvEst_gtmean),

               hstVar = (1/faEst^2) * (hstVar + (H_a_struct^2 * faVar) - 2 * H_a_struct * cvEst_hst),
               estVar = (1/faEst^2) * (estVar + (Eh_a_struct^2 * faVar) - 2 * Eh_a_struct * cvEst_est),
               sstVar = (1/faEst^2) * (sstVar + (S_a_struct^2 * faVar) - 2 * S_a_struct * cvEst_sst),
               hspVar = (1/faEst^2) * (hspVar + (H_a_species^2 * faVar) - 2 * H_a_species * cvEst_hsp),
               espVar = (1/faEst^2) * (espVar + (Eh_a_species^2 * faVar) - 2 * Eh_a_species * cvEst_esp),
               sspVar = (1/faEst^2) * (sspVar + (S_a_species^2 * faVar) - 2 * S_a_species * cvEst_ssp),
               stkVar = (1/faEst^2) * (stkVar + (STOCKING^2 * faVar) - 2 * STOCKING * cvEst_stk),

               ## RATIO SE
               TPA_RATE_SE = sqrt(ctVar) / abs(TPA_RATE) * 100,
               BAA_RATE_SE = sqrt(cbVar) / abs(BAA_RATE) * 100,
               SI_SE = sqrt(siVar) / abs(SI) * 100,
               TPA_MORT_SE = sqrt(tmortVar) / abs(TPA_MORT) * 100,
               TPA_RECR_SE = sqrt(trecrVar) / abs(TPA_RECR) * 100,
               BAA_MORT_SE = sqrt(bmortVar) / abs(BAA_MORT) * 100,
               BAA_RECR_SE = sqrt(brecrVar) / abs(BAA_RECR) * 100,
               BAA_GROW_SE = sqrt(bgrowVar) / abs(BAA_GROW) * 100,
               ELEV_SE = sqrt(elevVar) / abs(ELEV) * 100,
               INSECT_RATE_SE = sqrt(bugVar) / abs(INSECT_RATE) * 100,
               DISEASE_RATE_SE = sqrt(diseaseVar) / abs(DISEASE_RATE) * 100,
               FIRE_RATE_SE = sqrt(fireVar) / abs(FIRE_RATE) * 100,
               ANIMAL_RATE_SE = sqrt(animalVar) / abs(ANIMAL_RATE) * 100,
               WEATHER_RATE_SE = sqrt(weatherVar) / abs(WEATHER_RATE) * 100,
               VEG_RATE_SE = sqrt(vegVar) / abs(VEG_RATE) * 100,
               UNKNOWN_RATE_SE = sqrt(unVar) / abs(UNKNOWN_RATE) * 100,
               SILV_RATE_SE = sqrt(silvVar) / abs(SILV_RATE) * 100,
               MORT_RATE_SE = sqrt(mortVar) / abs(MORT_RATE) * 100,
               DROUGHT_SEV_SE = sqrt(dVar) / abs(DROUGHT_SEV) * 100,
               WET_SEV_SE = sqrt(wVar) / abs(WET_SEV) * 100,
               ALL_SEV_SE = sqrt(aVar) / abs(ALL_SEV) * 100,
               GROW_DROUGHT_SEV_SE = sqrt(gdVar) / abs(GROW_DROUGHT_SEV) * 100,
               GROW_WET_SEV_SE = sqrt(gwVar) / abs(GROW_WET_SEV) * 100,
               GROW_ALL_SEV_SE = sqrt(gaVar) / abs(GROW_WET_SEV) * 100,
               VPD_ANOM_SE = sqrt(vpdVar) / abs(VPD_ANOM) * 100,
               TMAX_ANOM_SE = sqrt(tmaxVar) / abs(TMAX_ANOM) * 100,
               TMEAN_ANOM_SE = sqrt(tmeanVar) / abs(TMEAN_ANOM) * 100,
               GROW_VPD_ANOM_SE = sqrt(gvpdVar) / abs(GROW_VPD_ANOM) * 100,
               GROW_TMAX_ANOM_SE = sqrt(gtmaxVar) / abs(GROW_TMAX_ANOM) * 100,
               GROW_TMEAN_ANOM_SE = sqrt(gtmeanVar) / abs(GROW_TMEAN_ANOM) * 100,

               H_a_species_SE = sqrt(hspVar) / abs(H_a_species) * 100,
               Eh_a_species_SE = sqrt(espVar) / abs(Eh_a_species) * 100,
               S_a_species_SE = sqrt(sspVar) / abs(S_a_species) * 100,
               H_a_struct_SE = sqrt(hstVar) / abs(H_a_struct) * 100,
               Eh_a_struct_SE = sqrt(estVar) / abs(Eh_a_struct) * 100,
               S_a_struct_SE = sqrt(sstVar) / abs(S_a_struct) * 100,
               STOCKING_SE = sqrt(stkVar) / abs(STOCKING) * 100,

               #SI_SE = sqrt(Mvar) / abs(SI) * 100,
               nPlots = plotIn_t,
               nTotal = nh,
               TPA_RATE_INT = abs(TPA_RATE) * 1.96 * TPA_RATE_SE / 100,
               BAA_RATE_INT = abs(BAA_RATE) * 1.96 * BAA_RATE_SE / 100,
               SI_INT = abs(SI) * 1.96 * SI_SE / 100,
               TPA_STATUS = case_when(
                 TPA_RATE < 0 & TPA_RATE + TPA_RATE_INT < 0 ~ 'Decline',
                 TPA_RATE < 0 & TPA_RATE + TPA_RATE_INT > 0 ~ 'Stable',
                 TPA_RATE > 0 & TPA_RATE - TPA_RATE_INT > 0 ~ 'Expand',
                 TRUE ~ 'Stable'
               ),
               BAA_STATUS = case_when(
                 BAA_RATE < 0 & BAA_RATE + BAA_RATE_INT < 0 ~ 'Decline',
                 BAA_RATE < 0 & BAA_RATE + BAA_RATE_INT > 0 ~ 'Stable',
                 BAA_RATE > 0 & BAA_RATE - BAA_RATE_INT > 0 ~ 'Expand',
                 TRUE ~ 'Stable'
               )
        ) %>%
        mutate(SI_STATUS = case_when(
          SI < 0 & SI + SI_INT < 0 ~ 'Decline',
          SI < 0 & SI + SI_INT > 0 ~ 'Stable',
          SI > 0 & SI - SI_INT > 0 ~ 'Expand',
          TRUE ~ 'Stable'
        ))
    })

    if (bySpecies == FALSE){
      ### REALLY DON'T LIKE WORKING WITH FULL TABLES, WORK ON A FIX
      ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
      fullArea <- full %>%
        left_join(select(pops, PLT_CN, YEAR), by = 'PLT_CN') %>%
        group_by(.dots = grpBy) %>%
        summarize(H_g_struct = divIndex(htClass, -BAA * tDI, index = 'H'),
                  Eh_g_struct = divIndex(htClass, -BAA * tDI, index = 'Eh'),
                  S_g_struct = divIndex(htClass, -BAA * tDI, index = 'S'),
                  H_g_species = divIndex(SPCD, -TPA_UNADJ * tDI, index = 'H'),
                  Eh_g_species = divIndex(SPCD, -TPA_UNADJ * tDI, index = 'Eh'),
                  S_g_species = divIndex(SPCD, -TPA_UNADJ * tDI, index = 'S'))

      tOut <- left_join(tOut, fullArea, by = grpBy) %>%
        mutate(H_b_struct = H_g_struct - H_a_struct,
               Eh_b_struct = Eh_g_struct - Eh_a_struct,
               S_b_struct = S_g_struct - S_a_struct,
               H_b_species = H_g_species - H_a_species,
               Eh_b_species = Eh_g_species - Eh_a_species,
               S_b_species = S_g_species - S_a_species)
    } else {
      tOut <- tOut %>%
        mutate(H_g_struct = NA,
               Eh_g_struct = NA,
               S_g_struct = NA,
               H_g_species = NA,
               Eh_g_species = NA,
               S_g_species = NA,
               H_b_struct = NA,
               Eh_b_struct = NA,
               S_b_struct = NA,
               H_b_species = NA,
               Eh_b_species = NA,
               S_b_species = NA)
    }




    if (totals) {
      tOut <- tOut %>%
        select(grpBy, SI, TPA_RATE, BAA_RATE, SI_STATUS, TPA_STATUS, BAA_STATUS,
               TPA_MORT, TPA_RECR, BAA_MORT, BAA_RECR, BAA_GROW, ELEV,
               SI_INT, TPA_RATE_INT, BAA_RATE_INT, CHNG_TPA, CHNG_BAA, PREV_TPA, PREV_BAA,
               INSECT_RATE, DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE, UNKNOWN_RATE, SILV_RATE,
               DROUGHT_SEV, WET_SEV, ALL_SEV, GROW_DROUGHT_SEV, GROW_WET_SEV, GROW_ALL_SEV,
               VPD_ANOM, TMAX_ANOM, TMEAN_ANOM, GROW_VPD_ANOM, GROW_TMAX_ANOM, GROW_TMEAN_ANOM,
               STOCKING, H_a_struct, H_b_struct, H_g_struct, S_a_struct, S_b_struct, S_g_struct, Eh_a_struct, Eh_b_struct, Eh_g_struct,
               H_a_species, H_b_species, H_g_species, S_a_species, S_b_species, S_g_species, Eh_a_species, Eh_b_species, Eh_g_species,
               SI_SE, TPA_RATE_SE, BAA_RATE_SE,
               TPA_MORT_SE, TPA_RECR_SE, BAA_MORT_SE, BAA_RECR_SE, BAA_GROW_SE, ELEV_SE,
               CHNG_TPA_SE, CHNG_BAA_SE, PREV_TPA_SE, PREV_BAA_SE,
               INSECT_RATE_SE, DISEASE_RATE_SE, FIRE_RATE_SE, ANIMAL_RATE_SE, WEATHER_RATE_SE, VEG_RATE_SE, UNKNOWN_RATE_SE, SILV_RATE_SE,
               DROUGHT_SEV_SE, WET_SEV_SE, ALL_SEV_SE, GROW_DROUGHT_SEV_SE, GROW_WET_SEV_SE, GROW_ALL_SEV_SE,
               VPD_ANOM_SE, TMAX_ANOM_SE, TMEAN_ANOM_SE, GROW_VPD_ANOM_SE, GROW_TMAX_ANOM_SE, GROW_TMEAN_ANOM_SE,
               STOCKING_SE, H_a_struct_SE, Eh_a_struct_SE, S_a_struct_SE, H_a_species_SE, Eh_a_species_SE, S_a_species_SE,
               nPlots)

    } else {
      tOut <- tOut %>%
        select(grpBy, SI, TPA_RATE, BAA_RATE, SI_STATUS, TPA_STATUS, BAA_STATUS,
               TPA_MORT, TPA_RECR, BAA_MORT, BAA_RECR, BAA_GROW, ELEV,
               INSECT_RATE, DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE, UNKNOWN_RATE, SILV_RATE,
               DROUGHT_SEV, WET_SEV, ALL_SEV, GROW_DROUGHT_SEV, GROW_WET_SEV, GROW_ALL_SEV,
               VPD_ANOM, TMAX_ANOM, TMEAN_ANOM, GROW_VPD_ANOM, GROW_TMAX_ANOM, GROW_TMEAN_ANOM,
               STOCKING, H_a_struct, H_b_struct, H_g_struct, S_a_struct, S_b_struct, S_g_struct, Eh_a_struct, Eh_b_struct, Eh_g_struct,
               H_a_species, H_b_species, H_g_species, S_a_species, S_b_species, S_g_species, Eh_a_species, Eh_b_species, Eh_g_species,
               SI_SE, TPA_RATE_SE, BAA_RATE_SE,
               TPA_MORT_SE, TPA_RECR_SE, BAA_MORT_SE, BAA_RECR_SE, BAA_GROW_SE, ELEV_SE,
               INSECT_RATE_SE, DISEASE_RATE_SE, FIRE_RATE_SE, ANIMAL_RATE_SE, WEATHER_RATE_SE, VEG_RATE_SE, UNKNOWN_RATE_SE, SILV_RATE_SE,
               DROUGHT_SEV_SE, WET_SEV_SE, ALL_SEV_SE, GROW_DROUGHT_SEV_SE, GROW_WET_SEV_SE, GROW_ALL_SEV_SE,
               VPD_ANOM_SE, TMAX_ANOM_SE, TMEAN_ANOM_SE, GROW_VPD_ANOM_SE, GROW_TMAX_ANOM_SE, GROW_TMEAN_ANOM_SE,
               STOCKING_SE, H_a_struct_SE, Eh_a_struct_SE, S_a_struct_SE, H_a_species_SE, Eh_a_species_SE, S_a_species_SE,
               nPlots)
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

  # Return a spatial object
  if (!is.null(polys) & byPlot == FALSE) {
    ## NO IMPLICIT NA
    nospGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE])
    nospSym <- syms(nospGrp)
    tOut <- complete(tOut, !!!nospSym)
    ## If species, we don't want unique combos of variables related to same species
    ## but we do want NAs in polys where species are present
    if (length(nospGrp) < length(grpBy)){
      spGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME')])
      spSym <- syms(spGrp)
      tOut <- complete(tOut, nesting(!!!nospSym))
    }

    suppressMessages({suppressWarnings({
      tOut <- left_join(tOut, polys, by = 'polyID') %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})

    ## Makes it horrible to work with as a dataframe
    if (returnSpatial == FALSE) tOut <- select(tOut, -c(geometry))
  } else if (!is.null(polys) & byPlot){
    polys <- as.data.frame(polys)
    tOut <- left_join(tOut, select(polys, -c(geometry)), by = 'polyID')
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
  return(tOut)
}

### RUNS CLIMATE DATA AS WELL -- Percent change -- asymmetrical
si_old <- function(db,
               grpBy = NULL,
               polys = NULL,
               returnSpatial = FALSE,
               bySpecies = FALSE,
               bySizeClass = FALSE,
               landType = 'forest',
               treeType = 'live',
               minLive = 0,
               method = 'annual',
               lambda = .5,
               treeDomain = NULL,
               areaDomain = NULL,
               totals = TRUE,
               byPlot = FALSE,
               nCores = 1) {

  if (!('grow_drought_sev' %in% names(db$PLOT))){
    stop("Need climate data?")
  }

  ## Need a plotCN
  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[10,] %>% # Just the first row
        left_join(select(db$COND, PLT_CN, names(db$COND)[names(db$COND) %in% names(db$PLOT) == FALSE]), by = 'PLT_CN') %>%
        inner_join(select(db$TREE, PLT_CN, names(db$TREE)[names(db$TREE) %in% c(names(db$PLOT), names(db$COND)) == FALSE]), by = 'PLT_CN') %>%
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

  reqTables <- c('PLOT', 'TREE', 'TREE_GRM_COMPONENT', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method, 'unknown. Defaulting to Temporally Indifferent (TI).'))
  }

  # ## No EXP_GROW available for Western States, make sure we warn that values will be returned as 0
  # # These states do not allow temporal queries. Things are extremely weird with their eval groups
  # noGrow <- c(02,03,04,07,08,11,14,15,16, 23, 30, 32, 35,43,49, 78)
  # if(any(unique(db$PLOT$STATECD) %in% noGrow)){
  #   vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
  #   fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
  #   warning(paste('Recruitment data unavailable for: ', toString(fancyName) , '. Returning 0 for all recruitment estimates which include these states.', sep = ''))
  # }
  # # These states do not allow change estimates.
  # if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
  #   vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
  #   fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
  #   stop(paste('Growth & Mortality Estimates unavailable for: ', paste(as.character(fancyName), collapse = ', '), sep = ''))
  # }



  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}
  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy

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
    #grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
    grpBy = c(grpBy, 'polyID')

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
        right_join(select(db$PLOT, PLT_CN, pltID), by = 'pltID') %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  }
  if (byPlot & returnSpatial){
    grpBy <- c(grpBy, 'LON', 'LAT')
  } # END AREAL

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
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

  # Same as above for tree (ex. trees > 20 ft tall)
  treeDomain <- substitute(treeDomain)
  tD <- eval(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$TREE$tD <- as.numeric(tD)


  ### Snag the EVALIDs that are needed
  db$POP_EVAL  <- db$POP_EVAL %>%
    #left_join(ga, by = 'END_INVYR') %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD', 'GROWTH_ACCT') %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP %in% c('EXPVOL')) %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)# %>%
  #group_by(END_INVYR) %>%
  #summarise(id = list(EVALID)

  ## Make an annual panel ID, associated with an INVYR

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'GROWTH_ACCT', 'END_INVYR', 'EVAL_TYP')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
    ## Join on REMPER PLOTS
    left_join(select(db$PLOT, PLT_CN, REMPER, PREV_PLT_CN, DESIGNCD), by = 'PLT_CN') %>%
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & DESIGNCD == 1) %>%
    mutate_if(is.factor,
              as.character)

  ### Which estimator to use?
  if (str_to_upper(method) %in% c('ANNUAL')){
    ## Want to use the year where plots are measured, no repeats
    ## Breaking this up into pre and post reporting becuase
    ## Estimation units get weird on us otherwise
    popOrig <- pops
    pops <- pops %>%
      group_by(STATECD) %>%
      filter(END_INVYR == INVYR) %>%
      ungroup()

    prePops <- popOrig %>%
      group_by(STATECD) %>%
      filter(INVYR < min(END_INVYR, na.rm = TRUE)) %>%
      distinct(PLT_CN, .keep_all = TRUE) %>%
      ungroup()

    pops <- bind_rows(pops, prePops) %>%
      mutate(YEAR = INVYR)

  } else {     # Otherwise temporally indifferent
    pops <- rename(pops, YEAR = END_INVYR)
  }

  ## P2POINTCNT column is NOT consistent for annnual estimates, plots
  ## within individual strata and est units are related to different INVYRs
  p2_INVYR <- pops %>%
    group_by(ESTN_UNIT_CN, STRATUM_CN, INVYR) %>%
    summarize(P2POINTCNT_INVYR = length(unique(PLT_CN)))
  ## Want a count of p2 points / eu, gets screwed up with grouping below
  p2eu_INVYR <- p2_INVYR %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, INVYR, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN, INVYR) %>%
    summarize(p2eu_INVYR = sum(P2POINTCNT_INVYR, na.rm = TRUE))
  p2eu <- pops %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN) %>%
    summarize(p2eu = sum(P2POINTCNT, na.rm = TRUE))

  ## Rejoin
  pops <- pops %>%
    left_join(p2_INVYR, by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'INVYR')) %>%
    left_join(p2eu_INVYR, by = c('ESTN_UNIT_CN', 'INVYR')) %>%
    left_join(p2eu, by = 'ESTN_UNIT_CN')


  ## Recode a few of the estimation methods to make things easier below
  pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')


  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
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


  ## Merging state and county codes
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
      out <- parLapply(cl, X = names(plts), fun = siHelper1_old, plts, db, grpBy, byPlot, minLive)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = siHelper1_old, plts, db, grpBy, byPlot, minLive, mc.cores = nCores)
    }
  })


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
    ## Population estimation
  } else {
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    a <- bind_rows(out[names(out) == 'a'])
    t <- bind_rows(out[names(out) == 't'])


    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)
    #aGrpBy <- c('YEAR', aGrpBy)


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
        out <- parLapply(cl, X = names(popState), fun = siHelper2_old, popState, a, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = siHelper2_old, popState, a, t, grpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])

    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
      ### ---- SIMPLE MOVING AVERAGE
      if (str_to_upper(method) == 'SMA'){
        ## Assuming a uniform weighting scheme
        wgts <- pops %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(wgt = 1 / length(unique(INVYR)))

        #aEst <- left_join(aEst, wgts, by = 'ESTN_UNIT_CN')
        tEst <- left_join(tEst, wgts, by = 'ESTN_UNIT_CN')

        #### ----- Linear MOVING AVERAGE
      } else if (str_to_upper(method) == 'LMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))

        ## Want a number of INVYRs per EU
        neu <- wgts %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(n = sum(rank, na.rm = TRUE))

        ## Rejoining and computing wgts
        wgts <- wgts %>%
          left_join(neu, by = 'ESTN_UNIT_CN') %>%
          mutate(wgt = rank / n) %>%
          ungroup() %>%
          select(ESTN_UNIT_CN, INVYR, wgt)

        #aEst <- left_join(aEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

        #### ----- EXPONENTIAL MOVING AVERAGE
      } else if (str_to_upper(method) == 'EMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))


        if (length(lambda) < 2){
          ## Want sum of weighitng functions
          neu <- wgts %>%
            mutate(l = lambda) %>%
            group_by(ESTN_UNIT_CN) %>%
            summarize(l = first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = 'ESTN_UNIT_CN') %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(ESTN_UNIT_CN, INVYR, wgt)
        } else {
          grpBy <- c('lambda', grpBy)
          #aGrpBy <- c('lambda', aGrpBy)
          ## Duplicate weights for each level of lambda
          yrWgts <- list()
          for (i in 1:length(unique(lambda))) {
            yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
          }
          wgts <- bind_rows(yrWgts)
          ## Want sum of weighitng functions
          neu <- wgts %>%
            group_by(lambda, ESTN_UNIT_CN) %>%
            summarize(l = first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = c('lambda', 'ESTN_UNIT_CN')) %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(lambda, ESTN_UNIT_CN, INVYR, wgt)
        }

        #aEst <- left_join(aEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

      }

      ### Applying the weights
      # Area
      # aEst <- aEst %>%
      #   mutate_at(vars(aEst), ~(.*wgt)) %>%
      #   mutate_at(vars(aVar), ~(.*(wgt^2))) %>%
      #   group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
      #   summarize_at(vars(aEst:plotIn_AREA), sum, na.rm = TRUE)


      tEst <- tEst %>%
        mutate_at(vars(ctEst:gaEst), ~(.*wgt)) %>%
        mutate_at(vars(ctVar:cvEst_ga), ~(.*(wgt^2))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize_at(vars(ctEst:plotIn_t), sum, na.rm = TRUE)

    }

    ##---------------------  TOTALS and RATIOS
    # Tree
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    ##---------------------  TOTALS and RATIOS
    suppressWarnings({
      tOut <- tTotal %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        mutate(CHNG_TPA = ctEst,
               CHNG_BAA = cbEst,
               PREV_TPA = ptEst,
               PREV_BAA = pbEst,
               TPA_RATE = ctEst / ptEst,
               BAA_RATE = cbEst / pbEst,
               INSECT_RATE = bugEst / pbEst,
               DISEASE_RATE = diseaseEst / pbEst,
               FIRE_RATE = fireEst / pbEst,
               ANIMAL_RATE = animalEst / pbEst,
               WEATHER_RATE = weatherEst / pbEst,
               VEG_RATE = vegEst / pbEst,
               UNKNOWN_RATE = unEst / pbEst,
               SILV_RATE = silvEst / pbEst,
               MORT_RATE = mortEst / ptEst,
               DROUGHT_SEV = dEst / faEst,
               WET_SEV = wEst / faEst,
               ALL_SEV = aEst / faEst,
               GROW_DROUGHT_SEV = gdEst / faEst,
               GROW_WET_SEV = gwEst / faEst,
               GROW_ALL_SEV = gaEst / faEst,
               VPD_ANOM = vpdEst / faEst,
               TMAX_ANOM = tmaxEst / faEst,
               TMEAN_ANOM = tmeanEst / faEst,
               GROW_VPD_ANOM = gvpdEst / faEst,
               GROW_TMAX_ANOM = gtmaxEst / faEst,
               GROW_TMEAN_ANOM = gtmeanEst / faEst,

               x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
               y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
               M = sqrt(x^2 + y^2),
               SUST_INDEX = if_else(x < 0, -M, M),
               ## TOTAL SE
               CHNG_TPA_SE = sqrt(ctVar) / abs(ctEst) * 100,
               CHNG_BAA_SE = sqrt(cbVar) / abs(cbEst) * 100,
               PREV_TPA_SE = sqrt(ptVar) / abs(ptEst) * 100,
               PREV_BAA_SE = sqrt(pbVar) / abs(pbEst) * 100,
               ## Ratio variance
               ctVar = (1/PREV_TPA^2) * (ctVar + (TPA_RATE^2 * ptVar) - 2 * TPA_RATE * cvEst_ct),
               cbVar = (1/PREV_BAA^2) * (cbVar + (BAA_RATE^2 * pbVar) - 2 * BAA_RATE * cvEst_cb),
               bugVar = (1/PREV_BAA^2) * (bugVar + (INSECT_RATE^2 * pbVar) - 2 * INSECT_RATE * cvEst_bug),
               diseaseVar = (1/PREV_BAA^2) * (diseaseVar + (DISEASE_RATE^2 * pbVar) - 2 * DISEASE_RATE * cvEst_disease),
               fireVar = (1/PREV_BAA^2) * (fireVar + (FIRE_RATE^2 * pbVar) - 2 * FIRE_RATE * cvEst_fire),
               animalVar = (1/PREV_BAA^2) * (animalVar + (ANIMAL_RATE^2 * pbVar) - 2 * ANIMAL_RATE * cvEst_animal),
               weatherVar = (1/PREV_BAA^2) * (weatherVar + (WEATHER_RATE^2 * pbVar) - 2 * WEATHER_RATE * cvEst_weather),
               vegVar = (1/PREV_BAA^2) * (vegVar + (VEG_RATE^2 * pbVar) - 2 * VEG_RATE * cvEst_veg),
               unVar = (1/PREV_BAA^2) * (unVar + (UNKNOWN_RATE^2 * pbVar) - 2 * UNKNOWN_RATE * cvEst_un),
               silvVar = (1/PREV_BAA^2) * (silvVar + (SILV_RATE^2 * pbVar) - 2 * SILV_RATE * cvEst_silv),
               mortVar = (1/PREV_TPA^2) * (mortVar + (MORT_RATE^2 * ptVar) - 2 * MORT_RATE * cvEst_mort),
               dVar = (1/faEst^2) * (dVar + (DROUGHT_SEV^2 * faVar) - 2 * DROUGHT_SEV * cvEst_d),
               wVar = (1/faEst^2) * (wVar + (WET_SEV^2 * faVar) - 2 * WET_SEV * cvEst_w),
               aVar = (1/faEst^2) * (aVar + (ALL_SEV^2 * faVar) - 2 * ALL_SEV * cvEst_a),
               gdVar = (1/faEst^2) * (gdVar + (GROW_DROUGHT_SEV^2 * faVar) - 2 * GROW_DROUGHT_SEV * cvEst_gd),
               gwVar = (1/faEst^2) * (gwVar + (GROW_WET_SEV^2 * faVar) - 2 * GROW_WET_SEV * cvEst_gw),
               gaVar = (1/faEst^2) * (gaVar + (GROW_ALL_SEV^2 * faVar) - 2 * GROW_ALL_SEV * cvEst_ga),
               vpdVar = (1/faEst^2) * (vpdVar + (VPD_ANOM^2 * faVar) - 2 * VPD_ANOM * cvEst_vpd),
               tmaxVar = (1/faEst^2) * (tmaxVar + (TMAX_ANOM^2 * faVar) - 2 * TMAX_ANOM * cvEst_tmax),
               tmeanVar = (1/faEst^2) * (tmeanVar + (TMEAN_ANOM^2 * faVar) - 2 * TMEAN_ANOM * cvEst_tmean),
               gvpdVar = (1/faEst^2) * (gvpdVar + (GROW_VPD_ANOM^2 * faVar) - 2 * GROW_VPD_ANOM * cvEst_gvpd),
               gtmaxVar = (1/faEst^2) * (gtmaxVar + (GROW_TMAX_ANOM^2 * faVar) - 2 * GROW_TMAX_ANOM * cvEst_gtmax),
               gtmeanVar = (1/faEst^2) * (gtmeanVar + (GROW_TMEAN_ANOM^2 * faVar) - 2 * GROW_TMEAN_ANOM * cvEst_gtmean),
               ## RATIO SE
               TPA_RATE_SE = sqrt(ctVar) / abs(TPA_RATE) * 100,
               BAA_RATE_SE = sqrt(cbVar) / abs(BAA_RATE) * 100,
               INSECT_RATE_SE = sqrt(bugVar) / abs(INSECT_RATE) * 100,
               DISEASE_RATE_SE = sqrt(diseaseVar) / abs(DISEASE_RATE) * 100,
               FIRE_RATE_SE = sqrt(fireVar) / abs(FIRE_RATE) * 100,
               ANIMAL_RATE_SE = sqrt(animalVar) / abs(ANIMAL_RATE) * 100,
               WEATHER_RATE_SE = sqrt(weatherVar) / abs(WEATHER_RATE) * 100,
               VEG_RATE_SE = sqrt(vegVar) / abs(VEG_RATE) * 100,
               UNKNOWN_RATE_SE = sqrt(unVar) / abs(UNKNOWN_RATE) * 100,
               SILV_RATE_SE = sqrt(silvVar) / abs(SILV_RATE) * 100,
               MORT_RATE_SE = sqrt(mortVar) / abs(MORT_RATE) * 100,
               DROUGHT_SEV_SE = sqrt(dVar) / abs(DROUGHT_SEV) * 100,
               WET_SEV_SE = sqrt(wVar) / abs(WET_SEV) * 100,
               ALL_SEV_SE = sqrt(aVar) / abs(ALL_SEV) * 100,
               GROW_DROUGHT_SEV_SE = sqrt(gdVar) / abs(GROW_DROUGHT_SEV) * 100,
               GROW_WET_SEV_SE = sqrt(gwVar) / abs(GROW_WET_SEV) * 100,
               GROW_ALL_SEV_SE = sqrt(gaVar) / abs(GROW_WET_SEV) * 100,
               VPD_ANOM_SE = sqrt(vpdVar) / abs(VPD_ANOM) * 100,
               TMAX_ANOM_SE = sqrt(tmaxVar) / abs(TMAX_ANOM) * 100,
               TMEAN_ANOM_SE = sqrt(tmeanVar) / abs(TMEAN_ANOM) * 100,
               GROW_VPD_ANOM_SE = sqrt(gvpdVar) / abs(GROW_VPD_ANOM) * 100,
               GROW_TMAX_ANOM_SE = sqrt(gtmaxVar) / abs(GROW_TMAX_ANOM) * 100,
               GROW_TMEAN_ANOM_SE = sqrt(gtmeanVar) / abs(GROW_TMEAN_ANOM) * 100,
               #SUST_INDEX_SE = sqrt(Mvar) / abs(SUST_INDEX) * 100,
               nPlots = plotIn_t,
               nTotal = nh,
               TPA_RATE_INT = abs(TPA_RATE) * 1.96 * TPA_RATE_SE / 100,
               BAA_RATE_INT = abs(BAA_RATE) * 1.96 * BAA_RATE_SE / 100,
               TPA_STATUS = case_when(
                 TPA_RATE < 0 & TPA_RATE + TPA_RATE_INT < 0 ~ 'Decline',
                 TPA_RATE < 0 & TPA_RATE + TPA_RATE_INT > 0 ~ 'Stable',
                 TPA_RATE > 0 & TPA_RATE - TPA_RATE_INT > 0 ~ 'Expand',
                 TRUE ~ 'Stable'
               ),
               BAA_STATUS = case_when(
                 BAA_RATE < 0 & BAA_RATE + BAA_RATE_INT < 0 ~ 'Decline',
                 BAA_RATE < 0 & BAA_RATE + BAA_RATE_INT > 0 ~ 'Stable',
                 BAA_RATE > 0 & BAA_RATE - BAA_RATE_INT > 0 ~ 'Expand',
                 TRUE ~ 'Stable'
               )
        ) %>%
        mutate(SI_STATUS = case_when(
          TPA_STATUS == 'Expand' & BAA_STATUS == 'Expand' ~ 'Expand',
          TPA_STATUS == 'Expand' & BAA_STATUS == 'Stable' ~ 'Marginal Expand',
          TPA_STATUS == 'Stable' & BAA_STATUS == 'Expand' ~ 'Marginal Expand',
          TPA_STATUS == 'Stable' & BAA_STATUS == 'Stable' ~ 'Stable',
          TPA_STATUS == 'Stable' & BAA_STATUS == 'Decline' ~ 'Marginal Decline',
          TPA_STATUS == 'Decline' & BAA_STATUS == 'Stable' ~ 'Marginal Decline',
          TPA_STATUS == 'Decline' & BAA_STATUS == 'Decline' ~ 'Decline',
          TPA_STATUS != BAA_STATUS & SUST_INDEX < 0  ~ 'Marginal Decline',
          TRUE ~ 'Marginal Expand'
        ))
    })



    if (totals) {
      tOut <- tOut %>%
        select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE, MORT_RATE, TPA_RATE_INT, BAA_RATE_INT,
               TPA_STATUS, SI_STATUS, BAA_STATUS, CHNG_TPA, CHNG_BAA, PREV_TPA, PREV_BAA,
               INSECT_RATE, DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE, UNKNOWN_RATE, SILV_RATE,
               DROUGHT_SEV, WET_SEV, ALL_SEV, GROW_DROUGHT_SEV, GROW_WET_SEV, GROW_ALL_SEV,
               VPD_ANOM, TMAX_ANOM, TMEAN_ANOM, GROW_VPD_ANOM, GROW_TMAX_ANOM, GROW_TMEAN_ANOM,
               TPA_RATE_SE, BAA_RATE_SE, MORT_RATE_SE, CHNG_TPA_SE, CHNG_BAA_SE, PREV_TPA_SE, PREV_BAA_SE,
               INSECT_RATE_SE, DISEASE_RATE_SE, FIRE_RATE_SE, ANIMAL_RATE_SE, WEATHER_RATE_SE, VEG_RATE_SE, UNKNOWN_RATE_SE, SILV_RATE_SE,
               DROUGHT_SEV_SE, WET_SEV_SE, ALL_SEV_SE, GROW_DROUGHT_SEV_SE, GROW_WET_SEV_SE, GROW_ALL_SEV_SE,
               VPD_ANOM_SE, TMAX_ANOM_SE, TMEAN_ANOM_SE, GROW_VPD_ANOM_SE, GROW_TMAX_ANOM_SE, GROW_TMEAN_ANOM_SE,

               nPlots)

    } else {
      tOut <- tOut %>%
        select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE, MORT_RATE, TPA_RATE_INT, BAA_RATE_INT,
               SI_STATUS, TPA_STATUS, BAA_STATUS,
               INSECT_RATE, DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE, UNKNOWN_RATE, SILV_RATE,
               DROUGHT_SEV, WET_SEV, ALL_SEV, GROW_DROUGHT_SEV, GROW_WET_SEV, GROW_ALL_SEV,
               VPD_ANOM, TMAX_ANOM, TMEAN_ANOM, GROW_VPD_ANOM, GROW_TMAX_ANOM, GROW_TMEAN_ANOM,
               TPA_RATE_SE, BAA_RATE_SE, MORT_RATE_SE,
               INSECT_RATE_SE, DISEASE_RATE_SE, FIRE_RATE_SE, ANIMAL_RATE_SE, WEATHER_RATE_SE, VEG_RATE_SE, UNKNOWN_RATE_SE, SILV_RATE_SE,
               DROUGHT_SEV_SE, WET_SEV_SE, ALL_SEV_SE, GROW_DROUGHT_SEV_SE, GROW_WET_SEV_SE, GROW_ALL_SEV_SE,
               VPD_ANOM_SE, TMAX_ANOM_SE, TMEAN_ANOM_SE, GROW_VPD_ANOM_SE, GROW_TMAX_ANOM_SE, GROW_TMEAN_ANOM_SE,
               nPlots)
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

  # Return a spatial object
  if (!is.null(polys) & byPlot == FALSE) {
    ## NO IMPLICIT NA
    nospGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE])
    nospSym <- syms(nospGrp)
    tOut <- complete(tOut, !!!nospSym)
    ## If species, we don't want unique combos of variables related to same species
    ## but we do want NAs in polys where species are present
    if (length(nospGrp) < length(grpBy)){
      spGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME')])
      spSym <- syms(spGrp)
      tOut <- complete(tOut, nesting(!!!nospSym))
    }

    suppressMessages({suppressWarnings({
      tOut <- left_join(tOut, polys, by = 'polyID') %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})

    ## Makes it horrible to work with as a dataframe
    if (returnSpatial == FALSE) tOut <- select(tOut, -c(geometry))
  } else if (!is.null(polys) & byPlot){
    polys <- as.data.frame(polys)
    tOut <- left_join(tOut, select(polys, -c(geometry)), by = 'polyID')
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)


  ## Add our scaling factors
  tOut$tpaMean <- tpaMean
  tOut$tpaSD <- tpaSD
  tOut$baaMean <- baaMean
  tOut$baaSD <- baaSD

  return(tOut)
}
