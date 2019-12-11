biomassNew <- function(db,
                       grpBy = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       bySpecies = FALSE,
                       bySizeClass = FALSE,
                       landType = 'forest',
                       treeType = 'live',
                       method = 'TI',
                       yrs = 5,
                       treeDomain = NULL,
                       areaDomain = NULL,
                       totals = FALSE,
                       byPlot = FALSE,
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
        inner_join(db$TREE, by = 'PLT_CN') %>%
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
  }
  # Tree Type domain indicator
  if (tolower(treeType) == 'live'){
    db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
  } else if (tolower(treeType) == 'dead'){
    db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 2 & db$TREE$STANDING_DEAD_CD == 1, 1, 0)
  } else if (tolower(treeType) == 'gs'){
    db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1 & db$TREE$DIA >= 5 & db$TREE$TREECLCD == 2, 1, 0)
  } else if (tolower(treeType) == 'all'){
    db$TREE$typeD <- 1
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
  db$POP_EVAL<- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
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

  # ### The population tables
  # pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'END_INVYR')) %>%
  #   rename(EVAL_CN = CN) %>%
  #   left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
  #   rename(ESTN_UNIT_CN = CN) %>%
  #   left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
  #   rename(STRATUM_CN = CN) %>%
  #   left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
  #   rename(YEAR = END_INVYR) %>%
  #   mutate_if(is.factor,
  #             as.character)
  #
  # ### Which estimator to use?
  # if (str_to_upper(method) %in% c('ANNUAL', "SMA")){
  #   popOrig <- pops
  #   ## Want to use the year where plots are measured, no repeats
  #   pops <- filter(pops, YEAR == INVYR)
  # }

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
    db$TREE$sizeClass <- makeClasses(db$TREE$DIA, interval = 2)
    db$TREE <- db$TREE[!is.na(db$TREE$sizeClass),]
  }


  # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
  if (!is.null(polys)){
    aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
  } else {
    aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
  }

  ## Only the necessary plots for EVAL of interest
  db$PLOT <- filter(db$PLOT, PLT_CN %in% pops$PLT_CN)

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
      out <- parLapply(cl, X = names(plts), fun = tpaHelper1, plts, db, grpBy, aGrpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = tpaHelper1, plts, db, grpBy, aGrpBy, byPlot, mc.cores = nCores)
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

    }
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
        ## Use the same cluster as above
        # cl <- makeCluster(nCores)
        # clusterEvalQ(cl, {
        #   library(dplyr)
        #   library(stringr)
        #   library(rFIA)
        # })
        out <- parLapply(cl, X = names(popState), fun = tpaHelper2, popState, a, t, grpBy, aGrpBy)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = tpaHelper2, popState, a, t, grpBy, aGrpBy, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    aEst <- bind_rows(out[names(out) == 'aEst'])
    tEst <- bind_rows(out[names(out) == 'tEst'])


    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA')){
      ## Need a STATECD on aEst and tEst to join wgts
      aEst <- left_join(aEst, select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN'))
      tEst <- left_join(tEst, select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN'))

      #### Summarizing to state level here to apply weights by panel
      # Area
      aEst <- aEst %>%
        group_by(STATECD, .dots = aGrpBy) %>%
        summarize(aEst = sum(aEst, na.rm = TRUE),
                  aVar = sum(aVar, na.rm = TRUE),
                  plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))
      # Tree
      tEst <- tEst %>%
        group_by(STATECD, .dots = grpBy) %>%
        #left_join(aTotal, by = c(aGrpBy)) %>%
        summarize(tEst = sum(tEst, na.rm = TRUE),
                  bEst = sum(bEst, na.rm = TRUE),
                  tTEst = sum(tTEst, na.rm = TRUE), ## Need to sum this
                  bTEst = sum(bTEst, na.rm = TRUE), ## Need to sum this
                  tTVar = sum(tTVar, na.rm = TRUE),
                  bTVar = sum(bTVar, na.rm = TRUE),
                  cvEst_tT = sum(cvEst_tT, na.rm = TRUE),
                  cvEst_bT = sum(cvEst_bT, na.rm = TRUE),
                  ## Variances
                  tVar = sum(tVar, na.rm = TRUE),
                  bVar = sum(bVar, na.rm = TRUE),
                  #aVar = first(aVar),
                  cvEst_t = sum(cvEst_t, na.rm = TRUE),
                  cvEst_b = sum(cvEst_b, na.rm = TRUE),
                  plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

      ### ---- SIMPLE MOVING AVERAGE
      if (method == 'SMA'){
        ## Assuming a uniform weighting scheme
        popOrig <- mutate(popOrig, YEAR = END_INVYR)
        wgts <- popOrig %>%
          group_by(YEAR, STATECD) %>%
          summarize(wgt = 1 / length(unique(INVYR))) %>%
          ## Expand it out again
          right_join(popOrig, by = c('YEAR', 'STATECD')) %>%
          distinct(YEAR, INVYR, STATECD, .keep_all = TRUE) %>%
          select(YEAR, INVYR, STATECD, wgt)

        # Area
        aEst <- left_join(wgts, aEst, by = c('INVYR' = 'YEAR', 'STATECD')) %>%
          group_by(STATECD, .dots = aGrpBy) %>%
          summarize(aEst = sum(aEst * wgt, na.rm = TRUE),
                    aVar = sum(aVar * wgt^2, na.rm = TRUE),
                    plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

        tEst <- left_join(wgts, tEst, by = c('INVYR' = 'YEAR', 'STATECD')) %>%
          group_by(STATECD, .dots = grpBy) %>%
          #left_join(aTotal, by = c(aGrpBy)) %>%
          summarize(tEst = sum(tEst * wgt, na.rm = TRUE),
                    bEst = sum(bEst* wgt, na.rm = TRUE),
                    tTEst = sum(tTEst* wgt, na.rm = TRUE), ## Need to sum this
                    bTEst = sum(bTEst* wgt, na.rm = TRUE), ## Need to sum this
                    tTVar = sum(tTVar* wgt^2, na.rm = TRUE),
                    bTVar = sum(bTVar* wgt^2, na.rm = TRUE),
                    cvEst_tT = sum(cvEst_tT* wgt^2, na.rm = TRUE),
                    cvEst_bT = sum(cvEst_bT* wgt^2, na.rm = TRUE),
                    ## Variances
                    tVar = sum(tVar* wgt^2, na.rm = TRUE),
                    bVar = sum(bVar* wgt^2, na.rm = TRUE),
                    #aVar = first(aVar),
                    cvEst_t = sum(cvEst_t* wgt^2, na.rm = TRUE),
                    cvEst_b = sum(cvEst_b* wgt^2, na.rm = TRUE),
                    plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))


        #### ----- EXPONENTIAL MOVING AVERAGE
      } else if (method == 'EMA'){
        ## Assuming a uniform weighting scheme
        popOrig <- mutate(popOrig, YEAR = END_INVYR)
        evals <- popOrig %>%
          distinct(YEAR, INVYR, STATECD, .keep_all = TRUE) %>%
          select(YEAR, INVYR, STATECD) %>%
          mutate(yrs = yrs)

        # Area
        aEst <- left_join(evals, aEst, by = c('INVYR' = 'YEAR', 'STATECD')) %>%
          #left_join(minEval, by = 'STATECD') %>%
          ## Make a previous year column, NA if earliest recorded
          #mutate(PREV_INVYR = ifelse(INVYR == minYear, NA_integer_, as.numeric(INVYR) -1))# %>%
          #%>%
          group_by(STATECD, .dots = aGrpBy) %>%
          summarize(aEst = ema(aEst, yrs),
                    aVar = ema(aVar, yrs, var = TRUE),
                    plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

        tEst <- left_join(evals, tEst, by = c('INVYR' = 'YEAR', 'STATECD')) %>%
          group_by(STATECD, .dots = grpBy) %>%
          #left_join(aTotal, by = c(aGrpBy)) %>%
          summarize(tEst = ema(tEst, yrs),
                    bEst = ema(bEst, yrs),
                    tTEst = ema(tTEst, yrs), ## Need to sum this
                    bTEst = ema(bTEst, yrs), ## Need to sum this
                    tTVar = ema(tTVar, yrs, var = TRUE),
                    bTVar = ema(bTVar, yrs, var = TRUE),
                    cvEst_tT = ema(cvEst_tT, yrs, var = TRUE),
                    cvEst_bT = ema(cvEst_bT, yrs, var = TRUE),
                    ## Variances
                    tVar = ema(tVar, yrs, var = TRUE),
                    bVar = ema(bVar, yrs, var = TRUE),
                    #aVar = first(aVar),
                    cvEst_t = ema(cvEst_t, yrs, var = TRUE),
                    cvEst_b = ema(cvEst_b, yrs, var = TRUE),
                    plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

      }

    }




    ##---------------------  TOTALS and RATIOS
    # Area
    aTotal <- aEst %>%
      group_by(.dots = aGrpBy) %>%
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                aVar = sum(aVar, na.rm = TRUE),
                AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
                nPlots_AREA = sum(plotIn_AREA, na.rm = TRUE))
    # Tree
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      #left_join(aTotal, by = c(aGrpBy)) %>%
      summarize(TREE_TOTAL = sum(tEst, na.rm = TRUE),
                BA_TOTAL = sum(bEst, na.rm = TRUE),
                ## Variances
                treeVar = sum(tVar, na.rm = TRUE),
                baVar = sum(bVar, na.rm = TRUE),
                #aVar = first(aVar),
                cvT = sum(cvEst_t, na.rm = TRUE),
                cvB = sum(cvEst_b, na.rm = TRUE),
                ## Sampling Errors
                TREE_SE = sqrt(treeVar) / TREE_TOTAL * 100,
                BA_SE = sqrt(baVar) / BA_TOTAL * 100,
                nPlots_TREE = sum(plotIn_TREE, na.rm = TRUE)) #%>%
    ## IF using polys, we treat each zone as a unique population
    if (!is.null(polys)){
      propGrp <- c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
    } else {
      propGrp <- 'YEAR'
    }
    ## Hand the proportions
    tpTotal <- tEst %>%
      group_by(.dots = propGrp) %>%
      summarize(TREE_TOTAL_full = sum(tTEst, na.rm = TRUE), ## Need to sum this
                BA_TOTAL_full = sum(bTEst, na.rm = TRUE), ## Need to sum this
                tTVar = sum(tTVar, na.rm = TRUE),
                bTVar = sum(bTVar, na.rm = TRUE),
                cvTT = sum(cvEst_tT, na.rm = TRUE),
                cvBT = sum(cvEst_bT, na.rm = TRUE))


    suppressWarnings({
      ## Bring them together
      tTotal <- tTotal %>%
        left_join(aTotal, by = aGrpBy) %>%
        left_join(tpTotal, by = propGrp) %>%
        mutate(TPA = TREE_TOTAL / AREA_TOTAL,
               BAA = BA_TOTAL / AREA_TOTAL,
               tpaVar = (1/AREA_TOTAL^2) * (treeVar + (TPA^2 * aVar) - 2 * TPA * cvT),
               baaVar = (1/AREA_TOTAL^2) * (baVar + (BAA^2 * aVar) - (2 * BAA * cvB)),
               TPA_SE = sqrt(tpaVar) / TPA * 100,
               BAA_SE = sqrt(baaVar) / BAA * 100,
               TPA_PERC = TREE_TOTAL / (TREE_TOTAL_full) * 100,
               BAA_PERC = BA_TOTAL / (BA_TOTAL_full) * 100,
               tpVar = (1/TREE_TOTAL_full^2) * (treeVar + (TPA_PERC^2 * tTVar) - 2 * TPA_PERC * cvTT),
               bpVar = (1/BA_TOTAL_full^2) * (baVar + (BAA_PERC^2 * bTVar) - (2 * BAA_PERC * cvBT)),
               TPA_PERC_SE = sqrt(tpVar) / TPA_PERC * 100,
               BAA_PERC_SE = sqrt(bpVar) / BAA_PERC * 100)
    })


    if (totals) {
      tOut <- tTotal %>%
        select(grpBy, TPA, BAA, TPA_PERC, BAA_PERC, TREE_TOTAL, BA_TOTAL, AREA_TOTAL, TPA_SE, BAA_SE,
               TPA_PERC_SE, BAA_PERC_SE, TREE_SE, BA_SE, AREA_TOTAL_SE, nPlots_TREE, nPlots_AREA)
    } else {
      tOut <- tTotal %>%
        select(grpBy, TPA, BAA, TPA_PERC, BAA_PERC,  TPA_SE, BAA_SE,
               TPA_PERC_SE, BAA_PERC_SE, nPlots_TREE, nPlots_AREA)
    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]


    # Return a spatial object
    if (!is.null(polys) & returnSpatial) {
      suppressMessages({suppressWarnings({tOut <- left_join(polys, tOut) %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})
    } else if (!is.null(polys) & returnSpatial == FALSE){
      tOut <- select(tOut, c('YEAR', grpByOrig, tNames, everything())) %>%
        filter(!is.na(polyID))
    }
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## Pretty output
  tOut <- drop_na(tOut, grpBy) %>%
    arrange(YEAR) %>%
    as_tibble()


  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
  return(tOut)
}















  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf') %>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    grpBy <- c(grpBy, 'LON', 'LAT')
  }


  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
  }
  # Tree Type domain indicator
  if (tolower(treeType) == 'live'){
    db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
  } else if (tolower(treeType) == 'dead'){
    db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 2 & db$TREE$STANDING_DEAD_CD == 1, 1, 0)
  } else if (tolower(treeType) == 'gs'){
    db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1 & db$TREE$DIA >= 5 & db$TREE$TREECLCD == 2, 1, 0)
  } else if (tolower(treeType) == 'all'){
    db$TREE$typeD <- 1
  }
  # update spatial domain indicator
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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

  ## Prep joins and filters
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)
    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))

  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'INVYR', 'MACRO_BREAKPOINT_DIA', grpP, 'sp', 'aD_p')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'landD', 'aD_c')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      left_join(select(db_clip$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', 'typeD', 'tD',
                                       'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG', 'DRYBIO_BG', 'CARBON_AG', 'CARBON_BG', grpT)), by = c('PLT_CN', 'CONDID')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM) %>%
      mutate_if(is.factor,
                as.character)%>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE)
    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }

    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
    data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp

    ## Add species to groups
    if (bySpecies) {
      data <- data %>%
        left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
        mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' '))%>%
        mutate_if(is.factor,
                  as.character)
      grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
      grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    }

    ## Break into size classes
    if (bySizeClass){
      grpBy <- c(grpBy, 'sizeClass')
      grpByOrig <- c(grpByOrig, 'sizeClass')
      data$sizeClass <- makeClasses(data$DIA, interval = 2)
      data <- data[!is.na(data$sizeClass),]
    }


    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {
      bOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
        # Compute estimates at plot level
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(NETVOL_ACRE = sum(VOLCFNET * TPA_UNADJ * tDI, na.rm = TRUE),
                  SAWVOL_ACRE = sum(VOLCSNET * TPA_UNADJ * tDI, na.rm = TRUE),
                  BIO_AG_ACRE = sum(DRYBIO_AG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                  BIO_BG_ACRE = sum(DRYBIO_BG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                  BIO_ACRE = sum(BIO_AG_ACRE, BIO_BG_ACRE, na.rm = TRUE),
                  CARB_AG_ACRE = sum(CARBON_AG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                  CARB_BG_ACRE = sum(CARBON_BG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                  CARB_ACRE = sum(CARB_AG_ACRE, CARB_BG_ACRE, na.rm = TRUE),
                  nStems = length(which(tDI == 1)))

      if (returnSpatial){
        bOut <- bOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
    } else {
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
      if (!is.null(polys)){
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND) | grpBy %in% names(pltSF)])
      } else {
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND)])
      }

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
            library(tidyr)
          })
          bOut <- parLapply(cl, X = names(combos), fun = biomassHelper, combos, data, grpBy, aGrpBy, totals, SE)
          stopCluster(cl)
        } else { # Unix systems
          bOut <- mclapply(X = names(combos), FUN = biomassHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
        }
      })

      if (SE){
        # Convert from list to dataframe
        bOut <- do.call(rbind,bOut)
      } else {
        # Pull out dataframe
        bOut <- bOut[[1]]
      }

      # Snag some names for below
      bNames <- names(bOut)[names(bOut) %in% grpBy == FALSE]

      # Return a spatial object
      if ('YEAR' %in% names(bOut)){
        # Return a spatial object
        if (!is.null(polys) & returnSpatial) {
          suppressMessages({suppressWarnings({bOut <- left_join(polys, bOut) %>%
            select(c(grpByOrig, bNames, names(polys))) %>%
            filter(!is.na(polyID))})})
        } else if (!is.null(polys) & returnSpatial == FALSE){
          bOut <- select(bOut, c(grpByOrig, bNames, everything())) %>%
            filter(!is.na(polyID))
        }
      } else { ## Function found no plots within the polygon, so it panics
        combos <- data %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize()
        bOut <- data.frame("YEAR" = combos$YEAR, "NETVOL_ACRE" = rep(NA, nrow(combos)),
                           "SAWVOL_ACRE" = rep(NA, nrow(combos)), "BIO_AG_ACRE" = rep(NA, nrow(combos)),
                           "BIO_BG_ACRE" = rep(NA, nrow(combos)), "BIO_ACRE" = rep(NA, nrow(combos)),
                           "CARB_AG_ACRE" = rep(NA, nrow(combos)),"CARB_BG_ACRE" = rep(NA, nrow(combos)),
                           "CARB_ACRE" = rep(NA, nrow(combos)), "NETVOL_ACRE_SE" = rep(NA, nrow(combos)),
                           "SAWVOL_ACRE_SE"  = rep(NA, nrow(combos)), "BIO_AG_ACRE_SE"  = rep(NA, nrow(combos)),
                           "BIO_BG_ACRE_SE"  = rep(NA, nrow(combos)), "BIO_ACRE_SE"  = rep(NA, nrow(combos)),
                           "CARB_AG_ACRE_SE" = rep(NA, nrow(combos)), "CARB_BG_ACRE_SE" = rep(NA, nrow(combos)),
                           "CARB_ACRE_SE"  = rep(NA, nrow(combos)), "nPlots_VOL" = rep(NA, nrow(combos)),
                           "nPlots_AREA" = rep(NA, nrow(combos)))
        if (!is.null(polys) & returnSpatial) {
          suppressMessages({suppressWarnings({
            polys = left_join(polys, combos)
            bOut <- left_join(polys, bOut)})})
        } else if (!is.null(polys) & returnSpatial == FALSE){
          bOut <- select(bOut, c(grpByOrig, everything()))
        }
      }
    } # End byPlot == FALSE
    out[[y]] <- bOut
    pb$tick()
  }
  bOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  bOut <- drop_na(bOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR) %>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) bOut <- st_sf(bOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) bOut <- unique(bOut)
  return(bOut)
}
