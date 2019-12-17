
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
                     lambda = .94,
                     stateVar = 'TPA',
                     treeDomain = NULL,
                     areaDomain = NULL,
                     totals = FALSE,
                     byPlot = FALSE,
                     nCores = 1) {

  ## Need a plotCN
  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Check if column exists
    allNames <- c(names(db$PLOT), names(db$COND), names(db$TREE))

    if (quo_name(grpBy_quo) %in% allNames){
      # Convert to character
      grpBy <- quo_name(grpBy_quo)
    } else {
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT, TREE, or COND tables. Did you accidentally quote the variables names? e.g. use grpBy = ECOSUBCD (correct) instead of grpBy = "ECOSUBCD". ', collapse = ', '))
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

  ## No EXP_GROW available for Western States, make sure we warn that values will be returned as 0
  # These states do not allow temporal queries. Things are extremely weird with their eval groups
  noGrow <- c(02,03,04,07,08,11,14,15,16, 23, 30, 32, 35,43,49, 78)
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


  ### HANDLE THE STATE VARIABLE, only applying to the midpoint table for consistency
  if (str_to_upper(stateVar) == 'TPA'){
   db$TREE_GRM_MIDPT$state <- 1
  } else if (str_to_upper(stateVar) == 'BAA'){
    db$TREE_GRM_MIDPT$state <- basalArea(db$TREE_GRM_MIDPT$DIA)
  } else if (str_to_upper(stateVar) == 'SAWVOL'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLCSNET
  } else if (str_to_upper(stateVar) == 'NETVOL'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$VOLCFNET
  } else if (str_to_upper(stateVar) == 'BIO_AG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_AG
  } else if (str_to_upper(stateVar) == 'BIO_BG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_BG
  } else if (str_to_upper(stateVar) == 'BIO'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYBIO_BG + db$TREE_GRM_MIDPT$DRYBIO_AG
  } else if (str_to_upper(stateVar) == 'CARB_AG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYCARB_AG * .5
  } else if (str_to_upper(stateVar) == 'CARB_BG'){
    db$TREE_GRM_MIDPT$state <- db$TREE_GRM_MIDPT$DRYCARB_BG * .5
  } else if (str_to_upper(stateVar) == 'CARB'){
    db$TREE_GRM_MIDPT$state <- (db$TREE_GRM_MIDPT$DRYCARB_BG + db$TREE_GRM_MIDPT$DRYCARB_AG) * .5
  } else {
    stop(paste('Method not known for stateVar:', stateVar))
  }

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'all'){
      db$TREE$typeD <- 1
      ## Rename some variables in grm
      db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                      TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_FOREST,
                                      TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_FOREST,
                                      TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_FOREST,
                                      SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_FOREST,
                                      COMPONENT = SUBP_COMPONENT_AL_FOREST)

    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
      db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                      TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_FOREST,
                                      TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_FOREST,
                                      TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_FOREST,
                                      SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_FOREST,
                                      COMPONENT = SUBP_COMPONENT_GS_FOREST)
    }
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'all'){
      db$TREE$typeD <- 1
      ## Rename some variables in grm
      db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                      TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_TIMBER,
                                      TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_TIMBER,
                                      TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_TIMBER,
                                      SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_TIMBER,
                                      COMPONENT = SUBP_COMPONENT_AL_TIMBER)

    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
      db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                      TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_TIMBER,
                                      TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_TIMBER,
                                      TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_TIMBER,
                                      SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_TIMBER,
                                      COMPONENT = SUBP_COMPONENT_GS_TIMBER)
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
  db$POP_EVAL<- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD', 'GROWTH_ACCT') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)# %>%
  #group_by(END_INVYR) %>%
  #summarise(id = list(EVALID)

  ## Make an annual panel ID, associated with an INVYR

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'GROWTH_ACCT', 'END_INVYR')) %>%
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
      filter(END_INVYR == INVYR) %>%
      ungroup()

    prePops <- popOrig %>%
      group_by(STATECD) %>%
      filter(INVYR < min(END_INVYR, na.rm = TRUE)) %>%
      distinct(PLT_CN, .keep_all = TRUE) %>%
      ungroup()

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
      out <- parLapply(cl, X = names(plts), fun = gmHelper1, plts, db, grpBy, aGrpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = gmHelper1, plts, db, grpBy, aGrpBy, byPlot, mc.cores = nCores)
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

      ## Modify some names if a different state variable was given
      names(tOut) <- str_replace(names(tOut), 'TPA', stateVar)

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
        out <- parLapply(cl, X = names(popState), fun = gmHelper2, popState, a, t, grpBy, aGrpBy)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = gmHelper2, popState, a, t, grpBy, aGrpBy, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    aEst <- bind_rows(out[names(out) == 'aEst'])
    tEst <- bind_rows(out[names(out) == 'tEst'])


    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA')){
      ## Need a STATECD on aEst and tEst to join wgts
      if ('STATECD' %in% names(tEst) == FALSE){
        ## Need a STATECD on aEst and tEst to join wgts
        tEst <- left_join(tEst, select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN'))
        aEst <- left_join(aEst, select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN'))
      }

      #### Summarizing to state level here to apply weights by panel
      #### Getting rid of ESTN_UNITS
      # Area
      aEst <- aEst %>%
        group_by(STATECD, .dots = aGrpBy) %>%
        summarize(aEst = sum(aEst, na.rm = TRUE),
                  aVar = sum(aVar, na.rm = TRUE),
                  plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))
      # Tree
      tEst <- tEst %>%
        group_by(STATECD, .dots = grpBy) %>%
        summarize_at(vars(tEst:plotIn_h),sum, na.rm = TRUE)

      ## Naming
      popOrig <- mutate(popOrig, YEAR = END_INVYR)

      ### ---- SIMPLE MOVING AVERAGE
      if (str_to_upper(method) == 'SMA'){
        ## Assuming a uniform weighting scheme
        wgts <- popOrig %>%
          group_by(YEAR, STATECD) %>%
          summarize(wgt = 1 / length(unique(INVYR))) %>%
          ## Expand it out again
          right_join(popOrig, by = c('YEAR', 'STATECD')) %>%
          distinct(YEAR, INVYR, STATECD, .keep_all = TRUE) %>%
          select(YEAR, INVYR, STATECD, wgt)

        #### ----- Linear MOVING AVERAGE
      } else if (str_to_upper(method) == 'LMA'){
        wgts <- popOrig %>%
          group_by(YEAR, STATECD) %>%
          summarize(n = length(unique(INVYR)),
                    minyr = min(INVYR, na.rm = TRUE)) %>%
          ## Expand it out again
          right_join(popOrig, by = c('YEAR', 'STATECD')) %>%
          distinct(YEAR, INVYR, STATECD, .keep_all = TRUE) %>%
          mutate(wgt = YEAR - minyr / sum(1:n, na.rm = TRUE)) %>%
          select(YEAR, INVYR, STATECD, wgt)


        #### ----- EXPONENTIAL MOVING AVERAGE
      } else if (str_to_upper(method) == 'EMA'){
        wgts <- popOrig %>%
          distinct(YEAR, INVYR, STATECD, .keep_all = TRUE) %>%
          select(YEAR, INVYR, STATECD)
        if (length(lambda) < 2){
          ## Weights based on temporal window
          wgts <- wgts %>%
            mutate(lambda = lambda,
                   l = lambda,
                   yrPrev = YEAR - INVYR,
                   wgt = l^yrPrev *(1-l))
        } else {
          grpBy <- c('lambda', grpBy)
          aGrpBy <- c('lambda', aGrpBy)
          ## Duplicate weights for each level of lambda
          yrWgts <- list()
          for (i in 1:length(unique(lambda))) {
            yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
          }
          wgts <- bind_rows(yrWgts) %>%
            mutate(l = lambda,
                   #l = 2 / (1 + lambda),
                   yrPrev = YEAR - INVYR,
                   wgt = l^yrPrev *(1-l))
        }

      }

      ### Applying the weights
      # Area
      aEst <- left_join(wgts, aEst, by = c('INVYR' = 'YEAR', 'STATECD')) %>%
        mutate_at(vars(aEst), ~(.*wgt)) %>%
        mutate_at(vars(aVar), ~(.*(wgt^2))) %>%
        group_by(STATECD, .dots = grpBy) %>%
        summarize_at(vars(aEst:plotIn_AREA), sum, na.rm = TRUE)


      tEst <- left_join(wgts, tEst, by = c('INVYR' = 'YEAR', 'STATECD')) %>%
        mutate_at(vars(tEst:hEst), ~(.*wgt)) %>%
        mutate_at(vars(tVar:cvEst_hp), ~(.*(wgt^2))) %>%
        group_by(STATECD, .dots = grpBy) %>%
        summarize_at(vars(tEst:plotIn_h), sum, na.rm = TRUE)

    }

    ##---------------------  TOTALS and RATIOS
    # Area
    # aTotal <- aEst %>%
    #   group_by(.dots = aGrpBy) %>%
    #   summarize(aEst = sum(aEst, na.rm = TRUE),
    #             aVar = sum(aVar, na.rm = TRUE),
    #             #AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
    #             plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))
    aTotal <- aEst %>%
      group_by(.dots = aGrpBy) %>%
      summarize_all(sum,na.rm = TRUE)
    # summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
    #           aVar = sum(aVar, na.rm = TRUE),
    #           AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
    #           nPlots_AREA = sum(plotIn_AREA, na.rm = TRUE))
    # Tree
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    suppressWarnings({
      ## Bring them together
      tOut <- tTotal %>%
        left_join(aTotal, by = aGrpBy) %>%
        # Renaming, computing ratios, and SE
        mutate(TREE_TOTAL = tEst,
               RECR_TOTAL = rEst,
               MORT_TOTAL = mEst,
               REMV_TOTAL = hEst,
               AREA_TOTAL = aEst,
               RECR_TPA = RECR_TOTAL / AREA_TOTAL,
               MORT_TPA = MORT_TOTAL / AREA_TOTAL,
               REMV_TPA = REMV_TOTAL / AREA_TOTAL,
               RECR_PERC = RECR_TOTAL / TREE_TOTAL * 100,
               MORT_PERC = MORT_TOTAL / TREE_TOTAL * 100,
               REMV_PERC = REMV_TOTAL / TREE_TOTAL * 100,
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
               ## SE TOTAL
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL *100,
               TREE_TOTAL_SE = sqrt(tVar) / TREE_TOTAL *100,
               RECR_TOTAL_SE = sqrt(rVar) / RECR_TOTAL *100,
               MORT_TOTAL_SE = sqrt(mVar) / MORT_TOTAL *100,
               REMV_TOTAL_SE = sqrt(hVar) / REMV_TOTAL *100,
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
      tOut <- tOut %>%
        select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
               TREE_TOTAL, RECR_TOTAL, MORT_TOTAL, REMV_TOTAL, AREA_TOTAL,
               RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE, RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE,
               TREE_TOTAL_SE, RECR_TOTAL_SE, MORT_TOTAL_SE, REMV_TOTAL_SE, AREA_TOTAL_SE,
               nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV, nPlots_AREA)
    } else {
      tOut <- tOut %>%
        select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
               RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE, RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE,
               nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV,nPlots_AREA)
    }


    ## Modify some names if a different state variable was given
    names(tOut) <- str_replace(names(tOut), 'TPA', stateVar)

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

  }

  ## Pretty output
  tOut <- drop_na(tOut, grpBy) %>%
    arrange(YEAR) %>%
    as_tibble()


  # Return a spatial object
  if (!is.null(polys)) {
    ## We don't like missing polygons through time, this makes them NA instead
    combos <- select(tOut, grpBy) %>%
      distinct() %>%
      mutate_all(as.factor) %>%
      group_by(.dots = grpBy, .drop = FALSE) %>%
      summarize() %>%
      ungroup() %>%
      mutate_all(as.character)

    combos <- matchColClasses(tOut, combos)
    tOut <- left_join(combos, tOut, by = grpBy)

    suppressMessages({suppressWarnings({
      tOut <- left_join(tOut, polys) %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
  return(tOut)
}










#
#   ## Which grpByNames are in which table? Helps us subset below
#   grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
#   grpC <- names(db$COND)[names(db$COND) %in% grpBy]
#   grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]
#
#
#   ### Snag the EVALIDs that are needed
#   ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
#   ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
#   if (!is.null(polys)){
#     ids <- pltSF %>%
#       left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
#       left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID', 'GROWTH_ACCT'), by = 'EVALID') %>%
#       inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
#       filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
#       filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
#       distinct(polyID, END_INVYR, EVALID, .keep_all = TRUE)
#
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
#       summarise(id = list(EVALID),
#                 ga = if_else(any(GROWTH_ACCT == 'Y'), 1, 0))
#   } else {
#     ids <- db$POP_EVAL %>%
#       select('CN', 'END_INVYR', 'EVALID', 'GROWTH_ACCT') %>%
#       inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
#       filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
#       filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
#       distinct(END_INVYR, EVALID, .keep_all = TRUE) %>%
#       group_by(END_INVYR) %>%
#       summarise(id = list(EVALID),
#                 ga = if_else(any(GROWTH_ACCT == 'Y'), 1, 0))
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
#     data <- select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')) %>%
#       left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
#       left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
#       right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
#       left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
#       left_join(select(db_clip$PLOT, c('PLT_CN', 'PREV_PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'REMPER', grpP, 'sp', 'aD_p')), by = 'PLT_CN') %>%
#       left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID','landD', 'aD_c', grpC)), by = 'PLT_CN') %>%
#       left_join(select(db_clip$TREE, c('TRE_CN', 'PREV_TRE_CN', 'TREE', 'PLT_CN', 'CONDID', 'PREVCOND', 'SPCD', grpT, 'typeD', 'tD')), by = c('PLT_CN', 'CONDID')) %>% ## ISSUE MAY BE HERE, SEE EVALIDATOR CODE
#       # GRM
#       left_join(select(db_clip$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
#       mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
#       mutate(tAdj = grmAdj(SUBPTYP_GRM, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
#       rename(YEAR = END_INVYR,
#              YEAR_RANGE = REPORT_YEAR_NM)
#     # Recode a few of the estimation methods to make things easier below
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
#
#     ## MODIFY FOR GROWTH ACCOUNTING
#     if (ids$ga[y]){
#       # Only subplots from cond change matrix
#       db_clip$SUBP_COND_CHNG_MTRX <- filter(db_clip$SUBP_COND_CHNG_MTRX, SUBPTYP == 1)
#
#       # Previous attributes
#       data <- data %>%
#         left_join(select(db_clip$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'CONDID'), suffix = c('', '.subp')) %>%
#         left_join(select(db_clip$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN.subp' = 'PLT_CN', 'PREVCOND.subp' = 'CONDID'), suffix = c('', '.chng')) %>%
#         left_join(select(db_clip$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
#         left_join(select(db_clip$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
#         left_join(select(db_clip$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
#         mutate_if(is.factor,
#                   as.character) %>%
#         mutate(aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.chng == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
#                tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0))
#
#       #If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
#       data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
#       data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
#       data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
#       data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
#       data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
#       data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)
#
#       ## Comprehensive indicator function
#       data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp * data$aChng
#       data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng
#       #data$aDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$sp.prev * data$aChng
#       # data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp * data$tChng
#
#       ## If growth accounting, summing proportions at subplot level, otherwise full plot
#       chngAdj <- .25
#
#       ## No growth accounting
#     } else {
#       # Previous attributes
#       data <- data %>%
#         left_join(select(db_clip$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
#         left_join(select(db_clip$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC)), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
#         left_join(select(db_clip$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
#         mutate_if(is.factor,
#                   as.character) %>%
#         ## Rename CONDPROP_UNADJ for consistency with above
#         rename(SUBPTYP_PROP_CHNG = CONDPROP_UNADJ)
#       ## Plot level proportion
#       chngAdj <- 1
#
#       # If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
#       data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
#       data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
#       data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
#       data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
#       data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
#       data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)
#
#       ## Comprehensive indicator function
#       data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
#       data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev
#
#     }
#
#
#
#     ## Add species to groups
#     if (bySpecies) {
#       data <- data %>%
#         left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
#         mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' ')) %>%
#         mutate_if(is.factor,
#                   as.character)
#       grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
#       grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
#     }
#
#     ## Break into size classes
#     if (bySizeClass){
#       grpBy <- c(grpBy, 'sizeClass')
#       grpByOrig <- c(grpByOrig, 'sizeClass')
#       data$sizeClass <- makeClasses(data$DIA, interval = 2)
#       data <- data[!is.na(data$sizeClass),]
#     }
#
#
#
#     ####################  COMPUTE ESTIMATES  ###########################
#     ### -- BYPLOT -- TPA Estimates at each plot location
#     if (byPlot) {
#       tOut <- data %>%
#         mutate(YEAR = INVYR) %>%
#         distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
#         #filter(EVALID %in% tID) %>%
#         # Compute estimates at plot level
#         group_by(.dots = grpBy, PLT_CN) %>%
#         summarize(TOTAL_TPA = sum(TPAGROW_UNADJ * tDI, na.rm = TRUE),
#                   RECR_TPA = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'] / REMPER[COMPONENT == 'INGROWTH'], na.rm = TRUE),
#                   MORT_TPA = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
#                   REMV_TPA = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
#                   RECR_PERC = RECR_TPA / TOTAL_TPA * 100,
#                   MORT_PERC = MORT_TPA / TOTAL_TPA * 100,
#                   REMV_PERC = REMV_TPA / TOTAL_TPA * 100,
#                   nStems = length(which(tDI == 1)))
#
#       if (returnSpatial){
#         tOut <- tOut %>%
#           filter(!is.na(LAT) & !is.na(LON)) %>%
#           st_as_sf(coords = c('LON', 'LAT'),
#                    crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
#
#       }
#
#     } else {
#       # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
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
#
#       # List of rows for lapply
#       combos <- split(combos, seq(nrow(combos)))
#
#       # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
#       if (!is.null(polys)){
#         aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND) | grpBy %in% names(pltSF)])
#       } else {
#         aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND)])
#       }
#
#       suppressWarnings({
#         ## Compute estimates in parallel -- Clusters in windows, forking otherwise
#         if (Sys.info()['sysname'] == 'Windows'){
#           cl <- makeCluster(nCores)
#           clusterEvalQ(cl, {
#             library(dplyr)
#             library(stringr)
#           })
#           tOut <- parLapply(cl, X = names(combos), fun = growMortHelper, combos, data, grpBy, aGrpBy, totals, SE, chngAdj)
#           stopCluster(cl)
#         } else { # Unix systems
#           tOut <- mclapply(X = names(combos), FUN = growMortHelper, combos, data, grpBy, aGrpBy, totals, SE, chngAdj, mc.cores = nCores)
#         }
#       })
#
#       if (SE){
#         # Convert from list to dataframe
#         tOut <- do.call(rbind,tOut)
#       } else {
#         # Pull out dataframe
#         tOut <- tOut[[1]]
#       }
#
#
#       # Snag some names for below
#       tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]
#
#       # Return a spatial object
#       if ('YEAR' %in% names(tOut)){
#         # Return a spatial object
#         if (!is.null(polys) & returnSpatial) {
#           suppressMessages({suppressWarnings({tOut <- left_join(polys, tOut) %>%
#             select(c(grpByOrig, tNames, names(polys))) %>%
#             filter(!is.na(polyID))})})
#         } else if (!is.null(polys) & returnSpatial == FALSE){
#           tOut <- select(tOut, c(grpByOrig, tNames, everything())) %>%
#             filter(!is.na(polyID))
#         }
#       } else { ## Function found no plots within the polygon, so it panics
#         combos <- data %>%
#           as.data.frame() %>%
#           group_by(.dots = grpBy) %>%
#           summarize()
#         tOut <- data.frame("YEAR" = combos$YEAR, "RECR_TPA" = rep(NA, nrow(combos)),
#                            "MORT_TPA" = rep(NA, nrow(combos)), "REMV_TPA" = rep(NA, nrow(combos)),
#                            "RECR_PERC" = rep(NA, nrow(combos)), "MORT_PERC" = rep(NA, nrow(combos)),
#                            "REMV_PERC" = rep(NA, nrow(combos)), "RECR_TPA_SE" = rep(NA, nrow(combos)),
#                            "MORT_TPA_SE" = rep(NA, nrow(combos)), "REMV_TPA_SE" = rep(NA, nrow(combos)),
#                            "RECR_PERC_SE" = rep(NA, nrow(combos)), "MORT_PERC_SE" = rep(NA, nrow(combos)),
#                            "REMV_PERC_SE" = rep(NA, nrow(combos)), "nPlots_TREE" = rep(NA, nrow(combos)),
#                            "nPlots_RECR"= rep(NA, nrow(combos)), "nPlots_MORT"= rep(NA, nrow(combos)),
#                            "nPlots_REMV"= rep(NA, nrow(combos)), "nPlots_AREA" = rep(NA, nrow(combos)))
#         if (!is.null(polys) & returnSpatial) {
#           suppressMessages({suppressWarnings({
#             polys = left_join(polys, combos)
#             tOut <- left_join(polys, tOut)})})
#         } else if (!is.null(polys) & returnSpatial == FALSE){
#           tOut <- data.frame(select(tOut, -c('YEAR')), combos)
#         }
#       }
#     } # End byPlot == FALSE
#     out[[y]] <- tOut
#     pb$tick()
#   }
#   tOut <- do.call(rbind, out)
#   ## For spatial plots
#   if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
#   tOut <- drop_na(tOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
#     arrange(YEAR)%>%
#     as_tibble()
#   ## Above converts to tibble
#   if (returnSpatial) tOut <- st_sf(tOut)
#   # ## remove any duplicates in byPlot (artifact of END_INYR loop)
#   if (byPlot) tOut <- unique(tOut)
#   return(tOut)
# }
