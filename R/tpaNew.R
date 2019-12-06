tpaNew <- function(db,
                grpBy = NULL,
                polys = NULL,
                returnSpatial = FALSE,
                bySpecies = FALSE,
                bySizeClass = FALSE,
                landType = 'forest',
                treeType = 'live',
                treeDomain = NULL,
                areaDomain = NULL,
                totals = FALSE,
                byPlot = FALSE,
                SE = TRUE,
                nCores = 1) {

  ## Need a plotCN
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN)

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
        return(0)
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

  # Save original grpBy for pretty return with spatial objects
  #grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
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

  ### Snag the EVALIDs that are needed
  db$POP_EVAL<- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)# %>%
    #group_by(END_INVYR) %>%
    #summarise(id = list(EVALID))

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'END_INVYR')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = 'STRATUM_CN') %>%
    rename(YEAR = END_INVYR) %>%
    mutate_if(is.factor,
              as.character) %>%
    filter(!is.na(YEAR))

  # ## How many plots per stratum
  # pltStrat <- pops %>%
  #   group_by(ESTN_UNIT_CN, STRATUM_CN) %>%
  #   summarize(plts = length(unique(PLT_CN)),
  #             p2 = first(P2POINTCNT))
  # ## Count total plts per EU, some get dropped in grouping below
  # p2eu <- pops %>%
  #   distinct(ESTN_UNIT_CN, STRATUM_CN, P2POINTCNT) %>%
  #   group_by(ESTN_UNIT_CN) %>%
  #   summarize(p2eu = sum(P2POINTCNT, na.rm = TRUE))
  #
  # ### Add columns from above
  # data <- data %>%
  #   left_join(pltStrat, by = c('ESTN_UNIT_CN', 'STRATUM_CN')) %>%
  #   left_join(p2eu, by = c('ESTN_UNIT_CN'))


  # ## Grab some adjustment factors from stratum
  # adj <- db$POP_PLOT_STRATUM_ASSGN %>%
  #   left_join(db$POP_STRATUM, by = c('STRATUM_CN' = 'CN')) #%>%
  #   group_by(PLT_CN) %>%
  #   summarize(ADJ_FACTOR_SUBP = first(ADJ_FACTOR_SUBP),
  #             ADJ_FACTOR_SUBPm = mean(ADJ_FACTOR_SUBP),
  #             ADJ_FACTOR_SUBPv = var(ADJ_FACTOR_SUBP),
  #             ADJ_FACTOR_MICR = first(ADJ_FACTOR_MICR),
  #             ADJ_FACTOR_MACR = first(ADJ_FACTOR_MACR))


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', grpP, 'aD_p', 'sp')) %>%
    filter(PLT_CN %in% pops$PLT_CN) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', grpT, 'tD', 'typeD')), by = c('PLT_CN', 'CONDID')) %>%
    ## Need a code that tells us where the tree was measured
    ## macroplot, microplot, subplot
    mutate(PLOT_BASIS = case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA) ~ NA_character_,
      ## When DIA is less than 5", use microplot value
      DIA < 5 ~ 'MICR',
      ## When DIA is greater than 5", use subplot value
      DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
      DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) #%>%
    #filter(!is.na(PLOT_BASIS))
    # rename(YEAR = INVYR) %>%
    # mutate_if(is.factor,
    #           as.character) %>%
    # filter(!is.na(YEAR))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
  data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp


  ## Recode a few of the estimation methods to make things easier below
  pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')

  # ### Make sure appropriate adjustment factors are used here
  # data <- data %>%
  #   mutate(
  #     ## AREA
  #     aAdj = case_when(
  #       ## When NA, stay NA
  #       is.na(PROP_BASIS) ~ NA_real_,
  #       ## If the proportion was measured for a macroplot,
  #       ## use the macroplot value
  #       PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
  #       ## Otherwise, use the subpplot value
  #       PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
  #     ## TREE
  #     tAdj = case_when(
  #       ## When DIA is na, adjustment is NA
  #       is.na(DIA) ~ ADJ_FACTOR_SUBP,
  #       ## When DIA is less than 5", use microplot value
  #       DIA < 5 ~ ADJ_FACTOR_MICR,
  #       ## When DIA is greater than 5", use subplot value
  #       DIA >= 5 ~ ADJ_FACTOR_SUBP
  #     ))

  ## Add species to groups
  if (bySpecies) {
    data <- data %>%
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
    data$sizeClass <- makeClasses(data$DIA, interval = 2)
    data <- data[!is.na(data$sizeClass),]
  }

  # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
  if (!is.null(polys)){
    aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
  } else {
    aGrpBy <- c(grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
  }


  ################# ------------ AREA ESTIMATES
  #### Having trouble filtering out plots that are used for a particular EVAL_TYP
  ## area -- plot level
  a <- data %>%
    ## Will be lots of trees here, so CONDPROP listed multiple times
    ## Adding PROP_BASIS so we can handle adjustment factors at strata level
    distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
    group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
    summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
              plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

  ## Adding YEAR to groups
  aGrpBy <- c('YEAR', aGrpBy)

  ## Strata level estimates
  aStrat <- a %>%
    ## Rejoin with population tables
    right_join(pops, by = 'PLT_CN') %>%
    mutate(
      ## AREA
      aAdj = case_when(
        ## When NA, stay NA
        is.na(PROP_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = fa * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = aGrpBy) %>%
    summarize(a_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              aStrat = mean(fa * a_t, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              ndif = nh - n,
              ## Strata level variances
              av = ifelse(first(ESTN_METHOD == 'simple'),
                          var(c(fa, numeric(ndif)) * first(a) / nh),
                          (sum((c(fa, numeric(ndif))^2)) - sum(nh * aStrat^2)) / (nh * (nh-1))))
  ## Estimation unit
  aEst <- aStrat %>%
    group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
              aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))


  ######## ------------------ TREE ESTIMATES + CV

  ## Tree plts
  t <- data %>%
    filter(!is.na(PLOT_BASIS)) %>%
    group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
    summarize(tPlot = sum(TPA_UNADJ * tDI, na.rm = TRUE),
              bPlot = sum(basalArea(DIA) * TPA_UNADJ * tDI, na.rm = TRUE),
              tTPlot = sum(TPA_UNADJ * pDI, na.rm = TRUE),
              bTPlot = sum(basalArea(DIA) * TPA_UNADJ * pDI, na.rm = TRUE),
              plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0))

  ## Adding YEAR to groups
  grpBy <- c('YEAR', grpBy)

  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    right_join(pops, by = 'PLT_CN') %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, aGrpBy[aGrpBy %in% 'YEAR' == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% 'YEAR' == FALSE])) %>%
    #Add adjustment factors
    mutate(
      ## AREA
      tAdj = case_when(
        ## When NA, stay NA
        is.na(PLOT_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PLOT_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PLOT_BASIS == 'SUBP' ~ as.numeric(ADJ_FACTOR_SUBP),
        PLOT_BASIS == 'MICR' ~ as.numeric(ADJ_FACTOR_MICR)),
      tPlot = tPlot * tAdj,
      bPlot = bPlot * tAdj,
      tTPlot = tTPlot * tAdj,
      bTPlot = bTPlot * tAdj) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(tPlot = sum(tPlot, na.rm = TRUE),
              bPlot = sum(bPlot, na.rm = TRUE),
              tTPlot = sum(tTPlot, na.rm = TRUE),
              bTPlot = sum(bTPlot, na.rm = TRUE),
              fa = first(fa),
              plotIn = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              tStrat = mean(tPlot * r_t, na.rm = TRUE),
              bStrat = mean(bPlot * r_t, na.rm = TRUE),
              tTStrat = mean(tTPlot * r_t, na.rm = TRUE),
              bTStrat = mean(bTPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_TREE = sum(plotIn, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              ndif = nh - n,
              ## Strata level variances
              #aVar = (sum(forArea^2) - sum(P2POINTCNT * aStrat^2)) / (P2POINTCNT * (P2POINTCNT-1)),
              tv = ifelse(first(ESTN_METHOD == 'simple'),
                          var(c(tPlot, numeric(ndif)) * first(a) / nh),
                          (sum(c(tPlot, numeric(ndif))^2) - sum(nh * tStrat^2)) / (nh * (nh-1))), # Stratified and double cases
              bv = ifelse(first(ESTN_METHOD == 'simple'),
                          var(c(tPlot, numeric(ndif))* first(a) / nh),
                          (sum(c(bPlot, numeric(ndif))^2) - sum(nh * bStrat^2)) / (nh * (nh-1))),
              tTv = ifelse(first(ESTN_METHOD == 'simple'),
                           var(c(tTPlot, numeric(ndif)) * first(a) / nh),
                           (sum(c(tTPlot, numeric(ndif))^2) - sum(nh * tTStrat^2)) / (nh * (nh-1))), # Stratified and double cases
              bTv = ifelse(first(ESTN_METHOD == 'simple'),
                           var(c(bTPlot, numeric(ndif)) * first(a) / nh),
                           (sum(c(bTPlot, numeric(ndif))^2) - sum(nh * bTStrat^2)) / (nh * (nh-1))),
              # Strata level covariances
              cvStrat_t = ifelse(first(ESTN_METHOD == 'simple'),
                                 cov(fa,tPlot),
                                 (sum(fa*tPlot, na.rm = TRUE) - sum(nh * aStrat *tStrat)) / (nh * (nh-1))), # Stratified and double cases
              cvStrat_b = ifelse(first(ESTN_METHOD == 'simple'),
                                 cov(fa,bPlot),
                                 (sum(fa*bPlot, na.rm = TRUE) - sum(nh * aStrat *bStrat)) / (nh * (nh-1))),
              cvStrat_tT = ifelse(first(ESTN_METHOD == 'simple'),
                                  cov(tTPlot,tPlot),
                                  (sum(tTPlot*tPlot, na.rm = TRUE) - sum(nh * tTStrat *tStrat)) / (nh * (nh-1))), # Stratified and double cases
              cvStrat_bT = ifelse(first(ESTN_METHOD == 'simple'),
                                  cov(bTPlot,bPlot),
                                  (sum(bTPlot*bPlot, na.rm = TRUE) - sum(nh * bTStrat *bStrat)) / (nh * (nh-1)))) %>%
    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
              bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
              tTEst = unitMean(ESTN_METHOD, a, nh, w, tTStrat),
              bTEst = unitMean(ESTN_METHOD, a, nh, w, bTStrat),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE),
              tVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, tv, tStrat, tEst),
              bVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bv, bStrat, bEst),
              tTVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, tTv, tTStrat, tTEst),
              bTVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bTv, bTStrat, bTEst),
              # Unit Covariance
              cvEst_t = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_t, tStrat, tEst, aStrat, aEst),
              cvEst_b = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_b, bStrat, bEst, aStrat, aEst),
              cvEst_tT = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_t, tStrat, tEst, tTStrat, tTEst),
              cvEst_bT = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_b, bStrat, bEst, bTStrat, bTEst))



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
              TREE_TOTAL_full = sum(tTEst, na.rm = TRUE),
              BA_TOTAL_full = sum(bTEst, na.rm = TRUE),
              ## Ratios
              #TPA = TREE_TOTAL / AREA_TOTAL,
              #BAA = BA_TOTAL / AREA_TOTAL,
              TPA_PERC = TREE_TOTAL / TREE_TOTAL_full * 100,
              BAA_PERC = BA_TOTAL / BA_TOTAL_full * 100,
              ## Variances
              treeVar = sum(tVar, na.rm = TRUE),
              baVar = sum(bVar, na.rm = TRUE),
              tTVar = sum(tTVar, na.rm = TRUE),
              bTVar = sum(bTVar, na.rm = TRUE),
              #aVar = first(aVar),
              cvT = sum(cvEst_t, na.rm = TRUE),
              cvB = sum(cvEst_b, na.rm = TRUE),
              cvTT = sum(cvEst_tT, na.rm = TRUE),
              cvBT = sum(cvEst_bT, na.rm = TRUE),
              #tpaVar = (1/AREA_TOTAL^2) * (treeVar + (TPA^2 * aVar) - 2 * TPA * cvT),
              #baaVar = (1/AREA_TOTAL^2) * (baVar + (BAA^2 * aVar) - (2 * BAA * cvB)),
              tpVar = (1/TREE_TOTAL_full^2) * (treeVar + (TPA_PERC^2 * tTVar) - 2 * TPA_PERC * cvTT),
              bpVar = (1/BA_TOTAL_full^2) * (baVar + (BAA_PERC^2 * bTVar) - (2 * BAA_PERC * cvBT)),
              ## Sampling Errors
              TREE_SE = sqrt(treeVar) / TREE_TOTAL * 100,
              BA_SE = sqrt(baVar) / BA_TOTAL * 100,
              #AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
              #TPA_SE = sqrt(tpaVar) / TPA * 100,
              #BAA_SE = sqrt(baaVar) / BAA * 100,
              TPA_PERC_SE = sqrt(tpVar) / TPA_PERC * 100,
              BAA_PERC_SE = sqrt(bpVar) / BAA_PERC * 100,
              nPlots_TREE = sum(plotIn_TREE, na.rm = TRUE))
              #nPlots_AREA = first(nPlots_AREA))

  ## Bring them together
  tTotal <- tTotal %>%
    left_join(aTotal, by = aGrpBy) %>%
    mutate(TPA = TREE_TOTAL / AREA_TOTAL,
           BAA = BA_TOTAL / AREA_TOTAL,
           tpaVar = (1/AREA_TOTAL^2) * (treeVar + (TPA^2 * aVar) - 2 * TPA * cvT),
           baaVar = (1/AREA_TOTAL^2) * (baVar + (BAA^2 * aVar) - (2 * BAA * cvB)),
           TPA_SE = sqrt(tpaVar) / TPA * 100,
           BAA_SE = sqrt(baaVar) / BAA * 100)

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


  # # Return a spatial object
  # if ('YEAR' %in% names(tOut)){
  #   # Return a spatial object
  #   if (!is.null(polys) & returnSpatial) {
  #     suppressMessages({suppressWarnings({tOut <- left_join(polys, tOut) %>%
  #       select(c(grpByOrig, tNames, names(polys))) %>%
  #       filter(!is.na(polyID))})})
  #   } else if (!is.null(polys) & returnSpatial == FALSE){
  #     tOut <- select(tOut, c(grpByOrig, tNames, everything())) %>%
  #       filter(!is.na(polyID))
  #     # Return spatial plots
  #   }
  # } else { ## Function found no plots within the polygon, so it panics
  #   combos <- data %>%
  #     as.data.frame() %>%
  #     group_by(.dots = grpBy) %>%
  #     summarize()
  #   tOut <- data.frame("YEAR" = combos$YEAR, "TPA" = rep(NA, nrow(combos)),
  #                      "BAA" = rep(NA, nrow(combos)), "TPA_PERC" = rep(NA, nrow(combos)),
  #                      "BAA_PERC" = rep(NA, nrow(combos)), "TPA_SE" = rep(NA, nrow(combos)),
  #                      "BAA_SE" = rep(NA, nrow(combos)), "TPA_PERC_SE" = rep(NA, nrow(combos)),
  #                      "BAA_PERC_SE" = rep(NA, nrow(combos)),"nPlots_TREE" = rep(NA, nrow(combos)),
  #                      "nPlots_AREA" = rep(NA, nrow(combos)))
  #   if (!is.null(polys) & returnSpatial) {
  #     suppressMessages({suppressWarnings({
  #       polys = left_join(polys, combos)
  #       tOut <- left_join(polys, tOut)})})
  #   } else if (!is.null(polys) & returnSpatial == FALSE){
  #     tOut <- data.frame(select(tOut, -c('YEAR')), combos)
  #   }
  # }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  tOut <- drop_na(tOut, grpBy) %>%
    arrange(YEAR) %>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
  return(tOut)
}
