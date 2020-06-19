dwmStarter <- function(x,
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
                       tidy = TRUE,
                       nCores = 1,
                       remote){

  reqTables <- c('PLOT', 'COND_DWM_CALC', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')

  if (remote){
    ## Store the original parameters here
    params <- db

    ## Read in one state at a time
    db <- readFIA(dir = db$dir, common = db$common,
                  tables = reqTables, states = x, ## x is the vector of state names
                  nCores = nCores)

    ## If a clip was specified, run it now
    if ('mostRecent' %in% names(params)){
      db <- clipFIA(db, mostRecent = params$mostRecent,
                    mask = params$mask, matchEval = params$matchEval,
                    evalid = params$evalid, designCD = params$designCD,
                    nCores = nCores)
    }

  } else {
    ## Really only want the required tables
    db <- db[names(db) %in% reqTables]
  }



  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))
  db$COND_DWM_CALC <- db[['COND_DWM_CALC']] %>% mutate(DWM_CN = CN)
  db$COND <- db[['COND']] %>% mutate(CND_CN = CN)



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
        select(!!grpBy_quo)
    )

    # If column doesnt exist, just returns 0, not a dataframe
    if (is.null(nrow(d_quo))){
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT or COND tables. Did you accidentally quote the variables names? e.g. use grpBy = ECOSUBCD (correct) instead of grpBy = "ECOSUBCD". ', collapse = ', '))
    } else {
      # Convert to character
      grpBy <- names(d_quo)
    }
  } else {
    grpBy <- NULL
  }

  reqTables <- c('PLOT', 'COND_DWM_CALC', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  ## Some warnings

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

  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}
  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy

  ## IF the object was clipped
  if ('prev' %in% names(db$PLOT)){
    ## Only want the current plots, no grm
    db$PLOT <- filter(db$PLOT, prev == 0)
  }

  ### DEAL WITH TEXAS
  if (any(db$POP_EVAL$STATECD %in% 48)){
    ## Will require manual updates, fix your shit texas
    txIDS <- db$POP_EVAL %>%
      filter(STATECD %in% 48) %>%
      filter(END_INVYR < 2017) %>%
      filter(END_INVYR > 2006) %>%
      ## Removing any inventory that references east or west, sorry
      filter(str_detect(str_to_upper(EVAL_DESCR), 'EAST', negate = TRUE) &
               str_detect(str_to_upper(EVAL_DESCR), 'WEST', negate = TRUE))
    db$POP_EVAL <- bind_rows(filter(db$POP_EVAL, !(STATECD %in% 48)), txIDS)
  }

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # # Convert polygons to an sf object
    # polys <- polys %>%
    #   as('sf')%>%
    #   mutate_if(is.factor,
    #             as.character)
    # ## A unique ID
    # polys$polyID <- 1:nrow(polys)
    #
    # # Add shapefile names to grpBy
    grpBy = c(grpBy, 'polyID')

    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('LON', 'LAT', pltID)) %>%
      filter(!is.na(LAT) & !is.na(LON)) %>%
      distinct(pltID, .keep_all = TRUE)
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys))

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
      left_join(select(pltSF, polyID, pltID), by = 'pltID')

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
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
  }

  # update spatial domain indicator
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(db$PLOT$pltID %in% pltSF$pltID, 1, 0)
  } else {
    db$PLOT$sp <- 1
  }

  # User defined domain indicator for area (ex. specific forest type)
  pcEval <- left_join(db$PLOT, select(db$COND, -c('STATECD', 'UNITCD', 'COUNTYCD', 'INVYR', 'PLOT')), by = 'PLT_CN')
  #areaDomain <- substitute(areaDomain)
  pcEval$aD <- rlang::eval_tidy(areaDomain, pcEval) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
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
    filter(EVAL_TYP == 'EXPDWM') %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)# %>%
  #group_by(END_INVYR) %>%
  #summarise(id = list(EVALID)

  ## Make an annual panel ID, associated with an INVYR

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'END_INVYR', 'EVAL_TYP')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
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


  ## Only the necessary plots for EVAL of interest
  db$PLOT <- filter(db$PLOT, PLT_CN %in% pops$PLT_CN)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  db$PLOT <- select(db$PLOT, c('PLT_CN', 'STATECD', 'COUNTYCD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', grpP, 'aD_p', 'sp'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD'))
  db$COND_DWM_CALC <- select(db$COND_DWM_CALC, -c( 'STATECD', 'COUNTYCD', 'UNITCD', 'INVYR', 'MEASYEAR', 'PLOT', 'EVALID'))
  #filter(DIA >= 5)

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
      out <- parLapply(cl, X = names(plts), fun = dwmHelper1, plts, db, grpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = dwmHelper1, plts, db, grpBy, byPlot, mc.cores = nCores)
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

    out <- list(tEst = tOut, grpBy = grpBy, grpByOrig = grpByOrig)
    ## Population estimation
  } else {
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    t <- bind_rows(out[names(out) == 't'])


    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)


    ## Splitting up by STATECD and groups of 25 ESTN_UNIT_CNs
    #estunit <- distinct(pops, ESTN_UNIT_CN) #%>%
    #mutate(estID)

    #estID <- seq(1, nrow(estunit), 50)
    #estunit$estID <- rep_len(estID, length.out = nrow(estunit))
    #pops <- pops %>%
    #  left_join(estunit, by = 'ESTN_UNIT_CN') #%>%
    #mutate(estBreaks = )

    #popState <- split(pops, as.factor(pops$estID))
    popState <- split(pops, as.factor(pops$STATECD))
    #
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
        out <- parLapply(cl, X = names(popState), fun = dwmHelper2, popState, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = dwmHelper2, popState, t, grpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])

    out <- list(tEst = tEst, grpBy = grpBy, grpByOrig = grpByOrig)
  }

  return(out)

}


#' @export
dwm <- function(db,
                           grpBy = NULL,
                           polys = NULL,
                           returnSpatial = FALSE,
                           landType = 'forest',
                           method = 'TI',
                           lambda = .5,
                           areaDomain = NULL,
                           byPlot = FALSE,
                           totals = FALSE,
                           tidy = TRUE,
                           nCores = 1) {

  ##  don't have to change original code
  grpBy_quo <- rlang::enquo(grpBy)
  areaDomain <- rlang::enquo(areaDomain)

  ### Is DB remote?
  remote <- ifelse(class(db) == 'Remote.FIA.Database', 1, 0)
  if (remote){

    iter <- db$states

    ## In memory
  } else {
    ## Some warnings
    if (class(db) != "FIA.Database"){
      stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
    }

    ## an iterator for remote
    iter <- 1

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
  }



  ## Run the main portion
  out <- lapply(X = iter, FUN = dwmStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                landType, method,
                lambda, areaDomain,
                byPlot, totals, tidy,
                nCores, remote)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  tEst <- bind_rows(out[names(out) == 'tEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  grpByOrig <- out[names(out) == 'grpByOrig'][[1]]




  if (byPlot){
    tOut <- tEst

    ## Population estimates
  } else {

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
      tEst <- tEst %>%
        mutate_at(vars(aEst:cEst), ~(.*wgt)) %>%
        mutate_at(vars(aVar:cvEst_c), ~(.*(wgt^2))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize_at(vars(aEst:cvEst_c), sum, na.rm = TRUE)
    }

    ##---------------------  TOTALS and RATIOS
    suppressWarnings({
      ## Bring them together
      tOut <- ungroup(tEst) %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        # Renaming, computing ratios, and SE
        mutate(VOL_DUFF = NA,
               VOL_LITTER = NA,VOL_1HR = vsmEst,
               VOL_10HR = vmdEst,
               VOL_100HR = vlgEst,
               VOL_1000HR =  vcEst,
               VOL_PILE =  vpEst,
               VOL =  vEst,
               BIO_DUFF = bdEst,
               BIO_LITTER = blEst,
               BIO_1HR = bsmEst,
               BIO_10HR = bmdEst,
               BIO_100HR = blgEst,
               BIO_1000HR = bcEst,
               BIO_PILE = bpEst,
               BIO = bEst,
               CARB_DUFF = cdEst,
               CARB_LITTER = clEst,
               CARB_1HR = csmEst,
               CARB_10HR = cmdEst,
               CARB_100HR = clgEst,
               CARB_1000HR = ccEst,
               CARB_PILE = cpEst,
               CARB = cEst,
               AREA_TOTAL = aEst,
               ## RATIOS
               VOL_DUFF_ACRE = NA,
               VOL_LITTER_ACRE = NA,
               VOL_1HR_ACRE = vsmEst / aEst,
               VOL_10HR_ACRE = vmdEst / aEst,
               VOL_100HR_ACRE = vlgEst / aEst,
               VOL_1000HR_ACRE =  vcEst / aEst,
               VOL_PILE_ACRE =  vpEst / aEst,
               VOL_ACRE =  vEst / aEst,
               BIO_DUFF_ACRE = bdEst / aEst,
               BIO_LITTER_ACRE = blEst / aEst,
               BIO_1HR_ACRE = bsmEst / aEst,
               BIO_10HR_ACRE = bmdEst / aEst,
               BIO_100HR_ACRE = blgEst / aEst,
               BIO_1000HR_ACRE = bcEst / aEst,
               BIO_PILE_ACRE = bpEst / aEst,
               BIO_ACRE = bEst / aEst,
               CARB_DUFF_ACRE = cdEst / aEst,
               CARB_LITTER_ACRE = clEst / aEst,
               CARB_1HR_ACRE = csmEst / aEst,
               CARB_10HR_ACRE = cmdEst / aEst,
               CARB_100HR_ACRE = clgEst / aEst,
               CARB_1000HR_ACRE = ccEst / aEst,
               CARB_PILE_ACRE = cpEst / aEst,
               CARB_ACRE = cEst / aEst,
               # Sampling Errors totals
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
               VOL_DUFF_SE = NA,
               VOL_LITTER_SE = NA,
               VOL_1HR_SE = sqrt(vsmVar) / VOL_1HR * 100,
               VOL_10HR_SE = sqrt(vmdVar) / VOL_10HR * 100,
               VOL_100HR_SE = sqrt(vlgVar) / VOL_100HR * 100,
               VOL_1000HR_SE = sqrt(vcVar) / VOL_1000HR * 100,
               VOL_PILE_SE = sqrt(vpVar) / VOL_PILE * 100,
               VOL_SE = sqrt(vVar) / VOL * 100,
               BIO_DUFF_SE = sqrt(bdVar) / BIO_DUFF * 100,
               BIO_LITTER_SE = sqrt(blVar) / BIO_LITTER * 100,
               BIO_1HR_SE = sqrt(bsmVar) / BIO_1HR * 100,
               BIO_10HR_SE = sqrt(bmdVar) / BIO_10HR * 100,
               BIO_100HR_SE = sqrt(blgVar) / BIO_100HR * 100,
               BIO_1000HR_SE = sqrt(bcVar) / BIO_1000HR * 100,
               BIO_PILE_SE = sqrt(bpVar) / BIO_PILE * 100,
               BIO_SE = sqrt(bVar) / BIO * 100,
               CARB_DUFF_SE = sqrt(cdVar) / CARB_DUFF * 100,
               CARB_LITTER_SE = sqrt(clVar) / CARB_LITTER * 100,
               CARB_1HR_SE = sqrt(csmVar) / CARB_1HR * 100,
               CARB_10HR_SE = sqrt(cmdVar) / CARB_10HR * 100,
               CARB_100HR_SE = sqrt(clgVar) / CARB_100HR * 100,
               CARB_1000HR_SE = sqrt(ccVar) / CARB_1000HR * 100,
               CARB_PILE_SE = sqrt(cpVar) / CARB_PILE * 100,
               CARB_SE = sqrt(cVar) / CARB * 100,
               # Per Acre variances
               vsmVar = (1/AREA_TOTAL^2) * (vsmVar + (VOL_1HR_ACRE^2 * aVar - 2 * VOL_1HR_ACRE * cvEst_vsm)),
               vmdVar = (1/AREA_TOTAL^2) * (vmdVar + (VOL_10HR_ACRE^2 * aVar - 2 * VOL_10HR_ACRE * cvEst_vmd)),
               vlgVar = (1/AREA_TOTAL^2) * (vlgVar + (VOL_100HR_ACRE^2 * aVar - 2 * VOL_100HR_ACRE * cvEst_vlg)),
               vcVar = (1/AREA_TOTAL^2) * (vcVar + (VOL_1000HR_ACRE^2 * aVar - 2 * VOL_1000HR_ACRE *cvEst_vc)),
               vpVar = (1/AREA_TOTAL^2) * (vpVar + (VOL_PILE_ACRE^2 * aVar - 2 * VOL_PILE_ACRE * cvEst_vp)),
               vVar = (1/AREA_TOTAL^2) * (vVar + (VOL_ACRE^2 * aVar - 2 * VOL_ACRE * cvEst_v)),
               bdVar = (1/AREA_TOTAL^2) * (bdVar + (BIO_DUFF_ACRE^2 * aVar - 2 * BIO_DUFF_ACRE * cvEst_bd)),
               blVar = (1/AREA_TOTAL^2) * (blVar + (BIO_LITTER_ACRE^2 * aVar - 2 * BIO_LITTER_ACRE * cvEst_bl)),
               bsmVar = (1/AREA_TOTAL^2) * (bsmVar + (BIO_1HR_ACRE^2 * aVar - 2 * BIO_1HR_ACRE * cvEst_bsm)),
               bmdVar = (1/AREA_TOTAL^2) * (bmdVar + (BIO_10HR_ACRE^2 * aVar - 2 * BIO_10HR_ACRE * cvEst_bmd)),
               blgVar = (1/AREA_TOTAL^2) * (blgVar + (BIO_100HR_ACRE^2 * aVar - 2 * BIO_100HR_ACRE * cvEst_blg)),
               bcVar = (1/AREA_TOTAL^2) * (bcVar + (BIO_1000HR_ACRE^2 * aVar - 2 * BIO_1000HR_ACRE * cvEst_bc)),
               bpVar = (1/AREA_TOTAL^2) * (bpVar + (BIO_PILE_ACRE^2 * aVar - 2 * BIO_PILE_ACRE * cvEst_bp)),
               bVar = (1/AREA_TOTAL^2) * (bVar + (BIO_ACRE^2 * aVar - 2 * BIO_ACRE * cvEst_b)),
               cdVar = (1/AREA_TOTAL^2) * (cdVar + (CARB_DUFF_ACRE^2 * aVar - 2 * CARB_DUFF_ACRE * cvEst_cd)),
               clVar = (1/AREA_TOTAL^2) * (clVar + (CARB_LITTER_ACRE^2 * aVar - 2 * CARB_LITTER_ACRE * cvEst_cl)),
               csmVar = (1/AREA_TOTAL^2) * (csmVar + (CARB_1HR_ACRE^2 * aVar - 2 * CARB_1HR_ACRE * cvEst_csm)),
               cmdVar = (1/AREA_TOTAL^2) * (cmdVar + (CARB_10HR_ACRE^2 * aVar - 2 * CARB_10HR_ACRE * cvEst_cmd)),
               clgVar = (1/AREA_TOTAL^2) * (clgVar + (CARB_100HR_ACRE^2 * aVar - 2 * CARB_100HR_ACRE * cvEst_clg)),
               ccVar = (1/AREA_TOTAL^2) * (ccVar + (CARB_1000HR_ACRE^2 * aVar - 2 * CARB_1000HR_ACRE * cvEst_cc)),
               cpVar = (1/AREA_TOTAL^2) * (cpVar + (CARB_PILE_ACRE^2 * aVar - 2 * CARB_PILE_ACRE * cvEst_cp)),
               cVar = (1/AREA_TOTAL^2) * (cVar + (CARB_ACRE^2 * aVar - 2 * CARB_ACRE * cvEst_c)),
               # Per acre sampling errors
               VOL_DUFF_ACRE_SE = NA,
               VOL_LITTER_ACRE_SE = NA,
               VOL_1HR_ACRE_SE = sqrt(vsmVar) / VOL_1HR_ACRE * 100,
               VOL_10HR_ACRE_SE = sqrt(vmdVar) / VOL_10HR_ACRE * 100,
               VOL_100HR_ACRE_SE = sqrt(vlgVar) / VOL_100HR_ACRE * 100,
               VOL_1000HR_ACRE_SE = sqrt(vcVar) / VOL_1000HR_ACRE * 100,
               VOL_PILE_ACRE_SE = sqrt(vpVar) / VOL_PILE_ACRE * 100,
               VOL_ACRE_SE = sqrt(vVar) / VOL_ACRE * 100,
               BIO_DUFF_ACRE_SE = sqrt(bdVar) / BIO_DUFF_ACRE * 100,
               BIO_LITTER_ACRE_SE = sqrt(blVar) / BIO_LITTER_ACRE * 100,
               BIO_1HR_ACRE_SE = sqrt(bsmVar) / BIO_1HR_ACRE * 100,
               BIO_10HR_ACRE_SE = sqrt(bmdVar) / BIO_10HR_ACRE * 100,
               BIO_100HR_ACRE_SE = sqrt(blgVar) / BIO_100HR_ACRE * 100,
               BIO_1000HR_ACRE_SE = sqrt(bcVar) / BIO_1000HR_ACRE * 100,
               BIO_PILE_ACRE_SE = sqrt(bpVar) / BIO_PILE_ACRE * 100,
               BIO_ACRE_SE = sqrt(bVar) / BIO_ACRE * 100,
               CARB_DUFF_ACRE_SE = sqrt(cdVar) / CARB_DUFF_ACRE * 100,
               CARB_LITTER_ACRE_SE = sqrt(clVar) / CARB_LITTER_ACRE * 100,
               CARB_1HR_ACRE_SE = sqrt(csmVar) / CARB_1HR_ACRE * 100,
               CARB_10HR_ACRE_SE = sqrt(cmdVar) / CARB_10HR_ACRE * 100,
               CARB_100HR_ACRE_SE = sqrt(clgVar) / CARB_100HR_ACRE * 100,
               CARB_1000HR_ACRE_SE = sqrt(ccVar) / CARB_1000HR_ACRE * 100,
               CARB_PILE_ACRE_SE = sqrt(cpVar) / CARB_PILE_ACRE * 100,
               CARB_ACRE_SE = sqrt(cVar) / CARB_ACRE * 100,
               nPlots_DWM = plotIn)
    })
    # Remove the total values if told to do so
    if (totals) {
      tOut <- tOut %>%
        select(grpBy, names(tOut)[str_detect(names(tOut), 'Var', negate = TRUE) &
                                    str_detect(names(tOut), 'cvEst', negate = TRUE) &
                                    str_detect(names(tOut), 'Est', negate = TRUE)], nPlots_DWM)
    } else {
      tOut <- tOut %>%
        select(grpBy, names(tOut)[str_detect(names(tOut), 'Var', negate = TRUE) &
                                    str_detect(names(tOut), 'cvEst', negate = TRUE) &
                                    str_detect(names(tOut), 'Est', negate = TRUE) &
                                    str_detect(names(tOut), 'ACRE')], nPlots_DWM)
    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

    ## Tidy things up if they didn't specify polys, returnSpatial
    if (tidy & is.null(polys) & returnSpatial == FALSE){
      ## pivot longer
      bio <- pivot_longer(select(tOut, grpBy, BIO_DUFF_ACRE:BIO_ACRE, nPlots_DWM), names_to = 'FUEL_TYPE', values_to = 'BIO_ACRE', cols = BIO_DUFF_ACRE:BIO_ACRE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      bio_SE <-pivot_longer(select(tOut, grpBy, BIO_DUFF_ACRE_SE:BIO_ACRE_SE), names_to = 'FUEL_TYPE', values_to = 'BIO_ACRE_SE', cols = BIO_DUFF_ACRE_SE:BIO_ACRE_SE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      vol <- pivot_longer(select(tOut, grpBy, VOL_DUFF_ACRE:VOL_ACRE), names_to = 'FUEL_TYPE', values_to = 'VOL_ACRE', cols = VOL_DUFF_ACRE:VOL_ACRE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      vol_SE <-pivot_longer(select(tOut, grpBy, VOL_DUFF_ACRE_SE:VOL_ACRE_SE), names_to = 'FUEL_TYPE', values_to = 'VOL_ACRE_SE', cols = VOL_DUFF_ACRE_SE:VOL_ACRE_SE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      carb <- pivot_longer(select(tOut, grpBy, CARB_DUFF_ACRE:CARB_ACRE), names_to = 'FUEL_TYPE', values_to = 'CARB_ACRE', cols = CARB_DUFF_ACRE:CARB_ACRE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      carb_SE <-pivot_longer(select(tOut, grpBy, CARB_DUFF_ACRE_SE:CARB_ACRE_SE), names_to = 'FUEL_TYPE', values_to = 'CARB_ACRE_SE', cols = CARB_DUFF_ACRE_SE:CARB_ACRE_SE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      ## rejoin
      fuel <- left_join(bio, bio_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
        left_join(vol, by = c(grpBy, 'FUEL_TYPE')) %>%
        left_join(vol_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
        left_join(carb, by = c(grpBy, 'FUEL_TYPE')) %>%
        left_join(carb_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
        select(grpBy, FUEL_TYPE, VOL_ACRE, BIO_ACRE, CARB_ACRE, VOL_ACRE_SE, BIO_ACRE_SE, CARB_ACRE_SE, nPlots_DWM) %>%
        filter(FUEL_TYPE %in% 'ACRE' == FALSE)

      if (totals){
        ## pivot longer
        bio <- pivot_longer(select(tOut, grpBy, BIO_DUFF:BIO), names_to = 'FUEL_TYPE', values_to = 'BIO_TOTAL', cols = BIO_DUFF:BIO) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        bio_SE <-pivot_longer(select(tOut, grpBy, BIO_DUFF_SE:BIO_SE), names_to = 'FUEL_TYPE', values_to = 'BIO_TOTAL_SE', cols = BIO_DUFF_SE:BIO_SE) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        vol <- pivot_longer(select(tOut, grpBy, VOL_DUFF:VOL), names_to = 'FUEL_TYPE', values_to = 'VOL_TOTAL', cols = VOL_DUFF:VOL) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        vol_SE <-pivot_longer(select(tOut, grpBy, VOL_DUFF_SE:VOL_SE), names_to = 'FUEL_TYPE', values_to = 'VOL_TOTAL_SE', cols = VOL_DUFF_SE:VOL_SE) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        carb <- pivot_longer(select(tOut, grpBy, CARB_DUFF:CARB), names_to = 'FUEL_TYPE', values_to = 'CARB_TOTAL', cols = CARB_DUFF:CARB) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        carb_SE <-pivot_longer(select(tOut, grpBy, CARB_DUFF_SE:CARB_SE), names_to = 'FUEL_TYPE', values_to = 'CARB_TOTAL_SE', cols = CARB_DUFF_SE:CARB_SE) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        ## Rejoin
        fuel <- fuel %>%
          left_join(bio, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(bio_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(vol, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(vol_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(carb, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(carb_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(select(tOut, AREA_TOTAL, AREA_TOTAL_SE, grpBy), by = c(grpBy))%>%
          select(grpBy, FUEL_TYPE, VOL_ACRE, BIO_ACRE, CARB_ACRE, VOL_TOTAL, BIO_TOTAL, CARB_TOTAL,
                 AREA_TOTAL, VOL_ACRE_SE, BIO_ACRE_SE, CARB_ACRE_SE,
                 VOL_TOTAL_SE, BIO_TOTAL_SE, CARB_TOTAL_SE,
                  AREA_TOTAL_SE, nPlots_DWM) %>%
          filter(FUEL_TYPE %in% 'ACRE' == FALSE)
      }
      tOut <- fuel

    }
  } # End byPlot
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

    suppressMessages({suppressWarnings({tOut <- left_join(tOut, polys, by = 'polyID') %>%
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




dwm_backup <- function(db,
                grpBy = NULL,
                polys = NULL,
                returnSpatial = FALSE,
                landType = 'forest',
                method = 'TI',
                lambda = .5,
                areaDomain = NULL,
                byPlot = FALSE,
                totals = FALSE,
                tidy = TRUE,
                nCores = 1) {

  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))
  db$COND_DWM_CALC <- db[['COND_DWM_CALC']] %>% mutate(DWM_CN = CN)
  db$COND <- db[['COND']] %>% mutate(CND_CN = CN)

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
        select(!!grpBy_quo)
    )

    # If column doesnt exist, just returns 0, not a dataframe
    if (is.null(nrow(d_quo))){
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT or COND tables. Did you accidentally quote the variables names? e.g. use grpBy = ECOSUBCD (correct) instead of grpBy = "ECOSUBCD". ', collapse = ', '))
    } else {
      # Convert to character
      grpBy <- names(d_quo)
    }
  }

  reqTables <- c('PLOT', 'COND_DWM_CALC', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
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

  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}
  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy

  ## IF the object was clipped
  if ('prev' %in% names(db$PLOT)){
    ## Only want the current plots, no grm
    db$PLOT <- filter(db$PLOT, prev == 0)
  }

  ### DEAL WITH TEXAS
  if (any(db$POP_EVAL$STATECD %in% 48)){
    ## Will require manual updates, fix your shit texas
    txIDS <- db$POP_EVAL %>%
      filter(STATECD %in% 48) %>%
      filter(END_INVYR < 2017) %>%
      filter(END_INVYR > 2006) %>%
      ## Removing any inventory that references east or west, sorry
      filter(str_detect(str_to_upper(EVAL_DESCR), 'EAST', negate = TRUE) &
               str_detect(str_to_upper(EVAL_DESCR), 'WEST', negate = TRUE))
    db$POP_EVAL <- bind_rows(filter(db$POP_EVAL, !(STATECD %in% 48)), txIDS)
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
    #grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
    grpBy = c(grpBy, 'polyID')

    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('LON', 'LAT', pltID)) %>%
      filter(!is.na(LAT) & !is.na(LON)) %>%
      distinct(pltID, .keep_all = TRUE)
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys))

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
      left_join(select(pltSF, polyID, pltID), by = 'pltID')

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
    filter(EVAL_TYP == 'EXPDWM') %>%
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
      out <- parLapply(cl, X = names(plts), fun = dwmHelper1, plts, db, grpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = dwmHelper1, plts, db, grpBy, byPlot, mc.cores = nCores)
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
        out <- parLapply(cl, X = names(popState), fun = dwmHelper2, popState, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = dwmHelper2, popState, t, grpBy, method, mc.cores = nCores)
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
      tEst <- tEst %>%
        mutate_at(vars(aEst:cEst), ~(.*wgt)) %>%
        mutate_at(vars(aVar:cvEst_c), ~(.*(wgt^2))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize_at(vars(aEst:cvEst_c), sum, na.rm = TRUE)
    }

    ##---------------------  TOTALS and RATIOS
    suppressWarnings({
      ## Bring them together
      tOut <- ungroup(tEst) %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        # Renaming, computing ratios, and SE
        mutate(VOL_DUFF = NA,
               VOL_LITTER = NA,VOL_1HR = vsmEst,
               VOL_10HR = vmdEst,
               VOL_100HR = vlgEst,
               VOL_1000HR =  vcEst,
               VOL_PILE =  vpEst,
               VOL =  vEst,
               BIO_DUFF = bdEst,
               BIO_LITTER = blEst,
               BIO_1HR = bsmEst,
               BIO_10HR = bmdEst,
               BIO_100HR = blgEst,
               BIO_1000HR = bcEst,
               BIO_PILE = bpEst,
               BIO = bEst,
               CARB_DUFF = cdEst,
               CARB_LITTER = clEst,
               CARB_1HR = csmEst,
               CARB_10HR = cmdEst,
               CARB_100HR = clgEst,
               CARB_1000HR = ccEst,
               CARB_PILE = cpEst,
               CARB = cEst,
               AREA_TOTAL = aEst,
               ## RATIOS
               VOL_DUFF_ACRE = NA,
               VOL_LITTER_ACRE = NA,
               VOL_1HR_ACRE = vsmEst / aEst,
               VOL_10HR_ACRE = vmdEst / aEst,
               VOL_100HR_ACRE = vlgEst / aEst,
               VOL_1000HR_ACRE =  vcEst / aEst,
               VOL_PILE_ACRE =  vpEst / aEst,
               VOL_ACRE =  vEst / aEst,
               BIO_DUFF_ACRE = bdEst / aEst,
               BIO_LITTER_ACRE = blEst / aEst,
               BIO_1HR_ACRE = bsmEst / aEst,
               BIO_10HR_ACRE = bmdEst / aEst,
               BIO_100HR_ACRE = blgEst / aEst,
               BIO_1000HR_ACRE = bcEst / aEst,
               BIO_PILE_ACRE = bpEst / aEst,
               BIO_ACRE = bEst / aEst,
               CARB_DUFF_ACRE = cdEst / aEst,
               CARB_LITTER_ACRE = clEst / aEst,
               CARB_1HR_ACRE = csmEst / aEst,
               CARB_10HR_ACRE = cmdEst / aEst,
               CARB_100HR_ACRE = clgEst / aEst,
               CARB_1000HR_ACRE = ccEst / aEst,
               CARB_PILE_ACRE = cpEst / aEst,
               CARB_ACRE = cEst / aEst,
               # Sampling Errors totals
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
               VOL_DUFF_SE = NA,
               VOL_LITTER_SE = NA,
               VOL_1HR_SE = sqrt(vsmVar) / VOL_1HR * 100,
               VOL_10HR_SE = sqrt(vmdVar) / VOL_10HR * 100,
               VOL_100HR_SE = sqrt(vlgVar) / VOL_100HR * 100,
               VOL_1000HR_SE = sqrt(vcVar) / VOL_1000HR * 100,
               VOL_PILE_SE = sqrt(vpVar) / VOL_PILE * 100,
               VOL_SE = sqrt(vVar) / VOL * 100,
               BIO_DUFF_SE = sqrt(bdVar) / BIO_DUFF * 100,
               BIO_LITTER_SE = sqrt(blVar) / BIO_LITTER * 100,
               BIO_1HR_SE = sqrt(bsmVar) / BIO_1HR * 100,
               BIO_10HR_SE = sqrt(bmdVar) / BIO_10HR * 100,
               BIO_100HR_SE = sqrt(blgVar) / BIO_100HR * 100,
               BIO_1000HR_SE = sqrt(bcVar) / BIO_1000HR * 100,
               BIO_PILE_SE = sqrt(bpVar) / BIO_PILE * 100,
               BIO_SE = sqrt(bVar) / BIO * 100,
               CARB_DUFF_SE = sqrt(cdVar) / CARB_DUFF * 100,
               CARB_LITTER_SE = sqrt(clVar) / CARB_LITTER * 100,
               CARB_1HR_SE = sqrt(csmVar) / CARB_1HR * 100,
               CARB_10HR_SE = sqrt(cmdVar) / CARB_10HR * 100,
               CARB_100HR_SE = sqrt(clgVar) / CARB_100HR * 100,
               CARB_1000HR_SE = sqrt(ccVar) / CARB_1000HR * 100,
               CARB_PILE_SE = sqrt(cpVar) / CARB_PILE * 100,
               CARB_SE = sqrt(cVar) / CARB * 100,
               # Per Acre variances
               vsmVar = (1/AREA_TOTAL^2) * (vsmVar + (VOL_1HR_ACRE^2 * aVar - 2 * VOL_1HR_ACRE * cvEst_vsm)),
               vmdVar = (1/AREA_TOTAL^2) * (vmdVar + (VOL_10HR_ACRE^2 * aVar - 2 * VOL_10HR_ACRE * cvEst_vmd)),
               vlgVar = (1/AREA_TOTAL^2) * (vlgVar + (VOL_100HR_ACRE^2 * aVar - 2 * VOL_100HR_ACRE * cvEst_vlg)),
               vcVar = (1/AREA_TOTAL^2) * (vcVar + (VOL_1000HR_ACRE^2 * aVar - 2 * VOL_1000HR_ACRE *cvEst_vc)),
               vpVar = (1/AREA_TOTAL^2) * (vpVar + (VOL_PILE_ACRE^2 * aVar - 2 * VOL_PILE_ACRE * cvEst_vp)),
               vVar = (1/AREA_TOTAL^2) * (vVar + (VOL_ACRE^2 * aVar - 2 * VOL_ACRE * cvEst_v)),
               bdVar = (1/AREA_TOTAL^2) * (bdVar + (BIO_DUFF_ACRE^2 * aVar - 2 * BIO_DUFF_ACRE * cvEst_bd)),
               blVar = (1/AREA_TOTAL^2) * (blVar + (BIO_LITTER_ACRE^2 * aVar - 2 * BIO_LITTER_ACRE * cvEst_bl)),
               bsmVar = (1/AREA_TOTAL^2) * (bsmVar + (BIO_1HR_ACRE^2 * aVar - 2 * BIO_1HR_ACRE * cvEst_bsm)),
               bmdVar = (1/AREA_TOTAL^2) * (bmdVar + (BIO_10HR_ACRE^2 * aVar - 2 * BIO_10HR_ACRE * cvEst_bmd)),
               blgVar = (1/AREA_TOTAL^2) * (blgVar + (BIO_100HR_ACRE^2 * aVar - 2 * BIO_100HR_ACRE * cvEst_blg)),
               bcVar = (1/AREA_TOTAL^2) * (bcVar + (BIO_1000HR_ACRE^2 * aVar - 2 * BIO_1000HR_ACRE * cvEst_bc)),
               bpVar = (1/AREA_TOTAL^2) * (bpVar + (BIO_PILE_ACRE^2 * aVar - 2 * BIO_PILE_ACRE * cvEst_bp)),
               bVar = (1/AREA_TOTAL^2) * (bVar + (BIO_ACRE^2 * aVar - 2 * BIO_ACRE * cvEst_b)),
               cdVar = (1/AREA_TOTAL^2) * (cdVar + (CARB_DUFF_ACRE^2 * aVar - 2 * CARB_DUFF_ACRE * cvEst_cd)),
               clVar = (1/AREA_TOTAL^2) * (clVar + (CARB_LITTER_ACRE^2 * aVar - 2 * CARB_LITTER_ACRE * cvEst_cl)),
               csmVar = (1/AREA_TOTAL^2) * (csmVar + (CARB_1HR_ACRE^2 * aVar - 2 * CARB_1HR_ACRE * cvEst_csm)),
               cmdVar = (1/AREA_TOTAL^2) * (cmdVar + (CARB_10HR_ACRE^2 * aVar - 2 * CARB_10HR_ACRE * cvEst_cmd)),
               clgVar = (1/AREA_TOTAL^2) * (clgVar + (CARB_100HR_ACRE^2 * aVar - 2 * CARB_100HR_ACRE * cvEst_clg)),
               ccVar = (1/AREA_TOTAL^2) * (ccVar + (CARB_1000HR_ACRE^2 * aVar - 2 * CARB_1000HR_ACRE * cvEst_cc)),
               cpVar = (1/AREA_TOTAL^2) * (cpVar + (CARB_PILE_ACRE^2 * aVar - 2 * CARB_PILE_ACRE * cvEst_cp)),
               cVar = (1/AREA_TOTAL^2) * (cVar + (CARB_ACRE^2 * aVar - 2 * CARB_ACRE * cvEst_c)),
               # Per acre sampling errors
               VOL_DUFF_ACRE_SE = NA,
               VOL_LITTER_ACRE_SE = NA,
               VOL_1HR_ACRE_SE = sqrt(vsmVar) / VOL_1HR_ACRE * 100,
               VOL_10HR_ACRE_SE = sqrt(vmdVar) / VOL_10HR_ACRE * 100,
               VOL_100HR_ACRE_SE = sqrt(vlgVar) / VOL_100HR_ACRE * 100,
               VOL_1000HR_ACRE_SE = sqrt(vcVar) / VOL_1000HR_ACRE * 100,
               VOL_PILE_ACRE_SE = sqrt(vpVar) / VOL_PILE_ACRE * 100,
               VOL_ACRE_SE = sqrt(vVar) / VOL_ACRE * 100,
               BIO_DUFF_ACRE_SE = sqrt(bdVar) / BIO_DUFF_ACRE * 100,
               BIO_LITTER_ACRE_SE = sqrt(blVar) / BIO_LITTER_ACRE * 100,
               BIO_1HR_ACRE_SE = sqrt(bsmVar) / BIO_1HR_ACRE * 100,
               BIO_10HR_ACRE_SE = sqrt(bmdVar) / BIO_10HR_ACRE * 100,
               BIO_100HR_ACRE_SE = sqrt(blgVar) / BIO_100HR_ACRE * 100,
               BIO_1000HR_ACRE_SE = sqrt(bcVar) / BIO_1000HR_ACRE * 100,
               BIO_PILE_ACRE_SE = sqrt(bpVar) / BIO_PILE_ACRE * 100,
               BIO_ACRE_SE = sqrt(bVar) / BIO_ACRE * 100,
               CARB_DUFF_ACRE_SE = sqrt(cdVar) / CARB_DUFF_ACRE * 100,
               CARB_LITTER_ACRE_SE = sqrt(clVar) / CARB_LITTER_ACRE * 100,
               CARB_1HR_ACRE_SE = sqrt(csmVar) / CARB_1HR_ACRE * 100,
               CARB_10HR_ACRE_SE = sqrt(cmdVar) / CARB_10HR_ACRE * 100,
               CARB_100HR_ACRE_SE = sqrt(clgVar) / CARB_100HR_ACRE * 100,
               CARB_1000HR_ACRE_SE = sqrt(ccVar) / CARB_1000HR_ACRE * 100,
               CARB_PILE_ACRE_SE = sqrt(cpVar) / CARB_PILE_ACRE * 100,
               CARB_ACRE_SE = sqrt(cVar) / CARB_ACRE * 100,
               nPlots_DWM = plotIn)
    })
    # Remove the total values if told to do so
    if (totals) {
      tOut <- tOut %>%
        select(grpBy, names(tOut)[str_detect(names(tOut), 'Var', negate = TRUE) &
                                    str_detect(names(tOut), 'cvEst', negate = TRUE) &
                                    str_detect(names(tOut), 'Est', negate = TRUE)], nPlots_DWM)
    } else {
      tOut <- tOut %>%
        select(grpBy, names(tOut)[str_detect(names(tOut), 'Var', negate = TRUE) &
                                    str_detect(names(tOut), 'cvEst', negate = TRUE) &
                                    str_detect(names(tOut), 'Est', negate = TRUE) &
                                    str_detect(names(tOut), 'ACRE')], nPlots_DWM)
    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

    ## Tidy things up if they didn't specify polys, returnSpatial
    if (tidy & is.null(polys) & returnSpatial == FALSE){
      ## pivot longer
      bio <- pivot_longer(select(tOut, grpBy, BIO_DUFF_ACRE:BIO_ACRE, nPlots_DWM), names_to = 'FUEL_TYPE', values_to = 'BIO_ACRE', cols = BIO_DUFF_ACRE:BIO_ACRE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      bio_SE <-pivot_longer(select(tOut, grpBy, BIO_DUFF_ACRE_SE:BIO_ACRE_SE), names_to = 'FUEL_TYPE', values_to = 'BIO_ACRE_SE', cols = BIO_DUFF_ACRE_SE:BIO_ACRE_SE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      vol <- pivot_longer(select(tOut, grpBy, VOL_DUFF_ACRE:VOL_ACRE), names_to = 'FUEL_TYPE', values_to = 'VOL_ACRE', cols = VOL_DUFF_ACRE:VOL_ACRE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      vol_SE <-pivot_longer(select(tOut, grpBy, VOL_DUFF_ACRE_SE:VOL_ACRE_SE), names_to = 'FUEL_TYPE', values_to = 'VOL_ACRE_SE', cols = VOL_DUFF_ACRE_SE:VOL_ACRE_SE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      carb <- pivot_longer(select(tOut, grpBy, CARB_DUFF_ACRE:CARB_ACRE), names_to = 'FUEL_TYPE', values_to = 'CARB_ACRE', cols = CARB_DUFF_ACRE:CARB_ACRE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      carb_SE <-pivot_longer(select(tOut, grpBy, CARB_DUFF_ACRE_SE:CARB_ACRE_SE), names_to = 'FUEL_TYPE', values_to = 'CARB_ACRE_SE', cols = CARB_DUFF_ACRE_SE:CARB_ACRE_SE) %>%
        mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
      ## rejoin
      fuel <- left_join(bio, bio_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
        left_join(vol, by = c(grpBy, 'FUEL_TYPE')) %>%
        left_join(vol_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
        left_join(carb, by = c(grpBy, 'FUEL_TYPE')) %>%
        left_join(carb_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
        select(grpBy, FUEL_TYPE, VOL_ACRE, BIO_ACRE, CARB_ACRE, VOL_ACRE_SE, BIO_ACRE_SE, CARB_ACRE_SE, nPlots_DWM) %>%
        filter(FUEL_TYPE %in% 'ACRE' == FALSE)

      if (totals){
        ## pivot longer
        bio <- pivot_longer(select(tOut, grpBy, BIO_DUFF:BIO), names_to = 'FUEL_TYPE', values_to = 'BIO_TOTAL', cols = BIO_DUFF:BIO) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        bio_SE <-pivot_longer(select(tOut, grpBy, BIO_DUFF_SE:BIO_SE), names_to = 'FUEL_TYPE', values_to = 'BIO_TOTAL_SE', cols = BIO_DUFF_SE:BIO_SE) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        vol <- pivot_longer(select(tOut, grpBy, VOL_DUFF:VOL), names_to = 'FUEL_TYPE', values_to = 'VOL_TOTAL', cols = VOL_DUFF:VOL) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        vol_SE <-pivot_longer(select(tOut, grpBy, VOL_DUFF_SE:VOL_SE), names_to = 'FUEL_TYPE', values_to = 'VOL_TOTAL_SE', cols = VOL_DUFF_SE:VOL_SE) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        carb <- pivot_longer(select(tOut, grpBy, CARB_DUFF:CARB), names_to = 'FUEL_TYPE', values_to = 'CARB_TOTAL', cols = CARB_DUFF:CARB) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        carb_SE <-pivot_longer(select(tOut, grpBy, CARB_DUFF_SE:CARB_SE), names_to = 'FUEL_TYPE', values_to = 'CARB_TOTAL_SE', cols = CARB_DUFF_SE:CARB_SE) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE,pattern= '_', simplify = TRUE,)[,2])
        ## Rejoin
        fuel <- fuel %>%
          left_join(bio, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(bio_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(vol, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(vol_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(carb, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(carb_SE, by = c(grpBy, 'FUEL_TYPE')) %>%
          left_join(select(tOut, AREA_TOTAL, AREA_TOTAL_SE, grpBy), by = c(grpBy))%>%
          select(grpBy, FUEL_TYPE, VOL_ACRE, BIO_ACRE, CARB_ACRE, VOL_TOTAL, BIO_TOTAL, CARB_TOTAL,
                 AREA_TOTAL, VOL_ACRE_SE, BIO_ACRE_SE, CARB_ACRE_SE,
                 VOL_TOTAL_SE, BIO_TOTAL_SE, CARB_TOTAL_SE,
                 AREA_TOTAL_SE, nPlots_DWM) %>%
          filter(FUEL_TYPE %in% 'ACRE' == FALSE)
      }
      tOut <- fuel

    }
  } # End byPlot
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

    suppressMessages({suppressWarnings({tOut <- left_join(tOut, polys, by = 'polyID') %>%
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


