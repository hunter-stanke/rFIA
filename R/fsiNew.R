fsiStarter <- function(x,
                       db,
                       grpBy_quo = NULL,
                       scaleBy_quo = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       bySpecies = FALSE,
                       bySizeClass = FALSE,
                       landType = 'forest',
                       treeType = 'live',
                       method = 'sma',
                       lambda = .5,
                       treeDomain = NULL,
                       areaDomain = NULL,
                       totals = FALSE,
                       byPlot = FALSE,
                       useLM = FALSE,
                       nCores = 1,
                       remote,
                       mr){

  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
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
  db$TREE$TRE_CN <- db$TREE$CN
  ##  don't have to change original code
  #grpBy_quo <- enquo(grpBy)

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
  } else {
    grpBy <- NULL
  }

  # Probably cheating, but it works
  if (quo_name(scaleBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[10,] %>% # Just the first row
        left_join(select(db$COND, PLT_CN, names(db$COND)[names(db$COND) %in% names(db$PLOT) == FALSE]), by = 'PLT_CN') %>%
        select(!!scaleBy_quo)
    )

    # If column doesnt exist, just returns 0, not a dataframe
    if (is.null(nrow(d_quo))){
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT or COND tables. Did you accidentally quote the variables names? e.g. use scaleBy = ECOSUBCD (correct) instead of scaleBy = "ECOSUBCD". ', collapse = ', '))
    } else {
      # Convert to character
      scaleBy <- names(d_quo)
    }
  } else {
    scaleBy <- NULL
  }


  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')

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
          library(sf)
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

  # Same as above for tree (ex. trees > 20 ft tall)
  #treeDomain <- substitute(treeDomain)
  tD <- rlang::eval_tidy(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$TREE$tD <- as.numeric(tD)



  ## We only want inventory/ population info from t2 plots, but we need the plot tree cond data
  ## for t1 and t2
  remPlts <- db$PLOT %>%
    select(PLT_CN, PREV_PLT_CN, DESIGNCD, REMPER, PLOT_STATUS_CD) %>%
    ## Has to have a remeasurement, be in the current sample, and of the national design
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & PLOT_STATUS_CD != 3 & DESIGNCD == 1) %>%
    left_join(select(db$PLOT, PLT_CN, DESIGNCD, PLOT_STATUS_CD), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    ## past emasurement must be in the previous sample and of national design
    filter(PLOT_STATUS_CD.prev != 3 & DESIGNCD.prev == 1)

  ### Snag the EVALIDs that are needed
  db$POP_EVAL<- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPVOL') %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)

  ## If a most-recent subset, make sure that we don't get two reporting years in
  ## western states
  if (mr) {
    db$POP_EVAL <- db$POP_EVAL %>%
      group_by(EVAL_TYP) %>%
      filter(END_INVYR == max(END_INVYR, na.rm = TRUE))
  }

  ## Make an annual panel ID, associated with an INVYR

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'END_INVYR', 'EVAL_TYP')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
    ## ONLY REMEASURED PLOTS MEETING CRITERIA ABOVE
    filter(PLT_CN %in% remPlts$PLT_CN) %>%
    ungroup() %>%
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
      ungroup() %>%
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

  ## Only the necessary plots for EVAL of interest
  remPltList <- unique(c(remPlts$PLT_CN, remPlts$PREV_PLT_CN))
  db$PLOT <- filter(db$PLOT, PLT_CN %in% remPltList)
  db$COND <- filter(db$COND, PLT_CN %in% remPltList)
  db$TREE <- filter(db$TREE, PLT_CN %in% remPltList)

  ## Tree basal area per acre
  db$TREE <- db$TREE %>%
    mutate(BAA = basalArea(DIA) * TPA_UNADJ)

  ## Narrow up the tables to the necessary variables, reduces memory load
  ## sent out the cores
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% c(grpBy, scaleBy)]
  grpC <- names(db$COND)[names(db$COND) %in% c(grpBy, scaleBy) & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% c(grpBy, scaleBy) & names(db$TREE) %in% c(grpP, grpC) == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  db$PLOT <- select(db$PLOT, c('PLT_CN', pltID, 'REMPER', 'DESIGNCD', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR',
                               'MEASYEAR', 'MEASMON', 'MEASDAY', 'PLOT_STATUS_CD', PREV_PLT_CN, grpP, 'aD_p', 'sp'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD'))
  db$TREE <- select(db$TREE, c('PLT_CN', 'TRE_CN', 'CONDID', 'DIA', 'TPA_UNADJ', 'BAA', 'SUBP', 'TREE', grpT, 'tD', 'typeD',
                               PREVCOND, PREV_TRE_CN, STATUSCD, SPCD))


  ## Merging state and county codes
  plts <- split(db$PLOT, as.factor(paste(db$PLOT$COUNTYCD, db$PLOT$STATECD, sep = '_')))


  ## Use the linear model procedures or strictly t2-t1 / remper?
  if (useLM){
    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        cl <- makeCluster(nCores)
        clusterEvalQ(cl, {
          library(dplyr)
          library(stringr)
          library(rFIA)
          library(tidyr)
          library(purrr)
        })
        out <- parLapply(cl, X = names(plts), fun = fsiHelper1_lm, plts, db, grpBy, scaleBy, byPlot)
        #stopCluster(cl) # Keep the cluster active for the next run
      } else { # Unix systems
        out <- mclapply(names(plts), FUN = fsiHelper1_lm, plts, db, grpBy, scaleBy, byPlot, mc.cores = nCores)
      }
    })
  } else {
    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        cl <- makeCluster(nCores)
        clusterEvalQ(cl, {
          library(dplyr)
          library(stringr)
          library(rFIA)
          library(tidyr)
        })
        out <- parLapply(cl, X = names(plts), fun = fsiHelper1, plts, db, grpBy, scaleBy, byPlot)
        #stopCluster(cl) # Keep the cluster active for the next run
      } else { # Unix systems
        out <- mclapply(names(plts), FUN = fsiHelper1, plts, db, grpBy, scaleBy, byPlot, mc.cores = nCores)
      }
    })
  }

  ## back to dataframes
  out <- unlist(out, recursive = FALSE)
  t <- bind_rows(out[names(out) == 't'])
  t1 <- bind_rows(out[names(out) == 't1'])
  a <- bind_rows(out[names(out) == 'a'])

  out <- list(t = t, t1 = t1, a = a, grpBy = grpBy, scaleBy = scaleBy, grpByOrig = grpByOrig, pops = pops)

}




#' @export
fsi <- function(db,
                   grpBy = NULL,
                   polys = NULL,
                   returnSpatial = FALSE,
                   bySpecies = FALSE,
                   bySizeClass = FALSE,
                   landType = 'forest',
                   treeType = 'live',
                   method = 'sma',
                   lambda = .5,
                   treeDomain = NULL,
                   areaDomain = NULL,
                   totals = TRUE,
                   byPlot = FALSE,
                   useLM = FALSE,
                   scaleBy = NULL,
                   nCores = 1) {

  ##  don't have to change original code
  grpBy_quo <- rlang::enquo(grpBy)
  scaleBy_quo <- rlang::enquo(scaleBy)
  areaDomain <- rlang::enquo(areaDomain)
  treeDomain <- rlang::enquo(treeDomain)

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

  ## Check for a most recent subset
  if (remote){
    if ('mostRecent' %in% names(db)){
      mr = db$mostRecent # logical
    } else {
      mr = FALSE
    }
    ## In-memory
  } else {
    if ('mostRecent' %in% names(db)){
      mr = TRUE
    } else {
      mr = FALSE
    }
  }

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf') %>%
      mutate_if(is.factor,
                as.character)
    ## A unique ID
    polys$polyID <- 1:nrow(polys)
  }


  ## Run the main portion
  out <- lapply(X = iter, FUN = fsiStarter, db,
                grpBy_quo = grpBy_quo, scaleBy_quo, polys, returnSpatial,
                bySpecies, bySizeClass,
                landType, treeType, method,
                lambda, treeDomain, areaDomain,
                totals, byPlot, useLM,
                nCores, remote, mr)

  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  t <- bind_rows(out[names(out) == 't'])
  t1 <- bind_rows(out[names(out) == 't1'])
  a <- bind_rows(out[names(out) == 'a'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  scaleBy <- out[names(out) == 'scaleBy'][[1]]
  grpByOrig <- out[names(out) == 'grpByOrig'][[1]]

  ## Have to update the PLT_CNs if useLM = TRUE
  if (useLM){
    t1 <- t1 %>%
      ungroup() %>%
      select(-c(PLT_CN)) %>%
      left_join(distinct(select(t, PLT_CN, pltID)), by = 'pltID')
    t$PREV_PLT_CN = NA
  }

  ## Prep the data for modeling the upper boundary
  ## of the size-density curve
  scaleSyms <- syms(scaleBy)

  ## Get groups prepped to fit model
  grpRates <- t1 %>%
    ungroup() %>%
    filter(TPA1 > 0) %>%
    mutate(t = log(TPA1),
           b = log(BA1)) %>%
    select(t, b, PLT_CN, !!!scaleSyms)

  if (!is.null(scaleBy)){
    ## group IDS
    grpRates <- mutate(grpRates, grps = paste(!!!scaleSyms))
  } else {
    grpRates$grps = 1
  }

  ## If more than one group use a mixed model
  if (length(unique(grpRates$grps)) > 1){

    ## Run lmm at the 99 percentile of the distribution
    mod <- lqmm(t ~ b, random = ~ b, group = grps,
                tau = .99, data = grpRates,
                control = list(method = "df"))

    suppressWarnings({
      ## Summarize results
      beta1 <- lqmm::coef.lqmm(mod)[1] + lqmm::ranef.lqmm(mod)[1]
      beta2 <- lqmm::coef.lqmm(mod)[2] + lqmm::ranef.lqmm(mod)[2]
      betas <- bind_cols(beta1, beta2) %>%
        mutate(grps = row.names(.))
      names(betas) <- c('int', 'rate', 'grps')
    })


    ## Rejoin and estimate slopes at each plot
    grpRates <- grpRates %>%
      left_join(betas, by = 'grps') %>%
      mutate_at(.vars = vars(t, b, int), .funs = exp) %>%
      mutate(slope = (int * rate) * (b^(rate-1))) %>%
      select(PLT_CN, slope, int, rate, !!!scaleSyms)

  } else {

    ## Run lmm at the 99 percentile of the distribution
    mod <- lqm(t ~ b,
                tau = .99, data = grpRates)

    ## Summarize results
    beta1 <- coef(mod)[1]
    beta2 <- coef(mod)[2]
    betas <- data.frame(beta1, beta2) %>%
      mutate(grps = 1)
    names(betas) <- c('int', 'rate', 'grps')

    ## Rejoin and estimate slopes at each plot
    grpRates <- grpRates %>%
      left_join(betas, by = 'grps') %>%
      mutate_at(.vars = vars(t, b, int), .funs = exp) %>%
      mutate(slope = (int * rate) * (b^(rate-1))) %>%
      select(PLT_CN, slope, int, rate, !!!scaleSyms)
  }



  if (byPlot){

    ## back to dataframes
    tOut <- t

    ## Scale the indices by group
    tOut <- tOut %>%
      left_join(grpRates, by = c('PLT_CN', scaleBy)) %>%
      ## When PREV_BA is 0, slope is infinite
      mutate(slope = case_when(PREV_BA == 0 ~ 100000,
                               TRUE ~ slope)) %>%
      mutate(dt = CHNG_TPA,
             db = CHNG_BA,
             t1 = PREV_TPA / REMPER,
             t2 = CURR_TPA / REMPER,
             b1 = PREV_BA / REMPER,
             b2 = CURR_BA / REMPER) %>%
      ## The FSI and % FSI
      mutate(FSI = projectPoints(db, dt, -(1/slope), 0, returnPoint = FALSE),
             ## can and will produce INF here due to division by zero, that's fine, just use the FSI if that matters to you
             FSI2 = projectPoints(b2, t2,-(1/slope), 0, returnPoint = FALSE),
             FSI1 = projectPoints(b1, t1, -(1/slope), 0, returnPoint = FALSE)) %>%
      ## Summing across scaleBy
      group_by(.dots = grpBy[!c(grpBy %in% 'YEAR')], YEAR, PLT_CN, PLOT_STATUS_CD, PREV_PLT_CN,
               REMPER) %>%
      summarize(PREV_TPA = sum(PREV_TPA, na.rm = TRUE),
                PREV_BAA = sum(PREV_BAA, na.rm = TRUE),
                PREV_BA = sum(PREV_BA, na.rm = TRUE),
                CHNG_TPA = sum(CHNG_TPA, na.rm = TRUE),
                CHNG_BAA = sum(CHNG_BAA, na.rm = TRUE),
                CHNG_BA = sum(CHNG_BA, na.rm = TRUE),
                CURR_TPA = sum(CURR_TPA, na.rm = TRUE),
                CURR_BAA = sum(CURR_BAA, na.rm = TRUE),
                CURR_BA = sum(CURR_BA, na.rm = TRUE),
                FSI = mean(FSI, na.rm = TRUE),
                FSI2 = mean(FSI2, na.rm = TRUE),
                FSI1 = mean(FSI1, na.rm = TRUE)) %>%
      mutate(PERC_FSI = (FSI2 - FSI1) / FSI1 * 100)




    ## Make it spatial
    if (returnSpatial){
      tOut <- tOut %>%
        filter(!is.na(LAT) & !is.na(LON)) %>%
        st_as_sf(coords = c('LON', 'LAT'),
                 crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }


    tOut <- select(tOut, YEAR, PLT_CN, any_of('PREV_PLT_CN'), PLOT_STATUS_CD, grpBy[grpBy != 'YEAR'],
                   FSI, PERC_FSI, REMPER, everything()) %>%
      select(-c(FSI2, FSI1))



    ## Population estimation
  } else {
    ## back to dataframes
    pops <- bind_rows(out[names(out) == 'pops'])

    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)
    popState <- split(pops, as.factor(pops$STATECD))


    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        cl <- makeCluster(nCores)
        clusterEvalQ(cl, {
          library(dplyr)
          library(stringr)
          library(rFIA)
          library(tidyr)
        })
        out <- parLapply(cl, X = names(popState), fun = fsiHelper2, popState, t, a, grpBy, scaleBy, method, grpRates)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = fsiHelper2, popState, t, a, grpBy, scaleBy, method, grpRates, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])
    tEst <- ungroup(tEst)

    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
      ### ---- SIMPLE MOVING AVERAGE
      if (str_to_upper(method) == 'SMA'){
        ## Assuming a uniform weighting scheme
        wgts <- pops %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(wgt = 1 / length(unique(INVYR)))

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
            summarize(l = 1 - first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = 'ESTN_UNIT_CN') %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(ESTN_UNIT_CN, INVYR, wgt)
        } else {
          grpBy <- c('lambda', grpBy)
          ## Duplicate weights for each level of lambda
          yrWgts <- list()
          for (i in 1:length(unique(lambda))) {
            yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
          }
          wgts <- bind_rows(yrWgts)
          ## Want sum of weighitng functions
          neu <- wgts %>%
            group_by(lambda, ESTN_UNIT_CN) %>%
            summarize(l = 1 - first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = c('lambda', 'ESTN_UNIT_CN')) %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(lambda, ESTN_UNIT_CN, INVYR, wgt)
        }

        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))

      }

      ### Applying the weights
      tEst <- tEst %>%
        mutate_at(vars(ctEst:faEst), ~(.*wgt)) %>%
        mutate_at(vars(ctVar:cvEst_psi), ~(.*(wgt^2))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize_at(vars(ctEst:plotIn_t), sum, na.rm = TRUE)
    }

    suppressMessages({suppressWarnings({
      ## If a clip was specified, handle the reporting years
      if (mr){
        ## If a most recent subset, ignore differences in reporting years across states
        ## instead combine most recent information from each state
        # ID mr years by group
        maxyearsT <- tEst %>%
          select(grpBy) %>%
          group_by(.dots = grpBy[!c(grpBy %in% 'YEAR')]) %>%
          summarise(YEAR = max(YEAR, na.rm = TRUE))

        # Combine estimates
        tEst <- tEst %>%
          ungroup() %>%
          select(-c(YEAR)) %>%
          left_join(maxyearsT, by = grpBy[!c(grpBy %in% 'YEAR')])
      }
    })})


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
        mutate(TPA_RATE = ctEst / ptEst,
               BA_RATE = cbEst / pbEst,
               FSI = siEst / faEst,
               PERC_FSI = siEst / si1Est,

               ## Ratio variance
               ctVar = (1/ptEst^2) * (ctVar + (TPA_RATE^2 * ptVar) - (2 * TPA_RATE * cvEst_ct)),
               cbVar = (1/pbEst^2) * (cbVar + (BA_RATE^2 * pbVar) - (2 * BA_RATE * cvEst_cb)),
               psiVar = (1/si1Est^2) * (siVar + (PERC_FSI^2 * si1Var) - (2 * PERC_FSI * cvEst_psi)),
               siVar = (1/faEst^2) * (siVar + (FSI^2 * faVar) - (2 * FSI * cvEst_si)),

               ## Make it a percent
               PERC_FSI = PERC_FSI * 100,
               psiVar = psiVar * (100^2),

               ## RATIO SE
               TPA_RATE_SE = sqrt(ctVar) / abs(TPA_RATE) * 100,
               BA_RATE_SE = sqrt(cbVar) / abs(BA_RATE) * 100,
               FSI_SE = sqrt(siVar) / abs(FSI) * 100,
               FSI_VAR = siVar,
               PERC_FSI_SE = sqrt(psiVar) / abs(PERC_FSI) * 100,
               PERC_FSI_VAR = psiVar,
               TPA_RATE_VAR = ctVar,
               BA_RATE_VAR = cbVar,
               nPlots = plotIn_t,
               N = nh,
               FSI_INT = qt(.975, df=N-1) * (sqrt(siVar)/sqrt(N)),
               PERC_FSI_INT = qt(.975, df=N-1) * (sqrt(psiVar)/sqrt(N))) %>%
        mutate(FSI_STATUS = case_when(
          FSI < 0 & FSI + FSI_INT < 0 ~ 'Decline',
          FSI < 0 & FSI + FSI_INT > 0 ~ 'Stable',
          FSI > 0 & FSI - FSI_INT > 0  ~ 'Expand',
          TRUE ~ 'Stable'
        ))
    })


    if (totals) {
      tOut <- tOut %>%
        select(grpBy, FSI, PERC_FSI, FSI_STATUS,
               FSI_INT, PERC_FSI_INT,
               TPA_RATE, BA_RATE,
               FSI_SE, PERC_FSI_SE, TPA_RATE_SE, BA_RATE_SE,
               FSI_VAR, PERC_FSI_VAR, TPA_RATE_VAR, BA_RATE_VAR,
               nPlots, N)

    } else {
      tOut <- tOut %>%
        select(grpBy, FSI, PERC_FSI, FSI_STATUS,
               FSI_INT, PERC_FSI_INT,
               TPA_RATE, BA_RATE,
               FSI_SE, PERC_FSI_SE, TPA_RATE_SE, BA_RATE_SE,
               FSI_VAR, PERC_FSI_VAR, TPA_RATE_VAR, BA_RATE_VAR,
               nPlots, N)
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





#### Pre July 28 re-write where we allows slopes to vary by stand age
# fsiStarter <- function(x,
#                        db,
#                        grpBy_quo = NULL,
#                        polys = NULL,
#                        returnSpatial = FALSE,
#                        bySpecies = FALSE,
#                        bySizeClass = FALSE,
#                        landType = 'forest',
#                        treeType = 'live',
#                        method = 'sma',
#                        lambda = .5,
#                        treeDomain = NULL,
#                        areaDomain = NULL,
#                        totals = FALSE,
#                        byPlot = FALSE,
#                        useLM = FALSE,
#                        nCores = 1,
#                        remote,
#                        mr){
#
#   reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
#                  'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
#
#   if (remote){
#     ## Store the original parameters here
#     params <- db
#
#     ## Read in one state at a time
#     db <- readFIA(dir = db$dir, common = db$common,
#                   tables = reqTables, states = x, ## x is the vector of state names
#                   nCores = nCores)
#
#     ## If a clip was specified, run it now
#     if ('mostRecent' %in% names(params)){
#       db <- clipFIA(db, mostRecent = params$mostRecent,
#                     mask = params$mask, matchEval = params$matchEval,
#                     evalid = params$evalid, designCD = params$designCD,
#                     nCores = nCores)
#     }
#
#   } else {
#     ## Really only want the required tables
#     db <- db[names(db) %in% reqTables]
#   }
#
#
#
#   ## Need a plotCN, and a new ID
#   db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
#                                 pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))
#   db$TREE$TRE_CN <- db$TREE$CN
#   ##  don't have to change original code
#   #grpBy_quo <- enquo(grpBy)
#
#   # Probably cheating, but it works
#   if (quo_name(grpBy_quo) != 'NULL'){
#     ## Have to join tables to run select with this object type
#     plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
#     ## We want a unique error message here to tell us when columns are not present in data
#     d_quo <- tryCatch(
#       error = function(cnd) {
#         return(0)
#       },
#       plt_quo[10,] %>% # Just the first row
#         left_join(select(db$COND, PLT_CN, names(db$COND)[names(db$COND) %in% names(db$PLOT) == FALSE]), by = 'PLT_CN') %>%
#         inner_join(select(db$TREE, PLT_CN, names(db$TREE)[names(db$TREE) %in% c(names(db$PLOT), names(db$COND)) == FALSE]), by = 'PLT_CN') %>%
#         select(!!grpBy_quo)
#     )
#
#     # If column doesnt exist, just returns 0, not a dataframe
#     if (is.null(nrow(d_quo))){
#       grpName <- quo_name(grpBy_quo)
#       stop(paste('Columns', grpName, 'not found in PLOT, TREE, or COND tables. Did you accidentally quote the variables names? e.g. use grpBy = ECOSUBCD (correct) instead of grpBy = "ECOSUBCD". ', collapse = ', '))
#     } else {
#       # Convert to character
#       grpBy <- names(d_quo)
#     }
#   } else {
#     grpBy <- NULL
#   }
#
#   reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
#                  'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
#
#   if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
#     stop('polys must be spatial polygons object of class sp or sf. ')
#   }
#   if (landType %in% c('timber', 'forest') == FALSE){
#     stop('landType must be one of: "forest" or "timber".')
#   }
#   if (any(reqTables %in% names(db) == FALSE)){
#     missT <- reqTables[reqTables %in% names(db) == FALSE]
#     stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
#   }
#   if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
#     warning(paste('Method', method, 'unknown. Defaulting to Temporally Indifferent (TI).'))
#   }
#
#   # I like a unique ID for a plot through time
#   if (byPlot) {grpBy <- c('pltID', grpBy)}
#   # Save original grpBy for pretty return with spatial objects
#   grpByOrig <- grpBy
#
#
#   ### DEAL WITH TEXAS
#   if (any(db$POP_EVAL$STATECD %in% 48)){
#     ## Will require manual updates, fix your shit texas
#     txIDS <- db$POP_EVAL %>%
#       filter(STATECD %in% 48) %>%
#       filter(END_INVYR < 2017) %>%
#       filter(END_INVYR > 2006) %>%
#       ## Removing any inventory that references east or west, sorry
#       filter(str_detect(str_to_upper(EVAL_DESCR), 'EAST', negate = TRUE) &
#                str_detect(str_to_upper(EVAL_DESCR), 'WEST', negate = TRUE))
#     db$POP_EVAL <- bind_rows(filter(db$POP_EVAL, !(STATECD %in% 48)), txIDS)
#   }
#
#   ### AREAL SUMMARY PREP
#   if(!is.null(polys)) {
#     # # Convert polygons to an sf object
#     # polys <- polys %>%
#     #   as('sf')%>%
#     #   mutate_if(is.factor,
#     #             as.character)
#     # ## A unique ID
#     # polys$polyID <- 1:nrow(polys)
#     #
#     # # Add shapefile names to grpBy
#     grpBy = c(grpBy, 'polyID')
#
#     ## Make plot data spatial, projected same as polygon layer
#     pltSF <- select(db$PLOT, c('LON', 'LAT', pltID)) %>%
#       filter(!is.na(LAT) & !is.na(LON)) %>%
#       distinct(pltID, .keep_all = TRUE)
#     coordinates(pltSF) <- ~LON+LAT
#     proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#     pltSF <- as(pltSF, 'sf') %>%
#       st_transform(crs = st_crs(polys))
#
#     ## Split up polys
#     polyList <- split(polys, as.factor(polys$polyID))
#     suppressWarnings({suppressMessages({
#       ## Compute estimates in parallel -- Clusters in windows, forking otherwise
#       if (Sys.info()['sysname'] == 'Windows'){
#         cl <- makeCluster(nCores)
#         clusterEvalQ(cl, {
#           library(dplyr)
#           library(stringr)
#           library(rFIA)
#           library(sf)
#         })
#         out <- parLapply(cl, X = names(polyList), fun = areal_par, pltSF, polyList)
#         #stopCluster(cl) # Keep the cluster active for the next run
#       } else { # Unix systems
#         out <- mclapply(names(polyList), FUN = areal_par, pltSF, polyList, mc.cores = nCores)
#       }
#     })})
#     pltSF <- bind_rows(out)
#
#     # A warning
#     if (length(unique(pltSF$pltID)) < 1){
#       stop('No plots in db overlap with polys.')
#     }
#     ## Add polygon names to PLOT
#     db$PLOT <- db$PLOT %>%
#       left_join(select(pltSF, polyID, pltID), by = 'pltID')
#
#
#     ## TO RETURN SPATIAL PLOTS
#   }
#   if (byPlot & returnSpatial){
#     grpBy <- c(grpBy, 'LON', 'LAT')
#   } # END AREAL
#
#   ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
#   # Land type domain indicator
#   if (tolower(landType) == 'forest'){
#     db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
#     # Tree Type domain indicator
#     if (tolower(treeType) == 'live'){
#       db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
#     } else if (tolower(treeType) == 'gs'){
#       db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
#     }
#   } else if (tolower(landType) == 'timber'){
#     db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
#     # Tree Type domain indicator
#     if (tolower(treeType) == 'live'){
#       db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
#     } else if (tolower(treeType) == 'gs'){
#       db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
#     }
#   }
#
#   # update spatial domain indicator
#   if(!is.null(polys)){
#     db$PLOT$sp <- ifelse(db$PLOT$pltID %in% pltSF$pltID, 1, 0)
#   } else {
#     db$PLOT$sp <- 1
#   }
#
#   # User defined domain indicator for area (ex. specific forest type)
#   pcEval <- left_join(db$PLOT, select(db$COND, -c('STATECD', 'UNITCD', 'COUNTYCD', 'INVYR', 'PLOT')), by = 'PLT_CN')
#   #areaDomain <- substitute(areaDomain)
#   pcEval$aD <- rlang::eval_tidy(areaDomain, pcEval) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
#   if(!is.null(pcEval$aD)) pcEval$aD[is.na(pcEval$aD)] <- 0 # Make NAs 0s. Causes bugs otherwise
#   if(is.null(pcEval$aD)) pcEval$aD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
#   pcEval$aD <- as.numeric(pcEval$aD)
#   db$COND <- left_join(db$COND, select(pcEval, c('PLT_CN', 'CONDID', 'aD')), by = c('PLT_CN', 'CONDID')) %>%
#     mutate(aD_c = aD)
#   aD_p <- pcEval %>%
#     group_by(PLT_CN) %>%
#     summarize(aD_p = as.numeric(any(aD > 0)))
#   db$PLOT <- left_join(db$PLOT, aD_p, by = 'PLT_CN')
#   rm(pcEval)
#
#   # Same as above for tree (ex. trees > 20 ft tall)
#   #treeDomain <- substitute(treeDomain)
#   tD <- rlang::eval_tidy(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
#   if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
#   if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
#   db$TREE$tD <- as.numeric(tD)
#
#
#
#   ## We only want inventory/ population info from t2 plots, but we need the plot tree cond data
#   ## for t1 and t2
#   remPlts <- db$PLOT %>%
#     select(PLT_CN, PREV_PLT_CN, DESIGNCD, REMPER, PLOT_STATUS_CD) %>%
#     ## Has to have a remeasurement, be in the current sample, and of the national design
#     filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & PLOT_STATUS_CD != 3 & DESIGNCD == 1) %>%
#     left_join(select(db$PLOT, PLT_CN, DESIGNCD, PLOT_STATUS_CD), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
#     ## past emasurement must be in the previous sample and of national design
#     filter(PLOT_STATUS_CD.prev != 3 & DESIGNCD.prev == 1)
#
#   ### Snag the EVALIDs that are needed
#   db$POP_EVAL<- db$POP_EVAL %>%
#     select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD') %>%
#     inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
#     filter(EVAL_TYP == 'EXPVOL') %>%
#     filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
#     distinct(END_INVYR, EVALID, .keep_all = TRUE)
#
#   ## If a most-recent subset, make sure that we don't get two reporting years in
#   ## western states
#   if (mr) {
#     db$POP_EVAL <- db$POP_EVAL %>%
#       group_by(EVAL_TYP) %>%
#       filter(END_INVYR == max(END_INVYR, na.rm = TRUE))
#   }
#
#   ## Make an annual panel ID, associated with an INVYR
#
#   ### The population tables
#   pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'END_INVYR', 'EVAL_TYP')) %>%
#     rename(EVAL_CN = CN) %>%
#     left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
#     rename(ESTN_UNIT_CN = CN) %>%
#     left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
#     rename(STRATUM_CN = CN) %>%
#     left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
#     ## ONLY REMEASURED PLOTS MEETING CRITERIA ABOVE
#     filter(PLT_CN %in% remPlts$PLT_CN) %>%
#     ungroup() %>%
#     mutate_if(is.factor,
#               as.character)
#
#   ### Which estimator to use?
#   if (str_to_upper(method) %in% c('ANNUAL')){
#     ## Want to use the year where plots are measured, no repeats
#     ## Breaking this up into pre and post reporting becuase
#     ## Estimation units get weird on us otherwise
#     popOrig <- pops
#     pops <- pops %>%
#       group_by(STATECD) %>%
#       filter(END_INVYR == INVYR) %>%
#       ungroup()
#
#     prePops <- popOrig %>%
#       group_by(STATECD) %>%
#       filter(INVYR < min(END_INVYR, na.rm = TRUE)) %>%
#       distinct(PLT_CN, .keep_all = TRUE) %>%
#       ungroup()
#
#     pops <- bind_rows(pops, prePops) %>%
#       mutate(YEAR = INVYR)
#
#   } else {     # Otherwise temporally indifferent
#     pops <- rename(pops, YEAR = END_INVYR)
#   }
#
#   ## P2POINTCNT column is NOT consistent for annnual estimates, plots
#   ## within individual strata and est units are related to different INVYRs
#   p2_INVYR <- pops %>%
#     group_by(ESTN_UNIT_CN, STRATUM_CN, INVYR) %>%
#     summarize(P2POINTCNT_INVYR = length(unique(PLT_CN)))
#   ## Want a count of p2 points / eu, gets screwed up with grouping below
#   p2eu_INVYR <- p2_INVYR %>%
#     distinct(ESTN_UNIT_CN, STRATUM_CN, INVYR, .keep_all = TRUE) %>%
#     group_by(ESTN_UNIT_CN, INVYR) %>%
#     summarize(p2eu_INVYR = sum(P2POINTCNT_INVYR, na.rm = TRUE))
#   p2eu <- pops %>%
#     distinct(ESTN_UNIT_CN, STRATUM_CN, .keep_all = TRUE) %>%
#     group_by(ESTN_UNIT_CN) %>%
#     summarize(p2eu = sum(P2POINTCNT, na.rm = TRUE))
#
#   ## Rejoin
#   pops <- pops %>%
#     left_join(p2_INVYR, by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'INVYR')) %>%
#     left_join(p2eu_INVYR, by = c('ESTN_UNIT_CN', 'INVYR')) %>%
#     left_join(p2eu, by = 'ESTN_UNIT_CN')
#
#
#   ## Recode a few of the estimation methods to make things easier below
#   pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
#                             `Post-Stratification` = 'strat',
#                             `Stratified random sampling` = 'strat',
#                             `Double sampling for stratification` = 'double',
#                             `Simple random sampling` = 'simple',
#                             `Subsampling units of unequal size` = 'simple')
#
#
#   ## Add species to groups
#   if (bySpecies) {
#     db$TREE <- db$TREE %>%
#       left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
#       mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' ')) %>%
#       ungroup() %>%
#       mutate_if(is.factor,
#                 as.character)
#     grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
#     grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
#   }
#
#   ## Break into size classes
#   if (bySizeClass){
#     grpBy <- c(grpBy, 'sizeClass')
#     grpByOrig <- c(grpByOrig, 'sizeClass')
#     db$TREE$sizeClass <- makeClasses(db$TREE$DIA, interval = 2, numLabs = TRUE)
#     db$TREE <- db$TREE[!is.na(db$TREE$sizeClass),]
#   }
#
#   ## Only the necessary plots for EVAL of interest
#   remPltList <- unique(c(remPlts$PLT_CN, remPlts$PREV_PLT_CN))
#   db$PLOT <- filter(db$PLOT, PLT_CN %in% remPltList)
#   db$COND <- filter(db$COND, PLT_CN %in% remPltList)
#   db$TREE <- filter(db$TREE, PLT_CN %in% remPltList)
#
#   ## Tree basal area per acre
#   db$TREE <- db$TREE %>%
#     mutate(BAA = basalArea(DIA) * TPA_UNADJ)
#
#   ## Narrow up the tables to the necessary variables, reduces memory load
#   ## sent out the cores
#   ## Which grpByNames are in which table? Helps us subset below
#   grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
#   grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
#   grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]
#
#   ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
#   db$PLOT <- select(db$PLOT, c('PLT_CN', pltID, 'REMPER', 'DESIGNCD', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR',
#                                'MEASYEAR', 'MEASMON', 'MEASDAY', 'PLOT_STATUS_CD', PREV_PLT_CN, grpP, 'aD_p', 'sp'))
#   db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD'))
#   db$TREE <- select(db$TREE, c('PLT_CN', 'TRE_CN', 'CONDID', 'DIA', 'TPA_UNADJ', 'BAA', 'SUBP', 'TREE', grpT, 'tD', 'typeD',
#                                PREVCOND, PREV_TRE_CN, STATUSCD, SPCD))
#
#
#   ## Merging state and county codes
#   plts <- split(db$PLOT, as.factor(paste(db$PLOT$COUNTYCD, db$PLOT$STATECD, sep = '_')))
#
#
#   ## Use the linear model procedures or strictly t2-t1 / remper?
#   if (useLM){
#     suppressWarnings({
#       ## Compute estimates in parallel -- Clusters in windows, forking otherwise
#       if (Sys.info()['sysname'] == 'Windows'){
#         cl <- makeCluster(nCores)
#         clusterEvalQ(cl, {
#           library(dplyr)
#           library(stringr)
#           library(rFIA)
#           library(tidyr)
#           library(purrr)
#         })
#         out <- parLapply(cl, X = names(plts), fun = fsiHelper1_lm, plts, db, grpBy, byPlot)
#         #stopCluster(cl) # Keep the cluster active for the next run
#       } else { # Unix systems
#         out <- mclapply(names(plts), FUN = fsiHelper1_lm, plts, db, grpBy, byPlot, mc.cores = nCores)
#       }
#     })
#   } else {
#     suppressWarnings({
#       ## Compute estimates in parallel -- Clusters in windows, forking otherwise
#       if (Sys.info()['sysname'] == 'Windows'){
#         cl <- makeCluster(nCores)
#         clusterEvalQ(cl, {
#           library(dplyr)
#           library(stringr)
#           library(rFIA)
#           library(tidyr)
#         })
#         out <- parLapply(cl, X = names(plts), fun = fsiHelper1, plts, db, grpBy, byPlot)
#         #stopCluster(cl) # Keep the cluster active for the next run
#       } else { # Unix systems
#         out <- mclapply(names(plts), FUN = fsiHelper1, plts, db, grpBy, byPlot, mc.cores = nCores)
#       }
#     })
#   }
#
#   ## back to dataframes
#   out <- unlist(out, recursive = FALSE)
#   t <- bind_rows(out[names(out) == 't'])
#   a <- bind_rows(out[names(out) == 'a'])
#
#   out <- list(t = t, a = a, grpBy = grpBy, grpByOrig = grpByOrig, pops = pops)
#
# }



# fsi <- function(db,
#                 grpBy = NULL,
#                 polys = NULL,
#                 returnSpatial = FALSE,
#                 bySpecies = FALSE,
#                 bySizeClass = FALSE,
#                 landType = 'forest',
#                 treeType = 'live',
#                 method = 'sma',
#                 lambda = .5,
#                 treeDomain = NULL,
#                 areaDomain = NULL,
#                 totals = TRUE,
#                 byPlot = FALSE,
#                 useLM = FALSE,
#                 nCores = 1) {
#
#   ##  don't have to change original code
#   grpBy_quo <- rlang::enquo(grpBy)
#   areaDomain <- rlang::enquo(areaDomain)
#   treeDomain <- rlang::enquo(treeDomain)
#
#   ### Is DB remote?
#   remote <- ifelse(class(db) == 'Remote.FIA.Database', 1, 0)
#   if (remote){
#
#     iter <- db$states
#
#     ## In memory
#   } else {
#     ## Some warnings
#     if (class(db) != "FIA.Database"){
#       stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
#     }
#
#     ## an iterator for remote
#     iter <- 1
#
#   }
#
#   ## Check for a most recent subset
#   if (remote){
#     if ('mostRecent' %in% names(db)){
#       mr = db$mostRecent # logical
#     } else {
#       mr = FALSE
#     }
#     ## In-memory
#   } else {
#     if ('mostRecent' %in% names(db)){
#       mr = TRUE
#     } else {
#       mr = FALSE
#     }
#   }
#
#   ### AREAL SUMMARY PREP
#   if(!is.null(polys)) {
#     # Convert polygons to an sf object
#     polys <- polys %>%
#       as('sf') %>%
#       mutate_if(is.factor,
#                 as.character)
#     ## A unique ID
#     polys$polyID <- 1:nrow(polys)
#   }
#
#
#   ## Run the main portion
#   out <- lapply(X = iter, FUN = fsiStarter, db,
#                 grpBy_quo = grpBy_quo, polys, returnSpatial,
#                 bySpecies, bySizeClass,
#                 landType, treeType, method,
#                 lambda, treeDomain, areaDomain,
#                 totals, byPlot, useLM,
#                 nCores, remote, mr)
#
#   ## Bring the results back
#   out <- unlist(out, recursive = FALSE)
#   t <- bind_rows(out[names(out) == 't'])
#   a <- bind_rows(out[names(out) == 'a'])
#   grpBy <- out[names(out) == 'grpBy'][[1]]
#   grpByOrig <- out[names(out) == 'grpByOrig'][[1]]
#
#
#
#   if (byPlot){
#
#     ## back to dataframes
#     tOut <- t
#
#
#     grpScale <- tOut %>%
#       filter(PLOT_STATUS_CD == 1) %>%
#       group_by(.dots = grpBy[grpBy %in% c('YEAR', 'pltID') == FALSE]) %>%
#       summarise(tpaSD = sd(c(PREV_TPA, CURR_TPA), na.rm = TRUE),
#                 baSD = sd(c(PREV_BA, CURR_BA), na.rm = TRUE))
#
#     ## Scale the indices by group
#     tOut <- tOut %>%
#       left_join(grpScale, by = grpBy[grpBy %in% c('YEAR', 'pltID') == FALSE]) %>%
#       mutate(dt = CHNG_TPA / tpaSD,
#              db = CHNG_BA / baSD,
#              t1 = PREV_TPA / tpaSD / REMPER,
#              t2 = CURR_TPA / tpaSD/ REMPER,
#              b1 = PREV_BA / baSD / REMPER,
#              b2 = CURR_BA / baSD / REMPER) %>%
#       ## The FSI and % FSI
#       mutate(FSI = projectPoints(dt, db, 0.8025, 0, returnPoint = FALSE),
#              ## can and will produce INF here due to division by zero, that's fine, just use the FSI if that matters to you
#              PERC_FSI = projectPoints(t2, b2, 0.8025, 0, returnPoint = FALSE) / projectPoints(t1, b1, 0.8025, 0, returnPoint = FALSE) -1,
#              TPA_RATE = dt,
#              BA_RATE = db)
#
#
#     ## Make it spatial
#     if (returnSpatial){
#       tOut <- tOut %>%
#         filter(!is.na(LAT) & !is.na(LON)) %>%
#         st_as_sf(coords = c('LON', 'LAT'),
#                  crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
#       grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
#
#     }
#
#
#     tOut <- select(tOut, YEAR, PLT_CN, any_of('PREV_PLT_CN'), PLOT_STATUS_CD, grpBy[grpBy != 'YEAR'], FSI, PERC_FSI, TPA_RATE, BA_RATE, REMPER, everything()) %>%
#       select(-c(dt, db, t1, t2, b1, b2))
#
#
#
#     ## Population estimation
#   } else {
#     ## back to dataframes
#     pops <- bind_rows(out[names(out) == 'pops'])
#
#     ### CANNOT USE vectors as they are to compute SD, because of PLOT_BASIS issues
#     ### SD is unnessarily high --> plot level first
#     pltRates <- t %>%
#       ungroup() %>%
#       filter(plotIn == 1) %>%
#       select(PLT_CN, REMPER, CHNG_TPA, CHNG_BAA, PREV_BAA, PREV_TPA, plotIn, grpBy) %>%
#       ## Replace any NAs with zeros
#       mutate(PREV_BAA = replace_na(PREV_BAA, 0),
#              PREV_TPA = replace_na(PREV_TPA, 0),
#              CHNG_BAA = replace_na(CHNG_BAA, 0),
#              CHNG_TPA = replace_na(CHNG_TPA, 0)) %>%
#       group_by(PLT_CN, plotIn, .dots = grpBy) %>%
#       ## Sum across micro, subp, and macroplots
#       summarize(PREV_BAA = sum(PREV_BAA, na.rm = TRUE),
#                 PREV_TPA = sum(PREV_TPA, na.rm = TRUE),
#                 CHNG_BAA = sum(CHNG_BAA, na.rm = TRUE),
#                 CHNG_TPA = sum(CHNG_TPA, na.rm = TRUE),
#                 REMPER = first(REMPER)) %>%
#       ## T2 attributes
#       mutate(CURR_BAA = PREV_BAA + CHNG_BAA,
#              CURR_TPA = PREV_TPA + CHNG_TPA) %>%
#       ## Change in average tree BA and QMD
#       mutate(PREV_BA = if_else(PREV_TPA != 0, PREV_BAA / PREV_TPA, 0),
#              CURR_BA = if_else(CURR_TPA != 0, CURR_BAA / CURR_TPA, 0),
#              CHNG_BA = CURR_BA - PREV_BA) %>%
#       ## SD by group
#       group_by(.dots = grpBy) %>%
#       summarise(tpaSD = sd(c(PREV_TPA, CURR_TPA), na.rm = TRUE),
#                 baSD = sd(c(PREV_BA, CURR_BA), na.rm = TRUE)) %>%
#       ungroup()
#
#
#     ## Adding YEAR to groups
#     grpBy <- c('YEAR', grpBy)
#     popState <- split(pops, as.factor(pops$STATECD))
#
#
#     suppressWarnings({
#       ## Compute estimates in parallel -- Clusters in windows, forking otherwise
#       if (Sys.info()['sysname'] == 'Windows'){
#         cl <- makeCluster(nCores)
#         clusterEvalQ(cl, {
#           library(dplyr)
#           library(stringr)
#           library(rFIA)
#           library(tidyr)
#         })
#         out <- parLapply(cl, X = names(popState), fun = fsiHelper2, popState, t, a, grpBy, method, pltRates)
#         stopCluster(cl)
#       } else { # Unix systems
#         out <- mclapply(names(popState), FUN = fsiHelper2, popState, t, a, grpBy, method, pltRates, mc.cores = nCores)
#       }
#     })
#     ## back to dataframes
#     out <- unlist(out, recursive = FALSE)
#     tEst <- bind_rows(out[names(out) == 'tEst'])
#     tEst <- ungroup(tEst)
#
#     ##### ----------------- MOVING AVERAGES
#     if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
#       ### ---- SIMPLE MOVING AVERAGE
#       if (str_to_upper(method) == 'SMA'){
#         ## Assuming a uniform weighting scheme
#         wgts <- pops %>%
#           group_by(ESTN_UNIT_CN) %>%
#           summarize(wgt = 1 / length(unique(INVYR)))
#
#         tEst <- left_join(tEst, wgts, by = 'ESTN_UNIT_CN')
#
#         #### ----- Linear MOVING AVERAGE
#       } else if (str_to_upper(method) == 'LMA'){
#         wgts <- pops %>%
#           distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
#           arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
#           group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
#           mutate(rank = min_rank(INVYR))
#
#         ## Want a number of INVYRs per EU
#         neu <- wgts %>%
#           group_by(ESTN_UNIT_CN) %>%
#           summarize(n = sum(rank, na.rm = TRUE))
#
#         ## Rejoining and computing wgts
#         wgts <- wgts %>%
#           left_join(neu, by = 'ESTN_UNIT_CN') %>%
#           mutate(wgt = rank / n) %>%
#           ungroup() %>%
#           select(ESTN_UNIT_CN, INVYR, wgt)
#
#         tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
#
#         #### ----- EXPONENTIAL MOVING AVERAGE
#       } else if (str_to_upper(method) == 'EMA'){
#         wgts <- pops %>%
#           distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
#           arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
#           group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
#           mutate(rank = min_rank(INVYR))
#
#
#         if (length(lambda) < 2){
#           ## Want sum of weighitng functions
#           neu <- wgts %>%
#             mutate(l = lambda) %>%
#             group_by(ESTN_UNIT_CN) %>%
#             summarize(l = 1 - first(lambda),
#                       sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))
#
#           ## Rejoining and computing wgts
#           wgts <- wgts %>%
#             left_join(neu, by = 'ESTN_UNIT_CN') %>%
#             mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
#             ungroup() %>%
#             select(ESTN_UNIT_CN, INVYR, wgt)
#         } else {
#           grpBy <- c('lambda', grpBy)
#           ## Duplicate weights for each level of lambda
#           yrWgts <- list()
#           for (i in 1:length(unique(lambda))) {
#             yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
#           }
#           wgts <- bind_rows(yrWgts)
#           ## Want sum of weighitng functions
#           neu <- wgts %>%
#             group_by(lambda, ESTN_UNIT_CN) %>%
#             summarize(l = 1 - first(lambda),
#                       sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))
#
#           ## Rejoining and computing wgts
#           wgts <- wgts %>%
#             left_join(neu, by = c('lambda', 'ESTN_UNIT_CN')) %>%
#             mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
#             ungroup() %>%
#             select(lambda, ESTN_UNIT_CN, INVYR, wgt)
#         }
#
#         tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
#
#       }
#
#       ### Applying the weights
#       tEst <- tEst %>%
#         mutate_at(vars(ctEst:faEst), ~(.*wgt)) %>%
#         mutate_at(vars(ctVar:cvEst_psi), ~(.*(wgt^2))) %>%
#         group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
#         summarize_at(vars(ctEst:plotIn_t), sum, na.rm = TRUE)
#     }
#
#     suppressMessages({suppressWarnings({
#       ## If a clip was specified, handle the reporting years
#       if (mr){
#         ## If a most recent subset, ignore differences in reporting years across states
#         ## instead combine most recent information from each state
#         # ID mr years by group
#         maxyearsT <- tEst %>%
#           select(grpBy) %>%
#           group_by(.dots = grpBy[!c(grpBy %in% 'YEAR')]) %>%
#           summarise(YEAR = max(YEAR, na.rm = TRUE))
#
#         # Combine estimates
#         tEst <- tEst %>%
#           ungroup() %>%
#           select(-c(YEAR)) %>%
#           left_join(maxyearsT, by = grpBy[!c(grpBy %in% 'YEAR')])
#       }
#     })})
#
#
#     ##---------------------  TOTALS and RATIOS
#     # Tree
#     tTotal <- tEst %>%
#       group_by(.dots = grpBy) %>%
#       summarize_all(sum,na.rm = TRUE)
#
#
#
#     ##---------------------  TOTALS and RATIOS
#     suppressWarnings({
#       tOut <- tTotal %>%
#         group_by(.dots = grpBy) %>%
#         summarize_all(sum,na.rm = TRUE) %>%
#         mutate(TPA_RATE = ctEst / ptEst,
#                BA_RATE = cbEst / pbEst,
#                FSI = siEst / faEst,
#                PERC_FSI = siEst / si1Est,
#
#                ## Ratio variance
#                ctVar = (1/ptEst^2) * (ctVar + (TPA_RATE^2 * ptVar) - (2 * TPA_RATE * cvEst_ct)),
#                cbVar = (1/pbEst^2) * (cbVar + (BA_RATE^2 * pbVar) - (2 * BA_RATE * cvEst_cb)),
#                psiVar = (1/si1Est^2) * (siVar + (PERC_FSI^2 * si1Var) - (2 * PERC_FSI * cvEst_psi)),
#                siVar = (1/faEst^2) * (siVar + (FSI^2 * faVar) - (2 * FSI * cvEst_si)),
#
#                ## Make it a percent
#                PERC_FSI = PERC_FSI * 100,
#                psiVar = psiVar * (100^2),
#
#                ## RATIO SE
#                TPA_RATE_SE = sqrt(ctVar) / abs(TPA_RATE) * 100,
#                BA_RATE_SE = sqrt(cbVar) / abs(BA_RATE) * 100,
#                FSI_SE = sqrt(siVar) / abs(FSI) * 100,
#                FSI_VAR = siVar,
#                PERC_FSI_SE = sqrt(psiVar) / abs(PERC_FSI) * 100,
#                PERC_FSI_VAR = psiVar,
#                TPA_RATE_VAR = ctVar,
#                BA_RATE_VAR = cbVar,
#                nPlots = plotIn_t,
#                N = nh,
#                FSI_INT = qt(.975, df=N-1) * (sqrt(siVar)/sqrt(N)),
#                PERC_FSI_INT = qt(.975, df=N-1) * (sqrt(psiVar)/sqrt(N))) %>%
#         mutate(FSI_STATUS = case_when(
#           FSI < 0 & FSI + FSI_INT < 0 ~ 'Decline',
#           FSI < 0 & FSI + FSI_INT > 0 ~ 'Stable',
#           FSI > 0 & FSI - FSI_INT > 0  ~ 'Expand',
#           TRUE ~ 'Stable'
#         ))
#     })
#
#
#     if (totals) {
#       tOut <- tOut %>%
#         select(grpBy, FSI, PERC_FSI, FSI_STATUS,
#                FSI_INT, PERC_FSI_INT,
#                TPA_RATE, BA_RATE,
#                FSI_SE, PERC_FSI_SE, TPA_RATE_SE, BA_RATE_SE,
#                FSI_VAR, PERC_FSI_VAR, TPA_RATE_VAR, BA_RATE_VAR,
#                nPlots, N)
#
#     } else {
#       tOut <- tOut %>%
#         select(grpBy, FSI, PERC_FSI, FSI_STATUS,
#                FSI_INT, PERC_FSI_INT,
#                TPA_RATE, BA_RATE,
#                FSI_SE, PERC_FSI_SE, TPA_RATE_SE, BA_RATE_SE,
#                FSI_VAR, PERC_FSI_VAR, TPA_RATE_VAR, BA_RATE_VAR,
#                nPlots, N)
#     }
#
#     # Snag the names
#     tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]
#
#   }
#
#   ## Pretty output
#   tOut <- tOut %>%
#     ungroup() %>%
#     mutate_if(is.factor, as.character) %>%
#     drop_na(grpBy) %>%
#     arrange(YEAR) %>%
#     as_tibble()
#
#   # Return a spatial object
#   if (!is.null(polys) & byPlot == FALSE) {
#     ## NO IMPLICIT NA
#     nospGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE])
#     nospSym <- syms(nospGrp)
#     tOut <- complete(tOut, !!!nospSym)
#     ## If species, we don't want unique combos of variables related to same species
#     ## but we do want NAs in polys where species are present
#     if (length(nospGrp) < length(grpBy)){
#       spGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME')])
#       spSym <- syms(spGrp)
#       tOut <- complete(tOut, nesting(!!!nospSym))
#     }
#
#     suppressMessages({suppressWarnings({
#       tOut <- left_join(tOut, polys, by = 'polyID') %>%
#         select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
#         filter(!is.na(polyID))})})
#
#     ## Makes it horrible to work with as a dataframe
#     if (returnSpatial == FALSE) tOut <- select(tOut, -c(geometry))
#   } else if (!is.null(polys) & byPlot){
#     polys <- as.data.frame(polys)
#     tOut <- left_join(tOut, select(polys, -c(geometry)), by = 'polyID')
#   }
#
#
#   ## For spatial plots
#   if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
#
#   ## Above converts to tibble
#   if (returnSpatial) tOut <- st_sf(tOut)
#   # ## remove any duplicates in byPlot (artifact of END_INYR loop)
#   if (byPlot) tOut <- unique(tOut)
#
#
#   return(tOut)
#
# }
