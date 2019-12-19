#' @export
biomass <- function(db,
                       grpBy = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       bySpecies = FALSE,
                       bySizeClass = FALSE,
                       landType = 'forest',
                       treeType = 'live',
                       method = 'TI',
                       lambda = .94,
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
  if (str_to_upper(method) %in% c('ANNUAL', "SMA", 'EMA', 'LMA')){
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
      out <- parLapply(cl, X = names(plts), fun = bioHelper1, plts, db, grpBy, aGrpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = bioHelper1, plts, db, grpBy, aGrpBy, byPlot, mc.cores = nCores)
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
        out <- parLapply(cl, X = names(popState), fun = bioHelper2, popState, a, t, grpBy, aGrpBy)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = bioHelper2, popState, a, t, grpBy, aGrpBy, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    aEst <- bind_rows(out[names(out) == 'aEst'])
    tEst <- bind_rows(out[names(out) == 'tEst'])


    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
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
        summarize_at(vars(nvEst:plotIn_TREE),sum, na.rm = TRUE)

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
          ungroup() %>%
          mutate(diffYear = (INVYR - minyr) + 1)
        ## Having trouble with dplyr here, hacky
        wgts$id <- 1:nrow(wgts)
        wgts <- wgts %>%
          group_by(YEAR, INVYR, STATECD, id) %>%
          summarize(wgt = diffYear / sum(1:n, na.rm = TRUE))


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
        group_by(STATECD, .dots = aGrpBy) %>%
        summarize_at(vars(aEst:plotIn_AREA), sum, na.rm = TRUE)


      tEst <- left_join(wgts, tEst, by = c('INVYR' = 'YEAR', 'STATECD')) %>%
        mutate_at(vars(nvEst:ctEst), ~(.*wgt)) %>%
        mutate_at(vars(nvVar:cvEst_ct), ~(.*(wgt^2))) %>%
        group_by(STATECD, .dots = grpBy) %>%
        summarize_at(vars(nvEst:plotIn_TREE), sum, na.rm = TRUE)

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
        mutate(AREA_TOTAL = aEst,
               NETVOL_TOTAL = nvEst,
               SAWVOL_TOTAL = svEst,
               BIO_AG_TOTAL = bagEst,
               BIO_BG_TOTAL = bbgEst,
               BIO_TOTAL = btEst,
               CARB_AG_TOTAL = cagEst,
               CARB_BG_TOTAL = cbgEst,
               CARB_TOTAL = ctEst,
               ## Ratios
               NETVOL_ACRE = NETVOL_TOTAL / AREA_TOTAL,
               SAWVOL_ACRE = SAWVOL_TOTAL / AREA_TOTAL,
               BIO_AG_ACRE = BIO_AG_TOTAL / AREA_TOTAL,
               BIO_BG_ACRE = BIO_BG_TOTAL / AREA_TOTAL,
               BIO_ACRE = BIO_TOTAL / AREA_TOTAL,
               CARB_AG_ACRE = CARB_AG_TOTAL / AREA_TOTAL,
               CARB_BG_ACRE = CARB_BG_TOTAL / AREA_TOTAL,
               CARB_ACRE = CARB_TOTAL / AREA_TOTAL,
               ## Ratio Var
               nvaVar = (1/AREA_TOTAL^2) * (nvVar + (NETVOL_ACRE^2 * aVar) - 2 * NETVOL_ACRE * cvEst_nv),
               svaVar = (1/AREA_TOTAL^2) * (svVar + (SAWVOL_ACRE^2 * aVar) - 2 * SAWVOL_ACRE * cvEst_sv),
               baaVar = (1/AREA_TOTAL^2) * (bagVar + (BIO_AG_ACRE^2 * aVar) - 2 * BIO_AG_ACRE * cvEst_bag),
               bbaVar = (1/AREA_TOTAL^2) * (bbgVar + (BIO_BG_ACRE^2 * aVar) - 2 * BIO_BG_ACRE * cvEst_bbg),
               btaVar = (1/AREA_TOTAL^2) * (btVar + (BIO_ACRE^2 * aVar) - 2 * BIO_ACRE * cvEst_bt),
               caaVar = (1/AREA_TOTAL^2) * (cagVar + (CARB_AG_ACRE^2 * aVar) - 2 * CARB_AG_ACRE * cvEst_cag),
               cbaVar = (1/AREA_TOTAL^2) * (bbgVar + (CARB_BG_ACRE^2 * aVar) - 2 * CARB_BG_ACRE * cvEst_bbg),
               ctaVar = (1/AREA_TOTAL^2) * (ctVar + (CARB_ACRE^2 * aVar) - 2 * CARB_ACRE * cvEst_ct),
               ## SE RATIO
               NETVOL_ACRE_SE = sqrt(nvaVar) / NETVOL_ACRE *100,
               SAWVOL_ACRE_SE = sqrt(svaVar) / SAWVOL_ACRE *100,
               BIO_AG_ACRE_SE = sqrt(baaVar) / BIO_AG_ACRE *100,
               BIO_BG_ACRE_SE = sqrt(bbaVar) / BIO_BG_ACRE *100,
               BIO_ACRE_SE = sqrt(btVar) / BIO_ACRE * 100,
               CARB_AG_ACRE_SE = sqrt(caaVar) / CARB_AG_ACRE *100,
               CARB_BG_ACRE_SE = sqrt(cbaVar) / CARB_BG_ACRE *100,
               CARB_ACRE_SE = sqrt(ctVar) / CARB_ACRE *100,
               ## SE TOTAL
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL *100,
               NETVOL_TOTAL_SE = sqrt(nvVar) / NETVOL_TOTAL *100,
               SAWVOL_TOTAL_SE = sqrt(svVar) / SAWVOL_TOTAL *100,
               BIO_AG_TOTAL_SE = sqrt(bagVar) / BIO_AG_TOTAL *100,
               BIO_BG_TOTAL_SE = sqrt(bbgVar) / BIO_BG_TOTAL *100,
               BIO_TOTAL_SE = sqrt(btVar) / BIO_TOTAL *100,
               CARB_AG_TOTAL_SE = sqrt(cagVar) / CARB_AG_TOTAL *100,
               CARB_BG_TOTAL_SE = sqrt(cbgVar) / CARB_BG_TOTAL *100,
               CARB_TOTAL_SE = sqrt(ctVar) / CARB_TOTAL *100,
               ## nPlots
               nPlots_TREE = plotIn_TREE,
               nPlots_AREA = plotIn_AREA)
    })


    if (totals) {

      tOut <- tOut %>%
        select(grpBy, "NETVOL_ACRE","SAWVOL_ACRE","BIO_AG_ACRE","BIO_BG_ACRE",
               "BIO_ACRE","CARB_AG_ACRE","CARB_BG_ACRE","CARB_ACRE","NETVOL_TOTAL",
               "SAWVOL_TOTAL","BIO_AG_TOTAL","BIO_BG_TOTAL","BIO_TOTAL","CARB_AG_TOTAL",
               "CARB_BG_TOTAL","CARB_TOTAL", "AREA_TOTAL","NETVOL_ACRE_SE",
               "SAWVOL_ACRE_SE","BIO_AG_ACRE_SE", "BIO_BG_ACRE_SE", "BIO_ACRE_SE",
               "CARB_AG_ACRE_SE","CARB_BG_ACRE_SE","CARB_ACRE_SE","NETVOL_TOTAL_SE",
               "SAWVOL_TOTAL_SE",  "BIO_AG_TOTAL_SE",  "BIO_BG_TOTAL_SE",
               "BIO_TOTAL_SE", "CARB_AG_TOTAL_SE", "CARB_BG_TOTAL_SE", "CARB_TOTAL_SE",
               "AREA_TOTAL_SE","nPlots_TREE","nPlots_AREA")

    } else {
      tOut <- tOut %>%
        select(grpBy, "NETVOL_ACRE","SAWVOL_ACRE","BIO_AG_ACRE","BIO_BG_ACRE",
               "BIO_ACRE","CARB_AG_ACRE","CARB_BG_ACRE","CARB_ACRE","NETVOL_ACRE_SE",
               "SAWVOL_ACRE_SE","BIO_AG_ACRE_SE", "BIO_BG_ACRE_SE", "BIO_ACRE_SE",
               "CARB_AG_ACRE_SE","CARB_BG_ACRE_SE","CARB_ACRE_SE","nPlots_TREE",
               "nPlots_AREA")
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
  if (!is.null(polys)) {
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
      tOut <- left_join(tOut, polys) %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})

    ## Makes it horrible to work with as a dataframe
    if (returnSpatial == FALSE) tOut <- select(tOut, -c(geometry))
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
  return(tOut)
}










