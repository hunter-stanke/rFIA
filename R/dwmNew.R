#' @export
dwm <- function(db,
                           grpBy = NULL,
                           polys = NULL,
                           returnSpatial = FALSE,
                           landType = 'forest',
                           method = 'TI',
                           lambda = .94,
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

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Check if column exists
    allNames <- c(names(db$PLOT), names(db$COND))

    if (quo_name(grpBy_quo) %in% allNames){
      # Convert to character
      grpBy <- quo_name(grpBy_quo)
    } else {
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT or COND tables. Did you accidentally quote the variables names? e.g. use grpBy = ECOSUBCD (correct) instead of grpBy = "ECOSUBCD". ', collapse = ', '))
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
      left_join(select(pltSF, polyID, pltID), by = 'pltID')

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
        out <- parLapply(cl, X = names(popState), fun = dwmHelper2, popState, t, grpBy)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = dwmHelper2, popState, t, grpBy, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])



    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
      if ('STATECD' %in% names(tEst) == FALSE){
        ## Need a STATECD on aEst and tEst to join wgts
        tEst <- left_join(tEst, select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN'))
      }



      #### Summarizing to state level here to apply weights by panel
      #### Getting rid of ESTN_UNITS
      # Tree
      tEst <- tEst %>%
        group_by(STATECD, .dots = grpBy) %>%
        summarize_at(vars(aEst:cvEst_c),sum, na.rm = TRUE)

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
          #aGrpBy <- c('lambda', aGrpBy)
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
      tEst <- left_join(wgts, tEst, by = c('INVYR' = 'YEAR', 'STATECD')) %>%
        mutate_at(vars(aEst:cEst), ~(.*wgt)) %>%
        mutate_at(vars(aVar:cvEst_c), ~(.*(wgt^2))) %>%
        group_by(STATECD, .dots = grpBy) %>%
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
  tOut <- drop_na(tOut, grpBy[grpBy %in% c('polyID') == FALSE]) %>%
    arrange(YEAR) %>%
    as_tibble()

  # Return a spatial object
  if (!is.null(polys)) {
    ### NO IMPLICIT NA
    grpSym <- syms(grpBy)
    combos <- tOut %>%
      expand(!!!grpSym)
    tOut <- left_join(combos, tOut, by = grpBy)
    suppressMessages({suppressWarnings({tOut <- left_join(tOut, polys, by = 'polyID') %>%
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


