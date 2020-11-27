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
                       remote,
                       mr){



  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'COND_DWM_CALC', 'COND', 'POP_PLOT_STRATUM_ASSGN',
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
  db$COND_DWM_CALC <- db[['COND_DWM_CALC']] %>% mutate(DWM_CN = CN)
  db$COND <- db[['COND']] %>% mutate(CND_CN = CN)

  ## Convert grpBy to character
  grpBy <- grpByToChar(db, grpBy_quo)

  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy


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

  ## Spatial boundary
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(!is.na(db$PLOT$polyID), 1, 0)
  } else {
    db$PLOT$sp <- 1
  }

  # User defined domain indicator for area (ex. specific forest type)
  db <- udAreaDomain(db, areaDomain)





  ## Handle population tables --------------------------------------------------
  ## Filtering out all inventories that are not relevant to the current estimation
  ## type. If using estimator other than TI, handle the differences in P2POINTCNT
  ## and in assigning YEAR column (YEAR = END_INVYR if method = 'TI')
  pops <- handlePops(db, evalType = c('EXPDWM'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (method != 'TI') {
    pops <- mergeSmallStrata(db, pops)
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

  ### Only joining tables necessary to produce plot level estimates
  db$PLOT <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA',
                               'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD',
                               all_of(grpP), 'aD_p', 'sp', 'COUNTYCD'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS',
                               'COND_STATUS_CD', 'CONDID',
                               all_of(grpC), 'aD_c', 'landD')) %>%
    filter(PLT_CN %in% db$PLOT$PLT_CN)
  db$COND_DWM_CALC <- select(db$COND_DWM_CALC, -c( 'STATECD', 'COUNTYCD',
                                                   'UNITCD', 'INVYR',
                                                   'MEASYEAR', 'PLOT',
                                                   'EVALID'))




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
      out <- parLapply(cl, X = names(plts), fun = dwmHelper1, plts,
                       db[names(db) %in% c('COND', 'COND_DWM_CALC')],
                       grpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = dwmHelper1, plts,
                      db[names(db) %in% c('COND', 'COND_DWM_CALC')],
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

    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)

    ## An iterator for population estimation
    popState <- split(pops, as.factor(pops$STATECD))
    #
    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        out <- parLapply(cl, X = names(popState), fun = dwmHelper2, popState, t, grpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = dwmHelper2, popState, t, grpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])





    ## Compute moving average weights if not TI ----------------------------------
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){

      ## Compute the weights
      wgts <- maWeights(pops, method, lambda)

      ## If moving average ribbons, add lambda to grpBy for easier summary
      if (str_to_upper(method) == 'EMA' & length(lambda) > 1){
        grpBy <- c('lambda', grpBy)
        aGrpBy <- c('lambda', aGrpBy)
      }


      ## Apply the weights
      if (str_to_upper(method) %in% c('LMA', 'EMA')){
        joinCols <- c('YEAR', 'STATECD', 'INVYR')
      } else {
        joinCols <- c('YEAR', 'STATECD')
      tEst <- tEst %>%
        left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
        left_join(wgts, by = joinCols) %>%
        mutate(across(aEst:cEst, ~(.*wgt))) %>%
        mutate(across(aVar:cvEst_c, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize(across(aEst:cvEst_c, sum, na.rm = TRUE))
      }




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
                           totals = FALSE,
                           variance = FALSE,
                           byPlot = FALSE,
                           tidy = TRUE,
                           nCores = 1) {

  ##  don't have to change original code
  grpBy_quo <- rlang::enquo(grpBy)
  areaDomain <- rlang::enquo(areaDomain)


  ## Handle iterator if db is remote
  remote <- ifelse(class(db) == 'Remote.FIA.Database', 1, 0)
  iter <- remoteIter(db, remote)

  ## Check for a most recent subset
  mr <- checkMR(db, remote)

  ## prep for areal summary
  polys <- arealSumPrep1(polys)



  ## Run the main portion
  out <- lapply(X = iter, FUN = dwmStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                landType, method,
                lambda, areaDomain,
                byPlot, totals, tidy,
                nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  tEst <- bind_rows(out[names(out) == 'tEst'])
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
               # variance totals
               AREA_TOTAL_VAR = aVar,
               VOL_DUFF_VAR = NA,
               VOL_LITTER_VAR = NA,
               VOL_1HR_VAR = vsmVar,
               VOL_10HR_VAR = vmdVar,
               VOL_100HR_VAR = vlgVar,
               VOL_1000HR_VAR = vcVar,
               VOL_PILE_VAR = vpVar,
               VOL_VAR = vVar,
               BIO_DUFF_VAR = bdVar,
               BIO_LITTER_VAR = blVar,
               BIO_1HR_VAR = bsmVar,
               BIO_10HR_VAR = bmdVar,
               BIO_100HR_VAR = blgVar,
               BIO_1000HR_VAR = bcVar,
               BIO_PILE_VAR = bpVar,
               BIO_VAR = bVar,
               CARB_DUFF_VAR = cdVar,
               CARB_LITTER_VAR = clVar,
               CARB_1HR_VAR = csmVar,
               CARB_10HR_VAR = cmdVar,
               CARB_100HR_VAR = clgVar,
               CARB_1000HR_VAR = ccVar,
               CARB_PILE_VAR = cpVar,
               CARB_VAR = cVar,

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
               # Per Acre variances -- for output
               VOL_DUFF_ACRE_VAR = NA,
               VOL_LITTER_ACRE_VAR = NA,
               VOL_1HR_ACRE_VAR = vsmVar,
               VOL_10HR_ACRE_VAR = vmdVar,
               VOL_100HR_ACRE_VAR = vlgVar,
               VOL_1000HR_ACRE_VAR = vcVar,
               VOL_PILE_ACRE_VAR = vpVar,
               VOL_ACRE_VAR = vVar,
               BIO_DUFF_ACRE_VAR = bdVar,
               BIO_LITTER_ACRE_VAR = blVar,
               BIO_1HR_ACRE_VAR = bsmVar,
               BIO_10HR_ACRE_VAR = bmdVar,
               BIO_100HR_ACRE_VAR = blgVar,
               BIO_1000HR_ACRE_VAR = bcVar,
               BIO_PILE_ACRE_VAR = bpVar,
               BIO_ACRE_VAR = bVar,
               CARB_DUFF_ACRE_VAR = cdVar,
               CARB_LITTER_ACRE_VAR = clVar,
               CARB_1HR_ACRE_VAR = csmVar,
               CARB_10HR_ACRE_VAR = cmdVar,
               CARB_100HR_ACRE_VAR = clgVar,
               CARB_1000HR_ACRE_VAR = ccVar,
               CARB_PILE_ACRE_VAR = cpVar,
               CARB_ACRE_VAR = cVar,

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
      if (variance){
        tOut <- tOut %>%
          select(grpBy, names(tOut)[str_detect(names(tOut), 'Var', negate = TRUE) &
                                      str_detect(names(tOut), 'cvEst', negate = TRUE) &
                                      str_detect(names(tOut), 'SE', negate = TRUE) &
                                      str_detect(names(tOut), 'Est', negate = TRUE)], nPlots_DWM, N)
      } else {
        tOut <- tOut %>%
          select(grpBy, names(tOut)[str_detect(names(tOut), 'Var', negate = TRUE) &
                                      str_detect(names(tOut), 'cvEst', negate = TRUE) &
                                      str_detect(names(tOut), 'VAR', negate = TRUE) &
                                      str_detect(names(tOut), 'Est', negate = TRUE)], nPlots_DWM, N)
      }

    } else {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, names(tOut)[str_detect(names(tOut), 'Var', negate = TRUE) &
                                      str_detect(names(tOut), 'cvEst', negate = TRUE) &
                                      str_detect(names(tOut), 'Est', negate = TRUE) &
                                      str_detect(names(tOut), 'SE', negate = TRUE) &
                                      str_detect(names(tOut), 'ACRE')], nPlots_DWM, N)
      } else {
        tOut <- tOut %>%
          select(grpBy, names(tOut)[str_detect(names(tOut), 'Var', negate = TRUE) &
                                      str_detect(names(tOut), 'cvEst', negate = TRUE) &
                                      str_detect(names(tOut), 'Est', negate = TRUE) &
                                      str_detect(names(tOut), 'VAR', negate = TRUE) &
                                      str_detect(names(tOut), 'ACRE')], nPlots_DWM, N)
      }
    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]


    ## Tidy things up if they didn't specify polys, returnSpatial
    if (tidy & is.null(polys) & returnSpatial == FALSE){

      ## Writing the variance options after the below, and don't want to change below
      ## So instead, do a temporary name swap
      if (variance) names(tOut) <- str_replace(names(tOut), '_VAR', '_SE')

      ## pivot longer
      bio <- pivot_longer(select(tOut, grpBy, BIO_DUFF_ACRE:BIO_ACRE, nPlots_DWM, N), names_to = 'FUEL_TYPE', values_to = 'BIO_ACRE', cols = BIO_DUFF_ACRE:BIO_ACRE) %>%
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
        select(grpBy, FUEL_TYPE, VOL_ACRE, BIO_ACRE, CARB_ACRE, VOL_ACRE_SE, BIO_ACRE_SE, CARB_ACRE_SE, nPlots_DWM, N) %>%
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
                  AREA_TOTAL_SE, nPlots_DWM, N) %>%
          filter(FUEL_TYPE %in% 'ACRE' == FALSE)
      }
      tOut <- fuel

      ## Writing the variance options after the below, and don't want to change below
      ## So instead, do a temporary name swap
      if (variance) names(tOut) <- str_replace(names(tOut), '_SE', '_VAR')

    }
  } # End byPlot

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

