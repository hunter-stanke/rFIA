vrStarter <- function(x,
                      db,
                      grpBy_quo = NULL,
                      polys = NULL,
                      returnSpatial = FALSE,
                      bySpecies = FALSE,
                      bySizeClass = FALSE,
                      landType = 'forest',
                      treeType = 'live',
                      method = 'TI',
                      lambda = .5,
                      treeDomain = NULL,
                      areaDomain = NULL,
                      totals = FALSE,
                      byPlot = FALSE,
                      nCores = 1,
                      remote,
                      mr){





  ## Read required data, prep the database -------------------------------------
  reqTables <- c('PLOT', 'TREE', 'COND',
                 'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT', 'TREE_GRM_BEGIN',
                 'SUBP_COND_CHNG_MTRX',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')


  ## If remote, read in state by state. Otherwise, drop all unnecessary tables
  db <- readRemoteHelper(db, remote, reqTables, nCores)

  # ## IF the object was clipped
  # if ('prev' %in% names(db$PLOT)){
  #   ## Only want the current plots, no grm
  #   db$PLOT <- filter(db$PLOT, prev == 0)
  # }

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
  if (treeType %in% c('live', 'dead', 'gs', 'all') == FALSE){
    stop('treeType must be one of: "live", "dead", "gs", or "all".')
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
  ## No EXP_GROW available for most Western States, make sure we warn that values will be returned as 0
  # These states do not allow temporal queries. Things are extremely weird with their eval groups
  noGrow <- c(02,03,04,07,08,11,14,15,16, 30, 32, 35,43,49, 78)
  if(any(unique(db$PLOT$STATECD) %in% noGrow)){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    warning(paste('Growth data unavailable for: ', toString(fancyName) , '. Returning 0 for all growth estimates which include these states.', sep = ''))
  }

  # These states do not allow change estimates.
  if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    stop(paste('Growth & Mortality Estimates unavailable for: ', as.character(fancyName), sep = ''))
  }


  ## Prep other variables ------------------------------------------------------
  ## Need a plotCN, and a new ID
  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

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

  ## Land type and tree type combined
  db <- typeDomain_grow(db, treeType, landType, type = 'vr')

  ## Spatial boundary
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(!is.na(db$PLOT$polyID), 1, 0)
  } else {
    db$PLOT$sp <- 1
  }

  # User defined domain indicator for area (ex. specific forest type)
  db <- udAreaDomain(db, areaDomain)

  # User defined domain indicator for tree (ex. trees > 20 ft tall)
  db <- udTreeDomain(db, treeDomain)




  ## Handle population tables --------------------------------------------------
  ## Filtering out all inventories that are not relevant to the current estimation
  ## type. If using estimator other than TI, handle the differences in P2POINTCNT
  ## and in assigning YEAR column (YEAR = END_INVYR if method = 'TI')
  pops <- handlePops(db, evalType = c('EXPGROW'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (method != 'TI') {
    pops <- mergeSmallStrata(db, pops)
  }




  ## Canned groups -------------------------------------------------------------
  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      left_join(select(intData$REF_SPECIES_2018,
                       c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
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





  ## Slim down the database for we hand it off to the estimators ---------------
  ## Reduces memory requirements and speeds up processing ----------------------

  ## Only the necessary plots for EVAL of interest
  db$PLOT <- filter(db$PLOT, PLT_CN %in% pops$PLT_CN)

  ## Narrow up the tables to the necessary variables
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy &
                           !c(names(db$COND) %in% grpP)]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy &
                           !c(names(db$TREE) %in% c(grpP, grpC))]

  ### Only joining tables necessary to produce plot level estimates
  db$PLOT <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA',
                               'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD',
                               all_of(grpP), 'aD_p', 'sp', 'COUNTYCD',
                               REMPER, PREV_PLT_CN))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS',
                               'COND_STATUS_CD', 'CONDID',
                               all_of(grpC), 'aD_c', 'landD')) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  db$TREE <- select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN',
                               'PREV_TRE_CN', 'SUBP', 'TREE', all_of(grpT), 'tD',
                               'typeD', 'DIA', 'DRYBIO_AG', 'VOLCFNET',
                               'VOLCSNET')) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))
  db$TREE_GRM_COMPONENT <- db$TREE_GRM_COMPONENT %>%
    select(c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ',
             'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')) %>%
    filter(TRE_CN %in% db$TREE$TRE_CN)
  db$TREE_GRM_MIDPT <- db$TREE_GRM_MIDPT %>%
    select(c('TRE_CN', 'DIA', 'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG')) %>%
    filter(TRE_CN %in% db$TREE$TRE_CN)
  db$TREE_GRM_BEGIN <- db$TREE_GRM_BEGIN %>%
    select(c('TRE_CN', 'DIA', 'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG')) %>%
    filter(TRE_CN %in% db$TREE$TRE_CN)
  db$SUBP_COND_CHNG_MTRX <- db$SUBP_COND_CHNG_MTRX %>%
    select(PLT_CN, PREV_PLT_CN, SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID) %>%
    filter(PLT_CN %in% c(db$PLOT$PLT_CN, db$PLOT$PREV_PLT_CN))


  # Separate area grouping names from tree grouping names
  if (!is.null(polys)){
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND), names(polys))]
  } else {
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND))]
  }




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
      out <- parLapply(cl, X = names(plts), fun = vrHelper1, plts,
                       db[names(db) %in% c('COND', 'TREE', 'TREE_GRM_COMPONENT',
                                           'TREE_GRM_MIDPT', 'TREE_GRM_BEGIN',
                                           'SUBP_COND_CHNG_MTRX')],
                       grpBy, aGrpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = vrHelper1, plts,
                      db[names(db) %in% c('COND', 'TREE', 'TREE_GRM_COMPONENT',
                                          'TREE_GRM_MIDPT', 'TREE_GRM_BEGIN',
                                          'SUBP_COND_CHNG_MTRX')],
                      grpBy, aGrpBy, byPlot, mc.cores = nCores)
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

    out <- list(tEst = tOut, grpBy = grpBy, aGrpBy = aGrpBy, grpByOrig = grpByOrig)
    ## Population estimation
  } else {
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    a <- bind_rows(out[names(out) == 'a'])
    t <- bind_rows(out[names(out) == 't'])


    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)
    aGrpBy <- c('YEAR', aGrpBy)


    ## An iterator for population estimation
    popState <- split(pops, as.factor(pops$STATECD))
    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        out <- parLapply(cl, X = names(popState), fun = vrHelper2, popState, a, t, grpBy, aGrpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = vrHelper2, popState, a, t, grpBy, aGrpBy, method, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    aEst <- bind_rows(out[names(out) == 'aEst'])
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
      }
      aEst <- aEst %>%
        left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
        left_join(wgts, by = joinCols) %>%
        mutate(across(aEst, ~(.*wgt))) %>%
        mutate(across(aVar, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
        summarize(across(aEst:plotIn_AREA, sum, na.rm = TRUE))
      tEst <- tEst %>%
        left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
        left_join(wgts, by = joinCols) %>%
        mutate(across(tEst:bioAEst, ~(.*wgt))) %>%
        mutate(across(tVar:cvEst_bioA, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize(across(tEst:plotIn_TREE, sum, na.rm = TRUE))


      ## If using an ANNUAL estimator --------------------------------------------
    } else if (str_to_upper(method) == 'ANNUAL') {

      ## ANNUAL ESTIMATOR is when END_INVYR = INVYR
      aEst <- aEst %>%
        group_by(INVYR, .dots = aGrpBy) %>%
        summarize(across(.cols = everything(),  sum, na.rm = TRUE)) %>%
        filter(YEAR == INVYR) %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = INVYR)
      tEst <- tEst %>%
        group_by(INVYR, .dots = grpBy) %>%
        summarize(across(.cols = everything(),  sum, na.rm = TRUE)) %>%
        filter(YEAR == INVYR)%>%
        select(-c(YEAR)) %>%
        mutate(YEAR = INVYR)


      ## Rather than choose the annual panel estimate when INVYR = END_INVYR,
      ## choose the END_INVYR that has the highest N for each INVYR. Doing this
      ## because maybe not all 2018 data had been entered by the time the 2018
      ## END_INVYR cycle was produced. Maybe 2019 has more info on 2018 plots.
      ## So, ideally we would choose the cycle with the most plots for a given
      ## panel. Doing that here, important distinction from previous.
      ## NOT USED CURRENTLY --------------------------------------------------

      # aEst <- aEst %>%
      #   left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      #   group_by(STATECD, INVYR, .dots = aGrpBy[aGrpBy != 'STATECD']) %>%
      #   summarize(across(.cols = everything(),  sum, na.rm = TRUE)) %>%
      #   group_by(STATECD, INVYR, .dots = aGrpBy[aGrpBy %in% c('STATECD', 'YEAR') == FALSE]) %>%
      #   filter(plotIn_AREA == max(plotIn_AREA, na.rm = TRUE)) %>%
      #   filter(YEAR == max(YEAR, na.rm = TRUE)) %>%
      #   select(-c(YEAR)) %>%
      #   mutate(YEAR = INVYR)



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


    out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy, grpByOrig = grpByOrig)
  }

  return(out)

}



#' @export
vitalRates <- function(db,
                       grpBy = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       bySpecies = FALSE,
                       bySizeClass = FALSE,
                       landType = 'forest',
                       treeType = 'live',
                       method = 'TI',
                       lambda = .5,
                       treeDomain = NULL,
                       areaDomain = NULL,
                       totals = FALSE,
                       variance = FALSE,
                       byPlot = FALSE,
                       nCores = 1) {

  ##  don't have to change original code
  grpBy_quo <- rlang::enquo(grpBy)
  areaDomain <- rlang::enquo(areaDomain)
  treeDomain <- rlang::enquo(treeDomain)

  ## Handle iterator if db is remote
  remote <- ifelse(class(db) == 'Remote.FIA.Database', 1, 0)
  iter <- remoteIter(db, remote)

  ## Check for a most recent subset
  mr <- checkMR(db, remote)

  ## prep for areal summary
  polys <- arealSumPrep1(polys)


  ## Run the main portion
  out <- lapply(X = iter, FUN = vrStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                bySpecies, bySizeClass,
                landType, treeType, method,
                lambda, treeDomain, areaDomain,
                totals, byPlot, nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  aEst <- bind_rows(out[names(out) == 'aEst'])
  tEst <- bind_rows(out[names(out) == 'tEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  aGrpBy <- out[names(out) == 'aGrpBy'][[1]]
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
      aEst <- combineMR(aEst, aGrpBy)
    }



    ## Totals and ratios -------------------------------------------------------
    aTotal <- aEst %>%
      group_by(.dots = aGrpBy) %>%
      summarize_all(sum,na.rm = TRUE)

    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    suppressWarnings({
      ## Bring them together
      tOut <- tTotal %>%
        left_join(aTotal, by = aGrpBy) %>%
        # Renaming, computing ratios, and SE
        mutate(TREE_TOTAL = tEst,
               DIA_TOTAL = dEst,
               BA_TOTAL = bEst,
               NETVOL_TOTAL = gEst,
               BIO_TOTAL = bioEst,
               BAA_TOTAL = baaEst,
               NETVOLA_TOTAL = gaEst,
               BIOA_TOTAL = bioAEst,
               AREA_TOTAL = aEst,
               DIA_GROW = DIA_TOTAL / TREE_TOTAL,
               BA_GROW = BA_TOTAL / TREE_TOTAL,
               NETVOL_GROW = NETVOL_TOTAL / TREE_TOTAL,
               BIO_GROW = BIO_TOTAL / TREE_TOTAL,
               BA_GROW_AC = BAA_TOTAL / AREA_TOTAL,
               NETVOL_GROW_AC = NETVOLA_TOTAL / AREA_TOTAL,
               BIO_GROW_AC = BIOA_TOTAL / AREA_TOTAL,
               ## Ratio Var
               dgVar = (1/TREE_TOTAL^2) * (dVar + (DIA_GROW^2 * tVar) - 2 * DIA_GROW * cvEst_d),
               bgVar = (1/TREE_TOTAL^2) * (bVar + (BA_GROW^2 * tVar) - 2 * BA_GROW * cvEst_b),
               ggVar = (1/TREE_TOTAL^2) * (gVar + (NETVOL_GROW^2 * tVar) - 2 * NETVOL_GROW * cvEst_g),
               biogVar = (1/TREE_TOTAL^2) * (bioVar + (BIO_GROW^2 * tVar) - 2 * BIO_GROW * cvEst_bio),
               baagVar = (1/AREA_TOTAL^2) * (baaVar + (BA_GROW_AC^2 * aVar) - 2 * BA_GROW_AC * cvEst_baa),
               gagVar = (1/AREA_TOTAL^2) * (gaVar + (NETVOL_GROW_AC^2 * aVar) - 2 * NETVOL_GROW_AC * cvEst_ga),
               bioAgVar = (1/AREA_TOTAL^2) * (bioAVar + (BIO_GROW_AC^2 * aVar) - 2 * BIO_GROW_AC * cvEst_bioA),

               ## SE RATIO
               DIA_GROW_SE = sqrt(dgVar) / DIA_GROW * 100,
               BA_GROW_SE = sqrt(bgVar) / BA_GROW * 100,
               NETVOL_GROW_SE = sqrt(ggVar) / NETVOL_GROW * 100,
               BIO_GROW_SE = sqrt(biogVar) / BIO_GROW * 100,
               BA_GROW_AC_SE = sqrt(baagVar) / BA_GROW_AC * 100,
               NETVOL_GROW_AC_SE = sqrt(gagVar) / NETVOL_GROW_AC * 100,
               BIO_GROW_AC_SE = sqrt(bioAgVar) / BIO_GROW_AC * 100,
               ## VAR RATIO
               DIA_GROW_VAR = dgVar,
               BA_GROW_VAR = bgVar,
               NETVOL_GROW_VAR = ggVar,
               BIO_GROW_VAR = biogVar,
               BA_GROW_AC_VAR = baagVar,
               NETVOL_GROW_AC_VAR = gagVar,
               BIO_GROW_AC_VAR = bioAgVar,
               ## SE TOTAL
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL *100,
               TREE_TOTAL_SE = sqrt(tVar) / TREE_TOTAL *100,
               DIA_TOTAL_SE = sqrt(dVar) / DIA_TOTAL *100,
               BA_TOTAL_SE = sqrt(bVar) / BA_TOTAL *100,
               NETVOL_TOTAL_SE = sqrt(gVar) / NETVOL_TOTAL *100,
               BIO_TOTAL_SE = sqrt(bioVar) / BIO_TOTAL *100,
               BAA_TOTAL_SE = sqrt(baaVar) / BAA_TOTAL *100,
               NETVOLA_TOTAL_SE = sqrt(gaVar) / NETVOLA_TOTAL *100,
               BIOA_TOTAL_SE = sqrt(bioAVar) / BIOA_TOTAL *100,
               ## VAR TOTAL
               AREA_TOTAL_VAR = aVar,
               TREE_TOTAL_VAR = tVar,
               DIA_TOTAL_VAR = dVar,
               BA_TOTAL_VAR = bVar,
               NETVOL_TOTAL_VAR = gVar,
               BIO_TOTAL_VAR = bioVar,
               BAA_TOTAL_VAR = baaVar,
               NETVOLA_TOTAL_VAR = gaVar,
               BIOA_TOTAL_VAR = bioAVar,
               ## nPlots
               # Non-zero plots
               nPlots_TREE = plotIn_TREE,
               nPlots_AREA = plotIn_AREA)
    })


    # Make some columns go away
    if (totals) {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, DIA_GROW, BA_GROW, NETVOL_GROW, BIO_GROW, BA_GROW_AC, NETVOL_GROW_AC,
                 BIO_GROW_AC, DIA_GROW_VAR, BA_GROW_VAR, NETVOL_GROW_VAR, BIO_GROW_VAR,
                 BA_GROW_AC_VAR, NETVOL_GROW_AC_VAR, BIO_GROW_AC_VAR,
                 TREE_TOTAL, DIA_TOTAL, BA_TOTAL, NETVOL_TOTAL, BIO_TOTAL, BAA_TOTAL,
                 NETVOLA_TOTAL,  BIOA_TOTAL, AREA_TOTAL,
                 TREE_TOTAL_VAR, DIA_TOTAL_VAR, BA_TOTAL_VAR, NETVOL_TOTAL_VAR, BIO_TOTAL_VAR, BAA_TOTAL_VAR,
                 NETVOLA_TOTAL_VAR,  BIOA_TOTAL_VAR, AREA_TOTAL_VAR, nPlots_TREE, nPlots_AREA, N)
      } else {
        tOut <- tOut %>%
          select(grpBy, DIA_GROW, BA_GROW, NETVOL_GROW, BIO_GROW, BA_GROW_AC, NETVOL_GROW_AC,
                 BIO_GROW_AC, DIA_GROW_SE, BA_GROW_SE, NETVOL_GROW_SE, BIO_GROW_SE,
                 BA_GROW_AC_SE, NETVOL_GROW_AC_SE, BIO_GROW_AC_SE,
                 TREE_TOTAL, DIA_TOTAL, BA_TOTAL, NETVOL_TOTAL, BIO_TOTAL, BAA_TOTAL,
                 NETVOLA_TOTAL,  BIOA_TOTAL, AREA_TOTAL,
                 TREE_TOTAL_SE, DIA_TOTAL_SE, BA_TOTAL_SE, NETVOL_TOTAL_SE, BIO_TOTAL_SE, BAA_TOTAL_SE,
                 NETVOLA_TOTAL_SE,  BIOA_TOTAL_SE, AREA_TOTAL_SE, nPlots_TREE, nPlots_AREA)
      }

    } else {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, DIA_GROW, BA_GROW, NETVOL_GROW, BIO_GROW, BA_GROW_AC, NETVOL_GROW_AC,
                 BIO_GROW_AC, DIA_GROW_VAR, BA_GROW_VAR, NETVOL_GROW_VAR, BIO_GROW_VAR,
                 BA_GROW_AC_VAR, NETVOL_GROW_AC_VAR, BIO_GROW_AC_VAR, nPlots_TREE, nPlots_AREA, N)
      } else {
        tOut <- tOut %>%
          select(grpBy, DIA_GROW, BA_GROW, NETVOL_GROW, BIO_GROW, BA_GROW_AC, NETVOL_GROW_AC,
                 BIO_GROW_AC, DIA_GROW_SE, BA_GROW_SE, NETVOL_GROW_SE, BIO_GROW_SE,
                 BA_GROW_AC_SE, NETVOL_GROW_AC_SE, BIO_GROW_AC_SE, nPlots_TREE, nPlots_AREA)
      }

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
  } else {
    ## Sometimes we see blanks for non EXPGROW years
    tOut <- filter(tOut, nPlots_AREA > 0)
  }

  return(tOut)
}


