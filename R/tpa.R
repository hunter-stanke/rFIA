tpaStarter <- function(x,
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
  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN',
                 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')

  ## If remote, read in state by state. Otherwise, drop all unnecessary tables
  db <- readRemoteHelper(x, db, remote, reqTables, nCores)

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


  ## Prep other variables ------------------------------------------------------
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>%
    mutate(PLT_CN = CN,
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
    db <- arealSumPrep2(db, grpBy, polys, nCores, remote)

    ## If there's nothing there, skip the state
    if (is.null(db)) return('no plots in polys')
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
  ## Tree type
  db$TREE$typeD <- treeTypeDomain(treeType,
                                  db$TREE$STATUSCD,
                                  db$TREE$DIA,
                                  db$TREE$TREECLCD)

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
  pops <- handlePops(db, evalType = c('EXPVOL', 'EXPCURR'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (str_to_upper(method) != 'TI') {
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
                               all_of(grpP), 'aD_p', 'sp', 'COUNTYCD'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS',
                               'COND_STATUS_CD', 'CONDID',
                               all_of(grpC), 'aD_c', 'landD')) %>%
    filter(PLT_CN %in% db$PLOT$PLT_CN)
  db$TREE <- select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ',
                               'SUBP', 'TREE', all_of(grpT), 'tD', 'typeD')) %>%
    filter(PLT_CN %in% db$PLOT$PLT_CN)

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
      out <- parLapply(cl, X = names(plts), fun = tpaHelper1, plts,
                       db[names(db) %in% c('COND', 'TREE')],
                       grpBy, aGrpBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = tpaHelper1, plts,
                      db[names(db) %in% c('COND', 'TREE')],
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
        out <- parLapply(cl, X = names(popState), fun = tpaHelper2, popState, a, t, grpBy, aGrpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = tpaHelper2, popState, a, t, grpBy, aGrpBy, method, mc.cores = nCores)
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
        mutate(across(tEst:bTEst, ~(.*wgt))) %>%
        mutate(across(tVar:cvEst_bT, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, A, .dots = grpBy) %>%
        summarize(across(tEst:cvEst_bT, sum, na.rm = TRUE))




    ## If using an ANNUAL estimator --------------------------------------------
    } else if (str_to_upper(method) == 'ANNUAL') {

      # If INVYR is in YEAR, choose the estimates when INVYR == YEAR
      # Otherwise, choose the estimates produced with the most plots
      aEst <- filterAnnual(aEst, aGrpBy, plotIn_AREA, db$POP_ESTN_UNIT)
      tEst <- filterAnnual(tEst, grpBy, plotIn_TREE, db$POP_ESTN_UNIT)
    }

    out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy, grpByOrig = grpByOrig)

  }

  return(out)

}














#' @export
tpa <- function(db,
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
  out <- lapply(X = iter, FUN = tpaStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                bySpecies, bySizeClass,
                landType, treeType, method,
                lambda, treeDomain, areaDomain,
                totals, byPlot, nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
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
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                aVar = sum(aVar, na.rm = TRUE),
                AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
                AREA_TOTAL_VAR = aVar,
                nPlots_AREA = sum(plotIn_AREA, na.rm = TRUE))

    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum, na.rm = TRUE) %>%
      mutate(TREE_TOTAL = tEst,
            BA_TOTAL = bEst,
            ## Variances
            treeVar = tVar,
            baVar = bVar,
            #aVar = first(aVar),
            cvT = cvEst_t,
            cvB = cvEst_b,
            ## Sampling Errors
            TREE_SE = sqrt(treeVar) / TREE_TOTAL * 100,
            TREE_VAR = treeVar,
            BA_SE = sqrt(baVar) / BA_TOTAL * 100,
            BA_VAR = baVar,
            #N = sum(N, na.rm = TRUE),
            nPlots_TREE = plotIn_TREE) %>%
      select(grpBy, TREE_TOTAL, BA_TOTAL, treeVar, baVar, cvT, cvB, TREE_SE, BA_SE, TREE_VAR, BA_VAR,
             nPlots_TREE, N, A)

    ## IF using polys, we treat each zone as a unique population
    if (!is.null(polys)){
      propGrp <- c('polyID', grpBy)
    } else {
      propGrp <- 'YEAR'
    }

    ## Hand the proportions
    tpTotal <- tEst %>%
      group_by(.dots = unique(propGrp)) %>%
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
        left_join(tpTotal, by = unique(propGrp)) %>%
        mutate(TPA = TREE_TOTAL / AREA_TOTAL,
               BAA = BA_TOTAL / AREA_TOTAL,
               tpaVar = (1/AREA_TOTAL^2) * (treeVar + (TPA^2 * aVar) - 2 * TPA * cvT),
               baaVar = (1/AREA_TOTAL^2) * (baVar + (BAA^2 * aVar) - (2 * BAA * cvB)),
               TPA_SE = sqrt(tpaVar) / TPA * 100,
               BAA_SE = sqrt(baaVar) / BAA * 100,
               TPA_VAR = tpaVar,
               BAA_VAR = baaVar,
               TPA_PERC = TREE_TOTAL / (TREE_TOTAL_full) * 100,
               BAA_PERC = BA_TOTAL / (BA_TOTAL_full) * 100,
               tpVar = (1/TREE_TOTAL_full^2) * (treeVar + (TPA_PERC^2 * tTVar) - 2 * TPA_PERC * cvTT),
               bpVar = (1/BA_TOTAL_full^2) * (baVar + (BAA_PERC^2 * bTVar) - (2 * BAA_PERC * cvBT)),
               TPA_PERC_SE = sqrt(tpVar) / TPA_PERC * 100,
               BAA_PERC_SE = sqrt(bpVar) / BAA_PERC * 100,
               TPA_PERC_VAR = tpVar,
               BAA_PERC_VAR = bpVar)


    })

    ## Naming conventions
    if (totals) {
      if (variance){
        tOut <- tTotal %>%
          select(grpBy, TPA, BAA, TPA_PERC, BAA_PERC, TREE_TOTAL, BA_TOTAL, AREA_TOTAL, TPA_VAR, BAA_VAR,
                 TPA_PERC_VAR, BAA_PERC_VAR, TREE_VAR, BA_VAR, AREA_TOTAL_VAR, nPlots_TREE, nPlots_AREA, N)
      } else {
        tOut <- tTotal %>%
          select(grpBy, TPA, BAA, TPA_PERC, BAA_PERC, TREE_TOTAL, BA_TOTAL, AREA_TOTAL, TPA_SE, BAA_SE,
                 TPA_PERC_SE, BAA_PERC_SE, TREE_SE, BA_SE, AREA_TOTAL_SE, nPlots_TREE, nPlots_AREA, N)
      }

    } else {
      if (variance) {
        tOut <- tTotal %>%
          select(grpBy, TPA, BAA, TPA_PERC, BAA_PERC,  TPA_VAR, BAA_VAR,
                 TPA_PERC_VAR, BAA_PERC_VAR, nPlots_TREE, nPlots_AREA, N)
      } else {
        tOut <- tTotal %>%
          select(grpBy, TPA, BAA, TPA_PERC, BAA_PERC,  TPA_SE, BAA_SE,
                 TPA_PERC_SE, BAA_PERC_SE, nPlots_TREE, nPlots_AREA, N)
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
  }

  return(tOut)
}





