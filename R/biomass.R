bioStarter <- function(x,
                       db,
                       grpBy_quo = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       bySpecies = FALSE,
                       bySizeClass = FALSE,
                       byComponent = FALSE,
                       landType = 'forest',
                       treeType = 'live',
                       method = 'TI',
                       lambda = .5,
                       treeDomain = NULL,
                       areaDomain = NULL,
                       totals = FALSE,
                       byPlot = FALSE,
                       component = 'AG',
                       bioMethod = 'CRM',
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
  if (any(str_to_upper(component) %in% c('AG', 'ROOTS', 'BOLE', "TOP", "FOLIAGE", "STUMP", "SAPLING", "WDLD_SPP", "TOTAL")) == FALSE) {
    stop('Unknown component. Must be a combination of: "AG", "ROOTS", "BOLE", "TOP", "FOLIAGE", "STUMP", "SAPLING", and "WDLD_SPP". Alternatively, use "TOTAL" for a sum of all components, and set "byComponent=TRUE" to estimate all components simultaneously.')
  }

  ## Biomass method warnings
  bioMethod <- stringr::str_to_upper(bioMethod)
  if (length(bioMethod) > 1) {
    stop('Only one "bioMethod" allowed at this time. Please choose one of "CRM" (Component Ratio Method) or "JENKINS" (Jenkins et al 2003).')
  }
  if (!c(bioMethod %in% c('CRM', 'JENKINS'))) {
    stop('"bioMethod" not recognized. Please choose one of "CRM" (Component Ratio Method) or "JENKINS" (Jenkins et al 2003).')
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


  ## When component = AG or total, replace with component names
  component = str_to_upper(unique(component))
  if ('TOTAL' %in% component | byComponent) {component <- c('ROOTS', 'BOLE', "TOP", "FOLIAGE", "STUMP", "SAPLING", "WDLD_SPP")}
  if ('AG' %in% component & byComponent == FALSE) {component <- unique(c(component[component != 'AG'], 'BOLE', "TOP", "STUMP", "SAPLING", "WDLD_SPP"))}



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



  ## estimate foliage biomass component here so we don't need to carry the extra columns
  db$TREE <- db$TREE %>%
    left_join(select(intData$REF_SPECIES_2018, -c(COMMON_NAME, GENUS, SPECIES)), by = 'SPCD') %>%
    mutate(jTotal = exp(JENKINS_TOTAL_B1 + JENKINS_TOTAL_B2 * log(DIA*2.54)) * 2.2046,
           stemRatio = exp(JENKINS_STEM_WOOD_RATIO_B1 + JENKINS_STEM_WOOD_RATIO_B2 / (DIA*2.54)),
           barkRatio = exp(JENKINS_STEM_BARK_RATIO_B1 + JENKINS_STEM_BARK_RATIO_B2 / (DIA*2.54)),
           leafRatio = exp(JENKINS_FOLIAGE_RATIO_B1 + JENKINS_FOLIAGE_RATIO_B2 / (DIA*2.54)),
           jBoleBio = (jTotal * stemRatio) + (jTotal * barkRatio),
           jLeafBio = jTotal * leafRatio,
           adj = case_when(is.na(DIA) ~ NA_real_,
                           !is.na(DRYBIO_WDLD_SPP) ~ DRYBIO_WDLD_SPP / (jTotal - jLeafBio),
                           DIA >= 5 ~ DRYBIO_BOLE / jBoleBio,
                           TRUE ~ JENKINS_SAPLING_ADJUSTMENT),
           DRYBIO_FOLIAGE = case_when(STATUSCD == 1 ~ jLeafBio * adj,
                                      STATUSCD == 2 ~ 0,
                                      TRUE ~ NA_real_))

  ## If we want to use Jenkin's methods instead of Component Ratio,
  ## swap out the bole and adjust other components
  if (bioMethod == 'JENKINS') {
    db$TREE <- db$TREE %>%
      ## Replacing component ratio biomass estimates w/ Jenkins
      mutate(DRYBIO_BOLE = jBoleBio,
             ## adj defined above - ratio of volume-based biomass estimates and
             ## diameter-based estimates for bole volume
             DRYBIO_TOP = DRYBIO_TOP / adj,
             DRYBIO_STUMP = DRYBIO_STUMP / adj,
             DRYBIO_BG = DRYBIO_BG / adj,
             DRYBIO_SAPLING = DRYBIO_SAPLING / adj,
             DRYBIO_WDLD_SPP = DRYBIO_WDLD_SPP / adj,
             DRYBIO_FOLIAGE = DRYBIO_WDLD_SPP / adj)
  }


  ## Add component to grpBy
  if (byComponent){
    grpBy <- c(grpBy, 'COMPONENT')
    grpByOrig <- c(grpByOrig, 'COMPONENT')
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
                               'SUBP', 'TREE', all_of(grpT), 'tD', 'typeD',
                               DRYBIO_TOP, DRYBIO_BOLE, DRYBIO_STUMP, DRYBIO_ROOTS = DRYBIO_BG,
                               DRYBIO_SAPLING, DRYBIO_WDLD_SPP, DRYBIO_FOLIAGE)) %>%
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
      out <- parLapply(cl, X = names(plts), fun = bioHelper1, plts,
                       db[names(db) %in% c('COND', 'TREE')],
                       grpBy, aGrpBy, byPlot, component)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = bioHelper1, plts,
                      db[names(db) %in% c('COND', 'TREE')],
                      grpBy, aGrpBy, byPlot, component, mc.cores = nCores)
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
        out <- parLapply(cl, X = names(popState), fun = bioHelper2, popState, a, t, grpBy, aGrpBy, method)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = bioHelper2, popState, a, t, grpBy, aGrpBy, method, mc.cores = nCores)
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
        group_by(ESTN_UNIT_CN, N, A, .dots = aGrpBy) %>%
        summarize(across(aEst:plotIn_AREA, sum, na.rm = TRUE))
      tEst <- tEst %>%
        left_join(select(db$POP_ESTN_UNIT, CN, STATECD), by = c('ESTN_UNIT_CN' = 'CN')) %>%
        left_join(wgts, by = joinCols) %>%
        mutate(across(bEst:cEst, ~(.*wgt))) %>%
        mutate(across(bVar:cvEst_c, ~(.*(wgt^2)))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize(across(bEst:plotIn_TREE, sum, na.rm = TRUE))



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
biomass <- function(db,
                    grpBy = NULL,
                    polys = NULL,
                    returnSpatial = FALSE,
                    bySpecies = FALSE,
                    bySizeClass = FALSE,
                    byComponent = FALSE,
                    landType = 'forest',
                    treeType = 'live',
                    method = 'TI',
                    lambda = .5,
                    treeDomain = NULL,
                    areaDomain = NULL,
                    totals = FALSE,
                    variance = FALSE,
                    byPlot = FALSE,
                    component = 'AG',
                    bioMethod = 'CRM',
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
  out <- lapply(X = iter, FUN = bioStarter, db,
                grpBy_quo = grpBy_quo, polys, returnSpatial,
                bySpecies, bySizeClass, byComponent,
                landType, treeType, method,
                lambda, treeDomain, areaDomain,
                totals, byPlot, component, bioMethod,
                nCores, remote, mr)
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
      summarize_all(sum,na.rm = TRUE)

    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)


    suppressWarnings({
      ## Bring them together
      tOut <- tTotal %>%
        left_join(aTotal, by = aGrpBy) %>%
        # Renaming, computing ratios, and SE
        mutate(AREA_TOTAL = aEst,
               BIO_TOTAL = bEst,
               CARB_TOTAL = cEst,
               ## Ratios
               BIO_ACRE = BIO_TOTAL / AREA_TOTAL,
               CARB_ACRE = CARB_TOTAL / AREA_TOTAL,
               ## Ratio Var
               baVar = (1/AREA_TOTAL^2) * (bVar + (BIO_ACRE^2 * aVar) - 2 * BIO_ACRE * cvEst_b),
               caVar = (1/AREA_TOTAL^2) * (cVar + (CARB_ACRE^2 * aVar) - 2 * CARB_ACRE * cvEst_c),
               ## SE RATIO
               BIO_ACRE_SE = sqrt(baVar) / BIO_ACRE * 100,
               CARB_ACRE_SE = sqrt(caVar) / CARB_ACRE *100,
               ## SE TOTAL
               AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL *100,
               BIO_TOTAL_SE = sqrt(bVar) / BIO_TOTAL *100,
               CARB_TOTAL_SE = sqrt(cVar) / CARB_TOTAL *100,

               ## VAR RATIO
               BIO_ACRE_VAR = baVar,
               CARB_ACRE_VAR = caVar,
               ## VAR TOTAL
               AREA_TOTAL_VAR = aVar,
               BIO_TOTAL_VAR = bVar,
               CARB_TOTAL_VAR = cVar,
               ## nPlots
               nPlots_TREE = plotIn_TREE,
               nPlots_AREA = plotIn_AREA)
    })


    if (totals) {

      if (variance){
        tOut <- tOut %>%
          select(grpBy, "BIO_ACRE","CARB_ACRE","BIO_TOTAL","CARB_TOTAL", "AREA_TOTAL", "BIO_ACRE_VAR",
                 "CARB_ACRE_VAR","BIO_TOTAL_VAR", "CARB_TOTAL_VAR",
                 "AREA_TOTAL_VAR","nPlots_TREE","nPlots_AREA", 'N')

      } else {
        tOut <- tOut %>%
          select(grpBy, "BIO_ACRE","CARB_ACRE","BIO_TOTAL","CARB_TOTAL", "AREA_TOTAL", "BIO_ACRE_SE",
                 "CARB_ACRE_SE","BIO_TOTAL_SE", "CARB_TOTAL_SE",
                 "AREA_TOTAL_SE","nPlots_TREE","nPlots_AREA", 'N')
      }



    } else {
      if (variance){
        tOut <- tOut %>%
          select(grpBy, "BIO_ACRE","CARB_ACRE", "BIO_ACRE_VAR",
                 "CARB_ACRE_VAR", "nPlots_TREE","nPlots_AREA", 'N')

      } else {
        tOut <- tOut %>%
          select(grpBy, "BIO_ACRE","CARB_ACRE", "BIO_ACRE_SE",
                 "CARB_ACRE_SE","nPlots_TREE","nPlots_AREA", 'N')
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









