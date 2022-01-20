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
                           treeList = FALSE,
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
    db$PLOT <- dplyr::filter(db$PLOT, prev == 0)
  }

  ## Handle TX issues - we only keep inventory years that are present in BOTH
  ## EAST AND WEST TX
  db <- handleTX(db)




  ## Some warnings if inputs are bogus -----------------------------------------
  if (!is.null(polys) &
      dplyr::first(class(polys)) %in%
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
  if (stringr::str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method,
                  'unknown. Defaulting to Temporally Indifferent (TI).'))
  }
  if (any(stringr::str_to_upper(component) %in% c('AG', 'ROOTS', 'BOLE', "TOP", "FOLIAGE", "STUMP", "SAPLING", "WDLD_SPP", "TOTAL")) == FALSE) {
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
    dplyr::mutate(PLT_CN = CN,
                  pltID = stringr::str_c(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))

  ## Convert grpBy to character
  grpBy <- grpByToChar(db, grpBy_quo)

  # I like a unique ID for a plot through time
  if (byPlot | treeList) {grpBy <- c('pltID', grpBy)}


  ## When component = AG or total, replace with component names
  component = stringr::str_to_upper(unique(component))
  if ('TOTAL' %in% component | byComponent) {component <- c('ROOTS', 'BOLE', "TOP", "FOLIAGE", "STUMP", "SAPLING", "WDLD_SPP")}
  if ('AG' %in% component & byComponent == FALSE) {component <- unique(c(component[component != 'AG'], 'BOLE', "TOP", "STUMP", "SAPLING", "WDLD_SPP"))}



  ## Intersect plots with polygons if polygons are given
  if (!is.null(polys)){

    ## Add shapefile names to grpBy
    grpBy = c(grpBy, names(polys)[names(polys) != 'geometry'])

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
  pops <- handlePops(db, evalType = c('VOL'), method, mr)

  ## A lot of states do their stratification in such a way that makes it impossible
  ## to estimate variance of annual panels w/ post-stratified estimator. That is,
  ## the number of plots within a panel within an stratum is less than 2. When
  ## this happens, merge strata so that all have at least two obs
  if (stringr::str_to_upper(method) != 'TI') {
    pops <- mergeSmallStrata(db, pops)
  }




  ## Canned groups -------------------------------------------------------------
  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      dplyr::left_join(dplyr::select(intData$REF_SPECIES_2018,
                                     c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      dplyr::mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' ')) %>%
      dplyr::mutate_if(is.factor,
                       as.character)
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }

  ## Break into size classes
  if (bySizeClass){
    grpBy <- c(grpBy, 'sizeClass')
    db$TREE$sizeClass <- makeClasses(db$TREE$DIA, interval = 2, numLabs = TRUE)
    db$TREE <- db$TREE[!is.na(db$TREE$sizeClass),]
  }



  ## estimate foliage biomass component here so we don't need to carry the extra columns
  db$TREE <- db$TREE %>%
    dtplyr::lazy_dt() %>%
    dplyr::left_join(dplyr::select(intData$REF_SPECIES_2018, -c(COMMON_NAME, GENUS, SPECIES)), by = 'SPCD') %>%
    dplyr::mutate(jTotal = exp(JENKINS_TOTAL_B1 + JENKINS_TOTAL_B2 * log(DIA*2.54)) * 2.2046,
                  stemRatio = exp(JENKINS_STEM_WOOD_RATIO_B1 + JENKINS_STEM_WOOD_RATIO_B2 / (DIA*2.54)),
                  barkRatio = exp(JENKINS_STEM_BARK_RATIO_B1 + JENKINS_STEM_BARK_RATIO_B2 / (DIA*2.54)),
                  leafRatio = exp(JENKINS_FOLIAGE_RATIO_B1 + JENKINS_FOLIAGE_RATIO_B2 / (DIA*2.54)),
                  jBoleBio = (jTotal * stemRatio) + (jTotal * barkRatio),
                  jLeafBio = jTotal * leafRatio,
                  adj = dplyr::case_when(is.na(DIA) ~ NA_real_,
                                         !is.na(DRYBIO_WDLD_SPP) ~ DRYBIO_WDLD_SPP / (jTotal - jLeafBio),
                                         DIA >= 5 ~ DRYBIO_BOLE / jBoleBio,
                                         TRUE ~ JENKINS_SAPLING_ADJUSTMENT),
                  DRYBIO_FOLIAGE = dplyr::case_when(STATUSCD == 1 ~ jLeafBio * adj,
                                                    STATUSCD == 2 ~ 0,
                                                    TRUE ~ NA_real_)) %>%
    as.data.frame()

  ## If we want to use Jenkin's methods instead of Component Ratio,
  ## swap out the bole and adjust other components
  if (bioMethod == 'JENKINS') {
    db$TREE <- db$TREE %>%
      ## Replacing component ratio biomass estimates w/ Jenkins
      dplyr::mutate(DRYBIO_BOLE = jBoleBio,
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
  }







  ## Prep the tree list --------------------------------------------------------

  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy &
                           !c(names(db$COND) %in% grpP)]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy &
                           !c(names(db$TREE) %in% c(grpP, grpC))]

  ## Dropping irrelevant rows and columns
  db$PLOT <- db$PLOT %>%
    dplyr::select(c(PLT_CN, STATECD, MACRO_BREAKPOINT_DIA,
                    INVYR, MEASYEAR, PLOT_STATUS_CD,
                    dplyr::all_of(grpP), sp, COUNTYCD)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    dplyr::filter(PLOT_STATUS_CD == 1 & sp == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)

  db$COND <- db$COND %>%
    dplyr::select(c(PLT_CN, CONDPROP_UNADJ, PROP_BASIS,
                    COND_STATUS_CD, CONDID,
                    dplyr::all_of(grpC), aD, landD)) %>%
    ## Drop non-forested plots, and those otherwise outside our domain of interest
    dplyr::filter(aD == 1 & landD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% pops$PLT_CN)

  db$TREE <- db$TREE %>%
    dplyr::select(c(PLT_CN, CONDID, DIA, SPCD, TPA_UNADJ,
                    SUBP, TREE, dplyr::all_of(grpT), tD, typeD,
                    DRYBIO_TOP, DRYBIO_BOLE, DRYBIO_STUMP, DRYBIO_ROOTS = DRYBIO_BG,
                    DRYBIO_SAPLING, DRYBIO_WDLD_SPP, DRYBIO_FOLIAGE)) %>%
    ## Drop plots outside our domain of interest
    dplyr::filter(!is.na(DIA) & TPA_UNADJ > 0 & tD == 1 & typeD == 1) %>%
    ## Drop visits not used in our eval of interest
    dplyr::filter(PLT_CN %in% db$PLOT$PLT_CN)

  # Separate area grouping names from tree grouping names
  if (!is.null(polys)){
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND), names(polys))]
  } else {
    aGrpBy <- grpBy[grpBy %in% c(names(db$PLOT), names(db$COND))]
  }



  ## Full tree list
  data <- db$PLOT %>%
    dplyr::left_join(db$COND, by = c('PLT_CN')) %>%
    dplyr::left_join(db$TREE, by = c('PLT_CN', 'CONDID'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp
  data$tDI <- data$landD * data$aD * data$tD * data$typeD * data$sp


  ## Convert to long format, where biomass component is the observation (multiple per tree)
  data <- data %>%
    tidyr::pivot_longer(cols = DRYBIO_TOP:DRYBIO_FOLIAGE,
                        names_to = c(".value", 'COMPONENT'),
                        names_sep = 7) %>%
    dplyr::rename(DRYBIO = DRYBIO_) %>%
    dplyr::filter(COMPONENT %in% component)



  ## Plot-level summaries ------------------------------------------------------
  if (byPlot & !treeList){

    grpBy <- c('YEAR', grpBy)
    grpSyms <- dplyr::syms(grpBy)
    aGrpSyms <- dplyr::syms(aGrpBy)

    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(PLT_CN, !!!aGrpSyms) %>%
      dplyr::summarize(PROP_FOREST = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE)) %>%
      as.data.frame()

    t <- data %>%
      dplyr::mutate(YEAR = MEASYEAR) %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, COMPONENT, .keep_all = TRUE) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(!!!grpSyms, PLT_CN) %>%
      dplyr::summarize(BIO_ACRE = sum(DRYBIO * TPA_UNADJ * tDI, na.rm = TRUE) / 2000) %>%
      dplyr::mutate(CARB_ACRE = BIO_ACRE * .5) %>%
      dplyr::left_join(a, by = c('PLT_CN', aGrpBy)) %>%
      as.data.frame()

    ## Make it spatial
    if (returnSpatial){
      t <- t %>%
        dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
        sf::st_as_sf(coords = c('LON', 'LAT'),
                     crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

    }

    out <- list(tEst = t, grpBy = grpBy, aGrpBy = aGrpBy)



    ## Population estimation -----------------------------------------------------
  } else {

    aGrpSyms <- dplyr::syms(aGrpBy)
    ### Condition list
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      dplyr::distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      dplyr::mutate(fa = CONDPROP_UNADJ * aDI) %>%
      dplyr::select(PLT_CN, AREA_BASIS = PROP_BASIS, CONDID, !!!aGrpSyms, PROP_BASIS, fa)

    grpSyms <- dplyr::syms(grpBy)
    ## Tree list
    t <- data %>%
      dplyr::distinct(PLT_CN, SUBP, TREE, COMPONENT, .keep_all = TRUE) %>%
      dplyr::mutate(bPlot = DRYBIO * TPA_UNADJ  * tDI / 2000,
                    cPlot = bPlot * .5) %>%
      ## Need a code that tells us where the tree was measured
      ## macroplot, microplot, subplot
      dplyr::mutate(
        TREE_BASIS = dplyr::case_when(
          ## When DIA is na, adjustment is NA
          is.na(DIA) ~ NA_character_,
          ## When DIA is less than 5", use microplot value
          DIA < 5 ~ 'MICR',
          ## When DIA is greater than 5", use subplot value
          DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
          DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
          DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) %>%
      dplyr::filter(!is.na(TREE_BASIS)) %>%
      dplyr::select(PLT_CN, TREE_BASIS, !!!grpSyms, SUBP, TREE, bPlot, cPlot)

    ## Return a tree/condition list ready to be handed to `customPSE`
    if (treeList) {

      tEst <- a %>%
        dplyr::left_join(t, by = c('PLT_CN', aGrpBy)) %>%
        ## Summarize over components so output isn't confusing to end user
        dtplyr::lazy_dt() %>%
        dplyr::mutate(EVAL_TYP = 'VOL') %>%
        dplyr::group_by(PLT_CN, EVAL_TYP, TREE_BASIS, AREA_BASIS,
                        !!!grpSyms, CONDID, SUBP, TREE, fa) %>%
        dplyr::summarise(dplyr::across(c(bPlot, cPlot), sum, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::select(PLT_CN, EVAL_TYP, TREE_BASIS, AREA_BASIS,
                      !!!grpSyms, CONDID, SUBP, TREE,
                      BIO_ACRE = bPlot,
                      CARB_ACRE = cPlot,
                      PROP_FOREST = fa) %>%
        as.data.frame()
      out <- list(tEst = tEst, aEst = NULL, grpBy = grpBy, aGrpBy = aGrpBy)

    ## Otherwise, proceed to population estimation
    } else {

      ## Sum variable(s) up to plot-level and adjust for non-response
      tPlt <- sumToPlot(t, pops, grpBy)
      aPlt <- sumToPlot(a, pops, aGrpBy)

      ## Adding YEAR to groups
      grpBy <- c('YEAR', grpBy)
      aGrpBy <- c('YEAR', aGrpBy)


      ## Sum variable(s) up to strata then estimation unit level
      eu.sums <- sumToEU(db, tPlt, aPlt, pops, grpBy, aGrpBy, method)
      tEst <- eu.sums$x
      aEst <- eu.sums$y


      out <- list(tEst = tEst, aEst = aEst, grpBy = grpBy, aGrpBy = aGrpBy)
    }



  }

  return(out)

}




#'@export
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
                    treeList = FALSE,
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
                totals, byPlot, treeList,
                component, bioMethod,
                nCores, remote, mr)
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  if (remote) out <- dropStatesOutsidePolys(out)
  aEst <- dplyr::bind_rows(out[names(out) == 'aEst'])
  tEst <- dplyr::bind_rows(out[names(out) == 'tEst'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  aGrpBy <- out[names(out) == 'aGrpBy'][[1]]
  grpSyms <- dplyr::syms(grpBy)
  aGrpSyms <- dplyr::syms(aGrpBy)

  ## Summarize population estimates across estimation units
  if (!byPlot & !treeList) {

    ## Combine most-recent population estimates across states with potentially
    ## different reporting schedules, i.e., if 2016 is most recent in MI and 2017 is
    ## most recent in WI, combine them and label as 2017
    if (mr) {
      tEst <- combineMR(tEst, grpBy)
      aEst <- combineMR(aEst, aGrpBy)
    }



    ## Totals and ratios -------------------------------------------------------
    aEst <- aEst %>%
      dplyr::group_by(!!!aGrpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!!aGrpSyms, fa_mean, fa_var, nPlots.y)

    tEst <- tEst %>%
      dplyr::group_by(!!!grpSyms) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), sum, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(aEst, by = aGrpBy) %>%
      # Renaming, computing ratios, and SE
      dplyr::mutate(BIO_TOTAL = bPlot_mean,
                    CARB_TOTAL = cPlot_mean,
                    AREA_TOTAL = fa_mean,

                    ## Ratios
                    BIO_ACRE = BIO_TOTAL / AREA_TOTAL,
                    CARB_ACRE = CARB_TOTAL / AREA_TOTAL,

                    ## Variances
                    BIO_TOTAL_VAR = bPlot_var,
                    CARB_TOTAL_VAR = cPlot_var,
                    AREA_TOTAL_VAR = fa_var,
                    BIO_ACRE_VAR = ratioVar(bPlot_mean, fa_mean, bPlot_var, fa_var, bPlot_cv),
                    CARB_ACRE_VAR = ratioVar(cPlot_mean, fa_mean, cPlot_var, fa_var, cPlot_cv),

                    ## Sampling errors
                    BIO_TOTAL_SE = sqrt(bPlot_var) / BIO_TOTAL * 100,
                    CARB_TOTAL_SE = sqrt(cPlot_var) / CARB_TOTAL * 100,
                    AREA_TOTAL_SE = sqrt(fa_var) / AREA_TOTAL * 100,
                    BIO_ACRE_SE = sqrt(BIO_ACRE_VAR) / BIO_ACRE * 100,
                    CARB_ACRE_SE = sqrt(CARB_ACRE_VAR) / CARB_ACRE * 100,

                    ## N plots
                    nPlots_TREE = nPlots.x,
                    nPlots_AREA = nPlots.y,
                    N = P2PNTCNT_EU) %>%
      dplyr::select(!!!grpSyms, BIO_ACRE, CARB_ACRE, BIO_TOTAL, CARB_TOTAL, AREA_TOTAL,
                    BIO_ACRE_VAR, CARB_ACRE_VAR, BIO_TOTAL_VAR, CARB_TOTAL_VAR, AREA_TOTAL_VAR,
                    BIO_ACRE_SE, CARB_ACRE_SE, BIO_TOTAL_SE, CARB_TOTAL_SE, AREA_TOTAL_SE,
                    nPlots_TREE, nPlots_AREA, N)

    ## Drop totals unless told not to
    if (!totals) {
      tEst <- tEst[,!stringr::str_detect(names(tEst), '_TOTAL')]
    }

    ## Select either variance or SE, depending on input
    if (variance) {
      tEst <- tEst[,!stringr::str_detect(names(tEst), '_SE')]
    } else {
      tEst <- tEst[,!stringr::str_detect(names(tEst), '_VAR')]
    }

  }


  ## Pretty output
  tEst <- tEst %>%
    dplyr::ungroup() %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    as_tibble()

  # We don't include YEAR in treeList output, and NA groups will be important
  # for retaining non-treed forestland
  if (!treeList) {
    tEst <- tEst %>%
      tidyr::drop_na(grpBy[!c(grpBy %in% names(polys))]) %>%
      dplyr::arrange(YEAR)
  }


  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]

  ## For spatial polygons
  if (returnSpatial & !byPlot) {
    tEst <- dplyr::left_join(tEst,
                             as.data.frame(dplyr::select(polys, polyID, geometry)),
                             by = 'polyID')
  }

  ## Above converts to tibble
  if (returnSpatial) tEst <- sf::st_sf(tEst)


  return(tEst)
}









