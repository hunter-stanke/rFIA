
#### A Start up Message ------------
.onAttach <- function(lib, pkg) {
  if(interactive() || getOption("verbose"))
    packageStartupMessage(sprintf("Package %s (%s) loaded. Check out our website at https://rfia.netlify.app/.\nType citation(\"%s\") for examples of how to cite rFIA.\n", pkg,
                                  packageDescription(pkg)$Version, pkg))
}



#' @import dtplyr
#' @import dplyr
#' @import methods
#' @import sf
#' @import stringr
#' @import bit64
#' @import tidyselect
#' @importFrom rlang eval_tidy enquo enquos quos quo
#' @importFrom data.table fread fwrite rbindlist fifelse
#' @importFrom parallel makeCluster detectCores mclapply parLapply stopCluster clusterEvalQ
#' @import tidyr
#' @import ggplot2
#' @importFrom stats cov var coef lm na.omit quantile
#' @importFrom utils adist object.size read.csv tail globalVariables type.convert download.file unzip packageDescription
NULL


## Iterator if db is remote
remoteIter <- function(db, remote){
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

  return(iter)
}

## Check most recent and handle remote dbs
checkMR <- function(db, remote){
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
  return(mr)
}

## Prep database for areal summary
arealSumPrep1 <- function(polys){

  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    ## A unique ID
    polys$polyID <- 1:nrow(polys)
  }

  return(polys)
}

## Do the spatial intersection of plots w/ polgyons
arealSumPrep2 <- function(db, grpBy, polys, nCores){

  ## Make plot data spatial, projected same as polygon layer
  pltSF <- select(db$PLOT, c('LON', 'LAT', pltID)) %>%
    filter(!is.na(LAT) & !is.na(LON)) %>%
    distinct(pltID, .keep_all = TRUE) %>%
    st_as_sf(coords = c('LON', 'LAT'))
  st_crs(pltSF) <- 4326
  pltSF <- pltSF %>%
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

  return(db)

}

## Helper function to read remote database by state within
## estimator functions
readRemoteHelper <- function(x, db, remote, reqTables, nCores){
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

  return(db)
}


## Convert grpBy from quos to character and make sure columns exist
grpByToChar <- function(db, grpBy_quo){

  ## If tree table exists
  if ('TREE' %in% names(db)){
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
  } else {
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
  }

  return(grpBy)

}

## Drop all inventories that are specific to east or west tx. Only retain evals
## that span the entire state
handleTX <- function(db){

  if (any(db$POP_EVAL$STATECD %in% 48)){

    badIDS <- db$POP_PLOT_STRATUM_ASSGN %>%
      filter(STATECD %in% 48) %>%
      distinct(EVALID, UNITCD) %>%
      group_by(EVALID) %>%
      summarise(n = n()) %>%
      filter(n != 7) ## i.e., only keep EVALIDs w/out all 7 units

    db$POP_EVAL <- filter(db$POP_EVAL, !c(EVALID %in% badIDS))

  }

  return(db)
}

## Land type domain indicator
landTypeDomain <- function(landType, COND_STATUS_CD, SITECLCD, RESERVCD) {
  if (tolower(landType) == 'forest'){
    landD <- ifelse(COND_STATUS_CD == 1, 1, 0)
  } else if (tolower(landType) == 'timber'){
    landD <- ifelse(COND_STATUS_CD == 1 & SITECLCD %in% c(1, 2, 3, 4, 5, 6) & RESERVCD == 0, 1, 0)
  } else if (tolower(landType) == 'non-forest'){
    landD <- ifelse(COND_STATUS_CD == 2)
  } else if (tolower(landType) == 'water'){
    landD <- ifelse(COND_STATUS_CD == 3 | COND_STATUS_CD == 4)
  } else if (tolower(landType) == 'census water'){
    landD <- ifelse(COND_STATUS_CD == 4)
  } else if (tolower(landType) == 'non-census water'){
    landD <- ifelse(COND_STATUS_CD == 4)
  } else if (tolower(landType) == 'all') {
    landD <- 1
  }

  return(landD)
}

## Tree type domain indicator
treeTypeDomain <- function(treeType, STATUSCD, DIA, TREECLCD) {
  if (tolower(treeType) == 'live'){
    typeD <- ifelse(STATUSCD == 1, 1, 0)
  } else if (tolower(treeType) == 'dead'){
    typeD <- ifelse(STATUSCD == 2, 1, 0)
  } else if (tolower(treeType) == 'gs'){
    typeD <- ifelse(STATUSCD == 1 & DIA >= 5 & TREECLCD == 2, 1, 0)
  } else if (tolower(treeType) == 'all'){
    typeD <- 1
  }
  return(typeD)
}

typeDomain_grow <- function(db, treeType, landType, type) {

  if (type == 'vr'){
    if (tolower(landType) == 'forest'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
      # Tree Type domain indicator
      if (tolower(treeType) == 'live'){
        db$TREE$typeD <- 1
        ## Rename some variables in grm
        db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                        TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_FOREST,
                                        TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_FOREST,
                                        TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_FOREST,
                                        SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_FOREST,
                                        COMPONENT = SUBP_COMPONENT_AL_FOREST)

      } else if (tolower(treeType) == 'gs'){
        db$TREE$typeD <- 1 #ifelse(db$TREE$DIA >= 5, 1, 0)
        db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                        TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_FOREST,
                                        TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_FOREST,
                                        TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_FOREST,
                                        SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_FOREST,
                                        COMPONENT = SUBP_COMPONENT_GS_FOREST)
      }
    } else if (tolower(landType) == 'timber'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
      # Tree Type domain indicator
      if (tolower(treeType) == 'live'){
        db$TREE$typeD <- 1
        ## Rename some variables in grm
        db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                        TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_TIMBER,
                                        TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_TIMBER,
                                        TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_TIMBER,
                                        SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_TIMBER,
                                        COMPONENT = SUBP_COMPONENT_AL_TIMBER)

      } else if (tolower(treeType) == 'gs'){
        db$TREE$typeD <- 1 #ifelse(db$TREE$DIA >= 5, 1, 0)
        db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                        TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_TIMBER,
                                        TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_TIMBER,
                                        TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_TIMBER,
                                        SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_TIMBER,
                                        COMPONENT = SUBP_COMPONENT_GS_TIMBER)
      }
    }
  } else if (type == 'gm') {
    ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
    # Land type domain indicator
    if (tolower(landType) == 'forest'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
      # Tree Type domain indicator
      if (tolower(treeType) == 'all'){
        db$TREE$typeD <- 1
        ## Rename some variables in grm
        db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                        TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_FOREST,
                                        TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_FOREST,
                                        TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_FOREST,
                                        SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_FOREST,
                                        COMPONENT = SUBP_COMPONENT_AL_FOREST) %>%
          mutate(TPARECR_UNADJ = case_when(
            is.na(COMPONENT) ~ NA_real_,
            COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2') ~ TPAGROW_UNADJ,
            TRUE ~ 0))

      } else if (tolower(treeType) == 'gs'){
        # db$TREE <- db$TREE %>%
        #   mutate(typeD = case_when(
        #     STATUSCD %in% 1:2 & DIA >=5 ~ 1,
        #     STATUSCD == 3 & PREVDIA >=5 ~ 1,
        #     TRUE ~ 0))
        db$TREE$typeD <- 1
        db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                        TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_FOREST,
                                        TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_FOREST,
                                        TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_FOREST,
                                        SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_FOREST,
                                        COMPONENT = SUBP_COMPONENT_GS_FOREST)%>%
          mutate(TPARECR_UNADJ = case_when(
            is.na(COMPONENT) ~ NA_real_,
            COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2') ~ TPAGROW_UNADJ,
            TRUE ~ 0))
      }
    } else if (tolower(landType) == 'timber'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
      # Tree Type domain indicator
      if (tolower(treeType) == 'all'){
        db$TREE$typeD <- 1
        ## Rename some variables in grm
        db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                        TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_TIMBER,
                                        TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_TIMBER,
                                        TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_TIMBER,
                                        SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_TIMBER,
                                        COMPONENT = SUBP_COMPONENT_AL_TIMBER)%>%
          mutate(TPARECR_UNADJ = case_when(
            is.na(COMPONENT) ~ NA_real_,
            COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2') ~ TPAGROW_UNADJ,
            TRUE ~ 0))

      } else if (tolower(treeType) == 'gs'){
        # db$TREE <- db$TREE %>%
        #   mutate(typeD = case_when(
        #     STATUSCD %in% 1:2 & DIA >=5 ~ 1,
        #     STATUSCD == 3 & PREVDIA >=5 ~ 1,
        #     TRUE ~ 0))
        db$TREE$typeD <- 1
        db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                        TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_TIMBER,
                                        TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_TIMBER,
                                        TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_TIMBER,
                                        SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_TIMBER,
                                        COMPONENT = SUBP_COMPONENT_GS_TIMBER)%>%
          mutate(TPARECR_UNADJ = case_when(
            is.na(COMPONENT) ~ NA_real_,
            COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2') ~ TPAGROW_UNADJ,
            TRUE ~ 0))
      }
    }
  }



  return(db)
}

## Build domain indicator for UD area domain
udAreaDomain <- function(db, areaDomain) {

  pcEval <- left_join(db$PLOT, select(db$COND, -c('STATECD', 'UNITCD', 'COUNTYCD', 'INVYR', 'PLOT')), by = 'PLT_CN')
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

  return(db)
}

## Build domain indicator for UD area domain
udTreeDomain <- function(db, treeDomain) {

  tD <- rlang::eval_tidy(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$TREE$tD <- as.numeric(tD)

  return(db)
}

## Build domain indicator for UD area domain
udSeedDomain <- function(db, treeDomain) {

  tD <- rlang::eval_tidy(treeDomain, db$SEEDLING) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$SEEDLING$tD <- as.numeric(tD)

  return(db)
}



## Handle population tables
handlePops <- function(db, evalType, method, mr, pltList = NULL, ga = FALSE){

  if (ga){
    ## What years are growth accounting years --> not all filled in
    ga <- db$POP_EVAL %>%
      group_by(END_INVYR) %>%
      summarize(ga = if_else(any(GROWTH_ACCT == 'Y'), 1, 0))

    ### Snag the EVALIDs that are needed
    db$POP_EVAL  <- db$POP_EVAL %>%
      left_join(ga, by = 'END_INVYR') %>%
      select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD', 'GROWTH_ACCT', 'ga', STATECD) %>%
      left_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID, .keep_all = TRUE)
    gaVars <- c('GROWTH_ACCT')
  } else {
    ### Snag the EVALIDs that are needed
    db$POP_EVAL<- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD', STATECD) %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP %in% evalType) %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID, .keep_all = TRUE)
    gaVars <- NULL
  }


  ## If a most-recent subset, make sure that we don't get two reporting years in
  ## western states
  if (mr) {
    db$POP_EVAL <- db$POP_EVAL %>%
      group_by(EVAL_TYP, STATECD) %>%
      filter(END_INVYR == max(END_INVYR, na.rm = TRUE)) %>%
      ungroup()
  }


  ## Cut STATECD
  db$POP_EVAL <- select(db$POP_EVAL, -c(STATECD))

  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'END_INVYR', 'EVAL_TYP', any_of(gaVars))) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    distinct(EVALID, ESTN_UNIT_CN, STRATUM_CN, .keep_all = TRUE) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
    distinct(EVALID, ESTN_UNIT_CN, STRATUM_CN, PLT_CN, .keep_all = TRUE) %>%
    ungroup() %>%
    mutate_if(is.factor,
              as.character) %>%
    rename(YEAR = END_INVYR)

  ## Dropping non-sampled plots for P3 variables
  if (!is.null(pltList)) {
    pops <- pops %>%
      filter(pops$PLT_CN %in% pltList)
  }

  ## P2POINTCNT column is NOT consistent for annual estimates, plots
  ## within individual strata and est units are related to different INVYRs
  pops <- pops %>%
    ## Count of plots by Stratum/INVYR
    group_by(EVALID, ESTN_UNIT_CN, STRATUM_CN, INVYR) %>%
    mutate(P2POINTCNT_INVYR = n()) %>%
    ## By unit / INVYR
    group_by(EVALID, ESTN_UNIT_CN, INVYR) %>%
    mutate(p2eu_INVYR = n()) %>%
    ## By unit for entire cycle
    group_by(EVALID, ESTN_UNIT_CN) %>%
    mutate(p2eu = n()) %>%
    ungroup()

  ## Recode a few of the estimation methods to make things easier below
  pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')

  return(pops)

}



## When something other than temporally indifferent is used, we may need to merge small strata
## There is no great way to go about this, but very important we do it for variance issues.
## So, what we will do is:
## (1) Identify strata/INVYR pairs with less than 2 ground plots
## (2) For each of those pairs, identify their most similar neighbor based on fuzzy string matching of STRATUM Descriptions.
##     Neighbors (other strata) must be from the same estimation unit and measured in the same year (INVYR). If no
##     neighbors exist, i.e., the small stratum was the only one measured in a given INVYR
## Important -- we are effectively allowing for stratification to vary across panels within an inventory cycle.
##              Also, we are adjusting strata weights by INVYR within the cycle, i.e., the same stratum may be
##              weighted differently in different years within the same cycle

mergeSmallStrata <- function(db, pops) {

  ## Make a unique ID for stratum / year pairs
  pops$stratID <- paste(pops$STRATUM_CN, pops$INVYR, sep = '_')
  ## We'll allow strata weights to vary by INVYR w/ same strata because not all
  ## strata are sampled in a given year
  pops$P1POINTCNT_INVYR <- pops$P1POINTCNT

  ## Stratum year pairs
  stratYr <- pops %>%
    left_join(select(db$POP_STRATUM, CN, STRATUM_DESCR), by = c('STRATUM_CN' = 'CN')) %>%
    distinct(STATECD, ESTN_UNIT_CN, STRATUM_CN, STRATUM_DESCR, INVYR,
             P2POINTCNT, P1POINTCNT, P2POINTCNT_INVYR, P1POINTCNT_INVYR, stratID) %>%
    ## If buffer is present in the name, then the stratum has a different intensity
    ## than other strata in the same estimation unit (PNW only).
    ## Only combine buffer w/ buffer
    mutate(buff = str_detect(STRATUM_DESCR, 'buff') & STATECD %in% c(53, 41, 6)) %>%
    mutate(wrong = P2POINTCNT_INVYR < 2) %>%
    group_by(ESTN_UNIT_CN, INVYR) %>%
    mutate(nStrata_INVYR = length(unique(STRATUM_CN))) %>%
    group_by(ESTN_UNIT_CN) %>%
    mutate(nStrata = length(unique(STRATUM_CN))) %>%
    ungroup() %>%
    arrange(P2POINTCNT_INVYR)

  ## Check if any fail
  warnMe <- c()

  ## If any are too small, i.e., only one plot --> do some merging
  if (sum(stratYr$wrong, na.rm = TRUE) > 0){

    for ( i in stratYr$stratID[stratYr$wrong == 1] ) {

      ## Subset the row
      dat <- filter(stratYr, stratID == i)


      ## Use fuzzy string matching if any other strata are available to merge with
      if (dat$nStrata_INVYR > 1) {
        ## Find its nearest neighbor of those in the same estimation
        ## unit and INVYR
        neighbors <- stratYr %>%
          filter(ESTN_UNIT_CN == dat$ESTN_UNIT_CN) %>%
          #filter(buff == dat$buff) %>%
          filter(INVYR == dat$INVYR) %>%
          filter(stratID != i)

        if (nrow(neighbors) < 1) {
          warnMe <- c(warnMe, TRUE)

        } else {
          warnMe <- c(warnMe, FALSE)

          ## Find the most similar neighbor in terms of stratum description
          msn <- adist(dat$STRATUM_DESCR, neighbors$STRATUM_DESCR)
          msnID <- neighbors$stratID[which.min(msn)]


          ## In pops, we want to update all rows of the giving and receiving strata
          ## Where giving gets a change in STRATUM_CN, P1POINTCNT, and P2POINTCNT

          ## Giving stratum ----------------------------------------------------
          pops[pops$stratID == i, 'STRATUM_CN'] <- unique(pops[pops$stratID == msnID, 'STRATUM_CN'])
          pops[pops$stratID == i, 'P1POINTCNT_INVYR'] <- unique(pops[pops$stratID == msnID, 'P1POINTCNT_INVYR']) + dat$P1POINTCNT_INVYR
          pops[pops$stratID == i, 'P2POINTCNT_INVYR'] <- unique(pops[pops$stratID == msnID, 'P2POINTCNT_INVYR']) + dat$P2POINTCNT_INVYR
          #pops[pops$stratID == i, 'stratID'] <- unique(pops[pops$stratID == msnID, 'stratID'])


          ## Receiving stratum -------------------------------------------------
          pops[pops$stratID == msnID, 'P1POINTCNT_INVYR'] <- unique(pops[pops$stratID == msnID, 'P1POINTCNT_INVYR']) + dat$P1POINTCNT_INVYR
          pops[pops$stratID == msnID, 'P2POINTCNT_INVYR'] <- unique(pops[pops$stratID == msnID, 'P2POINTCNT_INVYR']) + dat$P2POINTCNT_INVYR

        }




        ## If the small stratum is the only one available in the estimation unit in a given year,
        ## merge the INVYR within the same stratum
      } else {

        ## No other strata measured in the same year, so merge years instead
        neighbors <- stratYr %>%
          filter(STRATUM_CN == dat$STRATUM_CN) %>%
          filter(stratID != i)

        if (nrow(neighbors) > 0) {
          warnMe <- c(warnMe, FALSE)

          ## Find the most similar neighbor in terms of stratum description
          msn <- abs(neighbors$INVYR - (dat$INVYR + .01))
          msnID <- neighbors$stratID[which.min(msn)]

          ## Giving stratum ----------------------------------------------------
          pops[pops$stratID == i, 'P2POINTCNT_INVYR'] <- unique(pops[pops$stratID == msnID, 'P2POINTCNT_INVYR']) + dat$P2POINTCNT_INVYR
          pops[pops$stratID == i, 'INVYR'] <- unique(pops[pops$stratID == msnID, 'INVYR'])
          #pops[pops$stratID == i, 'stratID'] <- unique(pops[pops$stratID == msnID, 'stratID'])


          ## Receiving stratum -------------------------------------------------
          pops[pops$stratID == msnID, 'P2POINTCNT_INVYR'] <- unique(pops[pops$stratID == msnID, 'P2POINTCNT_INVYR']) + dat$P2POINTCNT_INVYR

        } else {
          warnMe <- c(warnMe, TRUE)
        }
      }
    }

    ## Keeps it from repeating above
    if (any(warnMe)) warning('Bad stratification, i.e., strata too small to compute variance of annual panels. If you are only interested in totals and/or ratio estimates, disregard this. However, if interested in variance (e.g., for confidence intervals) try using method = "TI".')


    ## Update adjustment factors in estimation unit
    pops <- pops %>%
      distinct(EVALID, ESTN_UNIT_CN, STRATUM_CN, INVYR, PLT_CN, .keep_all = TRUE) %>%
      ## Fix adjustment factors
      group_by(STRATUM_CN) %>%
      mutate(ADJ_FACTOR_MICR = mean(ADJ_FACTOR_MICR, na.rm = TRUE),
             ADJ_FACTOR_SUBP = mean(ADJ_FACTOR_SUBP, na.rm = TRUE),
             ADJ_FACTOR_MACR = mean(ADJ_FACTOR_MACR)) %>%
      ungroup()
  }


  ## Whether any strata are small are not, we will almost surely need to adjust
  ## strata weights for some INVYRs. Not all stratum are measured in each panel
  ## within a cycle, and when this happens, we need to weight the observed strata
  ## higher. If we don't, strata weights will not sum to 1 in some years and we
  ## will grossly underestimate means/totals in some cases
  ## Adjust stratum weights when not all strata are sampled in an INVYR
  stratYr <- pops %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, INVYR,
             P1POINTCNT_INVYR, P1PNTCNT_EU,
             P2POINTCNT_INVYR) %>%
    ## If multiple stratum were merged onto one, choose the maximum value (i.e.,
    ## the result of the last iteration)
    group_by(ESTN_UNIT_CN, STRATUM_CN, INVYR) %>%
    mutate(P1POINTCNT_INVYR = max(P1POINTCNT_INVYR),
           P2POINTCNT_INVYR = max(P2POINTCNT_INVYR)) %>%
    distinct() %>%
    mutate(stratWgt = P1POINTCNT_INVYR / P1PNTCNT_EU) %>%
    group_by(ESTN_UNIT_CN, INVYR) %>%
    mutate(propSampled = sum(stratWgt, na.rm = TRUE),
           stratWgt_INVYR = stratWgt / propSampled,
           P1POINTCNT_INVYR = P1PNTCNT_EU * stratWgt_INVYR,
           p2eu_INVYR = sum(P2POINTCNT_INVYR, na.rm = TRUE)) %>%
    ungroup() %>%
    select(STRATUM_CN, INVYR, P1POINTCNT_INVYR, P2POINTCNT_INVYR, p2eu_INVYR)

  pops <- pops %>%
    select(-c(P1POINTCNT_INVYR, stratID, P2POINTCNT_INVYR, p2eu_INVYR)) %>%
    left_join(stratYr, by = c('STRATUM_CN', 'INVYR')) %>%
    distinct(EVALID, ESTN_UNIT_CN, STRATUM_CN, PLT_CN, .keep_all = TRUE)


  return(pops)

}




## Moving average weights
maWeights <- function(pops, method, lambda){

  ### ---- SIMPLE MOVING AVERAGE
  if (str_to_upper(method) == 'SMA'){
    ## Assuming a uniform weighting scheme
    wgts <- pops %>%
      group_by(YEAR, STATECD) %>%
      summarize(wgt = 1 / length(unique(INVYR)))

    #### ----- Linear MOVING AVERAGE
  } else if (str_to_upper(method) == 'LMA'){

    wgts <- pops %>%
      distinct(YEAR, STATECD, INVYR, .keep_all = TRUE) %>%
      arrange(YEAR, STATECD, INVYR) %>%
      group_by(as.factor(YEAR), as.factor(STATECD)) %>%
      mutate(rank = min_rank(INVYR),
             nsum = sum(1:n())) %>%
      ungroup() %>%
      mutate(wgt = rank / nsum) %>%
      ungroup() %>%
      select(YEAR, STATECD, INVYR, wgt)


    #### ----- EXPONENTIAL MOVING AVERAGE
  } else if (str_to_upper(method) == 'EMA'){
    wgts <- pops %>%
      distinct(YEAR, STATECD, INVYR, .keep_all = TRUE) %>%
      arrange(YEAR, STATECD, INVYR) %>%
      group_by(as.factor(YEAR), as.factor(STATECD)) %>%
      mutate(rank = min_rank(INVYR))


    if (length(lambda) < 2){
      ## Want sum of weighitng functions
      neu <- wgts %>%
        mutate(l = lambda) %>%
        group_by(YEAR, STATECD) %>%
        summarize(l = 1 - first(lambda),
                  sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

      ## Rejoining and computing wgts
      wgts <- wgts %>%
        left_join(neu, by = c('YEAR', 'STATECD')) %>%
        mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
        ungroup() %>%
        select(YEAR, STATECD, INVYR, wgt)
    } else {
      ## Duplicate weights for each level of lambda
      yrWgts <- list()
      for (i in 1:length(unique(lambda))) {
        yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
      }
      wgts <- bind_rows(yrWgts)
      ## Want sum of weighitng functions
      neu <- wgts %>%
        group_by(lambda, YEAR, STATECD) %>%
        summarize(l = 1 - first(lambda),
                  sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))

      ## Rejoining and computing wgts
      wgts <- wgts %>%
        left_join(neu, by = c('lambda', 'YEAR', 'STATECD')) %>%
        mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
        ungroup() %>%
        select(lambda, YEAR, STATECD, INVYR, wgt)
    }
  }

  return(wgts)
}


## Combine most-recent population estimates across states with potentially
## different reporting schedules, i.e., if 2016 is most recent in MI and 2017 is
## most recent in WI, combine them and label as 2017
combineMR <- function(x, grpBy){

  suppressMessages({suppressWarnings({

    maxyears <- x %>%
      select(all_of(grpBy)) %>%
      group_by(.dots = grpBy[!c(grpBy %in% 'YEAR')]) %>%
      summarise(YEAR = max(YEAR, na.rm = TRUE))

    out <- x %>%
      ungroup() %>%
      select(-c(YEAR)) %>%
      left_join(maxyears, by = grpBy[!c(grpBy %in% 'YEAR')])

  })})

  return(out)
}

## Make implicit NA explicity for spatial summaries
prettyNamesSF <- function (tOut, polys, byPlot, grpBy, grpByOrig, tNames, returnSpatial) {

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
      tOut <- left_join(tOut, polys) %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})

    ## Makes it horrible to work with as a dataframe
    if (returnSpatial == FALSE) tOut <- select(tOut, -c(geometry))
  } else if (!is.null(polys) & byPlot){
    polys <- as.data.frame(polys)
    tOut <- left_join(tOut, select(polys, -c(geometry)), by = 'polyID')
  }

  return(tOut)
}

## Choose annual panels to return
filterAnnual <- function(x, grpBy, pltsVar) {
  pltquo <- rlang::enquo(pltsVar)
  x <- x %>%
    mutate(nplts = !!pltquo) %>%
    group_by(INVYR, .dots = grpBy) %>%
    summarize(across(.cols = everything(),  sum, na.rm = TRUE)) %>%
    ## Keep these
    group_by(INVYR) %>%
    mutate(keep = ifelse(INVYR %in% YEAR,
                         ifelse(YEAR == INVYR, 1, 0), ## When TRUE
                         ifelse(nplts == max(nplts, na.rm = TRUE), 1, 0))) %>% ## When INVYR not in YEAR, keep estimates from the inventory where panel has the most plots
    ungroup() %>%
    filter(keep == 1) %>%
    ## If there are multiple reporting years where a panel has the same number of plots
    ## then the estimate will be way too big, we fix this by taking the first row from each output group
    ## If the above worked it will have no effect. If the above failed, it will save our ass.
    mutate(YEAR = INVYR) %>%
    group_by(.dots = grpBy) %>%
    summarize(across(.cols = everything(), first)) %>%
    ungroup()

  return(x)
}




## Estimate skewness in a distribution of values
skewness <- function(x, na.rm = TRUE){

  ## Cut any NA
  if (na.rm) x <- x[!is.na(x)]

  ## Sample size
  n <- length(x)

  ## Estimate the skewness
  skew <- (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)

  return(skew)
}

projectPnts <- function(x, y, slope = NULL, yint = NULL){
  if (is.null(slope)){
    P = data.frame(xOrig = x, yOrig = y)
    P$x <- (P$yOrig+P$xOrig) / 2
    P$y <- P$x
  } else {
    P = data.frame(x, y)
    P$m <- slope
    P$n <- yint
    ## Perp Points
    P$x1 = P$x + -slope
    P$y1 = P$y + 1
    ## Perp Line
    P$m1 = (P$y1-P$y)/(P$x1-P$x)
    P$n1 = P$y - P$m1*P$x
    ## Line intersection
    P$x=(P$n1-P$n)/(P$m-P$m1)
    P$y=P$m*P$x+P$n
  }
  return(P)
}

projectPoints <- function(x, y, slope = 1, yint = 0, returnPoint = TRUE){
  ## Solve for 1:1 line by default

  ## So where does y = mx and y = -1/m * x + b converge
  perp_slope <-  - 1 / slope
  ## Solve for c given x and y
  perp_int <- -perp_slope*x + y

  ## Set equations equal to each other on y
  ## -1/m*x + b = mx
  xproj <- (perp_int - yint) / (slope + -perp_slope)
  yproj <- slope * xproj + yint

  if (returnPoint){
    out <- data.frame(x = xproj, y = yproj)
  } else {
    out <- sqrt((xproj^2) + (yproj^2))
    out <- if_else(xproj < 0, -out, out)
  }
  return(out)
}

matchColClasses <- function(df1, df2) {
  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)

  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)

  for (n in 1:length(sharedColNames)) {
    class(df2[, sharedColNames[n]]) <- sharedColTypes[n]
  }

  return(df2)
}

################ PREVIOUS FUNCTIONS ######################
#### SHANNON'S EVENESS INDEX (H)
#
# speciesObs: vector of observations (species or unique ID)
#
# Returns evenness score (numeric)
####
divIndex <- function(SPCD, TPA, index) {
  # Shannon's Index
  if(index == 'H'){
    species <- unique(SPCD[TPA > 0 & !is.na(TPA)])
    total <- sum(TPA, na.rm = TRUE)

    p <- c() # Empty vector to hold proportions
    for (i in 1:length(species)){
      p[i] <- sum(TPA[SPCD == species[i]], na.rm = TRUE) / total
    }
    value <- -sum(p*log(p), na.rm = TRUE)
  }
  if(index == 'Eh'){
    species <- unique(SPCD[TPA > 0 & !is.na(TPA)])
    total <- sum(TPA, na.rm = TRUE)
    p <- c() # Empty vector to hold proportions
    for (i in 1:length(species)){
      p[i] <- sum(TPA[SPCD == species[i]], na.rm = TRUE) / total
    }
    S <- length(unique(SPCD[TPA > 0 & !is.na(TPA)]))
    if(S == 0) S <- NA
    value <- -sum(p*log(p), na.rm = TRUE) / S
  }
  # Richness
  if(index == 'S'){
    value = length(unique(SPCD[TPA > 0 & !is.na(TPA)])) ## Assumes equal probabilty of detection, not true because of nested sampling design
    }
  # Simpsons Index
  if(index == 'D'){

    # species <- unique(SPCD[TPA > 0])
    # total <- sum(TPA, na.rm = TRUE)
    # p <- c() # Empty vector to hold proportions
    # for (i in 1:length(species)){
    #   p[i] <- sum(TPA[SPCD == species[i]], na.rm = TRUE) / total
    # }
    # value <- 1 /
    #
    #
    # total <- sum(TPA, na.rm = TRUE)
    # props <- data.frame(SPCD, TPA) %>%
    #   filter(TPA > 0) %>%
    #   group_by(SPCD) %>%
    #   summarize(tpa = sum(TPA, na.rm = TRUE)) %>%
    #   mutate(prop = tpa / total) %>%
    #   summarize(D = 1 / sum(prop^2, na.rm = TRUE))
    # props$D[is.infinite(props$D)] <- 0
    # value <- props$D
    # #value <- 1 - (sum(p*(p-1), na.rm = TRUE) / (total * (total-1)))
  }
  # Berger–Parker dominance
  if(index == 'BP'){
    species <- unique(SPCD[TPA > 0])
    total <- sum(TPA, na.rm = TRUE)
    p <- c() # Empty vector to hold proportions
    for (i in 1:length(species)){
      p[i] <- sum(TPA[SPCD == species[i]], na.rm = TRUE) / total
    }
    value <- max(p)
  }
  return(value)
}

areal_par <- function(x, pltSF, polys){
  pltSF <- st_intersection(pltSF, polys[[x]]) %>%
    as.data.frame() %>%
    select(-c('geometry')) # removes artifact of SF object
}

## Exponenetially weighted moving average
ema <- function(x, yrs, var = FALSE){
  l <- 2 / (1 +first(yrs))
  wgts <- c()
  for (i in 1:length(x)) wgts[i] <- l^(i-1)*(1-l)

  if (var){
    #out <- sum(wgts^2 * x,na.rm = TRUE)
    out <- wgts^2 * x
  } else {
    #out <- sum(wgts * x,na.rm = TRUE)
    out <- wgts * x
  }

  return(out)
}


#' @export
makeClasses <- function(x, interval = NULL, lower = NULL, upper = NULL, brks = NULL, numLabs = FALSE){
  if(is.null(brks)){
    ## If min & max isn't specified, then use the data to compute
    low <- ifelse(is.null(lower), min(x, na.rm = TRUE), lower)
    up <- ifelse(is.null(upper), max(x, na.rm = TRUE), upper)
    # Compute class intervals
    brks = c(low)
    while (as.numeric(tail(brks,1)) < as.numeric(up)){
      brks <- c(brks, tail(brks,1) + interval)
    }
  } else {
    low = lower
  }

  # Apply classes to data
  classes <- cut(x, breaks = brks, include.lowest = TRUE, right = FALSE)

  # Convert to numeric (lowest value of interval)
  if (numLabs){
    classes <- as.numeric(classes) * interval -interval + low
  }

  return(classes)
}

#### Basal Area Function (returns sq units of diameter)
basalArea <- function(diameter, DIA_MID = NULL){
  #ba = ((diameter/2)^2) * pi
  # if (!is.null(DIA_MID)){
  #   diameter[is.null(diameter)] <- DIA_MID[is.null(diameter)]
  # }
  # ba = diameter^2 * .005454 # SQ FT, consistency with FIA EVALIDator

  # ba <- case_when(
  #   is.na(diameter) ~ NA_real_,
  #   ## Growth accounting only
  #   diameter < 0 ~ -(diameter^2 * .005454),
  #   TRUE ~ diameter^2 * .005454)

  negative <- data.table::fifelse(diameter < 0, -1, 1)
  ba <- diameter^2 * .005454 * negative


  return(ba)
}

### MODE FUNCTIO

##### Classification of Stand Structural Stage ######
##
## Classifies stand structural stage as pole, mature, late-successional, or mosaic
##  based on relative basal area of live canopy trees within pole, mature & large classes
##
##  diameter: stem DBH (inches) (DIA)
##  crownClass: canopy position of stem, suppressed and open grown excluded (CCLCD)
structHelper <- function(dia, crownClass){

  # Exclude suppressed and open grown stems from analysis
  dia = dia[crownClass %in% c(2,3,4)]

  # Total basal area within plot
  totalBA = sum(basalArea(dia[dia >= 5]), na.rm = TRUE)

  # Calculate proportion of stems in each size class by basal area
  pole = sum(basalArea(dia[dia >= 5 & dia < 10.23622]), na.rm = TRUE) / totalBA
  mature = sum(basalArea(dia[dia >= 10.23622 & dia < 18.11024]), na.rm = TRUE) / totalBA
  large = sum(basalArea(dia[dia >=18.11024]), na.rm = TRUE) / totalBA

  # Series of conditionals to identify stand structural stage based on basal
  #   area proportions in each size class
  if(is.nan(pole) | is.nan(mature) | is.nan(large)){
    stage = 'mosaic'
  } else if ( ((pole + mature) > .67) &  (pole > mature)){
    stage = 'pole'
  } else if(((pole + mature) > .67) &  (pole < mature)){
    stage = 'mature'
  } else if(((mature + large) > .67) & (mature > large)){
    stage = 'mature'
  } else if(((mature + large) > .67) & (mature < large)){
    stage = 'late'
  } else{
    stage = 'mosaic'
  }

  return(as.factor(stage))
}


# Prop basis helper
adjHelper <- function(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR){
  # IF it doesnt exist make it massive
  MACRO_BREAKPOINT_DIA[is.na(MACRO_BREAKPOINT_DIA)] <- 10000
  # Replace DIA with adjustment factors
  adj <- DIA
  adj[is.na(DIA)] <- ADJ_FACTOR_SUBP[is.na(DIA)]
  adj[DIA < 5 & !is.na(DIA)] <- ADJ_FACTOR_MICR[DIA < 5 & !is.na(DIA)]
  adj[DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA & !is.na(DIA)] <- ADJ_FACTOR_SUBP[DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA & !is.na(DIA)]
  adj[DIA >= MACRO_BREAKPOINT_DIA & !is.na(DIA)] <- ADJ_FACTOR_MACR[DIA >= MACRO_BREAKPOINT_DIA & !is.na(DIA)]

  return(adj)

}

# GRM adjustment helper
grmAdj <- function(subtyp, adjMicr, adjSubp, adjMacr) {

  data <- data.frame(typ = as.numeric(subtyp), micr = as.numeric(adjMicr), subp =as.numeric(adjSubp), macr =as.numeric(adjMacr))

  data <- data %>%
    mutate(adj = case_when(
      typ == 0 ~ 0,
      typ == 1 ~ subp,
      typ == 2 ~ micr,
      typ == 3 ~ macr,
    ))

  return(data$adj)
}

# #stratVar <- function(x, a, p2, method, y = NULL){
#   p2 <- first(p2)
#   method <- first(method)
#   a <- first(a)
#   ## Variance Estimation
#   if (is.null(y)){
#     if (method == 'simple'){
#       out <- var(x * a / p2)
#     } else {
#       out <- (sum(x^2) - sum(p2 *  mean(x, na.rm = TRUE)^2)) / (p2 * (p2-1))
#     }
#     ## Covariance Estimation
#   } else {
#     if (method == 'simple'){
#       out <- cov(x,y)
#     } else {
#       out <- (sum(x*y) - sum(p2 * mean(x, na.rm = TRUE) * mean(y, na.rm = TRUE))) / (p2 * (p2-1))
#     }
#   }
#
# }
# stratVar <- function(ESTN_METHOD, x, xStrat, ndif, a, nh, y = NULL, yStrat = NULL){
#   ## Variance
#   if (is.null(y)){
#     v <- ifelse(first(ESTN_METHOD == 'simple'),
#                 var(c(x, numeric(ndif)) * first(a) / nh),
#                 (sum((c(x, numeric(ndif))^2), na.rm = TRUE) - nh * xStrat^2) / (nh * (nh-1)))
#     ## Covariance
#   } else {
#     v <- ifelse(first(ESTN_METHOD == 'simple'),
#                 cov(x,y),
#                 (sum(x*y, na.rm = TRUE) - sum(nh * xStrat *yStrat, na.rm = TRUE)) / (nh * (nh-1)))
#   }
# }

# Helper function to compute variance for estimation units (manages different estimation methods)
unitVarDT <- function(method, ESTN_METHOD, a, nh, w, v, stratMean, stratMean1 = NULL){
  unitM <- unitMean(ESTN_METHOD, a, nh, w, stratMean)
  unitM1 <- unitMean(ESTN_METHOD, a, nh, w, stratMean1)
  if(method == 'var'){
    uv = ifelse(first(ESTN_METHOD) == 'strat',
                ((first(a)^2)/sum(nh)) * (sum(w*nh*v) + sum((1-w)*(nh/sum(nh))*v)),
                ifelse(first(ESTN_METHOD) == 'double',
                       (first(a)^2) * (sum(((nh-1)/(sum(nh)-1))*(nh/sum(nh))*v) + ((1/(sum(nh)-1))*sum((nh/sum(nh))*(stratMean - (unitM/first(a)))^2))),
                       sum(v))) # Stratified random case
  } else { # Compute covariance
    cv = ifelse(first(ESTN_METHOD) == 'strat',
                ((first(a)^2)/sum(nh)) * (sum(w*nh*v) + sum((1-w)*(nh/sum(nh))*v)),
                ifelse(first(ESTN_METHOD) == 'double',
                       (first(a)^2) * (sum(((nh-1)/(sum(nh)-1))*(nh/sum(nh))*v) + ((1/(sum(nh)-1))*sum((nh/sum(nh))*(stratMean - unitM) * (stratMean1 - (unitM1/first(a)))))),
                       sum(v))) # Stratified random case (additive covariance)
  }
}

unitVar <- function(method, ESTN_METHOD, a, nh, w, v, stratMean, unitM, stratMean1 = NULL, unitM1 = NULL){
  if(method == 'var'){
    uv = ifelse(first(ESTN_METHOD) == 'strat',
                ((first(a)^2)/sum(nh)) * (sum(w*nh*v) + sum((1-w)*(nh/sum(nh))*v)),
                ifelse(first(ESTN_METHOD) == 'double',
                       (first(a)^2) * (sum(((nh-1)/(sum(nh)-1))*(nh/sum(nh))*v) + ((1/(sum(nh)-1))*sum((nh/sum(nh))*(stratMean - (unitM/first(a)))^2))),
                       sum(v))) # Stratified random case
  } else {
    cv = ifelse(first(ESTN_METHOD) == 'strat',
                ((first(a)^2)/sum(nh)) * (sum(w*nh*v) + sum((1-w)*(nh/sum(nh))*v)),
                ifelse(first(ESTN_METHOD) == 'double',
                       (first(a)^2) * (sum(((nh-1)/(sum(nh)-1))*(nh/sum(nh))*v) + ((1/(sum(nh)-1))*sum((nh/sum(nh))*(stratMean - unitM) * (stratMean1 - (unitM1/first(a)))))),
                       sum(v))) # Stratified random case (additive covariance)
  }
}

unitVarNew <- function(method, ESTN_METHOD, a, nh, n, w, v, stratMean, unitM, stratMean1 = NULL, unitM1 = NULL){
  if(method == 'var'){
    uv = ifelse(first(ESTN_METHOD) == 'strat',
                ((first(a)^2)/n) * (sum(w*nh*v, na.rm = TRUE) + sum((1-w)*(nh/n)*v, na.rm = TRUE)),
                ifelse(first(ESTN_METHOD) == 'double',
                       (first(a)^2) * (sum(((nh-1)/(n-1))*(nh/n)*v, na.rm = TRUE) + ((1/(n-1))*sum((nh/n)*(stratMean - (unitM/first(a)))^2, na.rm = TRUE))),
                       sum(v, na.rm = TRUE))) # Stratified random case
  } else {
    cv = ifelse(first(ESTN_METHOD) == 'strat',
                ((first(a)^2)/n) * (sum(w*nh*v, na.rm = TRUE) + sum((1-w)*(nh/n)*v, na.rm = TRUE)),
                ifelse(first(ESTN_METHOD) == 'double',
                       (first(a)^2) * (sum(((nh-1)/(n-1))*(nh/n)*v, na.rm = TRUE) + ((1/(n-1))*sum((nh/n)*(stratMean - unitM) * (stratMean1 - (unitM1/first(a))), na.rm = TRUE))),
                       sum(v))) # Stratified random case (additive covariance)
  }
}

## Compute ratio variances at the estimation unit level
rVar <- function(x, y, xVar, yVar, xyCov){
  ## Ratio
  r <- y / x
  ## Ratio variance
  rv <- (1 / x^2) * (yVar + (r^2 * xVar) - (2 * r * xyCov))

  return(rv)
}

# Helper function to compute variance for estimation units (manages different estimation methods)
unitMean <- function(ESTN_METHOD, a, nh, w, stratMean){
  um = ifelse(first(ESTN_METHOD) == 'strat',
              sum(stratMean * w, na.rm = TRUE) * first(a),
              ifelse(first(ESTN_METHOD) == 'double',
                     sum(stratMean * (nh / sum(nh)), na.rm = TRUE) * first(a),
                     mean(stratMean, na.rm = TRUE) * first(a))) # Simple random case
}

## Calculate change for VR
## Calculate change for VR
# vrNAHelper <- function(attribute2, attribute1){
#   ## IF one time is NA, then both must be NA
#   vals <- case_when(
#     is.na(attribute)
#   )
# }


## Replace current attributes with midpoint attributes depending on component
vrAttHelper <- function(attribute, attribute.prev, attribute.mid, attribute.beg, component, remper, oneortwo){

  ## ONLY WORKS FOR ATTRIBUTES DEFINED IN TRE_MIDPNT and TRE_BEGIN
  at <- case_when(
    oneortwo == 2 ~ case_when(
      str_detect(component, c('SURVIVOR|INGROWTH|REVERSION')) ~ attribute / remper,
      str_detect(component, c('CUT|DIVERSION')) ~ attribute.mid / remper),
    oneortwo == 1 ~ case_when(
      str_detect(component, c('SURVIVOR|CUT1|DIVERSION1|MORTALITY1')) ~ case_when(
        !is.na(attribute.beg) ~ - attribute.beg / remper,
        TRUE ~ - attribute.prev / remper)))

  return(at)
}

stratVar <- function(ESTN_METHOD, x, xStrat, ndif, a, nh, y = NULL, yStrat = NULL){
  ## Variance
  if (is.null(y)){
    v <- ifelse(first(ESTN_METHOD == 'simple'),
                var(c(x, numeric(ndif)) * first(a) / nh),
                (sum((c(x, numeric(ndif))^2), na.rm = TRUE) - nh * xStrat^2) / (nh * (nh-1)))
    ## Covariance
  } else {
    v <- ifelse(first(ESTN_METHOD == 'simple'),
                cov(x,y),
                (sum(x*y, na.rm = TRUE) - sum(nh * xStrat * yStrat, na.rm = TRUE)) / (nh * (nh-1)))
  }
}



















## Some base functions for the FIA Database Class
#' @export
summary.FIA.Database <- function(object, ...){
  cat('---- FIA Database Object -----', '\n')
  # Years available
  if (!is.null(object$POP_EVAL$END_INVYR)){
    cat('Reporting Years: ',
        unique(object$POP_EVAL$END_INVYR[order(object$POP_EVAL$END_INVYR)]), '\n')
  }
  # States Covered
  if (!is.null(object$PLOT$STATECD)){
    states <- unique(ifelse(str_length(object$PLOT$STATECD) < 2, paste(0, object$PLOT$STATECD, sep = ''), object$PLOT$STATECD))
    cat('States:          ',
        as.character(unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% states])), '\n')
  }
  # Number of Plots
  if (!is.null(object$POP_STRATUM)){
    eval <- rFIA::findEVALID(object, mostRecent = TRUE, type = 'CURR')
    nPlots <- object$POP_STRATUM %>%
      filter(EVALID %in% eval) %>%
      group_by(ESTN_UNIT_CN, CN) %>%
      summarise(n = first(P2POINTCNT)) %>%
      summarise(n = sum(n))
    cat('Total Plots:     ', sum(nPlots$n), '\n')
  }

  ## Memory Allocated
  mem <- object.size(object)
  cat('Memory Used:     ', format(mem, units = 'MB'), '\n')

  ## Tables included

  cat('Tables:          ', names(object), '\n')

}
#' @export
print.FIA.Database <- function(x, ...){
  cat('---- FIA Database Object -----', '\n')
  # Years available
  if (!is.null(x$POP_EVAL$END_INVYR)){
    cat('Reporting Years: ',
        unique(x$POP_EVAL$END_INVYR[order(x$POP_EVAL$END_INVYR)]), '\n')
  }
  # States Covered
  if (!is.null(x$PLOT$STATECD)){
    states <- unique(ifelse(str_length(x$PLOT$STATECD) < 2, paste(0, x$PLOT$STATECD, sep = ''), x$PLOT$STATECD))
    cat('States:          ',
        as.character(unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% states])), '\n')
  }
  # Number of Plots
  if (!is.null(x$POP_STRATUM)){
    eval <- rFIA::findEVALID(x, mostRecent = TRUE, type = 'CURR')
    nPlots <- x$POP_STRATUM %>%
      filter(EVALID %in% eval) %>%
      group_by(ESTN_UNIT_CN, CN) %>%
      summarise(n = first(P2POINTCNT)) %>%
      summarise(n = sum(n))
    cat('Total Plots:     ', sum(nPlots$n), '\n')
  }

  ## Memory Allocated
  mem <- object.size(x)
  cat('Memory Used:     ', format(mem, units = 'MB'), '\n')

  ## Tables included
  cat('Tables:          ', names(x), '\n', '\n')

  if (length(x) > 1){
    print(sapply(x, as_tibble))
  } else {
    print(as_tibble(x[1]))
  }
  cat('\n')
}

#' @export
summary.Remote.FIA.Database <- function(object, ...){
  cat('---- Remote FIA Database Object -----', '\n')

  # States Covered
  if (!is.null(object$states)){
    cat('States:          ',
        as.character(object$states), '\n')
  }

  ## Memory Allocated
  mem <- object.size(object)
  cat('Memory Used:     ', format(mem, units = 'MB'), '\n')


}
#' @export
print.Remote.FIA.Database <- function(x, ...){
  cat('---- Remote FIA Database Object -----', '\n')

  # States Covered
  if (!is.null(x$states)){
    cat('States:          ',
        as.character(x$states), '\n')
  }

  ## Memory Allocated
  mem <- object.size(x)
  cat('Memory Used:     ', format(mem, units = 'MB'), '\n')

}

#' @export
str.FIA.Database <- function(object, ...) {
  cat(paste('FIA.Database', "\n"))
}

#' @export
str.Remote.FIA.Database <- function(object, ...) {
  cat(paste('Remote.FIA.Database', "\n"))
}








#globalVariables(c('.'))

## Not exported
readFHM <- function(dir, tables = NULL, nCores = 1){
  # Add a slash to end of directory name if missing
  if (str_sub(dir,-1) != '/') dir <- paste(dir, '/', sep = "")
  # Grab all the file names in directory
  files <- list.files(dir)
  inTables <- list()

  # Some warnings
  if(!dir.exists(dir)) {
    stop(paste('Directory', dir, 'does not exist.'))
  }
  if(length(files[str_to_lower(str_sub(files,-4, -1)) == '.csv']) < 1){
    stop(paste('Directory', dir, 'contains no .csv files.'))
  }


  # Only read in the specified tables
  if (!is.null(tables)){
    if (any(str_sub(files, 3, 3) == '_')){
      files <- files[str_sub(files,4,-5) %in% tables]
    } else {
      files <- files[str_sub(files,1,-5) %in% tables]
    }
  }

  # Only csvs
  files <- files[str_to_lower(str_sub(files,-4,-1)) == '.csv']

  inTables <- list()
  for (n in 1:length(files)){
    # Read in and append each file to a list
    file <- fread(paste(dir, files[n], sep = ""), showProgress = FALSE, integer64 = 'double', logical01 = FALSE, nThread = nCores)
    # We don't want data.table formats
    #file <- as.data.frame(file)
    fileName <- str_sub(files[n], 1, -5)

    inTables[[fileName]] <- file
  }

  outTables <- list()
  names(inTables) <- str_sub(names(inTables), 4, -6)

  uniqueNames <- unique(names(inTables))
  ## Works regardless of whether or not there are duplicate names (multiple states)
  for (i in 1:length(uniqueNames)){
    outTables[[uniqueNames[i]]] <- rbindlist(inTables[names(inTables) == uniqueNames[i]], fill=TRUE)
  }

  # NEW CLASS NAME FOR FIA DATABASE OBJECTS
  outTables <- lapply(outTables, as_tibble)
  class(outTables) <- 'FIA.Database'

  ## If you are on windows, close explicitly
  #closeAllConnections()

  return(outTables)
}


## Not exported
getFHM <- function(states,
                   dir = NULL,
                   nCores = 1){

  if (!is.null(dir)){
    # Add a slash to end of directory name if missing
    if (str_sub(dir,-1) != '/'){
      dir <- paste(dir, '/', sep = "")
    }
    # Check to see directory exists, if not, make it
    if(!dir.exists(dir)) {
      dir.create(dir)
      message(paste('Creating directory:', dir))
    }
  }

  ## If dir is not specified, hold in a temporary directory
  if (is.null(dir)){tempDir <- tempdir()}

  #   ## Some warnings up front
  #   ## Do not try to merge ENTIRE with other states
  #   if (length(states) > 1 & any(str_detect(str_to_upper(states), 'ENTIRE'))){
  #     stop('Cannot merge ENITRE with other state tables. ENTIRE includes all state tables combined. Do you only need data for a particular region?')
  #   }
  #   ## Check to make sure states exist
  #   allStates <- c('AL', 'AK', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'HI', 'ID',
  #                  'IL', 'IN', 'IA', 'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 'MI', 'MN', 'MS',
  #                  'MO', 'MT', 'NE', 'NV', 'NH', 'NJ', 'NM', 'NY', 'NC', 'ND', 'OH', 'OK',
  #                  'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WA', 'WV',
  #                  'WI', 'WY', 'AS', 'FM', 'GU', 'MP', 'PW', 'PR', 'VI', 'ENTIRE', 'REF')
  #   if (any(str_to_upper(states) %in% allStates == FALSE)){
  #     missStates <- states[str_to_upper(states) %in% allStates == FALSE]
  #     stop(paste('Data unavailable for: ', paste(as.character(missStates),collapse = ', '), '. Did you use state/terrority abbreviations? e.g. use states = "AL" (correct) instead of states = "ALABAMA".'))
  #   }
  #
  #   ## Check to make sure tables exist
  #   allTables <- c("BOUNDARY", "COND_DWM_CALC", "COND","COUNTY","DWM_COARSE_WOODY_DEBRIS",
  #                  "DWM_DUFF_LITTER_FUEL","DWM_FINE_WOODY_DEBRIS","DWM_MICROPLOT_FUEL",
  #                  "DWM_RESIDUAL_PILE", "DWM_TRANSECT_SEGMENT", "DWM_VISIT","GRND_CVR",
  #                  "INVASIVE_SUBPLOT_SPP","LICHEN_LAB","LICHEN_PLOT_SUMMARY","LICHEN_VISIT",
  #                  "OZONE_BIOSITE_SUMMARY","OZONE_PLOT_SUMMARY","OZONE_PLOT","OZONE_SPECIES_SUMMARY",
  #                  "OZONE_VALIDATION","OZONE_VISIT", "P2VEG_SUBP_STRUCTURE","P2VEG_SUBPLOT_SPP",
  #                  "PLOT_REGEN","PLOT", "PLOTGEOM", "PLOTSNAP","POP_ESTN_UNIT","POP_EVAL_ATTRIBUTE",
  #                  "POP_EVAL_GRP","POP_EVAL_TYP","POP_EVAL","POP_PLOT_STRATUM_ASSGN","POP_STRATUM",
  #                  "SEEDLING_REGEN","SEEDLING","SITETREE","SOILS_EROSION","SOILS_LAB","SOILS_SAMPLE_LOC" ,
  #                  "SOILS_VISIT", "SUBP_COND_CHNG_MTRX","SUBP_COND","SUBPLOT_REGEN","SUBPLOT",
  #                  "SURVEY","TREE_GRM_BEGIN","TREE_GRM_COMPONENT","TREE_GRM_ESTN", "TREE_GRM_MIDPT",
  #                  "TREE_REGIONAL_BIOMASS", "TREE_WOODLAND_STEMS","TREE","VEG_PLOT_SPECIES",
  #                  "VEG_QUADRAT","VEG_SUBPLOT_SPP","VEG_SUBPLOT", "VEG_VISIT",
  #                  'CITATION', 'DIFFERENCE_TEST_PER_ACRE', 'DIFFERENCE_TEST_TOTALS',
  #                  'FIADB_VERSION', 'FOREST_TYPE', 'FOREST_TYPE_GROUP',
  #                  'GRM_TYPE', 'HABTYP_DESCRIPTION', 'HABTYP_PUBLICATION',
  #                  'INVASIVE_SPECIES', 'LICHEN_SPECIES', 'LICHEN_SPP_COMMENTS',
  #                  'NVCS_HEIRARCHY_STRCT', 'NVCS_LEVEL_1_CODES', 'NVCS_LEVEL_2_CODES',
  #                  'NVCS_LEVEL_3_CODES', 'NVCS_LEVEL_4_CODES', 'NVCS_LEVEL_5_CODES',
  #                  'NVCS_LEVEL_6_CODES', 'NVCS_LEVEL_7_CODES', 'NVCS_LEVEL_8_CODES',
  #                  'OWNGRPCD', 'PLANT_DICTIONARY', 'POP_ATTRIBUTE', 'POP_EVAL_TYP_DESCR',
  #                  'RESEARCH_STATION', 'SPECIES', 'SPECIES_GROUP', 'STATE_ELEV', 'UNIT')
  #   if (any(str_to_upper(tables) %in% allTables == FALSE)){
  #     missTables <- tables[str_to_upper(tables) %in% allTables == FALSE]
  #     stop(paste('Tables: ', paste(as.character(missTables),collapse = ', '), ' unavailble. Check the FIA Datamart at https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html for a list of available tables for each state. Alternatively, specify common = TRUE to download the most commonly used tables.
  #
  # Did you accidentally include the state abbreviation in front of the table name? e.g. tables = "AL_PLOT" (wrong) instead of tables = "PLOT" (correct).'))
  #   }
  #


  ## If individual tables are specified, then just grab those .csvs, otherwise download the .zip file, extract and read with fread. Should be quite a bit quicker.
  urlNames <- sapply(states, FUN = function(x){paste0(x,'.zip')})
  urlNames <- c(urlNames)
  ## Set up urls
  urls <- paste0('https://www.fia.fs.fed.us/tools-data/other_data/csv/', urlNames)

  # Make sure state Abbs are in right format
  states <- str_to_upper(states)


  ## If dir is not specified, hold in a temporary directory
  if (is.null(dir)){tempDir <- tempdir()}


  ## Download each state and extract to directory
  for (i in 1:length(states)){
    # Temporary directory to download to
    temp <- tempfile()
    ## Make the URL
    url <- paste0('https://www.fia.fs.fed.us/tools-data/other_data/csv/', states[i],'.zip')
    ## Download as temporary file
    download.file(url, temp)
    ## Extract
    if (is.null(dir)){
      unzip(temp, exdir = tempDir)
    } else {
      unzip(temp, exdir = str_sub(dir, 1, -2))
    }
    unlink(temp)
  }

  ## Read in the files w/ readFHM
  if (is.null(dir)){
    outTables <- readFHM(tempDir, nCores = nCores)

    unlink(tempDir)

    #unlink(tempDir, recursive = TRUE)

  } else {
    outTables <- readFHM(dir, nCores = nCores)
  }

  # NEW CLASS NAME FOR FIA DATABASE OBJECTS
  #outTables <- lapply(outTables, as.data.frame)
  class(outTables) <- 'FHM.Database'

  return(outTables)

}


### Connect to an SQLite3 backend
# connectFIA <- function(dir){
#   ## Connect to the database
#   db <- dbConnect(RSQLite::SQLite(), dir)
#
#   ## Grab the names and store object in list like those held in memory
#   tableNames <- dbListTables(db)
#   outList <- list()
#   for (i in 1:length(tableNames)){
#     outList[[tableNames[i]]] <- tbl(db, tableNames[i])
#   }
#
#   # NEW CLASS NAME FOR FIA DATABASE OBJECTS
#   #outTables <- lapply(outTables, as.data.frame)
#   class(outList) <- 'FIA.Database'
#
#   return(outList)
# }



#### THIS MAY NEED WORK. NOT ALL EVALIDs follow the same coding scheme (ex, CT 2005 --> 95322)
# Look up EVALID codes
#' @export
findEVALID <- function(db = NULL,
                       mostRecent = FALSE,
                       state = NULL,
                       year = NULL,
                       type = NULL){

  #### REWRITING FOR SIMPLICITY #####
  # Joing w/ evaltype code
  ids <- db$POP_EVAL %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_GRP_CN', 'EVAL_TYP')), by = 'EVAL_GRP_CN') %>%
    mutate(place = str_to_upper(LOCATION_NM))

  if (!is.null(state)){
    state <- str_to_upper(state)
    evalGrp <- intData$EVAL_GRP %>%
      select(STATECD, STATE) %>%
      mutate(STATECD = as.numeric(STATECD))
    ## Join state abbs with state codes in popeval
    ids <- left_join(ids, evalGrp, by = 'STATECD')
    # Check if any specified are missing from db
    if (any(unique(state) %in% unique(evalGrp$STATE) == FALSE)){
      missStates <- state[state %in% unique(evalGrp$STATE) == FALSE]
      fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% missStates])
      stop(paste('States: ', toString(fancyName) , 'not found in db.', sep = ''))
    }
    ids <- filter(ids, STATE %in% state)
  }
  if (!is.null(year)){
    #year <- ifelse(str_length(year) == 2, year, str_sub(year, -2,-1))
    ids <- filter(ids, END_INVYR %in% year)
  }
  if (!is.null(type)){
    ids <- filter(ids, EVAL_TYP %in% paste0('EXP', type))
  }
  if (mostRecent) {

    ## Grouped filter wasn't working as intended, use filtering join
    maxYear <- ids %>%
      ## Remove TX, do it seperately
      filter(!(STATECD %in% 48)) %>%
      mutate(place = str_to_upper(LOCATION_NM)) %>%
      group_by(place, EVAL_TYP) %>%
      summarize(END_INVYR = max(END_INVYR, na.rm = TRUE),
                LOCATION_NM = first(LOCATION_NM))

    ## Texas coding standards are very bad
    ## Name two different inventory units with 5 different names
    ## Due to that, only use inventories for the ENTIRE state, sorry
    if (any(ids$STATECD %in% 48)){
      # evalType <- c('EXP_ALL', 'EXP_VOL', '')
      # evalCode <- c('00', '01', '03', '07', '09', '29')
      #
      # txIDS <- ids %>%
      #   filter(STATECD %in% 48) %>%
      #   # ## Removing any inventory that references east or west, sorry
      #   # filter(str_detect(str_to_upper(EVAL_DESCR), 'EAST', negate = TRUE) &
      #   #          str_detect(str_to_upper(EVAL_DESCR), 'WEST', negate = TRUE)) %>%
      #   mutate(typeCode = str_sub(str_trim(EVALID), -2, -1))
      #
      #   mutate(place = str_to_upper(LOCATION_NM)) %>%
      #   group_by(place, EVAL_TYP) %>%
      #   summarize(END_INVYR = max(END_INVYR, na.rm = TRUE),
      #             LOCATION_NM = first(LOCATION_NM))

      ## Will require manual updates, fix your shit texas
      txIDS <- ids %>%
        filter(STATECD %in% 48) %>%
        filter(END_INVYR < 2017) %>%
        filter(END_INVYR > 2006) %>%
        ## Removing any inventory that references east or west, sorry
        filter(str_detect(str_to_upper(EVAL_DESCR), 'EAST', negate = TRUE) &
                 str_detect(str_to_upper(EVAL_DESCR), 'WEST', negate = TRUE)) %>%
        mutate(place = str_to_upper(LOCATION_NM)) %>%
        group_by(place, EVAL_TYP) %>%
        summarize(END_INVYR = max(END_INVYR, na.rm = TRUE),
                  LOCATION_NM = first(LOCATION_NM))

      maxYear <- bind_rows(maxYear, txIDS)
    }

    # ids <- ids %>%
    #   mutate(place = str_to_upper(LOCATION_NM)) %>%
    #   ### TEXAS IS REALLY ANNOYING LIKE THIS
    #   ### FOR NOW, ONLY THE ENTIRE STATE
    #   filter(place %in% c('TEXAS(EAST)', 'TEXAS(WEST)') == FALSE)


    ids <- left_join(maxYear, select(ids, c('place', 'EVAL_TYP', 'END_INVYR', 'EVALID')), by = c('place', 'EVAL_TYP', 'END_INVYR'))
  }

  # Output as vector
  ID <- unique(ids$EVALID)

  return(ID)
}

































