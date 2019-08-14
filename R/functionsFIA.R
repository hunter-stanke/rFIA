
#### A Start up Message ------------
.onAttach <- function(libname, pkgname){
  packageStartupMessage('Download FIA Data Here: https://apps.fs.usda.gov/fia/datamart/datamart.html')
}


# ### Check if columns could be converted to numeric
# numericcharacters <- function(x) {
#   x <- !any(is.na(suppressWarnings(as.numeric(x)))) & is.character(x)
#
#   return(x)
# }


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

    ### REALLY SLOW
    # props <- data.frame(SPCD, TPA) %>%
    #   filter(TPA > 0) %>%
    #   group_by(SPCD) %>%
    #   summarize(tpa = sum(TPA, na.rm = TRUE)) %>%
    #   mutate(prop = tpa / total) %>%
    #   summarize(H = -sum(prop * log(prop), na.rm = TRUE))
    # value <- props$H

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
  if (!is.null(DIA_MID)){
    diameter[is.null(diameter)] <- DIA_MID[is.null(diameter)]
  }
  ba = diameter^2 * .005454 # SQ FT, consistency with FIA EVALIDator
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
  totalBA = sum(basalArea(dia[dia > 5]), na.rm = TRUE)

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
  subtyp[is.na(subtyp)] <- 0

  adj <- subtyp
  # Replace with adj values
  adj[subtyp == 1] <- adjSubp[subtyp == 1]
  adj[subtyp == 2] <- adjMicr[subtyp == 2]
  adj[subtyp == 3] <- adjMacr[subtyp == 3]

  return(adj)
}

# Helper function to compute variance for estimation units (manages different estimation methods)
unitVar <- function(method, ESTN_METHOD, a, nh, w, v, stratMean, unitM, stratMean1 = NULL, unitM1 = NULL){
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

# Helper function to compute variance for estimation units (manages different estimation methods)
unitMean <- function(ESTN_METHOD, a, nh, w, stratMean){
  um = ifelse(first(ESTN_METHOD) == 'strat',
              sum(stratMean * w, na.rm = TRUE) * first(a),
              ifelse(first(ESTN_METHOD) == 'double',
                     sum(stratMean * (nh / sum(nh)), na.rm = TRUE) * first(a),
                     mean(stratMean, na.rm = TRUE) * first(a))) # Simple random case
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
  cat('Tables:          ', names(object))

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

  print(sapply(x, as_tibble))

}

#' @import dplyr
#' @import methods
#' @import sf
#' @import stringr
#' @import gganimate
#' @importFrom data.table fread
#' @importFrom ggplot2 ggplot geom_sf labs scale_fill_viridis_c theme_minimal theme aes element_blank element_text unit ggsave
#' @importFrom parallel makeCluster detectCores mclapply parLapply
#' @importFrom tidyr gather
#' @importFrom pbapply pblapply
#' @importFrom sp over proj4string<- coordinates<- spTransform proj4string
#' @importFrom stats cov var
#' @importFrom utils object.size read.csv tail globalVariables type.convert
NULL

#globalVariables(c('.'))


################### FIA FUNCTIONS #########################
# Read in FIA database files (.csv) from local directory
#' @export
readFIA <- function(dir,
                    common = TRUE,
                    tables = NULL,
                    nCores = 1,
                    ...){

  cat(sys.call()$dir)
  # Add a slash to end of directory name if missing
  if (str_sub(dir,-1) != '/') dir <- paste(dir, '/', sep = "")
  # Grab all the file names in directory
  files <- list.files(dir)
  inTables <- list()

  # Some warnings
  if(!dir.exists(dir)) {
    stop(paste('Directory', dir, 'does not exist.'))
  }
  if(length(files[str_sub(files,-4, -1) == '.csv']) < 1){
    stop(paste('Directory', dir, 'contains no .csv files.'))
  }



  # Only read in the specified tables
  if (!is.null(tables)){
    if (any(str_sub(files, 3, 3) == '_')){
      files <- files[str_sub(files,1,-5) %in% tables]
    } else {
      files <- files[str_sub(files,4,-5) %in% tables]
    }
  }

  # Only csvs
  files <- files[str_sub(files,-4,-1) == '.csv']

  # Only extract the tables needed to run functions in rFIA
  if (common){
    cFiles <- c('COND', 'COND_DWM_CALC', 'INVASIVE_SUBPLOT_SPP', 'PLOT', 'POP_ESTN_UNIT',
                'POP_EVAL', 'POP_EVAL_GRP', 'POP_EVAL_TYP', 'POP_PLOT_STRATUM_ASSGN', 'POP_STRATUM',
                'SUBPLOT', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_ESTN', 'SUBP_COND_CHNG_MTRX')
    if (any(str_sub(files, 3, 3) == '_')){
      files <- files[str_sub(files,4,-5) %in% cFiles]
    } else{
      files <- files[str_sub(files,1,-5) %in% cFiles]
    }
  }


  ## Compute estimates in parallel -- Clusters in windows, forking otherwise
  if (Sys.info()['sysname'] == 'Windows'){
    cl <- makeCluster(nCores) # Set up snow cluster
    inTables <- parLapply(cl, FUN = readFIAHelper1, X = files, dir)
  } else { # Unix systems
    inTables <- mclapply(files, FUN = readFIAHelper1, dir, mc.cores = nCores)
  }
  # Give them some names
  names(inTables) <- files

  # Check for corresponding tables (multiple states)
  # If they exist, loop through and merge corresponding tables
  tableNames <- names(inTables)
  outTables <- list()
  ## Check if the directory has the entire US naming convention or state naming convention
  if (any(str_sub(files, 3, 3) == '_')){ ## ENTIRE NAMING CONVENTION
    if (anyDuplicated(str_sub(tableNames, 4)) != 0){
      for (i in 1:length(unique(str_sub(tableNames, 4)))){
        subList <- inTables[str_sub(tableNames, 4) == unique(str_sub(tableNames,4))[i]]
        name <- unique(str_sub(tableNames, 4, -5))[i]
        # Give a ton of warnings about factors and characters, don't do that
        outTables[[name]] <- suppressWarnings(do.call(rbind, subList))
      }
    } else {
      outTables <- inTables
      names(outTables) <- unique(str_sub(tableNames, 4, -5))
    }
  } else{ ## STATE NAMING CONVENTION
    if (anyDuplicated(tableNames) != 0){
      for (i in 1:length(unique(tableNames))){
        subList <- inTables[tableNames == unique(tableNames)[i]]
        name <- unique(str_sub(tableNames, 1, -5))[i]
        # Give a ton of warnings about factors and characters, don't do that
        outTables[[name]] <- suppressWarnings(do.call(rbind, subList))
      }
    } else {
      outTables <- inTables
      names(outTables) <- unique(str_sub(tableNames, 1, -5))
    }
  }

  # NEW CLASS NAME FOR FIA DATABASE OBJECTS
  outTables <- lapply(outTables, as.data.frame)
  class(outTables) <- 'FIA.Database'

  return(outTables)
}


#### THIS MAY NEED WORK. NOT ALL EVALIDs follow the same coding scheme (ex, CT 2005 --> 95322)
# Look up EVALID codes
#' @export
findEVALID <- function(db = NULL,
                       mostRecent = TRUE,
                       state = NULL,
                       year = NULL,
                       type = NULL){
  # Naming deals
  # if (!is.null(db)){
  #   POP_EVAL = db[['POP_EVAL']]
  # }

  # # R does dumb things when reading strings w/ leading zeros
  # tempID <- c()
  # for (i in 1:length(POP_EVAL$EVALID)){
  #   tempID[i] <- ifelse(str_length(POP_EVAL$EVALID[i]) < 6, paste('0', db$POP_EVAL$EVALID[i], sep = ''), db$POP_EVAL$EVALID[i])
  # }
  #
  # POP_EVAL$EVALID <- tempID


  # Convert types into codes
  if (!is.null(type)){
    if (type == 'ALL'){
      type = '00'
    } else if (type == 'CURR'){
      type = '01'
    } else if(type == 'VOL'){
      type = '02'
    } else if(type == 'GROW'){
      type = '03'
    } else if(type == 'MORT'){
      type == '04'
    } else if(type == 'REMV'){
      type = '05'
    } else if(type == 'CHNG'){
      type = '06'
    } else if(type == 'DWM'){
      type = '07'
    } else if(type == 'REGEN'){
      type = '08'
    }
  }

  # If year is given as whole number, cut it down
  if (!is.null(year)){
    year <- ifelse(str_length(year) == 2, year, str_sub(year, -2,-1))
  }

  # User specifies state & year --> return all evalIDs within state and year
  if(!is.null(state) & !is.null(year) & !is.null(db)){
    base <- paste(first(intData$EVAL_GRP$STATECD[intData$EVAL_GRP$STATE == str_to_upper(state)]), year, type, sep = "")
    ID <- str_subset(db$POP_EVAL$EVALID, base)
  # User doesn't specifiy a year or state, and doesn't want most recent --> all codes present
  } else if (is.null(state) & is.null(year) & mostRecent == FALSE & !is.null(db)){
    if (is.null(type)) {
      ID <- unique(db$POP_EVAL$EVALID)
    } else {
      index <- str_which(str_sub(POP_EVAL$EVALID, -2,-1), type)
      ID <- db$POP_EVAL$EVALID[index]
    }
  # Most recent subset from given database or POP EVAL table
  } else if(is.null(state) & is.null(year) & mostRecent){
    if (is.null(type)) {
      yrMax <- max(db$POP_EVAL$END_INVYR)
      ID <- db$POP_EVAL$EVALID[db$POP_EVAL$END_INVYR == yrMax]
    } else {
      yrMax <- max(db$POP_EVAL$END_INVYR)
      ID <- db$POP_EVAL$EVALID[db$POP_EVAL$END_INVYR == yrMax & str_sub(db$POP_EVAL$EVALID,-2,-1) == type]
    }
  }

  return(ID)
}

## Spatiotemporal clip of FIADB
#' @export
clipFIA <- function(db,
                    mostRecent = TRUE,
                    mask = NULL,
                    matchEval = FALSE,
                    evalID = NULL,
                    designCD = NULL) {

  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Databse". Use readFIA() to load your FIA data.')
  }
  # These states do not allow temporal queries. Things are extremely weird with their eval groups
  if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 41, 53)) & mostRecent){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 41, 53)])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    warning(paste('Most recent subset unavailable for: ', as.character(fancyName) , '. All data required to compute estimates. Returning full database.', sep = ''))
  }
  if (!is.null(mask) & first(class(mask)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('mask must be spatial polygons object of class sp or sf. ')
  }


  if('GRM_COND' %in% names(db) == FALSE & mostRecent & nrow(db$SUBP_COND_CHNG_MTRX) > 1){
    suppressWarnings({
    ## Modify the Subplot Condition Change matrix to include current and previous condition status
    db$SUBP_COND_CHNG_MTRX <- db$SUBP_COND_CHNG_MTRX %>%
      filter(SUBPTYP == 1) %>%
      inner_join(select(db$COND, c('PLT_CN', 'CONDID', 'COND_STATUS_CD')), by = c('PLT_CN', 'CONDID')) %>%
      inner_join(select(db$COND, c('PLT_CN', 'CONDID', 'COND_STATUS_CD')),
                 by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'),
                 suffix = c('.CURR', '.PREV')) %>%
      group_by(PLT_CN, CONDID) %>%
      summarize(COND_CHANGE_CD = ifelse(any(COND_STATUS_CD.CURR == 1 & COND_STATUS_CD.PREV == 1), 1, 0))

    db$GRM_COND <- db$COND %>%
      inner_join(select(db$SUBP_COND_CHNG_MTRX, c('PLT_CN', 'CONDID', 'COND_CHANGE_CD')), by = c('PLT_CN', 'CONDID')) %>%
      distinct(PLT_CN, CONDID, .keep_all = TRUE)%>%
      mutate(CONDPROP_UNADJ = ifelse(COND_CHANGE_CD == 1, CONDPROP_UNADJ, 0)) # Has to be forested currently and at last measurment
    })
    }





  ######### IF USER SPECIES EVALID (OR MOST RECENT), EXTRACT APPROPRIATE PLOTS ##########
  if (!is.null(evalID)){
    # Prevent naming interference
    evalid <- str_sub(evalID, 2,-1)
    # Join appropriate tables and filter out specified EVALIDs
    tempData <- db$PLOT %>%
      mutate(PLT_CN = CN) %>%
      inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID')), by = 'PLT_CN') %>%
      #inner_join(select(POP_ESTN_UNIT, c('EVAL_GRP_CN', 'EVALID')), by = 'EVAL_GRP_CN') %>%
      filter(EVALID %in% evalid)
    # Extract plots which relate to specified EVALID
    db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]
  }

  ## Locate the most recent EVALID and subset plots
  if (mostRecent){
    suppressWarnings({
    tempData <- db$PLOT %>%
      mutate(PLT_CN = CN) %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db$POP_STRATUM, -c('STATECD', 'RSCD')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      left_join(select(db$POP_EVAL, c('EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM', 'LOCATION_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      mutate_if(is.factor,
                as.character) %>%
      filter(!is.na(END_INVYR))
    })
     if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 41, 53))){
      tempMR <- tempData %>%
        filter(STATECD %in% c(69, 72, 78, 41, 53) == FALSE) %>%
        group_by(LOCATION_NM, EVAL_TYP) %>%
        filter(END_INVYR == max(as.numeric(END_INVYR), na.rm = TRUE))
      tempFULL <- tempData %>%
        filter(STATECD %in% c(69, 72, 78, 41, 53))
      tempData <- bind_rows(tempMR, tempFULL)

     } else {
      tempData <- tempData %>%
        group_by(LOCATION_NM, EVAL_TYP) %>%
        filter(END_INVYR == max(as.numeric(END_INVYR), na.rm = TRUE))
     }

    # Both most recent and EVALID is specified, pulls most recent subset of the evalids given
    if (!is.null(evalID)){
      tempData <- tempData %>%
        filter(EVALID %in% evalid)
    }
    # # Locate the most recent EVALIDs present in data
    # yrMax <- max(tempData$END_INVYR, na.rm = TRUE)
    # tempData <- tempData %>%
    #   filter(END_INVYR == yrMax)
    # Extract appropriate PLOTS
    db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]
    #test = db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]

    # Write out evalids sot aht we don't have to repeat above later
    evalid <- unique(tempData$EVALID)

    ## If not most recent, but still want matching evals, go for it.
    } else if (matchEval){
    ### Keeping only years where all states represented are reported for
    if (length(unique(db$POP_EVAL$STATECD)) > 1){
      # Counting number of states measured by year, remove years which don't include all states
      numStates <- db$POP_EVAL %>%
        group_by(END_INVYR, STATECD) %>%
        summarize() %>%
        group_by(END_INVYR) %>%
        summarize(n = n()) %>%
        filter(n == length(unique(db$POP_EVAL$STATECD)))

      db$POP_EVAL <- db$POP_EVAL %>%
        filter(END_INVYR %in% numStates$END_INVYR)
    }

  }




  ################### MISC SUBSETS ON PLOT & OZONE PLOT #############################
  ## If no min-max dates specified, return entire thing
  # if(is.null(minDate)) minDate <- ymd('1900-01-01')
  # if(is.null(maxDate)) maxDate <- now()

  db$PLOT <- db$PLOT %>%
    mutate(PLT_CN = CN) %>%
    #mutate(date = ymd(paste(MEASYEAR, MEASMON, MEASDAY, sep = '-'))) %>%
    mutate(plotID = paste(STATECD, COUNTYCD, PLOT, sep = '')) # Make a unique ID for each plot, irrespective of time

  if (!is.null(designCD)) db$PLOT <- filter(db$PLOT, DESIGNCD == any(as.integer(designCD)))
  #filter(QA_STATUS == 1 | QA_STATUS == 7) %>% # Ensure all quality assurance data are present in FIA database


  ##########################  SPATIAL INTERSECTION ######################################
  ## Oringally did this based solely on PLT_CN, but this method will not work when we try to compute variance
  ##  estimates for a region, because estimates are computed at ESTN_UNIT level. Now the spatial functionality
  ##  of clipFIA is primarily intended to reduce the amount of data which must be held in RAM. As long as all
  ##  plot data is avaliable for any intersecting ESTN_UNIT, we can cut back on memory requirements and still
  ##  produce valid estimates
   if (!is.null(mask)){
     # Convert polygons to an sf object
     mask <- mask %>%
       as('sf')

     ## Make plot data spatial, projected same as polygon layer
     pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
     coordinates(pltSF) <- ~LON+LAT
     proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
     pltSF <- as(pltSF, 'sf') %>%
       st_transform(crs = st_crs(mask)$proj4string)

     # Intersect plot with polygons
     mask$polyID <- 1:nrow(mask)
     suppressMessages({suppressWarnings({
       pltSF <- st_intersection(pltSF, mask) %>%
         as.data.frame() %>%
         select(-c('geometry')) # removes artifact of SF object
     })})

     # Identify the estimation units to which plots within the mask belong to
     estUnits <- pltSF %>%
       inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'STRATUM_CN')), by = 'PLT_CN') %>%
       inner_join(select(db$POP_STRATUM, c('CN', 'ESTN_UNIT_CN')), by = c('STRATUM_CN' = 'CN')) %>%
       group_by(ESTN_UNIT_CN) %>%
       summarize()

     # Identify all the plots which fall inside the above estimation units
     plts <- db$PLOT %>%
       inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'STRATUM_CN')), by = 'PLT_CN') %>%
       inner_join(select(db$POP_STRATUM, c('CN', 'ESTN_UNIT_CN')), by = c('STRATUM_CN' = 'CN')) %>%
       filter(ESTN_UNIT_CN %in% estUnits$ESTN_UNIT_CN) %>%
       group_by(PLT_CN) %>%
       summarize()

     # Clip out the above plots from the full database, will reduce size by a shit pile
     db$PLOT <- filter(db$PLOT, db$PLOT$PLT_CN %in% plts$PLT_CN)

    ## IF ozone is specified, do a seperate intersection (PLOTs not colocated w/ veg PLOTs)
    if (!is.null(db$OZONE_PLOT)) {
      # Seperate spatial object
      ozoneSP <- db$OZONE_PLOT

      ## Make PLOT data spatial, projected same as mask layer
      coordinates(ozoneSP) <- ~LON+LAT
      proj4string(ozoneSP) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
      ozoneSP <- spTransform(ozoneSP, CRSobj = proj4string(mask))

      ## Spatial query (keep only PLOTs that fall within shell)
      db$OZONE_PLOT <- ozoneSP[mask,]

      ## Add coordinates back in to dataframe
      coords <- st_coordinates(db$OZONE_PLOT)
      db$OZONE_PLOT <- db$OZONE_PLOT %>%
        data.frame() %>%
        mutate(LAT = coords[,2]) %>%
        mutate(LON = coords[,1])
    }
   }

  # Starting w/ the most recent subset by state, just want to merge END_INVYR where we have overlapping polygons

  # #### TO deal with polygons which cross state boundaries with different most recent inventory years
  # if(matchEval & mostRecent & !is.null(mask)){
  #   # do a spatial intersection with plot to id poly, then merge the years in pop_eval
  #   evals <- pltSF %>%
  #     inner_join(db$POP_EVAL)
  #   # use pltSF from above
  # }



  #################  APPLY  QUERY TO REMAINING TABLES  ###########################
  ## User gives object names and not list object
  otherTables <- db

  # if (!is.list(otherTables) & !is.null(otherTables)){
  #   clippedData <- list()
  #   clippedData[['PLOT']] <- db$PLOT
  #   if (!is.null(db$OZONE_PLOT)){clippedData[['OZONE_PLOT']] <- db$OZONE_PLOT}
  #   # Spatial query for all tables directly related to PLOT
  #   for (i in 1:length(otherTables)){
  #     # Pull the object with the name listed from the global environment
  #     table <- get(otherTables[i], envir = .GlobalEnv)
  #     if ('PLT_CN' %in% colnames(table) & str_detect(otherTables[i], 'OZONE') == FALSE){
  #       clipTable <- table[table$PLT_CN %in% db$PLOT$CN,]
  #       clippedData[[otherTables[i]]] <- clipTable
  #
  #     } else if (str_detect(otherTables[i], 'OZONE') &'PLT_CN' %in% colnames(table)){
  #       clipTable <- table[table$PLT_CN %in% db$OZONE_PLOT$CN,]
  #       clippedData[[otherTables[i]]] <- clipTable
  #
  #     } else if ('CTY_CN' %in% colnames(table)){
  #       clipTable <- table[table$CTY_CN %in% db$PLOT$CTY_CN,]
  #       clippedData[[otherTables[i]]] <- clipTable
  #
  #     } else if('SRV_CN' %in% colnames(table)){
  #       clipTable <- table[table$SRV_CN %in% db$PLOT$SRV_CN,]
  #       clippedData[[otherTables[i]]] <- clipTable
  #
  #     } else{ # IF it doesn't connect in a way described above, return the whole thing
  #       clippedData[[otherTables[i]]] <- table
  #     }
  #   }
  #
  #   ### IF LIST OBJECT IS INGESTED FROM OTHER TABLES
  # } else

  if (is.list(otherTables)){
    tableNames <- names(otherTables)
    clippedData <- list()

    # Spatial query for all tables directly related to PLOT
    for (i in 1:length(otherTables)){
      # Pull the object with the name listed from the global environment
      table <- otherTables[[i]]
      name <- tableNames[i]

      if (name == "PLOT"){
        clippedData[['PLOT']] <- db$PLOT

      } else if (!is.null(db$OZONE_PLOT) & name == 'OZONE_PLOT'){
          clippedData[['OZONE_PLOT']] <- db$OZONE_PLOT

      } else if ('PLT_CN' %in% colnames(table) & str_detect(name, 'OZONE') == FALSE){
        clipTable <- table[table$PLT_CN %in% db$PLOT$CN,]
        clippedData[[name]] <- clipTable

      } else if ('PLT_CN' %in% colnames(table) & !is.null(db$OZONE_PLOT) & str_detect(name, 'OZONE')){
        clipTable <- table[table$PLT_CN %in% db$OZONE_PLOT$CN,]
        clippedData[[name]] <- clipTable

      } else if ('CTY_CN' %in% colnames(table)){
        clipTable <- table[table$CTY_CN %in% db$PLOT$CTY_CN,]
        clippedData[[name]] <- clipTable

      } else if('SRV_CN' %in% colnames(table)){
        clipTable <- table[table$SRV_CN %in% db$PLOT$SRV_CN,]
        clippedData[[name]] <- clipTable

      } else if(name == 'PLOTGEOM'){
        clipTable <- table[table$CN %in% db$PLOT$CN,]
        clippedData[[name]] <- clipTable

      } else if(name == 'SNAP'){
        clipTable <- table[table$CN %in% db$PLOT$CN,]
        clippedData[[name]] <- clipTable

      } else{ # IF it doesn't connect in a way described above, return the whole thing
        clippedData[[name]] <- table
      }
    }


    # Cascade spatial query for all other tables (indirectly related to PLOT)
    if ("TREE_REGIONAL_BIOMASS" %in% tableNames & 'TREE' %in% names(clippedData)){
      clippedData[["TREE_REGIONAL_BIOMASS"]] <- db$TREE_REGIONAL_BIOMASS %>%
        filter(TRE_CN %in% clippedData$TREE$CN)
    }
    # Deal with the POP Tables
    if(!is.null(evalID) | mostRecent){ # User specified EVALIDs
      # Links only to POP_EVAL
      if ('POP_EVAL' %in% tableNames){
        clippedData[['POP_EVAL']] <- db$POP_EVAL %>%
          filter(EVALID %in% evalid)
        if("POP_EVAL_GRP" %in% tableNames){
          clippedData[["POP_EVAL_GRP"]] <- db$POP_EVAL_GRP %>%
            filter(CN %in% clippedData$POP_EVAL$EVAL_GRP_CN)
        }
        if("POP_EVAL_ATTRIBUTE" %in% tableNames){
          clippedData[["POP_EVAL_ATTRIBUTE"]] <- db$POP_EVAL_ATTRIBUTE %>%
            filter(EVAL_CN %in% clippedData$POP_EVAL$CN)
        }
        if('POP_EVAL_TYP' %in% tableNames){
          clippedData[['POP_EVAL_TYP']] <- db$POP_EVAL_TYP %>%
            filter(EVAL_CN %in% clippedData$POP_EVAL$CN)
        }
      }
      # Links to both POP_EVAL & PLOT
      if ('POP_ESTN_UNIT' %in% tableNames){
        clippedData[['POP_ESTN_UNIT']] <- db$POP_ESTN_UNIT %>%
          filter(EVALID %in% evalid)
      }
      if ('POP_STRATUM' %in% tableNames){
        clippedData[['POP_STRATUM']] <- db$POP_STRATUM %>%
          filter(EVALID %in% evalid)
      }
      if ("POP_PLOT_STRATUM_ASSGN" %in% names(clippedData)){ # Should already be in here becuase contains PLT_CN
        clippedData[["POP_PLOT_STRATUM_ASSGN"]] <- clippedData[["POP_PLOT_STRATUM_ASSGN"]] %>%
          filter(EVALID %in% evalid)
      }

      # User did not specify evalids or most recent
    } else {
      if ("POP_PLOT_STRATUM_ASSGN_STRATUM_ASSGN" %in% names(clippedData) & 'POP_STRATUM' %in% tableNames){
        # Estimation units
        units <- otherTables$POP_STRATUM %>%
          filter(CN %in% clippedData$POP_PLOT_STRATUM_ASSGN$STRATUM_CN) %>%
          select(ESTN_UNIT_CN) %>%
          distinct(ESTN_UNIT_CN)
        # Returns all strata for an estimation unit in which plot falls
        clippedData[['POP_STRATUM']] <- otherTables$POP_STRATUM %>%
          filter(ESTN_UNIT_CN %in% units$ESTN_UNIT_CN)
        if ("POP_ESTN_UNIT" %in% tableNames){
          clippedData[['POP_ESTN_UNIT']] <- db$POP_ESTN_UNIT %>%
            filter(CN %in% units$ESTN_UNIT_CN)
        }
      }
      if ("POP_PLOT_STRATUM_ASSGN" %in% names(clippedData) & 'POP_EVAL' %in% tableNames){
        clippedData[['POP_EVAL']] <- db$POP_EVAL %>%
          filter(EVALID %in% clippedData$POP_PLOT_STRATUM_ASSGN$EVALID)
        if("POP_EVAL_GRP" %in% tableNames){
          clippedData[["POP_EVAL_GRP"]] <- db$POP_EVAL_GRP %>%
            filter(CN %in% clippedData$POP_EVAL$EVAL_GRP_CN)
        }
        if("POP_EVAL_ATTRIBUTE" %in% tableNames){
          clippedData[["POP_EVAL_ATTRIBUTE"]] <- db$POP_EVAL_ATTRIBUTE %>%
            filter(EVAL_CN %in% clippedData$POP_EVAL$CN)
        }
        if('POP_EVAL_TYP' %in% tableNames){
          clippedData[['POP_EVAL_TYP']] <- db$POP_EVAL_TYP %>%
            filter(EVAL_CN %in% clippedData$POP_EVAL$CN)
        }
      }
    }
    # User didn't specify other tables, skip above and just return plot
  } else if (is.null(otherTables) & is.null(db$OZONE_PLOT)){
    clippedData <- db$PLOT
  } else if (is.null(otherTables) & !is.null(db$OZONE_PLOT)){
    clippedData <- list(db$PLOT, db$OZONE_PLOT)
  }

  if (mostRecent) clippedData$mostRecent <- TRUE

  class(clippedData) <- 'FIA.Database'
  return(clippedData)
}

# Stand Structural Stage
#' @export
standStruct <- function(db,
                        grpBy = NULL,
                        polys = NULL,
                        returnSpatial = FALSE,
                        landType = 'forest',
                        areaDomain = NULL,
                        byPlot = FALSE,
                        totals = FALSE,
                        tidy = TRUE,
                        SE = TRUE,
                        progress = TRUE,
                        nCores = 1) {
  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  grpNames <- c(names(db$PLOT), names(db$COND))

  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(first(class(polys))) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (!is.null(grpBy) & class(grpBy) != 'character'){
    stop('grpBy must be of class character. Please specify variable names to group by as quoted list. Example: grpBy = c("OWNGRPCD", "SITECLCD")')
  }
  if (tidy & returnSpatial){
    warning('Returning multiple observations for each areal unit. If returnSpatial = TRUE, tidy = FALSE is recommended.')
  }
  if (landType %in% c('timber', 'forest', 'all') == FALSE){
    stop('landType must be one of: "forest", "timber", or "all".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (any(grpBy %in% grpNames == FALSE)){
    missG <- grpNames[grpNames %in% grpBy == FALSE]
    stop(paste('Columns', paste (as.character(missG), collapse = ', '), 'not found in PLOT or COND tables.'))

  }

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  # ## Pull out individual tables from database object
  # if (!is.null(db)){
  #   TREE <- db[['TREE']]
  #   COND <- db[['COND']]
  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)
  #   POP_PLOT_STRATUM_ASSGN <- db[['POP_PLOT_STRATUM_ASSGN_STRATUM_ASSGN']]
  #   PEU <- db$POP_ESTN_UNIT
  #   POP_EVAL <- db$POP_EVAL
  #   POP_STRATUM <- db$POP_STRATUM
  #   PET <- db$POP_EVAL_TYP
  #   PEG <- db$POP_EVAL_GRP
  # }

  message('Joining FIA Tables.....')


  ### Snag the EVALIDs that are needed & subset POP_EVAL to only include these
  ids <- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPCURR' | EVAL_TYP == 'EXPVOL') %>%
    distinct(EVALID)
  db$POP_EVAL <- db$POP_EVAL %>%
    filter(EVALID %in% ids$EVALID)


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }

    # # Convert back to dataframe
    # db$PLOT <- as.data.frame(db$PLOT) %>%
    #   select(-c('geometry')) # removes artifact of SF object

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
  }

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
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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

  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]

  ## Prep joins and filters
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', 'landD', 'aD_c', grpC)), by = c('PLT_CN')) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
    left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'STATUSCD', 'CCLCD', 'TREECLCD', 'STANDING_DEAD_CD', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE')), by = c('PLT_CN', 'CONDID')) %>%
    mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    rename(YEAR = END_INVYR,
           YEAR_RANGE = REPORT_YEAR_NM) %>%
    mutate_if(is.factor,
              as.character)

  ## Recode a few of the estimation methods to make things easier below
  data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')
  if(!is.null(polys)){
    data <- left_join(data, pltSF, by = 'PLT_CN')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        left_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))

      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }

  }

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$sp



  ####################  COMPUTE ESTIMATES  ###########################
  ### -- BYPLOT -- TPA Estimates at each plot location
  if (byPlot) {
    sOut <- data %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
                stage = structHelper(DIA, CCLCD),
                aAdj = first(aAdj),
                aDI = first(aDI)) #%>%
      # group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      # summarize(pole = sum(CONDPROP_UNADJ[stage == 'pole'] * aDI[stage == 'pole'] * aAdj[stage == 'pole'], na.rm = TRUE),
      #           mature = sum(CONDPROP_UNADJ[stage == 'mature'] * aDI[stage == 'mature'] * aAdj[stage == 'mature'], na.rm = TRUE),
      #           late = sum(CONDPROP_UNADJ[stage == 'late'] * aDI[stage == 'late'] * aAdj[stage == 'late'], na.rm = TRUE),
      #           mosaic = sum(CONDPROP_UNADJ[stage == 'mosaic'] * aDI[stage == 'mosaic'] * aAdj[stage == 'mosaic'], na.rm = TRUE),
      #           faFull = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
      #           plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
  } else {
      # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
      ## produce a unique domain indicator for each individual output observation (ex. Red Oak in Ingham County) to produce valid estimates (otherwise subsampling the
      ## estimation unit, and cause estimates to be inflated substantially
    combos <- select(data, c(grpBy)) %>%
      as.data.frame() %>%
      group_by(.dots = grpBy) %>%
      summarize() %>%
      filter(!is.na(YEAR))
    if(!is.null(polys)){
      combos <- filter(combos, !is.na(polyID))
    }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      message('Computing Summary Statistics.....')

      suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        cl <- makeCluster(nCores)
        if(progress){ # Include progress Bar
          sOut <- pblapply(names(combos), FUN = standStructHelper, combos, data, grpBy, totals, tidy, SE, cl = cl)
        } else { # No progress Bar
          sOut <- parLapply(cl, names(combos), FUN = standStructHelper, combos, data, grpBy, tidy, totals, SE)
        }
      } else { # Unix systems
        if(progress){
          sOut <- pblapply(names(combos), FUN = standStructHelper, combos, data, grpBy, totals, tidy, SE, cl = nCores)
          #sOut <- lapply(names(combos), FUN = standStructHelper, combos, data, grpBy, totals, tidy, SE)
        } else { # No progress Bar, much quicker
          sOut <- mclapply(names(combos), FUN = standStructHelper, combos, data, grpBy, totals, tidy, SE, mc.cores = nCores)
        }
      }
      })

      if (SE){
        # Convert from list to dataframe
        sOut <- do.call(rbind,sOut) %>% #bind_rows(sOut, .id = NULL) %>%
          as.data.frame()

        ## IF the user wants a tidy dataframe at the end, handle it for them
        if (tidy){
          # Gather up all those rando columns
          stage <- gather(sOut, key = 'STAGE', value = 'PERC_AREA', POLE_PERC:MOSAIC_PERC)
          stageSE <- gather(sOut, key = 'STAGE', value = 'PERC_AREA_SE', POLE_PERC_SE:MOSAIC_PERC_SE)
          # Join them back up all nice like
          sTidy <- bind_cols(select(stage, c(names(combos[[1]]), 'STAGE', 'PERC_AREA'), nPlots),
                             select(stageSE, PERC_AREA_SE))
          if(totals){
            stageT <- gather(sOut, key = 'STAGE', value = 'AREA', POLE_AREA:MOSAIC_AREA)
            stageTSE <- gather(sOut, key = 'STAGE', value = 'AREA_SE', POLE_AREA_SE:MOSAIC_AREA_SE)
            # Join them back up all nice like
            sTidy <- bind_cols(sTidy,
                               select(stageT, AREA),
                               select(stageTSE, AREA_SE))
          }
          sOut <- sTidy %>%
            select(-nPlots, nPlots) %>%
            mutate(STAGE = str_split(STAGE, "_", simplify = TRUE)[,1]) %>%
            arrange(YEAR)
        }
      } else {
        # Pull out dataframe
        sOut <- sOut[[1]] %>%
          ungroup() %>%
          as.data.frame()

        if (tidy){
          # Gather up all those rando columns
          stage <- gather(sOut, key = 'STAGE', value = 'PERC_AREA', POLE_PERC:MOSAIC_PERC)
          # Join them back up all nice like
          sTidy <- bind_cols(select(stage, c(names(combos[[1]]), 'STAGE', 'PERC_AREA'), nPlots))
          if(totals){
            stageT <- gather(sOut, key = 'STAGE', value = 'AREA', POLE_AREA:MOSAIC_AREA)
            # Join them back up all nice like
            sTidy <- bind_cols(sTidy,
                               select(stageT, AREA))
          }
          sOut <- sTidy %>%
            select(-nPlots, nPlots) %>%
            mutate(STAGE = str_split(STAGE, "_", simplify = TRUE)[,1]) %>%
            arrange(YEAR)
        }
      }

      # Names for below
      sNames <- names(sOut)[names(sOut) %in% grpBy == FALSE]

      # Return a spatial object
      if ('YEAR' %in% names(sOut)){
        # Return a spatial object
        if (!is.null(polys) & returnSpatial) {
          suppressMessages({suppressWarnings({sOut <- left_join(polys, sOut) %>%
            select(c(grpByOrig, sNames, names(polys))) %>%
            filter(!is.na(polyID))})})
        } else if (!is.null(polys) & returnSpatial == FALSE){
          sOut <- select(sOut, c(grpByOrig, sNames, everything())) %>%
            filter(!is.na(polyID))
        }
      } else { ## Function found no plots within the polygon, so it panics
        combos <- data %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize()
        sOut <- data.frame("YEAR" = combos$YEAR,
                           'POLE_PERC' = rep(NA, nrow(combos)), 'MATURE_PERC' = rep(NA, nrow(combos)),
                           'LATE_PERC' = rep(NA, nrow(combos)),
                           'MOSAIC_PERC' = rep(NA, nrow(combos)), 'POLE_PERC_SE' = rep(NA, nrow(combos)),
                           'MATURE_PERC_SE' = rep(NA, nrow(combos)),
                           'LATE_PERC_SE' = rep(NA, nrow(combos)),
                           'MOSAIC_PERC_SE' = rep(NA, nrow(combos)),
                           "nPlots" = rep(NA, nrow(combos)))
        if (!is.null(polys) & returnSpatial) {
          suppressMessages({suppressWarnings({
            polys = left_join(polys, combos)
            sOut <- left_join(polys, sOut)})})
        } else if (!is.null(polys) & returnSpatial == FALSE){
          sOut <- data.frame(select(sOut, -c('YEAR')), combos)
        }
      }


  } # End byPlot = FALSE
  remove(data)
  remove(db)
  sOut <- filter(sOut, !is.na(YEAR))
  gc()

  return(sOut)
}

##### Species diversity indices
#' @export
diversity <- function(db,
                      grpBy = NULL,
                      polys = NULL,
                      returnSpatial = FALSE,
                      bySizeClass = FALSE,
                      landType = 'forest',
                      treeType = 'live',
                      treeDomain = NULL,
                      areaDomain = NULL,
                      byPlot = FALSE,
                      SE = TRUE,
                      progress = TRUE,
                      nCores = 1) {

  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  grpNames <- c(names(db$PLOT), names(db$TREE), names(db$COND))

  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (!is.null(grpBy) & class(grpBy) != 'character'){
    stop('grpBy must be of class character. Please specify variable names to group by as quoted list. Example: grpBy = c("OWNGRPCD", "SITECLCD")')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (any(grpBy %in% grpNames == FALSE)){
    missG <- grpNames[grpNames %in% grpBy == FALSE]
    stop(paste('Columns', paste (as.character(missG), collapse = ', '), 'not found in PLOT, COND, or TREE tables.'))

  }


  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)


  message('Joining FIA Tables.....')

  ### Snag the EVALIDs that are needed & subset POP_EVAL to only include these
  ids <- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
    distinct(EVALID)
  db$POP_EVAL <- db$POP_EVAL %>%
    filter(EVALID %in% ids$EVALID)


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})

    # # Convert back to dataframe
    # db$PLOT <- as.data.frame(db$PLOT) %>%
    #   select(-c('geometry')) # removes artifact of SF object
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
  }

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
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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

  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ## Prep joins and filters
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
    left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA','SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', grpT, 'typeD', 'tD')), by = c('PLT_CN', 'CONDID')) %>%
    mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    rename(YEAR = END_INVYR,
           YEAR_RANGE = REPORT_YEAR_NM)%>%
    mutate_if(is.factor,
              as.character)

  ## Recode a few of the estimation methods to make things easier below
  data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')
  if(!is.null(polys)){
    data <- left_join(data, pltSF, by = 'PLT_CN')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        left_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))

      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }
  }

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
  data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp

  ## Break into size classes
  if (bySizeClass){
    grpBy <- c(grpBy, 'sizeClass')
    grpByOrig <- c(grpByOrig, 'sizeClass')
    data$sizeClass <- makeClasses(data$DIA, interval = 2)
  }



  ####################  COMPUTE ESTIMATES  ###########################
  ### -- BYPLOT -- TPA Estimates at each plot location
  if (byPlot) {
    dOut <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      summarize(H = divIndex(SPCD, TPA_UNADJ  * tDI, index = 'H'),
                S = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S'),
                Eh = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh'),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0))
      # group_by(.dots = grpBy, ESTN_UNIT_CN, STRATUM_CN, PLT_CN) %>%
      # summarize(H = sum(hCond, na.rm = TRUE),
      #           Eh = sum(EhCond, na.rm = TRUE),
      #           S = sum(sCond* plotIn, na.rm = TRUE),
      #           plotIn = sum(plotIn, na.rm = TRUE)) #%>%
    #filter(S > 0)

    ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
  } else {
    combos <- select(data, c(grpBy)) %>%
      as.data.frame() %>%
      group_by(.dots = grpBy) %>%
      summarize() %>%
      filter(!is.na(YEAR))
    if(!is.null(polys)){
      combos <- filter(combos, !is.na(polyID))
    }
    # List of rows for lapply
    combos <- split(combos, seq(nrow(combos)))

    message('Computing Summary Statistics.....')

    suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      if(progress){ # Include progress Bar
        dOut <- pblapply(names(combos), FUN = diversityHelper, combos, data, grpBy, SE, cl = cl)
      } else { # No progress Bar
        dOut <- parLapply(cl, names(combos), FUN = diversityHelper, combos, data, grpBy, SE)
      }
    } else { # Unix systems
      if(progress){
        dOut <- pblapply(names(combos), FUN = diversityHelper, combos, data, grpBy, SE, cl = nCores)
        #dOut <- pbapply(names(combos), FUN = diversityHelper, combos, data, grpBy, totals, SE)
      } else { # No progress Bar, much quicker
        dOut <- mclapply(names(combos), FUN = diversityHelper, combos, data, grpBy, SE, mc.cores = nCores)
      }
    }
    })

    if (SE){
      # Convert from list to dataframe
      #print(dOut)
      #names(dOut) <- 1:length(dOut)
      #print(names(dOut))
      dOut <- do.call(rbind,dOut)
    } else {
      # Pull out dataframe
      dOut <- dOut[[1]]
    }

    # Snag the names
    dNames <- names(dOut)[names(dOut) %in% grpBy == FALSE]

    # Return a spatial object
    if ('YEAR' %in% names(dOut)){
      # Return a spatial object
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({dOut <- left_join(polys, dOut) %>%
          select(c(grpByOrig, dNames, names(polys))) %>%
          filter(!is.na(polyID))})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        dOut <- select(dOut, c(grpByOrig, dNames, everything())) %>%
          filter(!is.na(polyID))
      }
    } else { ## Function found no plots within the polygon, so it panics
      combos <- data %>%
        as.data.frame() %>%
        group_by(.dots = grpBy) %>%
        summarize()
      dOut <- data.frame("YEAR" = combos$YEAR, "H_a" = rep(NA, nrow(combos)),
                         "H_b" = rep(NA, nrow(combos)), "H_g" = rep(NA, nrow(combos)),
                         "Eh_a" = rep(NA, nrow(combos)),
                         "Eh_b" = rep(NA, nrow(combos)), "Eh_g" = rep(NA, nrow(combos)),
                         "S_a" = rep(NA, nrow(combos)),
                         "S_b" = rep(NA, nrow(combos)), "S_g" = rep(NA, nrow(combos)),
                         "H_a_SE" = rep(NA, nrow(combos)),
                         "Eh_a_SE" = rep(NA, nrow(combos)), "S_a_SE" = rep(NA, nrow(combos)),
                         "nStands" = rep(NA, nrow(combos)))
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({
          polys = left_join(polys, combos)
          dOut <- left_join(polys, dOut)})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        dOut <- data.frame(select(dOut, -c('YEAR')), combos)
      }
    }


  } # End byPlot == FALSE
  dOut <- filter(dOut, !is.na(YEAR))
  return(dOut)
}

## TPA & BAA
#' @export
tpa <- function(db,
                grpBy = NULL,
                polys = NULL,
                returnSpatial = FALSE,
                bySpecies = FALSE,
                bySizeClass = FALSE,
                landType = 'forest',
                treeType = 'live',
                treeDomain = NULL,
                areaDomain = NULL,
                totals = FALSE,
                byPlot = FALSE,
                SE = TRUE,
                progress = TRUE,
                nCores = 1) {

  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  grpNames <- c(names(db$PLOT), names(db$TREE), names(db$COND))
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (!is.null(grpBy) & class(grpBy) != 'character'){
    stop('grpBy must be of class character. Please specify variable names to group by as quoted list. Example: grpBy = c("OWNGRPCD", "SITECLCD")')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (any(grpBy %in% grpNames == FALSE)){
    missG <- grpNames[grpNames %in% grpBy == FALSE]
    stop(paste('Columns', paste (as.character(missG), collapse = ', '), 'not found in PLOT, COND, or TREE tables.'))

  }

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  ## Need a plotCN
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN)


  message('Joining FIA Tables.....')


  ### Snag the EVALIDs that are needed & subset POP_EVAL to only include these
  ids <- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
    distinct(EVALID)
  db$POP_EVAL <- db$POP_EVAL %>%
    filter(EVALID %in% ids$EVALID)


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
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
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ## Prep joins and filters
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
    left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', grpT, 'tD', 'typeD')), by = c('PLT_CN', 'CONDID')) %>%
    mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    rename(YEAR = END_INVYR,
           YEAR_RANGE = REPORT_YEAR_NM) %>%
    mutate_if(is.factor,
              as.character)
  ## Recode a few of the estimation methods to make things easier below
  data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')

  if(!is.null(polys)){
    data <- left_join(data, pltSF, by = 'PLT_CN')
    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }
  }

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
  data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp

  ## Add species to groups
  if (bySpecies) {
    data <- data %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' '))
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }

  ## Break into size classes
  if (bySizeClass){
    grpBy <- c(grpBy, 'sizeClass')
    grpByOrig <- c(grpByOrig, 'sizeClass')
    data$sizeClass <- makeClasses(data$DIA, interval = 2)
  }


  ####################  COMPUTE ESTIMATES  ###########################
  ### -- BYPLOT -- TPA Estimates at each plot location
  if (byPlot) {
    tOut <- data %>%
      filter(EVAL_TYP == 'EXPVOL') %>%
      distinct(PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(TPA = sum(TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                BAA = sum(basalArea(DIA) * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                nStems = length(which(tDI == 1))) %>%
      filter(TPA > 0)

    ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
  } else {
    #if (SE){
    # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
    ## produce a unique domain indicator for each individual output observation (ex. Red Oak in Ingham County) to produce valid estimates (otherwise subsampling the
    ## estimation unit, and cause estimates to be inflated substantially)
    combos <- select(data, c(grpBy)) %>%
      as.data.frame() %>%
      group_by(.dots = grpBy) %>%
      summarize() %>%
      filter(!is.na(YEAR))
    if(!is.null(polys)){
      combos <- filter(combos, !is.na(polyID))
    }
    # List of rows for lapply
    combos <- split(combos, seq(nrow(combos)))

    # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
    if (!is.null(polys)){
      aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
    } else {
      aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
    }

    message('Computing Summary Statistics.....')
    #pb <- progress_bar$new(format = "[:bar] :percent eta: :eta", total = nrow(combos), clear = FALSE, width= 100)

    #for (i in 1:nrow(combos)){ #, .combine = 'rbind', .packages = 'dplyr', .export = c('data')
    #tOut <- foreach(i = 1:nrow(combos), .combine = 'rbind', .packages = 'dplyr', .export = c('data'))

    suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      if(progress){ # Include progress Bar
        tOut <- pblapply(names(combos), FUN = tpaHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = cl)
      } else { # No progress Bar
        tOut <- parLapply(cl, names(combos), FUN = tpaHelper, combos, data, grpBy, aGrpBy, totals, SE)
      }
    } else { # Unix systems
      if(progress){
        tOut <- pblapply(names(combos), FUN = tpaHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = nCores)
      } else { # No progress Bar, much quicker
        tOut <- mclapply(names(combos), FUN = tpaHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
      }
    }
    })

    if (SE){
      # Convert from list to dataframe
      tOut <- do.call(rbind,tOut)
    } else {
      # Pull out dataframe
      tOut <- tOut[[1]]
    }

    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

    # Return a spatial object
    if ('YEAR' %in% names(tOut)){
      # Return a spatial object
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({tOut <- left_join(polys, tOut) %>%
          select(c(grpByOrig, tNames, names(polys))) %>%
          filter(!is.na(polyID))})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        tOut <- select(tOut, c(grpByOrig, tNames, everything())) %>%
          filter(!is.na(polyID))
      }
    } else { ## Function found no plots within the polygon, so it panics
      combos <- data %>%
        as.data.frame() %>%
        group_by(.dots = grpBy) %>%
        summarize()
      tOut <- data.frame("YEAR" = combos$YEAR, "TPA" = rep(NA, nrow(combos)),
                         "BAA" = rep(NA, nrow(combos)), "TPA_PERC" = rep(NA, nrow(combos)),
                         "BAA_PERC" = rep(NA, nrow(combos)), "TPA_SE" = rep(NA, nrow(combos)),
                         "BAA_SE" = rep(NA, nrow(combos)), "TPA_PERC_SE" = rep(NA, nrow(combos)),
                         "BAA_PERC_SE" = rep(NA, nrow(combos)),"nPlots_TREE" = rep(NA, nrow(combos)),
                         "nPlots_AREA" = rep(NA, nrow(combos)))
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({
          polys = left_join(polys, combos)
          tOut <- left_join(polys, tOut)})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        tOut <- data.frame(select(tOut, -c('YEAR')), combos)
      }
    }

  } # End byPlot == FALSE
  remove(data)
  remove(db)
  tOut <- filter(tOut, !is.na(YEAR))
  return(tOut)
}

# Growth and Mortality
#' @export
growMort <- function(db,
                     grpBy = NULL,
                     polys = NULL,
                     returnSpatial = FALSE,
                     bySpecies = FALSE,
                     bySizeClass = FALSE,
                     landType = 'forest',
                     treeType = 'all',
                     treeDomain = NULL,
                     areaDomain = NULL,
                     totals = FALSE,
                     byPlot = FALSE,
                     SE = TRUE,
                     progress = TRUE,
                     nCores = 1) {
  reqTables <- c('PLOT', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_ESTN', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  grpNames <- c(names(db$PLOT), names(db$TREE), names(db$COND))
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (!is.null(grpBy) & class(grpBy) != 'character'){
    stop('grpBy must be of class character. Please specify variable names to group by as quoted list. Example: grpBy = c("OWNGRPCD", "SITECLCD")')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (any(grpBy %in% grpNames == FALSE)){
    missG <- grpNames[grpNames %in% grpBy == FALSE]
    stop(paste('Columns', paste (as.character(missG), collapse = ', '), 'not found in PLOT, COND, or TREE tables.'))

  }

  # These states do not allow change estimates.
  if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    stop(paste('Growth & Mortality Estimates unavailable for: ', as.character(fancyName), sep = ''))
  }

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  message('Joining FIA Tables.....')


  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)

  if ('GRM_COND' %in% names(db) == FALSE & nrow(db$SUBP_COND_CHNG_MTRX) > 1){
    ## Modify the Subplot Condition Change matrix to include current and previous condition status
    db$SUBP_COND_CHNG_MTRX <- db$SUBP_COND_CHNG_MTRX %>%
      filter(SUBPTYP == 1) %>%
      inner_join(select(db$COND, c('PLT_CN', 'CONDID', 'COND_STATUS_CD')), by = c('PLT_CN', 'CONDID')) %>%
      inner_join(select(db$COND, c('PLT_CN', 'CONDID', 'COND_STATUS_CD')),
                 by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'),
                 suffix = c('.CURR', '.PREV')) %>%
      group_by(PLT_CN, CONDID) %>%
      summarize(COND_CHANGE_CD = ifelse(any(COND_STATUS_CD.CURR == 1 & COND_STATUS_CD.PREV == 1), 1, 0))

    db$GRM_COND <- db$COND %>%
      inner_join(select(db$SUBP_COND_CHNG_MTRX, c('PLT_CN', 'CONDID', 'COND_CHANGE_CD')), by = c('PLT_CN', 'CONDID')) %>%
      distinct(PLT_CN, CONDID, .keep_all = TRUE)%>%
      mutate(CONDPROP_UNADJ = ifelse(COND_CHANGE_CD == 1, CONDPROP_UNADJ, 0)) # Has to be forested currently and at last measurment
  }


  ### Snag the EVALIDs that are needed & subset POP_EVAL to only include these
  ids <- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID') %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPGROW' | EVAL_TYP == 'EXPMORT' | EVAL_TYP == 'EXPREMV') %>%
    distinct(EVALID)
  db$POP_EVAL <- db$POP_EVAL %>%
    filter(EVALID %in% ids$EVALID)


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }

    # # Convert back to dataframe
    # db$PLOT <- as.data.frame(db$PLOT) %>%
    #   select(-c('geometry')) # removes artifact of SF object

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
  }

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
  }
  # Tree Type domain indicator
  if (tolower(treeType) == 'all'){
    db$TREE$typeD <- 1
  } else if (tolower(treeType) == 'gs'){
    db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
  }
  # update spatial domain indicator
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ## Prep joins and filters
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', grpP, 'sp', 'aD_p')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'landD', 'aD_c', 'CONDID', grpC)), by = c('PLT_CN')) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
    left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
    left_join(select(db$TREE, c('TRE_CN', 'PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', grpT, 'typeD', 'tD')), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$TREE_GRM_COMPONENT, c('PLT_CN', 'TRE_CN')), by = c('PLT_CN', 'TRE_CN')) %>%
    left_join(select(db$TREE_GRM_ESTN, c('PLT_CN', 'TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('PLT_CN', 'TRE_CN')) %>%
    mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    mutate(tAdj = grmAdj(SUBPTYP_GRM, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    rename(YEAR = END_INVYR,
           YEAR_RANGE = REPORT_YEAR_NM) %>%
    mutate(TPAMORT_UNADJ = ifelse(str_detect(COMPONENT, 'MORT'), TPAMORT_UNADJ, 0)) %>%
    mutate(TPAREMV_UNADJ = ifelse(str_detect(COMPONENT, 'CUT') | str_detect(COMPONENT, 'DIV'), TPAREMV_UNADJ, 0)) %>%
    mutate_if(is.factor,
              as.character)
    #mutate(CONDPROP_UNADJ = ifelse(is.na(TPAMORT_UNADJ), NA, CONDPROP_UNADJ))
  ## Recode a few of the estimation methods to make things easier below
  data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')
  if(!is.null(polys)){
    data <- left_join(data, pltSF, by = 'PLT_CN')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        left_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))

      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }
  }

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp

  ## Add species to groups
  if (bySpecies) {
    data <- data %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' '))
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }

  ## Break into size classes
  if (bySizeClass){
    grpBy <- c(grpBy, 'sizeClass')
    grpByOrig <- c(grpByOrig, 'sizeClass')
    data$sizeClass <- makeClasses(data$DIA, interval = 2)
  }



  ####################  COMPUTE ESTIMATES  ###########################
  ### -- BYPLOT -- TPA Estimates at each plot location
  if (byPlot) {
    tOut <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      #filter(EVALID %in% tID) %>%
      # Compute estimates at plot level
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(TOTAL_TPA = sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE),
                RECR_TPA = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'], na.rm = TRUE),
                MORT_TPA = sum(TPAMORT_UNADJ * tAdj * tDI, na.rm = TRUE),
                REMV_TPA = sum(TPAREMV_UNADJ * tAdj * tDI, na.rm = TRUE),
                LAMBDA = ((TOTAL_TPA / (TOTAL_TPA +MORT_TPA + REMV_TPA - RECR_TPA)) ^ (1/first(REMPER))) - 1)

  } else {
    # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
    ## produce a unique domain indicator for each individual output observation (ex. Red Oak in Ingham County) to produce valid estimates (otherwise subsampling the
    ## estimation unit, and cause estimates to be inflated substantially)
    combos <- select(data, c(grpBy)) %>%
      as.data.frame() %>%
      group_by(.dots = grpBy) %>%
      summarize() %>%
      filter(!is.na(YEAR))
    if(!is.null(polys)){
      combos <- filter(combos, !is.na(polyID))
    }
    # List of rows for lapply
    combos <- split(combos, seq(nrow(combos)))

    # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
    if (!is.null(polys)){
      aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
    } else {
      aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
    }

    message('Computing Summary Statistics.....')

    suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      if(progress){ # Include progress Bar
        tOut <- pblapply(names(combos), FUN = growMortHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = cl)
      } else { # No progress Bar
        tOut <- parLapply(cl, names(combos), FUN = growMortHelper, combos, data, grpBy, aGrpBy, totals, SE)
      }
    } else { # Unix systems
      if(progress){
        tOut <- pblapply(names(combos), FUN = growMortHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = nCores)
      } else { # No progress Bar, much quicker
        tOut <- mclapply(names(combos), FUN = growMortHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
      }
    }
    })

    if (SE){
      # Convert from list to dataframe
      tOut <- do.call(rbind,tOut)
    } else {
      # Pull out dataframe
      tOut <- tOut[[1]]
    }


    # Snag some names for below
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

    # Return a spatial object
    if ('YEAR' %in% names(tOut)){
      # Return a spatial object
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({tOut <- left_join(polys, tOut) %>%
          select(c(grpByOrig, tNames, names(polys))) %>%
          filter(!is.na(polyID))})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        tOut <- select(tOut, c(grpByOrig, tNames, everything())) %>%
          filter(!is.na(polyID))
      }
    } else { ## Function found no plots within the polygon, so it panics
      combos <- data %>%
        as.data.frame() %>%
        group_by(.dots = grpBy) %>%
        summarize()
      tOut <- data.frame("YEAR" = combos$YEAR, "RECR_TPA" = rep(NA, nrow(combos)),
                         "MORT_TPA" = rep(NA, nrow(combos)), "REMV_TPA" = rep(NA, nrow(combos)),
                         "RECR_PERC" = rep(NA, nrow(combos)), "MORT_PERC" = rep(NA, nrow(combos)),
                         "REMV_PERC" = rep(NA, nrow(combos)), "RECR_TPA_SE" = rep(NA, nrow(combos)),
                         "MORT_TPA_SE" = rep(NA, nrow(combos)), "REMV_TPA_SE" = rep(NA, nrow(combos)),
                         "RECR_PERC_SE" = rep(NA, nrow(combos)), "MORT_PERC_SE" = rep(NA, nrow(combos)),
                         "REMV_PERC_SE" = rep(NA, nrow(combos)), "nPlots_TREE" = rep(NA, nrow(combos)),
                         "nPlots_AREA" = rep(NA, nrow(combos)))
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({
          polys = left_join(polys, combos)
          tOut <- left_join(polys, tOut)})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        tOut <- data.frame(select(tOut, -c('YEAR')), combos)
      }
    }


  } # End byPlot == FALSE
  remove(data)
  remove(db)
  tOut <- filter(tOut, !is.na(YEAR))
  return(tOut)
}

# Growth Rates
#' @export
vitalRates <- function(db,
                       grpBy = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       bySpecies = FALSE,
                       bySizeClass = FALSE,
                       landType = 'forest',
                       treeType = 'live',
                       treeDomain = NULL,
                       areaDomain = NULL,
                       totals = FALSE,
                       byPlot = FALSE,
                       SE = TRUE,
                       progress = TRUE,
                       nCores = 1){
  reqTables <- c('PLOT', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_ESTN', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  grpNames <- c(names(db$PLOT), names(db$TREE), names(db$COND))
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (!is.null(grpBy) & class(grpBy) != 'character'){
    stop('grpBy must be of class character. Please specify variable names to group by as quoted list. Example: grpBy = c("OWNGRPCD", "SITECLCD")')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (any(grpBy %in% grpNames == FALSE)){
    missG <- grpNames[grpNames %in% grpBy == FALSE]
    stop(paste('Columns', paste (as.character(missG), collapse = ', '), 'not found in PLOT, COND, or TREE tables.'))

  }

  # These states do not allow change estimates.
  if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    stop(paste('Growth & Mortality Estimates unavailable for: ', as.character(fancyName), sep = ''))
  }

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  message('Joining FIA Tables.....')


  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)

  if ('GRM_COND' %in% names(db) == FALSE & nrow(db$SUBP_COND_CHNG_MTRX) > 1){
    ## Modify the Subplot Condition Change matrix to include current and previous condition status
    db$SUBP_COND_CHNG_MTRX <- db$SUBP_COND_CHNG_MTRX %>%
      filter(SUBPTYP == 1) %>%
      inner_join(select(db$COND, c('PLT_CN', 'CONDID', 'COND_STATUS_CD')), by = c('PLT_CN', 'CONDID')) %>%
      inner_join(select(db$COND, c('PLT_CN', 'CONDID', 'COND_STATUS_CD')),
                 by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'),
                 suffix = c('.CURR', '.PREV')) %>%
      group_by(PLT_CN, CONDID) %>%
      summarize(COND_CHANGE_CD = ifelse(any(COND_STATUS_CD.CURR == 1 & COND_STATUS_CD.PREV == 1), 1, 0))

    db$GRM_COND <- db$COND %>%
      inner_join(select(db$SUBP_COND_CHNG_MTRX, c('PLT_CN', 'CONDID', 'COND_CHANGE_CD')), by = c('PLT_CN', 'CONDID')) %>%
      distinct(PLT_CN, CONDID, .keep_all = TRUE)%>%
      mutate(CONDPROP_UNADJ = ifelse(COND_CHANGE_CD == 1, CONDPROP_UNADJ, 0)) # Has to be forested currently and at last measurment
  }


  ### Snag the EVALIDs that are needed & subset POP_EVAL to only include these
  ids <- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPGROW' | EVAL_TYP == 'EXPCURR') %>%
    distinct(EVALID)
  db$POP_EVAL <- db$POP_EVAL %>%
    filter(EVALID %in% ids$EVALID)


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }

    # # Convert back to dataframe
    # db$PLOT <- as.data.frame(db$PLOT) %>%
    #   select(-c('geometry')) # removes artifact of SF object

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
  }

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
  } else if (tolower(treeType) == 'gs'){
    db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1 & db$TREE$DIA >= 5 & db$TREE$TREECLCD == 2, 1, 0)
  }

  # update spatial domain indicator
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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

  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ## Prep joins and filters
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'REMPER', grpP, 'sp', 'aD_p')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'aD_c', 'landD', 'CONDID', grpC)), by = c('PLT_CN')) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
    left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
    left_join(select(db$TREE, c('TRE_CN', 'PLT_CN', 'CONDID', 'DIA', 'typeD', 'tD', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', grpT)), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$TREE_GRM_COMPONENT, c('PLT_CN', 'TRE_CN', 'ANN_DIA_GROWTH', 'ANN_HT_GROWTH', 'DIA_BEGIN', 'DIA_MIDPT', 'DIA_END',
                            'HT_BEGIN', 'HT_MIDPT', 'HT_END')), by = c('PLT_CN', 'TRE_CN')) %>%
    left_join(select(db$TREE_GRM_ESTN, c('PLT_CN', 'TRE_CN', 'TPAGROW_UNADJ', 'SUBPTYP_GRM', 'ANN_NET_GROWTH')), by = c('PLT_CN', 'TRE_CN')) %>%
    mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    mutate(tAdj = grmAdj(SUBPTYP_GRM, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    rename(YEAR = END_INVYR,
           YEAR_RANGE = REPORT_YEAR_NM) %>%
    mutate(ANN_HT_GROWTH = case_when(is.null(HT_END) ~ (HT_MIDPT - HT_BEGIN) / REMPER / 2,
                                     is.null(HT_BEGIN) ~ (HT_END - HT_MIDPT) / REMPER / 2,
                                     TRUE ~ (HT_END - HT_BEGIN) / REMPER)) %>%
    mutate(ANN_BA_GROWTH = case_when(is.null(DIA_END) ~ (basalArea(DIA_MIDPT) - basalArea(DIA_BEGIN)) / REMPER / 2,
                                     is.null(DIA_BEGIN) ~ (basalArea(DIA_END) - basalArea(DIA_MIDPT)) / REMPER / 2,
                                     TRUE ~ (basalArea(DIA_END) - basalArea(DIA_BEGIN)) / REMPER)) %>%
    mutate_if(is.factor,
              as.character)
  ## Recode a few of the estimation methods to make things easier below
  data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')
  if(!is.null(polys)){
    data <- left_join(data, pltSF, by = 'PLT_CN')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        left_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))

      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }

  }
  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp

  ## Add species to groups
  if (bySpecies) {
    data <- data %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' '))
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }

  ## Break into size classes
  if (bySizeClass){
    grpBy <- c(grpBy, 'sizeClass')
    grpByOrig <- c(grpByOrig, 'sizeClass')
    data$sizeClass <- makeClasses(data$DIA, interval = 2)
  }



  ####################  COMPUTE ESTIMATES  ###########################
  ### -- BYPLOT -- TPA Estimates at each plot location
  if (byPlot) {
    ### Compute total TREES in domain of interest
    tOut <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      #filter(EVALID %in% tID) %>%
      #filter(EVAL_TYP == 'EXPGROW') %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(TPA = sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE),
                DIA_GROW = mean(ANN_DIA_GROWTH * tAdj * tDI, na.rm = TRUE),
                BA_GROW = mean(ANN_BA_GROWTH * tAdj * tDI, na.rm = TRUE),
                BAA_GROW = sum(TPAGROW_UNADJ * ANN_BA_GROWTH * tAdj * tDI, na.rm = TRUE),
                HT_GROW = mean(ANN_HT_GROWTH * tAdj * tDI, na.rm = TRUE),
                NETVOL_GROW = mean(ANN_NET_GROWTH * tAdj * tDI, na.rm = TRUE),
                NETVOL_GROW_ACRE = sum(ANN_NET_GROWTH * TPAGROW_UNADJ *tAdj * tDI, na.rm = TRUE))

    ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
  } else {

    # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
    ## produce a unique domain indicator for each individual output observation (ex. Red Oak in Ingham County) to produce valid estimates (otherwise subsampling the
    ## estimation unit, and cause estimates to be inflated substantially)
    combos <- select(data, c(grpBy)) %>%
      as.data.frame() %>%
      group_by(.dots = grpBy) %>%
      summarize() %>%
      filter(!is.na(YEAR))
    if(!is.null(polys)){
      combos <- filter(combos, !is.na(polyID))
    }
    # List of rows for lapply
    combos <- split(combos, seq(nrow(combos)))

    # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
    if (!is.null(polys)){
      aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
    } else {
      aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
    }

    message('Computing Summary Statistics.....')

    suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      if(progress){ # Include progress Bar
        tOut <- pblapply(names(combos), FUN = vitalRatesHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = cl)
      } else { # No progress Bar
        tOut <- parLapply(cl, names(combos), FUN = vitalRatesHelper, combos, data, grpBy, aGrpBy, totals, SE)
      }
    } else { # Unix systems
      if(progress){
        tOut <- pblapply(names(combos), FUN = vitalRatesHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = nCores)
      } else { # No progress Bar, much quicker
        tOut <- mclapply(names(combos), FUN = vitalRatesHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
      }
    }
    })

    if (SE){
      # Convert from list to dataframe
      tOut <- do.call(rbind,tOut)
    } else {
      # Pull out dataframe
      tOut <- tOut[[1]]
    }


    # Snag some names for below
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]

    # Return a spatial object
    if ('YEAR' %in% names(tOut)){
      # Return a spatial object
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({tOut <- left_join(polys, tOut) %>%
          select(c(grpByOrig, tNames, names(polys))) %>%
          filter(!is.na(polyID))})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        tOut <- select(tOut, c(grpByOrig, tNames, everything())) %>%
          filter(!is.na(polyID))
      }
    } else { ## Function found no plots within the polygon, so it panics
      combos <- data %>%
        as.data.frame() %>%
        group_by(.dots = grpBy) %>%
        summarize()
      tOut <- data.frame("YEAR" = combos$YEAR, "DIA_GROW" = rep(NA, nrow(combos)),
                         "BA_GROW" = rep(NA, nrow(combos)), "BAA_GROW" = rep(NA, nrow(combos)),
                         "HT_GROW" = rep(NA, nrow(combos)), "NETVOL_GROW" = rep(NA, nrow(combos)),
                         "NETVOL_GROW_AC" = rep(NA, nrow(combos)),"DIA_GROW_SE" = rep(NA, nrow(combos)),
                         "BA_GROW_SE" = rep(NA, nrow(combos)), "BAA_GROW_SE" = rep(NA, nrow(combos)),
                         "HT_GROW_SE" = rep(NA, nrow(combos)), "NETVOL_GROW_SE" = rep(NA, nrow(combos)),
                         "NETVOL_GROW_AC_SE" = rep(NA, nrow(combos)),"nPlots_TREE" = rep(NA, nrow(combos)),
                         "nPlots_AREA" = rep(NA, nrow(combos)))
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({
          polys = left_join(polys, combos)
          tOut <- left_join(polys, tOut)})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        tOut <- select(tOut, c(grpByOrig, everything()))
      }
    }


  } # End byPlot == FALSE
  remove(data)
  remove(db)
  tOut <- filter(tOut, !is.na(YEAR))
  return(tOut)
}

# Volume Estimates
#' @export
biomass <- function(db,
                    grpBy = NULL,
                    polys = NULL,
                    returnSpatial = FALSE,
                    bySpecies = FALSE,
                    bySizeClass = FALSE,
                    landType = 'forest',
                    treeType = 'live',
                    treeDomain = NULL,
                    areaDomain = NULL,
                    totals = FALSE,
                    byPlot = FALSE,
                    SE = TRUE,
                    progress = TRUE,
                    nCores = 1) {
  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  grpNames <- c(names(db$PLOT), names(db$TREE), names(db$COND))
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (!is.null(grpBy) & class(grpBy) != 'character'){
    stop('grpBy must be of class character. Please specify variable names to group by as quoted list. Example: grpBy = c("OWNGRPCD", "SITECLCD")')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (any(grpBy %in% grpNames == FALSE)){
    missG <- grpNames[grpNames %in% grpBy == FALSE]
    stop(paste('Columns', paste (as.character(missG), collapse = ', '), 'not found in PLOT, COND, or TREE tables.'))

  }
  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  message('Joining FIA Tables.....')

  # ## Pull out individual tables from database object
  # if (!is.null(db)){
  #   TREE <- db[['TREE']]
  #   COND <- db[['COND']]
  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)
  #   POP_PLOT_STRATUM_ASSGN <- db[['POP_PLOT_STRATUM_ASSGN_STRATUM_ASSGN']]
  #   PEU <- db$POP_ESTN_UNIT
  #   POP_EVAL <- db$POP_EVAL
  #   POP_STRATUM <- db$POP_STRATUM
  #   PET <- db$POP_EVAL_TYP
  #   PEG <- db$POP_EVAL_GRP
  # }


  ### Snag the EVALIDs that are needed & subset POP_EVAL to only include these
  ids <- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
    distinct(EVALID)
  db$POP_EVAL <- db$POP_EVAL %>%
    filter(EVALID %in% ids$EVALID)


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf') %>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }

    # # Convert back to dataframe
    # db$PLOT <- as.data.frame(db$PLOT) %>%
    #   select(-c('geometry')) # removes artifact of SF object

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
  }


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
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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

  ## Prep joins and filters
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ## Prep joins and filters
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', grpP, 'sp', 'aD_p')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'landD', 'aD_c')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
    left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', 'typeD', 'tD',
                                'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG', 'DRYBIO_BG', 'CARBON_AG', 'CARBON_BG', grpT)), by = c('PLT_CN', 'CONDID')) %>%
    mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    rename(YEAR = END_INVYR,
           YEAR_RANGE = REPORT_YEAR_NM) %>%
    mutate_if(is.factor,
              as.character)
  ## Recode a few of the estimation methods to make things easier below
  data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')
  if(!is.null(polys)){
    data <- left_join(data, pltSF, by = 'PLT_CN')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        left_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
      }
  }

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp

  ## Add species to groups
  if (bySpecies) {
    data <- data %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' '))
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }

  ## Break into size classes
  if (bySizeClass){
    grpBy <- c(grpBy, 'sizeClass')
    grpByOrig <- c(grpByOrig, 'sizeClass')
    data$sizeClass <- makeClasses(data$DIA, interval = 2)
  }



  ####################  COMPUTE ESTIMATES  ###########################
  ### -- BYPLOT -- TPA Estimates at each plot location
  if (byPlot) {
    b <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(NETVOL = sum(VOLCFNET * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                SAWVOL = sum(VOLCSNET * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                BIO_AG = sum(DRYBIO_AG * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                BIO_BG = sum(DRYBIO_BG * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                BIO = sum(sum(DRYBIO_AG,DRYBIO_BG,na.rm = TRUE) * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                CARB_AG = sum(CARBON_AG * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                CARB_BG = sum(CARBON_BG * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                CARB = sum(sum(CARBON_AG,CARBON_BG,na.rm=TRUE) * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0))

    ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
  } else {
    combos <- select(data, c(grpBy)) %>%
      as.data.frame() %>%
      group_by(.dots = grpBy) %>%
      summarize() %>%
      filter(!is.na(YEAR))
    if(!is.null(polys)){
      combos <- filter(combos, !is.na(polyID))
    }
    # List of rows for lapply
    combos <- split(combos, seq(nrow(combos)))

    # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
    if (!is.null(polys)){
      aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
    } else {
      aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
    }

    message('Computing Summary Statistics.....')

    suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      if(progress){ # Include progress Bar
        bOut <- pblapply(names(combos), FUN = biomassHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = cl)
      } else { # No progress Bar
        bOut <- parLapply(cl, names(combos), FUN = biomassHelper, combos, data, grpBy, aGrpBy, totals, SE)
      }
    } else { # Unix systems
      if(progress){
        bOut <- pblapply(names(combos), FUN = biomassHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = nCores)
      } else { # No progress Bar, much quicker
        bOut <- mclapply(names(combos), FUN = biomassHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
      }
    }
    })

    if (SE){
      # Convert from list to dataframe
      bOut <- do.call(rbind,bOut)
    } else {
      # Pull out dataframe
      bOut <- bOut[[1]]
    }

    # Snag some names for below
    bNames <- names(bOut)[names(bOut) %in% grpBy == FALSE]

    # Return a spatial object
    if ('YEAR' %in% names(bOut)){
      # Return a spatial object
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({bOut <- left_join(polys, bOut) %>%
          select(c(grpByOrig, bNames, names(polys))) %>%
          filter(!is.na(polyID))})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        bOut <- select(bOut, c(grpByOrig, bNames, everything())) %>%
          filter(!is.na(polyID))
      }
    } else { ## Function found no plots within the polygon, so it panics
      combos <- data %>%
        as.data.frame() %>%
        group_by(.dots = grpBy) %>%
        summarize()
      bOut <- data.frame("YEAR" = combos$YEAR, "NETVOL_ACRE" = rep(NA, nrow(combos)),
                         "SAWVOL_ACRE" = rep(NA, nrow(combos)), "BIO_AG_ACRE" = rep(NA, nrow(combos)),
                         "BIO_BG_ACRE" = rep(NA, nrow(combos)), "BIO_ACRE" = rep(NA, nrow(combos)),
                         "CARB_AG_ACRE" = rep(NA, nrow(combos)),"CARB_BG_ACRE" = rep(NA, nrow(combos)),
                         "CARB_ACRE" = rep(NA, nrow(combos)), "NETVOL_ACRE_SE" = rep(NA, nrow(combos)),
                         "SAWVOL_ACRE_SE"  = rep(NA, nrow(combos)), "BIO_AG_ACRE_SE"  = rep(NA, nrow(combos)),
                         "BIO_BG_ACRE_SE"  = rep(NA, nrow(combos)), "BIO_ACRE_SE"  = rep(NA, nrow(combos)),
                         "CARB_AG_ACRE_SE" = rep(NA, nrow(combos)), "CARB_BG_ACRE_SE" = rep(NA, nrow(combos)),
                         "CARB_ACRE_SE"  = rep(NA, nrow(combos)), "nPlots_VOL" = rep(NA, nrow(combos)),
                         "nPlots_AREA" = rep(NA, nrow(combos)))
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({
          polys = left_join(polys, combos)
          bOut <- left_join(polys, bOut)})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        bOut <- select(bOut, c(grpByOrig, everything()))
      }
    }


  } # End byPlot == FALSE
  remove(data)
  remove(db)
  bOut <- filter(bOut, !is.na(YEAR))
  return(bOut)
}

# Coarse Woody Debris
#' @export
dwm <- function(db,
                grpBy = NULL,
                polys = NULL,
                returnSpatial = FALSE,
                landType = 'forest',
                areaDomain = NULL,
                byPlot = FALSE,
                totals = FALSE,
                tidy = TRUE,
                SE = TRUE,
                progress = TRUE,
                nCores = 1) {
  reqTables <- c('PLOT', 'COND_DWM_CALC', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  grpNames <- c(names(db$PLOT), names(db$COND))
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (!is.null(grpBy) & class(grpBy) != 'character'){
    stop('grpBy must be of class character. Please specify variable names to group by as quoted list. Example: grpBy = c("OWNGRPCD", "SITECLCD")')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (any(grpBy %in% grpNames == FALSE)){
    missG <- grpNames[grpNames %in% grpBy == FALSE]
    stop(paste('Columns', paste (as.character(missG), collapse = ', '), 'not found in PLOT or COND tables.'))

  }

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  message('Joining FIA Tables.....')

  db$COND_DWM_CALC <- db[['COND_DWM_CALC']] %>% mutate(DWM_CN = CN)
  db$COND <- db[['COND']] %>% mutate(CND_CN = CN)
  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)



  ### Snag the EVALIDs that are needed & subset POP_EVAL to only include these
  ids <- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPDWM') %>%
    distinct(EVALID)
  db$POP_EVAL <- db$POP_EVAL %>%
    filter(EVALID %in% ids$EVALID)


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf') %>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }

    # # Convert back to dataframe
    # db$PLOT <- as.data.frame(db$PLOT) %>%
    #   select(-c('geometry')) # removes artifact of SF object

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
  }

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
  }
  # update spatial domain indicator
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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

  ## Prep joins and filters
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]

  ## Prep joins and filters
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', grpP, 'sp', 'aD_p')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CND_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'aD_c', 'landD', 'CONDID', grpC)), by = c('PLT_CN')) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
    left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
    left_join(select(db$COND_DWM_CALC, -c( 'STATECD', 'COUNTYCD', 'UNITCD', 'INVYR', 'MEASYEAR', 'PLOT', 'CONDID', 'EVALID', 'STRATUM_CN')), by = c('PLT_CN', 'CND_CN')) %>%
    mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    #mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    rename(YEAR = END_INVYR,
           YEAR_RANGE = REPORT_YEAR_NM) %>%
    mutate_if(is.factor,
              as.character)
  ## Recode a few of the estimation methods to make things easier below
  data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')
  if(!is.null(polys)){
    data <- left_join(data, pltSF, by = 'PLT_CN')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        left_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))

      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }

  }

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp


  ####################  COMPUTE ESTIMATES  ###########################
  ### -- BYPLOT -- TPA Estimates at each plot location
  if (byPlot) {
    cOut <- data %>%
      #filter(EVAL_TYP == 'EXPDWM') %>%
      distinct(PLT_CN, CONDID, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
       summarize(#DEP_FUEL = mean(FUEL_DEPTH * aDI, na.rm = TRUE), # Mean for depth, because they are not per acre values
      #           DEP_LITTER = mean(LITTER_DEPTH * aDI, na.rm = TRUE),
      #           DEPTH_DUFF = mean(DUFF_DEPTH * aDI, na.rm = TRUE),
                VOL_1HR = sum(FWD_SM_VOLCF_ADJ * aDI, na.rm = TRUE),
                VOL_10HR = sum(FWD_MD_VOLCF_ADJ * aDI, na.rm = TRUE),
                VOL_100HR = sum(FWD_LG_VOLCF_ADJ * aDI, na.rm = TRUE),
                VOL_1000HR = sum(CWD_VOLCF_ADJ * aDI, na.rm = TRUE),
                VOL_PILE = sum(PILE_VOLCF_ADJ * aDI, na.rm = TRUE),
                VOL = sum(VOL_1HR, VOL_10HR, VOL_100HR, VOL_1000HR, VOL_PILE, na.rm = TRUE),
                # BIO_FUEL = sum(FUEL_BIOMASS * aDI / 2000, na.rm = TRUE),
                # BIO_LITTER = sum(LITTER_BIOMASS * aDI / 2000, na.rm = TRUE),
                # BIO_DUFF = sum(DUFF_BIOMASS* aDI / 2000, na.rm = TRUE),
                BIO_1HR = sum(FWD_SM_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                BIO_10HR = sum(FWD_MD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                BIO_100HR = sum(FWD_LG_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                BIO_1000HR = sum(CWD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                BIO_PILE = sum(PILE_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                BIO = sum(BIO_1HR, BIO_10HR, BIO_100HR, BIO_1000HR, BIO_PILE, na.rm = TRUE),
                CARB_1HR = sum(FWD_SM_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                CARB_10HR = sum(FWD_MD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                CARB_100HR = sum(FWD_LG_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                CARB_1000HR = sum(CWD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                CARB_PILE = sum(PILE_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                CARB = sum(CARB_1HR, CARB_10HR, CARB_100HR, CARB_1000HR, CARB_PILE, na.rm = TRUE))

    ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
  } else {
    # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
    ## produce a unique domain indicator for each individual output observation (ex. Red Oak in Ingham County) to produce valid estimates (otherwise subsampling the
    ## estimation unit, and cause estimates to be inflated substantially
    combos <- select(data, c(grpBy)) %>%
      as.data.frame() %>%
      group_by(.dots = grpBy) %>%
      summarize() %>%
      filter(!is.na(YEAR))
    if(!is.null(polys)){
      combos <- filter(combos, !is.na(polyID))
    }
    # List of rows for lapply
    combos <- split(combos, seq(nrow(combos)))

    message('Computing Summary Statistics.....')

    suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      if(progress){ # Include progress Bar
        cOut <- pblapply(names(combos), FUN = dwmHelper, combos, data, grpBy, totals, SE, cl = cl)
      } else { # No progress Bar
        cOut <- parLapply(cl, names(combos), FUN = dwmHelper, combos, data, grpBy, totals, SE)
      }
    } else { # Unix systems
      if(progress){
        cOut <- pblapply(names(combos), FUN = dwmHelper, combos, data, grpBy, totals, SE, cl = nCores)
        #sOut <- lapply(names(combos), FUN = standStructHelper, combos, data, grpBy, totals, tidy, SE)
      } else { # No progress Bar, much quicker
        cOut <- mclapply(names(combos), FUN = dwmHelper, combos, data, grpBy, totals, SE, mc.cores = nCores)
      }
    }
    })

    if (SE){
      # Convert from list to dataframe
      cOut <- do.call(rbind,cOut)  %>%
        as.data.frame()

      ## IF the user wants a tidy dataframe at the end, handle it for them
      if (tidy){
        # Gather up all those rando columns
        vol <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_ACRE', VOL_1HR_ACRE:VOL_PILE_ACRE)
        bio <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_ACRE', BIO_1HR_ACRE:BIO_PILE_ACRE)
        carb <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_ACRE', CARB_1HR_ACRE:CARB_PILE_ACRE)
        volSE <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_ACRE_SE', VOL_1HR_ACRE_SE:VOL_PILE_ACRE_SE)
        bioSE <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_ACRE_SE', BIO_1HR_ACRE_SE:BIO_PILE_ACRE_SE)
        carbSE <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_ACRE_SE', CARB_1HR_ACRE_SE:CARB_PILE_ACRE_SE)
        # Join them back up all nice like
        cTidy <- bind_cols(select(vol, c(names(combos[[1]]), 'FUEL_TYPE', 'VOL_ACRE', 'nPlots')),
                           select(bio, BIO_ACRE),
                           select(carb, CARB_ACRE),
                           select(volSE, VOL_ACRE_SE),
                           select(bioSE, BIO_ACRE_SE),
                           select(carbSE, CARB_ACRE_SE))
        if(totals){
          volT <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_TOTAL', VOL_1HR:VOL_PILE)
          bioT <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_TOTAL', BIO_1HR:BIO_PILE)
          carbT <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_TOTAL', CARB_1HR:CARB_PILE)
          volTSE <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_TOTAL_SE', VOL_1HR_SE:VOL_PILE_SE)
          bioTSE <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_TOTAL_SE', BIO_1HR_SE:BIO_PILE_SE)
          carbTSE <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_TOTAL_SE', CARB_1HR_SE:CARB_PILE_SE)
          cTidy <- bind_cols(cTidy, select(volT, VOL_TOTAL),
                             select(bioT, BIO_TOTAL),
                             select(carbT, CARB_TOTAL),
                             select(volTSE, VOL_TOTAL_SE),
                             select(bioTSE, BIO_TOTAL_SE),
                             select(carbTSE, CARB_TOTAL_SE))
        }
        cOut <- cTidy %>%
          select(-nPlots, nPlots) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE, "_", simplify = TRUE)[,2]) %>%
          arrange(YEAR)
      }
    } else {
      # Pull out dataframe
      cOut <- cOut[[1]]

      if (tidy){
        # Gather up all those rando columns
        vol <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_ACRE', VOL_1HR_ACRE:VOL_PILE_ACRE)
        bio <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_ACRE', BIO_1HR_ACRE:BIO_PILE_ACRE)
        carb <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_ACRE', CARB_1HR_ACRE:CARB_PILE_ACRE)
        # Join them back up all nice like
        cTidy <- bind_cols(select(vol, c(grpBy, 'FUEL_TYPE', 'VOL_ACRE', 'nPlots')),
                           select(bio, BIO_ACRE),
                           select(carb, CARB_ACRE))


        if(totals){
          volT <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_TOTAL', VOL_1HR:VOL_PILE)
          bioT <- gather(cOut, key = 'BIO_TYPE', value = 'BIO_TOTAL', BIO_1HR:BIO_PILE)
          carbT <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_TOTAL', CARB_1HR:CARB_PILE)
          cTidy <- bind_cols(cTidy, select(volT, VOL_TOTAL),
                             select(bioT, BIO_TOTAL),
                             select(carbT, CARB_TOTAL))
        }
        cOut <- cTidy %>%
          select(-nPlots, nPlots) %>%
          mutate(FUEL_TYPE = str_split(FUEL_TYPE, "_", simplify = TRUE)[,2]) %>%
          arrange(YEAR)
      }
    }


    # Names for below
    cNames <- names(cOut)[names(cOut) %in% grpBy == FALSE]

    # Return a spatial object
    if ('YEAR' %in% names(cOut)){
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({cOut <- left_join(polys, cOut) %>%
          select(c(grpByOrig, cNames, names(polys))) %>%
          filter(!is.na(polyID))})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        cOut <- select(cOut, c(grpByOrig, cNames, everything())) %>%
          filter(!is.na(polyID))
      }
    } else { ## Function found no plots within the polygon, so it panics
      combos <- data %>%
        as.data.frame() %>%
        group_by(.dots = grpBy) %>%
        summarize()


      cOut <- data.frame("YEAR" = combos$YEAR,
                         "VOL_1HR_ACRE" = rep(NA, nrow(combos)), "VOL_1HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "VOL_10HR_ACRE" = rep(NA, nrow(combos)), "VOL_10HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "VOL_100HR_ACRE" = rep(NA, nrow(combos)), "VOL_100HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "VOL_1000HR_ACRE" = rep(NA, nrow(combos)), "VOL_1000HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "VOL_PILE_ACRE" = rep(NA, nrow(combos)), "VOL_PILE_ACRE_SE" = rep(NA, nrow(combos)),
                         "VOL_ACRE" = rep(NA, nrow(combos)), "VOL_ACRE_SE" = rep(NA, nrow(combos)),
                         "BIO_1HR_ACRE" = rep(NA, nrow(combos)), "BIO_1HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "BIO_10HR_ACRE" = rep(NA, nrow(combos)), "BIO_10HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "BIO_100HR_ACRE" = rep(NA, nrow(combos)), "BIO_100HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "BIO_1000HR_ACRE" = rep(NA, nrow(combos)), "BIO_1000HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "BIO_PILE_ACRE" = rep(NA, nrow(combos)), "BIO_PILE_ACRE_SE" = rep(NA, nrow(combos)),
                         "BIO_ACRE" = rep(NA, nrow(combos)), "BIO_ACRE_SE" = rep(NA, nrow(combos)),
                         "CARB_1HR_ACRE" = rep(NA, nrow(combos)), "CARB_1HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "CARB_10HR_ACRE" = rep(NA, nrow(combos)), "CARB_10HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "CARB_100HR_ACRE" = rep(NA, nrow(combos)), "CARB_100HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "CARB_1000HR_ACRE" = rep(NA, nrow(combos)), "CARB_1000HR_ACRE_SE" = rep(NA, nrow(combos)),
                         "CARB_PILE_ACRE" = rep(NA, nrow(combos)), "CARB_PILE_ACRE_SE" = rep(NA, nrow(combos)),
                         "CARB_ACRE" = rep(NA, nrow(combos)), "CARB_ACRE_SE" = rep(NA, nrow(combos)),
                         "nPlots" =  rep(NA, nrow(combos)))
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({
          polys = left_join(polys, combos)
          cOut <- left_join(polys, cOut)})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        cOut <- select(cOut, c(grpByOrig, everything()))
      }
    }


  } # End byPlot = FALSE
  remove(data)
  remove(db)
  cOut <- filter(cOut, !is.na(YEAR))
  gc()

  return(cOut)
}

# Invasive coverage
#' @export
invasive <- function(db,
                     grpBy = NULL,
                     polys = NULL,
                     returnSpatial = FALSE,
                     landType = "forest",
                     areaDomain = NULL,
                     byPlot = FALSE,
                     totals = FALSE,
                     SE = TRUE,
                     progress = TRUE,
                     nCores = 1){

  reqTables <- c('PLOT', 'INVASIVE_SUBPLOT_SPP', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  grpNames <- c(names(db$PLOT), names(db$COND))
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (!is.null(grpBy) & class(grpBy) != 'character'){
    stop('grpBy must be of class character. Please specify variable names to group by as quoted list. Example: grpBy = c("OWNGRPCD", "SITECLCD")')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (any(grpBy %in% grpNames == FALSE)){
    missG <- grpNames[grpNames %in% grpBy == FALSE]
    stop(paste('Columns', paste (as.character(missG), collapse = ', '), 'not found in PLOT or COND tables.'))

  }


  ## Link up by species is required, cannot accurately estimate coverage of all species
  grpBy <- c("YEAR", grpBy, 'SYMBOL', 'SCIENTIFIC_NAME', 'COMMON_NAME')
  grpByOrig <- grpBy

  message('Joining FIA Tables.....')

  # if (!is.null(db)) {
  #   INV <- db[["INVASIVE_SUBPLOT_SPP"]]
  #   COND <- db[["COND"]]
  #   #SUBP <- db[['SUBPPLOT']]
  db$PLOT <- db[["PLOT"]] %>% mutate(PLT_CN = CN)
  #   POP_PLOT_STRATUM_ASSGN <- db[["POP_PLOT_STRATUM_ASSGN_STRATUM_ASSGN"]]
  #   PEU <- db$POP_ESTN_UNIT
  #   POP_EVAL <- db$POP_EVAL
  #   POP_STRATUM <- db$POP_STRATUM
  #   PET <- db$POP_EVAL_TYP
  #   PEG <- db$POP_EVAL_GRP
  # }


  ids <- db$POP_EVAL %>% select("CN", "END_INVYR", "EVALID") %>%
    inner_join(select(db$POP_EVAL_TYP, c("EVAL_CN", "EVAL_TYP")), by = c(CN = "EVAL_CN")) %>%
    filter(EVAL_TYP == "EXPVOL" | EVAL_TYP == "EXPCURR") %>%
    distinct(EVALID)
  db$POP_EVAL <- db$POP_EVAL %>% filter(EVALID %in% ids$EVALID)

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys)$proj4string)
    # Intersect plot with polygons
    polys$polyID <- 1:nrow(polys)
    suppressMessages({suppressWarnings({
      pltSF <- st_intersection(pltSF, polys) %>%
        as.data.frame() %>%
        select(-c('geometry')) # removes artifact of SF object
    })})
    # A warning
    if (length(unique(pltSF$PLT_CN)) < 1){
      stop('No plots in db overlap with polys.')
    }

    # # Convert back to dataframe
    # db$PLOT <- as.data.frame(db$PLOT) %>%
    #   select(-c('geometry')) # removes artifact of SF object

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
  }

  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
  }
  # update spatial domain indicator
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(db$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)
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

  ## Prep joins and filters
  ## Which grpByNames are in which table? Helps us subset below
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]

  ## Prep joins and filters
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVASIVE_SAMPLING_STATUS_CD', grpP, 'sp', 'aD_p')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'aD_c', 'landD', 'CONDID', grpC)), by = c('PLT_CN')) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
    left_join(select(db$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
    left_join(select(db$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
    full_join(select(db$INVASIVE_SUBPLOT_SPP, c('PLT_CN', 'COVER_PCT', 'VEG_SPCD', 'CONDID')), by = c("PLT_CN", "CONDID"))
  suppressWarnings({
  data <- data %>%
      left_join(intData$REF_PLANT_DICTIONARY, by = c('VEG_SPCD' = 'SYMBOL')) %>%
      mutate(SYMBOL = VEG_SPCD) %>%
      mutate(aAdj = ifelse(PROP_BASIS == "SUBP", ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
    #mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR, YEAR_RANGE = REPORT_YEAR_NM) %>%
      mutate_if(is.factor,
                as.character)
  })
  data$ESTN_METHOD = recode(.x = data$ESTN_METHOD, `Post-Stratification` = "strat",
                            `Stratified random sampling` = "strat", `Double sampling for stratification` = "double",
                            `Simple random sampling` = "simple", `Subsampling units of unequal size` = "simple")

  if(!is.null(polys)){
    data <- left_join(data, pltSF, by = 'PLT_CN')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        left_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))

      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }
  }

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp


  if (byPlot) {
    invOut <- data %>%
      filter(INVASIVE_SAMPLING_STATUS_CD == 1) %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, VEG_SPCD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(cover = sum(COVER_PCT/100 * CONDPROP_UNADJ * aAdj * aDI * 24^2*pi, na.rm = TRUE)) %>%
      filter(!is.na(SYMBOL))
  } else {
    # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
    ## produce a unique domain indicator for each individual output observation (ex. Red Oak in Ingham County) to produce valid estimates (otherwise subsampling the
    ## estimation unit, and cause estimates to be inflated substantially)
    combos <- select(data, c(grpBy)) %>%
      as.data.frame() %>%
      group_by(.dots = grpBy) %>%
      summarize() %>%
      filter(!is.na(YEAR))
    if(!is.null(polys)){
      combos <- filter(combos, !is.na(polyID))
    }
    # List of rows for lapply
    combos <- split(combos, seq(nrow(combos)))

    # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
    aGrpBy <- grpBy[grpBy %in% c('SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE]

    message('Computing Summary Statistics.....')
    #pb <- progress_bar$new(format = "[:bar] :percent eta: :eta", total = nrow(combos), clear = FALSE, width= 100)

    #for (i in 1:nrow(combos)){ #, .combine = 'rbind', .packages = 'dplyr', .export = c('data')
    #tOut <- foreach(i = 1:nrow(combos), .combine = 'rbind', .packages = 'dplyr', .export = c('data'))

    suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      if(progress){ # Include progress Bar
        invOut <- pblapply(names(combos), FUN = invasiveHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = cl)
      } else { # No progress Bar
        invOut <- parLapply(cl, names(combos), FUN = invasiveHelper, combos, data, grpBy, aGrpBy, totals, SE)
      }
    } else { # Unix systems
      if(progress){
        invOut <- pblapply(names(combos), FUN = invasiveHelper, combos, data, grpBy, aGrpBy, totals, SE, cl = nCores)
      } else { # No progress Bar, much quicker
        invOut <- mclapply(names(combos), FUN = invasiveHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
      }
    }

    if (SE){
      # Convert from list to dataframe
      invOut <- do.call(rbind,invOut)
    } else {
      # Pull out dataframe
      invOut <- invOut[[1]]
    }
    })


    # Snag the names
    invNames <- names(invOut)[names(invOut) %in% grpBy == FALSE]


    # Return a spatial object
    if ('YEAR' %in% names(invOut)){
      # Return a spatial object
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({invOut <- full_join(polys, invOut) %>%
          select(c(grpByOrig, invNames, names(polys))) %>%
          filter(!is.na(polyID))})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        invOut <- select(invOut, c(grpByOrig, invNames, everything())) %>%
          filter(!is.na(polyID))
      }
    } else { ## Function found no plots within the polygon, so it panics
      combos <- data %>%
        as.data.frame() %>%
        group_by(.dots = grpBy) %>%
        summarize()
      invOut <- data.frame("YEAR" = combos$YEAR, "COVER_PCT" = rep(NA, nrow(combos)),
                         "COVER_PCT_SE" = rep(NA, nrow(combos)), 'nPlots_INV' = rep(NA, nrow(combos)),
                         'nPlots_AREA' = rep(NA, nrow(combos)))
      if (!is.null(polys) & returnSpatial) {
        suppressMessages({suppressWarnings({
          polys = left_join(polys, combos)
          invOut <- left_join(polys, invOut)})})
      } else if (!is.null(polys) & returnSpatial == FALSE){
        invOut <- select(invOut, c(grpByOrig, everything()))
      }
    }




  } # End byPlot == FALSE
  remove(data)
  remove(db)
  invOut <- filter(invOut, !is.na(YEAR))
  return(invOut)

}

# # Area Estimates
# area <- function(db,
#                  landType = 'all',
#                  treeType = NULL,
#                  treeDomain = NULL,
#                  areaDomain = NULL,
#                  byLandType = FALSE,
#                  bySpecies = FALSE,
#                  bySizeClass = FALSE,
#                  grpBy = NULL,
#                  polys = NULL,
#                  returnSpatial = FALSE,
#                  byPlot = FALSE,
#                  SE = TRUE) {
# #   '''
# # Total forested area in CT (%) -->
# # Total area occupied by each species (%), on forestland --> forested area total
# # Total forested area by county, areal unit --> use land area total
# # Within a state, and county, total area of forestland --> need land area total
# # Total area occupied by a phyiographic class on forestland --> use forested area total
# # What if they want to group by state and county? which to use as denominator?
# #
# # ALWAYS DEFAULT PROPORTION ESTIMATES TO THE LOWEST GROUPING VALUE, ONLY WAY TO GO ABOUT IT
# #
# # not sure. Leaving area alone for now, play again when all other functions required by NPS area complete
# #
# # '''
#
#   # if(!is.null(eval(areaDomain))){
#   #   print('test')
#   # }
#
#   # Save original grpByfor pretty return with spatial objects
#   grpBy <- c('YEAR', grpBy)
#   grpByOrig <- grpBy
#
#   ## Pull out individual tables from database object
#   if (!is.null(db)){
#     TREE <- db[['TREE']]
#     COND <- db[['COND']]
#     PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)
#     POP_PLOT_STRATUM_ASSGN <- db[['POP_PLOT_STRATUM_ASSGN_STRATUM_ASSGN']]
#     PEU <- db$POP_ESTN_UNIT
#     POP_EVAL <- db$POP_EVAL
#     POP_STRATUM <- db$POP_STRATUM
#     PET <- db$POP_EVAL_TYP
#     PEG <- db$POP_EVAL_GRP
#   }
#
#
#
#   ### Keeping only years where all states represented are reported for
#   if (length(unique(POP_EVAL$STATECD)) > 1){
#     # Counting number of states measured by year, remove years which don't include all states
#     numStates <- POP_EVAL %>%
#       group_by(END_INVYR, STATECD) %>%
#       summarize() %>%
#       group_by(END_INVYR) %>%
#       summarize(n = n()) %>%
#       filter(n == length(unique(POP_EVAL$STATECD)))
#
#     POP_EVAL <- POP_EVAL %>%
#       filter(END_INVYR %in% numStates$END_INVYR)
#   }
#
#   ### Snag the EVALIDs that are needed & subset POP_EVAL to only include these
#   ids <- POP_EVAL %>%
#     select('CN', 'END_INVYR', 'EVALID') %>%
#     inner_join(select(PET, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
#     filter(EVAL_TYP == 'EXPCURR') %>%
#     distinct(EVALID)
#   POP_EVAL <- POP_EVAL %>%
#     filter(EVALID %in% ids$EVALID)
#
#
#   ### AREAL SUMMARY PREP
#   if(!is.null(polys)) {
#     # Convert polygons to an sf object
#     polys <- polys %>%
#       as('sf') %>%
#       mutate_if(is.factor,
#                 fct_explicit_na,
#                 na_level = "NA")
#     # Add shapefile names to grpBy
#     grpBy = c(names(polys)[1:ncol(polys)-1], 'polyID', grpBy)
#     ## Make plot data spatial, projected same as polygon layer
#     coordinates(PLOT) <- ~LON+LAT
#     proj4string(PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#     PLOT <- as(PLOT, 'sf') %>%
#       st_transform(plotSF, crs = st_crs(polys)$proj4string)
#     # Intersect plot with polygons
#     polys$polyID <- 1:nrow(polys)
#     suppressMessages({
#       PLOT <- st_intersection(PLOT, polys)
#     })
#
#     # Convert back to dataframe
#     PLOT <- as.data.frame(PLOT) %>%
#       select(-c('geometry')) # removes artifact of SF object
#
#   } else if (byPlot & returnSpatial){
#     ## Make plot data spatial, projected same as polygon layer
#     coordinates(PLOT) <- ~LON+LAT
#     proj4string(PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#     PLOT <- as(PLOT, 'sf')
#   }
#
#
#   ## Prep joins and filters
#   data <- PLOT %>%
#     inner_join(select(COND, -c('INVYR', 'STATECD', 'UNITCD', 'COUNTYCD', 'PLOT')), by = c('PLT_CN')) %>%
#     inner_join(select(POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
#     inner_join(select(POP_STRATUM, -c('STATECD', 'RSCD')), by = c('STRATUM_CN' = 'CN')) %>%
#     inner_join(select(PEU, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
#     inner_join(select(POP_EVAL, c('EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
#     inner_join(select(PET, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
#     inner_join(select(PEG, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
#     left_join(select(TREE, -c('INVYR', 'STATECD', 'UNITCD', 'COUNTYCD', 'PLOT')), by = c('PLT_CN', 'CONDID')) %>%
#     mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
#     distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, .keep_all = TRUE) %>%
#     mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
#     rename(YEAR = END_INVYR,
#            YEAR_RANGE = REPORT_YEAR_NM)
#   # # Make s the data huge, so only join up if the user wants to specify a criteria in the tree database
#   # if (!is.null(eval(treeDomain, data)) | bySpecies | bySizeClass){
#   #   data <- data %>%
#   #
#   # } else {
#   #   data <- data %>%
#   #     distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE)
#   # }
#
#   ## Recode a few of the estimation methods to make things easier below
#   data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
#                             `Post-Stratification` = 'strat',
#                             `Stratified random sampling` = 'strat',
#                             `Double sampling for stratification` = 'double',
#                             `Simple random sampling` = 'simple',
#                             `Subsampling units of unequal size` = 'simple')
#
#   if(byLandType == FALSE){
#     ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
#     # Land type domain indicator
#     if (tolower(landType) == 'forest'){
#       landD <- ifelse(data$COND_STATUS_CD == 1, 1, 0)
#     } else if (tolower(landType) == 'timber'){
#       landD <- ifelse(data$COND_STATUS_CD == 1 & data$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & data$RESERVCD == 0, 1, 0)
#     } else if (tolower(landType) == 'non-forest'){
#       landD <- ifelse(data$COND_STATUS_CD == 2)
#     } else if (tolower(landType) == 'water'){
#       landD <- ifelse(data$COND_STATUS_CD == 3 | data$COND_STATUS_CD == 4)
#     } else if (tolower(landType) == 'census water'){
#       landD <- ifelse(data$COND_STATUS_CD == 4)
#     } else if (tolower(landType) == 'non-census water'){
#       landD <- ifelse(data$COND_STATUS_CD == 4)
#     } else if (tolower(landType) == 'all') {
#       landD <- 1
#     }
#   } else {
#     landD <- 1
#   }
#   # Tree Type domain indicator
#   if(is.null(treeType)) {
#     typeD <- 1
#   } else if (tolower(treeType) == 'live'){
#     typeD <- ifelse(data$STATUSCD == 1, 1, 0)
#   } else if (tolower(treeType) == 'dead'){
#     typeD <- ifelse(data$STATUSCD == 2 & data$STANDING_DEAD_CD == 1, 1, 0)
#   } else if (tolower(treeType) == 'gs'){
#     typeD <- ifelse(data$STATUSCD == 1 & data$DIA >= 5 & data$TREECLCD == 2, 1, 0)
#   } else if (tolower(treeType) == 'all'){
#     typeD <- 1
#   }
#   # User defined domain indicator for area (ex. specific forest type)
#   areaDomain <- substitute(areaDomain)
#   aD <- eval(areaDomain, data) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
#   if(!is.null(aD)) aD[is.na(aD)] <- 0 # Make NAs 0s. Causes bugs otherwise
#   if(is.null(aD)) aD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
#   # Same as above for tree (ex. trees > 20 ft tall)
#   treeDomain <- substitute(treeDomain)
#   tD <- eval(treeDomain, data) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
#   if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
#   if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
#   ## Comprehensive indicator function
#   data$aDI <- landD * aD
#   data$tDI <- landD * aD * tD * typeD
#
#
#   ## Add species to groups
#   if (bySpecies) {
#     data <- data %>%
#       left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
#       mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' '))
#     # data$SPCD <- fct_explicit_na(data, SPCD)
#     # data$COMMON_NAME <- fct_explicit_na(data, COMMON_NAME)
#     # data$GENUS <- fct_explicit_na(data, GENUS)
#     # data$SPECIES <- fct_explicit_na(data, SPECIES)
#
#     grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
#     grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
#   }
#
#   ## Break into size classes
#   if (bySizeClass){
#     grpBy <- c(grpBy, 'sizeClass')
#     grpByOrig <- c(grpByOrig, 'sizeClass')
#     data$sizeClass <- makeClasses(data$DIA, interval = 2) %>%
#       fct_explicit_na()
#   }
#
#   # Make a new column that describes the land type
#   if (byLandType){
#     grpBy <- c(grpBy, 'landType')
#     grpByOrig <- c(grpByOrig, 'landType')
#     data <- data %>%
#       mutate(landType = case_when(
#         COND_STATUS_CD == 1 & SITECLCD %in% c(1:6) & RESERVCD ==0 ~ 'Timber',
#         COND_STATUS_CD == 1 ~ 'Non-Timber Forest',
#         COND_STATUS_CD == 2 ~ 'Non-Forest',
#         COND_STATUS_CD == 3 | COND_STATUS_CD == 4 ~ 'Water'))
#     #data$landType <- fct_explicit_na(data$landType)
#   }
#
#
#
#   ####################  COMPUTE ESTIMATES  ###########################
#   ### -- BYPLOT -- TPA Estimates at each plot location
#   if (byPlot) {
#     aOut <- data %>%
#       group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
#       summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
#                 tDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0),
#                 a = first(AREA_USED),
#                 p1EU = first(P1PNTCNT_EU),
#                 p1 = first(P1POINTCNT),
#                 p2 = first(P2POINTCNT),
#                 aAdj = first(aAdj),
#                 aDI = first(aDI)) %>%
#       group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#       summarize(propArea = sum(CONDPROP_UNADJ * tDI * aAdj, na.rm = TRUE),
#                 total = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
#                 plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))
#
#     ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
#   } else {
#     if (SE){
#       # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
#       ## produce a unique domain indicator for each individual output observation (ex. Red Oak in Ingham County) to produce valid estimates (otherwise subsampling the
#       ## estimation unit, and cause estimates to be inflated substantially
#       combos <- data %>%
#         as.data.frame() %>%
#         group_by(.dots = grpBy) %>%
#         summarize()
#
#       # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
#       aGrpBy <- c('YEAR', grpBy[grpBy %in% names(PLOT) | grpBy %in% names(COND)]) # Can't group by columns specific to TREE (or made from)
#
#
#
#       for (i in 1:nrow(combos)){
#         if (ncol(combos) > 0){
#           # Update domain indicator for each each column speficed in grpBy
#           td = 1 # Start both at 1, update as we iterate through
#           ad = 1
#           for (n in 1:ncol(combos)){
#             # Tree domain indicator for each column in
#             tObs <- combos[i,][[grpBy[n]]] == data[[grpBy[n]]]
#             td <- landD * aD * tD * typeD * tObs * td
#             # Area domain indicator for each column in
#             if(grpBy[n] %in% aGrpBy){
#               aObs <- combos[i,][[aGrpBy[n]]] == data[[aGrpBy[n]]]
#               aObs[is.na(aObs)] <- 0
#               ad <- landD * aD * aObs * ad
#             }
#           }
#           data$tDI <- td
#           data$tDI[is.na(data$tDI)] <- 0
#           data$fullDI <- ad
#           data$fullDI[is.na(data$aDI)] <- 0
#         } # Close DI loop & conditional
#
#         # Compute estimates
#         a <- data %>%
#           group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
#           # filter(tDI > 0) %>%
#           summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
#                     aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0),
#                     a = first(AREA_USED),
#                     p1EU = first(P1PNTCNT_EU),
#                     p1 = first(P1POINTCNT),
#                     p2 = first(P2POINTCNT),
#                     aAdj = first(aAdj),
#                     fullDI = first(fullDI)) %>%
#           group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#           summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
#                     faFull = sum(CONDPROP_UNADJ * fullDI * aAdj, na.rm = TRUE),
#                     plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
#                     p1EU = first(p1EU),
#                     a = first(a),
#                     p1 = first(p1),
#                     p2 = first(p2)) %>%
#           # Continue through totals
#           #d <- dInt %>%
#           #filter(plotIn > 0) %>%
#           group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
#           summarize(aStrat = mean(fa, na.rm = TRUE),
#                     fullStrat = mean(faFull, na.rm = TRUE),
#                     plotIn = sum(plotIn, na.rm = TRUE),
#                     a = first(a),
#                     w = first(p1) / first(p1EU), # Stratum weight
#                     nh = first(p2), # Number plots in stratum
#                     # Strata level variances
#                     av = ifelse(first(ESTN_METHOD == 'simple'),
#                                 var(fa * first(a) / nh),
#                                 (sum(fa^2, na.rm = TRUE) - sum(nh * aStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
#                     fullv = ifelse(first(ESTN_METHOD == 'simple'),
#                                    var(faFull * first(a) / nh),
#                                    (sum(faFull^2, na.rm = TRUE) - sum(nh * fullStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
#                     # cvStrat_t = ifelse(first(ESTN_METHOD == 'simple'),
#                     #                    cov(fa,tPlot),
#                     #                    (sum(fa*tPlot) - sum(nh * aStrat *tStrat)) / (nh * (nh-1))), # Stratified and double cases
#                     cvStrat = ifelse(first(ESTN_METHOD == 'simple'),
#                                        cov(faFull,fa),
#                                        (sum(faFull*fa, na.rm = TRUE) - sum(nh * aStrat *fullStrat, na.rm = TRUE)) / (nh * (nh-1)))) %>% # Stratified and double cases) %>% # Stratified and double cases
#           group_by(ESTN_UNIT_CN) %>%
#           summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
#                     fullEst = unitMean(ESTN_METHOD, a, nh, w, fullStrat),
#                     plotIn = sum(plotIn, na.rm = TRUE),
#                     # Estimation of unit variance,
#                     aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
#                     fullVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, fullv, fullStrat, fullEst),
#                     cvEst = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat, aStrat, aEst, fullStrat, fullEst)) %>%
#           # Compute totals
#           summarize(AREA = sum(aEst, na.rm = TRUE),
#                     TOTAL_AREA = sum(fullEst, na.rm = TRUE),
#                     PERC_AREA = AREA / TOTAL_AREA * 100,
#                     nPlots = sum(plotIn, na.rm = TRUE),
#                     areaVar = sum(aVar, na.rm = TRUE),
#                     fVar = sum(fullVar, na.rm = TRUE),
#                     cv = sum(cvEst, na.rm = TRUE),
#                     propVar = (1/TOTAL_AREA^2) * (areaVar + (PERC_AREA^2 * fVar) - 2 * PERC_AREA * cv),
#                     AREA_SE = ifelse(nPlots > 1, sqrt(areaVar) / AREA * 100,0),
#                     PERC_SE = ifelse(nPlots > 1, sqrt(propVar) / PERC_AREA * 100,0),
#                     TOTAL_SE = ifelse(nPlots > 1, sqrt(fVar) / TOTAL_AREA * 100,0)) %>%
#           select(AREA, PERC_AREA, TOTAL_AREA, AREA_SE, PERC_SE, TOTAL_SE, nPlots)
#
#
#         ## Add each row onto the output dataset
#         if(i > 1){
#           aOut <- bind_rows(aOut, a)
#         } else { # On the first iteration
#           aOut <- a
#         }
#       } # End Combo loop
#
#       # If grpBy was specified, rejoin with these names
#       if (ncol(combos) > 0){
#         aOut <- bind_cols(combos, aOut) #%>%
#         #filter(S_a > 0)
#       }
#
#       # Return a spatial object
#       if (!is.null(polys) & returnSpatial) {
#         suppressMessages({aOut <- left_join(polys, aOut) %>%
#           select(c(grpByOrig, 'AREA', 'PERC_AREA', 'TOTAL_AREA', 'AREA_SE', 'PERC_SE', 'TOTAL_SE', 'nPlots', names(polys)))})
#       } else if (!is.null(polys) & returnSpatial == FALSE){
#         aOut <- select(aOut, c(grpByOrig, 'AREA', 'PERC_AREA', 'TOTAL_AREA', 'AREA_SE', 'PERC_SE', 'TOTAL_SE', 'nPlots', everything()))
#       }
#     } else {# End SE == FALSE
#       ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
#       # IF they want tree level variables, have to switch it up a bit
#       # if(!is.null(treeDomain) | bySpecies | bySizeClass){
#       #   data <- data %>%
#       #     distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, EVALID, COND_STATUS_CD, .keep_all = TRUE)
#       # } else {
#       #   data <- data %>%
#       #     distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, EVALID, COND_STATUS_CD, .keep_all = TRUE)
#       # }
#       data$fullDI <- landD * aD
#
#       # Compute estimates
#       aSub <- data %>%
#         group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
#         # filter(tDI > 0) %>%
#         summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
#                   aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0),
#                   a = first(AREA_USED),
#                   p1EU = first(P1PNTCNT_EU),
#                   p1 = first(P1POINTCNT),
#                   p2 = first(P2POINTCNT),
#                   aAdj = first(aAdj),
#                   EXPNS = EXPNS) %>%
#         group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#         summarize(fa = sum(CONDPROP_UNADJ * aDI * EXPNS * aAdj, na.rm = TRUE),
#                   plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
#         group_by(.dots = grpBy) %>%
#         summarize(AREA = sum(fa, na.rm = TRUE),
#                   nPlots = sum(plotIn, na.rm = TRUE))
#
#       suppressMessages({
#         aOut <- data %>%
#           group_by(YEAR, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#           summarize(fa = sum(CONDPROP_UNADJ * fullDI * EXPNS * aAdj, na.rm = TRUE)) %>%
#           group_by(YEAR) %>%
#           summarize(TOTAL_AREA = sum(fa, na.rm = TRUE)) %>%
#           inner_join(aSub) %>%
#           mutate(PERC_AREA = AREA / TOTAL_AREA * 100) %>%
#           select(grpBy, AREA, PERC_AREA, TOTAL_AREA, nPlots)
#       })
#
#
#       # Return a spatial object
#       if (!is.null(polys) & returnSpatial) {
#         suppressMessages({aOut <- left_join(polys, aOut) %>%
#           select(c(grpByOrig, 'AREA', 'PERC_AREA', 'TOTAL_AREA', 'nPlots', names(polys)))})
#       } else if (!is.null(polys) & returnSpatial == FALSE){
#         aOut <- select(aOut, c(grpByOrig, 'AREA', 'PERC_AREA', 'TOTAL_AREA', 'nPlots', everything()))
#       }
#
#     }
#
#   }
#
#   return(aOut)
# }



# # Compute mean stand age attributes, add canopy cover too I guess
# age <- function(db = NULL,
#                      landType = 'Forest',
#                      treeType = 'Live',
#                      treeDomain = NULL,
#                      areaDomain = NULL,
#                      grpBy = NULL,
#                      polys = NULL,
#                      returnSpatial = FALSE,
#                      byPlot = FALSE,
#                      bySizeClass = FALSE,
#                      SE = TRUE){
# }


## make avStand & avTree functions, replace age above








