
#### A Start up Message ------------
.onAttach <- function(libname, pkgname){
  #packageStartupMessage('Download FIA Data Here: https://apps.fs.usda.gov/fia/datamart/datamart.html')
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

  ba <- case_when(
    is.na(diameter) ~ NA_real_,
    ## Growth accounting only
    diameter < 0 ~ -(diameter^2 * .005454),
    TRUE ~ diameter^2 * .005454)


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

stratVar <- function(ESTN_METHOD, x, xStrat, ndif, a, nh, y = NULL, yStrat = NULL){
  ## Variance
  if (is.null(y)){
    v <- ifelse(first(ESTN_METHOD == 'simple'),
                var(c(x, numeric(ndif)) * first(a) / nh),
                (sum(c(x, numeric(ndif))^2) - sum(nh * xStrat^2)) / (nh * (nh-1)))
    ## Covariance
  } else {
    v <- ifelse(first(ESTN_METHOD == 'simple'),
                cov(x,y),
                (sum(x*y, na.rm = TRUE) - sum(nh * xStrat *yStrat)) / (nh * (nh-1)))
  }
}

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

  if (length(x) > 1){
    print(sapply(x, as_tibble))
  } else {
    print(as_tibble(x[1]))
  }
}

#' @import dplyr
#' @import methods
#' @import sf
#' @import stringr
#' @import gganimate
#' @import ggplot2
#' @import bit64
#' @import progress
#' @importFrom data.table fread fwrite rbindlist
#' @importFrom parallel makeCluster detectCores mclapply parLapply stopCluster clusterEvalQ
#' @import tidyr
#' @importFrom sp over proj4string<- coordinates<- spTransform proj4string
#' @importFrom stats cov var
#' @importFrom utils object.size read.csv tail globalVariables type.convert download.file unzip
NULL

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


################### FIA FUNCTIONS #########################
# Read in FIA database files (.csv) from local directory
#' @export
readFIA <- function(dir,
                    common = TRUE,
                    tables = NULL,
                    nCores = 1,
                    ...){

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
      files <- files[str_sub(files,4,-5) %in% tables]
    } else {
      files <- files[str_sub(files,1,-5) %in% tables]
    }
  }

  # Only csvs
  files <- files[str_sub(files,-4,-1) == '.csv']

  # Only extract the tables needed to run functions in rFIA
  if (common){
    cFiles <- c('COND', 'COND_DWM_CALC', 'INVASIVE_SUBPLOT_SPP', 'PLOT', 'POP_ESTN_UNIT',
                'POP_EVAL', 'POP_EVAL_GRP', 'POP_EVAL_TYP', 'POP_PLOT_STRATUM_ASSGN', 'POP_STRATUM',
                'SUBPLOT', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT', 'TREE_GRM_BEGIN', 'SUBP_COND_CHNG_MTRX',
                'SEEDLING', 'SURVEY')
    if (any(str_sub(files, 3, 3) == '_')){
      files <- files[str_sub(files,4,-5) %in% cFiles]
    } else{
      files <- files[str_sub(files,1,-5) %in% cFiles]
    }
  }


  # ## Compute estimates in parallel -- Clusters in windows, forking otherwise
  # if (Sys.info()['sysname'] == 'Windows'){
  #   cl <- makeCluster(nCores) # Set up snow cluster
  #   inTables <- parLapply(cl, X = files, fun = readFIAHelper1, dir)
  # } else { # Unix systems
  #   inTables <- mclapply(files, FUN = readFIAHelper1, dir, mc.cores = nCores)
  # }
  inTables <- list()
  for (n in 1:length(files)){
    # Read in and append each file to a list
    file <- fread(paste(dir, files[n], sep = ""), showProgress = FALSE, integer64 = 'double', logical01 = FALSE, nThread = nCores, ...)
    # We don't want data.table formats
    #file <- as.data.frame(file)
    fileName <- str_sub(files[n], 1, -5)

    inTables[[fileName]] <- file
  }


  # Give them some names
  #names(inTables) <- files
  #inTables <- lapply(inTables, as.data.frame)

  outTables <- list()
  if (any(str_sub(names(inTables), 3, 3) == '_')){ ## STATE NAMING CONVENTION
    # Remove the state prefix
    names(inTables) <- str_sub(names(inTables), 4, -1)
  }
  uniqueNames <- unique(names(inTables))
  ## Works regardless of whether or not there are duplicate names (multiple states)
  for (i in 1:length(uniqueNames)){
    outTables[[uniqueNames[i]]] <- rbindlist(inTables[names(inTables) == uniqueNames[i]])
  }

  # NEW CLASS NAME FOR FIA DATABASE OBJECTS
  outTables <- lapply(outTables, as.data.frame)
  class(outTables) <- 'FIA.Database'

  ## If you are on windows, close explicitly
  #closeAllConnections()

  return(outTables)
}

## Access FIA Database files from the FIA Datamart
#' @export
getFIA <- function(states,
                   dir = NULL,
                   common = TRUE,
                   tables = NULL,
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

  ## Some warnings up front
  ## Do not try to merge ENTIRE with other states
  if (length(states) > 1 & any(str_detect(str_to_upper(states), 'ENTIRE'))){
    stop('Cannot merge ENITRE with other state tables. ENTIRE includes all state tables combined. Do you only need data for a particular region?')
  }
  ## Check to make sure states exist
  allStates <- c('AL', 'AK', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'HI', 'ID',
                 'IL', 'IN', 'IA', 'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 'MI', 'MN', 'MS',
                 'MO', 'MT', 'NE', 'NV', 'NH', 'NJ', 'NM', 'NY', 'NC', 'ND', 'OH', 'OK',
                 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WA', 'WV',
                 'WI', 'WY', 'AS', 'FM', 'GU', 'MP', 'PW', 'PR', 'VI', 'ENTIRE', 'REF')
  if (any(str_to_upper(states) %in% allStates == FALSE)){
    missStates <- states[str_to_upper(states) %in% allStates == FALSE]
    stop(paste('Data unavailable for: ', paste(as.character(missStates),collapse = ', '), '. Did you use state/terrority abbreviations? e.g. use states = "AL" (correct) instead of states = "ALABAMA".'))
  }

  ## Check to make sure tables exist
  allTables <- c("BOUNDARY", "COND_DWM_CALC", "COND","COUNTY","DWM_COARSE_WOODY_DEBRIS",
                 "DWM_DUFF_LITTER_FUEL","DWM_FINE_WOODY_DEBRIS","DWM_MICROPLOT_FUEL",
                 "DWM_RESIDUAL_PILE", "DWM_TRANSECT_SEGMENT", "DWM_VISIT","GRND_CVR",
                 "INVASIVE_SUBPLOT_SPP","LICHEN_LAB","LICHEN_PLOT_SUMMARY","LICHEN_VISIT",
                 "OZONE_BIOSITE_SUMMARY","OZONE_PLOT_SUMMARY","OZONE_PLOT","OZONE_SPECIES_SUMMARY",
                 "OZONE_VALIDATION","OZONE_VISIT", "P2VEG_SUBP_STRUCTURE","P2VEG_SUBPLOT_SPP",
                 "PLOT_REGEN","PLOT", "PLOTGEOM", "PLOTSNAP","POP_ESTN_UNIT","POP_EVAL_ATTRIBUTE",
                 "POP_EVAL_GRP","POP_EVAL_TYP","POP_EVAL","POP_PLOT_STRATUM_ASSGN","POP_STRATUM",
                 "SEEDLING_REGEN","SEEDLING","SITETREE","SOILS_EROSION","SOILS_LAB","SOILS_SAMPLE_LOC" ,
                 "SOILS_VISIT", "SUBP_COND_CHNG_MTRX","SUBP_COND","SUBPLOT_REGEN","SUBPLOT",
                 "SURVEY","TREE_GRM_BEGIN","TREE_GRM_COMPONENT","TREE_GRM_ESTN", "TREE_GRM_MIDPT",
                 "TREE_REGIONAL_BIOMASS", "TREE_WOODLAND_STEMS","TREE","VEG_PLOT_SPECIES",
                 "VEG_QUADRAT","VEG_SUBPLOT_SPP","VEG_SUBPLOT", "VEG_VISIT",
                 'CITATION', 'DIFFERENCE_TEST_PER_ACRE', 'DIFFERENCE_TEST_TOTALS',
                 'FIADB_VERSION', 'FOREST_TYPE', 'FOREST_TYPE_GROUP',
                 'GRM_TYPE', 'HABTYP_DESCRIPTION', 'HABTYP_PUBLICATION',
                 'INVASIVE_SPECIES', 'LICHEN_SPECIES', 'LICHEN_SPP_COMMENTS',
                 'NVCS_HEIRARCHY_STRCT', 'NVCS_LEVEL_1_CODES', 'NVCS_LEVEL_2_CODES',
                 'NVCS_LEVEL_3_CODES', 'NVCS_LEVEL_4_CODES', 'NVCS_LEVEL_5_CODES',
                 'NVCS_LEVEL_6_CODES', 'NVCS_LEVEL_7_CODES', 'NVCS_LEVEL_8_CODES',
                 'OWNGRPCD', 'PLANT_DICTIONARY', 'POP_ATTRIBUTE', 'POP_EVAL_TYP_DESCR',
                 'RESEARCH_STATION', 'SPECIES', 'SPECIES_GROUP', 'STATE_ELEV', 'UNIT')
  if (any(str_to_upper(tables) %in% allTables == FALSE)){
    missTables <- tables[str_to_upper(tables) %in% allTables == FALSE]
    stop(paste('Tables: ', paste(as.character(missTables),collapse = ', '), ' unavailble. Check the FIA Datamart at https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html for a list of available tables for each state. Alternatively, specify common = TRUE to download the most commonly used tables.

Did you accidentally include the state abbreviation in front of the table name? e.g. tables = "AL_PLOT" (wrong) instead of tables = "PLOT" (correct).'))
  }

  ## If individual tables are specified, then just grab those .csvs, otherwise download the .zip file, extract and read with fread. Should be quite a bit quicker.
  if (!is.null(tables)){
    ## Make a list of tables names to read in
    ## Append table names with state abbs and then add url link
    tables <- str_to_upper(tables)

    # Make sure state Abbs are in right format
    states <- str_to_upper(states)
    if ('ENTIRE' %in% states == FALSE) {
      states <- paste0(states, '_')
    } else {
      states <- ''
    }


    tblNames <- sapply(states, FUN = function(x, y){paste0(x,y,'.zip')}, y = tables)
    tblNames <- c(tblNames)
    urls <- paste0('https://apps.fs.usda.gov/fia/datamart/CSV/', tblNames)


    inTables = list()
    for (n in 1:length(urls)){

      newName <- paste0(str_sub(tblNames[n], 1, -5), '.csv')

      ## Download the zip to a temporary file
      temp <- tempfile()
      download.file(urls[n], temp)

      # Write the data out the directory they've chosen
      if(is.null(dir)){
        unzip(temp, exdir = tempDir)
        file <- fread(paste0(tempDir, '/', newName), showProgress = FALSE, logical01 = FALSE, integer64 = 'double', nThread = nCores)
      } else {
        #download.file(urls[n], paste0(dir, tblNames[n]))
        unzip(temp, exdir = str_sub(dir, 1, -2))
        file <- fread(paste0(dir, newName), showProgress = FALSE, logical01 = FALSE, integer64 = 'double', nThread = nCores)
      }

      unlink(temp)

      # We don't want data.table formats
      file <- as.data.frame(file)

      inTables[[str_sub(urls[n], 43, -5)]] <- file
    }

    # Check for corresponding tables (multiple states)
    # If they exist, loop through and merge corresponding tables
    ## Check if the directory has the entire US naming convention or state naming convention
    tableNames <- names(inTables)
    outTables <- list()
    if (any(str_sub(tableNames, 3, 3) == '_')){ ## STATE NAMING CONVENTION
      if (anyDuplicated(str_sub(tableNames, 4)) != 0){
        for (i in 1:length(unique(str_sub(tableNames, 4)))){
          subList <- inTables[str_sub(tableNames, 4) == unique(str_sub(tableNames,4))[i]]
          name <- unique(str_sub(tableNames, 4))[i]
          # Give a ton of warnings about factors and characters, don't do that
          outTables[[name]] <- suppressWarnings(do.call(rbind, subList))
        }
      } else {
        outTables <- inTables
        names(outTables) <- unique(str_sub(tableNames, 4))
      }
    } else{ ## ENTIRE NAMING CONVENTION
      if (anyDuplicated(tableNames) != 0){
        for (i in 1:length(unique(tableNames))){
          subList <- inTables[tableNames == unique(tableNames)[i]]
          name <- unique(str_sub(tableNames, 1))[i]
          # Give a ton of warnings about factors and characters, don't do that
          outTables[[name]] <- suppressWarnings(do.call(rbind, subList))
        }
      } else {
        outTables <- inTables
        names(outTables) <- unique(str_sub(tableNames, 1))
      }
    }
    # NEW CLASS NAME FOR FIA DATABASE OBJECTS
    #outTables <- lapply(outTables, as.data.frame)
    unlink(tempDir, recursive = TRUE)
    class(outTables) <- 'FIA.Database'
    #### DOWNLOADING THE WHOLE ZIP FILE
  } else {
    ## Download to a temporary file location, then extract to the permanent if wanted. Read extracted tables into R w/ readFIA.

    # Make sure state Abbs are in right format
    states <- str_to_upper(states)

    ## Download each state and extract to directory
    for (i in 1:length(states)){
      # Temporary directory to download to
      temp <- paste0(tempDir, '/', states[i],'.zip') #tempfile()
      ## Make the URL
      url <- paste0('https://apps.fs.usda.gov/fia/datamart/CSV/', states[i],'.zip')
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

    ## Read in the files w/ readFIA
    if (is.null(dir)){
      outTables <- readFIA(tempDir, nCores = nCores, common = common)
      #unlink(tempDir, recursive = TRUE)
    } else {
      outTables <- readFIA(dir, nCores = nCores, common = common)
    }

    ## If you are on windows, close explicitly
    #closeAllConnections()
    #unlink(temp)

  }
  #unlink(tempDir, recursive = TRUE)
  #closeAllConnections()

  if (is.null(dir)){
     tmp <- list.files(tempDir, full.names = TRUE, pattern = '.csv')
     invisible(file.remove(tmp))
    }

  return(outTables)

}

## Write out the raw FIA files
#' @export
writeFIA <- function(db,
                     dir,
                     nCores = 1,
                     ...){
  #cat(sys.call()$dir)
  if (!is.null(dir)){
    # Add a slash to end of directory name if missing
    if (str_sub(dir,-1) != '/'){
      dir <- paste(dir, '/', sep = "")
    }
    # Check to see directory exists
    if(!dir.exists(dir)) {
      stop(paste('Directory', dir, 'does not exist. Cannot create new directory.'))
    }
  }

  tableNames <- names(db)
  ## Write out tables
  ## Read/write tables in parallel -- Clusters in windows, forking otherwise
  # if (Sys.info()['sysname'] == 'Windows'){
  #   cl <- makeCluster(nCores) # Set up snow cluster
  #   parLapply(cl, X = tableNames, fun = writeFIAHelper, db, dir)
  # } else { # Unix systems
  #   mclapply(X = tableNames, FUN = writeFIAHelper, db, dir, mc.cores = nCores)
  # }
  for (i in 1:length(tableNames)){
    if (is.data.frame(db[[i]])){
      fwrite(x = db[[i]], file = paste0(dir, tableNames[i], '.csv'), showProgress = FALSE, nThread = nCores)
    }
  }

  ## If you are on windows, close explicitly
  closeAllConnections()

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
    left_join(select(db$POP_EVAL_TYP, c('EVAL_GRP_CN', 'EVAL_TYP')), by = 'EVAL_GRP_CN')

  if (!is.null(state)){
    state <- str_to_upper(state)
    ## Join state abbs with state codes in popeval
    ids <- left_join(ids, select(intData$EVAL_GRP, c('STATECD', 'STATE')), by = 'STATECD')
    # Check if any specified are missing from db
    if (any(unique(state) %in% unique(db$POP_EVAL_STATE) == FALSE)){
      missStates <- state[state %in% unique(db$POP_EVAL_STATE) == FALSE]
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
      mutate(place = str_to_upper(LOCATION_NM)) %>%
      group_by(place, EVAL_TYP) %>%
      summarize(END_INVYR = max(END_INVYR, na.rm = TRUE),
                LOCATION_NM = first(LOCATION_NM))
      #filter(END_INVYR == max(END_INVYR, na.rm = TRUE))

    ids <- left_join(maxYear, select(ids, c('LOCATION_NM', 'EVAL_TYP', 'END_INVYR', 'EVALID')), by = c('LOCATION_NM', 'EVAL_TYP', 'END_INVYR'))
  }

  # Output as vector
  ID <- unique(ids$EVALID)

  return(ID)
}

## Spatiotemporal clip of FIADB
#' @export
clipFIA <- function(db,
                    mostRecent = TRUE,
                    mask = NULL,
                    matchEval = FALSE,
                    evalid = NULL,
                    designCD = NULL) {

  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Databse". Use readFIA() to load your FIA data.')
  }
  if (!is.null(mask) & first(class(mask)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('mask must be spatial polygons object of class sp or sf. ')
  }

  ################### ADD UNIQUE ID TO PLOTS #############################
  db$PLOT <- db$PLOT %>%
    mutate(PLT_CN = CN) %>%
    #mutate(date = ymd(paste(MEASYEAR, MEASMON, MEASDAY, sep = '-'))) %>%
    mutate(pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_')) # Make a unique ID for each plot, irrespective of time

  if (!is.null(designCD)) db$PLOT <- filter(db$PLOT, DESIGNCD == any(as.integer(designCD)))


  ######### IF USER SPECIES EVALID (OR MOST RECENT), EXTRACT APPROPRIATE PLOTS ##########
  if (!is.null(evalid)){
    # Join appropriate tables and filter out specified EVALIDs
    tempData <- select(db$PLOT, CN, PREV_PLT_CN) %>%
      mutate(PLT_CN = CN) %>%
      inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID')), by = 'PLT_CN') %>%
      #inner_join(select(POP_ESTN_UNIT, c('EVAL_GRP_CN', 'EVALID')), by = 'EVAL_GRP_CN') %>%
      filter(EVALID %in% evalid)
    # Extract plots which relate to specified EVALID (previouy for change estimation)
    PPLOT <- db$PLOT[db$PLOT$CN %in% tempData$PREV_PLT_CN,]
    db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]
  }

  ## Locate the most recent EVALID and subset plots
  if (mostRecent){
    suppressWarnings({
    tempData <- select(db$PLOT, CN, PREV_PLT_CN) %>%
      mutate(PLT_CN = CN) %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'EVALID')), by = c('PLT_CN'))


    ## Most recent evals by
    mrids <- findEVALID(db, mostRecent = TRUE)

    tempData <- tempData %>%
      filter(EVALID %in% mrids)
    })

    # Extract appropriate PLOTS
    PPLOT <- db$PLOT[db$PLOT$CN %in% tempData$PREV_PLT_CN,]
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


  ##########################  SPATIAL-TEMPORAL INTERSECTION ######################################
  ## Oringally did this based solely on PLT_CN, but this method will not work when we try to compute variance
  ##  estimates for a region, because estimates are computed at ESTN_UNIT level. Now the spatial functionality
  ##  of clipFIA is primarily intended to reduce the amount of data which must be held in RAM. As long as all
  ##  plot data is avaliable for any intersecting ESTN_UNIT, we can cut back on memory requirements and still
  ##  produce valid estimates
   if (!is.null(mask)){
     # # Convert polygons to an sf object
     # mask <- mask %>%
     #   as('sf') %>%
     #   mutate(polyID = 1:nrow(mask))
     # ## Make plot data spatial, projected same as polygon layer
     # pltSF <- select(db$PLOT, c('LON', 'LAT', pltID)) %>%
     #   filter(!is.na(LAT) & !is.na(LON)) %>%
     #   distinct(pltID, .keep_all = TRUE)
     # coordinates(pltSF) <- ~LON+LAT
     # proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
     # pltSF <- as(pltSF, 'sf') %>%
     #   st_transform(crs = st_crs(mask)$proj4string)
     #
     # ## Split up mask
     # polyList <- split(mask, as.factor(mask$polyID))
     #
     # suppressMessages({suppressWarnings({
     #   ## Intersect it
     #   out <- lapply(names(polyList), FUN = areal_par, pltSF, polyList)
     # })})
     #
     #
     #
     # pltSF <- bind_rows(out) %>%
     #   left_join(select(db$PLOT, PLT_CN, PREV_PLT_CN, pltID), by = 'pltID')
     #
     # PPLOT <- filter(db$PLOT, db$PLOT$PLT_CN %in% pltSF$PREV_PLT_CN)
     # db$PLOT <- filter(db$PLOT, db$PLOT$PLT_CN %in% pltSF$PLT_CN)

     ###OLD
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
     plts <- select(db$PLOT, PLT_CN, PREV_PLT_CN) %>%
       inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'STRATUM_CN')), by = 'PLT_CN') %>%
       inner_join(select(db$POP_STRATUM, c('CN', 'ESTN_UNIT_CN')), by = c('STRATUM_CN' = 'CN')) %>%
       filter(ESTN_UNIT_CN %in% estUnits$ESTN_UNIT_CN) %>%
       group_by(PLT_CN, PREV_PLT_CN) %>%
       summarize()

     # Clip out the above plots from the full database, will reduce size by a shit pile
     PPLOT <- filter(db$PLOT, db$PLOT$PLT_CN %in% plts$PREV_PLT_CN)
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

  ## IF no spatial or temporal clip was specified, make PPLOT NULL
  if (mostRecent == FALSE & is.null(evalid) & is.null(mask)) PPLOT <- NULL

  #################  APPLY  QUERY TO REMAINING TABLES  ###########################
  ## User gives object names and not list object
  otherTables <- db

  if (is.list(otherTables)){
    tableNames <- names(otherTables)
    clippedData <- list()

    # Query for all tables directly related to PLOT
    for (i in 1:length(otherTables)){
      # Pull the object with the name listed from the global environment
      table <- otherTables[[i]]
      name <- tableNames[i]

      if (name == "PLOT"){
        db$PLOT$prev = 0
        if(nrow(PPLOT) > 0) PPLOT$prev = 1
        clippedData[['PLOT']] <- rbind(db$PLOT, PPLOT)

      } else if (!is.null(db$OZONE_PLOT) & name == 'OZONE_PLOT'){
          clippedData[['OZONE_PLOT']] <- db$OZONE_PLOT

      } else if (name == 'TREE' | name == 'COND'){ # Need previous attributes
        clipTable <- table[table$PLT_CN %in% c(db$PLOT$CN, PPLOT$PLT_CN),]
        clippedData[[name]] <- clipTable

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
    if(!is.null(evalid) | mostRecent){ # User specified EVALIDs
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
                        nCores = 1) {

  ## Need a plotCN
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        0
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
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
  }

  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  ## Some warnings
  if (class(db) != "FIA.Database"){
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
  }
  if (!is.null(polys) & first(first(class(polys))) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (tidy & returnSpatial & !is.null(polys)){
    warning('Returning multiple observations for each areal unit. If returnSpatial = TRUE, tidy = FALSE is recommended.')
  }
  if (landType %in% c('timber', 'forest', 'all') == FALSE){
    stop('landType must be one of: "forest", "timber", or "all".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    grpBy <- c(grpBy, 'LON', 'LAT')
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

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)
    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))
  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', grpP, 'aD_p', 'sp')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', 'landD', 'aD_c', grpC)), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      #left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      #left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      left_join(select(db_clip$TREE, c('PLT_CN', 'CONDID', 'DIA', 'STATUSCD', 'CCLCD', 'TREECLCD', 'STANDING_DEAD_CD', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE')), by = c('PLT_CN', 'CONDID')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM) %>%
      mutate_if(is.factor,
                as.character)%>%
      filter(!is.na(YEAR)) %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE)

    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }

    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
    data$tDI <- data$landD * data$aD_p * data$aD_c * data$sp



    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {
      sOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(stage = structHelper(DIA, CCLCD),
                  nStems = length(which(tDI == 1)))

      if (returnSpatial){
        sOut <- sOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
    } else {
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))
      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
            library(tidyr)
          })
          sOut <- parLapply(cl, X = names(combos), fun = standStructHelper, combos, data, grpBy, tidy, totals, SE)
          stopCluster(cl)
        } else { # Unix systems
          sOut <- mclapply(X = names(combos), FUN = standStructHelper, combos, data, grpBy, totals, tidy, SE, mc.cores = nCores)
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
    out[[y]] <- sOut
    pb$tick()
  }
  sOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  sOut <- drop_na(sOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR) %>%
    as_tibble()

  ## Above converts to tibble
  if (returnSpatial) sOut <- st_sf(sOut)

  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) sOut <- unique(sOut)

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
                      nCores = 1) {
  ## Need a plotCN
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        0
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        inner_join(db$TREE, by = 'PLT_CN') %>%
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

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    grpBy <- c(grpBy, 'LON', 'LAT')
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

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)
    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))
  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'INVYR', 'MACRO_BREAKPOINT_DIA', grpP, 'aD_p', 'sp')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      left_join(select(db_clip$TREE, c('PLT_CN', 'CONDID', 'DIA','SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', grpT, 'typeD', 'tD')), by = c('PLT_CN', 'CONDID')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM)%>%
      mutate_if(is.factor,
                as.character)%>%
      filter(!is.na(YEAR)) %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE)

    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
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
      data <- data[!is.na(data$sizeClass),]
    }



    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {
      dOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, SUBP, CONDID, TREE, .keep_all = TRUE) %>%
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(H = divIndex(SPCD, TPA_UNADJ  * tDI, index = 'H'),
                  S = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S'),
                  Eh = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh'),
                  nStems = length(which(tDI == 1)))

      if (returnSpatial){
        dOut <- dOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
    } else {
      # combos <- select(data, c(grpBy)) %>%
      #   as.data.frame() %>%
      #   group_by(.dots = grpBy) %>%
      #   summarize() %>%
      #   filter(!is.na(YEAR))
      # if(!is.null(polys)){
      #   combos <- filter(combos, !is.na(polyID))
      # }
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
          })
          dOut <- parLapply(cl, X = names(combos), fun = diversityHelper, combos, data, grpBy, SE)
          stopCluster(cl)
        } else { # Unix systems
          dOut <- mclapply(X = names(combos), FUN = diversityHelper, combos, data, grpBy, SE, mc.cores = nCores)
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
    out[[y]] <- dOut
    pb$tick()
  }
  dOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  dOut <- drop_na(dOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR) %>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) dOut <- st_sf(dOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) dOut <- unique(dOut)
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
                  nCores = 1) {

  ## Need a plotCN
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        inner_join(db$TREE, by = 'PLT_CN') %>%
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

  # Save original grpBy for pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
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

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)

    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))


  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', grpP, 'aD_p', 'sp')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
        right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      # left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      # left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      left_join(select(db_clip$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', grpT, 'tD', 'typeD')), by = c('PLT_CN', 'CONDID')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM) %>%
      mutate_if(is.factor,
                as.character) %>%
      filter(!is.na(YEAR)) %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE)

    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
    data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
    data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp

    ## Add species to groups
    if (bySpecies) {
      data <- data %>%
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
      data$sizeClass <- makeClasses(data$DIA, interval = 2)
      data <- data[!is.na(data$sizeClass),]
    }

    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {

      #grpBy_sym <- syms(grpBy)
      tOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
                  BAA = sum(basalArea(DIA) * TPA_UNADJ * tDI, na.rm = TRUE),
                  TPA_PERC = TPA / sum(TPA_UNADJ * pDI, na.rm = TRUE) * 100,
                  BAA_PERC = BAA / sum(basalArea(DIA) * TPA_UNADJ * pDI, na.rm = TRUE) * 100,
                  nStems = length(which(tDI == 1)))

      # tOut <- data %>%
      #   #lazy_dt() %>%
      #   #filter(EVAL_TYP == 'EXPVOL') %>%
      #   distinct(PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      #   group_by(!!grpBy_sym, PLT_CN) %>%
      #   summarize(TPA = sum(TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
      #             BAA = sum(basalArea(DIA) * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
      #             nStems = length(which(tDI == 1))) #%>%
      #   #as_tibble()



      if (returnSpatial){
        tOut <- tOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
    } else {
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
      if (!is.null(polys)){
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
      } else {
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
      }

      #message('Computing Summary Statistics.....')
      #pb <- progress_bar$new(format = "[:bar] :percent eta: :eta", total = nrow(combos), clear = FALSE, width= 100)

      #for (i in 1:nrow(combos)){ #, .combine = 'rbind', .packages = 'dplyr', .export = c('data')
      #tOut <- foreach(i = 1:nrow(combos), .combine = 'rbind', .packages = 'dplyr', .export = c('data'))

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
          })
          tOut <- parLapply(cl, X = names(combos), fun = tpaHelper, combos, data, grpBy, aGrpBy, totals, SE)
          stopCluster(cl)
        } else { # Unix systems
          tOut <- mclapply(X = names(combos), FUN = tpaHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
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
          # Return spatial plots
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
    out[[y]] <- tOut
    pb$tick()
  }

  tOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  tOut <- drop_na(tOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR) %>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
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
                     nCores = 1) {

  ## Need a plotCN
  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        0
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        inner_join(db$TREE, by = 'PLT_CN') %>%
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
  }

  reqTables <- c('PLOT', 'TREE', 'TREE_GRM_COMPONENT', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
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

  ## No EXP_GROW available for Western States, make sure we warn that values will be returned as 0
  # These states do not allow temporal queries. Things are extremely weird with their eval groups
  noGrow <- c(02,03,04,07,08,11,14,15,16, 23, 30, 32, 35,43,49, 78)
  if(any(unique(db$PLOT$STATECD) %in% noGrow)){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    warning(paste('Recruitment data unavailable for: ', toString(fancyName) , '. Returning 0 for all recruitment estimates which include these states.', sep = ''))
  }
  # These states do not allow change estimates.
  if(any(unique(db$PLOT$STATECD) %in% c(69, 72, 78, 15, 02))){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% c(69, 72, 78, 15, 02)])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    stop(paste('Growth & Mortality Estimates unavailable for: ', paste(as.character(fancyName), collapse = ', '), sep = ''))
  }

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    grpBy <- c(grpBy, 'LON', 'LAT')
  }

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
                                      COMPONENT = SUBP_COMPONENT_AL_FOREST)

    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
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
    if (tolower(treeType) == 'all'){
      db$TREE$typeD <- 1
      ## Rename some variables in grm
      db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                      TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_TIMBER,
                                      TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_TIMBER,
                                      TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_TIMBER,
                                      SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_TIMBER,
                                      COMPONENT = SUBP_COMPONENT_AL_TIMBER)

    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
      db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                      TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_TIMBER,
                                      TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_TIMBER,
                                      TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_TIMBER,
                                      SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_TIMBER,
                                      COMPONENT = SUBP_COMPONENT_GS_TIMBER)
    }
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


  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID', 'GROWTH_ACCT'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID, .keep_all = TRUE)

    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID),
                ga = if_else(any(GROWTH_ACCT == 'Y'), 1, 0))
  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID', 'GROWTH_ACCT') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID, .keep_all = TRUE) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID),
                ga = if_else(any(GROWTH_ACCT == 'Y'), 1, 0))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    data <- select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db_clip$PLOT, c('PLT_CN', 'PREV_PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'REMPER', grpP, 'sp', 'aD_p')), by = 'PLT_CN') %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID','landD', 'aD_c', grpC)), by = 'PLT_CN') %>%
      left_join(select(db_clip$TREE, c('TRE_CN', 'PREV_TRE_CN', 'TREE', 'PLT_CN', 'CONDID', 'PREVCOND', 'SPCD', grpT, 'typeD', 'tD')), by = c('PLT_CN', 'CONDID')) %>% ## ISSUE MAY BE HERE, SEE EVALIDATOR CODE
      # GRM
      left_join(select(db_clip$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      mutate(tAdj = grmAdj(SUBPTYP_GRM, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM)
    # Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }


    ## MODIFY FOR GROWTH ACCOUNTING
    if (ids$ga[y]){
      # Only subplots from cond change matrix
      db_clip$SUBP_COND_CHNG_MTRX <- filter(db_clip$SUBP_COND_CHNG_MTRX, SUBPTYP == 1)

      # Previous attributes
      data <- data %>%
        left_join(select(db_clip$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'CONDID'), suffix = c('', '.subp')) %>%
        left_join(select(db_clip$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN.subp' = 'PLT_CN', 'PREVCOND.subp' = 'CONDID'), suffix = c('', '.chng')) %>%
        left_join(select(db_clip$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
        left_join(select(db_clip$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
        left_join(select(db_clip$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
        mutate_if(is.factor,
                  as.character) %>%
        mutate(aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.chng == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
               tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0))

      #If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
      data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
      data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
      data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
      data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
      data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
      data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

      ## Comprehensive indicator function
      data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp * data$aChng
      data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng
      #data$aDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$sp.prev * data$aChng
      # data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp * data$tChng

      ## If growth accounting, summing proportions at subplot level, otherwise full plot
      chngAdj <- .25

      ## No growth accounting
    } else {
      # Previous attributes
      data <- data %>%
        left_join(select(db_clip$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
        left_join(select(db_clip$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC)), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
        left_join(select(db_clip$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
        mutate_if(is.factor,
                  as.character) %>%
        ## Rename CONDPROP_UNADJ for consistency with above
        rename(SUBPTYP_PROP_CHNG = CONDPROP_UNADJ)
      ## Plot level proportion
      chngAdj <- 1

      # If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
      data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
      data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
      data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
      data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
      data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
      data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

      ## Comprehensive indicator function
      data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
      data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev

    }



    ## Add species to groups
    if (bySpecies) {
      data <- data %>%
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
      data$sizeClass <- makeClasses(data$DIA, interval = 2)
      data <- data[!is.na(data$sizeClass),]
    }



    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {
      tOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
        #filter(EVALID %in% tID) %>%
        # Compute estimates at plot level
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(TOTAL_TPA = sum(TPAGROW_UNADJ * tDI, na.rm = TRUE),
                  RECR_TPA = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'] / REMPER[COMPONENT == 'INGROWTH'], na.rm = TRUE),
                  MORT_TPA = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
                  REMV_TPA = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
                  RECR_PERC = RECR_TPA / TOTAL_TPA * 100,
                  MORT_PERC = MORT_TPA / TOTAL_TPA * 100,
                  REMV_PERC = REMV_TPA / TOTAL_TPA * 100,
                  nStems = length(which(tDI == 1)))

      if (returnSpatial){
        tOut <- tOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

        }

    } else {
      # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }

      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
      if (!is.null(polys)){
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND) | grpBy %in% names(pltSF)])
      } else {
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND)])
      }

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
          })
          tOut <- parLapply(cl, X = names(combos), fun = growMortHelper, combos, data, grpBy, aGrpBy, totals, SE, chngAdj)
          stopCluster(cl)
        } else { # Unix systems
          tOut <- mclapply(X = names(combos), FUN = growMortHelper, combos, data, grpBy, aGrpBy, totals, SE, chngAdj, mc.cores = nCores)
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
                           "nPlots_RECR"= rep(NA, nrow(combos)), "nPlots_MORT"= rep(NA, nrow(combos)),
                           "nPlots_REMV"= rep(NA, nrow(combos)), "nPlots_AREA" = rep(NA, nrow(combos)))
        if (!is.null(polys) & returnSpatial) {
          suppressMessages({suppressWarnings({
            polys = left_join(polys, combos)
            tOut <- left_join(polys, tOut)})})
        } else if (!is.null(polys) & returnSpatial == FALSE){
          tOut <- data.frame(select(tOut, -c('YEAR')), combos)
        }
      }
    } # End byPlot == FALSE
    out[[y]] <- tOut
    pb$tick()
  }
  tOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  tOut <- drop_na(tOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR)%>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
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
                       nCores = 1){
  ## Need a plotCN
  db$TREE <- db[['TREE']] %>% mutate(TRE_CN = CN)
  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        0
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        inner_join(db$TREE, by = 'PLT_CN') %>%
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
  }

  reqTables <- c('PLOT', 'TREE', 'TREE_GRM_COMPONENT', 'TREE_GRM_MIDPT', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
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

  ## No EXP_GROW available for most Western States, make sure we warn that values will be returned as 0
  # These states do not allow temporal queries. Things are extremely weird with their eval groups
  noGrow <- c(02,03,04,07,08,11,14,15,16, 23, 30, 32, 35,43,49, 78)
  if(any(unique(db$PLOT$STATECD) %in% noGrow)){
    vState <- unique(db$PLOT$STATECD[db$PLOT$STATECD %in% noGrow])
    fancyName <- unique(intData$EVAL_GRP$STATE[intData$EVAL_GRP$STATECD %in% vState])
    warning(paste('Growth data unavailable for: ', toString(fancyName) , '. Returning 0 for all recruitment estimates which include these states.', sep = ''))
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

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    grpBy <- c(grpBy, 'LON', 'LAT')
  }


  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
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
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
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
    if (tolower(treeType) == 'all'){
      db$TREE$typeD <- 1
      ## Rename some variables in grm
      db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                      TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_AL_TIMBER,
                                      TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_AL_TIMBER,
                                      TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_AL_TIMBER,
                                      SUBPTYP_GRM = SUBP_SUBPTYP_GRM_AL_TIMBER,
                                      COMPONENT = SUBP_COMPONENT_AL_TIMBER)

    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5, 1, 0)
      db$TREE_GRM_COMPONENT <- rename(db$TREE_GRM_COMPONENT,
                                      TPAMORT_UNADJ = SUBP_TPAMORT_UNADJ_GS_TIMBER,
                                      TPAREMV_UNADJ = SUBP_TPAREMV_UNADJ_GS_TIMBER,
                                      TPAGROW_UNADJ = SUBP_TPAGROW_UNADJ_GS_TIMBER,
                                      SUBPTYP_GRM = SUBP_SUBPTYP_GRM_GS_TIMBER,
                                      COMPONENT = SUBP_COMPONENT_GS_TIMBER)
    }
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
  rm(aD_p)

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


  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID', 'GROWTH_ACCT'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPGROW') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID, .keep_all = TRUE)
    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID),
                ga = if_else(any(GROWTH_ACCT == 'Y'), 1, 0))%>%
      ## Need growth accounting
      filter(ga == 1)
  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID', 'GROWTH_ACCT') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPGROW') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID, .keep_all = TRUE) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID),
                ga = if_else(any(GROWTH_ACCT == 'Y'), 1, 0)) %>%
      ## Need growth accounting
      filter(ga == 1)
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    # Only subplots from cond change matrix
    db_clip$SUBP_COND_CHNG_MTRX <- filter(db_clip$SUBP_COND_CHNG_MTRX, SUBPTYP == 1)

    data <- select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db_clip$PLOT, c('PLT_CN', 'PREV_PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', grpP, 'sp', 'aD_p', 'REMPER')), by = 'PLT_CN') %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID','landD', 'aD_c', grpC)), by = 'PLT_CN') %>%
      left_join(select(db_clip$TREE, c('TRE_CN', 'PREV_TRE_CN', 'PLT_CN', 'TREE', 'DIA', 'DRYBIO_AG', 'VOLCFNET', 'VOLCSNET', 'CONDID', 'PREVCOND', 'SPCD', grpT, 'typeD', 'tD')), by = c('PLT_CN', 'CONDID')) %>%
      # GRM
      left_join(select(db_clip$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
      left_join(select(db_clip$TREE_GRM_MIDPT, c("TRE_CN", 'DIA', 'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG')), by = 'TRE_CN', suffix = c('', '.mid')) %>%
      left_join(select(db_clip$TREE_GRM_BEGIN, c("TRE_CN", 'DIA', 'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG')), by = 'TRE_CN', suffix = c('', '.beg')) %>%
      # Condition change attributes
      left_join(select(db_clip$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'CONDID'), suffix = c('', '.subp')) %>%
      left_join(select(db_clip$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN.subp' = 'PLT_CN', 'PREVCOND.subp' = 'CONDID'), suffix = c('', '.chng')) %>%
      # Previous attributes
      left_join(select(db_clip$PLOT, c('PLT_CN',  grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDID',  'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
      left_join(select(db_clip$TREE, c('TRE_CN', 'DIA',  'DRYBIO_AG', 'VOLCFNET', 'VOLCSNET', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
      ## Housekeeping and adjustment factors
      mutate_if(is.factor,
                as.character) %>%
      mutate(aChng = ifelse(COND_STATUS_CD.prev == 1 & COND_STATUS_CD == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
             tChng = ifelse(COND_STATUS_CD.prev == 1 & COND_STATUS_CD == 1, 1, 0)) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      mutate(tAdj = grmAdj(SUBPTYP_GRM, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM) %>%
      ## Remove duplicates
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, SUBP, CONDID, TRE_CN, .keep_all = TRUE)

    #If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
    data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
    data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
    data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
    data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
    data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
    data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp * data$aChng
    data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng

    # ## Comprehensive indicator function
    # data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp * data$aChng
    # data$tDI2 <- data$landD * data$aD_p * data$aD_c * data$sp * data$tD * data$typeD * data$tChng
    # data$tDI1 <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$sp.prev * data$tD.prev * data$typeD.prev * data$tChng

    xtraGrp = NULL
    ## Add species to groups
    if (bySpecies) {
      data <- data %>%
        left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
        mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' '))%>%
        mutate_if(is.factor,
                  as.character)
      grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
      grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
      xtraGrp <- c('SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    }

    ## Break into size classes
    if (bySizeClass){
      grpBy <- c(grpBy, 'sizeClass')
      grpByOrig <- c(grpByOrig, 'sizeClass')
      data$sizeClass <- makeClasses(data$DIA, interval = 2)
      data <- data[!is.na(data$sizeClass),]
      xtraGrp <- c(xtraGrp, 'sizeClass')
    }

    ## Modify  attributes depending on component (mortality uses midpoint)
    data <- data %>%
      mutate(DIA2 = vrAttHelper(DIA, DIA.prev, DIA.mid, DIA.beg, COMPONENT, REMPER, 2) * tDI,
             DIA1 = vrAttHelper(DIA, DIA.prev, DIA.mid, DIA.beg, COMPONENT, REMPER, 1) * tDI,
             BA2 = vrAttHelper(basalArea(DIA), basalArea(DIA.prev), basalArea(DIA.mid), basalArea(DIA.beg), COMPONENT, REMPER, 2) * tDI,
             BA1 = vrAttHelper(basalArea(DIA), basalArea(DIA.prev), basalArea(DIA.mid), basalArea(DIA.beg), COMPONENT, REMPER, 1) * tDI,
             VOLCFNET2 = vrAttHelper(VOLCFNET, VOLCFNET.prev, VOLCFNET.mid, VOLCFNET.beg, COMPONENT, REMPER, 2) * tDI,
             VOLCFNET1 = vrAttHelper(VOLCFNET, VOLCFNET.prev, VOLCFNET.mid, VOLCFNET.beg, COMPONENT, REMPER, 1) * tDI,
             DRYBIO_AG2 = vrAttHelper(DRYBIO_AG, DRYBIO_AG.prev, DRYBIO_AG.mid, DRYBIO_AG.beg, COMPONENT, REMPER, 2) * tDI,
             DRYBIO_AG1 = vrAttHelper(DRYBIO_AG, DRYBIO_AG.prev, DRYBIO_AG.mid, DRYBIO_AG.beg, COMPONENT, REMPER, 1) * tDI) %>%
      # ## Omitting columns where either time is an NA
      # mutate(DIA2 = case_when(is.na(DIA1) ~ NA_real_, TRUE ~ DIA2),
      #        DIA1 = case_when(is.na(DIA2) ~ NA_real_, TRUE ~ DIA1),
      #        BA2 = case_when(is.na(BA1) ~ NA_real_, TRUE ~ BA2),
      #        BA1 = case_when(is.na(BA2) ~ NA_real_, TRUE ~ BA1),
      #        VOLCFNET2 = case_when(is.na(VOLCFNET1) ~ NA_real_, TRUE ~ VOLCFNET2),
      #        VOLCFNET1 = case_when(is.na(VOLCFNET2) ~ NA_real_, TRUE ~ VOLCFNET1),
      #        DRYBIO_AG2 = case_when(is.na(DRYBIO_AG1) ~ NA_real_, TRUE ~ DRYBIO_AG2),
      #        DRYBIO_AG1 = case_when(is.na(DRYBIO_AG2) ~ NA_real_, TRUE ~ DRYBIO_AG1)) %>%
      select(-c(DIA.mid, VOLCFNET.mid, VOLCSNET.mid, DRYBIO_AG.mid,
                DIA.prev, VOLCFNET.prev, VOLCSNET.prev, DRYBIO_AG.prev,
                DIA.beg, VOLCFNET.beg, VOLCSNET.beg, DRYBIO_AG.beg,
                DIA, VOLCFNET, VOLCSNET, DRYBIO_AG))


    ## Just what we need
    data <- data %>%
      select(YEAR, ESTN_UNIT_CN, STRATUM_CN, PLT_CN, TRE_CN, SUBP, CONDID, ESTN_METHOD, EXPNS, INVYR, TREE, aDI, tDI, SUBPTYP_PROP_CHNG,
             grpP, grpC, grpT, xtraGrp, TPAGROW_UNADJ, tAdj, aAdj, AREA_USED, P1POINTCNT, P1PNTCNT_EU, P2POINTCNT,
             DIA2, DIA1, BA2, BA1, DRYBIO_AG2, DRYBIO_AG1, VOLCFNET2, VOLCFNET1) %>%
      ## Dropping NA columns
      #drop_na(SUBPTYP_PROP_CHNG) %>%
      ## Rearrange previous values as observations
      pivot_longer(cols = DIA2:VOLCFNET1,
                   names_to = c(".value", 'ONEORTWO'),
                   names_sep = -1)

    chngAdj <- .25

    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }

    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {
      ### Compute total TREES in domain of interest
      tOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
        group_by(.dots = grpBy, PLT_CN, SUBP, TREE) %>%
        summarize(d = sum(DIA * tDI, na.rm = TRUE),
                  ba = sum(BA * tDI, na.rm = TRUE),
                  baa = sum(TPAGROW_UNADJ * BA * tDI, na.rm = TRUE),
                  vol = sum(VOLCFNET * tDI, na.rm = TRUE),
                  volA = sum(TPAGROW_UNADJ * VOLCFNET * tDI, na.rm = TRUE),
                  bio = sum(DRYBIO_AG * tDI, na.rm = TRUE),
                  bioA = sum(TPAGROW_UNADJ * DRYBIO_AG * tDI, na.rm = TRUE)) %>%
        # Compute estimates at plot level
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(DIA_GROW = mean(d, na.rm = TRUE),
                  BA_GROW = mean(ba, na.rm = TRUE),
                  BAA_GROW = sum(baa, na.rm = TRUE),
                  NETVOL_GROW = mean(vol, na.rm = TRUE),
                  NETVOL_GROW_AC = sum(volA, na.rm = TRUE),
                  BIO_GROW = mean(vol, na.rm = TRUE),
                  BIO_GROW_AC = sum(volA, na.rm = TRUE),
                  nStems = length(which(!is.na(d))))

      if (returnSpatial){
        tOut <- tOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
    } else {

      # Unique combinations of specified grouping variables. Simply listing the grouping variables in estimation code below does not produce valid estimates. Have to
      ## produce a unique domain indicator for each individual output observation (ex. Red Oak in Ingham County) to produce valid estimates (otherwise subsampling the
      ## estimation unit, and cause estimates to be inflated substantially)

      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
      if (!is.null(polys)){
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND) | grpBy %in% names(pltSF)])
      } else {
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND)])
      }

      #suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
            library(tidyr)
          })
          tOut <- parLapply(cl, X = names(combos), fun = vitalRatesHelper, combos, data, grpBy, aGrpBy, totals, SE, chngAdj)
          stopCluster(cl)
        } else { # Unix systems
          tOut <- mclapply(X = names(combos), FUN = vitalRatesHelper, combos, data, grpBy, aGrpBy, totals, SE, chngAdj, mc.cores = nCores)
        }
      #})

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
    out[[y]] <- tOut
    pb$tick()
  }
  tOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  tOut <- drop_na(tOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR)%>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
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
                    nCores = 1) {
  ## Need a plotCN
  db$PLOT <- db[['PLOT']] %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        0
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        inner_join(db$TREE, by = 'PLT_CN') %>%
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

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf') %>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    grpBy <- c(grpBy, 'LON', 'LAT')
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

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)
    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))

  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'INVYR', 'MACRO_BREAKPOINT_DIA', grpP, 'sp', 'aD_p')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'landD', 'aD_c')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      left_join(select(db_clip$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', 'typeD', 'tD',
                                  'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG', 'DRYBIO_BG', 'CARBON_AG', 'CARBON_BG', grpT)), by = c('PLT_CN', 'CONDID')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM) %>%
      mutate_if(is.factor,
                as.character)%>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE)
    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }

    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
    data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp

    ## Add species to groups
    if (bySpecies) {
      data <- data %>%
        left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
        mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' '))%>%
        mutate_if(is.factor,
                  as.character)
      grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
      grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    }

    ## Break into size classes
    if (bySizeClass){
      grpBy <- c(grpBy, 'sizeClass')
      grpByOrig <- c(grpByOrig, 'sizeClass')
      data$sizeClass <- makeClasses(data$DIA, interval = 2)
      data <- data[!is.na(data$sizeClass),]
    }


    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {
      bOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
        # Compute estimates at plot level
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(NETVOL_ACRE = sum(VOLCFNET * TPA_UNADJ * tDI, na.rm = TRUE),
                  SAWVOL_ACRE = sum(VOLCSNET * TPA_UNADJ * tDI, na.rm = TRUE),
                  BIO_AG_ACRE = sum(DRYBIO_AG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                  BIO_BG_ACRE = sum(DRYBIO_BG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                  BIO_ACRE = sum(BIO_AG_ACRE, BIO_BG_ACRE, na.rm = TRUE),
                  CARB_AG_ACRE = sum(CARBON_AG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                  CARB_BG_ACRE = sum(CARBON_BG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                  CARB_ACRE = sum(CARB_AG_ACRE, CARB_BG_ACRE, na.rm = TRUE),
                  nStems = length(which(tDI == 1)))

      if (returnSpatial){
        bOut <- bOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
    } else {
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
      if (!is.null(polys)){
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND) | grpBy %in% names(pltSF)])
      } else {
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db_clip$PLOT) | grpBy %in% names(db_clip$COND)])
      }

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
            library(tidyr)
          })
          bOut <- parLapply(cl, X = names(combos), fun = biomassHelper, combos, data, grpBy, aGrpBy, totals, SE)
          stopCluster(cl)
        } else { # Unix systems
          bOut <- mclapply(X = names(combos), FUN = biomassHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
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
    out[[y]] <- bOut
    pb$tick()
  }
  bOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  bOut <- drop_na(bOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR) %>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) bOut <- st_sf(bOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) bOut <- unique(bOut)
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
                nCores = 1) {

  ## Need a plotCN
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        inner_join(db$TREE, by = 'PLT_CN') %>%
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
  }

  db$COND_DWM_CALC <- db[['COND_DWM_CALC']] %>% mutate(DWM_CN = CN)
  db$COND <- db[['COND']] %>% mutate(CND_CN = CN)

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

  # Save original grpByfor pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy

  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf') %>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    grpBy <- c(grpBy, 'LON', 'LAT')
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

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPDWM') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)
    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))
  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPDWM') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', grpP, 'sp', 'aD_p')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CND_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'aD_c', 'landD', 'CONDID', grpC)), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      #left_join(select(db_clip$COND_DWM_CALC, -c( 'STATECD', 'COUNTYCD', 'UNITCD', 'INVYR', 'MEASYEAR', 'PLOT', 'CONDID', 'EVALID', 'STRATUM_CN')), by = c('PLT_CN', 'CND_CN')) %>%
      left_join(select(db_clip$COND_DWM_CALC, -c( 'STATECD', 'COUNTYCD', 'UNITCD', 'INVYR', 'MEASYEAR', 'PLOT', 'EVALID')), by = c('PLT_CN', 'CONDID', 'STRATUM_CN')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      #mutate(tAdj = adjHelper(DIA, MACRO_BREAKPOINT_DIA, ADJ_FACTOR_MICR, ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM) %>%
      mutate_if(is.factor,
                as.character) %>%
      filter(EVAL_TYP == 'EXPDWM') %>%
      distinct(PLT_CN, CONDID, EVALID, COND_STATUS_CD, .keep_all = TRUE)

    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }


    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp


    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {
      cOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
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
          BIO_DUFF = sum(DUFF_BIOMASS* aDI / 2000, na.rm = TRUE),
          BIO_LITTER = sum(LITTER_BIOMASS * aDI / 2000, na.rm = TRUE),
          BIO_1HR = sum(FWD_SM_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
          BIO_10HR = sum(FWD_MD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
          BIO_100HR = sum(FWD_LG_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
          BIO_1000HR = sum(CWD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
          BIO_PILE = sum(PILE_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
          BIO = sum(BIO_LITTER, BIO_DUFF, BIO_1HR, BIO_10HR, BIO_100HR, BIO_1000HR, BIO_PILE, na.rm = TRUE),
          CARB_DUFF = sum(DUFF_CARBON* aDI / 2000, na.rm = TRUE),
          CARB_LITTER = sum(LITTER_CARBON * aDI / 2000, na.rm = TRUE),
          CARB_1HR = sum(FWD_SM_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
          CARB_10HR = sum(FWD_MD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
          CARB_100HR = sum(FWD_LG_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
          CARB_1000HR = sum(CWD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
          CARB_PILE = sum(PILE_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
          CARB = sum(CARB_LITTER, CARB_DUFF, CARB_1HR, CARB_10HR, CARB_100HR, CARB_1000HR, CARB_PILE, na.rm = TRUE))

      if (returnSpatial){
        cOut <- cOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN TPA -- Total number of trees in region & Mean TPA for the region
    } else {
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
            library(tidyr)
          })
          cOut <- parLapply(cl, X = names(combos), fun = dwmHelper, combos, data, grpBy, totals, SE)
          stopCluster(cl)
        } else { # Unix systems
          cOut <- mclapply(names(combos), FUN = dwmHelper, combos, data, grpBy, totals, SE, mc.cores = nCores)
        }
      })

      if (SE){
        # Convert from list to dataframe
        cOut <- do.call(rbind,cOut)  %>%
          as.data.frame()

        ## IF the user wants a tidy dataframe at the end, handle it for them
        if (tidy){
          # Gather up columns
          vol <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_ACRE', VOL_DUFF_ACRE:VOL_PILE_ACRE)
          bio <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_ACRE', BIO_DUFF_ACRE:BIO_PILE_ACRE)
          carb <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_ACRE', CARB_DUFF_ACRE:CARB_PILE_ACRE)
          volSE <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_ACRE_SE', VOL_DUFF_ACRE_SE:VOL_PILE_ACRE_SE)
          bioSE <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_ACRE_SE', BIO_DUFF_ACRE_SE:BIO_PILE_ACRE_SE)
          carbSE <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_ACRE_SE', CARB_DUFF_ACRE_SE:CARB_PILE_ACRE_SE)
          # Join them back up all nice like
          cTidy <- bind_cols(select(vol, c(names(combos[[1]]), 'FUEL_TYPE', 'VOL_ACRE', 'nPlots')),
                             select(bio, BIO_ACRE),
                             select(carb, CARB_ACRE),
                             select(volSE, VOL_ACRE_SE),
                             select(bioSE, BIO_ACRE_SE),
                             select(carbSE, CARB_ACRE_SE))
          if(totals){
            volT <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_TOTAL', VOL_DUFF:VOL_PILE)
            bioT <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_TOTAL', BIO_DUFF:BIO_PILE)
            carbT <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_TOTAL', CARB_DUFF:CARB_PILE)
            volTSE <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_TOTAL_SE', VOL_DUFF_SE:VOL_PILE_SE)
            bioTSE <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_TOTAL_SE', BIO_DUFF_SE:BIO_PILE_SE)
            carbTSE <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_TOTAL_SE', CARB_DUFF_SE:CARB_PILE_SE)
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
          # Gather up columns
          vol <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_ACRE', VOL_DUFF_ACRE:VOL_PILE_ACRE)
          bio <- gather(cOut, key = 'FUEL_TYPE', value = 'BIO_ACRE', BIO_DUFF_ACRE:BIO_PILE_ACRE)
          carb <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_ACRE', CARB_DUFF_ACRE:CARB_PILE_ACRE)
          # Join them back up all nice like
          cTidy <- bind_cols(select(vol, c(grpBy, 'FUEL_TYPE', 'VOL_ACRE', 'nPlots')),
                             select(bio, BIO_ACRE),
                             select(carb, CARB_ACRE))


          if(totals){
            volT <- gather(cOut, key = 'FUEL_TYPE', value = 'VOL_TOTAL', VOL_DUFF:VOL_PILE)
            bioT <- gather(cOut, key = 'BIO_TYPE', value = 'BIO_TOTAL', BIO_DUFF:BIO_PILE)
            carbT <- gather(cOut, key = 'FUEL_TYPE', value = 'CARB_TOTAL', CARB_DUFF:CARB_PILE)
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
                           "BIO_DUFF_ACRE" = rep(NA, nrow(combos)), "BIO_DUFF_ACRE_SE" = rep(NA, nrow(combos)),
                           "BIO_LITTER_ACRE" = rep(NA, nrow(combos)), "BIO_LITTER_ACRE_SE" = rep(NA, nrow(combos)),
                           "BIO_1HR_ACRE" = rep(NA, nrow(combos)), "BIO_1HR_ACRE_SE" = rep(NA, nrow(combos)),
                           "BIO_10HR_ACRE" = rep(NA, nrow(combos)), "BIO_10HR_ACRE_SE" = rep(NA, nrow(combos)),
                           "BIO_100HR_ACRE" = rep(NA, nrow(combos)), "BIO_100HR_ACRE_SE" = rep(NA, nrow(combos)),
                           "BIO_1000HR_ACRE" = rep(NA, nrow(combos)), "BIO_1000HR_ACRE_SE" = rep(NA, nrow(combos)),
                           "BIO_PILE_ACRE" = rep(NA, nrow(combos)), "BIO_PILE_ACRE_SE" = rep(NA, nrow(combos)),
                           "BIO_ACRE" = rep(NA, nrow(combos)), "BIO_ACRE_SE" = rep(NA, nrow(combos)),
                           "CARB_DUFF_ACRE" = rep(NA, nrow(combos)), "CARB_DUFF_ACRE_SE" = rep(NA, nrow(combos)),
                           "CARB_LITTER_ACRE" = rep(NA, nrow(combos)), "CARB_LITTER_ACRE_SE" = rep(NA, nrow(combos)),
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
    out[[y]] <- cOut
    pb$tick()
  }
  cOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  cOut <- drop_na(cOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR)%>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) cOut <- st_sf(cOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) {
    cOut <- unique(cOut)
  } else {
    cOut <- filter(cOut, nPlots > 0)
  }
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
                     nCores = 1){
  # Need a PLT_CN
  db$PLOT <- db[["PLOT"]] %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        0
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
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
  }

  reqTables <- c('PLOT', 'INVASIVE_SUBPLOT_SPP', 'COND',
                 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
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

  ## Link up by species is required, cannot accurately estimate coverage of all species
  grpBy <- c("YEAR", grpBy, 'SYMBOL', 'SCIENTIFIC_NAME', 'COMMON_NAME')
  grpByOrig <- grpBy


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS
  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    grpBy <- c(grpBy, 'LON', 'LAT')
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

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)
    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))
  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'INVASIVE_SAMPLING_STATUS_CD', grpP, 'sp', 'aD_p')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'aD_c', 'landD', 'CONDID', grpC)), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      full_join(select(db_clip$INVASIVE_SUBPLOT_SPP, c('PLT_CN', 'COVER_PCT', 'VEG_SPCD', 'SUBP', 'CONDID')), by = c("PLT_CN", "CONDID"))
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

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }


    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp


    if (byPlot) {
      invOut <- data %>%
        mutate(YEAR = INVYR) %>%
        filter(INVASIVE_SAMPLING_STATUS_CD == 1) %>%
        distinct(PLT_CN, CONDID, SUBP, VEG_SPCD, .keep_all = TRUE) %>%
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(cover = sum(COVER_PCT/100 * CONDPROP_UNADJ * aDI * 24^2*pi, na.rm = TRUE)) %>%
        filter(!is.na(SYMBOL))

      if (returnSpatial){
        invOut <- invOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }
    } else {
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
      aGrpBy <- grpBy[grpBy %in% c('SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE]

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
            library(tidyr)
          })
          invOut <- parLapply(cl, X = names(combos), fun = invasiveHelper, combos, data, grpBy, aGrpBy, totals, SE)
          stopCluster(cl)
        } else { # Unix systems
          invOut <- mclapply(X = names(combos), FUN = invasiveHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
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
    out[[y]] <- invOut
    pb$tick()
  }
  invOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  invOut <- drop_na(invOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR)%>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) invOut <- st_sf(invOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) invOut <- unique(invOut)
  return(invOut)

}



#' @export
area <- function(db,
                 grpBy = NULL,
                 polys = NULL,
                 returnSpatial = FALSE,
                 byLandType = FALSE,
                 landType = 'forest',
                 treeDomain = NULL,
                 areaDomain = NULL,
                 totals = FALSE,
                 byPlot = FALSE,
                 SE = TRUE,
                 nCores = 1) {

  ## Need a plotCN
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        inner_join(db$TREE, by = 'PLT_CN') %>%
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
  if (landType %in% c('timber', 'forest', 'non-forest', 'census water', 'non-census water', 'water', 'all') == FALSE){
    stop('landType must be one of: "forest", "non-forest", "census water", "non-census water", "water", "all".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }

  # Save original grpBy for pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
        inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID', 'STATECD')), by = 'PLT_CN') %>%
        inner_join(select(db$POP_EVAL, c('EVALID', 'END_INVYR')), by = 'EVALID') %>%
        group_by(polyID) %>%
        summarize(maxYear = max(END_INVYR, na.rm = TRUE))
    }

    ## TO RETURN SPATIAL PLOTS

  } else if (byPlot & returnSpatial){
    ## Make plot data spatial, projected same as polygon layer
    coordinates(db$PLOT) <- ~LON+LAT
    proj4string(db$PLOT) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    db$PLOT <- as(db$PLOT, 'sf')
  } # END AREAL


  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if(byLandType == FALSE){
    if (tolower(landType) == 'forest'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
    } else if (tolower(landType) == 'timber'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
    } else if (tolower(landType) == 'non-forest'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 2)
    } else if (tolower(landType) == 'water'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 3 | db$COND$COND_STATUS_CD == 4)
    } else if (tolower(landType) == 'census water'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 4)
    } else if (tolower(landType) == 'non-census water'){
      db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 4)
    } else if (tolower(landType) == 'all') {
      db$COND$landD <- 1
    }
  } else {
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

  # Same as above for tree (ex. trees > 20 ft tall)
  treeDomain <- substitute(treeDomain)
  tD <- eval(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$TREE$tD <- as.numeric(tD)

  # Make a new column that describes the land type and hold in COND
  if (byLandType){
    grpBy <- c(grpBy, 'landType')
    grpByOrig <- c(grpByOrig, 'landType')
    db$COND <- db$COND %>%
      mutate(landType = case_when(
        COND_STATUS_CD == 1 & SITECLCD %in% c(1:6) & RESERVCD ==0 ~ 'Timber',
        COND_STATUS_CD == 1 ~ 'Non-Timber Forest',
        COND_STATUS_CD == 2 ~ 'Non-Forest',
        COND_STATUS_CD == 3 | COND_STATUS_CD == 4 ~ 'Water'))
    db$COND <- db$COND[!is.na(db$COND$landType),]
  }

  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)
    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))
  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', grpP, 'aD_p', 'sp')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM) %>%
      mutate_if(is.factor,
                as.character) %>%
      filter(!is.na(YEAR))

    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    # Test if any polygons cross state boundaries w/ different recent inventory years
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1 & !is.null(polys)){
      # Replace YEAR from above w/ max year so that data is pooled across states
      data <- left_join(data, mergeYears, by = 'polyID') %>%
        select(-c(YEAR)) %>%
        mutate(YEAR = maxYear)
    }

    # If we need to, join the tree table
    if (length(grpT) > 0 | !is.null(treeDomain)){
      data <- data %>%
        left_join(select(db_clip$TREE, c('PLT_CN', 'CONDID', 'SUBP', 'TREE', grpT, 'tD')), by = c('PLT_CN', 'CONDID'))
    } else {
      data$tD <- 1
    }

    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp * data$tD




    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {
      aOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
        group_by(.dots = grpBy, PLT_CN, CONDID) %>%
        summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
                  aAdj = first(aAdj),
                  aDI = first(aDI)) %>%
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(area = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                  plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

      if (returnSpatial){
        invOut <- invOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN
    } else {
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
            library(tidyr)
          })
          aOut <- parLapply(cl, X = names(combos), fun = areaHelper, combos, data, grpBy, totals, SE)
          stopCluster(cl)
        } else { # Unix systems
          aOut <- mclapply(X = names(combos), FUN = areaHelper, combos, data, grpBy, totals, SE, mc.cores = nCores)
        }
      })

      if (SE){
        # Convert from list to dataframe
        aOut <- do.call(rbind,aOut)
      } else {
        # Pull out dataframe
        aOut <- aOut[[1]]
      }

      # Snag some names for below
      aNames <- names(aOut)[names(aOut) %in% grpBy == FALSE]

      # Return a spatial object
      if ('YEAR' %in% names(aOut)){
        # Return a spatial object
        if (!is.null(polys) & returnSpatial) {
          suppressMessages({suppressWarnings({aOut <- left_join(polys, aOut) %>%
            select(c(grpByOrig, aNames, names(polys))) %>%
            filter(!is.na(polyID))})})
        } else if (!is.null(polys) & returnSpatial == FALSE){
          aOut <- select(aOut, c(grpByOrig, aNames, everything())) %>%
            filter(!is.na(polyID))
        }
      } else { ## Function found no plots within the polygon, so it panics
        combos <- data %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize()
        aOut <- data.frame("YEAR" = combos$YEAR, "AREA" = rep(NA, nrow(combos)),
                           "AREA_SE" = rep(NA, nrow(combos)), "nPlots" = rep(NA, nrow(combos)))
        if (!is.null(polys) & returnSpatial) {
          suppressMessages({suppressWarnings({
            polys = left_join(polys, combos)
            aOut <- left_join(polys, aOut)})})
        } else if (!is.null(polys) & returnSpatial == FALSE){
          aOut <- select(aOut, c(grpByOrig, everything()))
        }
      }
    } # End byPlot == FALSE
    out[[y]] <- aOut
    pb$tick()
  }

  aOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  ## Remove NA values from groups
  aOut <- drop_na(aOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR)%>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) aOut <- st_sf(aOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) aOut <- unique(aOut)
  return(aOut)
}



## Seedling TPA
#' @export
seedling <- function(db,
                     grpBy = NULL,
                     polys = NULL,
                     returnSpatial = FALSE,
                     bySpecies = FALSE,
                     landType = 'forest',
                     treeDomain = NULL,
                     areaDomain = NULL,
                     totals = FALSE,
                     byPlot = FALSE,
                     SE = TRUE,
                     nCores = 1) {

  ## Need a plotCN
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN)

  ## Converting names given in grpBy to character vector (NSE to standard)
  ##  don't have to change original code
  grpBy_quo <- enquo(grpBy)

  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[1,] %>% # Just the first row
        inner_join(db$COND, by = 'PLT_CN') %>%
        inner_join(db$SEEDLING, by = 'PLT_CN') %>%
        select(!!grpBy_quo)
    )

    # If column doesnt exist, just returns 0, not a dataframe
    if (is.null(nrow(d_quo))){
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT, SEEDLING, or COND tables. Did you accidentally quote the variables names? e.g. use grpBy = ECOSUBCD (correct) instead of grpBy = "ECOSUBCD". ', collapse = ', '))
    } else {
      # Convert to character
      grpBy <- names(d_quo)
    }
  }

  reqTables <- c('PLOT', 'SEEDLING', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
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
  db$PLOT <- left_join(db$PLOT, select(db$SURVEY, STATECD, STATENM, ANN_INVENTORY, INVYR, RSCD), by = c("INVYR", "STATECD"))

  #22 and 23 is good to go for all manual years

  if(any(db$PLOT$RSCD %in% c(22,23) == FALSE)){
    if (any(db$PLOT$MANUAL < 2)) {
      ## Grab the states that are whacked out
      stYear <- db$PLOT %>%
        filter(RSCD %in% c(22,23) == FALSE & MANUAL < 2 & str_to_upper(ANN_INVENTORY) == 'Y') %>%
        distinct(RSCD, MANUAL, .keep_all = TRUE)

      warning(paste0('Seedling counts recorded as the actual seedling count up to six seedlings and then record 6+ if at least six seedlings were present in the following states/years: ', paste(as.character(stYear$STATENM), as.character(stYear$INVYR), collapse = ', '), '. Abundance is likely underestimated in these years.'))
    }
    }

  # Save original grpBy for pretty return with spatial objects
  grpBy <- c('YEAR', grpBy)
  grpByOrig <- grpBy


  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf')%>%
      mutate_if(is.factor,
                as.character)
    # Add shapefile names to grpBy
    grpBy = c(names(polys)[str_detect(names(polys), 'geometry') == FALSE], 'polyID', grpBy)
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
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(pltSF, by = 'PLT_CN')


    # Test if any polygons cross state boundaries w/ different recent inventory years (continued w/in loop)
    if ('mostRecent' %in% names(db) & length(unique(db$POP_EVAL$STATECD)) > 1){
      mergeYears <- pltSF %>%
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

  # Same as above for SEEDLING (ex. SEEDLINGs > 20 ft tall)
  treeDomain <- substitute(treeDomain)
  tD <- eval(treeDomain, db$SEEDLING) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$SEEDLING$tD <- as.numeric(tD)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$SEEDLING)[names(db$SEEDLING) %in% grpBy]

  ### Snag the EVALIDs that are needed
  ## To speed up processing time we will loop over reporting years and polygons and use clipFIA to reduce the number of rows of data
  ## Joining is quick, so we just do that on each core. Grouping and summarizing is slow with high n, so we focus on reducing n
  if (!is.null(polys)){
    ids <- pltSF %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, 'PLT_CN', 'EVALID'), by = 'PLT_CN') %>%
      left_join(select(db$POP_EVAL, 'CN', 'END_INVYR', 'EVALID'), by = 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(polyID, END_INVYR, EVALID)

    ## Must be a most recent subset for mergeYears to exists
    ## If so, we want to combine EVALIDs for poygons which straddle state boundaries
    ## w/ different most recent inventory years (ex. VA - 2016, NC - 2017)
    if (exists('mergeYears')){
      ids <- ids %>%
        left_join(mergeYears, by = 'polyID') %>%
        mutate(END_INVYR = maxYear)
    }
    ## Snag all the EVALIDs for each poly
    ids <- ids %>%
      group_by(polyID, END_INVYR) %>%
      summarise(id = list(EVALID))


  } else {
    ids <- db$POP_EVAL %>%
      select('CN', 'END_INVYR', 'EVALID') %>%
      inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
      filter(EVAL_TYP == 'EXPVOL' | EVAL_TYP == 'EXPCURR') %>%
      filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
      distinct(END_INVYR, EVALID) %>%
      group_by(END_INVYR) %>%
      summarise(id = list(EVALID))
  }

  ## Add a progress bar with ETA
  pb <- progress_bar$new(
    format = "  Computing Estimates [:bar] :percent  \t ETA: :eta \t Elapsed: :elapsed",
    total = nrow(ids), clear = FALSE, width= 100)

  # Looping over years (NOT PARALLEL, parallelization is applied to the groups to prevent spreading the entire db across cores)
  out <- list()
  for (y in 1:nrow(ids)){
    ## Clip out the necessary data
    if (!is.null(polys)){
      ## We do the spatial and temporal clip seperately because we can reuse the spatially clipped object
      ## in future temporal clips, assuming multiple reporting years exists within the polygons (not most recent clip)
      ## This prevents us from re-running st_intersection for each year - poly combo when poly is the same as previous.
      if (y == 1){
        ## First iteration, do both spatial and temporal
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else if (polys$polyID[polys$polyID == ids$polyID[y]] == polys$polyID[polys$polyID == ids$polyID[y-1]]) {
        ## Just rerun temporal clip with the same spatially clipped object
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      } else {
        ## Hit a new poly, make a new spatial clip
        db_poly <- clipFIA(db, mostRecent = FALSE, mask = polys[polys$polyID == ids$polyID[y],])
        db_clip <- clipFIA(db_poly, mostRecent = FALSE, evalid = ids$id[[y]])
      }
      # update spatial domain indicator
      db_clip$PLOT$sp <- ifelse(db_clip$PLOT$PLT_CN %in% pltSF$PLT_CN, 1, 0)

      ## Polys not specified, just temporal
    } else {
      db_clip <- clipFIA(db, mostRecent = FALSE, evalid = ids$id[[y]])
      # update spatial domain indicator
      db_clip$PLOT$sp <- 1
    }

    ## Prep joins and filters
    data <- select(db_clip$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', grpP, 'aD_p', 'sp')) %>%
      left_join(select(db_clip$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN')), by = c('PLT_CN')) %>%
      left_join(select(db_clip$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MACR', 'CN', 'P1POINTCNT')), by = c('STRATUM_CN' = 'CN')) %>%
      left_join(select(db_clip$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('ESTN_UNIT_CN' = 'CN')) %>%
      right_join(select(db_clip$POP_EVAL, c('EVALID', 'EVAL_GRP_CN', 'ESTN_METHOD', 'CN', 'END_INVYR', 'REPORT_YEAR_NM')), by = c('EVAL_CN' = 'CN')) %>%
      # left_join(select(db_clip$POP_EVAL_TYP, c('EVAL_TYP', 'EVAL_CN')), by = c('EVAL_CN')) %>%
      # left_join(select(db_clip$POP_EVAL_GRP, c('RSCD', 'CN', 'EVAL_GRP')), by = c('EVAL_GRP_CN' = 'CN')) %>%
      left_join(select(db_clip$SEEDLING, c('PLT_CN', 'CONDID', 'SPCD', 'TPA_UNADJ', 'SUBP', grpT, 'tD')), by = c('PLT_CN', 'CONDID')) %>%
      mutate(aAdj = ifelse(PROP_BASIS == 'SUBP', ADJ_FACTOR_SUBP, ADJ_FACTOR_MACR)) %>%
      mutate(tAdj = ADJ_FACTOR_MICR) %>%
      rename(YEAR = END_INVYR,
             YEAR_RANGE = REPORT_YEAR_NM) %>%
      mutate_if(is.factor,
                as.character) %>%
      filter(!is.na(YEAR)) %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, EVALID, COND_STATUS_CD, .keep_all = TRUE)

    ## Recode a few of the estimation methods to make things easier below
    data$ESTN_METHOD = recode(.x = data$ESTN_METHOD,
                              `Post-Stratification` = 'strat',
                              `Stratified random sampling` = 'strat',
                              `Double sampling for stratification` = 'double',
                              `Simple random sampling` = 'simple',
                              `Subsampling units of unequal size` = 'simple')

    ## Comprehensive indicator function
    data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
    data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp
    data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp

    ## Add species to groups
    if (bySpecies) {
      data <- data %>%
        left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
        mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' ')) %>%
        mutate_if(is.factor,
                  as.character)
      grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
      grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    }


    ####################  COMPUTE ESTIMATES  ###########################
    ### -- BYPLOT -- TPA Estimates at each plot location
    if (byPlot) {

      #grpBy_sym <- syms(grpBy)
      tOut <- data %>%
        mutate(YEAR = INVYR) %>%
        distinct(PLT_CN, SUBP, .keep_all = TRUE) %>%
        group_by(.dots = grpBy, PLT_CN) %>%
        summarize(TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
                  TPA_PERC = TPA / sum(TPA_UNADJ * pDI, na.rm = TRUE) * 100)

      # tOut <- data %>%
      #   #lazy_dt() %>%
      #   #filter(EVAL_TYP == 'EXPVOL') %>%
      #   distinct(PLT_CN, CONDID, SUBP, SEEDLING, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      #   group_by(!!grpBy_sym, PLT_CN) %>%
      #   summarize(TPA = sum(TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
      #             BAA = sum(basalArea(DIA) * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
      #             nStems = length(which(tDI == 1))) #%>%
      #   #as_tibble()



      if (returnSpatial){
        tOut <- tOut %>%
          filter(!is.na(LAT) & !is.na(LON)) %>%
          st_as_sf(coords = c('LON', 'LAT'),
                   crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

      }

      ### -- TOTALS & MEAN TPA -- Total number of SEEDLINGs in region & Mean TPA for the region
    } else {
      ## Duplicate and rbind so we have a unique poly key for the entire set of non-spatial combos
      if(is.null(polys)){
        combos <- select(data, c(grpBy)) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize() %>%
          filter(!is.na(YEAR))
      } else {
        ## Non spatial combos
        combosNS <- select(data, c(grpBy[grpBy %in% names(polys) == FALSE])) %>%
          as.data.frame() %>%
          group_by(.dots = grpBy[grpBy %in% names(polys) == FALSE]) %>%
          summarize() %>%
          filter(!is.na(YEAR))
        combosNSpoly <- combosNS %>%
          mutate(polyID = polys$polyID[polys$polyID == ids$polyID[y]])

        # New combos for spatial objects
        combos <- pltSF %>%
          select(-c(PLT_CN)) %>%
          distinct(polyID, .keep_all = TRUE) %>%
          right_join(combosNSpoly, by = 'polyID')
      }
      # List of rows for lapply
      combos <- split(combos, seq(nrow(combos)))

      # Seperate area grouping names, (ex. TPA red oak in total land area of ingham county, rather than only area where red oak occurs)
      if (!is.null(polys)){
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND) | grpBy %in% names(pltSF)])
      } else {
        aGrpBy <- c('YEAR', grpBy[grpBy %in% names(db$PLOT) | grpBy %in% names(db$COND)])
      }

      #message('Computing Summary Statistics.....')
      #pb <- progress_bar$new(format = "[:bar] :percent eta: :eta", total = nrow(combos), clear = FALSE, width= 100)

      #for (i in 1:nrow(combos)){ #, .combine = 'rbind', .packages = 'dplyr', .export = c('data')
      #tOut <- foreach(i = 1:nrow(combos), .combine = 'rbind', .packages = 'dplyr', .export = c('data'))

      suppressWarnings({
        ## Compute estimates in parallel -- Clusters in windows, forking otherwise
        if (Sys.info()['sysname'] == 'Windows'){
          cl <- makeCluster(nCores)
          clusterEvalQ(cl, {
            library(dplyr)
            library(stringr)
          })
          tOut <- parLapply(cl, X = names(combos), fun = seedHelper, combos, data, grpBy, aGrpBy, totals, SE)
          stopCluster(cl)
        } else { # Unix systems
          tOut <- mclapply(X = names(combos), FUN = seedHelper, combos, data, grpBy, aGrpBy, totals, SE, mc.cores = nCores)
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
          # Return spatial plots
        }
      } else { ## Function found no plots within the polygon, so it panics
        combos <- data %>%
          as.data.frame() %>%
          group_by(.dots = grpBy) %>%
          summarize()
        tOut <- data.frame("YEAR" = combos$YEAR,
                           "TPA" = rep(NA, nrow(combos)),
                           "TPA_PERC" = rep(NA, nrow(combos)),
                           "TPA_SE" = rep(NA, nrow(combos)),
                           "TPA_PERC_SE" = rep(NA, nrow(combos)),
                           "nPlots_SEEDLING" = rep(NA, nrow(combos)),
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
    out[[y]] <- tOut
    pb$tick()
  }

  tOut <- do.call(rbind, out)
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  tOut <- drop_na(tOut, grpBy[grpBy %in% names(polys) == FALSE]) %>%
    arrange(YEAR) %>%
    as_tibble()
  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
  return(tOut)
}



