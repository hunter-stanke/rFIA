
#### A Start up Message ------------
.onAttach <- function(libname, pkgname){
  #packageStartupMessage('Download FIA Data Here: https://apps.fs.usda.gov/fia/datamart/datamart.html')
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
                (sum((c(x, numeric(ndif))^2), na.rm = TRUE) - nh * xStrat^2) / (nh * (nh-1)))
    ## Covariance
  } else {
    v <- ifelse(first(ESTN_METHOD == 'simple'),
                cov(x,y),
                (sum(x*y, na.rm = TRUE) - sum(nh * xStrat *yStrat, na.rm = TRUE)) / (nh * (nh-1)))
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

#' @import dplyr
#' @import methods
#' @import sf
#' @import stringr
#' @import gganimate
#' @import ggplot2
#' @import bit64
#' @import progress
#' @import tidyselect
#' @import purrr
#' @importFrom rlang eval_tidy enquo enquos quos quo
#' @importFrom data.table fread fwrite rbindlist
#' @importFrom parallel makeCluster detectCores mclapply parLapply stopCluster clusterEvalQ
#' @import tidyr
#' @importFrom sp over proj4string<- coordinates<- spTransform proj4string
#' @importFrom stats cov var coef lm na.omit quantile
#' @importFrom utils object.size read.csv tail globalVariables type.convert download.file unzip
#' @import dtplyr
#' @importFrom R2jags autojags jags
#' @importFrom coda as.mcmc
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
                    states = NULL,
                    inMemory = TRUE,
                    nCores = 1,
                    ...){

  ## methods for reading the full database into memory
    # Add a slash to end of directory name if missing
    if (str_sub(dir,-1) != '/') dir <- paste(dir, '/', sep = "")
    # Grab all the file names in directory
    files <- list.files(dir)

  if (inMemory){

      inTables <- list()

    # Some warnings
    if(!dir.exists(dir)) {
      stop(paste('Directory', dir, 'does not exist.'))
    }
    if(length(files[str_sub(files,-4, -1) == '.csv']) < 1){
      stop(paste('Directory', dir, 'contains no .csv files.'))
    }

    ## Some warnings up front
    ## Do not try to merge ENTIRE with other states
    if (length(states) > 1 & any(str_detect(str_to_upper(states), 'ENTIRE'))){
      stop('Cannot merge ENITRE with other state tables. ENTIRE includes all state tables combined. Do you only need data for a particular region?')
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
                  'SEEDLING', 'SURVEY', 'SUBP_COND', 'P2VEG_SUBP_STRUCTURE')
      if (any(str_sub(files, 3, 3) == '_')){
        files <- files[str_sub(files,4,-5) %in% cFiles]
      } else{
        files <- files[str_sub(files,1,-5) %in% cFiles]
      }
    }


    ## If individual tables are specified, then just grab those .csvs, otherwise download the .zip file, extract and read with fread. Should be quite a bit quicker.
    if (!is.null(states)){
      ### I'm not very smart and like specify the name twice sometimes,
      ### --> making the function smarter than me
      states <- str_to_upper(unique(states))

      ## Check to make sure states exist
      allStates <- unique(str_to_upper(str_sub(files, 1, 2)))

      if (any(states %in% allStates == FALSE)){
        missStates <- states[states %in% allStates == FALSE]
        stop(paste('Data unavailable for: ', paste(as.character(missStates),collapse = ', '), '. States not found in specified directory.'))
      }

      files <- files[str_to_upper(str_sub(files, 1, 2)) %in% states]


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

    ### Methods for keeping data remote until they are needed
    ### Chunking up data into states
    ## inMemory = FALSE
  } else {

    ## IF states isn't given, default to all
    ## states in the directory
    if (is.null(states)){
      states <- unique(str_to_upper(str_sub(files, 1, 3)))
      states <- states[str_sub(states, 3,3) == '_']
      states <- str_sub(states, 1, 2)
      ## Don't fail if states have been merged
      if (length(states) < 1) states <- 1
    }


    ## Saving the call to readFIA, for eval later
    out <- list(dir = dir,
                common = common,
                tables = tables,
                states = states,
                ... = ...)

    class(out) <- 'Remote.FIA.Database'
    return(out)
  }

}




## Access FIA Database files from the FIA Datamart
#' @export
getFIA <- function(states,
                   dir = NULL,
                   common = TRUE,
                   tables = NULL,
                   load = TRUE,
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

    message(paste0('Saving to ', dir, '. NOTE: modifying FIA tables in Excel may corrupt csv files.'))
  }

  if (is.null(dir) & load == FALSE){
    stop('Must specify a directory ("dir") to save data when "load = FALSE".')
  }

  ## If dir is not specified, hold in a temporary directory
  #if (is.null(dir)){tempDir <- tempdir()}
  tempDir <- tempdir()
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

  ### I'm not very smart and like specify the name twice sometimes,
  ### --> making the function smarter than me
  states <- unique(states)

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
        ## Unlinking the directory is bad news, so just delete the file from tempdir
        file.remove(paste0(tempDir, '/', newName))
      } else {
        #download.file(urls[n], paste0(dir, tblNames[n]))
        unzip(temp, exdir = str_sub(dir, 1, -2))
        if (load) file <- fread(paste0(dir, newName), showProgress = FALSE, logical01 = FALSE, integer64 = 'double', nThread = nCores)
      }

      unlink(temp)

      if (load){
        # We don't want data.table formats
        file <- as.data.frame(file)
        inTables[[str_sub(urls[n], 43, -5)]] <- file
      }

    }

    if (load){
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
      class(outTables) <- 'FIA.Database'
    }



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
      #newName <- paste0(str_sub(url, 1, -4), 'csv')
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

    if (load){
      ## Read in the files w/ readFIA
      if (is.null(dir)){
        outTables <- readFIA(tempDir, nCores = nCores, common = common, states = states)
        #unlink(tempDir, recursive = TRUE)
      } else {
        outTables <- readFIA(dir, nCores = nCores, common = common, states = states)
      }
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

  if (load) return(outTables)
}



## Write out the raw FIA files
#' @export
writeFIA <- function(db,
                     dir,
                     byState = FALSE,
                     nCores = 1,
                     ...){
  if (class(db) == 'Remote.FIA.Database'){
    stop('Cannot write remote database.')
  }

  if (byState & !c('SURVEY' %in% names(db))){
    stop('Need survey table for state abbreviations.')
  }

  #cat(sys.call()$dir)
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
    message(paste0('Saving to ', dir, '. NOTE: modifying FIA tables in Excel may corrupt csv files.'))
  }

  ## Method to chunk up the database into states before writing it out
  if (byState){

    stateNames <- distinct(db$SURVEY, STATECD, STATEAB)
    db$PLOT <- db$PLOT %>%
      left_join(stateNames, by = 'STATECD')

    ## Unique state abbreviations
    states <- unique(db$PLOT$STATEAB)

    ## Chunk up plot
    pltList <- split(db$PLOT, as.factor(db$PLOT$STATEAB))

    ## Loop over states, do the writing
    for (s in 1:length(states)){
      db_clip <- db
      ## Overwriting plot with shortened table
      db_clip$PLOT <- pltList[[s]]
      ## Subsetting the remaining database based
      ## on plot table
      db_clip <- clipFIA(db_clip, mostRecent = FALSE)

      ## Write it all out
      tableNames <- names(db_clip)[names(db_clip) != 'mostRecent']
      tableNames <- paste(unique(db_clip$PLOT$STATEAB), tableNames, sep = '_')
      for (i in 1:length(tableNames)){
        if (is.data.frame(db_clip[[i]])){
          fwrite(x = db_clip[[i]], file = paste0(dir, tableNames[i], '.csv'), showProgress = FALSE, nThread = nCores)
        }
      }
    }


    ## Merge states together on writing
  } else {
    tableNames <- names(db)

    for (i in 1:length(tableNames)){
      if (is.data.frame(db[[i]])){
        fwrite(x = db[[i]], file = paste0(dir, tableNames[i], '.csv'), showProgress = FALSE, nThread = nCores)
      }
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

## Spatiotemporal clip of FIADB
#' @export
clipFIA <- function(db,
                    mostRecent = TRUE,
                    mask = NULL,
                    matchEval = FALSE,
                    evalid = NULL,
                    designCD = NULL,
                    nCores = 1) {

  ## Some warnings
  if (class(db) != "FIA.Database"){
    ## Add clipFIA arguments to remote database
    if (class(db) == 'Remote.FIA.Database'){

      db$mostRecent = mostRecent
      ## Lists won't accept NULL,
      ## but this is a hack
      db['mask'] = list(mask)

      db['evalid'] = list(evalid)

      db['designCD'] = list(designCD)

      return(db)
    }
    stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
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
    if (mostRecent) warning('Returning subset by EVALID only. For most recent subset, specify `evalid = NULL` and `mostRecent = TRUE`.')
    # Join appropriate tables and filter out specified EVALIDs
    tempData <- select(db$PLOT, CN, PREV_PLT_CN) %>%
      mutate(PLT_CN = CN) %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'EVALID')), by = 'PLT_CN') %>%
      #inner_join(select(POP_ESTN_UNIT, c('EVAL_GRP_CN', 'EVALID')), by = 'EVAL_GRP_CN') %>%
      filter(EVALID %in% evalid)
    # Extract plots which relate to specified EVALID (previous for change estimation)
    PPLOT <- db$PLOT[db$PLOT$CN %in% tempData$PREV_PLT_CN,]
    if (nrow(PPLOT) < 1) PPLOT <- NULL
    db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]
  }

  ## Locate the most recent EVALID and subset plots
  if (mostRecent & is.null(evalid)){
    suppressWarnings({
    tempData <- select(db$PLOT, CN, PREV_PLT_CN) %>%
      mutate(PLT_CN = CN) %>%
      left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'EVALID')), by = c('PLT_CN'))


    ## Most recent evals by
    mrids <- findEVALID(db, mostRecent = TRUE)

    tempData <- tempData %>%
      filter(EVALID %in% mrids)
    })

    # Write out evalids so that we don't have to repeat above later
    evalid <- unique(tempData$EVALID)

    # # Extract appropriate PLOTS
    # PPLOT <- db$PLOT[db$PLOT$CN %in% tempData$PREV_PLT_CN,]
    # if (nrow(PPLOT) < 1) PPLOT <- NULL
    # db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]
    # #test = db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]

    # Extract appropriate PLOTS
    PPLOT <- db$PLOT[db$PLOT$CN %in% tempData$PREV_PLT_CN,]
    if (nrow(PPLOT) < 1) {
      PPLOT <- NULL
    } else {
      ## Cutting duplicates
      PPLOT <- PPLOT %>%
        filter(!(PLT_CN %in% tempData$PLT_CN))
    }
    db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]
    #test = db$PLOT <- db$PLOT[db$PLOT$CN %in% tempData$PLT_CN,]


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


     # Convert polygons to an sf object
     mask <- mask %>%
       as('sf') %>%
       mutate(polyID = 1:nrow(mask))
     ## Make plot data spatial, projected same as polygon layer
     pltSF <- select(db$PLOT, c('LON', 'LAT', pltID)) %>%
       filter(!is.na(LAT) & !is.na(LON)) %>%
       distinct(pltID, .keep_all = TRUE)
     coordinates(pltSF) <- ~LON+LAT
     #proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
     #pltSF <- as(pltSF, 'sf') %>%
     #   st_transform(crs = st_crs(mask)$proj4string)
     pltSF <- as(pltSF, 'sf')
     st_crs(pltSF) <- 4326 #%>%
     pltSF <- st_transform(pltSF, crs = st_crs(mask))

     ## Split up polys
     polyList <- split(mask, as.factor(mask$polyID))
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

     pltSF <- bind_rows(out) %>%
       left_join(select(db$PLOT, PLT_CN, PREV_PLT_CN, pltID), by = 'pltID')

     #if(mostRecent == FALSE & is.null(evalid)) PPLOT <- filter(db$PLOT, db$PLOT$PLT_CN %in% pltSF$PREV_PLT_CN)
     if(mostRecent == FALSE & is.null(evalid)) {
       PPLOT <- NULL
     } else {
       PPLOT <- filter(PPLOT, pltID %in% pltSF$pltID)
     }
     db$PLOT <- filter(db$PLOT, db$PLOT$PLT_CN %in% pltSF$PLT_CN)

    #  ###OLD
    #  # Convert polygons to an sf object
    #  mask <- mask %>%
    #    as('sf')
    #
    #  ## Make plot data spatial, projected same as polygon layer
    #  pltSF <- select(db$PLOT, c('PLT_CN', 'LON', 'LAT'))
    #  coordinates(pltSF) <- ~LON+LAT
    #  proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    #  pltSF <- as(pltSF, 'sf') %>%
    #    st_transform(crs = st_crs(mask)$proj4string)
    #
    #  # Intersect plot with polygons
    #  mask$polyID <- 1:nrow(mask)
    #  suppressMessages({suppressWarnings({
    #    pltSF <- st_intersection(pltSF, mask) %>%
    #      as.data.frame() %>%
    #      select(-c('geometry')) # removes artifact of SF object
    #  })})
    #
    #  # Identify the estimation units to which plots within the mask belong to
    #  estUnits <- pltSF %>%
    #    inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'STRATUM_CN')), by = 'PLT_CN') %>%
    #    inner_join(select(db$POP_STRATUM, c('CN', 'ESTN_UNIT_CN')), by = c('STRATUM_CN' = 'CN')) %>%
    #    group_by(ESTN_UNIT_CN) %>%
    #    summarize()
    #
    #  # Identify all the plots which fall inside the above estimation units
    #  plts <- select(db$PLOT, PLT_CN, PREV_PLT_CN) %>%
    #    inner_join(select(db$POP_PLOT_STRATUM_ASSGN, c('PLT_CN', 'STRATUM_CN')), by = 'PLT_CN') %>%
    #    inner_join(select(db$POP_STRATUM, c('CN', 'ESTN_UNIT_CN')), by = c('STRATUM_CN' = 'CN')) %>%
    #    filter(ESTN_UNIT_CN %in% estUnits$ESTN_UNIT_CN) %>%
    #    group_by(PLT_CN, PREV_PLT_CN) %>%
    #    summarize()
    #
    #  # Clip out the above plots from the full database, will reduce size by a shit pile
    #  PPLOT <- filter(db$PLOT, db$PLOT$PLT_CN %in% plts$PREV_PLT_CN)
    #  db$PLOT <- filter(db$PLOT, db$PLOT$PLT_CN %in% plts$PLT_CN)
    #
    # ## IF ozone is specified, do a seperate intersection (PLOTs not colocated w/ veg PLOTs)
    # if (!is.null(db$OZONE_PLOT)) {
    #   # Seperate spatial object
    #   ozoneSP <- db$OZONE_PLOT
    #
    #   ## Make PLOT data spatial, projected same as mask layer
    #   coordinates(ozoneSP) <- ~LON+LAT
    #   proj4string(ozoneSP) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    #   ozoneSP <- spTransform(ozoneSP, CRSobj = proj4string(mask))
    #
    #   ## Spatial query (keep only PLOTs that fall within shell)
    #   db$OZONE_PLOT <- ozoneSP[mask,]
    #
    #   ## Add coordinates back in to dataframe
    #   coords <- st_coordinates(db$OZONE_PLOT)
    #   db$OZONE_PLOT <- db$OZONE_PLOT %>%
    #     data.frame() %>%
    #     mutate(LAT = coords[,2]) %>%
    #     mutate(LON = coords[,1])
    # }
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
        #if(nrow(PPLOT) > 0) PPLOT$prev = 1
        if(!is.null(PPLOT) > 0) PPLOT$prev = 1
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

clipFIA_old <- function(db,
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






























