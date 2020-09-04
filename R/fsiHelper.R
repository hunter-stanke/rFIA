

fsiHelper1 <- function(x, plts, db, grpBy, scaleBy, byPlot){

  ## Does not modify outside environment, just need scaleBy in here as well
  if (is.null(grpBy)){
    aGrps <- NULL
    grpBy <- scaleBy
  } else {
    aGrps <- grpBy
    grpBy <- unique(c(grpBy, scaleBy))
  }

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]


  ## Disturbance or treatment ever happen on the plot remeasurement? If so
  ## we want to cut it before we model the max size-density curve
  disturb <- select(db$PLOT, PLT_CN, pltID) %>%
    left_join(select(db$COND, PLT_CN, DSTRBCD1, TRTCD1), by = 'PLT_CN') %>%
    mutate(DSTRBCD1 = replace_na(DSTRBCD1, 0),
           TRTCD1 = replace_na(TRTCD1, 0)) %>%
    filter(DSTRBCD1 > 0) %>%
    ## Natural regen is ok
    filter(TRTCD1 > 0 & !c(TRTCD1 %in% 40)) %>%
    ## These plots were disturbed or treated
    distinct(PLT_CN, pltID)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% c(grpBy, scaleBy)]
  grpC <- names(db$COND)[names(db$COND) %in% c(grpBy, scaleBy) & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% c(grpBy, scaleBy) & names(db$TREE) %in% c(grpP, grpC) == FALSE]

  ## Making a treeID
  db$TREE$treID <- paste(db$TREE$SUBP, db$TREE$TREE, sep = '_')

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    filter(DESIGNCD == 1 & PLOT_STATUS_CD != 3 & !is.na(REMPER) & !is.na(PREV_PLT_CN)) %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$PLOT, c('PLT_CN', 'sp', 'aD_p', 'DESIGNCD', 'PLOT_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('2', '1')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('2', '1')) %>%
    left_join(select(db$TREE, c('TRE_CN', all_of(grpT), treID, 'typeD', 'tD', 'TPA_UNADJ', 'BAA', 'DIA', 'STATUSCD', SPCD)), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%
    filter(DESIGNCD1 == 1 & PLOT_STATUS_CD1 != 3) %>%
    mutate_if(is.factor,
              as.character)

  ## Comprehensive indicator function -- w/ growth accounting
  data$tDI2 <- data$landD2 * data$aD_p2 * data$aD_c2 * data$tD2 * data$typeD2 * data$sp2 *
    if_else(data$STATUSCD2 == 1, 1, 0)

  data$tDI1 <- data$landD1 * data$aD_p1 * data$aD_c1 * data$tD1 * data$typeD1 * data$sp1 *
    if_else(data$STATUSCD1 == 1, 1, 0)

  ## Comprehensive indicator function -- w/ growth accounting
  data$pDI2 <- data$landD2 * data$aD_p2 * data$aD_c2 * data$typeD2 * data$sp2 *
    if_else(data$STATUSCD2 == 1, 1, 0)

  data$pDI1 <- data$landD1 * data$aD_p1 * data$aD_c1 * data$typeD1 * data$sp1 *
    if_else(data$STATUSCD1 == 1, 1, 0)

  ## Save a copy for area calculations
  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  aData <- select(db$PLOT, c('PLT_CN', 'PREV_PLT_CN', 'pltID', 'DESIGNCD', 'REMPER', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'MEASMON', 'MEASDAY', 'PLOT_STATUS_CD', all_of(grpP), 'aD_p', 'sp')) %>%
    filter(DESIGNCD == 1 & PLOT_STATUS_CD != 3) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    mutate(aDI = landD * aD_p * aD_c * sp)


  ## PREVIOUS and CURRENT attributes
  data <- data %>%
    mutate(TPA_UNADJ1 = TPA_UNADJ1,
           TPA_UNADJ2 = TPA_UNADJ2,
           BAA1 = BAA1,
           BAA2 = BAA2,
           MORT = case_when(
             STATUSCD1 == 1 & STATUSCD2 == 2 ~ 1,
             STATUSCD1 == 1 & STATUSCD2 == 3 ~ 1,
             TRUE ~ 0),
           SURV = case_when(
             STATUSCD1 == 1 & STATUSCD2 == 1 ~ 1,
             TRUE ~ 0)
    )


  ## Just what we need
  data <- data %>%
    select(PLT_CN, PREV_PLT_CN, pltID, TRE_CN, SUBP, CONDID, TREE, CONDPROP_UNADJ,
           MEASYEAR, MACRO_BREAKPOINT_DIA, PROP_BASIS, grpP[grpP != 'PLOT_STATUS_CD'], grpC,
           REMPER, PLOT_STATUS_CD1, PLOT_STATUS_CD2,
           treID1, treID2,
           one_of(str_c(grpT,1),str_c(grpT,2)),
           tDI1, tDI2, pDI1, pDI2, STATUSCD1, STATUSCD2,
           DIA1, DIA2, BAA1, BAA2, TPA_UNADJ1, TPA_UNADJ2, SPCD1, SPCD2) %>%
    mutate(BAA1 = -(BAA1),
           TPA_UNADJ1 = -(TPA_UNADJ1)) %>%
    ## Rearrange previous values as observations
    pivot_longer(cols = -c(PLT_CN:REMPER),
                 names_to = c(".value", 'ONEORTWO'),
                 names_sep = -1) %>%
    mutate(PLOT_BASIS = case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA) ~ NA_character_,
      ## When DIA is less than 5", use microplot value
      DIA < 5 ~ 'MICR',
      ## When DIA is greater than 5", use subplot value
      DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
      DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR'))

  ## No zeros
  data <- data %>%
    mutate(TPA_UNADJ = replace_na(TPA_UNADJ, replace = 0),
           BAA = replace_na(BAA, replace = 0))

  ## Total trees for the size-density scaling
  t1 <- data %>%
    ## No disturbance/treatment plots
    filter(!c(pltID %in% disturb$pltID)) %>%
    distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
    filter(STATUSCD == 1) %>%
    filter(pDI == 1) %>%
    group_by(.dots = scaleBy, PLT_CN) %>%
    summarize(REMPER = first(REMPER),
              BAA1 = sum(-BAA[ONEORTWO == 1], na.rm = TRUE),
              TPA1 = sum(-TPA_UNADJ[ONEORTWO == 1], na.rm = TRUE),
              BAA2 = sum(BAA[ONEORTWO == 2], na.rm = TRUE),
              TPA2 = sum(TPA_UNADJ[ONEORTWO == 2], na.rm = TRUE),
              #times1 = round(TPA_UNADJ[ONEORTWO == 1]),
              skew1 = skewness(rep(DIA[ONEORTWO == 1], round(-TPA_UNADJ[ONEORTWO == 1]))),
              skew2 = skewness(rep(DIA[ONEORTWO == 2], round(TPA_UNADJ[ONEORTWO == 2])))
              ) %>%
    ## Mean BA
    mutate(BA1 = if_else(TPA1 != 0, BAA1 / TPA1, 0),
           BA2 = if_else(TPA2 != 0, BAA2 / TPA2, 0)) %>%
    ## Remove plots with high skewness
    filter(skew2 >= -1 & skew2 <= 1)


  if (byPlot){
    grpBy <- c('YEAR', grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      select(all_of(grpBy), PLT_CN, PREV_PLT_CN, PLOT_STATUS_CD, REMPER, SUBP, TREE,
             ONEORTWO, tDI, TPA_UNADJ, BAA)
      # # Compute estimates at plot level
      # group_by(.dots = grpBy, PLT_CN, PREV_PLT_CN) %>%
      # summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
      #           REMPER = first(REMPER),
      #           PREV_BAA = sum(-BAA[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1& STATUSCD == 1], na.rm = TRUE),
      #           PREV_TPA = sum(-TPA_UNADJ[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
      #           CHNG_TPA = sum(TPA_UNADJ[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
      #           ## Sum here to avoid issues w/ zeros
      #           CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
      #           PLOT_STATUS_CD = if_else(any(PLOT_STATUS_CD == 1), 1, 2)) %>%
      # ## Replace any NAs with zeros
      # mutate(PREV_BAA = replace_na(PREV_BAA, 0),
      #        PREV_TPA = replace_na(PREV_TPA, 0),
      #        CHNG_BAA = replace_na(CHNG_BAA, 0),
      #        CHNG_TPA = replace_na(CHNG_TPA, 0),
      #        ## T2 attributes
      #        CURR_BAA = PREV_BAA + CHNG_BAA,
      #        CURR_TPA = PREV_TPA + CHNG_TPA) %>%
      # ## Change in average tree BA and QMD
      # mutate(PREV_BA = if_else(PREV_TPA != 0, PREV_BAA / PREV_TPA, 0),
      #        CURR_BA = if_else(CURR_TPA != 0, CURR_BAA / CURR_TPA, 0),
      #        CHNG_BA = CURR_BA - PREV_BA,
      #        PREV_QMD = sqrt(PREV_BA / 0.005454154),
      #        CURR_QMD = sqrt(CURR_BA / 0.005454154),
      #        CHNG_QMD = CURR_QMD - PREV_QMD) %>%
      # ## All CHNG becomes an annual rate
      # mutate(CHNG_BAA = CHNG_BAA / REMPER,
      #        CHNG_TPA = CHNG_TPA / REMPER,
      #        CHNG_BA = CHNG_BA / REMPER,
      #        CHNG_QMD = CHNG_QMD / REMPER) %>%
      # #left_join(stems, by = c('PLT_CN', grpBy[!(grpBy %in% c('pltID', 'PLOT_STATUS_CD', 'YEAR'))])) %>%
      # ### Then divide by number of unique stems for an average
      # #mutate(CHNG_BAA = CHNG_BAA / n) %>%
      # ungroup() %>%
      # select(PLT_CN, PREV_PLT_CN, PLOT_STATUS_CD, REMPER, grpBy,
      #        PREV_TPA, PREV_BAA, PREV_BA, PREV_QMD,
      #        CHNG_TPA, CHNG_BAA, CHNG_BA, CHNG_QMD,
      #        CURR_TPA, CURR_BAA, CURR_BA, CURR_QMD,
      #        nLive)

    a = NULL

  } else {

    ### Plot-level estimates
    if (length(aGrps[aGrps %in% names(aData)]) < 1) {
      aGrps = NULL
    }  else {
      aGrps <- aGrps[aGrps %in% names(aData)]
    }


    ### Plot-level estimates
    a <- aData %>%
      ## date column
      mutate(date = paste(MEASYEAR, MEASMON, MEASDAY, sep = '-'),
             date = as.Date(date, "%Y-%m-%d")) %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, pltID, date, PROP_BASIS, CONDID, all_of(aGrps), .keep_all = TRUE) %>%
      group_by(PLT_CN, pltID, date, PROP_BASIS, CONDID, .dots = aGrps) %>%
      summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ * aDI)) %>%
      mutate(CONDPROP_UNADJ = replace_na(CONDPROP_UNADJ, 0)) %>%
      ## Average forested area between min and max date
      group_by(pltID, PROP_BASIS, .dots = aGrps) %>%
      summarize(minDate = min(date, na.rm = TRUE),
                maxDate = max(date, na.rm = TRUE),
                amin = sum(CONDPROP_UNADJ[date == minDate], na.rm = TRUE),
                amax = sum(CONDPROP_UNADJ[date == maxDate], na.rm = TRUE),
                fa = (amin + amax) / 2) %>%
      left_join(select(ungroup(db$PLOT), PLT_CN, pltID), by = c('pltID'))

    ### Compute total TREES in domain of interest
    t <- data %>%
      distinct(PLT_CN, TRE_CN, ONEORTWO, .keep_all = TRUE) %>%
      select(all_of(grpBy), PLT_CN, pltID, PLOT_STATUS_CD, PLOT_BASIS, MEASYEAR, REMPER, TRE_CN,
             ONEORTWO, STATUSCD, tDI, TPA_UNADJ, BAA)
      # # Compute estimates at plot level
      # group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
      # summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
      #           REMPER = first(REMPER),
      #           MEASYEAR = first(MEASYEAR),
      #           PREV_TPA = sum(-TPA_UNADJ[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
      #           PREV_BAA = sum(-BAA[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
      #           CHNG_TPA = sum(TPA_UNADJ[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
      #           CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
      #           CURR_TPA = PREV_TPA + CHNG_TPA,
      #           CURR_BAA = PREV_BAA + CHNG_BAA,
      #           plotIn = if_else(sum(tDI, na.rm = TRUE) > 0, 1, 0))
  }

  pltOut <- list(t = t, a = a, t1 = t1)
  return(pltOut)

}






fsiHelper2 <- function(x, popState, t, a, grpBy, scaleBy, method, useSeries){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  ######## ------------------ TREE ESTIMATES + CV
  aAdj <- a %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    mutate(
      ## AREA
      aAdj = case_when(
        ## When NA, stay NA
        is.na(PROP_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = fa * aAdj) %>%
    ungroup()

  ## Sometimes less specific area groups
  aGrps <- unique(grpBy[grpBy %in% names(aAdj)])


  ## Strata level estimates
  tEst <- t %>%
    ungroup() %>%
    ## Converting to average tree size
    mutate(BA = if_else(ONEORTWO == 1, -BAA / TPA_UNADJ, BAA / TPA_UNADJ)) %>%
    ## Summing within scaleBy
    group_by(.dots = unique(c(grpBy[!c(grpBy %in% c('YEAR', 'INVYR'))], scaleBy)), PLT_CN, pltID, MEASYEAR, REMPER, PLOT_BASIS) %>%
    summarize(PREV_RD = -sum(rd[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
              CURR_RD = sum(rd[ONEORTWO == 2] * tDI[ONEORTWO == 2], na.rm = TRUE),
              PREV_TPA = -sum(TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
              PREV_BA = -sum(BA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
              CHNG_TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE) / first(REMPER),
              CHNG_BA = sum(BA * tDI, na.rm = TRUE) / first(REMPER),
              plotIn_t = if_else(sum(tDI, na.rm = TRUE) > 0, 1, 0)) %>%
    ## Summing across scaleBy
    group_by(PLT_CN, pltID, MEASYEAR, PLOT_BASIS,
             REMPER, .dots = grpBy[!c(grpBy %in% c('YEAR', 'INVYR'))]) %>%
    summarize(CHNG_TPA = sum(CHNG_TPA, na.rm = TRUE),
              CHNG_BA = sum(CHNG_BA, na.rm = TRUE),
              PREV_TPA = sum(PREV_TPA, na.rm = TRUE),
              PREV_BA = sum(PREV_BA, na.rm = TRUE),
              PREV_RD = mean(PREV_RD, na.rm = TRUE),
              CURR_RD = mean(CURR_RD, na.rm = TRUE),
              plotIn_t = ifelse(sum(plotIn_t, na.rm = TRUE) >  0, 1,0)) %>%
    mutate(FSI = (CURR_RD - PREV_RD) / REMPER) %>%
    ungroup()

  ## If we want to use multiple remeasurements to estimate change,
  ## handle that here
  if (useSeries) {
    ## Get a unique ID for each remeasurement in the series
    nMeas <- t %>%
      distinct(pltID, PLT_CN, MEASYEAR, REMPER) %>%
      group_by(pltID) %>%
      mutate(n = length(unique(PLT_CN)),
             series = min_rank(MEASYEAR)) %>%
      ungroup() %>%
      select(pltID, PLT_CN, REMPER, n, series)

    ## Only if more than one remeasurement available
    if (any(nMeas$n > 1)){

        ## Now we loop over the unique values of n
        ## Basically have to chunk up the data each time
        ## in order to get intermediate estimates
        nRems <- unique(nMeas$n)
        remsList <- list()
        for (i in 1:length(nRems)){
          ## Temporal weights for each plot
          wgts <- nMeas %>%
            filter(series <= nRems[i] & n >= nRems[i]) %>%
            group_by(pltID) %>%
            ## Total remeasurement interval and weights for
            ## individual remeasurements
            mutate(fullRemp = sum(REMPER, na.rm = TRUE),
                   wgt = REMPER / fullRemp) %>%
            ungroup() %>%
            select(PLT_CN, n, series, wgt)

          dat <- tEst %>%
            left_join(wgts, by = c('PLT_CN')) %>%
            filter(series <= nRems[i] & n >= nRems[i]) %>%
            group_by(pltID, PLOT_BASIS, .dots = grpBy[!c(grpBy %in% c('YEAR', 'INVYR'))]) %>%
            summarize(FSI = sum(FSI*wgt, na.rm = TRUE),
                      CHNG_TPA = sum(CHNG_TPA*wgt, na.rm = TRUE),
                      CHNG_BA = sum(CHNG_BA*wgt, na.rm = TRUE),
                      PLT_CN = PLT_CN[which.max(series)],
                      CURR_RD = CURR_RD[which.max(series)],
                      PREV_RD = PREV_RD[which.min(series)],
                      PREV_TPA = PREV_TPA[which.min(series)],
                      PREV_BA = PREV_BA[which.min(series)],
                      plotIn_t = if_else(any(plotIn_t > 0), 1, 0)) %>%
            ungroup() %>%
            select(-c(pltID))
          remsList[[i]] <- dat
        }
        ## Bring it all back together
        dat <- bind_rows(remsList)

        ## Update columns in tEst
        tEst <- tEst %>%
          select(-c(CHNG_TPA:FSI)) %>%
          left_join(dat, by = c('PLT_CN', 'PLOT_BASIS', grpBy[!c(grpBy %in% c('YEAR', 'INVYR'))]))
    }
  }

  ## Now go to strata level and onward
  tEst <- tEst %>%
    ## Rejoin with population tables
    right_join(select(ungroup(popState[[x]]), -c(STATECD)), by = 'PLT_CN') %>%
    ungroup() %>%
    ## Need forest area to adjust SI indices
    #left_join(select(ungroup(aAdj), PLT_CN, aGrps, fa, aAdj), by = c('PLT_CN', aGrps)) %>%
    #Add adjustment factors
    mutate(tAdj = case_when(
      ## When NA, stay NA
      is.na(PLOT_BASIS) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      PLOT_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
      ## Otherwise, use the subpplot value
      PLOT_BASIS == 'SUBP' ~ as.numeric(ADJ_FACTOR_SUBP),
      PLOT_BASIS == 'MICR' ~ as.numeric(ADJ_FACTOR_MICR)),
      CHNG_TPA = CHNG_TPA * tAdj,
      CHNG_BA = CHNG_BA * tAdj,
      CURR_RD = CURR_RD * tAdj,
      PREV_TPA = PREV_TPA * tAdj,
      PREV_BA = PREV_BA * tAdj,
      PREV_RD = PREV_RD * tAdj,
      FSI = FSI * tAdj) %>%
    ## Extra step for variance issues - summing micro, subp, and macr components
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(CHNG_TPA = sum(CHNG_TPA, na.rm = TRUE),
              CHNG_BA = sum(CHNG_BA, na.rm = TRUE),
              PREV_TPA = sum(PREV_TPA, na.rm = TRUE),
              PREV_BA = sum(PREV_BA, na.rm = TRUE),
              PREV_RD = sum(PREV_RD, na.rm = TRUE),
              CURR_RD = sum(CURR_RD, na.rm = TRUE),
              FSI = sum(FSI, na.rm = TRUE),
              plotIn_t = ifelse(sum(plotIn_t, na.rm = TRUE) >  0, 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    ## Add on area
    left_join(select(aAdj, ESTN_UNIT_CN, STRATUM_CN, PLT_CN, aGrps, fa),
              by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'PLT_CN', aGrps)) %>%
    ## FSI is area adjusted
    mutate(si = FSI * fa,
           ra1 = PREV_RD * fa,
           ra2 = CURR_RD * fa) %>%
    ## Replace NAN w/ zeros
    mutate(si = replace_na(si, 0),
           ra1 = replace_na(ra1, 0),
           ra2 = replace_na(ra2, 0)) %>%
    ungroup() %>%
    ## Strata-level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              ctStrat = mean(CHNG_TPA * r_t, na.rm = TRUE),
              cbStrat = mean(CHNG_BA * r_t, na.rm = TRUE),
              ptStrat = mean(PREV_TPA * r_t, na.rm = TRUE),
              pbStrat = mean(PREV_BA * r_t, na.rm = TRUE),
              siStrat = mean(si * r_t, na.rm = TRUE),
              ra1Strat = mean(ra1 * r_t, na.rm = TRUE),
              ra2Strat = mean(ra2 * r_t, na.rm = TRUE),
              faStrat = mean(fa * r_t, na.rm = TRUE),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              ctv = stratVar(ESTN_METHOD, CHNG_TPA, ctStrat, ndif, a, nh),
              cbv = stratVar(ESTN_METHOD, CHNG_BA, cbStrat, ndif, a, nh),
              ptv = stratVar(ESTN_METHOD, PREV_TPA, ptStrat, ndif, a, nh),
              pbv = stratVar(ESTN_METHOD, PREV_BA, pbStrat, ndif, a, nh),
              siv = stratVar(ESTN_METHOD, si, siStrat, ndif, a, nh),
              ra1v = stratVar(ESTN_METHOD, ra1, ra1Strat, ndif, a, nh),
              ra2v = stratVar(ESTN_METHOD, ra2, ra2Strat, ndif, a, nh),
              fav = stratVar(ESTN_METHOD, fa, faStrat, ndif, a, nh),

              # Strata level covariances
              cvStrat_ct = stratVar(ESTN_METHOD, CHNG_TPA, ctStrat, ndif, a, nh, PREV_TPA, ptStrat),
              cvStrat_cb = stratVar(ESTN_METHOD, CHNG_BA, cbStrat, ndif, a, nh, PREV_BA, pbStrat),
              cvStrat_si = stratVar(ESTN_METHOD, si, siStrat, ndif, a, nh, fa, faStrat),
              cvStrat_ra1 = stratVar(ESTN_METHOD, ra1, ra1Strat, ndif, a, nh, fa, faStrat),
              cvStrat_ra2 = stratVar(ESTN_METHOD, ra2, ra2Strat, ndif, a, nh, fa, faStrat),
              cvStrat_psi = stratVar(ESTN_METHOD, si, siStrat, ndif, a, nh, ra1, ra1Strat)) %>%

    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
              cbEst = unitMean(ESTN_METHOD, a, nh, w, cbStrat),
              ptEst = unitMean(ESTN_METHOD, a, nh, w, ptStrat),
              pbEst = unitMean(ESTN_METHOD, a, nh, w, pbStrat),
              siEst = unitMean(ESTN_METHOD, a, nh, w, siStrat),
              ra1Est = unitMean(ESTN_METHOD, a, nh, w, ra1Strat),
              ra2Est = unitMean(ESTN_METHOD, a, nh, w, ra2Strat),
              faEst = unitMean(ESTN_METHOD, a, nh, w, faStrat),

              nh = first(nh),
              # Estimation of unit variance
              ctVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ctv, ctStrat, ctEst),
              cbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cbv, cbStrat, cbEst),
              ptVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ptv, ptStrat, ptEst),
              pbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pbv, pbStrat, pbEst),
              siVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, siv, siStrat, siEst),
              ra1Var = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ra1v, ra1Strat, ra1Est),
              ra2Var = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ra2v, ra2Strat, ra2Est),
              faVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, fav, faStrat, faEst),

              ## Covariances
              cvEst_ct = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ct, ctStrat, ctEst, ptStrat, ptEst),
              cvEst_cb = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_cb, cbStrat, cbEst, pbStrat, pbEst),
              cvEst_si = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_si, siStrat, siEst, faStrat, faEst),
              cvEst_ra1 = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ra1, ra1Strat, ra1Est, faStrat, faEst),
              cvEst_ra2 = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ra2, ra2Strat, ra2Est, faStrat, faEst),
              cvEst_psi = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_psi, siStrat, siEst, ra1Strat, ra1Est),

              plotIn_t = sum(plotIn_t, na.rm = TRUE)) %>%
    ungroup()

  out <- list(tEst = tEst, aEst = NULL)

  return(out)
}



# fsiHelper1_lm <- function(x, plts, db, grpBy, scaleBy, byPlot){
#
#   ## Does not modify outside environment, just need scaleBy in here as well
#   if (is.null(grpBy)){
#     aGrps <- NULL
#     grpBy <- scaleBy
#   } else {
#     aGrps <- grpBy
#     grpBy <- unique(c(grpBy, scaleBy))
#   }
#
#   ## Selecting the plots for one county
#   db$PLOT <- plts[[x]]
#   ## Carrying out filter across all tables
#   #db <- clipFIA(db, mostRecent = FALSE)
#
#   # Only subplots from cond change matrix
#   #db$SUBP_COND_CHNG_MTRX <- filter(db$SUBP_COND_CHNG_MTRX, SUBPTYP == 1)
#
#
#   ## Which grpByNames are in which table? Helps us subset below
#   grpP <- names(db$PLOT)[names(db$PLOT) %in% c(grpBy, scaleBy)]
#   grpC <- names(db$COND)[names(db$COND) %in% c(grpBy, scaleBy) & names(db$COND) %in% grpP == FALSE]
#   grpT <- names(db$TREE)[names(db$TREE) %in% c(grpBy, scaleBy) & names(db$TREE) %in% c(grpP, grpC) == FALSE]
#
#   ## Making a treeID
#   db$TREE$treID <- paste(db$TREE$SUBP, db$TREE$TREE, sep = '_')
#
#   ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
#   data <- db$PLOT %>%
#     filter(DESIGNCD == 1 & PLOT_STATUS_CD != 3) %>%
#     left_join(db$COND, by = c('PLT_CN')) %>%
#     left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
#
#     ## Need a code that tells us where the tree was measured
#     ## macroplot, microplot, subplot
#     mutate(PLOT_BASIS = case_when(
#       ## When DIA is na, adjustment is NA
#       ## NAs occur at harvest and/or recruitement events
#       ## Doesnt matter if we choose SUBP or MICR, adj will by * 0, then additive
#       is.na(DIA) ~ 'SUBP',
#       #is.na(DIA) ~ NA_character_,
#       ## When DIA is less than 5", use microplot value
#       DIA < 5 ~ 'MICR',
#       ## When DIA is greater than 5", use subplot value
#       DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
#       DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
#       DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) %>%
#     ## No NAs allowed as absence values, must be zero
#     mutate(TPA_UNADJ = replace_na(TPA_UNADJ, 0),
#            BAA = replace_na(BAA, 0))
#
#   ## Domain indicator
#   data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
#   data$pDI <- data$landD * data$aD_p * data$aD_c * data$typeD * data$sp
#   data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
#
#
#   ## TOTAL number of observations on each plot
#   obs <- data %>%
#     distinct(PLT_CN, pltID) %>%
#     group_by(pltID) %>%
#     summarize(obs = n()) %>%
#     filter(obs > 1)
#
#   ## TOTAL number of unique stems on each plot through time
#   stems <- data %>%
#     filter(pltID %in% obs$pltID) %>%
#     distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
#     filter(!is.na(treID)) %>%
#     filter(tDI == 1) %>%
#     group_by(pltID, .dots = grpBy[grpBy != 'pltID']) %>%
#     summarize(n = length(unique(treID)))
#
#
#   ## Series of remeasurements for each plot
#   remSeries <- data %>%
#     filter(pltID %in% obs$pltID) %>%
#     distinct(PLT_CN, pltID, MEASYEAR) %>%
#     arrange(pltID, MEASYEAR) %>%
#     group_by(pltID) %>%
#     mutate(meas = min_rank(MEASYEAR))
#   ## Every plot will have a one and a two
#   ## We need to iterate through each of these
#   ## possible 1 - n options
#
#   ## Total trees for the size-density scaling
#   t1 <- data %>%
#     left_join(remSeries, by = c('PLT_CN', 'pltID', 'MEASYEAR')) %>%
#     filter(!c(pltID %in% disturb$pltID)) %>%
#     ## Previous only
#     filter(meas == 1) %>%
#     distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
#     group_by(.dots = scaleBy, PLT_CN, pltID) %>%
#     summarize(REMPER = first(REMPER),
#               BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
#               TPA1 = sum(TPA_UNADJ[STATUSCD == 1] * pDI[STATUSCD == 1], na.rm = TRUE)) %>%
#     ## Change in average tree BA and QMD
#     mutate(BA1 = if_else(TPA1 != 0, BAA / TPA1, 0))
#
#
#   ## TPA and BAA model functions for map
#   ## Returns slope of each series in units of
#   ## annual change
#   t_lm <- function(df){
#     coef(lm(TPA ~ date, data = df))[2]*365
#   }
#   b_lm <- function(df){
#     coef(lm(BAA ~ date, data = df))[2]*365
#   }
#   remper <- function(df){
#     round(as.numeric(max(df$date) - min(df$date)) / 365, 1)
#   }
#   maxYear <- function(df){
#     max(df$MEASYEAR, na.rm = TRUE)
#   }
#   ptpa <- function(df){
#     df$TPA[which.min(df$date)]
#   }
#   pbaa <- function(df){
#     df$BAA[which.min(df$date)]
#   }
#   pstatus <- function(df){
#     if_else(any(df$PLOT_STATUS_CD == 1), 1, 2)
#   }
#   pstatusPop <- function(df){
#     if_else(any(df$plotIn == 1), 1, 0)
#   }
#   completePop <- function(df, grps){
#     df %>%
#       ungroup() %>%
#       complete(nesting(date, MEASYEAR), nesting(PLOT_BASIS, !!!grps)) %>%
#       mutate(TPA = replace_na(TPA, 0),
#              BAA = replace_na(BAA, 0))
#   }
#
#   completePlot <- function(df, grps){
#     df %>%
#       ungroup() %>%
#       complete(nesting(date, MEASYEAR), nesting(!!!grps)) %>%
#       mutate(TPA = replace_na(TPA, 0),
#              BAA = replace_na(BAA, 0))
#   }
#
#   ## Convert to syms for unquoting in dplyr chain
#   grps <- grpBy[!(grpBy %in% c('pltID', 'PLOT_STATUS_CD'))]
#   if (length(grps) > 0) {
#     grps <- syms(grps)
#   } else{
#     grps <- NULL
#   }
#
#
#   if (byPlot){
#     #grpBy <- c('YEAR', grpBy)
#
#     t <- data %>%
#       filter(pltID %in% obs$pltID) %>%
#       mutate(YEAR = MEASYEAR) %>%
#       distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
#       group_by(.dots = grpBy, PLT_CN, PLOT_STATUS_CD, MEASYEAR, MEASMON, MEASDAY) %>%
#       summarize(TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
#                 BAA = sum(BAA * tDI, na.rm = TRUE)) %>%
#       mutate(date = paste(MEASYEAR, MEASMON, MEASDAY, sep = '-'),
#              date = as.Date(date, "%Y-%m-%d")) %>%
#       ungroup() %>%
#       left_join(select(ungroup(remSeries), PLT_CN, meas), by = 'PLT_CN') %>%
#       left_join(select(ungroup(obs), pltID, obs), by = 'pltID')
#
#     tList <- list()
#     ## If plot measured more than twice, we want to return
#     ## intermediate estimates as well
#     for (i in 2:max(t$meas, na.rm = TRUE)){
#       t_int <- t %>%
#         ungroup() %>%
#         filter(meas <= i) %>%
#         filter(obs >= i) %>%
#         #filter(!(PLT_CN %in% cns)) %>%
#         select(-c(PLT_CN, MEASMON, MEASDAY)) %>%
#         ## This step completes the observations - implicit NA becomes explicit 0
#         nest(df = c(MEASYEAR, date, PLOT_STATUS_CD, TPA, BAA, meas, obs, !!!grps)) %>%
#         mutate(df = map(df, completePlot, grps)) %>%
#         unnest(df) %>%
#         ## Nest again, by group, to run the models
#         nest(df = c(MEASYEAR, date, PLOT_STATUS_CD, TPA, BAA, meas, obs)) %>%
#         #group_by(.dots = grpBy) %>%
#         #nest(df = c(MEASYEAR, date, TPA, BAA, meas, obs, grpBy[grpBy != 'pltID'])) %>%
#         mutate(t_rate = map_dbl(df, t_lm),
#                b_rate = map_dbl(df, b_lm),
#                REMPER = map_dbl(df, remper),
#                MEASYEAR = map_dbl(df, maxYear),
#                PREV_TPA = map_dbl(df, ptpa),
#                PREV_BAA = map_dbl(df, pbaa),
#                PLOT_STATUS_CD = map_dbl(df, pstatus),
#                ## For consistency with other function (direct remeasurements)
#                CHNG_TPA = t_rate * REMPER,
#                CHNG_BAA = b_rate * REMPER) %>%
#         ## Replace any NAs with zeros
#         mutate(PREV_BAA = replace_na(PREV_BAA, 0),
#                PREV_TPA = replace_na(PREV_TPA, 0),
#                CHNG_BAA = replace_na(CHNG_BAA, 0),
#                CHNG_TPA = replace_na(CHNG_TPA, 0),
#                ## T2 attributes
#                CURR_BAA = PREV_BAA + CHNG_BAA,
#                CURR_TPA = PREV_TPA + CHNG_TPA) %>%
#         ## Change in average tree BA and QMD
#         mutate(PREV_BA = if_else(PREV_TPA != 0, PREV_BAA / PREV_TPA, 0),
#                CURR_BA = if_else(CURR_TPA != 0, CURR_BAA / CURR_TPA, 0),
#                CHNG_BA = CURR_BA - PREV_BA,
#                PREV_QMD = sqrt(PREV_BA / 0.005454154),
#                CURR_QMD = sqrt(CURR_BA / 0.005454154),
#                CHNG_QMD = CURR_QMD - PREV_QMD) %>%
#         ## All CHNG becomes an annual rate
#         mutate(CHNG_BAA = CHNG_BAA / REMPER,
#                CHNG_TPA = CHNG_TPA / REMPER,
#                CHNG_BA = CHNG_BA / REMPER,
#                CHNG_QMD = CHNG_QMD / REMPER) %>%
#         #left_join(stems, by = c('PLT_CN', grpBy[!(grpBy %in% c('pltID', 'PLOT_STATUS_CD', 'YEAR'))])) %>%
#         ### Then divide by number of unique stems for an average
#         #mutate(CHNG_BAA = CHNG_BAA / n) %>%
#         ungroup() %>%
#         select(MEASYEAR, PLOT_STATUS_CD, REMPER, grpBy,
#                PREV_TPA, PREV_BAA, PREV_BA, PREV_QMD,
#                CHNG_TPA, CHNG_BAA, CHNG_BA, CHNG_QMD,
#                CURR_TPA, CURR_BAA, CURR_BA, CURR_QMD)
#
#       ## update lists
#       #cns <- c(cns, unique(t_int$PLT_CN))
#       tList[[i]] <- t_int
#     }
#
#     ## Back to dataframe
#     t <- bind_rows(tList) %>%
#       left_join(select(ungroup(db$PLOT), PLT_CN, MEASYEAR, pltID), by = c('pltID', 'MEASYEAR')) %>%
#       rename(YEAR = MEASYEAR)
#
#     a = NULL
#
#
#   } else {
#
#     # ### Plot-level estimates -- growth accounting
#     ### Plot-level estimates
#     a <- data %>%
#       filter(pltID %in% obs$pltID) %>%
#       ## date column
#       mutate(date = paste(MEASYEAR, MEASMON, MEASDAY, sep = '-'),
#              date = as.Date(date, "%Y-%m-%d")) %>%
#       ## Will be lots of trees here, so CONDPROP listed multiple times
#       ## Adding PROP_BASIS so we can handle adjustment factors at strata level
#       distinct(PLT_CN, pltID, date, PROP_BASIS, CONDID, .keep_all = TRUE) %>%
#       group_by(PLT_CN, pltID, date, PROP_BASIS, CONDID, .dots = aGrps) %>%
#       summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ * aDI)) %>%
#       mutate(CONDPROP_UNADJ = replace_na(CONDPROP_UNADJ, 0)) %>%
#       ## Average forested area between min and max date
#       group_by(pltID, PROP_BASIS, .dots = aGrps) %>%
#       summarize(minDate = min(date, na.rm = TRUE),
#                 maxDate = max(date, na.rm = TRUE),
#                 amin = sum(CONDPROP_UNADJ[date == minDate], na.rm = TRUE),
#                 amax = sum(CONDPROP_UNADJ[date == maxDate], na.rm = TRUE),
#                 fa = (amin + amax) / 2)
#
#     tStart <- data %>%
#       filter(pltID %in% obs$pltID) %>%
#       distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
#       group_by(.dots = grpBy, PLT_CN, pltID, PLOT_BASIS, MEASYEAR, MEASMON, MEASDAY) %>%
#       summarize(TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
#                 BAA = sum(BAA * tDI, na.rm = TRUE),
#                 plotIn = if_else(sum(tDI, na.rm = TRUE) > 0, 1, 0)) %>%
#       mutate(date = paste(MEASYEAR, MEASMON, MEASDAY, sep = '-'),
#              date = as.Date(date, "%Y-%m-%d")) %>%
#       #filter(!is.na(PLOT_BASIS)) %>%
#       #mutate(PLOT_BASIS = replace_na(PLOT_BASIS, 'SUBP')) %>%
#       ungroup() %>%
#       left_join(select(ungroup(remSeries), PLT_CN, meas), by = 'PLT_CN') %>%
#       left_join(select(ungroup(obs), pltID, obs), by = 'pltID') #%>%
#     #tidyr::complete(nesting(PLT_CN, pltID, date, MEASYEAR, meas, obs, !!!grps), PLOT_BASIS) %>%
#     #mutate(BAA = replace_na(BAA, 0),
#     #       TPA = replace_na(TPA, 0))
#
#     if (nrow(tStart) > 0) {
#       tList <- list()
#       ## If plot measured more than twice, we want to return
#       ## intermediate estimates as well
#       for (i in 2:max(tStart$meas, na.rm = TRUE)){
#         t_int <- tStart %>%
#           filter(meas <= i) %>%
#           filter(obs >= i) %>%
#           #filter(!(PLT_CN %in% cns)) %>%
#           select(-c(PLT_CN, MEASMON, MEASDAY)) %>%
#           #group_by(pltID, PLOT_BASIS, .dots = grpBy) %>%
#           ## This step completes the observations - implicit NA becomes explicit 0
#           nest(df = c(MEASYEAR, date, PLOT_BASIS, TPA, BAA, plotIn, meas, obs, !!!grps)) %>%
#           mutate(df = map(df, completePop, grps)) %>%
#           unnest(df) %>%
#           ungroup() %>%
#           ## Nest again, by group, to run the models
#           #group_by(pltID, PLOT_BASIS, .dots = grpBy) %>%
#           nest(df = c(MEASYEAR, date, TPA, BAA, plotIn, meas, obs)) %>%
#           mutate(t_rate = map_dbl(df, t_lm),
#                  b_rate = map_dbl(df, b_lm),
#                  REMPER = map_dbl(df, remper),
#                  MEASYEAR = map_dbl(df, maxYear),
#                  PREV_TPA = map_dbl(df, ptpa),
#                  PREV_BAA = map_dbl(df, pbaa),
#                  plotIn = map_dbl(df, pstatusPop),
#                  ## For consistency with other function (direct remeasurements)
#                  CHNG_TPA = t_rate * REMPER,
#                  CHNG_BAA = b_rate * REMPER,
#                  CURR_TPA = PREV_TPA + CHNG_TPA,
#                  CURR_BAA = PREV_BAA + CHNG_BAA) %>%
#           select(-c(df))
#         ## update lists
#         #cns <- c(cns, unique(t_int$PLT_CN))
#         tList[[i]] <- t_int
#       }
#
#       ## Back to dataframe
#       t <- bind_rows(tList)  %>%
#         left_join(select(db$PLOT, PLT_CN, MEASYEAR, pltID), by = c('pltID', 'MEASYEAR')) %>%
#         #rename(YEAR = MEASYEAR) %>%
#         mutate(CHNG_BAA = replace_na(CHNG_BAA, 0),
#                CHNG_TPA = replace_na(CHNG_TPA, 0)) %>%
#         filter(!is.na(t_rate)) %>%
#         left_join(select(ungroup(tStart), PLT_CN, PLOT_BASIS, all_of(grpBy)), by = c('PLT_CN', 'PLOT_BASIS', grpBy)) %>%
#         mutate(PREV_PLT_CN = NA)
#     } else {
#       t = NULL
#     }
#
#     ## Need a PLT_CN on a
#     a <- a %>%
#       left_join(select(ungroup(db$PLOT), PLT_CN, pltID), by = c('pltID'))
#
#   }
#
#
#
#   pltOut <- list(t = t, a = a, t1 = t1)
#   return(pltOut)
#
# }






#
#
# fsiHelper1 <- function(x, plts, db, grpBy, byPlot){
#
#   ## Selecting the plots for one county
#   db$PLOT <- plts[[x]]
#
#   ## Making a treeID
#   db$TREE$treID <- paste(db$TREE$SUBP, db$TREE$TREE, sep = '_')
#
#   ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
#   data <- db$PLOT %>%
#     left_join(db$COND, by = c('PLT_CN')) %>%
#     ## AGENTCD at remeasurement, died during the measurement interval
#     left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
#     left_join(db$PLOT, by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('2', '1')) %>%
#     left_join(db$COND, suffix = c('2', '1')) %>%
#     left_join(db$TREE, by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%
#     mutate_if(is.factor,
#               as.character)
#
#   ## Comprehensive indicator function -- w/ growth accounting
#   data$tDI2 <- data$landD2 * data$aD_p2 * data$aD_c2 * data$tD2 * data$typeD2 * data$sp2 *
#     if_else(data$STATUSCD2 == 1, 1, 0) # Live tree
#
#   data$tDI1 <- data$landD1 * data$aD_p1 * data$aD_c1 * data$tD1 * data$typeD1 * data$sp1 *
#     if_else(data$STATUSCD1 == 1, 1, 0) # Live tree
#
#
#   ## PREVIOUS and CURRENT attributes
#   data <- data %>%
#     mutate(TPA_UNADJ1 = TPA_UNADJ1,
#            TPA_UNADJ2 = TPA_UNADJ2,
#            BAA1 = BAA1,
#            BAA2 = BAA2,
#            MORT = case_when(
#              STATUSCD1 == 1 & STATUSCD2 == 2 ~ 1,
#              STATUSCD1 == 1 & STATUSCD2 == 3 ~ 1,
#              TRUE ~ 0),
#            SURV = case_when(
#              STATUSCD1 == 1 & STATUSCD2 == 1 ~ 1,
#              TRUE ~ 0)
#     )
#
#
#   ## Just what we need
#   data <- data %>%
#     select(PLT_CN, PREV_PLT_CN, pltID, TRE_CN, SUBP, CONDID, TREE, CONDPROP_UNADJ,
#            MEASYEAR, MACRO_BREAKPOINT_DIA, PROP_BASIS, grpP[grpP != 'PLOT_STATUS_CD'], grpC,
#            REMPER, PLOT_STATUS_CD1, PLOT_STATUS_CD2,
#            treID1, treID2,
#            one_of(str_c(grpT,1),str_c(grpT,2)),
#            tDI1, tDI2, STATUSCD1, STATUSCD2,
#            DIA1, DIA2, BAA1, BAA2, TPA_UNADJ1, TPA_UNADJ2) %>%
#     ## Negate previous measurements
#     mutate(BAA1 = -(BAA1),
#            TPA_UNADJ1 = -(TPA_UNADJ1)) %>%
#     ## Rearrange previous values as observations
#     pivot_longer(cols = -c(PLT_CN:REMPER),
#                  names_to = c(".value", 'ONEORTWO'),
#                  names_sep = -1) %>%
#     mutate(PLOT_BASIS = case_when(
#       ## When DIA is na, adjustment is NA
#       is.na(DIA) ~ NA_character_,
#       ## When DIA is less than 5", use microplot value
#       DIA < 5 ~ 'MICR',
#       ## When DIA is greater than 5", use subplot value
#       DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
#       DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
#       DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR'))
#
#   ## Under Daves reformulation, we cannot have NAs as absence values
#   ## Change all NA's to zeros for TPA and BAA
#   data <- data %>%
#     mutate(TPA_UNADJ = replace_na(TPA_UNADJ, replace = 0),
#            BAA = replace_na(BAA, replace = 0))
#
#   ## TOTAL number of unique stems on each plot through time
#   stems <- data %>%
#     filter(!is.na(treID)) %>%
#     filter(tDI == 1) %>%
#     group_by(PLT_CN, .dots = grpBy[!(grpBy %in% c('pltID', 'PLOT_STATUS_CD', 'YEAR'))]) %>%
#     summarize(n = length(unique(treID)),
#               n1 = length(unique(treID[ONEORTWO == 1])))
#
#
#   if (byPlot){
#     grpBy <- c('YEAR', grpBy)
#
#     t <- data %>%
#       mutate(YEAR = MEASYEAR) %>%
#       distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
#       # Compute estimates at plot level
#       group_by(.dots = grpBy, PLT_CN, PREV_PLT_CN) %>%
#       summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
#                 REMPER = first(REMPER),
#                 PREV_BAA = sum(-BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
#                 PREV_TPA = sum(-TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
#                 CHNG_TPA = sum(TPA_UNADJ[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
#                 #CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
#                 #CHNG_BAA = mean((BAA[ONEORTWO == 2] * tDI[ONEORTWO == 2]) + (BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1]), na.rm = TRUE),
#                 ## Sum here to avoid issues w/ zeros
#                 CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
#                 CURR_TPA = PREV_TPA + CHNG_TPA,
#                 CURR_BAA = PREV_BAA + CHNG_BAA) %>%
#       left_join(stems, by = c('PLT_CN', grpBy[!(grpBy %in% c('pltID', 'PLOT_STATUS_CD', 'YEAR'))])) %>%
#       ## Then divide by number of unique stems for an average
#       mutate(CHNG_BAA = CHNG_BAA / n) %>%
#       ungroup() %>%
#       select(PLT_CN, PREV_PLT_CN, REMPER, grpBy,
#              PREV_TPA, PREV_BAA, CHNG_TPA, CHNG_BAA, CURR_TPA, CURR_BAA,
#              n, nLive)
#
#   } else {
#
#     ### Compute total TREES in domain of interest
#     t <- data %>%
#       ## Need to count number of trees in each
#       mutate(treID = paste(pltID, SUBP, TREE)) %>%
#       distinct(PLT_CN, TRE_CN, ONEORTWO, .keep_all = TRUE) %>%
#       #group_by(PLT_CN, PLOT_BASIS, TRE_CN) %>%
#       #summarize()
#       # Compute estimates at plot level
#       group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
#       summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
#                 REMPER = first(REMPER),
#                 MEASYEAR = first(MEASYEAR),
#                 PREV_TPA = sum(-TPA_UNADJ[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
#                 PREV_BAA = sum(-BAA[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
#                 CHNG_TPA = sum(TPA_UNADJ[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
#                 CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
#                 plotIn = if_else(sum(tDI, na.rm = TRUE) > 0, 1, 0)) %>%
#       left_join(stems, by = c('PLT_CN', grpBy))
#
#   }
#
#   pltOut <- list(t = t)
#   return(pltOut)
#
# }
#
#
# fsiHelper2 <- function(x, popState, t, grpBy, method){
#
#   ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
#   if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
#     grpBy <- c(grpBy, 'INVYR')
#     popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
#     popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
#
#   }
#
#
#   ## Strata level estimates
#   tEst <- t %>%
#     ## Rejoin with population tables
#     right_join(select(popState[[x]], -c(STATECD, REMPER)), by = 'PLT_CN') %>%
#     ungroup() %>%
#     #Add adjustment factors
#     mutate(tAdj = case_when(
#       ## When NA, stay NA
#       is.na(PLOT_BASIS) ~ NA_real_,
#       ## If the proportion was measured for a macroplot,
#       ## use the macroplot value
#       PLOT_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
#       ## Otherwise, use the subpplot value
#       PLOT_BASIS == 'SUBP' ~ as.numeric(ADJ_FACTOR_SUBP),
#       PLOT_BASIS == 'MICR' ~ as.numeric(ADJ_FACTOR_MICR)),
#       ct = CHNG_TPA * tAdj,
#       cb = CHNG_BAA * tAdj,
#       pt = PREV_TPA * tAdj,
#       pb = PREV_BAA * tAdj) %>%
#     ## Computing change
#     mutate(ct = (ct) / REMPER,
#            cb = (cb) / REMPER) %>%
#     ## Extra step for variance issues
#     ## Summing across micro, subp, macro
#     group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
#     summarize(ctPlot = sum(ct, na.rm = TRUE),
#               cbPlot = sum(cb, na.rm = TRUE),
#               ptPlot = sum(pt, na.rm = TRUE),
#               pbPlot = sum(pb, na.rm = TRUE),
#               plotIn = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
#               nh = first(P2POINTCNT),
#               p2eu = first(p2eu),
#               a = first(AREA_USED),
#               w = first(P1POINTCNT) / first(P1PNTCNT_EU),
#               REMPER = first(REMPER),
#               ## Total unique number of trees
#               n = first(n)) %>%
#
#
#     ## Do not want to compute SI for micro and subp seperately, handle it here
#     mutate(TPA_RATE = TPA_RATE / REMPER,
#            BAA_RATE = BAA_RATE / REMPER / n,
#            #CURR_TPA = ptPlot + (ctPlot * REMPER),
#            #CURR_BAA = pbPlot + (cbPlot * REMPER),
#            #TPA_RATE = ctPlot, #/ (CURR_TPA + ptPlot) * 2,
#            #BAA_RATE = cbPlot, #/ (CURR_BAA + pbPlot) * 2,
#            x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
#            #y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
#            siPlot = sqrt(x^2 + x^2),
#            siPlot = if_else(x < 0, -siPlot, siPlot),
#            siPlot = case_when(
#              is.na(siPlot) ~ 0,
#              TRUE ~ siPlot)
#     ) %>%
#     ungroup() %>%
#     left_join(select(aAdj, PLT_CN, aGrps, fa), by = c('PLT_CN', aGrps)) %>%
#
#     group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
#     summarize(r_t = length(unique(PLT_CN)) / first(nh),
#               ctStrat = mean(ctPlot * r_t, na.rm = TRUE),
#               cbStrat = mean(cbPlot * r_t, na.rm = TRUE),
#               ptStrat = mean(ptPlot * r_t, na.rm = TRUE),
#               pbStrat = mean(pbPlot * r_t, na.rm = TRUE),
#               siStrat = mean(siPlot * fa * r_t, na.rm = TRUE),
#               faStrat = mean(fa * r_t, na.rm = TRUE),
#               plotIn_t = sum(plotIn_t, na.rm = TRUE),
#               n = n(),
#               ## We don't want a vector of these values, since they are repeated
#               nh = first(nh),
#               a = first(a),
#               w = first(w),
#               p2eu = first(p2eu),
#               ndif = nh - n,
#               # ## Strata level variances
#               ctv = stratVar(ESTN_METHOD, ctPlot, ctStrat, ndif, a, nh),
#               cbv = stratVar(ESTN_METHOD, cbPlot, cbStrat, ndif, a, nh),
#               ptv = stratVar(ESTN_METHOD, ptPlot, ptStrat, ndif, a, nh),
#               pbv = stratVar(ESTN_METHOD, pbPlot, pbStrat, ndif, a, nh),
#               siv = stratVar(ESTN_METHOD, siPlot * fa, siStrat, ndif, a, nh),
#               fav = stratVar(ESTN_METHOD, fa, faStrat, ndif, a, nh),
#
#
#               # Strata level covariances
#               cvStrat_ct = stratVar(ESTN_METHOD, ctPlot, ctStrat, ndif, a, nh, ptPlot, ptStrat),
#               cvStrat_cb = stratVar(ESTN_METHOD, cbPlot, cbStrat, ndif, a, nh, pbPlot, pbStrat),
#               cvStrat_si = stratVar(ESTN_METHOD, siPlot * fa, siStrat, ndif, a, nh, fa, faStrat)) %>%
#
#     ## Estimation unit
#     group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
#     summarize(ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
#               cbEst = unitMean(ESTN_METHOD, a, nh, w, cbStrat),
#               ptEst = unitMean(ESTN_METHOD, a, nh, w, ptStrat),
#               pbEst = unitMean(ESTN_METHOD, a, nh, w, pbStrat),
#               siEst = unitMean(ESTN_METHOD, a, nh, w, siStrat),
#               faEst = unitMean(ESTN_METHOD, a, nh, w, faStrat),
#
#               nh = first(nh),
#               # Estimation of unit variance
#               ctVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ctv, ctStrat, ctEst),
#               cbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cbv, cbStrat, cbEst),
#               ptVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ptv, ptStrat, ptEst),
#               pbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pbv, pbStrat, pbEst),
#               siVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, siv, siStrat, siEst),
#               faVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, fav, faStrat, faEst),
#
#               ## Covariances
#               cvEst_ct = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ct, ctStrat, ctEst, ptStrat, ptEst),
#               cvEst_cb = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_cb, cbStrat, cbEst, pbStrat, pbEst),
#               cvEst_si = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_si, siStrat, siEst, faStrat, faEst),
#
#               plotIn_t = sum(plotIn_t, na.rm = TRUE)) %>%
#     ungroup()
#
#   out <- list(tEst = tEst, aEst = NULL)
#
#   return(out)
# }
