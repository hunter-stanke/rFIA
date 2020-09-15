invHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)

  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$SUBP_COND, by = c('PLT_CN', 'CONDID')) %>%
    left_join(db$INVASIVE_SUBPLOT_SPP, by = c("PLT_CN", "CONDID", 'SUBP'))
  #filter(DIA >= 5)

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp



  if (byPlot){
    ## Return proportion of plot area covered by invasive
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      filter(!is.na(SYMBOL)) %>%
      distinct(PLT_CN, CONDID, SUBP, VEG_SPCD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, SUBP) %>%
      summarize(cover = sum(COVER_PCT/100 * SUBPCOND_PROP * aDI, na.rm = TRUE)) %>%
      group_by(PLT_CN, .dots = grpBy) %>%
      summarize(cover = mean(cover, na.rm = TRUE))
    a = NULL

  } else {
    # Compute estimates
    t <- data %>%
      distinct(PLT_CN, CONDID, SUBP, VEG_SPCD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, SUBP) %>%
      summarize(iPlot = sum(COVER_PCT/100 * SUBPCOND_PROP * aDI, na.rm = TRUE),
                plotIn_INV = ifelse(sum(aDI, na.rm = TRUE) >  0 & first(INVASIVE_SAMPLING_STATUS_CD) == 1, 1,0)) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(iPlot = mean(iPlot, na.rm = TRUE),
                plotIn_INV = ifelse(any(plotIn_INV > 0), 1, 0))

    a <- data %>%
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, PLT_CN, PROP_BASIS) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn_AREA = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))
  }

  pltOut <- list(t = t, a = a)
  return(pltOut)

}



invHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    aGrpBy <- c(aGrpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR

  }

  ## Strata level estimates
  aStrat <- a %>%
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
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = aGrpBy) %>%
    summarize(a_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              aStrat = mean(fa * a_t, na.rm = TRUE),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              p2eu = first(p2eu),
              ndif = nh - n,
              ## Strata level variances
              av = stratVar(ESTN_METHOD, fa, aStrat, ndif, a, nh))
  ## Estimation unit
  aEst <- aStrat %>%
    group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, av, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))


  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
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
      fa = fa * aAdj,
      iPlot = iPlot * ADJ_FACTOR_SUBP) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              iStrat = mean(iPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_INV = sum(plotIn_INV, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              p2eu = first(p2eu),
              ndif = nh - n,
              ## Strata level variances
              ## Strata level variances
              iv = stratVar(ESTN_METHOD, iPlot, iStrat, ndif, a, nh),
              cvStrat_i = stratVar(ESTN_METHOD, iPlot, iStrat, ndif, a, nh, fa, aStrat)) %>%
    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(iEst = unitMean(ESTN_METHOD, a, nh,  w, iStrat),
              N = first(p2eu),
              plotIn_INV = sum(plotIn_INV, na.rm = TRUE),
              iVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, iv, iStrat, iEst),
              # Unit Covariance
              cvEst_i = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_i, iStrat, iEst, aStrat, aEst))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}








invasiveHelper <- function(x, combos, data, grpBy, aGrpBy, totals, SE){

  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  ad = 1

  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(tObs))) == length(tObs)) tObs <- 1
    td <- data$aDI * tObs * td
    # Area domain indicator for each column in
    if(grpBy[n] %in% aGrpBy){
      aObs <- as.character(combos[[x]][[aGrpBy[n]]]) == as.character(data[[aGrpBy[n]]])
      if (length(which(is.na(aObs))) == length(aObs)) aObs <- 1
      #aObs[is.na(aObs)] <- 0
      ad <- data$aDI * aObs * ad

    }
  }


  if(SE){
    data$tDI <- td
    data$tDI[is.na(data$tDI)] <- 0
    data$aDI <- ad
    data$aDI[is.na(data$aDI)] <- 0
    invPrep <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SYMBOL, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(iPlot = ifelse(any(INVASIVE_SAMPLING_STATUS_CD == 1),
                               sum(COVER_PCT/100 * CONDPROP_UNADJ * aAdj * tDI, na.rm = TRUE),
                               NA),
                iPlots = ifelse(sum(tDI > 0 & iPlot > 0, na.rm = TRUE), 1, 0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT)) #%>%
    inv <- invPrep %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(iStrat = mean(iPlot, na.rm = TRUE),
                iPlots = sum(iPlots, na.rm = TRUE),
                a = first(a),
                w = first(p1)/first(p1EU),
                nh = first(p2),
                iv = ifelse(first(ESTN_METHOD == "simple"),
                            var(iPlot * first(a)/iPlots),
                            (sum(iPlot^2, na.rm = TRUE) - sum(iPlots * iStrat^2))/(iPlots * (iPlots - 1)))) %>%
      group_by(ESTN_UNIT_CN) %>%
      summarize(iEst = unitMean(ESTN_METHOD, a, iPlots, w, iStrat),
                iPlots = sum(iPlots, na.rm = TRUE),
                iVar = unitVar(method = "var", ESTN_METHOD, a, iPlots, w, iv, iStrat, iEst)) %>%
      summarize(COVER_AREA = sum(iEst, na.rm = TRUE),
                iVar = sum(iVar, na.rm = TRUE),
                COVER_AREA_SE = sqrt(iVar)/COVER_AREA * 100,
                nPlots_INV = sum(iPlots, na.rm = TRUE))

    aPrep <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
                plotIn = ifelse(sum(aDI > 0, na.rm = TRUE), 1, 0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT)) #%>%
    a <- aPrep %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = first(a),
                w = first(p1)/first(p1EU),
                nh = first(p2),
                av = ifelse(first(ESTN_METHOD == "simple"),
                            var(fa * first(a)/nh),
                            (sum(fa^2) - sum(nh * aStrat^2))/(nh * (nh - 1)))) %>%
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                plotIn = sum(plotIn, na.rm = TRUE),
                aVar = unitVar(method = "var", ESTN_METHOD, a, nh, w, av, aStrat, aEst)) %>%
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                aVar = sum(aVar, na.rm = TRUE),
                AREA_TOTAL_SE = sqrt(aVar)/AREA_TOTAL * 100,
                nPlots_AREA = sum(plotIn, na.rm = TRUE))


    ## Compute COVARIANCE between numerator and denominator (for ratio estimates of variance)
    covar <- invPrep %>%
      inner_join(aPrep, by = c('PLT_CN', 'ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN'), suffix = c('_i', '_a')) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                iStrat = mean(iPlot, na.rm = TRUE),
                a = first(a_i),
                w = first(p1_i) / first(p1EU_a), # Stratum weight
                nh = first(p2_i), # Number plots in stratum
                # Strata level covariances
                cvStrat = ifelse(first(ESTN_METHOD == "simple"),
                                 cov(fa, iPlot),
                                 (sum(fa * iPlot, na.rm = TRUE) - sum(nh * aStrat * iStrat))/(nh * (nh - 1)))) %>% # Stratified and double cases)
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                iEst = unitMean(ESTN_METHOD, a, iPlots, w, iStrat),
                cvEst = unitVar(method = "cov", ESTN_METHOD, a, nh, w, cvStrat, iStrat, iEst, aStrat, aEst)) %>%
      summarize(cv = sum(cvEst, na.rm = TRUE))

    ## Join up tree, area, and covariance tables
    inv <- data.frame(inv, a, covar)

    ## Make TPA, BAA and SE
    inv <- inv %>%
      mutate(COVER_PCT = COVER_AREA / AREA_TOTAL * 100,
             coverVar = (1/AREA_TOTAL^2) * (iVar + (COVER_PCT^2 * aVar) - 2 * COVER_PCT * cv),
             COVER_PCT_SE = sqrt(coverVar)/COVER_PCT * 100)

    if (totals) {
      inv <- inv %>%
        select(COVER_PCT, COVER_AREA, AREA_TOTAL, COVER_PCT_SE, COVER_AREA_SE, AREA_TOTAL_SE, nPlots_INV, nPlots_AREA)
    }
    else {
      inv <- inv %>%
        select(COVER_PCT, COVER_PCT_SE, nPlots_INV, nPlots_AREA)
    }

    inv <- data.frame(combos[[x]], inv) %>%
    filter(!is.na(SYMBOL))

  } else { # IF SE is FALSE
    a <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn = ifelse(sum(aDI > 0, na.rm = TRUE), 1, 0)) %>%
      group_by(.dots = aGrpBy) %>%
      summarize(AREA_TOTAL = sum(fa, na.rm = TRUE),
                nPlots_AREA = sum(plotIn, na.rm = TRUE))
    # Used to make a new EXPNS value below, number of plots where invasives were sampled
    nPlots <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, STRATUM_CN, PLT_CN) %>%
      summarize(plotIn = ifelse(sum(INVASIVE_SAMPLING_STATUS_CD == 1, na.rm = TRUE), 1, 0)) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, STRATUM_CN) %>%
      summarize(plotsInv = sum(plotIn, na.rm = TRUE))

    inv <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SYMBOL, .keep_all = TRUE) %>%
      inner_join(nPlots, by = c(aGrpBy, 'ESTN_UNIT_CN', 'STRATUM_CN')) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(expInv = first(AREA_USED)*first(P1POINTCNT) / first(P1PNTCNT_EU) / first(plotsInv),
                iPlot = ifelse(any(INVASIVE_SAMPLING_STATUS_CD == 1),
                               sum(COVER_PCT/100 * CONDPROP_UNADJ * aAdj * expInv * aDI, na.rm = TRUE),
                               NA)) %>%
      group_by(.dots = grpBy) %>%
      summarize(COVER_AREA = sum(iPlot, na.rm = TRUE),
                nPlots_INV = length(which(!is.na(iPlot))))


    # Rejoin
    suppressMessages({
      inv <- inner_join(inv, a) %>%
        filter(!is.na(SYMBOL)) %>%
        mutate(COVER_PCT = COVER_AREA / AREA_TOTAL * 100)
    })



    if (totals) {
      inv <- inv %>%
        select(grpBy, COVER_PCT, COVER_AREA, AREA_TOTAL, nPlots_INV, nPlots_AREA)
    }
    else {
      inv <- inv %>%
        select(grpBy, COVER_PCT, nPlots_INV, nPlots_AREA)
    }

  } # End SE conditional

  # Do some cleanup
  #gc()

  #Return a dataframe
  inv
}
