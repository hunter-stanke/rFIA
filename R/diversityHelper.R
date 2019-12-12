divHelper1 <- function(x, plts, db, grpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'grp', 'state', 'SUBP', 'TREE', grpT, 'tD', 'typeD')), by = c('PLT_CN', 'CONDID')) %>%
    ## Need a code that tells us where the tree was measured
    ## macroplot, microplot, subplot
    mutate(PLOT_BASIS = case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA) ~ NA_character_,
      ## When DIA is less than 5", use microplot value
      DIA < 5 ~ 'MICR',
      ## When DIA is greater than 5", use subplot value
      DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
      DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) #%>%


  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
  data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, CONDID, TREE, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(H = divIndex(grp, state  * tDI, index = 'H'),
                S = divIndex(grp, state * tDI, index = 'S'),
                Eh = divIndex(grp, state * tDI, index = 'Eh'),
                nStems = length(which(tDI == 1)))
    a = NULL

  } else {
    # Diversity is computed at the stand (condition level), and we continue to use the ratio of means estimator to get at
    #  average of the attribute of interest weighted by the area in which it occurs.
    t <- data %>%
      distinct(PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, PROP_BASIS, CONDID) %>%
      # filter(tDI > 0) %>%
      summarize(condArea = first(CONDPROP_UNADJ),
                hCond = divIndex(grp, state * tDI, index = 'H') * condArea,
                sCond = divIndex(grp, state * tDI, index = 'S') * condArea,
                EhCond = divIndex(grp, state * tDI, index = 'Eh') * condArea,
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0))  %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(hPlot = sum(hCond, na.rm = TRUE),
                EhPlot = sum(EhCond, na.rm = TRUE),
                sPlot = sum(sCond * plotIn, na.rm = TRUE),
                fa = sum(condArea * aDI, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE))
  }

  pltOut <- list(t = t)
  return(pltOut)
}



bioHelper2 <- function(x, popState, t, grpBy){

  ## Strata level estimates
  tEst <- t %>%
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
      fa = fa * aAdj,
      hPlot = hPlot * aAdj,
      EhPlot = EhPlot * aAdj,
      sPlot = sPlot * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(a_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              aStrat = mean(fa * a_t, na.rm = TRUE),
              hStrat = mean(hPlot * a_t, na.rm = TRUE),
              ehStrat = mean(EhPlot * a_t, na.rm = TRUE),
              sStrat = mean(sPlot * a_t, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              p2eu = first(p2eu),
              ndif = nh - n,
              ## Strata level variances
              av = ifelse(first(ESTN_METHOD == 'simple'),
                          var(c(fa, numeric(ndif)) * first(a) / nh),
                          (sum((c(fa, numeric(ndif))^2)) - nh * aStrat^2) / (nh * (nh-1))))
  ## Estimation unit
  aEst <- aStrat %>%
    group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, av, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

  ######## ------------------ TREE ESTIMATES + CV

  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
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
      nvPlot = nvPlot * tAdj,
      svPlot = svPlot * tAdj,
      bagPlot = bagPlot* tAdj,
      bbgPlot = bbgPlot* tAdj,
      btPlot = btPlot* tAdj,
      cagPlot = cagPlot* tAdj,
      cbgPlot = cbgPlot* tAdj,
      ctPlot = ctPlot* tAdj) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(nvPlot = sum(nvPlot, na.rm = TRUE),
              svPlot = sum(svPlot = svPlot, na.rm = TRUE),
              bagPlot = sum(bagPlot = bagPlot, na.rm = TRUE),
              bbgPlot = sum(bbgPlot = bbgPlot, na.rm = TRUE),
              btPlot = sum(btPlot, na.rm = TRUE),
              cagPlot = sum(cagPlot, na.rm = TRUE),
              cbgPlot = sum(cbgPlot, na.rm = TRUE),
              ctPlot = sum(ctPlot, na.rm = TRUE),
              fa = first(fa),
              plotIn = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              nvStrat = mean(nvPlot * r_t, na.rm = TRUE),
              svStrat = mean(svPlot * r_t, na.rm = TRUE),
              bagStrat = mean(bagPlot * r_t, na.rm = TRUE),
              bbgStrat = mean(bbgPlot * r_t, na.rm = TRUE),
              btStrat = mean(btPlot * r_t, na.rm = TRUE),
              cagStrat = mean(cagPlot * r_t, na.rm = TRUE),
              cbgStrat = mean(cbgPlot * r_t, na.rm = TRUE),
              ctStrat = mean(ctPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_TREE = sum(plotIn, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              nvv = stratVar(ESTN_METHOD, nvPlot, nvStrat, ndif, a, nh),
              svv = stratVar(ESTN_METHOD, svPlot, svStrat, ndif, a, nh),
              bagv = stratVar(ESTN_METHOD, bagPlot, bagStrat, ndif, a, nh),
              bbgv = stratVar(ESTN_METHOD, bbgPlot, bbgStrat, ndif, a, nh),
              btv = stratVar(ESTN_METHOD, btPlot, btStrat, ndif, a, nh),
              cagv = stratVar(ESTN_METHOD, cagPlot, cagStrat, ndif, a, nh),
              cbgv = stratVar(ESTN_METHOD, cbgPlot, cbgStrat, ndif, a, nh),
              ctv = stratVar(ESTN_METHOD, ctPlot, ctStrat, ndif, a, nh),
              # Strata level covariances
              cvStrat_nv = stratVar(ESTN_METHOD, nvPlot, nvStrat, ndif, a, nh, fa, aStrat),
              cvStrat_sv = stratVar(ESTN_METHOD, svPlot, svStrat, ndif, a, nh, fa, aStrat),
              cvStrat_bag = stratVar(ESTN_METHOD, bagPlot, bagStrat, ndif, a, nh, fa, aStrat),
              cvStrat_bbg = stratVar(ESTN_METHOD, bbgPlot, bbgStrat, ndif, a, nh, fa, aStrat),
              cvStrat_bt = stratVar(ESTN_METHOD, btPlot, btStrat, ndif, a, nh, fa, aStrat),
              cvStrat_cag = stratVar(ESTN_METHOD, cagPlot, cagStrat, ndif, a, nh, fa, aStrat),
              cvStrat_cbg = stratVar(ESTN_METHOD, cagPlot, cagStrat, ndif, a, nh, fa, aStrat),
              cvStrat_ct = stratVar(ESTN_METHOD, cagPlot, cagStrat, ndif, a, nh, fa, aStrat)
    ) %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(nvEst = unitMean(ESTN_METHOD, a, nh, w, nvStrat),
              svEst = unitMean(ESTN_METHOD, a, nh, w, svStrat),
              bagEst = unitMean(ESTN_METHOD, a, nh, w, bagStrat),
              bbgEst = unitMean(ESTN_METHOD, a, nh, w, bbgStrat),
              btEst = unitMean(ESTN_METHOD, a, nh, w, btStrat),
              cagEst = unitMean(ESTN_METHOD, a, nh, w, cagStrat),
              cbgEst = unitMean(ESTN_METHOD, a, nh, w, cbgStrat),
              ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
              # Estimation of unit variance
              nvVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, nvv, nvStrat, nvEst),
              svVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, svv, svStrat, svEst),
              bagVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bagv, bagStrat, bagEst),
              bbgVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bbgv, bbgStrat, bbgEst),
              btVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, btv, btStrat, btEst),
              cagVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cagv, cagStrat, cagEst),
              cbgVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cbgv, cbgStrat, cbgEst),
              ctVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ctv, ctStrat, ctEst),
              cvEst_nv = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_nv, nvStrat, nvEst, aStrat, aEst),
              cvEst_sv = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_sv, svStrat, svEst, aStrat, aEst),
              cvEst_bag = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_bag, bagStrat, bagEst, aStrat, aEst),
              cvEst_bbg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_bbg, bbgStrat, bbgEst, aStrat, aEst),
              cvEst_bt = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_bt, btStrat, btEst, aStrat, aEst),
              cvEst_cag = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_cag, cagStrat, cagEst, aStrat, aEst),
              cvEst_cbg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_cbg, cbgStrat, cbgEst, aStrat, aEst),
              cvEst_ct = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ct, ctStrat, ctEst, aStrat, aEst),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}






diversityHelper <- function(x, combos, data, grpBy, SE){
  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  #ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(tObs))) == length(tObs)) tObs <- 1
    td <- data$tDI * tObs * td
  }


  # IF we want sampling errors returned
  if(SE){
    data$tDI <- td
    data$tDI[is.na(data$tDI)] <- 0
    ## We produce an intermediate object in this chain as it is needed to compute the ratio of means variance
    ## Numerator and denominator are in different domains of interest, and may be grouped by different variables
    ## see covariance estimation below

    # Diversity is computed at the stand (condition level), and we continue to use the ratio of means estimator to get at
    #  average of the attribute of interest weighted by the area in which it occurs.
    d <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      # filter(tDI > 0) %>%
      summarize(condArea = first(CONDPROP_UNADJ),
                hCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'H') * condArea,
                sCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S') * condArea,
                EhCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh') * condArea,
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT),
                aAdj = first(aAdj),
                tAdj = first(tAdj),
                test = length(unique(SPCD)),
                nstems = sum(tDI, na.rm = TRUE))  %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(hPlot = sum(hCond * tAdj, na.rm = TRUE),
                EhPlot = sum(EhCond * tAdj, na.rm = TRUE),
                sPlot = sum(sCond * tAdj * plotIn, na.rm = TRUE),
                fa = sum(condArea * aDI * aAdj, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = first(a),
                p1EU = first(p1EU),
                p1 = first(p1),
                p2 = first(p2)) %>%
      # Continue through totals
      #d <- dInt %>%
      #filter(plotIn > 0) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                hStrat = mean(hPlot, na.rm = TRUE),
                EhStrat = mean(EhPlot, na.rm = TRUE),
                sStrat = mean(sPlot, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = first(a),
                w = first(p1) / first(p1EU), # Stratum weight
                nh = first(p2), # Number plots in stratum
                # Strata level variances
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * first(a) / nh),
                            (sum(fa^2, na.rm = TRUE) - sum(nh * aStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
                hv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(hPlot * first(a) / nh),
                            (sum(hPlot^2, na.rm = TRUE) - sum(nh * hStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                Ehv = ifelse(first(ESTN_METHOD == 'simple'),
                             var(EhPlot * first(a) / nh),
                             (sum(EhPlot^2, na.rm = TRUE) - sum(nh * EhStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                sv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(sPlot * first(a) / nh),
                            (sum(sPlot^2, na.rm = TRUE) - sum(nh * sStrat^2)) / (nh * (nh-1))),
                cvStrat_h = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(fa,hPlot),
                                   (sum(fa*hPlot) - sum(nh * aStrat *hStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_Eh = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(fa,EhPlot),
                                    (sum(fa*EhPlot) - sum(nh * aStrat *EhStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_s = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(fa,sPlot),
                                   (sum(fa*sPlot) - sum(nh * aStrat *sStrat)) / (nh * (nh-1))))  %>% # Stratified and double cases
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                h = unitMean(ESTN_METHOD, a, nh, w, hStrat),
                eh = unitMean(ESTN_METHOD, a, nh, w, EhStrat),
                s = unitMean(ESTN_METHOD, a, nh, w, sStrat),
                plotIn = sum(plotIn, na.rm = TRUE),
                # Estimation of unit variance
                hVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, hv, hStrat, h),
                ehVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, Ehv, EhStrat, eh),
                sVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, sv, sStrat, s),
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
                cvEst_h = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_h, hStrat, h, aStrat, aEst),
                cvEst_eh = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_Eh, EhStrat, eh, aStrat, aEst),
                cvEst_s = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_s, sStrat, s, aStrat, aEst)) %>%
      # Compute totals
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                H_a = sum(h, na.rm = TRUE) / AREA_TOTAL,
                Eh_a = sum(eh, na.rm = TRUE) / AREA_TOTAL,
                S_a = sum(s, na.rm = TRUE) / AREA_TOTAL,
                hVar = sum(hVar, na.rm = TRUE),
                ehVar = sum(ehVar, na.rm = TRUE),
                sVar = sum(sVar, na.rm = TRUE),
                nStands = sum(plotIn, na.rm = TRUE),
                areaVar = sum(aVar, na.rm = TRUE),
                cvH = sum(cvEst_h, na.rm = TRUE),
                cveH = sum(cvEst_eh, na.rm = TRUE),
                cvS = sum(cvEst_s, na.rm = TRUE),
                hVar = (1/AREA_TOTAL^2) * (hVar + (H_a^2 * areaVar) - 2 * H_a * cvH),
                ehVar = (1/AREA_TOTAL^2) * (ehVar + (Eh_a^2 * areaVar) - 2 * Eh_a * cveH),
                sVar = (1/AREA_TOTAL^2) * (sVar + (S_a^2 * areaVar) - 2 * S_a * cvS),
                H_a_SE = sqrt(hVar) / H_a * 100,
                Eh_a_SE = sqrt(ehVar) / Eh_a * 100,
                S_a_SE = sqrt(sVar) / S_a * 100) %>%
      select(H_a, Eh_a, S_a, H_a_SE, Eh_a_SE, S_a_SE, nStands)

    # Beta & gamma diversity indices
    dbg <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by() %>%
      summarize(H_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'H'),
                H_b = H_g - d$H_a,
                Eh_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh'),
                Eh_b = Eh_g - d$Eh_a,
                S_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S'),
                S_b = S_g - d$S_a)

    # Join up the alpha beta gamma
    d <- data.frame(d, dbg) %>%
      select(H_a, H_b, H_g, Eh_a, Eh_b, Eh_g, S_a, S_b, S_g, H_a_SE, Eh_a_SE, S_a_SE, nStands)

    # Rejoin with groupby
    d <- data.frame(combos[[x]], d)

  } else { # No sampling errors
    ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
    d <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      summarize(aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0),
                condArea = first(CONDPROP_UNADJ),
                hCond = divIndex(SPCD, TPA_UNADJ  * tDI, index = 'H') * condArea,
                sCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S') * condArea,
                EhCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh') * condArea,
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                EXPNS = first(EXPNS),
                tAdj = first(tAdj),
                aAdj = first(aAdj)) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(condArea * aDI * aAdj * EXPNS, na.rm = TRUE),
                hPlot = sum(hCond * EXPNS * tAdj, na.rm = TRUE),
                EhPlot = sum(EhCond * EXPNS * tAdj, na.rm = TRUE),
                sPlot = sum(sCond * EXPNS * tAdj * plotIn, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE)) %>%
      group_by(.dots = grpBy) %>%
      summarize(AREA_TOTAL = sum(fa, na.rm = TRUE),
                H_a = sum(hPlot, na.rm = TRUE) / AREA_TOTAL,
                Eh_a = sum(EhPlot, na.rm = TRUE) / AREA_TOTAL,
                S_a = sum(sPlot, na.rm = TRUE) / AREA_TOTAL,
                nStands = sum(plotIn, na.rm = TRUE)) #%>%
    #filter(S > 0) #%>%
    #select(c(grpByOrig, 'H_a', 'Eh_a', 'S_a', 'nStands'))

    # Beta & gamma diversity indices
    dbg <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy) %>%
      summarize(H_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'H'),
                Eh_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh'),
                S_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S'))
    # Join up the alpha beta gamma
    suppressMessages({
      d <- inner_join(d, dbg) %>%
        mutate(H_b = H_g - H_a,
               Eh_b = Eh_g - Eh_a,
               S_b = S_g - S_a) %>%
        ungroup() %>%
        select(grpBy, H_a, H_b, H_g, Eh_a, Eh_b, Eh_g, S_a, S_b, S_g, nStands)
    })
  }

  # Do some cleanup
  #gc()

  #Return a dataframe
  d

}
