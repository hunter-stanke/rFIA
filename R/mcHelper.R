mcHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)

  # Only subplots from cond change matrix
  #db$SUBP_COND_CHNG_MTRX <- filter(db$SUBP_COND_CHNG_MTRX, SUBPTYP == 1)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', 'STATUSCD', 'DIA', 'PREVDIA', 'TPA_UNADJ', grpT, 'tD', 'typeD')), by = c('PLT_CN', 'CONDID')) %>%
    #left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPARECR_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
    #left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
    #left_join(select(db$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'CONDID'), suffix = c('', '.subp')) %>%
    #left_join(select(db$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN.subp' = 'PLT_CN', 'PREVCOND.subp' = 'CONDID'), suffix = c('', '.chng')) %>%
    left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
    left_join(select(db$TREE, c('TRE_CN', 'DIA', grpT, 'typeD', 'tD', 'STATUSCD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
    mutate_if(is.factor,
              as.character) %>%
    ## Identify Mortality & Recruitment
    mutate(mort = if_else(STATUSCD %in% c(2,3) & STATUSCD.prev == 1 & DIA >=5, 1, 0),
           ## Ingrowth on subplot
           recr1 = if_else(!is.na(PREV_PLT_CN) & is.na(PREV_TRE_CN) & DIA >=5 & STATUSCD == 1, 1, 0),
           recr2 = if_else(DIA >= 5 & PREVDIA < 5 & STATUSCD == 1, 1, 0)) %>%
    ## Make some columns for easy use later
    mutate(M_TPA = case_when(
              mort == 1 ~ TPA_UNADJ,
              TRUE ~ 0),
           R_TPA = case_when(
             recr1 == 1 | recr2 == 1 ~ TPA_UNADJ,
             TRUE ~ 0),
           T_TPA = case_when(
             DIA >=5 & STATUSCD == 1 ~ TPA_UNADJ,
             TRUE ~ 0),
           BA1 = case_when(
             STATUSCD.prev == 1 & DIA >=5 & !is.na(PREVDIA) ~ basalArea(PREVDIA),
             TRUE ~ 0),
           BA2 = case_when(
             STATUSCD == 1 & DIA >=5 ~ basalArea(DIA),
             TRUE ~ 0))

    # mutate(#SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25,
    #        #TPAREMV_UNADJ = TPAREMV_UNADJ * REMPER,
    #        #TPAMORT_UNADJ = TPAMORT_UNADJ * REMPER,
    #        #TPARECR_UNADJ = TPARECR_UNADJ,
    #        #BA2 = basalArea(DIA) * TPAGROW_UNADJ, ## CURRENT
    #        ## Previous totals now rather than current totals
    #        #TPAGROW_UNADJ = case_when(
    #       #   COMPONENT %in% c('SURVIVOR', 'MORTALITY1', 'CUT1') ~ TPAGROW_UNADJ,
    #       #   TRUE ~ 0),
    #       # BA1 = basalArea(DIA.prev) * TPAGROW_UNADJ, ## PREVIOUS
    #        #BA_GROW = BA2- BA1,
    #        #aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.chng == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
    #        #tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0)
    #   )

  # If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
  data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
  data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
  data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
  data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
  data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
  data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

  ## Comprehensive indicator function -- w/ growth accounting
  #data$aDI_ga <- data$landD * data$aD_p * data$aD_c * data$sp * data$aChng
  #data$tDI_ga <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng
  #data$tDI_ga_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp #* data$tChng

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev
  data$tDI_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp




  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(TOTAL_TPA = sum(T_TPA * tDI, na.rm = TRUE),
                RECR_TPA = sum(R_TPA * tDI, na.rm = TRUE),
                MORT_TPA = sum(M_TPA * tDI, na.rm = TRUE),
                BAA1 = sum(BA1 * tDI, na.rm = TRUE),
                BAA2 = sum(BA2 * tDI, na.rm = TRUE),
                BAA_GROW = (BAA2 - BAA1) / BAA1,
                REMPER = first(REMPER),
                ## Demographic rates
                RECR_RATE = ((1+(RECR_TPA / TOTAL_TPA))^(1/REMPER)) - 1,
                MORT_RATE = 1 - ((1-(MORT_TPA / TOTAL_TPA))^(1/REMPER)),
                S = RECR_RATE - MORT_RATE,
                ## Basal area growth on survivors
                BAA_RATE = ((1+(BAA_GROW))^(1/REMPER)) - 1,
                nStems = length(which(tDI == 1)),
                nLive = length(which(STATUSCD == 1 & DIA >= 5))
                )

    a = NULL

  } else {
    ### Plot-level estimates -- growth accounting
    a_ga <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa_ga = sum(SUBPTYP_PROP_CHNG * aDI_ga, na.rm = TRUE),
                plotIn_ga = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))
    ### Plot-level estimates
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      left_join(select(a_ga, PLT_CN, PROP_BASIS, aGrpBy, fa_ga, plotIn_ga), by = c('PLT_CN', 'PROP_BASIS', aGrpBy))


    ### Compute total TREES in domain of interest
    t <- data %>%
      distinct(PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(PLT_CN, SUBPTYP_GRM, .dots = grpBy) %>%
      summarize(tPlot_ga = sum(TPAGROW_UNADJ * tDI_ga, na.rm = TRUE),
                rPlot_ga = sum(TPARECR_UNADJ * tDI_ga_r, na.rm = TRUE),
                mPlot_ga = sum(TPAMORT_UNADJ * tDI_ga, na.rm = TRUE),
                hPlot_ga = sum(TPAREMV_UNADJ * tDI_ga, na.rm = TRUE),
                plotIn_t_ga = ifelse(tPlot_ga >  0, 1,0),
                plotIn_r_ga = ifelse(rPlot_ga >  0, 1,0),
                plotIn_m_ga = ifelse(mPlot_ga > 0, 1,0),
                plotIn_h_ga = ifelse(hPlot_ga >  0, 1,0),
                ## No growth accoutning
                tPlot = sum(TPAGROW_UNADJ * tDI, na.rm = TRUE),
                rPlot = sum(TPARECR_UNADJ * tDI_r, na.rm = TRUE),
                mPlot = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
                hPlot = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
                plotIn_t = ifelse(tPlot >  0, 1,0),
                plotIn_r = ifelse(rPlot >  0, 1,0),
                plotIn_m = ifelse(mPlot > 0, 1,0),
                plotIn_h = ifelse(hPlot >  0, 1,0))
  }




  pltOut <- list(a = a, t = t)
  return(pltOut)

}



mcHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

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
    filter(EVAL_TYP %in% c('EXPGROW')) %>%
    #filter(EVAL_TYP %in% c('EXPCURR')) %>%
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
      fa = case_when(
        GROWTH_ACCT == 'Y' ~ fa_ga * aAdj,
        TRUE ~ fa * aAdj),
      plotIn = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_ga,
        TRUE ~ plotIn)) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = aGrpBy) %>%
    summarize(a_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              aStrat = mean(fa * a_t, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
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

  ######## ------------------ TREE ESTIMATES + CV

  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
    #filter(EVAL_TYP %in% c('EXPGROW')) %>%
    ## Need this for covariance later on
    left_join(select(a, fa, fa_ga, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
    mutate(tAdj = case_when(
      ## When NA, stay NA
      is.na(SUBPTYP_GRM) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      SUBPTYP_GRM == 0 ~ 0,
      SUBPTYP_GRM == 1 ~ as.numeric(ADJ_FACTOR_SUBP),
      SUBPTYP_GRM == 2 ~ as.numeric(ADJ_FACTOR_MICR),
      SUBPTYP_GRM == 3 ~ as.numeric(ADJ_FACTOR_MACR)),
      ## AREA
      aAdj = case_when(
        ## When NA, stay NA
        is.na(PROP_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = case_when(
        GROWTH_ACCT == 'Y' ~ fa_ga * aAdj,
        TRUE ~ fa * aAdj),
      tPlot = case_when(
        GROWTH_ACCT == 'Y' ~ tPlot_ga * tAdj,
        TRUE ~ tPlot * tAdj),
      rPlot = case_when(
        GROWTH_ACCT == 'Y' ~ rPlot_ga * tAdj,
        TRUE ~ rPlot * tAdj),
      mPlot = case_when(
        GROWTH_ACCT == 'Y' ~ mPlot_ga * tAdj,
        TRUE ~ mPlot * tAdj),
      hPlot = case_when(
        GROWTH_ACCT == 'Y' ~ hPlot_ga * tAdj,
        TRUE ~ hPlot * tAdj),
      plotIn_t = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_t_ga,
        TRUE ~ plotIn_t),
      plotIn_r = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_r_ga,
        TRUE ~ plotIn_r),
      plotIn_m = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_m_ga,
        TRUE ~ plotIn_m),
      plotIn_h = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_h_ga,
        TRUE ~ plotIn_h)) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(tPlot = sum(tPlot, na.rm = TRUE),
              rPlot = sum(rPlot, na.rm = TRUE),
              mPlot = sum(mPlot, na.rm = TRUE),
              hPlot = sum(hPlot, na.rm = TRUE),
              fa = first(fa),
              plotIn_t = ifelse(sum(plotIn_t >  0, na.rm = TRUE), 1,0),
              plotIn_r = ifelse(sum(plotIn_r >  0, na.rm = TRUE), 1,0),
              plotIn_m = ifelse(sum(plotIn_m >  0, na.rm = TRUE), 1,0),
              plotIn_h = ifelse(sum(plotIn_h >  0, na.rm = TRUE), 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              tStrat = mean(tPlot * r_t, na.rm = TRUE),
              rStrat = mean(rPlot * r_t, na.rm = TRUE),
              mStrat = mean(mPlot * r_t, na.rm = TRUE),
              hStrat = mean(hPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              plotIn_r = sum(plotIn_r, na.rm = TRUE),
              plotIn_m = sum(plotIn_m, na.rm = TRUE),
              plotIn_h = sum(plotIn_h, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              tv = stratVar(ESTN_METHOD, tPlot, tStrat, ndif, a, nh),
              rv = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh),
              mv = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh),
              hv = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh),
              # Strata level covariances
              cvStrat_r = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh, fa, aStrat),
              cvStrat_m = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh, fa, aStrat),
              cvStrat_h = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh, fa, aStrat),
              cvStrat_rp = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh, tPlot, tStrat),
              cvStrat_mp = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh, tPlot, tStrat),
              cvStrat_hp = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh, tPlot, tStrat)
    ) %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
              rEst = unitMean(ESTN_METHOD, a, nh, w, rStrat),
              mEst = unitMean(ESTN_METHOD, a, nh, w, mStrat),
              hEst = unitMean(ESTN_METHOD, a, nh, w, hStrat),
              #aEst = first(aEst),
              # Estimation of unit variance
              tVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, tv, tStrat, tEst),
              rVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, rv, rStrat, rEst),
              mVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mv, mStrat, mEst),
              hVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, hv, hStrat, hEst),
              cvEst_r = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_r, rStrat, rEst, aStrat, aEst),
              cvEst_m = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_m, mStrat, mEst, aStrat, aEst),
              cvEst_h = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_h, hStrat, hEst, aStrat, aEst),
              cvEst_rp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_rp, rStrat, rEst, tStrat, tEst),
              cvEst_mp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_mp, mStrat, mEst, tStrat, tEst),
              cvEst_hp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_hp, hStrat, hEst, tStrat, tEst),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              plotIn_r = sum(plotIn_r, na.rm = TRUE),
              plotIn_m = sum(plotIn_m, na.rm = TRUE),
              plotIn_h = sum(plotIn_h, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}


