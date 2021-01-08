gmHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]


  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
      left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'state_recr', TPA_UNADJ, STATUSCD, DIA)), by = c('PLT_CN', 'CONDID')) %>%
      left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPARECR_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
      left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA', 'state')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
      left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
      left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
      left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
    mutate_if(is.factor,
              as.character) %>%
    mutate(TPAGROW_UNADJ = TPAGROW_UNADJ * state,
           TPAREMV_UNADJ = TPAREMV_UNADJ * state,
           TPAMORT_UNADJ = TPAMORT_UNADJ * state,
           TPARECR_UNADJ = TPARECR_UNADJ * state_recr / REMPER,
           ## State recruit is the state variable adjustment for ALL TREES at T2,
           ## So we can estimate live TPA at t2 (t1 unavailable w/out growth accounting) with:
           TPA_UNADJ = TPA_UNADJ * state_recr * if_else(STATUSCD == 1 & DIA >= 5, 1, 0),
           aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
           tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0),
           test = if_else(COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2'), 1, 0))


  # If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
  data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
  data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
  data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
  ## Issue is here for pre-growth accounting
  #data$landD.prev <- ifelse(data$landD == 1 & data$landD.prev == 1, 1, 0)
  data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
  data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
  data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

  ## Comprehensive indicator function -- w/ growth accounting
  data$tDI_ga <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng
  data$tDI_ga_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp * data$tChng

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev
  data$tDI_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp

  ### DOING AREA SEPARATELY NOW FOR GROWTH ACCOUNTING PLOTS
  aData <- select(db$PLOT,c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$SUBP_COND_CHNG_MTRX, PLT_CN, PREV_PLT_CN, SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID), by = c('PLT_CN', 'PREV_PLT_CN')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$COND, c('PLT_CN', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
    mutate(aChng = if_else(COND_STATUS_CD == 1 &
                             COND_STATUS_CD.prev == 1 &
                             !is.null(CONDPROP_UNADJ) &
                             SUBPTYP == 1,
                           1, 0),
           SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25)

  aData$landD <- ifelse(aData$landD == 1 & aData$landD.prev == 1, 1, 0)
  aData$aDI_ga <- aData$landD * aData$aD_p * aData$aD_c * aData$sp * aData$aChng


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)
    t <- data %>%
      lazy_dt() %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(RECR_TPA = sum(TPARECR_UNADJ * tDI, na.rm = TRUE),
                MORT_TPA = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
                REMV_TPA = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
                TOTAL_TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
                nStems = length(which(tDI == 1))) %>%
      mutate(RECR_PERC = RECR_TPA / TOTAL_TPA * 100,
             MORT_PERC = MORT_TPA / TOTAL_TPA * 100,
             REMV_PERC = REMV_TPA / TOTAL_TPA * 100) %>%
      as.data.frame()

    a = NULL

  } else {
    ### Plot-level estimates -- growth accounting
    a_ga <- aData %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      #distinct(PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa_ga = sum(SUBPTYP_PROP_CHNG * aDI_ga, na.rm = TRUE),
                plotIn_ga = ifelse(sum(aDI_ga >  0, na.rm = TRUE), 1,0))
    ### Plot-level estimates
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      left_join(select(a_ga, PLT_CN, PROP_BASIS, aGrpBy, fa_ga, plotIn_ga), by = c('PLT_CN', 'PROP_BASIS', aGrpBy))

    grpSyms <- syms(grpBy)

    ### Compute total TREES in domain of interest
    t <- data %>%
      distinct(PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
      lazy_dt() %>%
      # Compute estimates at plot level
      group_by(PLT_CN, SUBPTYP_GRM, !!! grpSyms) %>%
      summarize(rPlot_ga = sum(TPARECR_UNADJ * tDI_ga_r, na.rm = TRUE),
                mPlot_ga = sum(TPAMORT_UNADJ * tDI_ga, na.rm = TRUE),
                hPlot_ga = sum(TPAREMV_UNADJ * tDI_ga, na.rm = TRUE),
                tPlot_ga = sum(TPA_UNADJ * tDI_ga, na.rm = TRUE),
                ## No growth accoutning
                rPlot = sum(TPARECR_UNADJ * tDI_r, na.rm = TRUE),
                mPlot = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
                hPlot = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
                tPlot = sum(TPA_UNADJ * tDI, na.rm = TRUE)) %>%
      mutate(plotIn_t_ga = ifelse(tPlot_ga >  0, 1,0),
             plotIn_r_ga = ifelse(rPlot_ga >  0, 1,0),
             plotIn_m_ga = ifelse(mPlot_ga > 0, 1,0),
             plotIn_h_ga = ifelse(hPlot_ga >  0, 1,0),
             plotIn_t = ifelse(tPlot >  0, 1,0),
             plotIn_r = ifelse(rPlot >  0, 1,0),
             plotIn_m = ifelse(mPlot > 0, 1,0),
             plotIn_h = ifelse(hPlot >  0, 1,0)) %>%
      as.data.frame()
  }




  pltOut <- list(a = a, t = t)
  return(pltOut)

}



gmHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    aGrpBy <- c(aGrpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$P1POINTCNT <- popState[[x]]$P1POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  aGrpSyms <- syms(aGrpBy)

  ## Strata level estimates
  aStrat <- a %>%
    lazy_dt() %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    mutate(
      ## AREA
      aAdj = dplyr::case_when(
        ## When NA, stay NA
        is.na(PROP_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ fa_ga * aAdj,
        TRUE ~ fa * aAdj),
      plotIn = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_ga,
        TRUE ~ plotIn)) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!aGrpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              aStrat = sum(fa, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              p2eu = dplyr::first(p2eu),
              av = sum(fa^2, na.rm = TRUE)) %>%
    mutate(aStrat = aStrat / nh,
           av = (av - (nh * aStrat^2)) / (nh * (nh-1))) %>%
    as.data.frame()
  ## Estimation unit
  aEst <- aStrat %>%
    group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, av, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

  grpSyms <- syms(grpBy)

  ######## ------------------ TREE ESTIMATES + CV

  ## Strata level estimates
  tEst <- t %>%
    lazy_dt() %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
    #filter(EVAL_TYP %in% c('EXPGROW')) %>%
    ## Need this for covariance later on
    left_join(select(a, fa, fa_ga, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
    mutate(tAdj = dplyr::case_when(
      ## When NA, stay NA
      is.na(SUBPTYP_GRM) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      SUBPTYP_GRM == 0 ~ 0,
      SUBPTYP_GRM == 1 ~ as.numeric(ADJ_FACTOR_SUBP),
      SUBPTYP_GRM == 2 ~ as.numeric(ADJ_FACTOR_MICR),
      SUBPTYP_GRM == 3 ~ as.numeric(ADJ_FACTOR_MACR)),
      ## AREA
      aAdj = dplyr::case_when(
        ## When NA, stay NA
        is.na(PROP_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ fa_ga * aAdj,
        TRUE ~ fa * aAdj),
      tPlot = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ tPlot_ga * tAdj,
        TRUE ~ tPlot * tAdj),
      rPlot = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ rPlot_ga * tAdj,
        TRUE ~ rPlot * tAdj),
      mPlot = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ mPlot_ga * tAdj,
        TRUE ~ mPlot * tAdj),
      hPlot = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ hPlot_ga * tAdj,
        TRUE ~ hPlot * tAdj),
      plotIn_t = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_t_ga,
        TRUE ~ plotIn_t),
      plotIn_r = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_r_ga,
        TRUE ~ plotIn_r),
      plotIn_m = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_m_ga,
        TRUE ~ plotIn_m),
      plotIn_h = dplyr::case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_h_ga,
        TRUE ~ plotIn_h)) %>%
    filter(!is.na(SUBPTYP_GRM)) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, !!!grpSyms) %>%
    summarize(tPlot = sum(tPlot, na.rm = TRUE),
              rPlot = sum(rPlot, na.rm = TRUE),
              mPlot = sum(mPlot, na.rm = TRUE),
              hPlot = sum(hPlot, na.rm = TRUE),
              aZero = sum(fa, na.rm = TRUE),
              fa = dplyr::first(fa),
              plotIn_t = ifelse(sum(plotIn_t >  0, na.rm = TRUE), 1,0),
              plotIn_r = ifelse(sum(plotIn_r >  0, na.rm = TRUE), 1,0),
              plotIn_m = ifelse(sum(plotIn_m >  0, na.rm = TRUE), 1,0),
              plotIn_h = ifelse(sum(plotIn_h >  0, na.rm = TRUE), 1,0),
              nh = dplyr::first(P2POINTCNT),
              p2eu = dplyr::first(p2eu),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU)) %>%
    # If area is 0, so is numerator
    mutate(aZero = ifelse(aZero > 0, 1, 0), # Binary, zero if 0 and 1 otherwise
           tPlot = tPlot * aZero,
           rPlot = rPlot * aZero,
           mPlot = mPlot * aZero,
           hPlot = hPlot * aZero,
           plotIn_t = plotIn_t * aZero,
           plotIn_r = plotIn_r * aZero,
           plotIn_m = plotIn_m * aZero,
           plotIn_h = plotIn_h * aZero) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(nh),
              a = dplyr::first(a),
              w = dplyr::first(w),
              p2eu = dplyr::first(p2eu),

              ## dtplyr is fast, but requires a few extra steps, so we'll finish
              ## means and variances in subseqent mutate step

              ## Strata means
              tStrat = sum(tPlot, na.rm = TRUE),
              rStrat = sum(rPlot, na.rm = TRUE),
              mStrat = sum(mPlot, na.rm = TRUE),
              hStrat = sum(hPlot, na.rm = TRUE),
              aStrat = dplyr::first(aStrat),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              plotIn_r = sum(plotIn_r, na.rm = TRUE),
              plotIn_m = sum(plotIn_m, na.rm = TRUE),
              plotIn_h = sum(plotIn_h, na.rm = TRUE),
              # ## Strata level variances
              tv = sum(tPlot^2, na.rm = TRUE),
              rv = sum(rPlot^2, na.rm = TRUE),
              mv = sum(mPlot^2, na.rm = TRUE),
              hv = sum(hPlot^2, na.rm = TRUE),
              # Strata level covariances
              cvStrat_r = sum(rPlot*fa, na.rm = TRUE),
              cvStrat_m = sum(mPlot*fa, na.rm = TRUE),
              cvStrat_h = sum(hPlot*fa, na.rm = TRUE),
              cvStrat_rp = sum(rPlot*tPlot, na.rm = TRUE),
              cvStrat_mp = sum(mPlot*tPlot, na.rm = TRUE),
              cvStrat_hp = sum(hPlot*tPlot, na.rm = TRUE)) %>%
    mutate(tStrat = tStrat / nh,
           rStrat = rStrat / nh,
           mStrat = mStrat / nh,
           hStrat = hStrat / nh,
           adj = nh * (nh-1),
           tv = (tv - (nh*tStrat^2)) / adj,
           rv = (rv - (nh*rStrat^2)) / adj,
           mv = (mv - (nh*mStrat^2)) / adj,
           hv = (hv - (nh*hStrat^2)) / adj,
           cvStrat_r = (cvStrat_r - (nh * rStrat * aStrat)) / adj,
           cvStrat_m = (cvStrat_m - (nh * mStrat * aStrat)) / adj,
           cvStrat_h = (cvStrat_h - (nh * hStrat * aStrat)) / adj,
           cvStrat_rp = (cvStrat_rp - (nh * rStrat * tStrat)) / adj,
           cvStrat_mp = (cvStrat_mp - (nh * mStrat * tStrat)) / adj,
           cvStrat_hp = (cvStrat_hp - (nh * hStrat * tStrat)) / adj) %>%
    as.data.frame() %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
              rEst = unitMean(ESTN_METHOD, a, nh, w, rStrat),
              mEst = unitMean(ESTN_METHOD, a, nh, w, mStrat),
              hEst = unitMean(ESTN_METHOD, a, nh, w, hStrat),
              N = dplyr::first(p2eu),
              #aEst = dplyr::first(aEst),
              # Estimation of unit variance
              tVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, tv, tStrat, tEst),
              rVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, rv, rStrat, rEst),
              mVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, mv, mStrat, mEst),
              hVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, hv, hStrat, hEst),
              cvEst_r = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_r, rStrat, rEst, aStrat, aEst),
              cvEst_m = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_m, mStrat, mEst, aStrat, aEst),
              cvEst_h = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_h, hStrat, hEst, aStrat, aEst),
              cvEst_rp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_rp, rStrat, rEst, tStrat, tEst),
              cvEst_mp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_mp, mStrat, mEst, tStrat, tEst),
              cvEst_hp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_hp, hStrat, hEst, tStrat, tEst),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              plotIn_r = sum(plotIn_r, na.rm = TRUE),
              plotIn_m = sum(plotIn_m, na.rm = TRUE),
              plotIn_h = sum(plotIn_h, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}





