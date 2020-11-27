bioHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
    ## Need a code that tells us where the tree was measured
    ## macroplot, microplot, subplot
    mutate(PLOT_BASIS = dplyr::case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA) ~ NA_character_,
      ## When DIA is less than 5", use microplot value
      DIA < 5 ~ 'MICR',
      ## When DIA is greater than 5", use subplot value
      DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
      DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp


  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(NETVOL_ACRE = sum(VOLCFNET * TPA_UNADJ * tDI, na.rm = TRUE),
                SAWVOL_ACRE = sum(VOLCSNET * TPA_UNADJ * tDI, na.rm = TRUE),
                BIO_AG_ACRE = sum(DRYBIO_AG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                BIO_BG_ACRE = sum(DRYBIO_BG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                CARB_AG_ACRE = sum(CARBON_AG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                CARB_BG_ACRE = sum(CARBON_BG * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                nStems = length(which(tDI == 1))) %>%
      mutate(BIO_ACRE = BIO_AG_ACRE + BIO_BG_ACRE,
             CARB_ACRE = CARB_AG_ACRE + CARB_BG_ACRE) %>%
      as.data.frame()

    a = NULL

  } else {

    grpSyms <- syms(grpBy)

    ### Plot-level estimates
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    ## Tree plts
    t <- data %>%
      lazy_dt() %>%
      filter(!is.na(PLOT_BASIS)) %>%
      group_by(PLT_CN, PLOT_BASIS, !!!grpSyms) %>%
      summarize(nvPlot = sum(VOLCFNET * TPA_UNADJ * tDI, na.rm = TRUE),
                svPlot = sum(VOLCSNET * TPA_UNADJ *  tDI, na.rm = TRUE),
                bagPlot = sum(DRYBIO_AG * TPA_UNADJ  * tDI  / 2000, na.rm = TRUE),
                bbgPlot = sum(DRYBIO_BG * TPA_UNADJ * tDI  / 2000, na.rm = TRUE),
                cagPlot = sum(CARBON_AG * TPA_UNADJ * tDI / 2000, na.rm = TRUE),
                cbgPlot = sum(CARBON_BG * TPA_UNADJ * tDI / 2000, na.rm = TRUE),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0)) %>%
      mutate(btPlot = bagPlot + bbgPlot,
             ctPlot = cagPlot + cbgPlot) %>%
      as.data.frame()
  }




  pltOut <- list(a = a, t = t)
  return(pltOut)

}



bioHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

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
        ## Otherwise, use the subplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = fa * aAdj) %>%
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
    group_by(ESTN_UNIT_CN, !!!aGrpSyms) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, av, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))



  ######## ------------------ TREE ESTIMATES + CV
  grpSyms <- syms(grpBy)
  ## Strata level estimates
  tEst <- t %>%
    lazy_dt() %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
    mutate(tAdj = dplyr::case_when(
        ## When NA, stay NA
        is.na(PLOT_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PLOT_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PLOT_BASIS == 'SUBP' ~ as.numeric(ADJ_FACTOR_SUBP),
        PLOT_BASIS == 'MICR' ~ as.numeric(ADJ_FACTOR_MICR)),
      ## AREA
      aAdj = dplyr::case_when(
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
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, !!!grpSyms) %>%
    summarize(nvPlot = sum(nvPlot, na.rm = TRUE),
              svPlot = sum(svPlot, na.rm = TRUE),
              bagPlot = sum(bagPlot, na.rm = TRUE),
              bbgPlot = sum(bbgPlot, na.rm = TRUE),
              btPlot = sum(btPlot, na.rm = TRUE),
              cagPlot = sum(cagPlot, na.rm = TRUE),
              cbgPlot = sum(cbgPlot, na.rm = TRUE),
              ctPlot = sum(ctPlot, na.rm = TRUE),
              fa = dplyr::first(fa),
              plotIn = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
              nh = dplyr::first(P2POINTCNT),
              p2eu = dplyr::first(p2eu),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU)) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(nh),
              a = dplyr::first(a),
              w = dplyr::first(w),
              p2eu = dplyr::first(p2eu),

              ## dtplyr is fast, but requires a few extra steps, so we'll finish
              ## means and variances in subsequent mutate step
              ## Strata sums
              nvStrat = sum(nvPlot, na.rm = TRUE),
              svStrat = sum(svPlot, na.rm = TRUE),
              bagStrat = sum(bagPlot, na.rm = TRUE),
              bbgStrat = sum(bbgPlot, na.rm = TRUE),
              btStrat = sum(btPlot, na.rm = TRUE),
              cagStrat = sum(cagPlot, na.rm = TRUE),
              cbgStrat = sum(cbgPlot, na.rm = TRUE),
              ctStrat = sum(ctPlot, na.rm = TRUE),
              aStrat = dplyr::first(aStrat),
              plotIn_TREE = sum(plotIn, na.rm = TRUE),

              ## Strata level variances
              nvv = sum(nvPlot^2, na.rm = TRUE),
              svv = sum(svPlot^2, na.rm = TRUE),
              bagv = sum(bagPlot^2, na.rm = TRUE),
              bbgv = sum(bbgPlot^2, na.rm = TRUE),
              btv = sum(btPlot^2, na.rm = TRUE),
              cagv = sum(cagPlot^2, na.rm = TRUE),
              cbgv = sum(cbgPlot^2, na.rm = TRUE),
              ctv = sum(ctPlot^2, na.rm = TRUE),

              ## Strata level covariances
              cvStrat_nv =sum(fa* nvPlot, na.rm = TRUE),
              cvStrat_sv =sum(fa* svPlot, na.rm = TRUE),
              cvStrat_bag =sum(fa* bagPlot, na.rm = TRUE),
              cvStrat_bbg =sum(fa* bbgPlot, na.rm = TRUE),
              cvStrat_bt =sum(fa* btPlot, na.rm = TRUE),
              cvStrat_cag =sum(fa* cagPlot, na.rm = TRUE),
              cvStrat_cbg =sum(fa* cbgPlot, na.rm = TRUE),
              cvStrat_ct =sum(fa* ctPlot, na.rm = TRUE)) %>%
    mutate(nvStrat = nvStrat / nh,
           svStrat = svStrat / nh,
           bagStrat = bagStrat / nh,
           bbgStrat = bbgStrat / nh,
           btStrat = btStrat / nh,
           cagStrat = cagStrat / nh,
           cbgStrat = cbgStrat / nh,
           ctStrat = ctStrat / nh,
           ## Variances
           adj = nh * (nh-1),
           nvv = (nvv - (nh*nvStrat^2)) / adj,
           svv = (svv - (nh*svStrat^2)) / adj,
           bagv = (bagv - (nh*bagStrat^2)) / adj,
           bbgv = (bbgv - (nh*bbgStrat^2)) / adj,
           btv = (btv - (nh*btStrat^2)) / adj,
           cagv = (cagv - (nh*cagStrat^2)) / adj,
           cbgv = (cbgv - (nh*cbgStrat^2)) / adj,
           ctv = (ctv - (nh*ctStrat^2)) / adj,
           ## Covariances
           cvStrat_nv = (cvStrat_nv - (nh * nvStrat * aStrat)) / adj,
           cvStrat_sv = (cvStrat_sv - (nh * svStrat * aStrat)) / adj,
           cvStrat_bag =  (cvStrat_bag - (nh * bagStrat * aStrat)) / adj,
           cvStrat_bbg = (cvStrat_bbg - (nh * bbgStrat * aStrat)) / adj,
           cvStrat_bt = (cvStrat_bt - (nh * btStrat * aStrat)) / adj,
           cvStrat_cag = (cvStrat_cag - (nh * cagStrat * aStrat)) / adj,
           cvStrat_cbg = (cvStrat_cbg - (nh * cbgStrat * aStrat)) / adj,
           cvStrat_ct = (cvStrat_ct - (nh * ctStrat * aStrat)) / adj) %>%
    as.data.frame() %>%

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
              N = dplyr::first(p2eu),
              # Estimation of unit variance
              nvVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, nvv, nvStrat, nvEst),
              svVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, svv, svStrat, svEst),
              bagVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bagv, bagStrat, bagEst),
              bbgVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bbgv, bbgStrat, bbgEst),
              btVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, btv, btStrat, btEst),
              cagVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cagv, cagStrat, cagEst),
              cbgVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cbgv, cbgStrat, cbgEst),
              ctVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, ctv, ctStrat, ctEst),
              cvEst_nv = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_nv, nvStrat, nvEst, aStrat, aEst),
              cvEst_sv = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_sv, svStrat, svEst, aStrat, aEst),
              cvEst_bag = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_bag, bagStrat, bagEst, aStrat, aEst),
              cvEst_bbg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_bbg, bbgStrat, bbgEst, aStrat, aEst),
              cvEst_bt = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_bt, btStrat, btEst, aStrat, aEst),
              cvEst_cag = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_cag, cagStrat, cagEst, aStrat, aEst),
              cvEst_cbg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_cbg, cbgStrat, cbgEst, aStrat, aEst),
              cvEst_ct = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_ct, ctStrat, ctEst, aStrat, aEst),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}







biomassHelper <- function(x, combos, data, grpBy, aGrpBy, totals, SE){
  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  ad = 1

  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(tObs))) == length(tObs)) tObs <- 1
    td <- data$tDI * tObs * td
    # Area domain indicator for each column in
    if(grpBy[n] %in% aGrpBy){
      aObs <- as.character(combos[[x]][[aGrpBy[n]]]) == as.character(data[[aGrpBy[n]]])
      if (length(which(is.na(aObs))) == length(aObs)) aObs <- 1
      aObs[is.na(aObs)] <- 0
      ad <- data$aDI * aObs * ad
    }
  }

  if(SE){
    data$tDI <- td
    data$tDI[is.na(data$tDI)] <- 0
    data$aDI <- ad
    data$aDI[is.na(data$aDI)] <- 0
    ## We produce an intermediate object in this chain as it is needed to compute the ratio of means variance
    ## Numerator and denominator are in different domains of interest, and may be grouped by different variables
    ## see covariance estimation below

    ### Compute total TREES in domain of interest
    bInt <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      #filter(EVALID %in% tID) %>%
      # Compute estimates at plot level
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(nvPlot = sum(VOLCFNET * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                svPlot = sum(VOLCSNET * TPA_UNADJ * tAdj * tDI, na.rm = TRUE),
                bagPlot = sum(DRYBIO_AG * TPA_UNADJ * tAdj * tDI  / 2000, na.rm = TRUE),
                bbgPlot = sum(DRYBIO_BG * TPA_UNADJ * tAdj * tDI  / 2000, na.rm = TRUE),
                btPlot = sum(bagPlot, bbgPlot, na.rm = TRUE),
                cagPlot = sum(CARBON_AG * TPA_UNADJ * tAdj * tDI / 2000, na.rm = TRUE),
                cbgPlot = sum(CARBON_BG * TPA_UNADJ * tAdj * tDI / 2000, na.rm = TRUE),
                ctPlot = sum(cagPlot, cbgPlot, na.rm = TRUE),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                a = dplyr::first(AREA_USED),
                p1EU = dplyr::first(P1PNTCNT_EU),
                p1 = dplyr::first(P1POINTCNT),
                p2 = dplyr::first(P2POINTCNT))
    # Continue through totals
    b <- bInt %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(nvStrat = mean(nvPlot, na.rm = TRUE),
                svStrat = mean(svPlot, na.rm = TRUE),
                bagStrat = mean(bagPlot, na.rm = TRUE),
                bbgStrat = mean(bbgPlot, na.rm = TRUE),
                btStrat = mean(btPlot, na.rm = TRUE),
                cagStrat = mean(cagPlot, na.rm = TRUE),
                cbgStrat = mean(cbgPlot, na.rm = TRUE),
                ctStrat = mean(ctPlot, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = dplyr::first(a),
                w = dplyr::first(p1) / dplyr::first(p1EU), # Stratum weight
                nh = dplyr::first(p2), # Number plots in stratum
                # Strata level variances
                nvv = ifelse(first(ESTN_METHOD == 'simple'),
                             var(nvPlot * dplyr::first(a) / nh),
                             (sum(nvPlot^2) - sum(nh * nvStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                svv = ifelse(first(ESTN_METHOD == 'simple'),
                             var(svPlot * dplyr::first(a) / nh),
                             (sum(svPlot^2) - sum(nh * svStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bagv = ifelse(first(ESTN_METHOD == 'simple'),
                              var(bagPlot * dplyr::first(a) / nh),
                              (sum(bagPlot^2) - sum(nh * bagStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bbgv = ifelse(first(ESTN_METHOD == 'simple'),
                              var(bbgPlot * dplyr::first(a) / nh),
                              (sum(bbgPlot^2) - sum(nh * bbgStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                btv = ifelse(first(ESTN_METHOD == 'simple'),
                             var(btPlot * dplyr::first(a) / nh),
                             (sum(btPlot^2) - sum(nh * btStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cagv = ifelse(first(ESTN_METHOD == 'simple'),
                              var(cagPlot * dplyr::first(a) / nh),
                              (sum(cagPlot^2) - sum(nh * cagStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cbgv = ifelse(first(ESTN_METHOD == 'simple'),
                              var(cbgPlot * dplyr::first(a) / nh),
                              (sum(cbgPlot^2) - sum(nh * cbgStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                ctv = ifelse(first(ESTN_METHOD == 'simple'),
                             var(ctPlot * dplyr::first(a) / nh),
                             (sum(ctPlot^2) - sum(nh * ctStrat^2)) / (nh * (nh-1)))) %>% # Stratified and double cases
      group_by(ESTN_UNIT_CN) %>%
      summarize(nvEst = unitMean(ESTN_METHOD, a, nh, w, nvStrat),
                svEst = unitMean(ESTN_METHOD, a, nh, w, svStrat),
                bagEst = unitMean(ESTN_METHOD, a, nh, w, bagStrat),
                bbgEst = unitMean(ESTN_METHOD, a, nh, w, bbgStrat),
                btEst = unitMean(ESTN_METHOD, a, nh, w, btStrat),
                cagEst = unitMean(ESTN_METHOD, a, nh, w, cagStrat),
                cbgEst = unitMean(ESTN_METHOD, a, nh, w, cbgStrat),
                ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
                plotIn = sum(plotIn, na.rm = TRUE),
                # Estimation of unit variance
                nvVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, nvv, nvStrat, nvEst),
                svVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, svv, svStrat, svEst),
                bagVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bagv, bagStrat, bagEst),
                bbgVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bbgv, bbgStrat, bbgEst),
                btVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, btv, btStrat, btEst),
                cagVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, cagv, cagStrat, cagEst),
                cbgVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, cbgv, cbgStrat, cbgEst),
                ctVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, ctv, ctStrat, ctEst)) %>%
      # Compute totals
      summarize(NETVOL_TOTAL = sum(nvEst, na.rm = TRUE),
                SAWVOL_TOTAL = sum(svEst, na.rm = TRUE),
                BIO_AG_TOTAL = sum(bagEst, na.rm = TRUE),
                BIO_BG_TOTAL = sum(bbgEst, na.rm = TRUE),
                BIO_TOTAL = sum(btEst, na.rm = TRUE),
                CARB_AG_TOTAL = sum(cagEst, na.rm = TRUE),
                CARB_BG_TOTAL = sum(cbgEst, na.rm = TRUE),
                CARB_TOTAL = sum(ctEst, na.rm = TRUE),
                nvVar = sum(nvVar, na.rm = TRUE),
                NETVOL_TOT_SE = sqrt(nvVar) / NETVOL_TOTAL * 100,
                svVar = sum(svVar, na.rm = TRUE),
                SAWVOL_TOT_SE = sqrt(svVar) / SAWVOL_TOTAL * 100,
                bagVar = sum(bagVar, na.rm = TRUE),
                BIO_AG_TOT_SE = sqrt(bagVar) / BIO_AG_TOTAL * 100,
                bbgVar = sum(bbgVar, na.rm = TRUE),
                BIO_BG_TOT_SE = sqrt(bbgVar) / BIO_BG_TOTAL * 100,
                btVar = sum(btVar, na.rm = TRUE),
                BIO_TOT_SE = sqrt(btVar) / BIO_TOTAL * 100,
                cagVar = sum(cagVar, na.rm = TRUE),
                CARB_AG_TOT_SE = sqrt(cagVar) / CARB_AG_TOTAL * 100,
                cbgVar = sum(cbgVar, na.rm = TRUE),
                CARB_BG_TOT_SE = sqrt(cbgVar) / CARB_BG_TOTAL * 100,
                ctVar = sum(ctVar, na.rm = TRUE),
                CARB_TOT_SE = sqrt(ctVar) / CARB_TOTAL * 100,
                nPlots_VOL = sum(plotIn, na.rm = TRUE))


    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      #filter(EVALID %in% aID) %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
                a = dplyr::first(AREA_USED),
                p1EU = dplyr::first(P1PNTCNT_EU),
                p1 = dplyr::first(P1POINTCNT),
                p2 = dplyr::first(P2POINTCNT)) #%>%
    #distinct(PLT_CN, .keep_all = TRUE)
    a <- aInt %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = dplyr::first(a),
                w = dplyr::first(p1) / dplyr::first(p1EU), # Stratum weight
                nh = dplyr::first(p2), # Number plots in stratum
                # Strata level variance
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * dplyr::first(a) / nh),
                            (sum(fa^2) - sum(nh * aStrat^2)) / (nh * (nh-1)))) %>% # Stratified and double cases
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                plotIn = sum(plotIn, na.rm = TRUE),
                # Estimation unit variance
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst)) %>%
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                areaVar = sum(aVar, na.rm = TRUE),
                AREA_TOTAL_SE = sqrt(areaVar) / AREA_TOTAL * 100,
                nPlots_AREA = sum(plotIn, na.rm = TRUE))


    ## Compute COVARIANCE between numerator and denominator (for ratio estimates of variance)
    covar <- bInt %>%
      inner_join(aInt, by = c('PLT_CN', 'ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN'), suffix = c('_b', '_a')) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                nvStrat = mean(nvPlot, na.rm = TRUE),
                svStrat = mean(svPlot, na.rm = TRUE),
                bagStrat = mean(bagPlot, na.rm = TRUE),
                bbgStrat = mean(bbgPlot, na.rm = TRUE),
                btStrat = mean(btPlot, na.rm = TRUE),
                cagStrat = mean(cagPlot, na.rm = TRUE),
                cbgStrat = mean(cbgPlot, na.rm = TRUE),
                ctStrat = mean(ctPlot, na.rm = TRUE),
                a = dplyr::first(a_b),
                w = dplyr::first(p1_b) / dplyr::first(p1EU_a), # Stratum weight
                nh = dplyr::first(p2_b), # Number plots in stratum
                # Strata level covariances
                cvStrat_nv = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(fa,nvPlot),
                                    (sum(fa*nvPlot) - sum(nh * aStrat *nvStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_sv = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(fa,svPlot),
                                    (sum(fa*svPlot) - sum(nh * aStrat *svStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_bag = ifelse(first(ESTN_METHOD == 'simple'),
                                     cov(fa,bagPlot),
                                     (sum(fa*bagPlot) - sum(nh * aStrat *bagStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_bbg = ifelse(first(ESTN_METHOD == 'simple'),
                                     cov(fa,bbgPlot),
                                     (sum(fa*bbgPlot) - sum(nh * aStrat *bbgStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_bt = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(fa,btPlot),
                                    (sum(fa*btPlot) - sum(nh * aStrat *btStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_cag = ifelse(first(ESTN_METHOD == 'simple'),
                                     cov(fa,cagPlot),
                                     (sum(fa*cagPlot) - sum(nh * aStrat *cagStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_cbg = ifelse(first(ESTN_METHOD == 'simple'),
                                     cov(fa,cbgPlot),
                                     (sum(fa*cbgPlot) - sum(nh * aStrat *cbgStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_ct = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(fa,ctPlot),
                                    (sum(fa*ctPlot) - sum(nh * aStrat *ctStrat)) / (nh * (nh-1)))) %>% # Stratified and double cases)
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                nvEst = unitMean(ESTN_METHOD, a, nh, w, nvStrat),
                svEst = unitMean(ESTN_METHOD, a, nh, w, svStrat),
                bagEst = unitMean(ESTN_METHOD, a, nh, w, bagStrat),
                bbgEst = unitMean(ESTN_METHOD, a, nh, w, bbgStrat),
                btEst = unitMean(ESTN_METHOD, a, nh, w, btStrat),
                cagEst = unitMean(ESTN_METHOD, a, nh, w, cagStrat),
                cbgEst = unitMean(ESTN_METHOD, a, nh, w, cbgStrat),
                ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
                cvEst_nv = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_nv, nvStrat, nvEst, aStrat, aEst),
                cvEst_sv = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_sv, svStrat, svEst, aStrat, aEst),
                cvEst_bag = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_bag, bagStrat, bagEst, aStrat, aEst),
                cvEst_bbg = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_bbg, bbgStrat, bbgEst, aStrat, aEst),
                cvEst_bt = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_bt, btStrat, btEst, aStrat, aEst),
                cvEst_cag = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_cag, cagStrat, cagEst, aStrat, aEst),
                cvEst_cbg = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_cbg, cbgStrat, cbgEst, aStrat, aEst),
                cvEst_ct = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_ct, ctStrat, ctEst, aStrat, aEst)) %>%
      summarize(nvCV = sum(cvEst_nv, na.rm = TRUE),
                svCV = sum(cvEst_sv, na.rm = TRUE),
                bagCV = sum(cvEst_bag, na.rm = TRUE),
                bbgCV = sum(cvEst_bbg, na.rm = TRUE),
                btCV = sum(cvEst_bt, na.rm = TRUE),
                cagCV = sum(cvEst_cag, na.rm = TRUE),
                cbgCV = sum(cvEst_cbg, na.rm = TRUE),
                ctCV = sum(cvEst_ct, na.rm = TRUE))

    ## Join up tree, area, and covariance tables
    b <- data.frame(b, a, covar)

    ## Make TPA, BAA and SE
    b <- b %>%
      mutate(NETVOL_ACRE = NETVOL_TOTAL / AREA_TOTAL,
             SAWVOL_ACRE = SAWVOL_TOTAL / AREA_TOTAL,
             BIO_AG_ACRE = BIO_AG_TOTAL / AREA_TOTAL,
             BIO_BG_ACRE = BIO_BG_TOTAL / AREA_TOTAL,
             BIO_ACRE = BIO_TOTAL / AREA_TOTAL,
             CARB_AG_ACRE = CARB_AG_TOTAL / AREA_TOTAL,
             CARB_BG_ACRE = CARB_BG_TOTAL / AREA_TOTAL,
             CARB_ACRE = CARB_TOTAL / AREA_TOTAL,
             nvaVar = (1/AREA_TOTAL^2) * (nvVar + (NETVOL_ACRE^2 * areaVar) - 2 * NETVOL_ACRE * nvCV),
             svaVar = (1/AREA_TOTAL^2) * (svVar + (SAWVOL_ACRE^2 * areaVar) - 2 * SAWVOL_ACRE * svCV),
             bagaVar = (1/AREA_TOTAL^2) * (bagVar + (BIO_AG_ACRE^2 * areaVar) - 2 * BIO_AG_ACRE * bagCV),
             bbgaVar = (1/AREA_TOTAL^2) * (bbgVar + (BIO_BG_ACRE^2 * areaVar) - 2 * BIO_BG_ACRE * bbgCV),
             btaVar = (1/AREA_TOTAL^2) * (btVar + (BIO_ACRE^2 * areaVar) - 2 * BIO_ACRE * btCV),
             cagaVar = (1/AREA_TOTAL^2) * (cagVar + (CARB_AG_ACRE^2 * areaVar) - 2 * CARB_AG_ACRE * cagCV),
             cbgaVar = (1/AREA_TOTAL^2) * (cbgVar + (CARB_BG_ACRE^2 * areaVar) - 2 * CARB_BG_ACRE * cbgCV),
             ctaVar = (1/AREA_TOTAL^2) * (ctVar + (CARB_ACRE^2 * areaVar) - 2 * CARB_ACRE * ctCV),
             NETVOL_ACRE_SE = sqrt(nvaVar) / NETVOL_ACRE * 100,
             SAWVOL_ACRE_SE = sqrt(svaVar) / SAWVOL_ACRE * 100,
             BIO_AG_ACRE_SE = sqrt(bagaVar) / BIO_AG_ACRE * 100,
             BIO_BG_ACRE_SE = sqrt(bbgaVar) / BIO_BG_ACRE * 100,
             BIO_ACRE_SE = sqrt(btaVar) / BIO_ACRE * 100,
             CARB_AG_ACRE_SE = sqrt(cagaVar) / CARB_AG_ACRE * 100,
             CARB_BG_ACRE_SE = sqrt(cbgaVar) / CARB_BG_ACRE * 100,
             CARB_ACRE_SE = sqrt(ctaVar) / CARB_ACRE * 100)

    if (totals) {
      b <- b %>%
        select(NETVOL_ACRE, SAWVOL_ACRE,  BIO_AG_ACRE,  BIO_BG_ACRE, BIO_ACRE, CARB_AG_ACRE, CARB_BG_ACRE, CARB_ACRE,
               NETVOL_ACRE_SE, SAWVOL_ACRE_SE,  BIO_AG_ACRE_SE,  BIO_BG_ACRE_SE, BIO_ACRE_SE, CARB_AG_ACRE_SE, CARB_BG_ACRE_SE, CARB_ACRE_SE,
               NETVOL_TOTAL, SAWVOL_TOTAL,  BIO_AG_TOTAL,  BIO_BG_TOTAL, BIO_TOTAL, CARB_AG_TOTAL, CARB_BG_TOTAL, CARB_TOTAL, AREA_TOTAL,
               nPlots_VOL, nPlots_AREA)
    } else {
      b <- b %>%
        select(names(b)[str_detect(names(b), 'Var', negate = TRUE) & str_detect(names(b), 'ACRE')], nPlots_VOL, nPlots_AREA)
    }
    # Rejoin with some grpBy Names
    b <- data.frame(combos[[x]], b)

  } else { # No sampling errors
    ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
    b <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(nvPlot = sum(VOLCFNET * TPA_UNADJ * tAdj * tDI * EXPNS, na.rm = TRUE),
                svPlot = sum(VOLCSNET * TPA_UNADJ * tAdj * tDI * EXPNS, na.rm = TRUE),
                bagPlot = sum(DRYBIO_AG * TPA_UNADJ * tAdj * tDI * EXPNS  / 2000, na.rm = TRUE),
                bbgPlot = sum(DRYBIO_BG * TPA_UNADJ * tAdj * tDI * EXPNS  / 2000, na.rm = TRUE),
                cagPlot = sum(CARBON_AG * TPA_UNADJ * tAdj * tDI * EXPNS  / 2000, na.rm = TRUE),
                cbgPlot = sum(CARBON_BG * TPA_UNADJ * tAdj * tDI * EXPNS  / 2000, na.rm = TRUE),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0)) %>%
      group_by(.dots = grpBy) %>%
      summarize(NETVOL_TOTAL = sum(nvPlot, na.rm = TRUE),
                SAWVOL_TOTAL = sum(svPlot, na.rm = TRUE),
                BIO_AG_TOTAL = sum(bagPlot, na.rm = TRUE),
                BIO_BG_TOTAL = sum(bbgPlot, na.rm = TRUE),
                BIO_TOTAL = sum(BIO_AG_TOTAL, BIO_BG_TOTAL, na.rm = TRUE),
                CARB_AG_TOTAL = sum(cagPlot, na.rm = TRUE),
                CARB_BG_TOTAL = sum(cbgPlot, na.rm = TRUE),
                CARB_TOTAL = sum(CARB_AG_TOTAL, CARB_BG_TOTAL, na.rm = TRUE),
                nPlots_VOL = sum(plotIn, na.rm = TRUE))

    a <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      group_by(.dots = aGrpBy) %>%
      summarize(AREA_TOTAL = sum(fa, na.rm = TRUE),
                nPlots_AREA = sum(plotIn, na.rm = TRUE))

    suppressMessages({
      b <- inner_join(b, a) %>%
        mutate(NETVOL_ACRE = NETVOL_TOTAL / AREA_TOTAL,
               SAWVOL_ACRE = SAWVOL_TOTAL / AREA_TOTAL,
               BIO_AG_ACRE = BIO_AG_TOTAL / AREA_TOTAL,
               BIO_BG_ACRE = BIO_BG_TOTAL / AREA_TOTAL,
               BIO_ACRE = BIO_TOTAL / AREA_TOTAL,
               CARB_AG_ACRE = CARB_AG_TOTAL / AREA_TOTAL,
               CARB_BG_ACRE = CARB_BG_TOTAL / AREA_TOTAL,
               CARB_ACRE = CARB_TOTAL / AREA_TOTAL)

      if (totals) {
        b <- b %>%
          select(grpBy, names(b)[str_detect(names(bOut), 'Var', negate = TRUE)], nPlots_VOL, nPlots_AREA)
      } else {
        b <- b %>%
          select(grpBy, names(b)[str_detect(names(b), 'Var', negate = TRUE) & str_detect(names(b), 'ACRE')], nPlots_VOL, nPlots_AREA)
      }
    })


  } # End SE Conditional
  #gc()

  # Return t
  b
}



