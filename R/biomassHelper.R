bioHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot, component){

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


  ## Convert to long format, where biomass component is the observation (multiple per tree)
  data <- data %>%
    pivot_longer(cols = DRYBIO_TOP:DRYBIO_FOLIAGE,
                 names_to = c(".value", 'COMPONENT'),
                 names_sep = 7) %>%
    rename(DRYBIO = DRYBIO_) %>%
    filter(COMPONENT %in% component) #%>%
  #filter(!is.na(DRYBIO))


  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, COMPONENT, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(BIO_ACRE = sum(DRYBIO * TPA_UNADJ * tDI, na.rm = TRUE) / 2000,
                nStems = length(which(tDI == 1))) %>%
      mutate(CARB_ACRE = BIO_ACRE * .5) %>%
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
      summarize(bPlot = sum(DRYBIO * TPA_UNADJ  * tDI / 2000, na.rm = TRUE),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0)) %>%
      mutate(cPlot = bPlot * .5) %>%
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
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE),
              N = dplyr::first(p2eu),
              A = dplyr::first(a))



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
      bPlot = bPlot* tAdj,
      cPlot = cPlot* tAdj) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, !!!grpSyms) %>%
    summarize(bPlot = sum(bPlot, na.rm = TRUE),
              cPlot = sum(cPlot, na.rm = TRUE),
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
              bStrat = sum(bPlot, na.rm = TRUE),
              cStrat = sum(cPlot, na.rm = TRUE),
              aStrat = dplyr::first(aStrat),
              plotIn_TREE = sum(plotIn, na.rm = TRUE),

              ## Strata level variances
              bv = sum(bPlot^2, na.rm = TRUE),
              cv = sum(cPlot^2, na.rm = TRUE),

              ## Strata level covariances
              cvStrat_b =sum(fa* bPlot, na.rm = TRUE),
              cvStrat_c =sum(fa* cPlot, na.rm = TRUE)) %>%
    mutate(bStrat = bStrat / nh,
           cStrat = cStrat / nh,

           ## Variances
           adj = nh * (nh-1),
           bv = (bv - (nh*bStrat^2)) / adj,
           cv = (cv - (nh*cStrat^2)) / adj,

           ## Covariances
           cvStrat_b = (cvStrat_b - (nh * bStrat * aStrat)) / adj,
           cvStrat_c = (cvStrat_c - (nh * cStrat * aStrat)) / adj) %>%
    as.data.frame() %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
              cEst = unitMean(ESTN_METHOD, a, nh, w, cStrat),
              # Estimation of unit variance
              bVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bv, bStrat, bEst),
              cVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cv, cStrat, cEst),
              cvEst_b = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_b, bStrat, bEst, aStrat, aEst),
              cvEst_c = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_c, cStrat, cEst, aStrat, aEst),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}







## DEPRECATED, will be cut in future release
bioHelper1_old <- function(x, plts, db, grpBy, aGrpBy, byPlot){

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

bioHelper2_old <- function(x, popState, a, t, grpBy, aGrpBy, method){

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





