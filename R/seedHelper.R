seedHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$SEEDLING, by = c('PLT_CN', 'CONDID'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp
  data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp


  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
                TPA_PERC = sum(TPA_UNADJ * pDI, na.rm = TRUE) * 100,
                nStems = sum(TREECOUNT_CALC,na.rm = TRUE)) %>%
      mutate(TPA_PERC = TPA / TPA_PERC) %>%
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
      group_by(PLT_CN, !!!grpSyms) %>%
      summarize(tPlot = sum(TPA_UNADJ * tDI, na.rm = TRUE),
                tTPlot = sum(TPA_UNADJ * pDI, na.rm = TRUE),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0)) %>%
      as.data.frame()
  }




  pltOut <- list(a = a, t = t)
  return(pltOut)

}



seedHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    aGrpBy <- c(aGrpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$P1POINTCNT <- popState[[x]]$P1POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  aGrpSyms <- syms(aGrpBy)
  grpSyms <- syms(grpBy)

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

  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
    mutate(
      ## AREA
      tAdj = as.numeric(ADJ_FACTOR_MICR),
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
      tPlot = tPlot * tAdj,
      tTPlot = tTPlot * tAdj) %>%
    # ## Extra step for variance issues
    # group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    # summarize(tPlot = sum(tPlot, na.rm = TRUE),
    #           tTPlot = sum(tTPlot, na.rm = TRUE),
    #           fa = dplyr::first(fa),
    #           plotIn = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
    #           nh = dplyr::first(P2POINTCNT),
    #           p2eu = dplyr::first(p2eu),
    #           a = dplyr::first(AREA_USED),
    #           w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU)) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              p2eu = dplyr::first(p2eu),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              tStrat = sum(tPlot, na.rm = TRUE),
              tTStrat = sum(tTPlot, na.rm = TRUE),
              aStrat = dplyr::first(aStrat),
              plotIn_TREE = sum(plotIn, na.rm = TRUE),

              ## Strata level variances
              tv = sum(tPlot^2, na.rm = TRUE),
              tTv = sum(tTPlot^2, na.rm = TRUE),

              ## Strata level covariances
              cvStrat_t = sum(fa*tPlot, na.rm = TRUE),
              cvStrat_tT = sum(tPlot*tTPlot,na.rm = TRUE)) %>%
    mutate(tStrat = tStrat / nh,
           tTStrat = tTStrat / nh,
           adj = nh * (nh-1),
           tv = (tv - (nh*tStrat^2)) / adj,
           tTv = (tTv - (nh*tTStrat^2)) / adj,
           cvStrat_t = (cvStrat_t - (nh * tStrat * aStrat)) / adj,
           cvStrat_tT = (cvStrat_tT - (nh * tStrat * tTStrat)) / adj) %>%
    as.data.frame() %>%
    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(tEst = unitMean(ESTN_METHOD, a, nh,  w, tStrat),
              tTEst = unitMean(ESTN_METHOD, a, nh,  w, tTStrat),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE),
              N = dplyr::first(p2eu),
              A = dplyr::first(a),

              tVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, tv, tStrat, tEst),
              tTVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, tTv, tTStrat, tTEst),
              # Unit Covariance
              cvEst_t = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_t, tStrat, tEst, aStrat, aEst),
              cvEst_tT = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_t, tStrat, tEst, tTStrat, tTEst))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}







