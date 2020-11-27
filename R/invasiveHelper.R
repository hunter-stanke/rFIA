invHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]

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
    grpSyms <- syms(grpBy)
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      filter(!is.na(SYMBOL)) %>%
      distinct(PLT_CN, CONDID, SUBP, VEG_SPCD, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN, SUBP) %>%
      summarize(cover = sum(COVER_PCT/100 * SUBPCOND_PROP * aDI, na.rm = TRUE)) %>%
      group_by(PLT_CN, !!!grpSyms) %>%
      summarize(cover = mean(cover, na.rm = TRUE)) %>%
      as.data.frame()
    a = NULL

  } else {
    grpSyms <- syms(grpBy)
    aGrpSyms <- syms(aGrpBy)
    # Compute estimates
    t <- data %>%
      distinct(PLT_CN, CONDID, SUBP, VEG_SPCD, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN, SUBP) %>%
      summarize(iPlot = sum(COVER_PCT/100 * SUBPCOND_PROP * aDI, na.rm = TRUE),
                plotIn_INV = ifelse(sum(aDI, na.rm = TRUE) >  0 &dplyr::first(INVASIVE_SAMPLING_STATUS_CD) == 1, 1,0)) %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(iPlot = mean(iPlot, na.rm = TRUE),
                plotIn_INV = ifelse(any(plotIn_INV > 0), 1, 0)) %>%
      as.data.frame()

    a <- data %>%
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!aGrpSyms, PLT_CN, PROP_BASIS) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn_AREA = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      as.data.frame()
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
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE),
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

  grpSyms <- syms(grpBy)

  ## Strata level estimates
  tEst <- t %>%
    lazy_dt() %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
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
      fa = fa * aAdj,
      iPlot = iPlot * ADJ_FACTOR_SUBP) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh =dplyr::first(P2POINTCNT),
              a =dplyr::first(AREA_USED),
              w =dplyr::first(P1POINTCNT) /dplyr::first(P1PNTCNT_EU),
              p2eu =dplyr::first(p2eu),
              iStrat = sum(iPlot, na.rm = TRUE),
              aStrat =dplyr::first(aStrat),
              plotIn_INV = sum(plotIn_INV, na.rm = TRUE),
              ## Strata level variances
              iv = sum(iPlot^2, na.rm = TRUE),
              cvStrat_i = sum(iPlot * fa, na.rm = TRUE)) %>%
    mutate(iStrat = iStrat / nh,
           adj = nh * (nh-1),
           iv = (iv - (nh*iStrat^2)) / adj,
           cvStrat_i = (cvStrat_i - (nh * iStrat * aStrat)) / adj) %>%
    as.data.frame() %>%
    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(iEst = unitMean(ESTN_METHOD, a, nh,  w, iStrat),
              N =dplyr::first(p2eu),
              plotIn_INV = sum(plotIn_INV, na.rm = TRUE),
              iVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh,dplyr::first(p2eu), w, iv, iStrat, iEst),
              # Unit Covariance
              cvEst_i = unitVarNew(method = 'cov', ESTN_METHOD, a, nh,dplyr::first(p2eu), w, cvStrat_i, iStrat, iEst, aStrat, aEst))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}






