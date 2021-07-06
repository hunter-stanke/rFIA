vegStructHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){


  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$SUBP_COND, by = c('PLT_CN', 'CONDID')) %>%
    left_join(db$P2VEG_SUBP_STRUCTURE, by = c("PLT_CN", "CONDID", 'SUBP'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp



  if (byPlot){
    ## Return proportion of plot area covered by each growth/habit layer
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, CONDID, SUBP, LAYER, GROWTH_HABIT, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN, SUBP) %>%
      summarize(cover = sum(COVER_PCT/100 * SUBPCOND_PROP * aDI, na.rm = TRUE)) %>%
      group_by(PLT_CN, !!!grpSyms) %>%
      summarize(cover = mean(cover, na.rm = TRUE)) %>%
      as.data.frame()
    a = NULL

  } else {
    grpSyms <- syms(grpBy)


    # Compute estimates
    t <- data %>%
      distinct(PLT_CN, CONDID, SUBP, LAYER, GROWTH_HABIT, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN, SUBP) %>%
      summarize(cPlot = sum(COVER_PCT/100 * SUBPCOND_PROP * aDI, na.rm = TRUE),
                plotIn_VEG = if_else(sum(aDI, na.rm = TRUE) > 0, 1, 0)) %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(cPlot = mean(cPlot, na.rm = TRUE),
                plotIn_VEG = ifelse(any(plotIn_VEG > 0), 1, 0)) %>%
      as.data.frame()

    a <- data %>%
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, PLT_CN, PROP_BASIS) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn_AREA = if_else(sum(aDI, na.rm = TRUE)  >  0, 1,0))
  }

  pltOut <- list(t = t, a = a)
  return(pltOut)

}



vegStructHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

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
      fa = fa * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!aGrpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              p2eu = dplyr::first(p2eu),
              aStrat = sum(fa, na.rm = TRUE),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE),
              ## Strata level variances
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
      cPlot = cPlot * ADJ_FACTOR_SUBP) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              p2eu = dplyr::first(p2eu),
              cStrat = sum(cPlot, na.rm = TRUE),
              aStrat = dplyr::first(aStrat),
              plotIn_VEG = sum(plotIn_VEG, na.rm = TRUE),
              ## Strata level variances
              cv = sum(cPlot^2, na.rm = TRUE),
              cvStrat_c = sum(cPlot * fa, na.rm = TRUE)) %>%
    mutate(cStrat = cStrat / nh,
           adj = nh * (nh-1),
           cv = (cv - (nh * cStrat^2)) / (nh * (nh-1)),
           cvStrat_c = (cvStrat_c - (nh * cStrat * aStrat)) / adj) %>%
    as.data.frame() %>%
    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, !!!grpSyms) %>%
    summarize(cEst = unitMean(ESTN_METHOD, a, nh,  w, cStrat),
              plotIn_VEG = sum(plotIn_VEG, na.rm = TRUE),
              N = dplyr::first(p2eu),
              A = dplyr::first(a),
              cVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cv, aStrat, aEst),
              # Unit Covariance
              cvEst_c = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_c, cStrat, cEst, cTStrat, cTEst))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}

