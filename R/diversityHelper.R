divHelper1 <- function(x, plts, db, grpBy, byPlot){

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
  data$aDI <- data$landD * data$aD * data$sp
  data$tDI <- data$landD * data$aD * data$tD * data$typeD * data$sp


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, CONDID, TREE, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(H = divIndex(grp, state  * tDI, index = 'H'),
                S = divIndex(grp, state * tDI, index = 'S'),
                Eh = divIndex(grp, state * tDI, index = 'Eh'),
                nStems = length(which(tDI == 1))) %>%
      as.data.frame()
    a = NULL
    full = NULL

  } else {

    grpSyms <- syms(grpBy)

    ## Using this to return a tree list for gamma and beta
    full <- data %>%
      filter(PLOT_STATUS_CD == 1) %>%
      mutate(state = state * tDI) %>%
      distinct(PLT_CN, !!!grpSyms, TRE_CN, grp, state)

    # Diversity is computed at the stand (condition level), and we continue to use the ratio of means estimator to get at
    #  average of the attribute of interest weighted by the area in which it occurs.
    t <- data %>%
      distinct(PLT_CN, CONDID, TREE, COND_STATUS_CD, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN, PROP_BASIS, CONDID) %>%
      # filter(tDI > 0) %>%
      summarize(condArea = dplyr::first(CONDPROP_UNADJ),
                hCond = divIndex(grp, state * tDI, index = 'H') * dplyr::first(CONDPROP_UNADJ),
                sCond = divIndex(grp, state * tDI, index = 'S') * dplyr::first(CONDPROP_UNADJ),
                EhCond = divIndex(grp, state * tDI, index = 'Eh') * dplyr::first(CONDPROP_UNADJ),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0))  %>%
      group_by(!!!grpSyms, PROP_BASIS, PLT_CN) %>%
      summarize(hPlot = sum(hCond, na.rm = TRUE),
                EhPlot = sum(EhCond, na.rm = TRUE),
                sPlot = sum(sCond * plotIn, na.rm = TRUE),
                fa = sum(condArea * aDI, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE)) %>%
      as.data.frame()
  }

  pltOut <- list(t = t, full = full)
  return(pltOut)
}



divHelper2 <- function(x, popState, t, full, grpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$P1POINTCNT <- popState[[x]]$P1POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  grpSyms <- syms(grpBy)

  ## Strata level estimates
  tEst <- t %>%
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
      fa = fa * aAdj,
      hPlot = hPlot * aAdj,
      EhPlot = EhPlot * aAdj,
      sPlot = sPlot * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              p2eu = dplyr::first(p2eu),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),

              ## dtplyr is fast, but requires a few extra steps, so we'll finish
              ## means and variances in subsequent mutate step
              aStrat = sum(fa, na.rm = TRUE),
              hStrat = sum(hPlot, na.rm = TRUE),
              ehStrat = sum(EhPlot, na.rm = TRUE),
              sStrat = sum(sPlot, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),

              ## Strata level variances
              av = sum(fa^2, na.rm = TRUE),
              hv = sum(hPlot^2, na.rm = TRUE),
              ehv = sum(EhPlot^2, na.rm = TRUE),
              sv = sum(sPlot^2, na.rm = TRUE),
              # Strata level covariances
              cvStrat_h = sum(hPlot * fa, na.rm = TRUE),
              cvStrat_eh = sum(EhPlot * fa, na.rm = TRUE),
              cvStrat_s = sum(sPlot * fa, na.rm = TRUE)) %>%
    mutate(aStrat = aStrat / nh,
           hStrat = hStrat / nh,
           ehStrat = ehStrat / nh,
           sStrat = sStrat / nh,
           adj = nh * (nh-1),
           av = (av - (nh*aStrat^2)) / adj,
           hv = (hv - (nh*hStrat^2)) / adj,
           sv = (sv - (nh*sStrat^2)) / adj,
           ehv = (ehv - (nh*ehStrat^2)) / adj,
           cvStrat_h = (cvStrat_h - (nh * hStrat * aStrat)) / adj,
           cvStrat_eh = (cvStrat_eh - (nh * ehStrat * aStrat)) / adj,
           cvStrat_s = (cvStrat_s - (nh * sStrat * aStrat)) / adj) %>%
    as.data.frame() %>%
  ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              hEst = unitMean(ESTN_METHOD, a, nh,  w, hStrat),
              ehEst = unitMean(ESTN_METHOD, a, nh,  w, ehStrat),
              sEst = unitMean(ESTN_METHOD, a, nh,  w, sStrat),
              N = dplyr::first(p2eu),
              A = dplyr::first(a),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, av, aStrat, aEst),
              hVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, hv, hStrat, hEst),
              ehVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, ehv, ehStrat, ehEst),
              sVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, sv, sStrat, sEst),
              cvEst_h = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_h, hStrat, hEst, aStrat, aEst),
              cvEst_eh = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_eh, ehStrat, ehEst, aStrat, aEst),
              cvEst_s = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_s, sStrat, sEst, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

  ## For gamma and beta, annoying but can't think of another way
  full <- full %>%
    inner_join(select(popState[[x]], c(YEAR, PLT_CN)), by = 'PLT_CN') %>%
    filter(!is.na(YEAR) & !is.na(state) & !is.na(grp))


  out <- list(tEst = tEst, full = full)

  return(out)
}



