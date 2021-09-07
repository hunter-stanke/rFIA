ssHelper1 <- function(x, plts, db, grpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$TREE, by = c('PLT_CN', 'CONDID'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp



  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = INVYR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(stage = structHelper(DIA, CCLCD),
                nStems = length(which(aDI == 1))) %>%
      as.data.frame() %>%
      ungroup() %>%
      mutate_if(is.factor, as.character)

  } else {
    grpSyms <- syms(grpBy)
    # Compute estimates
    t <- data %>%
      distinct(PLT_CN, SUBP, CONDID, TREE, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN, PROP_BASIS, CONDID) %>%
      summarize(CONDPROP_UNADJ = dplyr::first(CONDPROP_UNADJ), # Area
                stage = structHelper(DIA, CCLCD),
                aDI = dplyr::first(aDI)) %>%
      group_by(!!!grpSyms, PLT_CN, PROP_BASIS) %>%
      summarize(p = sum(CONDPROP_UNADJ[stage == 'pole'] * aDI[stage == 'pole'], na.rm = TRUE),
                ma = sum(CONDPROP_UNADJ[stage == 'mature'] * aDI[stage == 'mature'], na.rm = TRUE),
                l = sum(CONDPROP_UNADJ[stage == 'late'] * aDI[stage == 'late'], na.rm = TRUE),
                mo = sum(CONDPROP_UNADJ[stage == 'mosaic'] * aDI[stage == 'mosaic'], na.rm = TRUE),
                fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      as.data.frame()
  }

  pltOut <- list(t = t)
  return(pltOut)

}



ssHelper2 <- function(x, popState, t, grpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    #aGrpBy <- c(aGrpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$P1POINTCNT <- popState[[x]]$P1POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  grpSyms <- syms(grpBy)

  ## Strata level estimates
  tEst <- t %>%
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
      p = p * aAdj,
      ma = ma * aAdj,
      l = l * aAdj,
      mo = mo * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              p2eu = dplyr::first(p2eu),
              aStrat = sum(fa, na.rm = TRUE),
              pStrat = sum(p, na.rm = TRUE),
              maStrat = sum(ma, na.rm = TRUE),
              lStrat = sum(l, na.rm = TRUE),
              moStrat = sum(mo, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),

              # ## Strata level variances
              av = sum(fa^2, na.rm = TRUE),
              pv = sum(p^2, na.rm = TRUE),
              mav = sum(ma^2, na.rm = TRUE),
              lv = sum(l^2, na.rm = TRUE),
              mov = sum(mo^2, na.rm = TRUE),
              # Strata level covariances
              cvStrat_p = sum(p*fa, na.rm = TRUE),
              cvStrat_ma = sum(ma*fa, na.rm = TRUE),
              cvStrat_l = sum(l*fa, na.rm = TRUE),
              cvStrat_mo = sum(mo*fa, na.rm = TRUE)) %>%
    mutate(aStrat = aStrat / nh,
           pStrat = pStrat / nh,
           maStrat = maStrat / nh,
           lStrat = lStrat / nh,
           moStrat = moStrat / nh,
           adj = nh * (nh-1),
           av = (av - (nh*aStrat^2)) / adj,
           pv = (pv - (nh*pStrat^2)) / adj,
           mav = (mav - (nh*maStrat^2)) / adj,
           lv = (lv - (nh*lStrat^2)) / adj,
           mov = (mov - (nh*moStrat^2)) / adj,
           cvStrat_p = (cvStrat_p - (nh * pStrat * aStrat)) / adj,
           cvStrat_ma = (cvStrat_ma - (nh * maStrat * aStrat)) / adj,
           cvStrat_l = (cvStrat_l - (nh * lStrat * aStrat)) / adj,
           cvStrat_mo = (cvStrat_mo - (nh * moStrat * aStrat)) / adj) %>%
    as.data.frame() %>%
    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              pEst = unitMean(ESTN_METHOD, a, nh,  w, pStrat),
              maEst = unitMean(ESTN_METHOD, a, nh,  w, maStrat),
              lEst = unitMean(ESTN_METHOD, a, nh,  w, lStrat),
              moEst = unitMean(ESTN_METHOD, a, nh,  w, moStrat),
              N = dplyr::first(p2eu),
              A = dplyr::first(a),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, av, aStrat, aEst),
              pVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, pv, pStrat, pEst),
              maVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, mav, maStrat, maEst),
              lVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, lv, lStrat, lEst),
              moVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, mov, moStrat, moEst),
              cvEst_p = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_p, pStrat, pEst, aStrat, aEst),
              cvEst_ma = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_ma, maStrat, maEst, aStrat, aEst),
              cvEst_l = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_l, lStrat, lEst, aStrat, aEst),
              cvEst_mo = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_mo, moStrat, moEst, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

  out <- list(tEst = tEst)

  return(out)
}





