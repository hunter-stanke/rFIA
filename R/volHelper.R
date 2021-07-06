volHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

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
      summarize(BOLE_CF_ACRE = sum(bcf * TPA_UNADJ * tDI, na.rm = TRUE),
                SAW_CF_ACRE = sum(scf * TPA_UNADJ * tDI, na.rm = TRUE),
                SAW_MBF_ACRE = sum(sbf * TPA_UNADJ * tDI, na.rm = TRUE) / 1000,
                nStems = length(which(tDI == 1))) %>%
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
      summarize(bcfPlot = sum(bcf * TPA_UNADJ * tDI, na.rm = TRUE),
                scfPlot = sum(scf * TPA_UNADJ * tDI, na.rm = TRUE),
                sbfPlot = sum(sbf * TPA_UNADJ * tDI, na.rm = TRUE) / 1000,
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0)) %>%
      as.data.frame()
  }




  pltOut <- list(a = a, t = t)
  return(pltOut)

}



volHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

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
      bcfPlot = bcfPlot * tAdj,
      scfPlot = scfPlot * tAdj,
      sbfPlot = sbfPlot* tAdj) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, !!!grpSyms) %>%
    summarize(bcfPlot = sum(bcfPlot, na.rm = TRUE),
              scfPlot = sum(scfPlot, na.rm = TRUE),
              sbfPlot = sum(sbfPlot, na.rm = TRUE),
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
              bcfStrat = sum(bcfPlot, na.rm = TRUE),
              scfStrat = sum(scfPlot, na.rm = TRUE),
              sbfStrat = sum(sbfPlot, na.rm = TRUE),
              aStrat = dplyr::first(aStrat),
              plotIn_TREE = sum(plotIn, na.rm = TRUE),

              ## Strata level variances
              bcfv = sum(bcfPlot^2, na.rm = TRUE),
              scfv = sum(scfPlot^2, na.rm = TRUE),
              sbfv = sum(sbfPlot^2, na.rm = TRUE),

              ## Strata level covariances
              cvStrat_bcf =sum(fa* bcfPlot, na.rm = TRUE),
              cvStrat_scf =sum(fa* scfPlot, na.rm = TRUE),
              cvStrat_sbf =sum(fa* sbfPlot, na.rm = TRUE)) %>%
    mutate(bcfStrat = bcfStrat / nh,
           scfStrat = scfStrat / nh,
           sbfStrat = sbfStrat / nh,
           ## Variances
           adj = nh * (nh-1),
           bcfv = (bcfv - (nh*bcfStrat^2)) / adj,
           scfv = (scfv - (nh*scfStrat^2)) / adj,
           sbfv = (sbfv - (nh*sbfStrat^2)) / adj,
           ## Covariances
           cvStrat_bcf = (cvStrat_bcf - (nh * bcfStrat * aStrat)) / adj,
           cvStrat_scf = (cvStrat_scf - (nh * scfStrat * aStrat)) / adj,
           cvStrat_sbf =  (cvStrat_sbf - (nh * sbfStrat * aStrat)) / adj) %>%
    as.data.frame() %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(bcfEst = unitMean(ESTN_METHOD, a, nh, w, bcfStrat),
              scfEst = unitMean(ESTN_METHOD, a, nh, w, scfStrat),
              sbfEst = unitMean(ESTN_METHOD, a, nh, w, sbfStrat),
              N = dplyr::first(p2eu),
              A = dplyr::first(a),

              # Estimation of unit variance
              bcfVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bcfv, bcfStrat, bcfEst),
              scfVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, scfv, scfStrat, scfEst),
              sbfVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, sbfv, sbfStrat, sbfEst),
              cvEst_bcf = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_bcf, bcfStrat, bcfEst, aStrat, aEst),
              cvEst_scf = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_scf, scfStrat, scfEst, aStrat, aEst),
              cvEst_sbf = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_sbf, sbfStrat, sbfEst, aStrat, aEst),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}



