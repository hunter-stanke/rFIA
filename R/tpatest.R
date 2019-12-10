plotsum <- function(x, plts, db, grpBy, aGrpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SPCD', 'TPA_UNADJ', 'SUBP', 'TREE', grpT, 'tD', 'typeD')), by = c('PLT_CN', 'CONDID')) %>%
    ## Need a code that tells us where the tree was measured
    ## macroplot, microplot, subplot
    mutate(PLOT_BASIS = case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA) ~ NA_character_,
      ## When DIA is less than 5", use microplot value
      DIA < 5 ~ 'MICR',
      ## When DIA is greater than 5", use subplot value
      DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
      DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) #%>%
  #filter(!is.na(PLOT_BASIS))
  # rename(YEAR = INVYR) %>%
  # mutate_if(is.factor,
  #           as.character) %>%
  # filter(!is.na(YEAR))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
  data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
                BAA = sum(basalArea(DIA) * TPA_UNADJ * tDI, na.rm = TRUE),
                TPA_PERC = TPA / sum(TPA_UNADJ * pDI, na.rm = TRUE) * 100,
                BAA_PERC = BAA / sum(basalArea(DIA) * TPA_UNADJ * pDI, na.rm = TRUE) * 100,
                nStems = length(which(tDI == 1)))

    a = NULL

  } else {
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
      filter(!is.na(PLOT_BASIS)) %>%
      group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
      summarize(tPlot = sum(TPA_UNADJ * tDI, na.rm = TRUE),
                bPlot = sum(basalArea(DIA) * TPA_UNADJ * tDI, na.rm = TRUE),
                tTPlot = sum(TPA_UNADJ * pDI, na.rm = TRUE),
                bTPlot = sum(basalArea(DIA) * TPA_UNADJ * pDI, na.rm = TRUE),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0))
  }




  pltOut <- list(a = a, t = t)
  return(pltOut)

}



estsum <- function(x, popState, a, t, grpBy, aGrpBy){

  ## Strata level estimates
  aStrat <- a %>%
    ## Rejoin with population tables
    right_join(popState[[x]], by = 'PLT_CN') %>%
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
      fa = fa * aAdj) %>%
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
              av = ifelse(first(ESTN_METHOD == 'simple'),
                          var(c(fa, numeric(ndif)) * first(a) / nh),
                          (sum((c(fa, numeric(ndif))^2)) - nh * aStrat^2) / (nh * (nh-1))))
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
    right_join(popState[[x]], by = 'PLT_CN') %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% 'YEAR' == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% 'YEAR' == FALSE])) %>%
    #Add adjustment factors
    mutate(
      ## AREA
      tAdj = case_when(
        ## When NA, stay NA
        is.na(PLOT_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PLOT_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PLOT_BASIS == 'SUBP' ~ as.numeric(ADJ_FACTOR_SUBP),
        PLOT_BASIS == 'MICR' ~ as.numeric(ADJ_FACTOR_MICR)),
      ## AREA
      aAdj = case_when(
        ## When NA, stay NA
        is.na(PROP_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = fa * aAdj,
      tPlot = tPlot * tAdj,
      bPlot = bPlot * tAdj,
      tTPlot = tTPlot * tAdj,
      bTPlot = bTPlot * tAdj) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(tPlot = sum(tPlot, na.rm = TRUE),
              bPlot = sum(bPlot, na.rm = TRUE),
              tTPlot = sum(tTPlot, na.rm = TRUE),
              bTPlot = sum(bTPlot, na.rm = TRUE),
              fa = first(fa),
              plotIn = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
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
              bStrat = mean(bPlot * r_t, na.rm = TRUE),
              tTStrat = mean(tTPlot * r_t, na.rm = TRUE),
              bTStrat = mean(bTPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_TREE = sum(plotIn, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              ## Strata level variances
              #aVar = (sum(forArea^2) - sum(P2POINTCNT * aStrat^2)) / (P2POINTCNT * (P2POINTCNT-1)),
              tv = ifelse(first(ESTN_METHOD == 'simple'),
                          var(c(tPlot, numeric(ndif)) * first(a) / nh),
                          (sum(c(tPlot, numeric(ndif))^2) - sum(nh * tStrat^2)) / (nh * (nh-1))), # Stratified and double cases
              bv = ifelse(first(ESTN_METHOD == 'simple'),
                          var(c(tPlot, numeric(ndif))* first(a) / nh),
                          (sum(c(bPlot, numeric(ndif))^2) - sum(nh * bStrat^2)) / (nh * (nh-1))),
              tTv = ifelse(first(ESTN_METHOD == 'simple'),
                           var(c(tTPlot, numeric(ndif)) * first(a) / nh),
                           (sum(c(tTPlot, numeric(ndif))^2) - sum(nh * tTStrat^2)) / (nh * (nh-1))), # Stratified and double cases
              bTv = ifelse(first(ESTN_METHOD == 'simple'),
                           var(c(bTPlot, numeric(ndif)) * first(a) / nh),
                           (sum(c(bTPlot, numeric(ndif))^2) - sum(nh * bTStrat^2)) / (nh * (nh-1))),
              # Strata level covariances
              cvStrat_t = ifelse(first(ESTN_METHOD == 'simple'),
                                 cov(fa,tPlot),
                                 (sum(fa*tPlot, na.rm = TRUE) - sum(nh * aStrat *tStrat)) / (nh * (nh-1))), # Stratified and double cases
              cvStrat_b = ifelse(first(ESTN_METHOD == 'simple'),
                                 cov(fa,bPlot),
                                 (sum(fa*bPlot, na.rm = TRUE) - sum(nh * aStrat *bStrat)) / (nh * (nh-1))),
              cvStrat_tT = ifelse(first(ESTN_METHOD == 'simple'),
                                  cov(tTPlot,tPlot),
                                  (sum(tTPlot*tPlot, na.rm = TRUE) - sum(nh * tTStrat *tStrat)) / (nh * (nh-1))), # Stratified and double cases
              cvStrat_bT = ifelse(first(ESTN_METHOD == 'simple'),
                                  cov(bTPlot,bPlot),
                                  (sum(bTPlot*bPlot, na.rm = TRUE) - sum(nh * bTStrat *bStrat)) / (nh * (nh-1)))) %>%
    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(tEst = unitMean(ESTN_METHOD, a, nh,  w, tStrat),
              bEst = unitMean(ESTN_METHOD, a, nh,  w, bStrat),
              tTEst = unitMean(ESTN_METHOD, a, nh,  w, tTStrat),
              bTEst = unitMean(ESTN_METHOD, a, nh,  w, bTStrat),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE),
              tVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, tv, tStrat, tEst),
              bVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bv, bStrat, bEst),
              tTVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, tTv, tTStrat, tTEst),
              bTVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bTv, bTStrat, bTEst),
              # Unit Covariance
              cvEst_t = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_t, tStrat, tEst, aStrat, aEst),
              cvEst_b = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_b, bStrat, bEst, aStrat, aEst),
              cvEst_tT = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_t, tStrat, tEst, tTStrat, tTEst),
              cvEst_bT = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_b, bStrat, bEst, bTStrat, bTEst))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}

