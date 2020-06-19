vegStructHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)

  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$SUBP_COND, by = c('PLT_CN', 'CONDID')) %>%
    left_join(db$P2VEG_SUBP_STRUCTURE, by = c("PLT_CN", "CONDID", 'SUBP'))
  #filter(DIA >= 5)

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp



  if (byPlot){
    ## Return proportion of plot area covered by each growth/habit layer
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, CONDID, SUBP, LAYER, GROWTH_HABIT, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, SUBP) %>%
      summarize(cover = sum(COVER_PCT/100 * SUBPCOND_PROP * aDI, na.rm = TRUE)) %>%
      group_by(PLT_CN, .dots = grpBy) %>%
      summarize(cover = mean(cover, na.rm = TRUE))
    a = NULL

  } else {
    # Compute estimates
    t <- data %>%
      distinct(PLT_CN, CONDID, SUBP, LAYER, GROWTH_HABIT, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, SUBP) %>%
      summarize(cPlot = sum(COVER_PCT/100 * SUBPCOND_PROP * aDI, na.rm = TRUE),
                plotIn_VEG = if_else(sum(aDI, na.rm = TRUE) > 0, 1, 0)) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(cPlot = mean(cPlot, na.rm = TRUE),
                plotIn_VEG = ifelse(any(plotIn_VEG > 0), 1, 0))

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
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR

  }

  ## Strata level estimates
  aStrat <- a %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
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
      fa = fa * aAdj) #%>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = aGrpBy) %>%
    summarize(a_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              aStrat = mean(fa * a_t, na.rm = TRUE),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              p2eu = first(p2eu),
              ndif = nh - n,
              ## Strata level variances
              av = stratVar(ESTN_METHOD, fa, aStrat, ndif, a, nh))
  ## Estimation unit
  aEst <- aStrat %>%
    group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, av, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))


  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
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
      fa = fa * aAdj,
      cPlot = cPlot * ADJ_FACTOR_SUBP) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              cStrat = mean(cPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_VEG = sum(plotIn_VEG, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              p2eu = first(p2eu),
              ndif = nh - n,
              ## Strata level variances
              ## Strata level variances
              cv = stratVar(ESTN_METHOD, cPlot, cStrat, ndif, a, nh),
              cvStrat_c = stratVar(ESTN_METHOD, cPlot, cStrat, ndif, a, nh, fa, aStrat)) %>%
    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(cEst = unitMean(ESTN_METHOD, a, nh,  w, cStrat),
              plotIn_VEG = sum(plotIn_VEG, na.rm = TRUE),
              cVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cv, aStrat, aEst),
              # Unit Covariance
              cvEst_c = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_c, cStrat, cEst, cTStrat, cTEst))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}

