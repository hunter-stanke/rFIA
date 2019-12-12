ssHelper1 <- function(x, plts, db, grpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'TPA_UNADJ', 'CCLCD', 'SUBP', 'TREE')), by = c('PLT_CN', 'CONDID')) %>%
    filter(DIA >= 5)

   ## Comprehensive indicator function
   data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
   data$tDI <- data$landD * data$aD_p * data$aD_c * data$sp



  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    t <- data %>%
      mutate(YEAR = INVYR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(stage = structHelper(DIA, CCLCD),
                nStems = length(which(tDI == 1)))

  } else {
    # Compute estimates
    t <- data %>%
      distinct(PLT_CN, SUBP, CONDID, TREE, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, PROP_BASIS, CONDID) %>%
      summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ), # Area
                stage = structHelper(DIA, CCLCD),
                aDI = first(aDI)) %>%
      group_by(.dots = grpBy, PLT_CN, PROP_BASIS) %>%
      summarize(p = sum(CONDPROP_UNADJ[stage == 'pole'] * aDI[stage == 'pole'], na.rm = TRUE),
                ma = sum(CONDPROP_UNADJ[stage == 'mature'] * aDI[stage == 'mature'], na.rm = TRUE),
                l = sum(CONDPROP_UNADJ[stage == 'late'] * aDI[stage == 'late'], na.rm = TRUE),
                mo = sum(CONDPROP_UNADJ[stage == 'mosaic'] * aDI[stage == 'mosaic'], na.rm = TRUE),
                fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))
  }

  pltOut <- list(t = t)
  return(pltOut)

}



ssHelper2 <- function(x, popState, t, grpBy){

  ## Strata level estimates
  tEst <- t %>%
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
      fa = fa * aAdj,
      p = p * aAdj,
      ma = ma * aAdj,
      l = l * aAdj,
      mo = mo * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(a_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              aStrat = mean(fa * a_t, na.rm = TRUE),
              pStrat = mean(p * a_t, na.rm = TRUE),
              maStrat = mean(ma * a_t, na.rm = TRUE),
              lStrat = mean(l * a_t, na.rm = TRUE),
              moStrat = mean(mo * a_t, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              av = stratVar(ESTN_METHOD, fa, aStrat, ndif, a, nh),
              pv = stratVar(ESTN_METHOD, p, pStrat, ndif, a, nh),
              mav = stratVar(ESTN_METHOD, ma, maStrat, ndif, a, nh),
              lv = stratVar(ESTN_METHOD, l, lStrat, ndif, a, nh),
              mov = stratVar(ESTN_METHOD, mo, moStrat, ndif, a, nh),
              # Strata level covariances
              cvStrat_p = stratVar(ESTN_METHOD, p, pStrat, ndif, a, nh, fa, aStrat),
              cvStrat_ma = stratVar(ESTN_METHOD, ma, maStrat, ndif, a, nh, fa, aStrat),
              cvStrat_l = stratVar(ESTN_METHOD, l, lStrat, ndif, a, nh, fa, aStrat),
              cvStrat_mo = stratVar(ESTN_METHOD, mo, moStrat, ndif, a, nh, fa, aStrat)) %>%
    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              pEst = unitMean(ESTN_METHOD, a, nh,  w, pStrat),
              maEst = unitMean(ESTN_METHOD, a, nh,  w, maStrat),
              lEst = unitMean(ESTN_METHOD, a, nh,  w, lStrat),
              moEst = unitMean(ESTN_METHOD, a, nh,  w, moStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, av, aStrat, aEst),
              pVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pv, pStrat, pEst),
              maVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mav, maStrat, maEst),
              lVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, lv, lStrat, lEst),
              moVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mov, moStrat, moEst),
              cvEst_p = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pv, pStrat, pEst, aStrat, aEst),
              cvEst_ma = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mav, maStrat, maEst, aStrat, aEst),
              cvEst_l = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, lv, lStrat, lEst, aStrat, aEst),
              cvEst_mo = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mov, moStrat, moEst, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

  out <- list(tEst = tEst)

  return(out)
}



























standStructHelper <- function(x, combos, data, grpBy, totals, tidy, SE){
  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Area domain indicator for each column in
    aObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(aObs))) == length(aObs)) aObs <- 1
    ad <- data$aDI * aObs * ad
  }

  if(SE){
    data$aDI <- ad
    data$aDI[is.na(data$aDI)] <- 0
    # Compute estimates
    s <- data %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      # filter(tDI > 0) %>%
      summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
                stage = structHelper(DIA, CCLCD),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT),
                aAdj = first(aAdj),
                aDI = first(aDI)) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(p = sum(CONDPROP_UNADJ[stage == 'pole'] * aDI[stage == 'pole'] * aAdj[stage == 'pole'], na.rm = TRUE),
                ma = sum(CONDPROP_UNADJ[stage == 'mature'] * aDI[stage == 'mature'] * aAdj[stage == 'mature'], na.rm = TRUE),
                l = sum(CONDPROP_UNADJ[stage == 'late'] * aDI[stage == 'late'] * aAdj[stage == 'late'], na.rm = TRUE),
                mo = sum(CONDPROP_UNADJ[stage == 'mosaic'] * aDI[stage == 'mosaic'] * aAdj[stage == 'mosaic'], na.rm = TRUE),
                faFull = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
                p1EU = first(p1EU),
                a = first(a),
                p1 = first(p1),
                p2 = first(p2)) %>%
      # Continue through totals
      #d <- dInt %>%
      #filter(plotIn > 0) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(pStrat = mean(p, na.rm = TRUE),
                maStrat = mean(ma, na.rm = TRUE),
                lStrat = mean(l, na.rm = TRUE),
                moStrat = mean(mo, na.rm = TRUE),
                fullStrat = mean(faFull, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                #vPlots = ifelse(plotIn > 1, plotIn, NA),
                a = first(a),
                w = first(p1) / first(p1EU), # Stratum weight
                nh = first(p2), # Number plots in stratum
                # Strata level variances
                pv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(p * first(a) / nh),
                            (sum(p^2, na.rm = TRUE) - sum(nh * pStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
                mav = ifelse(first(ESTN_METHOD == 'simple'),
                             var(ma * first(a) / nh),
                             (sum(ma^2, na.rm = TRUE) - sum(nh * maStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
                lv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(l * first(a) / nh),
                            (sum(l^2, na.rm = TRUE) - sum(nh * lStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
                mov = ifelse(first(ESTN_METHOD == 'simple'),
                             var(mo * first(a) / nh),
                             (sum(mo^2, na.rm = TRUE) - sum(nh * moStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
                fullv = ifelse(first(ESTN_METHOD == 'simple'),
                               var(faFull * first(a) / nh),
                               (sum(faFull^2, na.rm = TRUE) - sum(nh * fullStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
                # cvStrat_t = ifelse(first(ESTN_METHOD == 'simple'),
                #                    cov(fa,tPlot),
                #                    (sum(fa*tPlot) - sum(nh * aStrat *tStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_p = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(faFull,p),
                                   (sum(faFull*p, na.rm = TRUE) - sum(nh * pStrat *fullStrat, na.rm = TRUE)) / (nh * (nh-1))),
                cvStrat_ma = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(faFull,ma),
                                    (sum(faFull*ma, na.rm = TRUE) - sum(nh * maStrat *fullStrat, na.rm = TRUE)) / (nh * (nh-1))),
                cvStrat_l = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(faFull,l),
                                   (sum(faFull*l, na.rm = TRUE) - sum(nh * lStrat *fullStrat, na.rm = TRUE)) / (nh * (nh-1))),
                cvStrat_mo = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(faFullmo),
                                    (sum(faFull*mo, na.rm = TRUE) - sum(nh * moStrat *fullStrat, na.rm = TRUE)) / (nh * (nh-1)))) %>% # Stratified and double cases) %>% # Stratified and double cases
      group_by(ESTN_UNIT_CN) %>%
      summarize(plotIn = sum(plotIn, na.rm = TRUE),
                #vPlots = ifelse(plotIn > 1, plotIn, NA),
                pEst = unitMean(ESTN_METHOD, a, plotIn, w, pStrat),
                maEst = unitMean(ESTN_METHOD, a, plotIn, w, maStrat),
                lEst = unitMean(ESTN_METHOD, a, plotIn, w, lStrat),
                moEst = unitMean(ESTN_METHOD, a, plotIn, w, moStrat),
                fullEst = unitMean(ESTN_METHOD, a, plotIn, w, fullStrat),
                # Estimation of unit variance,
                pVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, pv, pStrat, pEst),
                maVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, mav, maStrat, maEst),
                lVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, lv, lStrat, lEst),
                moVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, mov, moStrat, moEst),
                fullVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, fullv, fullStrat, fullEst),
                cvEst_p = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_p, pStrat, pEst, fullStrat, fullEst),
                cvEst_ma = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_ma, maStrat, maEst, fullStrat, fullEst),
                cvEst_l = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_l, lStrat, lEst, fullStrat, fullEst),
                cvEst_mo = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_mo, moStrat, moEst, fullStrat, fullEst)) %>%
      # Compute totals
      summarize(POLE_AREA = sum(pEst, na.rm = TRUE),
                MATURE_AREA = sum(maEst, na.rm = TRUE),
                LATE_AREA = sum(lEst, na.rm = TRUE),
                MOSAIC_AREA = sum(moEst, na.rm = TRUE),
                TOTAL_AREA = sum(fullEst, na.rm = TRUE),
                POLE_PERC = POLE_AREA / TOTAL_AREA * 100,
                MATURE_PERC = MATURE_AREA / TOTAL_AREA * 100,
                LATE_PERC = LATE_AREA / TOTAL_AREA * 100,
                MOSAIC_PERC = MOSAIC_AREA / TOTAL_AREA * 100,
                nPlots = sum(plotIn, na.rm = TRUE),
                pVar = sum(pVar, na.rm = TRUE),
                maVar = sum(maVar, na.rm = TRUE),
                lVar = sum(lVar, na.rm = TRUE),
                moVar = sum(moVar, na.rm = TRUE),
                fVar = sum(fullVar, na.rm = TRUE),
                cvP = sum(cvEst_p, na.rm = TRUE),
                cvMa = sum(cvEst_ma, na.rm = TRUE),
                cvL = sum(cvEst_l, na.rm = TRUE),
                cvMo = sum(cvEst_mo, na.rm = TRUE),
                ppVar = (1/TOTAL_AREA^2) * (pVar + (POLE_PERC^2 * fVar) - 2 * POLE_PERC * cvP),
                mapVar = (1/TOTAL_AREA^2) * (maVar + (MATURE_PERC^2 * fVar) - 2 * MATURE_PERC * cvMa),
                lpVar = (1/TOTAL_AREA^2) * (lVar + (LATE_PERC^2 * fVar) - 2 * LATE_PERC * cvL),
                mopVar = (1/TOTAL_AREA^2) * (moVar + (MOSAIC_PERC^2 * fVar) - 2 * MOSAIC_PERC * cvMo),
                POLE_AREA_SE = sqrt(pVar) / POLE_AREA * 100,
                MATURE_AREA_SE = sqrt(maVar) / MATURE_AREA * 100,
                LATE_AREA_SE = sqrt(lVar) / LATE_AREA * 100,
                MOSAIC_AREA_SE = sqrt(moVar) / MOSAIC_AREA * 100,
                POLE_PERC_SE = sqrt(ppVar) / POLE_PERC * 100,
                MATURE_PERC_SE = sqrt(mapVar) / MATURE_PERC * 100,
                LATE_PERC_SE = sqrt(lpVar) / LATE_PERC * 100,
                MOSAIC_PERC_SE = sqrt(mopVar) / MOSAIC_PERC * 100,
                TOTAL_AREA_SE = sqrt(fVar) / TOTAL_AREA * 100) %>%
      select(POLE_PERC, MATURE_PERC, LATE_PERC, MOSAIC_PERC,
             POLE_AREA, MATURE_AREA, LATE_AREA, MOSAIC_AREA, TOTAL_AREA,
             POLE_PERC_SE, MATURE_PERC_SE, LATE_PERC_SE, MOSAIC_PERC_SE,
             POLE_AREA_SE, MATURE_AREA_SE, LATE_AREA_SE, MOSAIC_AREA_SE,
             TOTAL_AREA_SE, nPlots)
    if(totals == FALSE){
      s <- s %>%
        select(POLE_PERC, MATURE_PERC, LATE_PERC, MOSAIC_PERC,
               POLE_PERC_SE, MATURE_PERC_SE, LATE_PERC_SE, MOSAIC_PERC_SE, nPlots)
    }

    # Rejoin w/ groupby names
    s <- data.frame(combos[[x]], s)

  } else {
    # Compute estimates
    s <- data %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      # filter(tDI > 0) %>%
      summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
                stage = structHelper(DIA, CCLCD),
                aAdj = first(aAdj),
                aDI = first(aDI),
                EXPNS = first(EXPNS)) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(p = sum(CONDPROP_UNADJ[stage == 'pole'] * aDI[stage == 'pole'] * aAdj[stage == 'pole'] * EXPNS[stage == 'pole'], na.rm = TRUE),
                ma = sum(CONDPROP_UNADJ[stage == 'mature'] * aDI[stage == 'mature'] * aAdj[stage == 'mature'] * EXPNS[stage == 'mature'], na.rm = TRUE),
                l = sum(CONDPROP_UNADJ[stage == 'late'] * aDI[stage == 'late'] * aAdj[stage == 'late'] * EXPNS[stage == 'late'], na.rm = TRUE),
                mo = sum(CONDPROP_UNADJ[stage == 'mosaic'] * aDI[stage == 'mosaic'] * aAdj[stage == 'mosaic'] * EXPNS[stage == 'mosaic'], na.rm = TRUE),
                faFull = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      group_by(.dots = grpBy) %>%
      summarize(POLE_AREA = sum(p, na.rm = TRUE),
                MATURE_AREA = sum(ma, na.rm = TRUE),
                LATE_AREA = sum(l, na.rm = TRUE),
                MOSAIC_AREA = sum(mo, na.rm = TRUE),
                TOTAL_AREA = sum(faFull, na.rm = TRUE),
                POLE_PERC = POLE_AREA / TOTAL_AREA * 100,
                MATURE_PERC = MATURE_AREA / TOTAL_AREA * 100,
                LATE_PERC = LATE_AREA / TOTAL_AREA * 100,
                MOSAIC_PERC = MOSAIC_AREA / TOTAL_AREA * 100,
                nPlots = sum(plotIn, na.rm = TRUE)) %>%
      select(grpBy, POLE_PERC, MATURE_PERC, LATE_PERC, MOSAIC_PERC,
             POLE_AREA, MATURE_AREA, LATE_AREA, MOSAIC_AREA, TOTAL_AREA,
             nPlots)

    if(totals == FALSE){
      s <- s %>%
        select(grpBy, POLE_PERC, MATURE_PERC, LATE_PERC, MOSAIC_PERC, nPlots)
    }


  } # End SE Conditional

  # Do some cleanup
  #gc()

  #Return a dataframe
  s


}
