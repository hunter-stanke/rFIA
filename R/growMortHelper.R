growMortHelper <- function(x, combos, data, grpBy, aGrpBy, totals, SE){
  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    td <- data$tDI * tObs * td
    # Area domain indicator for each column in
    if(grpBy[n] %in% aGrpBy){
      aObs <- as.character(combos[[x]][[aGrpBy[n]]]) == as.character(data[[aGrpBy[n]]])
      aObs[is.na(aObs)] <- 0
      ad <- data$aDI * aObs * ad
    }
  }

  if(SE){
    data$tDI <- td
    data$tDI[is.na(data$tDI)] <- 0
    data$aDI <- ad
    data$aDI[is.na(data$aDI)] <- 0
    ## We produce an intermediate object in this chain as it is needed to compute the ratio of means variance
    ## Numerator and denominator are in different domains of interest, and may be grouped by different variables
    ## see covariance estimation below

    # ### Compute total TREES in domain of interest
    # tInt <- data %>%
    #   distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, COMPONENT, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
    #   #filter(EVALID %in% tID) %>%
    #   # Compute estimates at plot level
    #   group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%

      ### Compute total TREES in domain of interest
    tInt <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      #filter(EVALID %in% tID) %>%
      # Compute estimates at plot level
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE),
                rPlot = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'], na.rm = TRUE),
                mPlot = sum(TPAMORT_UNADJ * tAdj * tDI, na.rm = TRUE),
                hPlot = sum(TPAREMV_UNADJ * tAdj * tDI, na.rm = TRUE),
                #prevPop = tPlot + mPlot * first(REMPER) + hPlot * first(REMPER) - rPlot * first(REMPER),
                #lPlot = (tPlot / ifelse(prevPop < 1, NA, prevPop)) ^ (1/first(REMPER)) - 1,
                #REMPER = first(REMPER),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT))
    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      #filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT))


    ## Compute COVARIANCE between numerator and denominator (for ratio estimates of variance)
    t <- tInt %>%
      inner_join(aInt, by = c('PLT_CN', 'ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN'), suffix = c('_t', '_a')) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(tStrat = mean(tPlot, na.rm = TRUE),
                rStrat = mean(rPlot, na.rm = TRUE),
                mStrat = mean(mPlot, na.rm = TRUE),
                hStrat = mean(hPlot, na.rm = TRUE),
                #lStrat = mean(lPlot, na.rm = TRUE),
                aStrat = mean(fa, na.rm = TRUE),
                a = first(a_t),
                w = first(p1_t) / first(p1EU_a), # Stratum weight
                nh = first(p2_t), # Number plots in stratum
                nPlots_TREE = sum(plotIn_t, na.rm = TRUE),
                nPlots_AREA = sum(plotIn_a, na.rm = TRUE),
                # Strata level variances
                tv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(tPlot * first(a) / nh),
                            (sum(tPlot^2) - sum(nh * tStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                rv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(rPlot * first(a) / nh),
                            (sum(rPlot^2) - sum(nh * rStrat^2)) / (nh * (nh-1))),
                mv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(mPlot * first(a) / nh),
                            (sum(mPlot^2) - sum(nh * mStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                hv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(hPlot * first(a) / nh),
                            (sum(hPlot^2) - sum(nh * hStrat^2)) / (nh * (nh-1))),
                #lv = ifelse(first(ESTN_METHOD == 'simple'),
                #             var(lPlot * first(a) / nh),
                #             (sum(lPlot^2) - sum(nh * lStrat^2)) / (nh * (nh-1))),
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * first(a) / nh),
                            (sum(fa^2) - sum(nh * aStrat^2)) / (nh * (nh-1))),
                # Strata level covariances
                cvStrat_r = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(fa,rPlot),
                                   (sum(fa*rPlot) - sum(nh * aStrat *rStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_m = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(fa,mPlot),
                                   (sum(fa*mPlot) - sum(nh * aStrat *mStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_h = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(fa,hPlot),
                                   (sum(fa*hPlot) - sum(nh * aStrat *hStrat)) / (nh * (nh-1))),
                #cvStrat_l = ifelse(first(ESTN_METHOD == 'simple'),
               #                     cov(fa,lPlot),
               #                     (sum(fa*lPlot) - sum(nh * aStrat *lStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_rT = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(tPlot,rPlot),
                                    (sum(tPlot*rPlot) - sum(nh * tStrat *rStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_mT = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(tPlot,mPlot),
                                    (sum(tPlot*mPlot) - sum(nh * tStrat *mStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_hT = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(tPlot,hPlot),
                                    (sum(tPlot*hPlot) - sum(nh * tStrat *hStrat)) / (nh * (nh-1)))) %>% # Stratified and double casesc) %>%
      # Estimation Unit
      group_by(ESTN_UNIT_CN) %>%
      summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
                rEst = unitMean(ESTN_METHOD, a, nh, w, rStrat),
                mEst = unitMean(ESTN_METHOD, a, nh, w, mStrat),
                hEst = unitMean(ESTN_METHOD, a, nh, w, hStrat),
                #lEst = unitMean(ESTN_METHOD, a, nh, w, lStrat),
                aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
                nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE),
                #Variance estimates
                tVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, tv, tStrat, tEst),
                rVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, rv, rStrat, rEst),
                mVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, mv, mStrat, mEst),
                hVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, hv, hStrat, hEst),
                #lVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, lv, lStrat, lEst),
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
                # Covariance estimates
                cvEst_r = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_r, rStrat, rEst, aStrat, aEst),
                cvEst_m = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_m, mStrat, mEst, aStrat, aEst),
                cvEst_h = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_h, hStrat, hEst, aStrat, aEst),
                #cvEst_l = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_l, lStrat, lEst, aStrat, aEst),
                cvEst_rT = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_rT, rStrat, rEst, tStrat, tEst),
                cvEst_mT = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_mT, mStrat, mEst, tStrat, tEst),
                cvEst_hT = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_hT, hStrat, hEst, tStrat, tEst)) %>%
      ## Full region
      summarize(TREE_TOTAL = sum(tEst, na.rm = TRUE),
                RECR_TOTAL = sum(rEst, na.rm = TRUE),
                MORT_TOTAL = sum(mEst, na.rm = TRUE),
                REMV_TOTAL = sum(hEst, na.rm = TRUE),
                AREA_TOTAL = sum(aEst, na.rm = TRUE),
                RECR_TPA = RECR_TOTAL / AREA_TOTAL,
                MORT_TPA = MORT_TOTAL / AREA_TOTAL,
                REMV_TPA = REMV_TOTAL / AREA_TOTAL,
                #LAMBDA = sum(lEst, na.rm = TRUE), # / AREA_TOTAL,
                RECR_PERC = RECR_TOTAL / TREE_TOTAL * 100,
                MORT_PERC = MORT_TOTAL / TREE_TOTAL * 100,
                REMV_PERC = REMV_TOTAL / TREE_TOTAL * 100,
                # Variance/covariance
                tVar = sum(tVar, na.rm = TRUE),
                rVar = sum(rVar, na.rm = TRUE),
                mVar = sum(mVar, na.rm = TRUE),
                hVar = sum(hVar, na.rm = TRUE),
                #lVar = sum(lVar, na.rm = TRUE),
                aVar = sum(aVar, na.rm = TRUE),
                cvR = sum(cvEst_r, na.rm = TRUE),
                cvM = sum(cvEst_m, na.rm = TRUE),
                cvH = sum(cvEst_h, na.rm = TRUE),
                #cvL = sum(cvEst_l, na.rm = TRUE),
                cvRT = sum(cvEst_rT, na.rm = TRUE),
                cvMT = sum(cvEst_mT, na.rm = TRUE),
                cvHT = sum(cvEst_hT, na.rm = TRUE),
                raVar = (1/AREA_TOTAL^2) * (rVar + (RECR_TPA^2 * aVar) - 2 * RECR_TPA * cvR),
                maVar = (1/AREA_TOTAL^2) * (mVar + (MORT_TPA^2 * aVar) - 2 * MORT_TPA * cvM),
                haVar = (1/AREA_TOTAL^2) * (hVar + (REMV_TPA^2 * aVar) - 2 * REMV_TPA * cvH),
                #laVar = (1/AREA_TOTAL^2) * (lVar + (LAMBDA^2 * aVar) - 2 * LAMBDA * cvL),
                rtVar = (1/TREE_TOTAL^2) * (rVar + (RECR_PERC^2 * tVar) - 2 * RECR_PERC * cvRT),
                mtVar = (1/TREE_TOTAL^2) * (mVar + (MORT_PERC^2 * tVar) - 2 * MORT_PERC * cvMT),
                htVar = (1/TREE_TOTAL^2) * (hVar + (REMV_PERC^2 * tVar) - 2 * REMV_PERC * cvHT),
                # Sampling Errors
                TREE_TOTAL_SE = sqrt(tVar) / TREE_TOTAL * 100,
                RECR_TOTAL_SE = sqrt(rVar) / RECR_TOTAL * 100,
                MORT_TOTAL_SE = sqrt(mVar) / MORT_TOTAL * 100,
                REMV_TOTAL_SE = sqrt(hVar) / REMV_TOTAL * 100,
                #LAMBDA_SE = sqrt(laVar) / LAMBDA * 100,
                AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
                RECR_TPA_SE = sqrt(raVar) / RECR_TPA * 100,
                MORT_TPA_SE = sqrt(maVar) / MORT_TPA * 100,
                REMV_TPA_SE = sqrt(haVar) / REMV_TPA * 100,
                RECR_PERC_SE = sqrt(rtVar) / RECR_TPA * 100,
                MORT_PERC_SE = sqrt(mtVar) / MORT_TPA * 100,
                REMV_PERC_SE = sqrt(htVar) / REMV_TPA * 100,
                # Non-zero plots
                nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
                nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE))

    # Make some columns go away
    if (totals) {
      t <- t %>%
        select(RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC, #LAMBDA,
               TREE_TOTAL, RECR_TOTAL, MORT_TOTAL, REMV_TOTAL, AREA_TOTAL,
               RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE,
               RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE, #LAMBDA_SE,
               TREE_TOTAL_SE, RECR_TOTAL_SE, MORT_TOTAL_SE, REMV_TOTAL_SE, AREA_TOTAL_SE,
               nPlots_TREE, nPlots_AREA)
    } else {
      t <- t %>%
        select(RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC, #LAMBDA,
               RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE,
               RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE, #LAMBDA_SE,
               nPlots_TREE, nPlots_AREA)
    }
    # Rejoin with some grpBy Names
    t <- data.frame(combos[[x]], t)


  } else { # No sampling errors
    ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
    ### Compute total TREES in domain of interest
    tInt <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      #filter(EVALID %in% tID) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE),
                rPlot = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'], na.rm = TRUE),
                mPlot = sum(TPAMORT_UNADJ * tAdj * tDI, na.rm = TRUE),
                hPlot = sum(TPAREMV_UNADJ * tAdj * tDI, na.rm = TRUE),
                #prevPop = tPlot + mPlot * first(REMPER) + hPlot * first(REMPER) - rPlot * first(REMPER),
                #lPlot = (tPlot / ifelse(prevPop < 1, NA, prevPop)) ^ (1/first(REMPER)) - 1,
                EXPNS = first(EXPNS),
                #lPlot = (((tPlot / (tPlot + mPlot * first(REMPER) + hPlot * first(REMPER) - rPlot * first(REMPER))) ^ (1/first(REMPER))) - 1) * first(EXPNS),
                # lPlot = ((sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE) /
                #             ((sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE) +
                #                 sum(TPAMORT_UNADJ * tAdj * tDI, na.rm = TRUE) * first(REMPER) +
                #                 sum(TPAMORT_UNADJ * tAdj * tDI, na.rm = TRUE) * first(REMPER) -
                #                 sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'], na.rm = TRUE) * first(REMPER))) ^
                #             (1/first(REMPER))) - 1),
                plotIn_t = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0))
    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      #filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn_a = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    suppressMessages({
      t <- tInt %>%
        inner_join(aInt) %>%
        # group_by(.dots = grpBy, ESTN_UNIT_CN, STRATUM_CN) %>%
        # summarize(lStrat = sum(lPlot, na.rm = TRUE) * first(EXPNS))# %>%
        # group_by(.dots = grpBy, ESTN_UNIT_CN) %>%
        # summarize(lEst = sum(lStrat, na.rm = TRUE))
        ## Full region
        group_by(.dots = grpBy) %>%
        summarize(TREE_TOTAL = sum(tPlot * EXPNS, na.rm = TRUE),
                  RECR_TOTAL = sum(rPlot* EXPNS, na.rm = TRUE),
                  MORT_TOTAL = sum(mPlot* EXPNS, na.rm = TRUE),
                  REMV_TOTAL = sum(hPlot* EXPNS, na.rm = TRUE),
                  AREA_TOTAL = sum(fa, na.rm = TRUE),
                  #LAMBDA = sum(lPlot * EXPNS, na.rm = TRUE),
                  RECR_TPA = RECR_TOTAL / AREA_TOTAL,
                  MORT_TPA = MORT_TOTAL / AREA_TOTAL,
                  REMV_TPA = REMV_TOTAL / AREA_TOTAL,
                  RECR_PERC = RECR_TOTAL / TREE_TOTAL * 100,
                  MORT_PERC = MORT_TOTAL / TREE_TOTAL * 100,
                  REMV_PERC = REMV_TOTAL / TREE_TOTAL * 100,
                  # Non-zero plots
                  nPlots_TREE = sum(plotIn_t, na.rm = TRUE),
                  nPlots_AREA = sum(plotIn_a, na.rm = TRUE))
    })

    # Make some columns go away
    if (totals) {
      t <- t %>%
        select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC, #LAMBDA,
               TREE_TOTAL, RECR_TOTAL, MORT_TOTAL, REMV_TOTAL, AREA_TOTAL,
               nPlots_TREE, nPlots_AREA)
    } else {
      t <- t %>%
        select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC, #LAMBDA,
               nPlots_TREE, nPlots_AREA)
    }
  } # End SE Conditional

  # Some cleanup
  gc()

  # Return t
  t
}
