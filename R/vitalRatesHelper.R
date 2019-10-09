vitalRatesHelper <- function(x, combos, data, grpBy, aGrpBy, totals, SE, chngAdj){
  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(tObs))) == length(tObs)) tObs <- 1
    td <-  tObs * td * data$tDI
    # Area domain indicator for each column in
    if(grpBy[n] %in% aGrpBy){
      aObs <- as.character(combos[[x]][[aGrpBy[n]]]) == as.character(data[[aGrpBy[n]]])
      if (length(which(is.na(aObs))) == length(aObs)) aObs <- 1
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
    ### Compute total TREES in domain of interest

    tInt <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, TRE_CN, ONEORTWO, .keep_all = TRUE) %>%
      ## Plot level estimates
      # ## Omitting columns where either time is an NA
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE),
                dPlot = sum(DIA * tAdj * tDI, na.rm = TRUE),
                bPlot = sum(BA * tAdj * tDI, na.rm = TRUE),
                baaPlot = sum(TPAGROW_UNADJ * BA * tAdj * tDI, na.rm = TRUE),
                gPlot = sum(VOLCFNET * tAdj * tDI, na.rm = TRUE),
                gaPlot = sum(TPAGROW_UNADJ *VOLCFNET * tAdj * tDI, na.rm = TRUE),
                bioPlot = sum(DRYBIO_AG* tAdj * tDI / 2000, na.rm = TRUE),
                bioAPlot = sum(TPAGROW_UNADJ * DRYBIO_AG* tAdj * tDI / 2000, na.rm = TRUE),
                plotIn_t = ifelse(sum(tDI, na.rm = TRUE) > 0, 1,0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT))

    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(SUBPTYP_PROP_CHNG * chngAdj * aDI * aAdj, na.rm = TRUE),
                plotIn_a = ifelse(sum(aDI, na.rm = TRUE) >  0, 1,0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT))

    ## Compute COVARIANCE between numerator and denominator (for ratio estimates of variance)
    t <- tInt %>%
      #inner_join(aInt, by = c('PLT_CN', 'ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN'), suffix = c('_t', '_a')) %>%
      left_join(aInt, by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', 'PLT_CN'), suffix = c('_t', '_a'))  %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(tStrat = mean(tPlot, na.rm = TRUE),
                dStrat = mean(dPlot, na.rm = TRUE),
                bStrat = mean(bPlot, na.rm = TRUE),
                baaStrat = mean(baaPlot, na.rm = TRUE),
                #htStrat = mean(htPlot, na.rm = TRUE),
                gStrat = mean(gPlot, na.rm = TRUE),
                gaStrat = mean(gaPlot, na.rm = TRUE),
                #sStrat = mean(sPlot, na.rm = TRUE),
                #saStrat = mean(saPlot, na.rm = TRUE),
                bioStrat = mean(bioPlot, na.rm = TRUE),
                bioaStrat = mean(bioAPlot, na.rm = TRUE),
                #carbStrat = mean(carbPlot, na.rm = TRUE),
                #carbaStrat = mean(carbAPlot, na.rm = TRUE),
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
                dv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(dPlot * first(a) / nh),
                            (sum(dPlot^2) - sum(nh * dStrat^2)) / (nh * (nh-1))),
                bv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(bPlot * first(a) / nh),
                            (sum(bPlot^2) - sum(nh * bStrat^2)) / (nh * (nh-1))),
                baav = ifelse(first(ESTN_METHOD == 'simple'),
                              var(baaPlot * first(a) / nh),
                              (sum(baaPlot^2) - sum(nh * baaStrat^2)) / (nh * (nh-1))),
                # htv = ifelse(first(ESTN_METHOD == 'simple'),
                #              var(htPlot * first(a) / nh),
                #              (sum(htPlot^2) - sum(nh * htStrat^2)) / (nh * (nh-1))),
                gv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(gPlot * first(a) / nh),
                            (sum(gPlot^2) - sum(nh * gStrat^2)) / (nh * (nh-1))),
                gav = ifelse(first(ESTN_METHOD == 'simple'),
                             var(gaPlot * first(a) / nh),
                             (sum(gaPlot^2) - sum(nh * gaStrat^2)) / (nh * (nh-1))),
                # sv = ifelse(first(ESTN_METHOD == 'simple'),
                #             var(sPlot * first(a) / nh),
                #             (sum(sPlot^2) - sum(nh * sStrat^2)) / (nh * (nh-1))),
                # sav = ifelse(first(ESTN_METHOD == 'simple'),
                #              var(saPlot * first(a) / nh),
                #              (sum(saPlot^2) - sum(nh * saStrat^2)) / (nh * (nh-1))),
                biov = ifelse(first(ESTN_METHOD == 'simple'),
                              var(bioPlot * first(a) / nh),
                              (sum(bioPlot^2) - sum(nh * bioStrat^2)) / (nh * (nh-1))),
                bioav = ifelse(first(ESTN_METHOD == 'simple'),
                               var(bioAPlot * first(a) / nh),
                               (sum(bioAPlot^2) - sum(nh * bioaStrat^2)) / (nh * (nh-1))),
                # carbv = ifelse(first(ESTN_METHOD == 'simple'),
                #             var(carbPlot * first(a) / nh),
                #             (sum(carbPlot^2) - sum(nh * carbStrat^2)) / (nh * (nh-1))),
                # carbav = ifelse(first(ESTN_METHOD == 'simple'),
                #              var(carbAPlot * first(a) / nh),
                #              (sum(carbAPlot^2) - sum(nh * carbaStrat^2)) / (nh * (nh-1))),
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * first(a) / nh),
                            (sum(fa^2) - sum(nh * aStrat^2)) / (nh * (nh-1))),
                # Strata level covariances
                cvStrat_d = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(tPlot,dPlot),
                                   (sum(tPlot*dPlot) - sum(nh * tStrat *dStrat)) / (nh * (nh-1))), # Stratified and double casesc
                cvStrat_b = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(tPlot,bPlot),
                                   (sum(tPlot*bPlot) - sum(nh * tStrat *bStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_baa = ifelse(first(ESTN_METHOD == 'simple'),
                                     cov(fa,baaPlot),
                                     (sum(fa*baaPlot) - sum(nh * aStrat *baaStrat)) / (nh * (nh-1))),
                # cvStrat_ht = ifelse(first(ESTN_METHOD == 'simple'),
                #                     cov(tPlot,htPlot),
                #                     (sum(tPlot*htPlot) - sum(nh * tStrat *htStrat)) / (nh * (nh-1))), # Stratified and double case
                cvStrat_g = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(tPlot,gPlot),
                                   (sum(tPlot*gPlot) - sum(nh * tStrat *gStrat)) / (nh * (nh-1))),
                cvStrat_ga = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(fa,gaPlot),
                                    (sum(fa*gaPlot) - sum(nh * aStrat *gaStrat)) / (nh * (nh-1))),
                # cvStrat_s = ifelse(first(ESTN_METHOD == 'simple'),
                #                    cov(tPlot,sPlot),
                #                    (sum(tPlot*sPlot) - sum(nh * tStrat *sStrat)) / (nh * (nh-1))),
                # cvStrat_sa = ifelse(first(ESTN_METHOD == 'simple'),
                #                     cov(fa,saPlot),
                #                     (sum(fa*saPlot) - sum(nh * aStrat *saStrat)) / (nh * (nh-1))),
                cvStrat_bio = ifelse(first(ESTN_METHOD == 'simple'),
                                     cov(tPlot,bioPlot),
                                     (sum(tPlot*bioPlot) - sum(nh * tStrat *bioStrat)) / (nh * (nh-1))),
                cvStrat_bioa = ifelse(first(ESTN_METHOD == 'simple'),
                                      cov(fa,bioAPlot),
                                      (sum(fa*bioAPlot) - sum(nh * aStrat *bioaStrat)) / (nh * (nh-1))),
                # cvStrat_carb = ifelse(first(ESTN_METHOD == 'simple'),
                #                    cov(tPlot,carbPlot),
                #                    (sum(tPlot*carbPlot) - sum(nh * tStrat *carbStrat)) / (nh * (nh-1))),
                # cvStrat_carba = ifelse(first(ESTN_METHOD == 'simple'),
                #                     cov(fa,carbAPlot),
                #                     (sum(fa*carbAPlot) - sum(nh * aStrat *carbaStrat)) / (nh * (nh-1)))
      ) %>%
      # Estimation Unit
      group_by(ESTN_UNIT_CN) %>%
      summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
                dEst = unitMean(ESTN_METHOD, a, nh, w, dStrat),
                bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
                baaEst = unitMean(ESTN_METHOD, a, nh, w, baaStrat),
                #htEst = unitMean(ESTN_METHOD, a, nh, w, htStrat),
                gEst = unitMean(ESTN_METHOD, a, nh, w, gStrat),
                gaEst = unitMean(ESTN_METHOD, a, nh, w, gaStrat),
                # sEst = unitMean(ESTN_METHOD, a, nh, w, sStrat),
                # saEst = unitMean(ESTN_METHOD, a, nh, w, saStrat),
                bioEst = unitMean(ESTN_METHOD, a, nh, w, bioStrat),
                bioaEst = unitMean(ESTN_METHOD, a, nh, w, bioaStrat),
                # carbEst = unitMean(ESTN_METHOD, a, nh, w, carbStrat),
                # carbaEst = unitMean(ESTN_METHOD, a, nh, w, carbaStrat),
                aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
                nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE),
                # Variance estimates
                tVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, tv, tStrat, tEst),
                dVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, dv, dStrat, dEst),
                bVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bv, bStrat, bEst),
                baaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, baav, baaStrat, baaEst),
                #htVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, htv, htStrat, htEst),
                gVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, gv, gStrat, gEst),
                gaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, gav, gaStrat, gaEst),
                # sVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, sv, sStrat, sEst),
                # saVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, sav, sStrat, saEst),
                bioVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, biov, bioStrat, bioEst),
                bioaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bioav, bioaStrat, bioaEst),
                # carbVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, carbv, carbStrat, carbEst),
                # carbaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, carbav, carbaStrat, carbaEst),
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
                # Covariance estimates
                cvEst_d = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_d, dStrat, dEst, tStrat, tEst),
                cvEst_b = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_b, bStrat, bEst, tStrat, tEst),
                cvEst_baa = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_baa, baaStrat, baaEst, aStrat, aEst),
                #cvEst_ht = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_ht, htStrat, htEst, tStrat, tEst),
                cvEst_g = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_g, gStrat, gEst, tStrat, tEst),
                cvEst_ga = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_ga, gaStrat, gaEst, aStrat, aEst),
                # cvEst_s = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_s, sStrat, sEst, tStrat, tEst),
                # cvEst_sa = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_sa, saStrat, saEst, aStrat, aEst),
                cvEst_bio = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_bio, bioStrat, bioEst, tStrat, tEst),
                cvEst_bioa = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_bioa, bioaStrat, bioaEst, aStrat, aEst),
                # cvEst_carb = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_carb, carbStrat, carbEst, tStrat, tEst),
                # cvEst_carba = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_carba, carbaStrat, carbaEst, aStrat, aEst)
      ) %>%
      ## Full region
      summarize(TREE_TOTAL = sum(tEst, na.rm = TRUE),
                DIA_TOTAL = sum(dEst, na.rm = TRUE),
                BA_TOTAL = sum(bEst, na.rm = TRUE),
                BAA_TOTAL = sum(baaEst, na.rm = TRUE),
                #HT_TOTAL = sum(htEst, na.rm = TRUE),
                NETVOL_TOTAL = sum(gEst, na.rm = TRUE),
                NETVOL_AC_TOT = sum(gaEst, na.rm = TRUE),
                # SAWVOL_TOTAL = sum(sEst, na.rm = TRUE),
                # SAWVOL_AC_TOT = sum(saEst, na.rm = TRUE),
                BIO_TOTAL = sum(bioEst, na.rm = TRUE),
                BIO_AC_TOT = sum(bioaEst, na.rm = TRUE),
                # CARB_TOTAL = sum(carbEst, na.rm = TRUE),
                # CARB_AC_TOT = sum(carbaEst, na.rm = TRUE),
                AREA_TOTAL = sum(aEst, na.rm = TRUE),
                DIA_GROW = DIA_TOTAL / TREE_TOTAL,
                BA_GROW = BA_TOTAL / TREE_TOTAL,
                BAA_GROW = BAA_TOTAL / AREA_TOTAL,
                #HT_GROW = HT_TOTAL / TREE_TOTAL,
                NETVOL_GROW = NETVOL_TOTAL / TREE_TOTAL,
                NETVOL_GROW_AC = NETVOL_AC_TOT / AREA_TOTAL,
                # SAWVOL_GROW = SAWVOL_TOTAL / TREE_TOTAL,
                # SAWVOL_GROW_AC = SAWVOL_AC_TOT / AREA_TOTAL,
                BIO_GROW = BIO_TOTAL / TREE_TOTAL,
                BIO_GROW_AC = BIO_AC_TOT / AREA_TOTAL,
                # CARB_GROW = CARB_TOTAL / TREE_TOTAL,
                # CARB_GROW_AC = CARB_AC_TOT / AREA_TOTAL,
                # Variance/covariance
                treeVar = sum(tVar, na.rm = TRUE),
                dVar = sum(dVar, na.rm = TRUE),
                bVar = sum(bVar, na.rm = TRUE),
                baaVar = sum(baaVar, na.rm = TRUE),
                #htVar = sum(htVar, na.rm = TRUE),
                gVar = sum(gVar, na.rm = TRUE),
                gaVar = sum(gaVar, na.rm = TRUE),
                # sVar = sum(sVar, na.rm = TRUE),
                # saVar = sum(saVar, na.rm = TRUE),
                bioVar = sum(bioVar, na.rm = TRUE),
                bioaVar = sum(bioaVar, na.rm = TRUE),
                # carbVar = sum(carbVar, na.rm = TRUE),
                # carbaVar = sum(carbaVar, na.rm = TRUE),
                aVar = sum(aVar, na.rm = TRUE),
                cvD = sum(cvEst_d, na.rm = TRUE),
                cvB = sum(cvEst_b, na.rm = TRUE),
                cvBAA = sum(cvEst_baa, na.rm = TRUE),
                #cvHT = sum(cvEst_ht, na.rm = TRUE),
                cvG = sum(cvEst_g, na.rm = TRUE),
                cvGA = sum(cvEst_ga, na.rm = TRUE),
                # cvS = sum(cvEst_s, na.rm = TRUE),
                # cvSA = sum(cvEst_sa, na.rm = TRUE),
                cvBIO = sum(cvEst_bio, na.rm = TRUE),
                cvBIOA = sum(cvEst_bioa, na.rm = TRUE),
                # cvCARB = sum(cvEst_carb, na.rm = TRUE),
                # cvCARBA = sum(cvEst_carba, na.rm = TRUE),
                dgVar = (1/TREE_TOTAL^2) * (dVar + (DIA_GROW^2 * treeVar) - 2 * DIA_GROW * cvD),
                bgVar = (1/TREE_TOTAL^2) * (bVar + (BA_GROW^2 * treeVar) - 2 * BA_GROW * cvB),
                baagVar = (1/AREA_TOTAL^2) * (baaVar + (BAA_GROW^2 * aVar) - 2 * DIA_GROW * cvBAA),
                # htgVar = (1/TREE_TOTAL^2) * (htVar + (HT_GROW^2 * treeVar) - 2 * HT_GROW * cvHT),
                ggVar = (1/TREE_TOTAL^2) * (gVar + (NETVOL_GROW^2 * treeVar) - 2 * NETVOL_GROW * cvG),
                gagVar = (1/AREA_TOTAL^2) * (gaVar + (NETVOL_GROW_AC^2 * aVar) - 2 * NETVOL_GROW_AC * cvGA),
                # sgVar = (1/TREE_TOTAL^2) * (sVar + (SAWVOL_GROW^2 * treeVar) - 2 * SAWVOL_GROW * cvS),
                # sagVar = (1/AREA_TOTAL^2) * (saVar + (SAWVOL_GROW_AC^2 * aVar) - 2 * SAWVOL_GROW_AC * cvSA),
                biogVar = (1/TREE_TOTAL^2) * (bioVar + (BIO_GROW^2 * treeVar) - 2 * BIO_GROW * cvBIO),
                bioagVar = (1/AREA_TOTAL^2) * (bioaVar + (BIO_GROW_AC^2 * aVar) - 2 * BIO_GROW_AC * cvBIOA),
                # carbgVar = (1/TREE_TOTAL^2) * (carbVar + (CARB_GROW^2 * treeVar) - 2 * CARB_GROW * cvCARB),
                # carbagVar = (1/AREA_TOTAL^2) * (carbaVar + (CARB_GROW_AC^2 * aVar) - 2 * CARB_GROW_AC * cvCARBA),
                # Sampling Errors
                TREE_TOTAL_SE = sqrt(treeVar) / TREE_TOTAL * 100,
                DIA_TOTAL_SE = sqrt(dVar) / DIA_TOTAL * 100,
                BA_TOTAL_SE = sqrt(bVar) / BA_TOTAL * 100,
                BAA_TOTAL_SE = sqrt(baaVar) / BA_TOTAL * 100,
                # HT_TOTAL_SE = sqrt(htVar) / HT_TOTAL * 100,
                NETVOL_TOTAL_SE = sqrt(gVar) / NETVOL_TOTAL * 100,
                NETVOL_AC_TOT_SE = sqrt(gaVar) / NETVOL_AC_TOT * 100,
                AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
                DIA_GROW_SE = sqrt(dgVar) / DIA_GROW * 100,
                BA_GROW_SE = sqrt(bgVar) / BA_GROW * 100,
                BAA_GROW_SE = sqrt(baagVar) / BAA_GROW * 100,
                # HT_GROW_SE = sqrt(htgVar) / HT_GROW * 100,
                NETVOL_GROW_SE = sqrt(ggVar) / NETVOL_GROW * 100,
                NETVOL_GROW_AC_SE = sqrt(gagVar) / NETVOL_GROW_AC * 100,
                # SAWVOL_GROW_SE = sqrt(sgVar) / SAWVOL_GROW * 100,
                # SAWVOL_GROW_AC_SE = sqrt(sagVar) / SAWVOL_GROW_AC * 100,
                BIO_GROW_SE = sqrt(biogVar) / BIO_GROW * 100,
                BIO_GROW_AC_SE = sqrt(bioagVar) / BIO_GROW_AC * 100,
                # CARB_GROW_SE = sqrt(carbgVar) / CARB_GROW * 100,
                # CARB_GROW_AC_SE = sqrt(carbagVar) / CARB_GROW_AC * 100,
                # Non-zero plots
                nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
                nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE))

    # Make some columns go away
    if (totals) {
      t <- t %>%
        select(DIA_GROW, BA_GROW, BAA_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               BIO_GROW, BIO_GROW_AC,
               TREE_TOTAL, DIA_TOTAL, BA_TOTAL, BAA_TOTAL, NETVOL_TOTAL,
               NETVOL_AC_TOT, BIO_TOTAL, BIO_AC_TOT, AREA_TOTAL, DIA_GROW_SE, BA_GROW_SE, BAA_GROW_SE, NETVOL_GROW_SE,
               NETVOL_GROW_AC_SE, BIO_GROW_SE, BIO_GROW_AC_SE,
               nPlots_TREE, nPlots_AREA)
    } else {
      t <- t %>%
        select(DIA_GROW, BA_GROW, BAA_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               BIO_GROW, BIO_GROW_AC,
               DIA_GROW_SE, BA_GROW_SE, BAA_GROW_SE, NETVOL_GROW_SE,
               BIO_GROW_SE, BIO_GROW_AC_SE,
               NETVOL_GROW_AC_SE, nPlots_TREE, nPlots_AREA)
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
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * EXPNS * tDI, na.rm = TRUE),
                dPlot = sum(ANN_DIA_GROWTH * tAdj * EXPNS * tDI, na.rm = TRUE),
                bPlot = sum(ANN_BA_GROWTH * tAdj *  EXPNS *tDI, na.rm = TRUE),
                baaPlot = sum(TPAGROW_UNADJ * ANN_BA_GROWTH * EXPNS * tAdj * tDI, na.rm = TRUE),
                htPlot = sum(ANN_HT_GROWTH * tAdj *  EXPNS *tDI, na.rm = TRUE),
                gPlot = sum(ANN_NET_GROWTH * tAdj * EXPNS * tDI, na.rm = TRUE),
                gaPlot = sum(ANN_NET_GROWTH * TPAGROW_UNADJ *tAdj * EXPNS * tDI, na.rm = TRUE),
                plotIn_t = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0))
    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj *EXPNS, na.rm = TRUE),
                plotIn_a = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    suppressMessages({
      t <- tInt %>%
        inner_join(aInt) %>%
        ## Full region
        group_by(.dots = grpBy) %>%
        summarize(TREE_TOTAL = sum(tPlot, na.rm = TRUE),
                  DIA_TOTAL = sum(dPlot, na.rm = TRUE),
                  BA_TOTAL = sum(bPlot, na.rm = TRUE),
                  BAA_TOTAL = sum(baaPlot, na.rm = TRUE),
                  HT_TOTAL = sum(htPlot, na.rm = TRUE),
                  NETVOL_TOTAL = sum(gPlot, na.rm = TRUE),
                  NETVOL_AC_TOT = sum(gaPlot, na.rm = TRUE),
                  AREA_TOTAL = sum(fa, na.rm = TRUE),
                  DIA_GROW = DIA_TOTAL / TREE_TOTAL,
                  BA_GROW = BA_TOTAL / TREE_TOTAL,
                  BAA_GROW = BAA_TOTAL / AREA_TOTAL,
                  HT_GROW = HT_TOTAL / TREE_TOTAL,
                  NETVOL_GROW = NETVOL_TOTAL / TREE_TOTAL,
                  NETVOL_GROW_AC = NETVOL_AC_TOT / AREA_TOTAL,
                  # Non-zero plots
                  nPlots_TREE = sum(plotIn_t, na.rm = TRUE),
                  nPlots_AREA = sum(plotIn_a, na.rm = TRUE)) %>%
        ungroup()
    })

    if (totals) {
      t <- t %>%
        select(grpBy, DIA_GROW, BA_GROW, BAA_GROW, HT_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               DIA_TOTAL, BA_TOTAL, BAA_TOTAL, HT_TOTAL, NETVOL_TOTAL, NETVOL_AC_TOT,
               TREE_TOTAL, AREA_TOTAL, nPlots_TREE, nPlots_AREA)
    } else {
      t <- t %>%
        select(grpBy, DIA_GROW, BA_GROW, BAA_GROW, HT_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               nPlots_TREE, nPlots_AREA)
    }

  } # End SE Conditional

  # Some cleanup
  gc()

  # Return t
  t
}

vitalRatesHelper_old <- function(x, combos, data, grpBy, aGrpBy, totals, SE){
  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(tObs))) == length(tObs)) tObs <- 1
    td <- data$tDI * tObs * td
    # Area domain indicator for each column in
    if(grpBy[n] %in% aGrpBy){
      aObs <- as.character(combos[[x]][[aGrpBy[n]]]) == as.character(data[[aGrpBy[n]]])
      if (length(which(is.na(aObs))) == length(aObs)) aObs <- 1
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
    ### Compute total TREES in domain of interest
    tInt <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      #filter(EVALID %in% tID) %>%
      #filter(EVAL_TYP == 'EXPGROW') %>%
      # Compute estimates at plot level
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE),
                dPlot = sum(ANN_DIA_GROWTH * tAdj * tDI, na.rm = TRUE),
                bPlot = sum(ANN_BA_GROWTH * tAdj * tDI, na.rm = TRUE),
                baaPlot = sum(TPAGROW_UNADJ * ANN_BA_GROWTH * tAdj * tDI, na.rm = TRUE),
                htPlot = sum(ANN_HT_GROWTH * tAdj * tDI, na.rm = TRUE),
                gPlot = sum(ANN_NET_GROWTH * tAdj * tDI, na.rm = TRUE),
                gaPlot = sum(ANN_NET_GROWTH * TPAGROW_UNADJ *tAdj * tDI, na.rm = TRUE),
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
                dStrat = mean(dPlot, na.rm = TRUE),
                bStrat = mean(bPlot, na.rm = TRUE),
                baaStrat = mean(baaPlot, na.rm = TRUE),
                htStrat = mean(htPlot, na.rm = TRUE),
                gStrat = mean(gPlot, na.rm = TRUE),
                gaStrat = mean(gaPlot, na.rm = TRUE),
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
                dv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(dPlot * first(a) / nh),
                            (sum(dPlot^2) - sum(nh * dStrat^2)) / (nh * (nh-1))),
                bv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(bPlot * first(a) / nh),
                            (sum(bPlot^2) - sum(nh * bStrat^2)) / (nh * (nh-1))),
                baav = ifelse(first(ESTN_METHOD == 'simple'),
                              var(baaPlot * first(a) / nh),
                              (sum(baaPlot^2) - sum(nh * baaStrat^2)) / (nh * (nh-1))),
                htv = ifelse(first(ESTN_METHOD == 'simple'),
                             var(htPlot * first(a) / nh),
                             (sum(htPlot^2) - sum(nh * htStrat^2)) / (nh * (nh-1))),
                gv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(gPlot * first(a) / nh),
                            (sum(gPlot^2) - sum(nh * gStrat^2)) / (nh * (nh-1))),
                gav = ifelse(first(ESTN_METHOD == 'simple'),
                             var(gaPlot * first(a) / nh),
                             (sum(gaPlot^2) - sum(nh * gaStrat^2)) / (nh * (nh-1))),
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * first(a) / nh),
                            (sum(fa^2) - sum(nh * aStrat^2)) / (nh * (nh-1))),
                # Strata level covariances
                cvStrat_d = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(tPlot,dPlot),
                                   (sum(tPlot*dPlot) - sum(nh * tStrat *dStrat)) / (nh * (nh-1))), # Stratified and double casesc
                cvStrat_b = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(tPlot,bPlot),
                                   (sum(tPlot*bPlot) - sum(nh * tStrat *bStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_baa = ifelse(first(ESTN_METHOD == 'simple'),
                                     cov(fa,baaPlot),
                                     (sum(fa*baaPlot) - sum(nh * aStrat *baaStrat)) / (nh * (nh-1))),
                cvStrat_ht = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(tPlot,htPlot),
                                    (sum(tPlot*htPlot) - sum(nh * tStrat *htStrat)) / (nh * (nh-1))), # Stratified and double case
                cvStrat_g = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(tPlot,gPlot),
                                   (sum(tPlot*gPlot) - sum(nh * tStrat *gStrat)) / (nh * (nh-1))),
                cvStrat_ga = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(fa,gaPlot),
                                    (sum(fa*gaPlot) - sum(nh * aStrat *gaStrat)) / (nh * (nh-1)))) %>%
      # Estimation Unit
      group_by(ESTN_UNIT_CN) %>%
      summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
                dEst = unitMean(ESTN_METHOD, a, nh, w, dStrat),
                bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
                baaEst = unitMean(ESTN_METHOD, a, nh, w, baaStrat),
                htEst = unitMean(ESTN_METHOD, a, nh, w, htStrat),
                gEst = unitMean(ESTN_METHOD, a, nh, w, gStrat),
                gaEst = unitMean(ESTN_METHOD, a, nh, w, gaStrat),
                aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
                nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE),
                # Variance estimates
                tVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, tv, tStrat, tEst),
                dVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, dv, dStrat, dEst),
                bVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bv, bStrat, bEst),
                baaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, baav, baaStrat, baaEst),
                htVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, htv, htStrat, htEst),
                gVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, gv, gStrat, gEst),
                gaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, gav, tStrat, gaEst),
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
                # Covariance estimates
                cvEst_d = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_d, dStrat, dEst, tStrat, tEst),
                cvEst_b = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_b, bStrat, bEst, tStrat, tEst),
                cvEst_baa = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_baa, baaStrat, baaEst, aStrat, aEst),
                cvEst_ht = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_ht, htStrat, htEst, tStrat, tEst),
                cvEst_g = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_g, gStrat, gEst, tStrat, tEst),
                cvEst_ga = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_ga, gaStrat, gaEst, tStrat, tEst)) %>%
      ## Full region
      summarize(TREE_TOTAL = sum(tEst, na.rm = TRUE),
                DIA_TOTAL = sum(dEst, na.rm = TRUE),
                BA_TOTAL = sum(bEst, na.rm = TRUE),
                BAA_TOTAL = sum(baaEst, na.rm = TRUE),
                HT_TOTAL = sum(htEst, na.rm = TRUE),
                NETVOL_TOTAL = sum(gEst, na.rm = TRUE),
                NETVOL_AC_TOT = sum(gaEst, na.rm = TRUE),
                AREA_TOTAL = sum(aEst, na.rm = TRUE),
                DIA_GROW = DIA_TOTAL / TREE_TOTAL,
                BA_GROW = BA_TOTAL / TREE_TOTAL,
                BAA_GROW = BAA_TOTAL / AREA_TOTAL,
                HT_GROW = HT_TOTAL / TREE_TOTAL,
                NETVOL_GROW = NETVOL_TOTAL / TREE_TOTAL,
                NETVOL_GROW_AC = NETVOL_AC_TOT / AREA_TOTAL,
                # Variance/covariance
                treeVar = sum(tVar, na.rm = TRUE),
                dVar = sum(dVar, na.rm = TRUE),
                bVar = sum(bVar, na.rm = TRUE),
                baaVar = sum(baaVar, na.rm = TRUE),
                htVar = sum(htVar, na.rm = TRUE),
                gVar = sum(gVar, na.rm = TRUE),
                gaVar = sum(gaVar, na.rm = TRUE),
                aVar = sum(aVar, na.rm = TRUE),
                cvD = sum(cvEst_d, na.rm = TRUE),
                cvB = sum(cvEst_b, na.rm = TRUE),
                cvBAA = sum(cvEst_baa, na.rm = TRUE),
                cvHT = sum(cvEst_ht, na.rm = TRUE),
                cvG = sum(cvEst_g, na.rm = TRUE),
                cvGA = sum(cvEst_ga, na.rm = TRUE),
                dgVar = (1/TREE_TOTAL^2) * (dVar + (DIA_GROW^2 * treeVar) - 2 * DIA_GROW * cvD),
                bgVar = (1/TREE_TOTAL^2) * (bVar + (BA_GROW^2 * treeVar) - 2 * BA_GROW * cvB),
                baagVar = (1/AREA_TOTAL^2) * (baaVar + (BAA_GROW^2 * aVar) - 2 * DIA_GROW * cvBAA),
                htgVar = (1/TREE_TOTAL^2) * (htVar + (HT_GROW^2 * treeVar) - 2 * HT_GROW * cvHT),
                ggVar = (1/TREE_TOTAL^2) * (gVar + (NETVOL_GROW^2 * treeVar) - 2 * NETVOL_GROW * cvG),
                gagVar = (1/AREA_TOTAL^2) * (gaVar + (NETVOL_GROW_AC^2 * aVar) - 2 * NETVOL_GROW_AC * cvGA),
                # Sampling Errors
                TREE_TOTAL_SE = sqrt(treeVar) / TREE_TOTAL * 100,
                DIA_TOTAL_SE = sqrt(dVar) / DIA_TOTAL * 100,
                BA_TOTAL_SE = sqrt(bVar) / BA_TOTAL * 100,
                BAA_TOTAL_SE = sqrt(baaVar) / BA_TOTAL * 100,
                HT_TOTAL_SE = sqrt(htVar) / HT_TOTAL * 100,
                NETVOL_TOTAL_SE = sqrt(gVar) / NETVOL_TOTAL * 100,
                NETVOL_AC_TOT_SE = sqrt(gaVar) / NETVOL_AC_TOT * 100,
                AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
                DIA_GROW_SE = sqrt(dgVar) / DIA_GROW * 100,
                BA_GROW_SE = sqrt(bgVar) / BA_GROW * 100,
                BAA_GROW_SE = sqrt(baagVar) / BAA_GROW * 100,
                HT_GROW_SE = sqrt(htgVar) / HT_GROW * 100,
                NETVOL_GROW_SE = sqrt(ggVar) / NETVOL_GROW * 100,
                NETVOL_GROW_AC_SE = sqrt(gagVar) / NETVOL_GROW_AC * 100,
                # Non-zero plots
                nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
                nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE))

    # Make some columns go away
    if (totals) {
      t <- t %>%
        select(-c(names(t)[str_detect(names(t), 'Var') & str_detect(names(t), 'cv')]))
    } else {
      t <- t %>%
        select(DIA_GROW, BA_GROW, BAA_GROW, HT_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               DIA_GROW_SE, BA_GROW_SE, BAA_GROW_SE, HT_GROW_SE, NETVOL_GROW_SE,
               NETVOL_GROW_AC_SE, nPlots_TREE, nPlots_AREA)
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
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * EXPNS * tDI, na.rm = TRUE),
                dPlot = sum(ANN_DIA_GROWTH * tAdj * EXPNS * tDI, na.rm = TRUE),
                bPlot = sum(ANN_BA_GROWTH * tAdj *  EXPNS *tDI, na.rm = TRUE),
                baaPlot = sum(TPAGROW_UNADJ * ANN_BA_GROWTH * EXPNS * tAdj * tDI, na.rm = TRUE),
                htPlot = sum(ANN_HT_GROWTH * tAdj *  EXPNS *tDI, na.rm = TRUE),
                gPlot = sum(ANN_NET_GROWTH * tAdj * EXPNS * tDI, na.rm = TRUE),
                gaPlot = sum(ANN_NET_GROWTH * TPAGROW_UNADJ *tAdj * EXPNS * tDI, na.rm = TRUE),
                plotIn_t = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0))
    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj *EXPNS, na.rm = TRUE),
                plotIn_a = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    suppressMessages({
      t <- tInt %>%
        inner_join(aInt) %>%
        ## Full region
        group_by(.dots = grpBy) %>%
        summarize(TREE_TOTAL = sum(tPlot, na.rm = TRUE),
                  DIA_TOTAL = sum(dPlot, na.rm = TRUE),
                  BA_TOTAL = sum(bPlot, na.rm = TRUE),
                  BAA_TOTAL = sum(baaPlot, na.rm = TRUE),
                  HT_TOTAL = sum(htPlot, na.rm = TRUE),
                  NETVOL_TOTAL = sum(gPlot, na.rm = TRUE),
                  NETVOL_AC_TOT = sum(gaPlot, na.rm = TRUE),
                  AREA_TOTAL = sum(fa, na.rm = TRUE),
                  DIA_GROW = DIA_TOTAL / TREE_TOTAL,
                  BA_GROW = BA_TOTAL / TREE_TOTAL,
                  BAA_GROW = BAA_TOTAL / AREA_TOTAL,
                  HT_GROW = HT_TOTAL / TREE_TOTAL,
                  NETVOL_GROW = NETVOL_TOTAL / TREE_TOTAL,
                  NETVOL_GROW_AC = NETVOL_AC_TOT / AREA_TOTAL,
                  # Non-zero plots
                  nPlots_TREE = sum(plotIn_t, na.rm = TRUE),
                  nPlots_AREA = sum(plotIn_a, na.rm = TRUE)) %>%
        ungroup()
    })

    if (totals) {
      t <- t %>%
        select(grpBy, DIA_GROW, BA_GROW, BAA_GROW, HT_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               DIA_TOTAL, BA_TOTAL, BAA_TOTAL, HT_TOTAL, NETVOL_TOTAL, NETVOL_AC_TOT,
               TREE_TOTAL, AREA_TOTAL, nPlots_TREE, nPlots_AREA)
    } else {
      t <- t %>%
        select(grpBy, DIA_GROW, BA_GROW, BAA_GROW, HT_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               nPlots_TREE, nPlots_AREA)
    }

  } # End SE Conditional

  # Some cleanup
  gc()

  # Return t
  t
}
