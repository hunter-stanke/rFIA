diversityHelper <- function(x, combos, data, grpBy, SE){
  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  #ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(tObs))) == length(tObs)) tObs <- 1
    td <- data$tDI * tObs * td
  }


  # IF we want sampling errors returned
  if(SE){
    data$tDI <- td
    data$tDI[is.na(data$tDI)] <- 0
    ## We produce an intermediate object in this chain as it is needed to compute the ratio of means variance
    ## Numerator and denominator are in different domains of interest, and may be grouped by different variables
    ## see covariance estimation below

    # Diversity is computed at the stand (condition level), and we continue to use the ratio of means estimator to get at
    #  average of the attribute of interest weighted by the area in which it occurs.
    d <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      # filter(tDI > 0) %>%
      summarize(condArea = first(CONDPROP_UNADJ),
                hCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'H') * condArea,
                sCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S') * condArea,
                EhCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh') * condArea,
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT),
                aAdj = first(aAdj),
                tAdj = first(tAdj),
                test = length(unique(SPCD)),
                nstems = sum(tDI, na.rm = TRUE))  %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(hPlot = sum(hCond * tAdj, na.rm = TRUE),
                EhPlot = sum(EhCond * tAdj, na.rm = TRUE),
                sPlot = sum(sCond * tAdj * plotIn, na.rm = TRUE),
                fa = sum(condArea * aDI * aAdj, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = first(a),
                p1EU = first(p1EU),
                p1 = first(p1),
                p2 = first(p2)) %>%
      # Continue through totals
      #d <- dInt %>%
      #filter(plotIn > 0) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                hStrat = mean(hPlot, na.rm = TRUE),
                EhStrat = mean(EhPlot, na.rm = TRUE),
                sStrat = mean(sPlot, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = first(a),
                w = first(p1) / first(p1EU), # Stratum weight
                nh = first(p2), # Number plots in stratum
                # Strata level variances
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * first(a) / nh),
                            (sum(fa^2, na.rm = TRUE) - sum(nh * aStrat^2, na.rm = TRUE)) / (nh * (nh-1))),
                hv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(hPlot * first(a) / nh),
                            (sum(hPlot^2, na.rm = TRUE) - sum(nh * hStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                Ehv = ifelse(first(ESTN_METHOD == 'simple'),
                             var(EhPlot * first(a) / nh),
                             (sum(EhPlot^2, na.rm = TRUE) - sum(nh * EhStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                sv = ifelse(first(ESTN_METHOD == 'simple'),
                            var(sPlot * first(a) / nh),
                            (sum(sPlot^2, na.rm = TRUE) - sum(nh * sStrat^2)) / (nh * (nh-1))),
                cvStrat_h = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(fa,hPlot),
                                   (sum(fa*hPlot) - sum(nh * aStrat *hStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_Eh = ifelse(first(ESTN_METHOD == 'simple'),
                                    cov(fa,EhPlot),
                                    (sum(fa*EhPlot) - sum(nh * aStrat *EhStrat)) / (nh * (nh-1))), # Stratified and double cases
                cvStrat_s = ifelse(first(ESTN_METHOD == 'simple'),
                                   cov(fa,sPlot),
                                   (sum(fa*sPlot) - sum(nh * aStrat *sStrat)) / (nh * (nh-1))))  %>% # Stratified and double cases
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                h = unitMean(ESTN_METHOD, a, nh, w, hStrat),
                eh = unitMean(ESTN_METHOD, a, nh, w, EhStrat),
                s = unitMean(ESTN_METHOD, a, nh, w, sStrat),
                plotIn = sum(plotIn, na.rm = TRUE),
                # Estimation of unit variance
                hVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, hv, hStrat, h),
                ehVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, Ehv, EhStrat, eh),
                sVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, sv, sStrat, s),
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
                cvEst_h = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_h, hStrat, h, aStrat, aEst),
                cvEst_eh = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_Eh, EhStrat, eh, aStrat, aEst),
                cvEst_s = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_s, sStrat, s, aStrat, aEst)) %>%
      # Compute totals
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                H_a = sum(h, na.rm = TRUE) / AREA_TOTAL,
                Eh_a = sum(eh, na.rm = TRUE) / AREA_TOTAL,
                S_a = sum(s, na.rm = TRUE) / AREA_TOTAL,
                hVar = sum(hVar, na.rm = TRUE),
                ehVar = sum(ehVar, na.rm = TRUE),
                sVar = sum(sVar, na.rm = TRUE),
                nStands = sum(plotIn, na.rm = TRUE),
                areaVar = sum(aVar, na.rm = TRUE),
                cvH = sum(cvEst_h, na.rm = TRUE),
                cveH = sum(cvEst_eh, na.rm = TRUE),
                cvS = sum(cvEst_s, na.rm = TRUE),
                hVar = (1/AREA_TOTAL^2) * (hVar + (H_a^2 * areaVar) - 2 * H_a * cvH),
                ehVar = (1/AREA_TOTAL^2) * (ehVar + (Eh_a^2 * areaVar) - 2 * Eh_a * cveH),
                sVar = (1/AREA_TOTAL^2) * (sVar + (S_a^2 * areaVar) - 2 * S_a * cvS),
                H_a_SE = sqrt(hVar) / H_a * 100,
                Eh_a_SE = sqrt(ehVar) / Eh_a * 100,
                S_a_SE = sqrt(sVar) / S_a * 100) %>%
      select(H_a, Eh_a, S_a, H_a_SE, Eh_a_SE, S_a_SE, nStands)

    # Beta & gamma diversity indices
    dbg <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by() %>%
      summarize(H_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'H'),
                H_b = H_g - d$H_a,
                Eh_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh'),
                Eh_b = Eh_g - d$Eh_a,
                S_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S'),
                S_b = S_g - d$S_a)

    # Join up the alpha beta gamma
    d <- data.frame(d, dbg) %>%
      select(H_a, H_b, H_g, Eh_a, Eh_b, Eh_g, S_a, S_b, S_g, H_a_SE, Eh_a_SE, S_a_SE, nStands)

    # Rejoin with groupby
    d <- data.frame(combos[[x]], d)

  } else { # No sampling errors
    ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
    d <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      summarize(aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0),
                condArea = first(CONDPROP_UNADJ),
                hCond = divIndex(SPCD, TPA_UNADJ  * tDI, index = 'H') * condArea,
                sCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S') * condArea,
                EhCond = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh') * condArea,
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                EXPNS = first(EXPNS),
                tAdj = first(tAdj),
                aAdj = first(aAdj)) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(condArea * aDI * aAdj * EXPNS, na.rm = TRUE),
                hPlot = sum(hCond * EXPNS * tAdj, na.rm = TRUE),
                EhPlot = sum(EhCond * EXPNS * tAdj, na.rm = TRUE),
                sPlot = sum(sCond * EXPNS * tAdj * plotIn, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE)) %>%
      group_by(.dots = grpBy) %>%
      summarize(AREA_TOTAL = sum(fa, na.rm = TRUE),
                H_a = sum(hPlot, na.rm = TRUE),
                Eh_a = sum(EhPlot, na.rm = TRUE),
                S_a = sum(sPlot, na.rm = TRUE),
                nStands = sum(plotIn, na.rm = TRUE)) #%>%
    #filter(S > 0) #%>%
    #select(c(grpByOrig, 'H_a', 'Eh_a', 'S_a', 'nStands'))

    # Beta & gamma diversity indices
    dbg <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy) %>%
      summarize(H_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'H'),
                Eh_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'Eh'),
                S_g = divIndex(SPCD, TPA_UNADJ * tDI, index = 'S'))
    # Join up the alpha beta gamma
    suppressMessages({
      d <- inner_join(d, dbg) %>%
        mutate(H_b = H_g - H_a,
               Eh_b = Eh_g - Eh_a,
               S_b = S_g - S_a) %>%
        ungroup() %>%
        select(grpBy, H_a, H_b, H_g, Eh_a, Eh_b, Eh_g, S_a, S_b, S_g, nStands)
    })
  }

  # Do some cleanup
  gc()

  #Return a dataframe
  d

}
