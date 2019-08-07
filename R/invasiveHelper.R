invasiveHelper <- function(x, combos, data, grpBy, aGrpBy, totals, SE){

  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    td <- data$aDI * tObs * td
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
    invPrep <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SYMBOL, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(iPlot = ifelse(any(INVASIVE_SAMPLING_STATUS_CD == 1),
                               sum(COVER_PCT/100 * CONDPROP_UNADJ * aAdj * tDI, na.rm = TRUE),
                               NA),
                iPlots = ifelse(sum(tDI > 0 & iPlot > 0, na.rm = TRUE), 1, 0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT)) #%>%
    inv <- invPrep %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(iStrat = mean(iPlot, na.rm = TRUE),
                iPlots = sum(iPlots, na.rm = TRUE),
                a = first(a),
                w = first(p1)/first(p1EU),
                nh = first(p2),
                iv = ifelse(first(ESTN_METHOD == "simple"),
                            var(iPlot * first(a)/iPlots),
                            (sum(iPlot^2, na.rm = TRUE) - sum(iPlots * iStrat^2))/(iPlots * (iPlots - 1)))) %>%
      group_by(ESTN_UNIT_CN) %>%
      summarize(iEst = unitMean(ESTN_METHOD, a, iPlots, w, iStrat),
                iPlots = sum(iPlots, na.rm = TRUE),
                iVar = unitVar(method = "var", ESTN_METHOD, a, iPlots, w, iv, iStrat, iEst)) %>%
      summarize(COVER_AREA = sum(iEst, na.rm = TRUE),
                iVar = sum(iVar, na.rm = TRUE),
                COVER_AREA_SE = sqrt(iVar)/COVER_AREA * 100,
                nPlots_INV = sum(iPlots, na.rm = TRUE))

    aPrep <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
                plotIn = ifelse(sum(aDI > 0, na.rm = TRUE), 1, 0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT)) #%>%
    a <- aPrep %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = first(a),
                w = first(p1)/first(p1EU),
                nh = first(p2),
                av = ifelse(first(ESTN_METHOD == "simple"),
                            var(fa * first(a)/nh),
                            (sum(fa^2) - sum(nh * aStrat^2))/(nh * (nh - 1)))) %>%
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                plotIn = sum(plotIn, na.rm = TRUE),
                aVar = unitVar(method = "var", ESTN_METHOD, a, nh, w, av, aStrat, aEst)) %>%
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                aVar = sum(aVar, na.rm = TRUE),
                AREA_TOTAL_SE = sqrt(aVar)/AREA_TOTAL * 100,
                nPlots_AREA = sum(plotIn, na.rm = TRUE))


    ## Compute COVARIANCE between numerator and denominator (for ratio estimates of variance)
    covar <- invPrep %>%
      inner_join(aPrep, by = c('PLT_CN', 'ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN'), suffix = c('_i', '_a')) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                iStrat = mean(iPlot, na.rm = TRUE),
                a = first(a_i),
                w = first(p1_i) / first(p1EU_a), # Stratum weight
                nh = first(p2_i), # Number plots in stratum
                # Strata level covariances
                cvStrat = ifelse(first(ESTN_METHOD == "simple"),
                                 cov(fa, iPlot),
                                 (sum(fa * iPlot, na.rm = TRUE) - sum(nh * aStrat * iStrat))/(nh * (nh - 1)))) %>% # Stratified and double cases)
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                iEst = unitMean(ESTN_METHOD, a, iPlots, w, iStrat),
                cvEst = unitVar(method = "cov", ESTN_METHOD, a, nh, w, cvStrat, iStrat, iEst, aStrat, aEst)) %>%
      summarize(cv = sum(cvEst, na.rm = TRUE))

    ## Join up tree, area, and covariance tables
    inv <- data.frame(inv, a, covar)

    ## Make TPA, BAA and SE
    inv <- inv %>%
      mutate(COVER_PCT = COVER_AREA / AREA_TOTAL * 100,
             coverVar = (1/AREA_TOTAL^2) * (iVar + (COVER_PCT^2 * aVar) - 2 * COVER_PCT * cv),
             COVER_PCT_SE = sqrt(coverVar)/COVER_PCT * 100)

    if (totals) {
      inv <- inv %>%
        select(COVER_PCT, COVER_AREA, AREA_TOTAL, COVER_PCT_SE, COVER_AREA_SE, AREA_TOTAL_SE, nPlots_INV, nPlots_AREA)
    }
    else {
      inv <- inv %>%
        select(COVER_PCT, COVER_PCT_SE, nPlots_INV, nPlots_AREA)
    }

    inv <- data.frame(combos[[x]], inv) %>%
    filter(!is.na(SYMBOL))

  } else { # IF SE is FALSE
    a <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn = ifelse(sum(aDI > 0, na.rm = TRUE), 1, 0)) %>%
      group_by(.dots = aGrpBy) %>%
      summarize(AREA_TOTAL = sum(fa, na.rm = TRUE),
                nPlots_AREA = sum(plotIn, na.rm = TRUE))
    # Used to make a new EXPNS value below, number of plots where invasives were sampled
    nPlots <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, STRATUM_CN, PLT_CN) %>%
      summarize(plotIn = ifelse(sum(INVASIVE_SAMPLING_STATUS_CD == 1, na.rm = TRUE), 1, 0)) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, STRATUM_CN) %>%
      summarize(plotsInv = sum(plotIn, na.rm = TRUE))

    inv <- data %>%
      filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SYMBOL, .keep_all = TRUE) %>%
      inner_join(nPlots, by = c(aGrpBy, 'ESTN_UNIT_CN', 'STRATUM_CN')) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(expInv = first(AREA_USED)*first(P1POINTCNT) / first(P1PNTCNT_EU) / first(plotsInv),
                iPlot = ifelse(any(INVASIVE_SAMPLING_STATUS_CD == 1),
                               sum(COVER_PCT/100 * CONDPROP_UNADJ * aAdj * expInv * aDI, na.rm = TRUE),
                               NA)) %>%
      group_by(.dots = grpBy) %>%
      summarize(COVER_AREA = sum(iPlot, na.rm = TRUE),
                nPlots_INV = length(which(!is.na(iPlot))))


    # Rejoin
    suppressMessages({
      inv <- inner_join(inv, a) %>%
        filter(!is.na(SYMBOL)) %>%
        mutate(COVER_PCT = COVER_AREA / AREA_TOTAL * 100)
    })



    if (totals) {
      inv <- inv %>%
        select(grpBy, COVER_PCT, COVER_AREA, AREA_TOTAL, nPlots_INV, nPlots_AREA)
    }
    else {
      inv <- inv %>%
        select(grpBy, COVER_PCT, nPlots_INV, nPlots_AREA)
    }

  } # End SE conditional

  # Do some cleanup
  gc()

  #Return a dataframe
  inv
}
