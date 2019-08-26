areaHelper <- function(x, combos, data, grpBy, totals, SE){
  # Update domain indicator for each each column speficed in grpBy
  td = 1 # Start both at 1, update as we iterate through
  ad = 1

  for (n in 1:ncol(combos[[x]])){
    # Tree domain indicator for each column in
    tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(tObs))) == length(tObs)) tObs <- 1
    td <- data$aDI * tObs * td
  }

  if(SE){
    data$aDI <- td
    data$aDI[is.na(data$aDI)] <- 0

    # Compute estimates
    a <- data %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, CONDID) %>%
      summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
                aDI = ifelse(sum(aDI > 0, na.rm = TRUE), 1, 0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT),
                aAdj = first(aAdj)) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
                p1EU = first(p1EU),
                a = first(a),
                p1 = first(p1),
                p2 = first(p2)) %>%
      # Continue through totals
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = first(a),
                w = first(p1) / first(p1EU), # Stratum weight
                nh = first(p2), # Number plots in stratum
                # Strata level variances
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * first(a) / nh),
                            (sum(fa^2, na.rm = TRUE) - sum(nh * aStrat^2, na.rm = TRUE)) / (nh * (nh-1)))) %>% # Stratified and double cases) %>% # Stratified and double cases
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                plotIn = sum(plotIn, na.rm = TRUE),
                # Estimation of unit variance,
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst)) %>%
      # Compute totals
      summarize(AREA = sum(aEst, na.rm = TRUE),
                nPlots = sum(plotIn, na.rm = TRUE),
                areaVar = sum(aVar, na.rm = TRUE),
                AREA_SE = ifelse(nPlots > 1, sqrt(areaVar) / AREA * 100,0)) %>%
      select(AREA, AREA_SE, nPlots)

    # Rejoin with some grpBy Names
    a <- data.frame(combos[[x]], a)

  } else { # No sampling errors
    ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
    # Compute estimates
    a <- data %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, EXPNS, PLT_CN, CONDID) %>%
      summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ),
                aDI = ifelse(sum(aDI > 0, na.rm = TRUE), 1, 0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT),
                aAdj = first(aAdj)) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      group_by(.dots = grpBy) %>%
      summarize(AREA = sum(fa, na.rm = TRUE),
                nPlots = sum(plotIn, na.rm = TRUE)) %>%
      select(grpBy, AREA, nPlots)

  } # End SE Conditional

  # Some cleanup
  gc()

  # Return t
  a
}
