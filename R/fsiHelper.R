
fsiHelper1 <- function(x, plts, db, grpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)

  # Only subplots from cond change matrix
  #db$SUBP_COND_CHNG_MTRX <- filter(db$SUBP_COND_CHNG_MTRX, SUBPTYP == 1)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]

  ## Making a replica condition ID so that we can group by attributes of previous conditions
  ## PREV COND is an issue because of recruitment values (PREVCOND = NA)
  #db$COND <- mutate(db$COND, PREV_CONDID = CONDID)

  ## Making a treeID
  # db$TREE$treID <- paste(db$TREE$SUBP, db$TREE$TREE, sep = '_')

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'pltID', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR',
                            'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp', 'DESIGNCD')) %>%
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & DESIGNCD == 1 & PLOT_STATUS_CD != 3) %>%

    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    ## AGENTCD at remeasurement, died during the measurement interval
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'TPA_UNADJ', 'BAA', 'DIA', 'AGENTCD', 'MORTYR', 'STATUSCD', 'SPCD')), by = c('PLT_CN', 'CONDID')) %>%
    #left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', DIA_BEGIN, DIA_END)), by = c('TRE_CN')) %>%

    left_join(select(db$PLOT, c('PLT_CN', 'sp', 'aD_p', 'DESIGNCD', 'PLOT_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('2', '1')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('2', '1')) %>%
    left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD', 'TPA_UNADJ', 'BAA', 'DIA', 'STATUSCD', htClass, SPCD)), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%
    #left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%

    mutate_if(is.factor,
              as.character)

  ## Comprehensive indicator function -- w/ growth accounting
  data$tDI2 <- data$landD2 * data$aD_p2 * data$aD_c2 * data$tD2 * data$typeD2 * data$sp2 *
    if_else(data$PLOT_STATUS_CD1 == 1 & data$PLOT_STATUS_CD2 == 1, 1, 0) *
    if_else(data$STATUSCD2 == 1, 1, 0) *
    ## If it is NA, then it is a recruit, otherwise, make sure it was not
    ## flagged incorrectly previously
    if_else(is.na(data$STATUSCD1), 1, if_else(data$STATUSCD1 == 0 | data$STATUSCD2 == 0, 0, 1))

  data$tDI1 <- data$landD1 * data$aD_p1 * data$aD_c1 * data$tD1 * data$typeD1 * data$sp1 *
    if_else(data$PLOT_STATUS_CD1 == 1 & data$PLOT_STATUS_CD2 == 1, 1, 0) *
    if_else(data$STATUSCD1 == 1, 1, 0) *
    if_else(is.na(data$STATUSCD1), NA_real_, if_else(data$STATUSCD1 == 0 | data$STATUSCD2 == 0, 0, 1))

  ## Doesn't allow recruitment to happen, becuase no previous status
  ## Adjust later, after the pivot
  #if_else(data$STATUSCD1 == 0 | data$STATUSCD2 == 0, 0, 1)



  ## PREVIOUS and CURRENT attributes
  data <- data %>%
    mutate(TPA_UNADJ1 = TPA_UNADJ1,
           TPA_UNADJ2 = TPA_UNADJ2,
           BAA1 = BAA1,
           BAA2 = BAA2,
           MORT = case_when(
             STATUSCD1 == 1 & STATUSCD2 == 2 ~ 1,
             STATUSCD1 == 1 & STATUSCD2 == 3 ~ 1,
             TRUE ~ 0),
           SURV = case_when(
             STATUSCD1 == 1 & STATUSCD2 == 1 ~ 1,
             TRUE ~ 0)
           #DISTURB = if_else()
    ) #%>%

  ## Just what we need
  data <- data %>%
    select(PLT_CN, PREV_PLT_CN, pltID, TRE_CN, SUBP, CONDID, TREE, CONDPROP_UNADJ,
           MEASYEAR, MACRO_BREAKPOINT_DIA, PROP_BASIS, grpP[grpP != 'PLOT_STATUS_CD'], grpC,
           REMPER, PLOT_STATUS_CD1, PLOT_STATUS_CD2,
           #one_of(str_c(grpP,1), str_c(grpC,1), str_c(grpT,1),
           #        str_c(grpP,2), str_c(grpC,2), str_c(grpT,2)),
           one_of(str_c(grpT,1),str_c(grpT,2)),
           tDI1, tDI2, STATUSCD1, STATUSCD2,
           DIA1, DIA2, BAA1, BAA2, TPA_UNADJ1, TPA_UNADJ2, SPCD1, SPCD2) %>%
    mutate(BAA1 = -(BAA1),
           TPA_UNADJ1 = -(TPA_UNADJ1)) %>%
    ## Rearrange previous values as observations
    pivot_longer(cols = -c(PLT_CN:REMPER),
                 names_to = c(".value", 'ONEORTWO'),
                 names_sep = -1) %>%
    mutate(PLOT_BASIS = case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA) ~ NA_character_,
      ## When DIA is less than 5", use microplot value
      DIA < 5 ~ 'MICR',
      ## When DIA is greater than 5", use subplot value
      DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
      DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) #%>%
  #mutate(tDI = tDI * if_else())



  if (byPlot){
    grpBy <- c('YEAR', grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, PLT_CN, PREV_PLT_CN) %>%
      summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
                REMPER = first(REMPER),
                PREV_BAA = sum(-BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
                PREV_TPA = sum(-TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
                CHNG_TPA = sum(TPA_UNADJ[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
                #CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
                CHNG_BAA = mean((BAA[ONEORTWO == 2] * tDI[ONEORTWO == 2]) + (BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1]), na.rm = TRUE),
                CURR_TPA = PREV_TPA + CHNG_TPA,
                CURR_BAA = PREV_BAA + CHNG_BAA,
                #DIFF_BAA = mean((BAA[ONEORTWO == 2] * tDI[ONEORTWO == 2]) + (BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1]), na.rm = TRUE),
                nStems = length(which(tDI == 1))) %>%
      ungroup() %>%
      select(PLT_CN, PREV_PLT_CN, REMPER, grpBy,
             PREV_TPA, PREV_BAA, CHNG_TPA, CHNG_BAA, CURR_TPA, CURR_BAA,
             nStems, nLive)
    a = NULL

  } else {
    # ### Plot-level estimates -- growth accounting
    # a <- data %>%
    #   ## Will be lots of trees here, so CONDPROP listed multiple times
    #   ## Adding PROP_BASIS so we can handle adjustment factors at strata level
    #   distinct(PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
    #   group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
    #   summarize(fa = sum(SUBPTYP_PROP_CHNG * aDI, na.rm = TRUE),
    #             plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))
    ### Plot-level estimates
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      #distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, CONDID, .dots = grpBy) %>%
      summarize(aDI = if_else(sum(tDI, na.rm = TRUE) > 0, 1, 0),
                CONDPROP_UNADJ = first(CONDPROP_UNADJ)) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = grpBy) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE))

    ### Compute total TREES in domain of interest
    t <- data %>%
      ## Need to count number of trees in each
      mutate(treID = paste(pltID, SUBP, TREE)) %>%
      distinct(PLT_CN, TRE_CN, ONEORTWO, .keep_all = TRUE) %>%
      #group_by(PLT_CN, PLOT_BASIS, TRE_CN) %>%
      #summarize()
      # Compute estimates at plot level
      group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
      summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
                REMPER = first(REMPER),
                MEASYEAR = first(MEASYEAR),
                PREV_TPA = sum(-TPA_UNADJ[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
                PREV_BAA = sum(-BAA[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
                CHNG_TPA = sum(TPA_UNADJ[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
                CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
                #CHNG_BAA1 = sum((BAA[ONEORTWO == 2] * tDI[ONEORTWO == 2]) + (BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1]), na.rm = TRUE),
                plotIn = if_else(sum(tDI, na.rm = TRUE) > 0, 1, 0),
                n = length(unique(TRE_CN[tDI == 1])))
                #n1 = length(unique(treID[tDI == 1])))
  }

  pltOut <- list(a = a, t = t)
  return(pltOut)

}


fsiHelper2 <- function(x, popState, a, t, grpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR

  }

  ######## ------------------ TREE ESTIMATES + CV
  aAdj <- a %>%
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
      fa = fa * aAdj) %>%
    ungroup()

  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD, REMPER)), by = 'PLT_CN') %>%
    ## Need forest area to adjust SI indices
    left_join(select(aAdj, PLT_CN, grpBy, fa, aAdj), by = c('PLT_CN', grpBy)) %>%
    #Add adjustment factors
    mutate(tAdj = case_when(
      ## When NA, stay NA
      is.na(PLOT_BASIS) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      PLOT_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
      ## Otherwise, use the subpplot value
      PLOT_BASIS == 'SUBP' ~ as.numeric(ADJ_FACTOR_SUBP),
      PLOT_BASIS == 'MICR' ~ as.numeric(ADJ_FACTOR_MICR)),
      ct = CHNG_TPA * tAdj,
      cb = CHNG_BAA * tAdj,
      pt = PREV_TPA * tAdj,
      pb = PREV_BAA * tAdj) %>%
    ## Computing change
    mutate(ct = (ct) / REMPER,
           cb = (cb) / REMPER) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(ctPlot = sum(ct, na.rm = TRUE),
              cbPlot = sum(cb, na.rm = TRUE),
              ptPlot = sum(pt, na.rm = TRUE),
              pbPlot = sum(pb, na.rm = TRUE),
              plotIn_t = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              REMPER = first(REMPER),
              ## Summing change across subp and micr
              TPA_RATE = sum(TPA_RATE, na.rm = TRUE),
              BAA_RATE = sum(BAA_RATE, na.rm = TRUE),
              ## Total unique number of trees
              n = sum(n, na.rm = TRUE)) %>%
    ## Do not want to compute SI for micro and subp seperately, handle it here
    mutate(TPA_RATE = TPA_RATE / REMPER,
           BAA_RATE = BAA_RATE / REMPER / n,
           #CURR_TPA = ptPlot + (ctPlot * REMPER),
           #CURR_BAA = pbPlot + (cbPlot * REMPER),
           #TPA_RATE = ctPlot, #/ (CURR_TPA + ptPlot) * 2,
           #BAA_RATE = cbPlot, #/ (CURR_BAA + pbPlot) * 2,
           x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
           #y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
           siPlot = sqrt(x^2 + x^2),
           siPlot = if_else(x < 0, -siPlot, siPlot),
           siPlot = case_when(
             is.na(siPlot) ~ 0,
             TRUE ~ siPlot)) %>%

    left_join(select(aAdj, PLT_CN, grpBy, fa), by = c('PLT_CN', grpBy)) %>%

    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              ctStrat = mean(ctPlot * r_t, na.rm = TRUE),
              cbStrat = mean(cbPlot * r_t, na.rm = TRUE),
              ptStrat = mean(ptPlot * r_t, na.rm = TRUE),
              pbStrat = mean(pbPlot * r_t, na.rm = TRUE),
              siStrat = mean(siPlot * fa * r_t, na.rm = TRUE),
              faStrat = mean(fa * r_t, na.rm = TRUE),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              ctv = stratVar(ESTN_METHOD, ctPlot, ctStrat, ndif, a, nh),
              cbv = stratVar(ESTN_METHOD, cbPlot, cbStrat, ndif, a, nh),
              ptv = stratVar(ESTN_METHOD, ptPlot, ptStrat, ndif, a, nh),
              pbv = stratVar(ESTN_METHOD, pbPlot, pbStrat, ndif, a, nh),
              siv = stratVar(ESTN_METHOD, siPlot * fa, siStrat, ndif, a, nh),
              fav = stratVar(ESTN_METHOD, fa, faStrat, ndif, a, nh),


              # Strata level covariances
              cvStrat_ct = stratVar(ESTN_METHOD, ctPlot, ctStrat, ndif, a, nh, ptPlot, ptStrat),
              cvStrat_cb = stratVar(ESTN_METHOD, cbPlot, cbStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_si = stratVar(ESTN_METHOD, siPlot * fa, siStrat, ndif, a, nh, fa, faStrat)) %>%

    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
              cbEst = unitMean(ESTN_METHOD, a, nh, w, cbStrat),
              ptEst = unitMean(ESTN_METHOD, a, nh, w, ptStrat),
              pbEst = unitMean(ESTN_METHOD, a, nh, w, pbStrat),
              siEst = unitMean(ESTN_METHOD, a, nh, w, siStrat),
              faEst = unitMean(ESTN_METHOD, a, nh, w, faStrat),

              nh = first(nh),
              # Estimation of unit variance
              ctVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ctv, ctStrat, ctEst),
              cbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cbv, cbStrat, cbEst),
              ptVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ptv, ptStrat, ptEst),
              pbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pbv, pbStrat, pbEst),
              siVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, siv, siStrat, siEst),
              faVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, fav, faStrat, faEst),

              ## Covariances
              cvEst_ct = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ct, ctStrat, ctEst, ptStrat, ptEst),
              cvEst_cb = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_cb, cbStrat, cbEst, pbStrat, pbEst),
              cvEst_si = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_si, siStrat, siEst, faStrat, faEst),

              plotIn_t = sum(plotIn_t, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = NULL)

  return(out)
}
