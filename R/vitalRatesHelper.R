vrHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

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

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'DIA', 'DRYBIO_AG', 'VOLCFNET', 'VOLCSNET')), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
    left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA', 'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
    left_join(select(db$TREE_GRM_BEGIN, c('TRE_CN', 'DIA', 'VOLCFNET', 'VOLCSNET', 'DRYBIO_AG')), by = c('TRE_CN'), suffix = c('', '.beg')) %>%
    #left_join(select(db$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'CONDID'), suffix = c('', '.subp')) %>%
    #left_join(select(db$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN.subp' = 'PLT_CN', 'PREVCOND.subp' = 'CONDID'), suffix = c('', '.chng')) %>%
    left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
    left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD', 'DIA',  'DRYBIO_AG', 'VOLCFNET', 'VOLCSNET')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
    mutate_if(is.factor,
              as.character) %>%
    mutate(aChng = ifelse(COND_STATUS_CD.prev == 1 & COND_STATUS_CD == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
           tChng = ifelse(COND_STATUS_CD.prev == 1 & COND_STATUS_CD == 1, 1, 0))

  #If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
  data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
  data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
  data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
  data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
  data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
  data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp #* data$aChng
  data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng


  ## Modify  attributes depending on component (mortality uses midpoint)
  data <- data %>%
    mutate(DIA2 = vrAttHelper(DIA, DIA.prev, DIA.mid, DIA.beg, COMPONENT, REMPER, 2) * tDI,
           DIA1 = vrAttHelper(DIA, DIA.prev, DIA.mid, DIA.beg, COMPONENT, REMPER, 1) * tDI,
           BA2 = vrAttHelper(basalArea(DIA), basalArea(DIA.prev), basalArea(DIA.mid), basalArea(DIA.beg), COMPONENT, REMPER, 2) * tDI,
           BA1 = vrAttHelper(basalArea(DIA), basalArea(DIA.prev), basalArea(DIA.mid), basalArea(DIA.beg), COMPONENT, REMPER, 1) * tDI,
           VOLCFNET2 = vrAttHelper(VOLCFNET, VOLCFNET.prev, VOLCFNET.mid, VOLCFNET.beg, COMPONENT, REMPER, 2) * tDI,
           VOLCFNET1 = vrAttHelper(VOLCFNET, VOLCFNET.prev, VOLCFNET.mid, VOLCFNET.beg, COMPONENT, REMPER, 1) * tDI,
           DRYBIO_AG2 = vrAttHelper(DRYBIO_AG, DRYBIO_AG.prev, DRYBIO_AG.mid, DRYBIO_AG.beg, COMPONENT, REMPER, 2) * tDI,
           DRYBIO_AG1 = vrAttHelper(DRYBIO_AG, DRYBIO_AG.prev, DRYBIO_AG.mid, DRYBIO_AG.beg, COMPONENT, REMPER, 1) * tDI) %>%
    # ## Omitting columns where either time is an NA
    # mutate(DIA2 = case_when(is.na(DIA1) ~ NA_real_, TRUE ~ DIA2),
    #        DIA1 = case_when(is.na(DIA2) ~ NA_real_, TRUE ~ DIA1),
    #        BA2 = case_when(is.na(BA1) ~ NA_real_, TRUE ~ BA2),
    #        BA1 = case_when(is.na(BA2) ~ NA_real_, TRUE ~ BA1),
    #        VOLCFNET2 = case_when(is.na(VOLCFNET1) ~ NA_real_, TRUE ~ VOLCFNET2),
    #        VOLCFNET1 = case_when(is.na(VOLCFNET2) ~ NA_real_, TRUE ~ VOLCFNET1),
    #        DRYBIO_AG2 = case_when(is.na(DRYBIO_AG1) ~ NA_real_, TRUE ~ DRYBIO_AG2),
    #        DRYBIO_AG1 = case_when(is.na(DRYBIO_AG2) ~ NA_real_, TRUE ~ DRYBIO_AG1)) %>%
    select(-c(DIA.mid, VOLCFNET.mid, VOLCSNET.mid, DRYBIO_AG.mid,
              DIA.prev, VOLCFNET.prev, VOLCSNET.prev, DRYBIO_AG.prev,
              DIA.beg, VOLCFNET.beg, VOLCSNET.beg, DRYBIO_AG.beg,
              DIA, VOLCFNET, VOLCSNET, DRYBIO_AG))

  ## Just what we need
  data <- data %>%
    select(PLT_CN, TRE_CN, SUBP, CONDID, TREE, tDI,
           grpP, grpC, grpT, TPAGROW_UNADJ, PROP_BASIS, SUBPTYP_GRM, PLOT_STATUS_CD,
           DIA2, DIA1, BA2, BA1, DRYBIO_AG2, DRYBIO_AG1, VOLCFNET2, VOLCFNET1, MEASYEAR) %>%
    ## Dropping NA columns
    #drop_na(SUBPTYP_PROP_CHNG) %>%
    ## Rearrange previous values as observations
    pivot_longer(cols = DIA2:VOLCFNET1,
                 names_to = c(".value", 'ONEORTWO'),
                 names_sep = -1)

  ### DOING AREA SEPARATELY NOW FOR GROWTH ACCOUNTING PLOTS
  aData <- select(db$PLOT,c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$SUBP_COND_CHNG_MTRX, PLT_CN, PREV_PLT_CN, SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID), by = c('PLT_CN', 'PREV_PLT_CN')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$COND, c('PLT_CN', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
    #left_join(select(db$POP_PLOT_STRATUM_ASSGN, by = 'PLT_CN')) %>%
    #left_join(select(db$POP_STRATUM, by = c('STRATUM_CN' = 'CN'))) %>%
    mutate(aChng = if_else(COND_STATUS_CD == 1 &
                             COND_STATUS_CD.prev == 1 &
                             !is.null(CONDPROP_UNADJ) &
                             SUBPTYP == 1,
                           1, 0),
           SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25)

  aData$aDI <- aData$landD * aData$aD_p * aData$aD_c * aData$sp * aData$aChng


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, SUBP, TREE) %>%
      summarize(d = sum(DIA * tDI, na.rm = TRUE),
                ba = sum(BA * tDI, na.rm = TRUE),
                baa = sum(TPAGROW_UNADJ * BA * tDI, na.rm = TRUE),
                vol = sum(VOLCFNET * tDI, na.rm = TRUE),
                volA = sum(TPAGROW_UNADJ * VOLCFNET * tDI, na.rm = TRUE),
                bio = sum(DRYBIO_AG * tDI, na.rm = TRUE),
                bioA = sum(TPAGROW_UNADJ * DRYBIO_AG * tDI, na.rm = TRUE)) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(DIA_GROW = mean(d, na.rm = TRUE),
                BA_GROW = mean(ba, na.rm = TRUE),
                BAA_GROW = sum(baa, na.rm = TRUE),
                NETVOL_GROW = mean(vol, na.rm = TRUE),
                NETVOL_GROW_AC = sum(volA, na.rm = TRUE),
                BIO_GROW = mean(vol, na.rm = TRUE),
                BIO_GROW_AC = sum(volA, na.rm = TRUE),
                nStems = length(which(!is.na(d))))

    a = NULL

  } else {
    ### Plot-level estimates -- growth accounting
    a <- aData %>%
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa = sum(SUBPTYP_PROP_CHNG * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))
    # ### Plot-level estimates
    # a <- data %>%
    #   ## Will be lots of trees here, so CONDPROP listed multiple times
    #   ## Adding PROP_BASIS so we can handle adjustment factors at strata level
    #   distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
    #   group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
    #   summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
    #             plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
    #   left_join(select(a_ga, PLT_CN, PROP_BASIS, aGrpBy, fa_ga, plotIn_ga), by = c('PLT_CN', 'PROP_BASIS', aGrpBy))

    t <- data %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, SUBPTYP_GRM) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ  * tDI, na.rm = TRUE),
                dPlot = sum(DIA  * tDI, na.rm = TRUE),
                bPlot = sum(BA  * tDI, na.rm = TRUE),
                baaPlot = sum(TPAGROW_UNADJ * BA  * tDI, na.rm = TRUE),
                gPlot = sum(VOLCFNET  * tDI, na.rm = TRUE),
                gaPlot = sum(TPAGROW_UNADJ *VOLCFNET  * tDI, na.rm = TRUE),
                bioPlot = sum(DRYBIO_AG * tDI / 2000, na.rm = TRUE),
                bioAPlot = sum(TPAGROW_UNADJ * DRYBIO_AG * tDI / 2000, na.rm = TRUE),
                plotIn_t = ifelse(sum(tDI, na.rm = TRUE) > 0, 1,0))
  }




  pltOut <- list(a = a, t = t)
  return(pltOut)

}



vrHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

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
    filter(EVAL_TYP %in% c('EXPGROW')) %>%
    #filter(EVAL_TYP %in% c('EXPCURR')) %>%
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
              av = stratVar(ESTN_METHOD, fa, aStrat, ndif, a, nh))
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
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    filter(EVAL_TYP %in% c('EXPGROW')) %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
    mutate(tAdj = case_when(
      ## When NA, stay NA
      is.na(SUBPTYP_GRM) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      SUBPTYP_GRM == 0 ~ 0,
      SUBPTYP_GRM == 1 ~ as.numeric(ADJ_FACTOR_SUBP),
      SUBPTYP_GRM == 2 ~ as.numeric(ADJ_FACTOR_MICR),
      SUBPTYP_GRM == 3 ~ as.numeric(ADJ_FACTOR_MACR)),
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
      dPlot = dPlot * tAdj,
      bPlot = bPlot * tAdj,
      baaPlot = baaPlot * tAdj,
      gPlot = gPlot * tAdj,
      gaPlot = gaPlot * tAdj,
      bioPlot = bioPlot * tAdj,
      bioAPlot = bioAPlot * tAdj) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(tPlot = sum(tPlot, na.rm = TRUE),
              dPlot = sum(dPlot, na.rm = TRUE),
              bPlot = sum(bPlot, na.rm = TRUE),
              baaPlot = sum(baaPlot, na.rm = TRUE),
              gPlot = sum(gPlot, na.rm = TRUE),
              gaPlot = sum(gaPlot, na.rm = TRUE),
              bioPlot = sum(bioPlot, na.rm = TRUE),
              bioAPlot = sum(bioAPlot, na.rm = TRUE),
              fa = first(fa),
              plotIn_TREE = ifelse(sum(plotIn_t >  0, na.rm = TRUE), 1,0),
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
              dStrat = mean(dPlot * r_t, na.rm = TRUE),
              bStrat = mean(bPlot * r_t, na.rm = TRUE),
              baaStrat = mean(baaPlot * r_t, na.rm = TRUE),
              gStrat = mean(gPlot * r_t, na.rm = TRUE),
              gaStrat = mean(gaPlot * r_t, na.rm = TRUE),
              bioStrat = mean(bioPlot * r_t, na.rm = TRUE),
              bioAStrat = mean(bioAPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              tv = stratVar(ESTN_METHOD, tPlot, tStrat, ndif, a, nh),
              dv = stratVar(ESTN_METHOD, dPlot, dStrat, ndif, a, nh),
              bv = stratVar(ESTN_METHOD, bPlot, bStrat, ndif, a, nh),
              baav = stratVar(ESTN_METHOD, baaPlot, baaStrat, ndif, a, nh),
              gv = stratVar(ESTN_METHOD, gPlot, gStrat, ndif, a, nh),
              gav = stratVar(ESTN_METHOD, gaPlot, gaStrat, ndif, a, nh),
              biov = stratVar(ESTN_METHOD, bioPlot, bioStrat, ndif, a, nh),
              bioAv = stratVar(ESTN_METHOD, bioAPlot, bioAStrat, ndif, a, nh),
              # Strata level covariances
              cvStrat_d = stratVar(ESTN_METHOD, dPlot, dStrat, ndif, a, nh, tPlot, tStrat),
              cvStrat_b = stratVar(ESTN_METHOD, bPlot, bStrat, ndif, a, nh, tPlot, tStrat),
              cvStrat_baa = stratVar(ESTN_METHOD, baaPlot, baaStrat, ndif, a, nh, fa, aStrat),
              cvStrat_g = stratVar(ESTN_METHOD, gPlot, gStrat, ndif, a, nh, tPlot, tStrat),
              cvStrat_ga = stratVar(ESTN_METHOD, gaPlot, gaStrat, ndif, a, nh, fa, aStrat),
              cvStrat_bio = stratVar(ESTN_METHOD, bioPlot, bioStrat, ndif, a, nh, tPlot, tStrat),
              cvStrat_bioA = stratVar(ESTN_METHOD, bioAPlot, bioAStrat, ndif, a, nh, fa, aStrat)

    ) %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
              dEst = unitMean(ESTN_METHOD, a, nh, w, dStrat),
              bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
              baaEst = unitMean(ESTN_METHOD, a, nh, w, baaStrat),
              gEst = unitMean(ESTN_METHOD, a, nh, w, gStrat),
              gaEst = unitMean(ESTN_METHOD, a, nh, w, gaStrat),
              bioEst = unitMean(ESTN_METHOD, a, nh, w, bioStrat),
              bioAEst = unitMean(ESTN_METHOD, a, nh, w, bioAStrat),
              #aEst = first(aEst),
              # Estimation of unit variance
              tVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, tv, tStrat, tEst),
              dVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, dv, dStrat, dEst),
              bVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bv, bStrat, bEst),
              baaVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, baav, baaStrat, baaEst),
              gVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, gv, gStrat, gEst),
              gaVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, gav, gaStrat, gaEst),
              bioVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, biov, bioStrat, bioEst),
              bioAVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bioAv, bioAStrat, bioAEst),

              cvEst_d = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_d, dStrat, dEst, tStrat, tEst),
              cvEst_b = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_b, bStrat, bEst, tStrat, tEst),
              cvEst_g = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_g, gStrat, gEst, tStrat, tEst),
              cvEst_bio = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_bio, bioStrat, bioEst, tStrat, tEst),
              cvEst_baa = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_baa, baaStrat, baaEst, aStrat, aEst),
              cvEst_ga = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ga, gaStrat, gaEst, aStrat, aEst),
              cvEst_bioA = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_bioA, bioAStrat, bioAEst, aStrat, aEst),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}


























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
      summarize(## Totals
                tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
                dEst = unitMean(ESTN_METHOD, a, nh, w, dStrat),
                bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
                baaEst = unitMean(ESTN_METHOD, a, nh, w, baaStrat),
                gEst = unitMean(ESTN_METHOD, a, nh, w, gStrat),
                gaEst = unitMean(ESTN_METHOD, a, nh, w, gaStrat),
                bioEst = unitMean(ESTN_METHOD, a, nh, w, bioStrat),
                bioaEst = unitMean(ESTN_METHOD, a, nh, w, bioaStrat),
                aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
                nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE),
                ## Variance estimates -- totals
                tVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, tv, tStrat, tEst),
                dVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, dv, dStrat, dEst),
                bVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bv, bStrat, bEst),
                baaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, baav, baaStrat, baaEst),
                gVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, gv, gStrat, gEst),
                gaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, gav, gaStrat, gaEst),
                bioVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, biov, bioStrat, bioEst),
                bioaVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bioav, bioaStrat, bioaEst),
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
                ## Covariance estimates -- ratios
                cvEst_d = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_d, dStrat, dEst, tStrat, tEst),
                cvEst_b = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_b, bStrat, bEst, tStrat, tEst),
                cvEst_baa = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_baa, baaStrat, baaEst, aStrat, aEst),
                cvEst_g = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_g, gStrat, gEst, tStrat, tEst),
                cvEst_ga = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_ga, gaStrat, gaEst, aStrat, aEst),
                cvEst_bio = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_bio, bioStrat, bioEst, tStrat, tEst),
                cvEst_bioa = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_bioa, bioaStrat, bioaEst, aStrat, aEst)
      ) %>%
      ## Full region
      summarize(TREE_TOTAL = sum(tEst, na.rm = TRUE),
                DIA_TOTAL = sum(dEst, na.rm = TRUE),
                BA_TOTAL = sum(bEst, na.rm = TRUE),
                BAA_TOTAL = sum(baaEst, na.rm = TRUE),
                NETVOL_TOTAL = sum(gEst, na.rm = TRUE),
                NETVOL_AC_TOT = sum(gaEst, na.rm = TRUE),
                BIO_TOTAL = sum(bioEst, na.rm = TRUE),
                BIO_AC_TOT = sum(bioaEst, na.rm = TRUE),
                AREA_TOTAL = sum(aEst, na.rm = TRUE),
                DIA_GROW = DIA_TOTAL / TREE_TOTAL,
                BA_GROW = BA_TOTAL / TREE_TOTAL,
                BAA_GROW = BAA_TOTAL / AREA_TOTAL,
                NETVOL_GROW = NETVOL_TOTAL / TREE_TOTAL,
                NETVOL_GROW_AC = NETVOL_AC_TOT / AREA_TOTAL,
                BIO_GROW = BIO_TOTAL / TREE_TOTAL,
                BIO_GROW_AC = BIO_AC_TOT / AREA_TOTAL,
                # Variance/covariance
                treeVar = sum(tVar, na.rm = TRUE),
                dVar = sum(dVar, na.rm = TRUE),
                bVar = sum(bVar, na.rm = TRUE),
                baaVar = sum(baaVar, na.rm = TRUE),
                gVar = sum(gVar, na.rm = TRUE),
                gaVar = sum(gaVar, na.rm = TRUE),
                bioVar = sum(bioVar, na.rm = TRUE),
                bioaVar = sum(bioaVar, na.rm = TRUE),
                aVar = sum(aVar, na.rm = TRUE),
                cvD = sum(cvEst_d, na.rm = TRUE),
                cvB = sum(cvEst_b, na.rm = TRUE),
                cvBAA = sum(cvEst_baa, na.rm = TRUE),
                cvG = sum(cvEst_g, na.rm = TRUE),
                cvGA = sum(cvEst_ga, na.rm = TRUE),
                cvBIO = sum(cvEst_bio, na.rm = TRUE),
                cvBIOA = sum(cvEst_bioa, na.rm = TRUE),
                dgVar = (1/TREE_TOTAL^2) * (dVar + (DIA_GROW^2 * treeVar) - 2 * DIA_GROW * cvD),
                bgVar = (1/TREE_TOTAL^2) * (bVar + (BA_GROW^2 * treeVar) - 2 * BA_GROW * cvB),
                baagVar = (1/AREA_TOTAL^2) * (baaVar + (BAA_GROW^2 * aVar) - 2 * DIA_GROW * cvBAA),
                ggVar = (1/TREE_TOTAL^2) * (gVar + (NETVOL_GROW^2 * treeVar) - 2 * NETVOL_GROW * cvG),
                gagVar = (1/AREA_TOTAL^2) * (gaVar + (NETVOL_GROW_AC^2 * aVar) - 2 * NETVOL_GROW_AC * cvGA),
                biogVar = (1/TREE_TOTAL^2) * (bioVar + (BIO_GROW^2 * treeVar) - 2 * BIO_GROW * cvBIO),
                bioagVar = (1/AREA_TOTAL^2) * (bioaVar + (BIO_GROW_AC^2 * aVar) - 2 * BIO_GROW_AC * cvBIOA),
                # Sampling Errors
                TREE_TOTAL_SE = sqrt(treeVar) / TREE_TOTAL * 100,
                DIA_TOTAL_SE = sqrt(dVar) / abs(DIA_TOTAL) * 100,
                BA_TOTAL_SE = sqrt(bVar) / abs(BA_TOTAL) * 100,
                BAA_TOTAL_SE = sqrt(baaVar) / abs(BA_TOTAL) * 100,
                # HT_TOTAL_SE = sqrt(htVar) / HT_TOTAL * 100,
                NETVOL_TOTAL_SE = sqrt(gVar) / abs(NETVOL_TOTAL) * 100,
                NETVOL_AC_TOT_SE = sqrt(gaVar) / abs(NETVOL_AC_TOT) * 100,
                AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
                DIA_GROW_SE = sqrt(dgVar) / abs(DIA_GROW) * 100,
                BA_GROW_SE = sqrt(bgVar) / abs(BA_GROW) * 100,
                BAA_GROW_SE = sqrt(baagVar) / abs(BAA_GROW) * 100,
                # HT_GROW_SE = sqrt(htgVar) / HT_GROW * 100,
                NETVOL_GROW_SE = sqrt(ggVar) / abs(NETVOL_GROW) * 100,
                NETVOL_GROW_AC_SE = sqrt(gagVar) / abs(NETVOL_GROW_AC) * 100,
                # SAWVOL_GROW_SE = sqrt(sgVar) / SAWVOL_GROW * 100,
                # SAWVOL_GROW_AC_SE = sqrt(sagVar) / SAWVOL_GROW_AC * 100,
                BIO_GROW_SE = sqrt(biogVar) / abs(BIO_GROW) * 100,
                BIO_GROW_AC_SE = sqrt(bioagVar) / abs(BIO_GROW_AC) * 100,
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
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, TRE_CN, ONEORTWO, .keep_all = TRUE) %>%
      #filter(EVALID %in% tID) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI * EXPNS, na.rm = TRUE),
                dPlot = sum(DIA * tAdj * tDI * EXPNS, na.rm = TRUE),
                bPlot = sum(BA * tAdj * tDI * EXPNS, na.rm = TRUE),
                baaPlot = sum(TPAGROW_UNADJ * BA * tAdj * tDI * EXPNS, na.rm = TRUE),
                gPlot = sum(VOLCFNET * tAdj * tDI * EXPNS, na.rm = TRUE),
                gaPlot = sum(TPAGROW_UNADJ *VOLCFNET * tAdj * tDI * EXPNS, na.rm = TRUE),
                bioPlot = sum(DRYBIO_AG* tAdj * tDI  * EXPNS/ 2000, na.rm = TRUE),
                bioAPlot = sum(TPAGROW_UNADJ * DRYBIO_AG* tAdj * tDI * EXPNS / 2000, na.rm = TRUE),
                plotIn_t = ifelse(sum(tDI, na.rm = TRUE) > 0, 1,0))
    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(SUBPTYP_PROP_CHNG * chngAdj * aDI * aAdj * EXPNS, na.rm = TRUE),
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
                  NETVOL_TOTAL = sum(gPlot, na.rm = TRUE),
                  NETVOL_AC_TOT = sum(gaPlot, na.rm = TRUE),
                  BIO_TOTAL = sum(bioPlot, na.rm = TRUE),
                  BIO_AC_TOT = sum(bioAPlot, na.rm = TRUE),
                  AREA_TOTAL = sum(fa, na.rm = TRUE),
                  DIA_GROW = DIA_TOTAL / TREE_TOTAL,
                  BA_GROW = BA_TOTAL / TREE_TOTAL,
                  BAA_GROW = BAA_TOTAL / AREA_TOTAL,
                  NETVOL_GROW = NETVOL_TOTAL / TREE_TOTAL,
                  NETVOL_GROW_AC = NETVOL_AC_TOT / AREA_TOTAL,
                  BIO_GROW = BIO_TOTAL / TREE_TOTAL,
                  BIO_GROW_AC = BIO_TOTAL / AREA_TOTAL,
                  # Non-zero plots
                  nPlots_TREE = sum(plotIn_t, na.rm = TRUE),
                  nPlots_AREA = sum(plotIn_a, na.rm = TRUE)) %>%
        ungroup()
    })

    # Make some columns go away
    if (totals) {
      t <- t %>%
        select(grpBy, DIA_GROW, BA_GROW, BAA_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               BIO_GROW, BIO_GROW_AC,
               TREE_TOTAL, DIA_TOTAL, BA_TOTAL, BAA_TOTAL, NETVOL_TOTAL,
               NETVOL_AC_TOT, BIO_TOTAL, BIO_AC_TOT, AREA_TOTAL,
               nPlots_TREE, nPlots_AREA)
    } else {
      t <- t %>%
        select(grpBy, DIA_GROW, BA_GROW, BAA_GROW, NETVOL_GROW, NETVOL_GROW_AC,
               BIO_GROW, BIO_GROW_AC, nPlots_TREE, nPlots_AREA)
    }

  } # End SE Conditional

  # Some cleanup
  #gc()

  # Return t
  t
}

