vrHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot, treeType){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'DIA', 'DRYBIO_AG', 'VOLCFNET', 'VOLBFNET', 'STATUSCD')), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
    left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA', 'VOLCFNET', 'VOLBFNET', 'DRYBIO_AG')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
    left_join(select(db$TREE_GRM_BEGIN, c('TRE_CN', 'DIA', 'VOLCFNET', 'VOLBFNET', 'DRYBIO_AG')), by = c('TRE_CN'), suffix = c('', '.beg')) %>%
    left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
    left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD', 'DIA',  'DRYBIO_AG', 'VOLCFNET',  'VOLBFNET', 'STATUSCD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
    mutate_if(is.factor,
              as.character) %>%
    mutate(aChng = ifelse(COND_STATUS_CD.prev == 1 & COND_STATUS_CD == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
           tChng = ifelse(COND_STATUS_CD.prev == 1 & COND_STATUS_CD == 1, 1, 0),
           status = case_when(STATUSCD == 1 & STATUSCD.prev == 1 ~ 1,
                                                            TRUE ~ 0))

  #If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
  data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
  data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
  data$landD.prev <- ifelse(data$landD == 1 & data$landD.prev == 1, 1, 0)
    #ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
  data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
  data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
  data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

  ## Comprehensive indicator function
  data$aDI <- data$landD.prev * data$aD_p * data$aD_c * data$sp * data$aChng
  data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng

  ## Only if we're only considering live trees that stayed live
  if (tolower(treeType) == 'live') {data$tDI <- data$tDI * data$status}

  ## Modify  attributes depending on component (mortality uses midpoint)
  data <- data %>%
    mutate(DIA2 = vrAttHelper(DIA, DIA.prev, DIA.mid, DIA.beg, COMPONENT, REMPER, 2) * tDI,
           DIA1 = vrAttHelper(DIA, DIA.prev, DIA.mid, DIA.beg, COMPONENT, REMPER, 1) * tDI,
           BA2 = vrAttHelper(basalArea(DIA), basalArea(DIA.prev), basalArea(DIA.mid), basalArea(DIA.beg), COMPONENT, REMPER, 2) * tDI,
           BA1 = vrAttHelper(basalArea(DIA), basalArea(DIA.prev), basalArea(DIA.mid), basalArea(DIA.beg), COMPONENT, REMPER, 1) * tDI,
           VOLCFNET2 = vrAttHelper(VOLCFNET, VOLCFNET.prev, VOLCFNET.mid, VOLCFNET.beg, COMPONENT, REMPER, 2) * tDI,
           VOLCFNET1 = vrAttHelper(VOLCFNET, VOLCFNET.prev, VOLCFNET.mid, VOLCFNET.beg, COMPONENT, REMPER, 1) * tDI,
           VOLBFNET2 = vrAttHelper(VOLBFNET, VOLBFNET.prev, VOLBFNET.mid, VOLBFNET.beg, COMPONENT, REMPER, 2) * tDI,
           VOLBFNET1 = vrAttHelper(VOLBFNET, VOLBFNET.prev, VOLBFNET.mid, VOLBFNET.beg, COMPONENT, REMPER, 1) * tDI,
           DRYBIO_AG2 = vrAttHelper(DRYBIO_AG, DRYBIO_AG.prev, DRYBIO_AG.mid, DRYBIO_AG.beg, COMPONENT, REMPER, 2) * tDI,
           DRYBIO_AG1 = vrAttHelper(DRYBIO_AG, DRYBIO_AG.prev, DRYBIO_AG.mid, DRYBIO_AG.beg, COMPONENT, REMPER, 1) * tDI) %>%
    select(-c(DIA.mid, VOLCFNET.mid, VOLBFNET.mid, DRYBIO_AG.mid,
              DIA.prev, VOLCFNET.prev, VOLBFNET.prev, DRYBIO_AG.prev,
              DIA.beg, VOLCFNET.beg, VOLBFNET.beg, DRYBIO_AG.beg,
              DIA, VOLCFNET, VOLBFNET, DRYBIO_AG))

  ## Just what we need
  data <- data %>%
    select(PLT_CN, TRE_CN, SUBP, CONDID, TREE, tDI,
           grpP, grpC, grpT, TPAGROW_UNADJ, PROP_BASIS, SUBPTYP_GRM, PLOT_STATUS_CD,
           DIA2, DIA1, BA2, BA1, DRYBIO_AG2, DRYBIO_AG1, VOLCFNET2, VOLCFNET1, VOLBFNET2, VOLBFNET1, MEASYEAR) %>%
    ## Rearrange previous values as observations
    pivot_longer(cols = DIA2:VOLBFNET1,
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

  aData$aDI <- aData$landD * aData$landD.prev * aData$aD_p * aData$aD_c * aData$sp * aData$aChng


  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    #grpSyms1 <- syms(grpBy[grpBy %in% c('SUBP', 'TREE') == FALSE])
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(t = sum(TPAGROW_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), ## Previous only
                d = sum(DIA * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                ba = sum(BA * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                vol = sum(VOLCFNET * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                svol = sum(VOLBFNET * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                bio = sum(DRYBIO_AG * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                nStems = length(unique(TRE_CN))) %>%
      mutate(DIA_GROW = d / t,
             BA_GROW = ba / t,
             NETVOL_GROW = vol / t,
             SAWVOL_GROW = svol / t,
             BIO_GROW = bio / t,
             BAA_GROW = ba,
             NETVOL_GROW_AC = vol,
             SAWVOL_GROW_AC = svol,
             BIO_GROW_AC = bio,
             PREV_TPA = t) %>%
      select(-c(t:bio)) %>%
      as.data.frame() %>%
      relocate(nStems, .after = PREV_TPA)


    a = NULL

  } else {
    ### Plot-level estimates -- growth accounting
    a <- aData %>%
      #filter(SUBPTYP_GRM == 1) %>%
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa = sum(SUBPTYP_PROP_CHNG * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    grpSyms <- syms(grpBy)
    t <- data %>%
      lazy_dt() %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      group_by(!!!grpSyms, PLT_CN, SUBPTYP_GRM) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), ## Previous only
                dPlot = sum(DIA * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                bPlot = sum(BA * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                gPlot = sum(VOLCFNET * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                sPlot = sum(VOLBFNET * TPAGROW_UNADJ * tDI, na.rm = TRUE),
                bioPlot = sum(DRYBIO_AG * TPAGROW_UNADJ * tDI / 2000, na.rm = TRUE),
                plotIn_t = ifelse(sum(tDI, na.rm = TRUE) > 0, 1,0)) %>%
      as.data.frame()
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
    popState[[x]]$P1POINTCNT <- popState[[x]]$P1POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  aGrpSyms <- syms(aGrpBy)

  ## Strata level estimates
  aStrat <- a %>%
    lazy_dt() %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    filter(EVAL_TYP %in% c('EXPGROW')) %>%
    mutate(
      ## AREA
      aAdj = dplyr::case_when(
        ## When NA, stay NA
        is.na(PROP_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = fa * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!aGrpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              p2eu = dplyr::first(p2eu),
              aStrat = sum(fa, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
              ## Strata level variances
              av = sum(fa^2, na.rm = TRUE)) %>%
    mutate(aStrat = aStrat / nh,
           av = (av - (nh * aStrat^2)) / (nh * (nh-1))) %>%
    as.data.frame()
  ## Estimation unit
  aEst <- aStrat %>%
    group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, av, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

  ######## ------------------ TREE ESTIMATES + CV

  grpSyms <- syms(grpBy)
  ## Strata level estimates
  tEst <- t %>%
    lazy_dt() %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    filter(EVAL_TYP %in% c('EXPGROW')) %>%
    ## Need this for covariance later on
    left_join(select(a, fa, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
    #Add adjustment factors
    mutate(tAdj = dplyr::case_when(
      ## When NA, stay NA
      is.na(SUBPTYP_GRM) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      SUBPTYP_GRM == 0 ~ 0,
      SUBPTYP_GRM == 1 ~ as.numeric(ADJ_FACTOR_SUBP),
      SUBPTYP_GRM == 2 ~ as.numeric(ADJ_FACTOR_MICR),
      SUBPTYP_GRM == 3 ~ as.numeric(ADJ_FACTOR_MACR)),
      ## AREA
      aAdj = dplyr::case_when(
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
      gPlot = gPlot * tAdj,
      sPlot = sPlot * tAdj,
      bioPlot = bioPlot * tAdj) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, !!!grpSyms) %>%
    summarize(tPlot = sum(tPlot, na.rm = TRUE),
              dPlot = sum(dPlot, na.rm = TRUE),
              bPlot = sum(bPlot, na.rm = TRUE),
              gPlot = sum(gPlot, na.rm = TRUE),
              sPlot = sum(sPlot, na.rm = TRUE),
              bioPlot = sum(bioPlot, na.rm = TRUE),
              fa = dplyr::first(fa),
              plotIn_TREE = ifelse(sum(plotIn_t >  0, na.rm = TRUE), 1,0),
              nh = dplyr::first(P2POINTCNT),
              p2eu = dplyr::first(p2eu),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU)) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(nh),
              a = dplyr::first(a),
              w = dplyr::first(w),
              p2eu = dplyr::first(p2eu),

              ## dtplyr is fast, but requires a few extra steps, so we'll finish
              ## means and variances in subseqent mutate step

              tStrat = sum(tPlot , na.rm = TRUE),
              dStrat = sum(dPlot , na.rm = TRUE),
              bStrat = sum(bPlot , na.rm = TRUE),
              gStrat = sum(gPlot , na.rm = TRUE),
              sStrat = sum(sPlot , na.rm = TRUE),
              bioStrat = sum(bioPlot , na.rm = TRUE),
              aStrat = dplyr::first(aStrat),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE),
              # ## Strata level variances
              tv = sum(tPlot^2, na.rm = TRUE),
              dv = sum(dPlot^2, na.rm = TRUE),
              bv = sum(bPlot^2, na.rm = TRUE),
              gv = sum(gPlot^2, na.rm = TRUE),
              sv = sum(sPlot^2, na.rm = TRUE),
              biov = sum(bioPlot^2, na.rm = TRUE),
              # Strata level covariances
              cvStrat_d = sum(dPlot*tPlot, na.rm = TRUE),
              cvStrat_b = sum(bPlot*tPlot, na.rm = TRUE),
              cvStrat_baa = sum(bPlot*fa, na.rm = TRUE),
              cvStrat_g = sum(gPlot*tPlot, na.rm = TRUE),
              cvStrat_ga = sum(gPlot*fa, na.rm = TRUE),
              cvStrat_s = sum(sPlot*tPlot, na.rm = TRUE),
              cvStrat_sa = sum(sPlot*fa, na.rm = TRUE),
              cvStrat_bio = sum(bioPlot*tPlot, na.rm = TRUE),
              cvStrat_bioA = sum(bioPlot*fa, na.rm = TRUE)) %>%
    mutate(tStrat = tStrat / nh,
           dStrat = dStrat / nh,
           bStrat = bStrat / nh,
           gStrat = gStrat / nh,
           sStrat = sStrat / nh,
           bioStrat = bioStrat / nh,
           adj = nh * (nh-1),
           tv = (tv - (nh*tStrat^2)) / adj,
           dv = (dv - (nh*dStrat^2)) / adj,
           bv = (bv - (nh*bStrat^2)) / adj,
           gv = (gv - (nh*gStrat^2)) / adj,
           sv = (sv - (nh*sStrat^2)) / adj,
           biov = (biov - (nh*bioStrat^2)) / adj,
           cvStrat_d = (cvStrat_d - (nh * dStrat * tStrat)) / adj,
           cvStrat_b = (cvStrat_b - (nh * bStrat * tStrat)) / adj,
           cvStrat_baa = (cvStrat_baa - (nh * bStrat * aStrat)) / adj,
           cvStrat_g = (cvStrat_g - (nh * gStrat * tStrat)) / adj,
           cvStrat_ga = (cvStrat_ga - (nh * gStrat * aStrat)) / adj,
           cvStrat_s = (cvStrat_s - (nh * sStrat * tStrat)) / adj,
           cvStrat_sa = (cvStrat_sa - (nh * sStrat * aStrat)) / adj,
           cvStrat_bio = (cvStrat_bio - (nh * bioStrat * tStrat)) / adj,
           cvStrat_bioA = (cvStrat_bioA - (nh * bioStrat * aStrat)) / adj) %>%
    as.data.frame() %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
              dEst = unitMean(ESTN_METHOD, a, nh, w, dStrat),
              bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
              gEst = unitMean(ESTN_METHOD, a, nh, w, gStrat),
              sEst = unitMean(ESTN_METHOD, a, nh, w, sStrat),
              bioEst = unitMean(ESTN_METHOD, a, nh, w, bioStrat),

              N = dplyr::first(p2eu),
              #aEst = dplyr::first(aEst),
              # Estimation of unit variance
              tVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, tv, tStrat, tEst),
              dVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, dv, dStrat, dEst),
              bVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bv, bStrat, bEst),
              gVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, gv, gStrat, gEst),
              sVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, sv, sStrat, sEst),
              bioVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, biov, bioStrat, bioEst),

              cvEst_d = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_d, dStrat, dEst, tStrat, tEst),
              cvEst_b = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_b, bStrat, bEst, tStrat, tEst),
              cvEst_g = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_g, gStrat, gEst, tStrat, tEst),
              cvEst_s = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_s, sStrat, sEst, tStrat, tEst),
              cvEst_bio = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_bio, bioStrat, bioEst, tStrat, tEst),
              cvEst_baa = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_baa, bStrat, bEst, aStrat, aEst),
              cvEst_ga = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_ga, gStrat, gEst, aStrat, aEst),
              cvEst_sa = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_sa, sStrat, sEst, aStrat, aEst),
              cvEst_bioA = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvStrat_bioA, bioStrat, bioEst, aStrat, aEst),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}










