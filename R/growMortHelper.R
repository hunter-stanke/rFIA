gmHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){

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
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'state_recr')), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPARECR_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
    left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA', 'state')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
    #left_join(select(db$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'SUBP', 'CONDID'), suffix = c('', '.subp')) %>%
    #left_join(select(db$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN.subp' = 'PLT_CN', 'PREVCOND.subp' = 'CONDID'), suffix = c('', '.chng')) %>%
    left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
    left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
    mutate_if(is.factor,
              as.character) %>%
    mutate(TPAGROW_UNADJ = TPAGROW_UNADJ * state,
           TPAREMV_UNADJ = TPAREMV_UNADJ * state,
           TPAMORT_UNADJ = TPAMORT_UNADJ * state,
           TPARECR_UNADJ = TPARECR_UNADJ * state_recr / REMPER,
           # mutate(SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25,
           #        TPAGROW_UNADJ = TPAGROW_UNADJ,
           #        TPAMORT_UNADJ1 = TPAMORT_UNADJ,
           #        TPAREMV_UNADJ = TPAREMV_UNADJ,
           #        TPARECR_UNADJ = TPARECR_UNADJ / REMPER,
           #aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.chng == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
           aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
           tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0),
           test = if_else(COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2'), 1, 0))

  # ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  # data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
  #   left_join(select(db$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'PREV_PLT_CN')) %>%
  #   left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN', 'CONDID')) %>%
  #   left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
  #   left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'state_recr')), by = c('PLT_CN', 'CONDID', 'SUBP')) %>%
  #   left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPARECR_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
  #   left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA', 'state')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
  #   #left_join(select(db$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.chng')) %>%
  #   left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
  #   left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
  #   mutate_if(is.factor,
  #             as.character) %>%
  #   mutate(SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25,
  #          TPAGROW_UNADJ = TPAGROW_UNADJ * state,
  #          TPAREMV_UNADJ = TPAREMV_UNADJ * state,
  #          TPAMORT_UNADJ = TPAMORT_UNADJ * state,
  #          TPARECR_UNADJ = TPARECR_UNADJ * state_recr / REMPER,
  #          # mutate(SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25,
  #          #        TPAGROW_UNADJ = TPAGROW_UNADJ,
  #          #        TPAMORT_UNADJ1 = TPAMORT_UNADJ,
  #          #        TPAREMV_UNADJ = TPAREMV_UNADJ,
  #          #        TPARECR_UNADJ = TPARECR_UNADJ / REMPER,
  #          #aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.chng == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
  #          aChng = if_else(COND_STATUS_CD == 1 &
  #                          COND_STATUS_CD.prev == 1 &
  #                          !is.null(CONDPROP_UNADJ) &
  #                          SUBPTYP == 1,
  #                        1, 0),
  #          tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0))

  # If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
  data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
  data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
  data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
  ## Issue is here for pre-growth accounting
  #data$landD.prev <- ifelse(data$landD == 1 & data$landD.prev == 1, 1, 0)
  data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
  data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
  data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

  ## Comprehensive indicator function -- w/ growth accounting
  #data$aDI_ga <- data$landD * data$aD_p * data$aD_c * data$sp * data$aChng
  data$tDI_ga <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev #* data$tChng
  data$tDI_ga_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp #* data$tChng

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev
  data$tDI_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp

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

  aData$landD <- ifelse(aData$landD == 1 & aData$landD.prev == 1, 1, 0)
  aData$aDI_ga <- aData$landD * aData$aD_p * aData$aD_c * aData$sp * aData$aChng


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(RECR_TPA = sum(TPARECR_UNADJ * tDI, na.rm = TRUE),
                MORT_TPA = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
                REMV_TPA = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
                TOTAL_TPA = sum(TPAGROW_UNADJ * tDI, na.rm = TRUE) + MORT_TPA + REMV_TPA - RECR_TPA,
                RECR_PERC = RECR_TPA / TOTAL_TPA * 100,
                MORT_PERC = MORT_TPA / TOTAL_TPA * 100,
                REMV_PERC = REMV_TPA / TOTAL_TPA * 100,
                nStems = length(which(tDI == 1)))

    a = NULL

  } else {
    ### Plot-level estimates -- growth accounting
    a_ga <- aData %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      #distinct(PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa_ga = sum(SUBPTYP_PROP_CHNG * aDI_ga, na.rm = TRUE),
                plotIn_ga = ifelse(sum(aDI_ga >  0, na.rm = TRUE), 1,0))
    ### Plot-level estimates
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      left_join(select(a_ga, PLT_CN, PROP_BASIS, aGrpBy, fa_ga, plotIn_ga), by = c('PLT_CN', 'PROP_BASIS', aGrpBy))


    ### Compute total TREES in domain of interest
    t <- data %>%
      distinct(PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(PLT_CN, SUBPTYP_GRM, .dots = grpBy) %>%
      summarize(rPlot_ga = sum(TPARECR_UNADJ * tDI_ga_r, na.rm = TRUE),
                mPlot_ga = sum(TPAMORT_UNADJ * tDI_ga, na.rm = TRUE),
                hPlot_ga = sum(TPAREMV_UNADJ * tDI_ga, na.rm = TRUE),
                tPlot_ga = sum(TPAGROW_UNADJ * tDI_ga, na.rm = TRUE) + mPlot_ga + hPlot_ga - rPlot_ga,
                plotIn_t_ga = ifelse(tPlot_ga >  0, 1,0),
                plotIn_r_ga = ifelse(rPlot_ga >  0, 1,0),
                plotIn_m_ga = ifelse(mPlot_ga > 0, 1,0),
                plotIn_h_ga = ifelse(hPlot_ga >  0, 1,0),
                ## No growth accoutning
                rPlot = sum(TPARECR_UNADJ * tDI_r, na.rm = TRUE),
                mPlot = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
                hPlot = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
                tPlot = sum(TPAGROW_UNADJ * tDI, na.rm = TRUE) + mPlot + hPlot - rPlot,
                plotIn_t = ifelse(tPlot >  0, 1,0),
                plotIn_r = ifelse(rPlot >  0, 1,0),
                plotIn_m = ifelse(mPlot > 0, 1,0),
                plotIn_h = ifelse(hPlot >  0, 1,0))
  }




  pltOut <- list(a = a, t = t)
  return(pltOut)

}



gmHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){

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
    #filter(EVAL_TYP %in% c('EXPMORT')) %>%
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
      fa = case_when(
        GROWTH_ACCT == 'Y' ~ fa_ga * aAdj,
        TRUE ~ fa * aAdj),
      plotIn = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_ga,
        TRUE ~ plotIn)) %>%
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
    filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
    #filter(EVAL_TYP %in% c('EXPGROW')) %>%
    ## Need this for covariance later on
    left_join(select(a, fa, fa_ga, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
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
      fa = case_when(
        GROWTH_ACCT == 'Y' ~ fa_ga * aAdj,
        TRUE ~ fa * aAdj),
      tPlot = case_when(
        GROWTH_ACCT == 'Y' ~ tPlot_ga * tAdj,
        TRUE ~ tPlot * tAdj),
      rPlot = case_when(
        GROWTH_ACCT == 'Y' ~ rPlot_ga * tAdj,
        TRUE ~ rPlot * tAdj),
      mPlot = case_when(
        GROWTH_ACCT == 'Y' ~ mPlot_ga * tAdj,
        TRUE ~ mPlot * tAdj),
      hPlot = case_when(
        GROWTH_ACCT == 'Y' ~ hPlot_ga * tAdj,
        TRUE ~ hPlot * tAdj),
      plotIn_t = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_t_ga,
        TRUE ~ plotIn_t),
      plotIn_r = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_r_ga,
        TRUE ~ plotIn_r),
      plotIn_m = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_m_ga,
        TRUE ~ plotIn_m),
      plotIn_h = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn_h_ga,
        TRUE ~ plotIn_h)) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(tPlot = sum(tPlot, na.rm = TRUE),
              rPlot = sum(rPlot, na.rm = TRUE),
              mPlot = sum(mPlot, na.rm = TRUE),
              hPlot = sum(hPlot, na.rm = TRUE),
              fa = first(fa),
              plotIn_t = ifelse(sum(plotIn_t >  0, na.rm = TRUE), 1,0),
              plotIn_r = ifelse(sum(plotIn_r >  0, na.rm = TRUE), 1,0),
              plotIn_m = ifelse(sum(plotIn_m >  0, na.rm = TRUE), 1,0),
              plotIn_h = ifelse(sum(plotIn_h >  0, na.rm = TRUE), 1,0),
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
              rStrat = mean(rPlot * r_t, na.rm = TRUE),
              mStrat = mean(mPlot * r_t, na.rm = TRUE),
              hStrat = mean(hPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              plotIn_r = sum(plotIn_r, na.rm = TRUE),
              plotIn_m = sum(plotIn_m, na.rm = TRUE),
              plotIn_h = sum(plotIn_h, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              tv = stratVar(ESTN_METHOD, tPlot, tStrat, ndif, a, nh),
              rv = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh),
              mv = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh),
              hv = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh),
              # Strata level covariances
              cvStrat_r = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh, fa, aStrat),
              cvStrat_m = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh, fa, aStrat),
              cvStrat_h = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh, fa, aStrat),
              cvStrat_rp = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh, tPlot, tStrat),
              cvStrat_mp = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh, tPlot, tStrat),
              cvStrat_hp = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh, tPlot, tStrat)
    ) %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
              rEst = unitMean(ESTN_METHOD, a, nh, w, rStrat),
              mEst = unitMean(ESTN_METHOD, a, nh, w, mStrat),
              hEst = unitMean(ESTN_METHOD, a, nh, w, hStrat),
              #aEst = first(aEst),
              # Estimation of unit variance
              tVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, tv, tStrat, tEst),
              rVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, rv, rStrat, rEst),
              mVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mv, mStrat, mEst),
              hVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, hv, hStrat, hEst),
              cvEst_r = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_r, rStrat, rEst, aStrat, aEst),
              cvEst_m = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_m, mStrat, mEst, aStrat, aEst),
              cvEst_h = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_h, hStrat, hEst, aStrat, aEst),
              cvEst_rp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_rp, rStrat, rEst, tStrat, tEst),
              cvEst_mp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_mp, mStrat, mEst, tStrat, tEst),
              cvEst_hp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_hp, hStrat, hEst, tStrat, tEst),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              plotIn_r = sum(plotIn_r, na.rm = TRUE),
              plotIn_m = sum(plotIn_m, na.rm = TRUE),
              plotIn_h = sum(plotIn_h, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}








growMortHelper <- function(x, combos, data, grpBy, aGrpBy, totals, SE, chngAdj){
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

    # ### Compute total TREES in domain of interest
    # tInt <- data %>%
    #   distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, COMPONENT, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
    #   #filter(EVALID %in% tID) %>%
    #   # Compute estimates at plot level
    #   group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%

    ### Compute total TREES in domain of interest
    tInt <- data %>%
      filter(EVAL_TYP %in% c('EXPGROW','EXPMORT', 'EXPREMV')) %>%
      #distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE),
                rPlot = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'] / REMPER[COMPONENT == 'INGROWTH'], na.rm = TRUE),
                mPlot = sum(TPAMORT_UNADJ* tAdj * tDI, na.rm = TRUE),
                hPlot = sum(TPAREMV_UNADJ * tAdj * tDI, na.rm = TRUE),
                #prevPop = tPlot + mPlot * first(REMPER) + hPlot * first(REMPER) - rPlot * first(REMPER),
                #lPlot = (tPlot / ifelse(prevPop < 1, NA, prevPop)) ^ (1/first(REMPER)) - 1,
                #REMPER = first(REMPER),
                plotIn_g = ifelse(tPlot >  0, 1,0),
                plotIn_r = ifelse(rPlot >  0, 1,0),
                plotIn_m = ifelse(mPlot > 0, 1,0),
                plotIn_h = ifelse(hPlot >  0, 1,0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT))

    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      #filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(SUBPTYP_PROP_CHNG * chngAdj * aDI * aAdj, na.rm = TRUE),
                plotIn_a = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT))


    ## Compute COVARIANCE between numerator and denominator (for ratio estimates of variance)
    t <- tInt %>%
      left_join(aInt, by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', 'PLT_CN'), suffix = c('_t', '_a'))  %>%
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
                nPlots_TREE = sum(plotIn_g, na.rm = TRUE),
                nPlots_RECR = sum(plotIn_r, na.rm = TRUE),
                nPlots_MORT = sum(plotIn_m, na.rm = TRUE),
                nPlots_REMV = sum(plotIn_h, na.rm = TRUE),
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
                nPlots_RECR = sum(nPlots_RECR, na.rm = TRUE),
                nPlots_MORT = sum(nPlots_MORT, na.rm = TRUE),
                nPlots_REMV = sum(nPlots_REMV, na.rm = TRUE),
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
                nPlots_RECR = sum(nPlots_RECR, na.rm = TRUE),
                nPlots_MORT = sum(nPlots_MORT, na.rm = TRUE),
                nPlots_REMV = sum(nPlots_REMV, na.rm = TRUE),
                nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE))

    # Make some columns go away
    if (totals) {
      t <- t %>%
        select(RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
               TREE_TOTAL, RECR_TOTAL, MORT_TOTAL, REMV_TOTAL, AREA_TOTAL,
               RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE,
               RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE,
               TREE_TOTAL_SE, RECR_TOTAL_SE, MORT_TOTAL_SE, REMV_TOTAL_SE, AREA_TOTAL_SE,
               nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV, nPlots_AREA)
    } else {
      t <- t %>%
        select(RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
               RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE,
               RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE,
               nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV, nPlots_AREA)
    }
    # Rejoin with some grpBy Names
    t <- data.frame(combos[[x]], t)


  } else { # No sampling errors
    ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
    ### Compute total TREES in domain of interest
    tInt <- data %>%
      filter(EVAL_TYP %in% c('EXPGROW','EXPMORT', 'EXPREMV')) %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI * EXPNS, na.rm = TRUE),
                rPlot = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'] * EXPNS[COMPONENT == 'INGROWTH'], na.rm = TRUE),
                mPlot = sum(TPAMORT_UNADJ* tAdj * tDI * EXPNS, na.rm = TRUE),
                hPlot = sum(TPAREMV_UNADJ * tAdj * tDI * EXPNS, na.rm = TRUE),
                plotIn_g = ifelse(tPlot >  0, 1,0),
                plotIn_r = ifelse(rPlot >  0, 1,0),
                plotIn_m = ifelse(mPlot > 0, 1,0),
                plotIn_h = ifelse(hPlot >  0, 1,0))

    ### Compute total AREA in the domain of interest
    aInt <- data %>%
      #filter(EVAL_TYP == 'EXPCURR') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
      group_by(.dots = aGrpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(fa = sum(SUBPTYP_PROP_CHNG * chngAdj * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn_a = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    suppressMessages({
      t <- tInt %>%
        inner_join(aInt) %>%
        ## Full region
        group_by(.dots = grpBy) %>%
        summarize(TREE_TOTAL = sum(tPlot, na.rm = TRUE),
                  RECR_TOTAL = sum(rPlot, na.rm = TRUE),
                  MORT_TOTAL = sum(mPlot, na.rm = TRUE),
                  REMV_TOTAL = sum(hPlot, na.rm = TRUE),
                  AREA_TOTAL = sum(fa, na.rm = TRUE),
                  RECR_TPA = RECR_TOTAL / AREA_TOTAL,
                  MORT_TPA = MORT_TOTAL / AREA_TOTAL,
                  REMV_TPA = REMV_TOTAL / AREA_TOTAL,
                  RECR_PERC = RECR_TOTAL / TREE_TOTAL * 100,
                  MORT_PERC = MORT_TOTAL / TREE_TOTAL * 100,
                  REMV_PERC = REMV_TOTAL / TREE_TOTAL * 100,
                  # Non-zero plots
                  nPlots_TREE = sum(plotIn_g, na.rm = TRUE),
                  nPlots_RECR = sum(plotIn_r, na.rm = TRUE),
                  nPlots_MORT = sum(plotIn_m, na.rm = TRUE),
                  nPlots_REMV = sum(plotIn_h, na.rm = TRUE),
                  nPlots_AREA = sum(plotIn_a, na.rm = TRUE))
    })

    # Make some columns go away
    if (totals) {
      t <- t %>%
        select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
               TREE_TOTAL, RECR_TOTAL, MORT_TOTAL, REMV_TOTAL, AREA_TOTAL,
               nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV,nPlots_AREA)
    } else {
      t <- t %>%
        select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
               nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV,nPlots_AREA)
    }
  } # End SE Conditional

  # Some cleanup
  #gc()

  # Return t
  t
}



#
#
# gmHelper1 <- function(x, plts, db, grpBy, aGrpBy, byPlot){
#
#   ## Selecting the plots for one county
#   db$PLOT <- plts[[x]]
#   ## Carrying out filter across all tables
#   #db <- clipFIA(db, mostRecent = FALSE)
#
#   # Only subplots from cond change matrix
#   #db$SUBP_COND_CHNG_MTRX <- filter(db$SUBP_COND_CHNG_MTRX, SUBPTYP == 1)
#
#
#   ## Which grpByNames are in which table? Helps us subset below
#   grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
#   grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
#   grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]
#
#   ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
#   data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
#     left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
#     left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'state_recr')), by = c('PLT_CN', 'CONDID')) %>%
#     left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPARECR_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
#     left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA', 'state')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
#     #left_join(select(db$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'SUBP', 'CONDID'), suffix = c('', '.subp')) %>%
#     #left_join(select(db$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN.subp' = 'PLT_CN', 'PREVCOND.subp' = 'CONDID'), suffix = c('', '.chng')) %>%
#     left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
#     left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
#     left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
#     mutate_if(is.factor,
#               as.character) %>%
#     mutate(TPAGROW_UNADJ = TPAGROW_UNADJ * state,
#            TPAREMV_UNADJ = TPAREMV_UNADJ * state,
#            TPAMORT_UNADJ = TPAMORT_UNADJ * state,
#            TPARECR_UNADJ = TPARECR_UNADJ * state_recr / REMPER,
#     # mutate(SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25,
#     #        TPAGROW_UNADJ = TPAGROW_UNADJ,
#     #        TPAMORT_UNADJ1 = TPAMORT_UNADJ,
#     #        TPAREMV_UNADJ = TPAREMV_UNADJ,
#     #        TPARECR_UNADJ = TPARECR_UNADJ / REMPER,
#            #aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.chng == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
#            aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
#            tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0),
#            test = if_else(COMPONENT %in% c('INGROWTH', 'CUT2', 'MORTALITY2'), 1, 0))
#
#   # ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
#   # data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
#   #   left_join(select(db$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'PREV_PLT_CN')) %>%
#   #   left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN', 'CONDID')) %>%
#   #   left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
#   #   left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'state_recr')), by = c('PLT_CN', 'CONDID', 'SUBP')) %>%
#   #   left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPARECR_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT')), by = c('TRE_CN')) %>%
#   #   left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA', 'state')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
#   #   #left_join(select(db$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.chng')) %>%
#   #   left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
#   #   left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
#   #   mutate_if(is.factor,
#   #             as.character) %>%
#   #   mutate(SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25,
#   #          TPAGROW_UNADJ = TPAGROW_UNADJ * state,
#   #          TPAREMV_UNADJ = TPAREMV_UNADJ * state,
#   #          TPAMORT_UNADJ = TPAMORT_UNADJ * state,
#   #          TPARECR_UNADJ = TPARECR_UNADJ * state_recr / REMPER,
#   #          # mutate(SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25,
#   #          #        TPAGROW_UNADJ = TPAGROW_UNADJ,
#   #          #        TPAMORT_UNADJ1 = TPAMORT_UNADJ,
#   #          #        TPAREMV_UNADJ = TPAREMV_UNADJ,
#   #          #        TPARECR_UNADJ = TPARECR_UNADJ / REMPER,
#   #          #aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.chng == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
#   #          aChng = if_else(COND_STATUS_CD == 1 &
#   #                          COND_STATUS_CD.prev == 1 &
#   #                          !is.null(CONDPROP_UNADJ) &
#   #                          SUBPTYP == 1,
#   #                        1, 0),
#   #          tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0))
#
#   #If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
#   data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
#   data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
#   data$landD.prev <- ifelse(data$landD == 1 & data$landD.prev == 1, 1, 0)
#   data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
#   data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
#   data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)
#
#   ## Comprehensive indicator function -- w/ growth accounting
#   #data$aDI_ga <- data$landD * data$aD_p * data$aD_c * data$sp * data$aChng
#   data$tDI_ga <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng
#   data$tDI_ga_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp #* data$tChng
#
#   ## Comprehensive indicator function
#   data$aDI <- data$landD.prev * data$aD_p * data$aD_c * data$sp
#   data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev
#   data$tDI_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
#
#   ### DOING AREA SEPARATELY NOW FOR GROWTH ACCOUNTING PLOTS
#   aData <- select(db$PLOT,c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
#     left_join(select(db$SUBP_COND_CHNG_MTRX, PLT_CN, PREV_PLT_CN, SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID), by = c('PLT_CN', 'PREV_PLT_CN')) %>%
#     left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN', 'CONDID')) %>%
#     left_join(select(db$COND, c('PLT_CN', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
#     #left_join(select(db$POP_PLOT_STRATUM_ASSGN, by = 'PLT_CN')) %>%
#     #left_join(select(db$POP_STRATUM, by = c('STRATUM_CN' = 'CN'))) %>%
#     mutate(aChng = if_else(COND_STATUS_CD == 1 &
#                              COND_STATUS_CD.prev == 1 &
#                              !is.null(CONDPROP_UNADJ) &
#                              SUBPTYP == 1,
#                            1, 0),
#            SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25)
#
#   aData$aDI_ga <- aData$landD * aData$landD.prev * aData$aD_p * aData$aD_c * aData$sp * aData$aChng
#   aData$aDI <- aData$landD * aData$aD_p * aData$aD_c * aData$sp
#
#   # ### DOING AREA SEPARATELY NOW FOR GROWTH ACCOUNTING PLOTS
#   # aData <- select(db$PLOT,c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
#   #   left_join(select(db$SUBP_COND_CHNG_MTRX, PLT_CN, PREV_PLT_CN, SUBPTYP, SUBPTYP_PROP_CHNG, PREVCOND, CONDID), by = c('PLT_CN', 'PREV_PLT_CN')) %>%
#   #   left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN', 'CONDID')) %>%
#   #   left_join(select(db$COND, c('PLT_CN', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
#   #   #left_join(select(db$POP_PLOT_STRATUM_ASSGN, by = 'PLT_CN')) %>%
#   #   #left_join(select(db$POP_STRATUM, by = c('STRATUM_CN' = 'CN'))) %>%
#   #   mutate(aChng = if_else(COND_STATUS_CD == 1 &
#   #                          COND_STATUS_CD.prev == 1 &
#   #                          !is.null(CONDPROP_UNADJ) &
#   #                          SUBPTYP == 1,
#   #                        1, 0),
#   #          SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25)
#   #
#   # aData$landD.prev <- ifelse(aData$landD == 1 & aData$landD.prev == 1, 1, 0)
#   # aData$aD_p.prev <- ifelse(is.na(aData$aD_p.prev), aData$aD_p, aData$aD_p.prev)
#   # aData$aD_c.prev <- ifelse(is.na(aData$aD_c.prev), aData$aD_c, aData$aD_c.prev)
#   # aData$sp.prev <- ifelse(is.na(aData$sp.prev), aData$sp, aData$sp.prev)
#   #
#   # aData$aDI_ga <- aData$landD.prev * aData$aD_p * aData$aD_c * aData$sp * aData$aChng
#
#
#   if (byPlot){
#     grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
#     t <- data %>%
#       mutate(YEAR = MEASYEAR) %>%
#       distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
#       # Compute estimates at plot level
#       group_by(.dots = grpBy, PLT_CN) %>%
#       summarize(RECR_TPA = sum(TPARECR_UNADJ * tDI, na.rm = TRUE),
#                 MORT_TPA = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
#                 REMV_TPA = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
#                 TOTAL_TPA = sum(TPAGROW_UNADJ * tDI, na.rm = TRUE) + MORT_TPA + REMV_TPA - RECR_TPA,
#                 RECR_PERC = RECR_TPA / TOTAL_TPA * 100,
#                 MORT_PERC = MORT_TPA / TOTAL_TPA * 100,
#                 REMV_PERC = REMV_TPA / TOTAL_TPA * 100,
#                 nStems = length(which(tDI == 1)))
#
#     a = NULL
#
#   } else {
#     ### Plot-level estimates -- growth accounting
#     a_ga <- aData %>%
#       filter(!is.na(PROP_BASIS)) %>%
#       ## Will be lots of trees here, so CONDPROP listed multiple times
#       ## Adding PROP_BASIS so we can handle adjustment factors at strata level
#       #distinct(PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
#       group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
#       summarize(fa_ga = sum(SUBPTYP_PROP_CHNG * aDI_ga, na.rm = TRUE),
#                 plotIn_ga = ifelse(sum(aDI_ga > 0, na.rm = TRUE), 1,0))
#     ### Plot-level estimates
#     a <- data %>%
#       ## Will be lots of trees here, so CONDPROP listed multiple times
#       ## Adding PROP_BASIS so we can handle adjustment factors at strata level
#       distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
#       group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
#       summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
#                 plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
#       left_join(select(a_ga, PLT_CN, PROP_BASIS, aGrpBy, fa_ga, plotIn_ga), by = c('PLT_CN', 'PROP_BASIS', aGrpBy))
#
#
#     ### Compute total TREES in domain of interest
#     t <- data %>%
#       distinct(PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
#       # Compute estimates at plot level
#       group_by(PLT_CN, SUBPTYP_GRM, .dots = grpBy) %>%
#       summarize(rPlot_ga = sum(TPARECR_UNADJ * tDI_ga_r, na.rm = TRUE),
#                 mPlot_ga = sum(TPAMORT_UNADJ * tDI_ga, na.rm = TRUE),
#                 hPlot_ga = sum(TPAREMV_UNADJ * tDI_ga, na.rm = TRUE),
#                 tPlot_ga = sum(TPAGROW_UNADJ * tDI_ga, na.rm = TRUE) + mPlot_ga + hPlot_ga - rPlot_ga,
#                 plotIn_t_ga = ifelse(tPlot_ga >  0, 1,0),
#                 plotIn_r_ga = ifelse(rPlot_ga >  0, 1,0),
#                 plotIn_m_ga = ifelse(mPlot_ga > 0, 1,0),
#                 plotIn_h_ga = ifelse(hPlot_ga >  0, 1,0),
#                 ## No growth accoutning
#                 rPlot = sum(TPARECR_UNADJ * tDI_r, na.rm = TRUE),
#                 mPlot = sum(TPAMORT_UNADJ * tDI, na.rm = TRUE),
#                 hPlot = sum(TPAREMV_UNADJ * tDI, na.rm = TRUE),
#                 tPlot = sum(TPAGROW_UNADJ * tDI, na.rm = TRUE) + mPlot + hPlot - rPlot,
#                 plotIn_t = ifelse(tPlot >  0, 1,0),
#                 plotIn_r = ifelse(rPlot >  0, 1,0),
#                 plotIn_m = ifelse(mPlot > 0, 1,0),
#                 plotIn_h = ifelse(hPlot >  0, 1,0))
#   }
#
#
#
#
#   pltOut <- list(a = a, t = t)
#   return(pltOut)
#
# }
#
#
#
# gmHelper2 <- function(x, popState, a, t, grpBy, aGrpBy, method){
#
#   ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
#   if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
#     grpBy <- c(grpBy, 'INVYR')
#     aGrpBy <- c(aGrpBy, 'INVYR')
#     popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
#     popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
#
#   }
#
#   ## Strata level estimates
#   aStrat <- a %>%
#     ## Rejoin with population tables
#     right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
#     filter(EVAL_TYP %in% c('EXPGROW')) %>%
#     #filter(EVAL_TYP %in% c('EXPCURR')) %>%
#     mutate(
#       ## AREA
#       aAdj = case_when(
#         ## When NA, stay NA
#         is.na(PROP_BASIS) ~ NA_real_,
#         ## If the proportion was measured for a macroplot,
#         ## use the macroplot value
#         PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
#         ## Otherwise, use the subpplot value
#         PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
#       fa = case_when(
#         GROWTH_ACCT == 'Y' ~ fa_ga * aAdj,
#         TRUE ~ fa * aAdj),
#       plotIn = case_when(
#         GROWTH_ACCT == 'Y' ~ plotIn_ga,
#         TRUE ~ plotIn)) %>%
#     group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = aGrpBy) %>%
#     summarize(a_t = length(unique(PLT_CN)) / first(P2POINTCNT),
#               aStrat = mean(fa * a_t, na.rm = TRUE),
#               plotIn_AREA = sum(plotIn, na.rm = TRUE),
#               n = n(),
#               ## We don't want a vector of these values, since they are repeated
#               nh = first(P2POINTCNT),
#               a = first(AREA_USED),
#               w = first(P1POINTCNT) / first(P1PNTCNT_EU),
#               p2eu = first(p2eu),
#               ndif = nh - n,
#               ## Strata level variances
#               av = stratVar(ESTN_METHOD, fa, aStrat, ndif, a, nh))
#   ## Estimation unit
#   aEst <- aStrat %>%
#     group_by(ESTN_UNIT_CN, .dots = aGrpBy) %>%
#     summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
#               aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, av, aStrat, aEst),
#               plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))
#
#   ######## ------------------ TREE ESTIMATES + CV
#
#   ## Strata level estimates
#   tEst <- t %>%
#     ## Rejoin with population tables
#     right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
#     filter(EVAL_TYP %in% c('EXPGROW', 'EXPMORT', 'EXPREMV')) %>%
#     #filter(EVAL_TYP %in% c('EXPGROW')) %>%
#     ## Need this for covariance later on
#     left_join(select(a, fa, fa_ga, PLT_CN, PROP_BASIS, aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE]), by = c('PLT_CN', aGrpBy[aGrpBy %in% c('YEAR', 'INVYR') == FALSE])) %>%
#     #Add adjustment factors
#     mutate(tAdj = case_when(
#       ## When NA, stay NA
#       is.na(SUBPTYP_GRM) ~ NA_real_,
#       ## If the proportion was measured for a macroplot,
#       ## use the macroplot value
#       SUBPTYP_GRM == 0 ~ 0,
#       SUBPTYP_GRM == 1 ~ as.numeric(ADJ_FACTOR_SUBP),
#       SUBPTYP_GRM == 2 ~ as.numeric(ADJ_FACTOR_MICR),
#       SUBPTYP_GRM == 3 ~ as.numeric(ADJ_FACTOR_MACR)),
#       ## AREA
#       aAdj = case_when(
#         ## When NA, stay NA
#         is.na(PROP_BASIS) ~ NA_real_,
#         ## If the proportion was measured for a macroplot,
#         ## use the macroplot value
#         PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
#         ## Otherwise, use the subpplot value
#         PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
#       fa = case_when(
#         GROWTH_ACCT == 'Y' ~ fa_ga * aAdj,
#         TRUE ~ fa * aAdj),
#       tPlot = case_when(
#         GROWTH_ACCT == 'Y' ~ tPlot_ga * tAdj,
#         TRUE ~ tPlot * tAdj),
#       rPlot = case_when(
#         GROWTH_ACCT == 'Y' ~ rPlot_ga * tAdj,
#         TRUE ~ rPlot * tAdj),
#       mPlot = case_when(
#         GROWTH_ACCT == 'Y' ~ mPlot_ga * tAdj,
#         TRUE ~ mPlot * tAdj),
#       hPlot = case_when(
#         GROWTH_ACCT == 'Y' ~ hPlot_ga * tAdj,
#         TRUE ~ hPlot * tAdj),
#       plotIn_t = case_when(
#         GROWTH_ACCT == 'Y' ~ plotIn_t_ga,
#         TRUE ~ plotIn_t),
#       plotIn_r = case_when(
#         GROWTH_ACCT == 'Y' ~ plotIn_r_ga,
#         TRUE ~ plotIn_r),
#       plotIn_m = case_when(
#         GROWTH_ACCT == 'Y' ~ plotIn_m_ga,
#         TRUE ~ plotIn_m),
#       plotIn_h = case_when(
#         GROWTH_ACCT == 'Y' ~ plotIn_h_ga,
#         TRUE ~ plotIn_h)) %>%
#     ## Extra step for variance issues
#     group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
#     summarize(tPlot = sum(tPlot, na.rm = TRUE),
#               rPlot = sum(rPlot, na.rm = TRUE),
#               mPlot = sum(mPlot, na.rm = TRUE),
#               hPlot = sum(hPlot, na.rm = TRUE),
#               fa = first(fa),
#               plotIn_t = ifelse(sum(plotIn_t >  0, na.rm = TRUE), 1,0),
#               plotIn_r = ifelse(sum(plotIn_r >  0, na.rm = TRUE), 1,0),
#               plotIn_m = ifelse(sum(plotIn_m >  0, na.rm = TRUE), 1,0),
#               plotIn_h = ifelse(sum(plotIn_h >  0, na.rm = TRUE), 1,0),
#               nh = first(P2POINTCNT),
#               p2eu = first(p2eu),
#               a = first(AREA_USED),
#               w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
#     ## Joining area data so we can compute ratio variances
#     left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
#     ## Strata level
#     group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
#     summarize(r_t = length(unique(PLT_CN)) / first(nh),
#               tStrat = mean(tPlot * r_t, na.rm = TRUE),
#               rStrat = mean(rPlot * r_t, na.rm = TRUE),
#               mStrat = mean(mPlot * r_t, na.rm = TRUE),
#               hStrat = mean(hPlot * r_t, na.rm = TRUE),
#               aStrat = first(aStrat),
#               plotIn_t = sum(plotIn_t, na.rm = TRUE),
#               plotIn_r = sum(plotIn_r, na.rm = TRUE),
#               plotIn_m = sum(plotIn_m, na.rm = TRUE),
#               plotIn_h = sum(plotIn_h, na.rm = TRUE),
#               n = n(),
#               ## We don't want a vector of these values, since they are repeated
#               nh = first(nh),
#               a = first(a),
#               w = first(w),
#               p2eu = first(p2eu),
#               ndif = nh - n,
#               # ## Strata level variances
#               tv = stratVar(ESTN_METHOD, tPlot, tStrat, ndif, a, nh),
#               rv = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh),
#               mv = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh),
#               hv = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh),
#               # Strata level covariances
#               cvStrat_r = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh, fa, aStrat),
#               cvStrat_m = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh, fa, aStrat),
#               cvStrat_h = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh, fa, aStrat),
#               cvStrat_rp = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh, tPlot, tStrat),
#               cvStrat_mp = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh, tPlot, tStrat),
#               cvStrat_hp = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh, tPlot, tStrat)
#     ) %>%
#
#     ## Estimation unit
#     left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
#     group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
#     summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
#               rEst = unitMean(ESTN_METHOD, a, nh, w, rStrat),
#               mEst = unitMean(ESTN_METHOD, a, nh, w, mStrat),
#               hEst = unitMean(ESTN_METHOD, a, nh, w, hStrat),
#               #aEst = first(aEst),
#               # Estimation of unit variance
#               tVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, tv, tStrat, tEst),
#               rVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, rv, rStrat, rEst),
#               mVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mv, mStrat, mEst),
#               hVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, hv, hStrat, hEst),
#               cvEst_r = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_r, rStrat, rEst, aStrat, aEst),
#               cvEst_m = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_m, mStrat, mEst, aStrat, aEst),
#               cvEst_h = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_h, hStrat, hEst, aStrat, aEst),
#               cvEst_rp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_rp, rStrat, rEst, tStrat, tEst),
#               cvEst_mp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_mp, mStrat, mEst, tStrat, tEst),
#               cvEst_hp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_hp, hStrat, hEst, tStrat, tEst),
#               plotIn_t = sum(plotIn_t, na.rm = TRUE),
#               plotIn_r = sum(plotIn_r, na.rm = TRUE),
#               plotIn_m = sum(plotIn_m, na.rm = TRUE),
#               plotIn_h = sum(plotIn_h, na.rm = TRUE))
#
#   out <- list(tEst = tEst, aEst = aEst)
#
#   return(out)
# }
#
#
#
#
#
#
#
#
# growMortHelper <- function(x, combos, data, grpBy, aGrpBy, totals, SE, chngAdj){
#   # Update domain indicator for each each column speficed in grpBy
#   td = 1 # Start both at 1, update as we iterate through
#   ad = 1
#   for (n in 1:ncol(combos[[x]])){
#     # Tree domain indicator for each column in
#     tObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
#     if (length(which(is.na(tObs))) == length(tObs)) tObs <- 1
#     td <- data$tDI * tObs * td
#     # Area domain indicator for each column in
#     if(grpBy[n] %in% aGrpBy){
#       aObs <- as.character(combos[[x]][[aGrpBy[n]]]) == as.character(data[[aGrpBy[n]]])
#       if (length(which(is.na(aObs))) == length(aObs)) aObs <- 1
#       aObs[is.na(aObs)] <- 0
#       ad <- data$aDI * aObs * ad
#     }
#   }
#
#
#   if(SE){
#     data$tDI <- td
#     data$tDI[is.na(data$tDI)] <- 0
#     data$aDI <- ad
#     data$aDI[is.na(data$aDI)] <- 0
#     ## We produce an intermediate object in this chain as it is needed to compute the ratio of means variance
#     ## Numerator and denominator are in different domains of interest, and may be grouped by different variables
#     ## see covariance estimation below
#
#     # ### Compute total TREES in domain of interest
#     # tInt <- data %>%
#     #   distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, COMPONENT, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
#     #   #filter(EVALID %in% tID) %>%
#     #   # Compute estimates at plot level
#     #   group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#
#       ### Compute total TREES in domain of interest
#     tInt <- data %>%
#       filter(EVAL_TYP %in% c('EXPGROW','EXPMORT', 'EXPREMV')) %>%
#       #distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CONDID, SUBP, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
#       distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
#       # Compute estimates at plot level
#       group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#       summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI, na.rm = TRUE),
#                 rPlot = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'] / REMPER[COMPONENT == 'INGROWTH'], na.rm = TRUE),
#                 mPlot = sum(TPAMORT_UNADJ* tAdj * tDI, na.rm = TRUE),
#                 hPlot = sum(TPAREMV_UNADJ * tAdj * tDI, na.rm = TRUE),
#                 #prevPop = tPlot + mPlot * first(REMPER) + hPlot * first(REMPER) - rPlot * first(REMPER),
#                 #lPlot = (tPlot / ifelse(prevPop < 1, NA, prevPop)) ^ (1/first(REMPER)) - 1,
#                 #REMPER = first(REMPER),
#                 plotIn_g = ifelse(tPlot >  0, 1,0),
#                 plotIn_r = ifelse(rPlot >  0, 1,0),
#                 plotIn_m = ifelse(mPlot > 0, 1,0),
#                 plotIn_h = ifelse(hPlot >  0, 1,0),
#                 a = first(AREA_USED),
#                 p1EU = first(P1PNTCNT_EU),
#                 p1 = first(P1POINTCNT),
#                 p2 = first(P2POINTCNT))
#
#     ### Compute total AREA in the domain of interest
#     aInt <- data %>%
#       #filter(EVAL_TYP == 'EXPCURR') %>%
#       distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
#       group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#       summarize(fa = sum(SUBPTYP_PROP_CHNG * chngAdj * aDI * aAdj, na.rm = TRUE),
#                 plotIn_a = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
#                 a = first(AREA_USED),
#                 p1EU = first(P1PNTCNT_EU),
#                 p1 = first(P1POINTCNT),
#                 p2 = first(P2POINTCNT))
#
#
#     ## Compute COVARIANCE between numerator and denominator (for ratio estimates of variance)
#     t <- tInt %>%
#       left_join(aInt, by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', 'PLT_CN'), suffix = c('_t', '_a'))  %>%
#       group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
#       summarize(tStrat = mean(tPlot, na.rm = TRUE),
#                 rStrat = mean(rPlot, na.rm = TRUE),
#                 mStrat = mean(mPlot, na.rm = TRUE),
#                 hStrat = mean(hPlot, na.rm = TRUE),
#                 #lStrat = mean(lPlot, na.rm = TRUE),
#                 aStrat = mean(fa, na.rm = TRUE),
#                 a = first(a_t),
#                 w = first(p1_t) / first(p1EU_a), # Stratum weight
#                 nh = first(p2_t), # Number plots in stratum
#                 nPlots_TREE = sum(plotIn_g, na.rm = TRUE),
#                 nPlots_RECR = sum(plotIn_r, na.rm = TRUE),
#                 nPlots_MORT = sum(plotIn_m, na.rm = TRUE),
#                 nPlots_REMV = sum(plotIn_h, na.rm = TRUE),
#                 nPlots_AREA = sum(plotIn_a, na.rm = TRUE),
#                 # Strata level variances
#                 tv = ifelse(first(ESTN_METHOD == 'simple'),
#                             var(tPlot * first(a) / nh),
#                             (sum(tPlot^2) - sum(nh * tStrat^2)) / (nh * (nh-1))), # Stratified and double cases
#                 rv = ifelse(first(ESTN_METHOD == 'simple'),
#                             var(rPlot * first(a) / nh),
#                             (sum(rPlot^2) - sum(nh * rStrat^2)) / (nh * (nh-1))),
#                 mv = ifelse(first(ESTN_METHOD == 'simple'),
#                             var(mPlot * first(a) / nh),
#                             (sum(mPlot^2) - sum(nh * mStrat^2)) / (nh * (nh-1))), # Stratified and double cases
#                 hv = ifelse(first(ESTN_METHOD == 'simple'),
#                             var(hPlot * first(a) / nh),
#                             (sum(hPlot^2) - sum(nh * hStrat^2)) / (nh * (nh-1))),
#                 #lv = ifelse(first(ESTN_METHOD == 'simple'),
#                 #             var(lPlot * first(a) / nh),
#                 #             (sum(lPlot^2) - sum(nh * lStrat^2)) / (nh * (nh-1))),
#                 av = ifelse(first(ESTN_METHOD == 'simple'),
#                             var(fa * first(a) / nh),
#                             (sum(fa^2) - sum(nh * aStrat^2)) / (nh * (nh-1))),
#                 # Strata level covariances
#                 cvStrat_r = ifelse(first(ESTN_METHOD == 'simple'),
#                                    cov(fa,rPlot),
#                                    (sum(fa*rPlot) - sum(nh * aStrat *rStrat)) / (nh * (nh-1))), # Stratified and double cases
#                 cvStrat_m = ifelse(first(ESTN_METHOD == 'simple'),
#                                    cov(fa,mPlot),
#                                    (sum(fa*mPlot) - sum(nh * aStrat *mStrat)) / (nh * (nh-1))), # Stratified and double cases
#                 cvStrat_h = ifelse(first(ESTN_METHOD == 'simple'),
#                                    cov(fa,hPlot),
#                                    (sum(fa*hPlot) - sum(nh * aStrat *hStrat)) / (nh * (nh-1))),
#                 #cvStrat_l = ifelse(first(ESTN_METHOD == 'simple'),
#                #                     cov(fa,lPlot),
#                #                     (sum(fa*lPlot) - sum(nh * aStrat *lStrat)) / (nh * (nh-1))), # Stratified and double cases
#                 cvStrat_rT = ifelse(first(ESTN_METHOD == 'simple'),
#                                     cov(tPlot,rPlot),
#                                     (sum(tPlot*rPlot) - sum(nh * tStrat *rStrat)) / (nh * (nh-1))), # Stratified and double cases
#                 cvStrat_mT = ifelse(first(ESTN_METHOD == 'simple'),
#                                     cov(tPlot,mPlot),
#                                     (sum(tPlot*mPlot) - sum(nh * tStrat *mStrat)) / (nh * (nh-1))), # Stratified and double cases
#                 cvStrat_hT = ifelse(first(ESTN_METHOD == 'simple'),
#                                     cov(tPlot,hPlot),
#                                     (sum(tPlot*hPlot) - sum(nh * tStrat *hStrat)) / (nh * (nh-1)))) %>% # Stratified and double casesc) %>%
#       # Estimation Unit
#       group_by(ESTN_UNIT_CN) %>%
#       summarize(tEst = unitMean(ESTN_METHOD, a, nh, w, tStrat),
#                 rEst = unitMean(ESTN_METHOD, a, nh, w, rStrat),
#                 mEst = unitMean(ESTN_METHOD, a, nh, w, mStrat),
#                 hEst = unitMean(ESTN_METHOD, a, nh, w, hStrat),
#                 #lEst = unitMean(ESTN_METHOD, a, nh, w, lStrat),
#                 aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
#                 nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
#                 nPlots_RECR = sum(nPlots_RECR, na.rm = TRUE),
#                 nPlots_MORT = sum(nPlots_MORT, na.rm = TRUE),
#                 nPlots_REMV = sum(nPlots_REMV, na.rm = TRUE),
#                 nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE),
#                 #Variance estimates
#                 tVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, tv, tStrat, tEst),
#                 rVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, rv, rStrat, rEst),
#                 mVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, mv, mStrat, mEst),
#                 hVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, hv, hStrat, hEst),
#                 #lVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, lv, lStrat, lEst),
#                 aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
#                 # Covariance estimates
#                 cvEst_r = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_r, rStrat, rEst, aStrat, aEst),
#                 cvEst_m = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_m, mStrat, mEst, aStrat, aEst),
#                 cvEst_h = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_h, hStrat, hEst, aStrat, aEst),
#                 #cvEst_l = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_l, lStrat, lEst, aStrat, aEst),
#                 cvEst_rT = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_rT, rStrat, rEst, tStrat, tEst),
#                 cvEst_mT = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_mT, mStrat, mEst, tStrat, tEst),
#                 cvEst_hT = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cvStrat_hT, hStrat, hEst, tStrat, tEst)) %>%
#       ## Full region
#       summarize(TREE_TOTAL = sum(tEst, na.rm = TRUE),
#                 RECR_TOTAL = sum(rEst, na.rm = TRUE),
#                 MORT_TOTAL = sum(mEst, na.rm = TRUE),
#                 REMV_TOTAL = sum(hEst, na.rm = TRUE),
#                 AREA_TOTAL = sum(aEst, na.rm = TRUE),
#                 RECR_TPA = RECR_TOTAL / AREA_TOTAL,
#                 MORT_TPA = MORT_TOTAL / AREA_TOTAL,
#                 REMV_TPA = REMV_TOTAL / AREA_TOTAL,
#                 #LAMBDA = sum(lEst, na.rm = TRUE), # / AREA_TOTAL,
#                 RECR_PERC = RECR_TOTAL / TREE_TOTAL * 100,
#                 MORT_PERC = MORT_TOTAL / TREE_TOTAL * 100,
#                 REMV_PERC = REMV_TOTAL / TREE_TOTAL * 100,
#                 # Variance/covariance
#                 tVar = sum(tVar, na.rm = TRUE),
#                 rVar = sum(rVar, na.rm = TRUE),
#                 mVar = sum(mVar, na.rm = TRUE),
#                 hVar = sum(hVar, na.rm = TRUE),
#                 #lVar = sum(lVar, na.rm = TRUE),
#                 aVar = sum(aVar, na.rm = TRUE),
#                 cvR = sum(cvEst_r, na.rm = TRUE),
#                 cvM = sum(cvEst_m, na.rm = TRUE),
#                 cvH = sum(cvEst_h, na.rm = TRUE),
#                 #cvL = sum(cvEst_l, na.rm = TRUE),
#                 cvRT = sum(cvEst_rT, na.rm = TRUE),
#                 cvMT = sum(cvEst_mT, na.rm = TRUE),
#                 cvHT = sum(cvEst_hT, na.rm = TRUE),
#                 raVar = (1/AREA_TOTAL^2) * (rVar + (RECR_TPA^2 * aVar) - 2 * RECR_TPA * cvR),
#                 maVar = (1/AREA_TOTAL^2) * (mVar + (MORT_TPA^2 * aVar) - 2 * MORT_TPA * cvM),
#                 haVar = (1/AREA_TOTAL^2) * (hVar + (REMV_TPA^2 * aVar) - 2 * REMV_TPA * cvH),
#                 #laVar = (1/AREA_TOTAL^2) * (lVar + (LAMBDA^2 * aVar) - 2 * LAMBDA * cvL),
#                 rtVar = (1/TREE_TOTAL^2) * (rVar + (RECR_PERC^2 * tVar) - 2 * RECR_PERC * cvRT),
#                 mtVar = (1/TREE_TOTAL^2) * (mVar + (MORT_PERC^2 * tVar) - 2 * MORT_PERC * cvMT),
#                 htVar = (1/TREE_TOTAL^2) * (hVar + (REMV_PERC^2 * tVar) - 2 * REMV_PERC * cvHT),
#                 # Sampling Errors
#                 TREE_TOTAL_SE = sqrt(tVar) / TREE_TOTAL * 100,
#                 RECR_TOTAL_SE = sqrt(rVar) / RECR_TOTAL * 100,
#                 MORT_TOTAL_SE = sqrt(mVar) / MORT_TOTAL * 100,
#                 REMV_TOTAL_SE = sqrt(hVar) / REMV_TOTAL * 100,
#                 #LAMBDA_SE = sqrt(laVar) / LAMBDA * 100,
#                 AREA_TOTAL_SE = sqrt(aVar) / AREA_TOTAL * 100,
#                 RECR_TPA_SE = sqrt(raVar) / RECR_TPA * 100,
#                 MORT_TPA_SE = sqrt(maVar) / MORT_TPA * 100,
#                 REMV_TPA_SE = sqrt(haVar) / REMV_TPA * 100,
#                 RECR_PERC_SE = sqrt(rtVar) / RECR_TPA * 100,
#                 MORT_PERC_SE = sqrt(mtVar) / MORT_TPA * 100,
#                 REMV_PERC_SE = sqrt(htVar) / REMV_TPA * 100,
#                 # Non-zero plots
#                 nPlots_TREE = sum(nPlots_TREE, na.rm = TRUE),
#                 nPlots_RECR = sum(nPlots_RECR, na.rm = TRUE),
#                 nPlots_MORT = sum(nPlots_MORT, na.rm = TRUE),
#                 nPlots_REMV = sum(nPlots_REMV, na.rm = TRUE),
#                 nPlots_AREA = sum(nPlots_AREA, na.rm = TRUE))
#
#     # Make some columns go away
#     if (totals) {
#       t <- t %>%
#         select(RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
#                TREE_TOTAL, RECR_TOTAL, MORT_TOTAL, REMV_TOTAL, AREA_TOTAL,
#                RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE,
#                RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE,
#                TREE_TOTAL_SE, RECR_TOTAL_SE, MORT_TOTAL_SE, REMV_TOTAL_SE, AREA_TOTAL_SE,
#                nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV, nPlots_AREA)
#     } else {
#       t <- t %>%
#         select(RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
#                RECR_TPA_SE, MORT_TPA_SE, REMV_TPA_SE,
#                RECR_PERC_SE, MORT_PERC_SE, REMV_PERC_SE,
#                nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV, nPlots_AREA)
#     }
#     # Rejoin with some grpBy Names
#     t <- data.frame(combos[[x]], t)
#
#
#   } else { # No sampling errors
#     ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
#     ### Compute total TREES in domain of interest
#     tInt <- data %>%
#       filter(EVAL_TYP %in% c('EXPGROW','EXPMORT', 'EXPREMV')) %>%
#       distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, TRE_CN, COMPONENT, .keep_all = TRUE) %>%
#       # Compute estimates at plot level
#       group_by(.dots = grpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#       summarize(tPlot = sum(TPAGROW_UNADJ * tAdj * tDI * EXPNS, na.rm = TRUE),
#                 rPlot = sum(TPAGROW_UNADJ[COMPONENT == 'INGROWTH'] * tAdj[COMPONENT == 'INGROWTH'] * tDI[COMPONENT == 'INGROWTH'] * EXPNS[COMPONENT == 'INGROWTH'], na.rm = TRUE),
#                 mPlot = sum(TPAMORT_UNADJ* tAdj * tDI * EXPNS, na.rm = TRUE),
#                 hPlot = sum(TPAREMV_UNADJ * tAdj * tDI * EXPNS, na.rm = TRUE),
#                 plotIn_g = ifelse(tPlot >  0, 1,0),
#                 plotIn_r = ifelse(rPlot >  0, 1,0),
#                 plotIn_m = ifelse(mPlot > 0, 1,0),
#                 plotIn_h = ifelse(hPlot >  0, 1,0))
#
#     ### Compute total AREA in the domain of interest
#     aInt <- data %>%
#       #filter(EVAL_TYP == 'EXPCURR') %>%
#       distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
#       group_by(.dots = aGrpBy, ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
#       summarize(fa = sum(SUBPTYP_PROP_CHNG * chngAdj * aDI * aAdj * EXPNS, na.rm = TRUE),
#                 plotIn_a = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))
#
#     suppressMessages({
#       t <- tInt %>%
#         inner_join(aInt) %>%
#         ## Full region
#         group_by(.dots = grpBy) %>%
#         summarize(TREE_TOTAL = sum(tPlot, na.rm = TRUE),
#                   RECR_TOTAL = sum(rPlot, na.rm = TRUE),
#                   MORT_TOTAL = sum(mPlot, na.rm = TRUE),
#                   REMV_TOTAL = sum(hPlot, na.rm = TRUE),
#                   AREA_TOTAL = sum(fa, na.rm = TRUE),
#                   RECR_TPA = RECR_TOTAL / AREA_TOTAL,
#                   MORT_TPA = MORT_TOTAL / AREA_TOTAL,
#                   REMV_TPA = REMV_TOTAL / AREA_TOTAL,
#                   RECR_PERC = RECR_TOTAL / TREE_TOTAL * 100,
#                   MORT_PERC = MORT_TOTAL / TREE_TOTAL * 100,
#                   REMV_PERC = REMV_TOTAL / TREE_TOTAL * 100,
#                   # Non-zero plots
#                   nPlots_TREE = sum(plotIn_g, na.rm = TRUE),
#                   nPlots_RECR = sum(plotIn_r, na.rm = TRUE),
#                   nPlots_MORT = sum(plotIn_m, na.rm = TRUE),
#                   nPlots_REMV = sum(plotIn_h, na.rm = TRUE),
#                   nPlots_AREA = sum(plotIn_a, na.rm = TRUE))
#     })
#
#     # Make some columns go away
#     if (totals) {
#       t <- t %>%
#         select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
#                TREE_TOTAL, RECR_TOTAL, MORT_TOTAL, REMV_TOTAL, AREA_TOTAL,
#                nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV,nPlots_AREA)
#     } else {
#       t <- t %>%
#         select(grpBy, RECR_TPA, MORT_TPA, REMV_TPA, RECR_PERC, MORT_PERC, REMV_PERC,
#                nPlots_TREE, nPlots_RECR, nPlots_MORT, nPlots_REMV,nPlots_AREA)
#     }
#   } # End SE Conditional
#
#   # Some cleanup
#   #gc()
#
#   # Return t
#   t
# }
