acHelper1 <- function(x, plts, db, grpBy, byPlot, keepThese, chngType,
                      areaDomain, treeDomain) {

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]

  grpSyms <- syms(grpBy)


  ## Narrow up the tables to the necessary variables
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy &
                           !c(names(db$COND) %in% grpP)]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy &
                           !c(names(db$TREE) %in% c(grpP, grpC))]


  ## First, make a domain indicator for tree attributes, if any are 1s within the condition,
  ## then the condition level domain indicator is 1

  ## was treeDomain NULL? If so, replace NAs w/ 1 below
  treeD <- ifelse(mean(db$TREE$tD, na.rm = TRUE) == 1, 1, 0)

  tree <- db$TREE %>%
    lazy_dt() %>%
    group_by(PLT_CN, CONDID) %>%
    summarize(tD = sum(tD)) %>%
    mutate(tD = case_when(tD > 0 ~ 1,
                          TRUE ~ 0)) %>%
    as.data.frame()


  ## Current and previous groups
  grp1 <- if (length(grpBy) > 0 ) { paste0(grpBy, '1') } else { character(0) }
  grp2 <- if (length(grpBy) > 0 ) { paste0(grpBy, '2') } else { character(0) }
  grp1Syms <- syms(grp1)
  grp2Syms <- syms(grp2)


  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    filter(PLT_CN %in% keepThese$PLT_CN) %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(tree, by = c('PLT_CN', 'CONDID')) %>%
    filter(!is.na(PROP_BASIS) & !is.na(CONDPROP_UNADJ)) %>%
    left_join(db$SUBP_COND_CHNG_MTRX, by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$PLOT, PLT_CN, sp, all_of(grpP), aD_p), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('2', '1')) %>%
    left_join(select(db$COND, PLT_CN, landD, all_of(grpC), aD_c, CONDPROP_UNADJ, CONDID), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('2', '1')) %>%
    left_join(tree, by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('2', '1')) %>%
    rename(CONDID1 = PREVCOND,
           CONDID2 = CONDID) %>%
    ## Don't want to drop non-treed forestland
    mutate(tD1 = replace_na(tD1, treeD),
           tD2 = replace_na(tD2, treeD)) %>%
    ## Drop all microplot proportions
    filter(SUBPTYP != 2) %>%
    ## Check if macroplot is available
    group_by(PLT_CN) %>%
    mutate(macro = ifelse(3 %in% unique(SUBPTYP), 1, 0)) %>%
    ungroup() %>%
    ## If so, drop the subplot entries
    mutate(dropThese = case_when(macro == 1 & SUBPTYP == 1 ~ 0,
                                 TRUE ~ 1)) %>%
    filter(dropThese == 1) %>%
    select(-c(dropThese))





  # Use growth accounting for net change
  if (chngType == 'net') {

    data <- data %>%
      # Clean up domain indicator and drop unnecessary cols
      mutate(tDI1 = landD1 * aD_c1 * sp1 * aD_p1 * tD1,
             tDI2 = landD2 * aD_c2 * sp2 * aD_p2 * tD2) %>%
      select(PLT_CN, REMPER, SUBP, MEASYEAR, PLOT_STATUS_CD, PROP_BASIS,
             tDI1, tDI2, all_of(grp1), all_of(grp2),
             CONDID1, CONDID2,
             SUBPTYP_PROP_CHNG) %>%
      pivot_longer(cols = tDI1:CONDID2,
                   names_to = c(".value", 'ONEORTWO'),
                   names_sep = -1) %>%
      ## Adjust for subplot, negate previous values, and apply domain indicator
      mutate(PREV_CONDPROP = case_when(ONEORTWO == 1 ~ SUBPTYP_PROP_CHNG * tDI / 4,
                                       TRUE ~ 0),
             CONDPROP_CHNG = case_when(ONEORTWO == 1 ~ -SUBPTYP_PROP_CHNG * tDI / 4 / REMPER,
                                       TRUE ~ SUBPTYP_PROP_CHNG * tDI / 4 / REMPER))


  ## Change by component
  } else {

    data <- data %>%
      # Clean up domain indicator and drop unnecessary cols
      mutate(TREE_DOMAIN1 = tD1,
             TREE_DOMAIN2 = tD2,
             AREA_DOMAIN1 = landD1 * aD_c1 * sp1 * aD_p1,
             AREA_DOMAIN2 = landD2 * aD_c2 * sp2 * aD_p2) %>%
      select(PLT_CN, REMPER, SUBP, MEASYEAR, PLOT_STATUS_CD, PROP_BASIS,
             TREE_DOMAIN1, TREE_DOMAIN2, AREA_DOMAIN1, AREA_DOMAIN2,
             all_of(grp1), all_of(grp2),
             CONDID1, CONDID2,
             SUBPTYP_PROP_CHNG) %>%
      ## Did the groups change? If so, track them and drop the rest
      mutate(chng = case_when(paste(TREE_DOMAIN1, AREA_DOMAIN1, !!!grp1Syms) != paste(TREE_DOMAIN2, AREA_DOMAIN2, !!!grp2Syms) ~ 1,
                              TRUE ~ 0)) %>%
      ## Adjust for subplot, negate previous values, and apply domain indicator
      mutate(tDI1 = TREE_DOMAIN1 * AREA_DOMAIN1,
             tDI2 = TREE_DOMAIN2 * AREA_DOMAIN2,
             tDI = ifelse(tDI1 + tDI2 > 0, 1, 0), # for nPlots
             PREV_CONDPROP = SUBPTYP_PROP_CHNG * tDI / 4,
             CONDPROP_CHNG = case_when(chng == 1 ~ SUBPTYP_PROP_CHNG / 4 / REMPER,
                                       TRUE ~ 0))

    ## Update grpBy
    grpBy = sort(c(grp1, grp2))

    if (quo_name(treeDomain) != 'NULL') grpBy <- c('TREE_DOMAIN1', 'TREE_DOMAIN2', grpBy)
    if (quo_name(areaDomain) != 'NULL') grpBy <- c('AREA_DOMAIN1', 'AREA_DOMAIN2', grpBy)


  }

  ## Summarize
  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN, REMPER) %>%
      summarize(PROP_CHNG = sum(CONDPROP_CHNG, na.rm = TRUE) * 100) %>%
      as.data.frame()


  } else {

    grpSyms <- syms(grpBy)

    ## Plot-level estimates
    t <- data %>%
      lazy_dt() %>%
      group_by(PLT_CN, PROP_BASIS, !!!grpSyms) %>%
      summarize(ac = sum(CONDPROP_CHNG, na.rm = TRUE),
                prev = sum(PREV_CONDPROP, na.rm = TRUE),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0)) %>%
      as.data.frame()
  }




  pltOut <- list(t = t, grpBy = grpBy)
  return(pltOut)


}



acHelper2 <- function(x, popState, t, grpBy, method) {

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$P1POINTCNT <- popState[[x]]$P1POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  grpSyms <- syms(grpBy)


  ## Totals
  tEst <- t %>%
    lazy_dt() %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
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
      ac = ac * aAdj,
      prev = prev * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              pStrat = sum(prev, na.rm = TRUE),
              cStrat = sum(ac, na.rm = TRUE),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              p2eu = dplyr::first(p2eu),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
              pv = sum(prev^2, na.rm = TRUE),
              cv = sum(ac^2, na.rm = TRUE),
              cov_cv = sum(prev * ac, na.rm = TRUE)) %>%
    mutate(pStrat = pStrat / nh, # Strata mean
           cStrat = cStrat / nh,
           pv = (pv - (nh * pStrat^2)) / (nh * (nh-1)),
           cv = (cv - (nh * cStrat^2)) / (nh * (nh-1)),
           cov_cv = (cov_cv - (nh * cStrat * pStrat)) /  (nh * (nh-1))) %>% # Strata variance
    as.data.frame() %>%
    ## Estimation unit
    group_by(ESTN_UNIT_CN, !!!grpSyms) %>%
    summarize(pEst = unitMean(ESTN_METHOD, a, nh, w, pStrat),
              cEst = unitMean(ESTN_METHOD, a, nh, w, cStrat),
              pVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, pv, pStrat, pEst),
              cVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cv, cStrat, cEst),
              cCV = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cov_cv, cStrat, cEst, pStrat, pEst),
              N = dplyr::first(p2eu),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE)) %>%
    distinct()


  out <- list(tEst = tEst)

}
