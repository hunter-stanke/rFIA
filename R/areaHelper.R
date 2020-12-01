areaHelper1 <- function(x, plts, db, grpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]

  ## was treeDomain NULL? If so, replace NAs w/ 1
  treeD <- ifelse(mean(db$TREE$tD, na.rm = TRUE) == 1, 1, 0)

  grpSyms <- syms(grpBy)

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
    lazy_dt() %>%
    mutate(tD = tidyr::replace_na(tD, treeD)) %>%
    group_by(PLT_CN, PROP_BASIS, CONDID, !!!grpSyms) %>%
    mutate(tD = ifelse(sum(tD, na.rm = TRUE) > 0, 1, 0)) %>%
    ungroup() %>%
    as.data.frame()


  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp * data$tD

  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = INVYR) %>%
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN, CONDID) %>%
      summarize(CONDPROP_UNADJ = dplyr::first(CONDPROP_UNADJ),
                aDI =dplyr::first(aDI)) %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(PERC_AREA = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                nCond = sum(aDI >  0, na.rm = TRUE)) %>%
      as.data.frame()


  } else {

    grpSyms <- syms(grpBy)

    ### Plot-level estimates
    t <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at stratum level
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(PLT_CN, PROP_BASIS, !!!grpSyms) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      as.data.frame()

  }

  pltOut <- list(t = t)
  return(pltOut)

}



areaHelper2 <- function(x, popState, t, grpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$P1POINTCNT <- popState[[x]]$P1POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  grpSyms <- syms(grpBy)

  ## Strata level estimates
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
      fa = fa * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, !!!grpSyms) %>%
    summarize(nh = dplyr::first(P2POINTCNT),
              aStrat = sum(fa, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              p2eu = dplyr::first(p2eu),
              av = sum(fa^2, na.rm = TRUE)) %>%
    mutate(aStrat = aStrat / nh, # Strata mean
           av = (av - (nh * aStrat^2)) / (nh * (nh-1))) %>% # Strata variance
    as.data.frame() %>%
    ## Estimation unit
    group_by(ESTN_UNIT_CN, !!!grpSyms) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, av, aStrat, aEst),
              N = dplyr::first(p2eu),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE)) %>%
    distinct()

  out <- list(tEst = tEst)

  return(out)
}













