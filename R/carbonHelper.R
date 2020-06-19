carbonHelper1 <- function(x, plts, db, grpBy, byPlot, byPool, byComponent, modelSnag){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
    ## Need a code that tells us where the tree was measured
    ## macroplot, microplot, subplot
    mutate(PLOT_BASIS = case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA) ~ NA_character_,
      ## When DIA is less than 5", use microplot value
      DIA < 5 ~ 'MICR',
      ## When DIA is greater than 5", use subplot value
      DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
      DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) %>%
    mutate(live = case_when(STATUSCD == 1 ~ 1,
                            is.na(DIA) ~ NA_real_,
                            TRUE ~ 0),
           dead = case_when(STATUSCD == 2 ~ 1,
                            is.na(DIA) ~ NA_real_,
                            TRUE ~ 0))
  #filter(!is.na(PLOT_BASIS))
  # rename(YEAR = INVYR) %>%
  # mutate_if(is.factor,
  #           as.character) %>%
  # filter(!is.na(YEAR))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp


  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')

    ### Plot-level estimates
    a <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, .dots = grpBy) %>%
      summarize(AG_UNDER_LIVE = sum(CONDPROP_UNADJ *CARBON_UNDERSTORY_AG * aDI, na.rm = TRUE),
                BG_UNDER_LIVE = sum(CONDPROP_UNADJ *CARBON_UNDERSTORY_BG * aDI, na.rm = TRUE),
                DOWN_DEAD = sum(CONDPROP_UNADJ *CARBON_DOWN_DEAD * aDI, na.rm = TRUE),
                STAND_DEAD_MOD = sum(CONDPROP_UNADJ *CARBON_STANDING_DEAD * aDI, na.rm = TRUE),
                LITTER = sum(CONDPROP_UNADJ *CARBON_LITTER * aDI, na.rm = TRUE),
                SOIL_ORG = sum(CONDPROP_UNADJ *CARBON_SOIL_ORG * aDI, na.rm = TRUE))

    ## Tree plts
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      group_by(PLT_CN, .dots = grpBy) %>%
      summarize(AG_OVER_LIVE = sum(CARBON_AG * TPA_UNADJ * live * aDI, na.rm = TRUE) / 2000,
                BG_OVER_LIVE = sum(CARBON_BG * TPA_UNADJ * live * aDI, na.rm = TRUE) / 2000,
                AG_OVER_DEAD = sum(CARBON_AG * TPA_UNADJ * dead * aDI, na.rm = TRUE) / 2000,
                BG_OVER_DEAD = sum(CARBON_BG * TPA_UNADJ * dead * aDI, na.rm = TRUE) / 2000)

    ## Join these back together
    t <- t %>%
      left_join(a, by = c('PLT_CN', grpBy))

    ## Decide which estimate to use for snags
    if (modelSnag){
      ## Model based
      t <- t %>%
        mutate(STAND_DEAD = STAND_DEAD_MOD) %>%
        select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))

    } else {
      ## Design based
      t <- t %>%
        mutate(STAND_DEAD = AG_OVER_DEAD + BG_OVER_DEAD) %>%
        select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))
    }


    ## Can't think of a better way to go about this so..
    if (byPool & byComponent == FALSE){
      ## Pool only
      t <- t %>%
        mutate(AG_LIVE = AG_UNDER_LIVE + AG_OVER_LIVE,
               BG_LIVE = BG_OVER_LIVE + BG_UNDER_LIVE,
               DEAD_WOOD = STAND_DEAD + DOWN_DEAD,
               LITTER = LITTER,
               SOIL_ORG = SOIL_ORG) %>%
        select(all_of(grpBy), PLT_CN, contains('PLOT_BASIS'), contains('PROP_BASIS'),
               contains('fa'), AG_LIVE, BG_LIVE, DEAD_WOOD, LITTER, SOIL_ORG,
               contains('plotIn_TREE'), contains('plotIn_AREA')) %>%
        pivot_longer(cols = AG_LIVE:SOIL_ORG, names_to = 'POOL', values_to = 'CARB_ACRE') %>%
        mutate(CARB_ACRE = CARB_ACRE * 0.90718474)

    } else if (byComponent){
      ## Get both pool and component here
      t <- t %>%
        select(all_of(grpBy), PLT_CN, contains('PLOT_BASIS'), contains('PROP_BASIS'),
               contains('fa'),
               AG_OVER_LIVE, AG_UNDER_LIVE, BG_OVER_LIVE, BG_UNDER_LIVE,
               DOWN_DEAD, STAND_DEAD, LITTER, SOIL_ORG, contains('plotIn_TREE'), contains('plotIn_AREA')) %>%
        pivot_longer(cols = AG_OVER_LIVE:SOIL_ORG,
                     names_to = 'COMPONENT', values_to = 'CARB_ACRE') %>%
        mutate(CARB_ACRE = CARB_ACRE * 0.90718474)

    } else{
      ## Totals
      t <- t %>%
        mutate(CARB_ACRE = AG_OVER_LIVE + AG_UNDER_LIVE +
                 BG_OVER_LIVE + BG_UNDER_LIVE +
                 DOWN_DEAD + STAND_DEAD + LITTER + SOIL_ORG) %>%
        mutate(CARB_ACRE = CARB_ACRE * 0.90718474) %>%
        select(grpBy, PLT_CN, contains('PLOT_BASIS'), contains('PROP_BASIS'),
               contains('fa'), CARB_ACRE, contains('plotIn_TREE'), contains('plotIn_AREA'))
    }

    a = NULL


  } else {
    ### Plot-level estimates
    a <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = grpBy) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                AG_UNDER_LIVE = sum(CONDPROP_UNADJ *CARBON_UNDERSTORY_AG * aDI, na.rm = TRUE),
                BG_UNDER_LIVE = sum(CONDPROP_UNADJ *CARBON_UNDERSTORY_BG * aDI, na.rm = TRUE),
                DOWN_DEAD = sum(CONDPROP_UNADJ *CARBON_DOWN_DEAD * aDI, na.rm = TRUE),
                STAND_DEAD_MOD = sum(CONDPROP_UNADJ *CARBON_STANDING_DEAD * aDI, na.rm = TRUE),
                LITTER = sum(CONDPROP_UNADJ *CARBON_LITTER * aDI, na.rm = TRUE),
                SOIL_ORG = sum(CONDPROP_UNADJ *CARBON_SOIL_ORG * aDI, na.rm = TRUE),
                plotIn_AREA = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    ## Tree plts
    t <- data %>%
      filter(!is.na(PLOT_BASIS)) %>%
      group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
      summarize(AG_OVER_LIVE = sum(CARBON_AG * TPA_UNADJ * live * aDI, na.rm = TRUE) / 2000,
                BG_OVER_LIVE = sum(CARBON_BG * TPA_UNADJ * live * aDI, na.rm = TRUE) / 2000,
                AG_OVER_DEAD = sum(CARBON_AG * TPA_UNADJ * dead * aDI, na.rm = TRUE) / 2000,
                BG_OVER_DEAD = sum(CARBON_BG * TPA_UNADJ * dead * aDI, na.rm = TRUE) / 2000,
                plotIn_TREE = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))

    # ## Join these back together
    # t <- t %>%
    #   left_join(a, by = c('PLT_CN', grpBy))

  }


  pltOut <- list(t = t, a = a)
  return(pltOut)

}



carbonHelper2 <- function(x, popState, t, a, grpBy, method, byPool, byComponent, modelSnag){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR

  }

  ######## ------------------ ADJUSTMENT FACTORS
  aAdjusted <- a %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    mutate(aAdj = case_when(
      ## When NA, stay NA
      is.na(PROP_BASIS) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
      ## Otherwise, use the subpplot value
      PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP)) %>%
    mutate(fa = fa * aAdj,
           AG_UNDER_LIVE = AG_UNDER_LIVE * aAdj,
           BG_UNDER_LIVE = BG_UNDER_LIVE * aAdj,
           DOWN_DEAD = DOWN_DEAD * aAdj,
           STAND_DEAD_MOD = STAND_DEAD_MOD * aAdj,
           LITTER = LITTER * aAdj,
           SOIL_ORG = SOIL_ORG * aAdj) %>%
    select(YEAR, PLT_CN:plotIn_AREA)

  t <- t %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    #Add adjustment factors
    mutate(tAdj = case_when(
      ## When NA, stay NA
      is.na(PLOT_BASIS) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      PLOT_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
      ## Otherwise, use the subpplot value
      PLOT_BASIS == 'SUBP' ~ as.numeric(ADJ_FACTOR_SUBP),
      PLOT_BASIS == 'MICR' ~ as.numeric(ADJ_FACTOR_MICR))) %>%
    mutate(AG_OVER_LIVE = AG_OVER_LIVE * tAdj,
           BG_OVER_LIVE = BG_OVER_LIVE * tAdj,
           AG_OVER_DEAD = AG_OVER_DEAD * tAdj,
           BG_OVER_DEAD = BG_OVER_DEAD * tAdj) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy[!(grpBy %in% c('POOL', 'COMPONENT'))]) %>%
    summarize(AG_OVER_LIVE = sum(AG_OVER_LIVE, na.rm = TRUE),
              BG_OVER_LIVE = sum(BG_OVER_LIVE, na.rm = TRUE),
              AG_OVER_DEAD = sum(AG_OVER_DEAD, na.rm = TRUE),
              BG_OVER_DEAD = sum(BG_OVER_DEAD, na.rm = TRUE),
              plotIn_TREE = ifelse(sum(plotIn_TREE >  0, na.rm = TRUE), 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    left_join(aAdjusted, by = c('PLT_CN', grpBy[!(grpBy %in% c('POOL', 'COMPONENT', 'INVYR'))]))

  ######## ------------------ SWAP TO LONG FORMAT
  ## Decide which estimate to use for snags
  if (modelSnag){
    ## Model based
    t <- t %>%
      mutate(STAND_DEAD = STAND_DEAD_MOD) %>%
      select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))

  } else {
    ## Design based
    t <- t %>%
      mutate(STAND_DEAD = AG_OVER_DEAD + BG_OVER_DEAD) %>%
      select(-c(AG_OVER_DEAD, BG_OVER_DEAD, STAND_DEAD_MOD))
  }


  ## Can't think of a better way to go about this so..
  if (byPool & byComponent == FALSE){
    ## Pool only
    t <- t %>%
      mutate(AG_LIVE = AG_UNDER_LIVE + AG_OVER_LIVE,
             BG_LIVE = BG_OVER_LIVE + BG_UNDER_LIVE,
             DEAD_WOOD = STAND_DEAD + DOWN_DEAD,
             LITTER = LITTER,
             SOIL_ORG = SOIL_ORG) %>%
      select(all_of(grpBy[!(grpBy %in% c('POOL', 'COMPONENT'))]),
             ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN,
             fa, AG_LIVE, BG_LIVE, DEAD_WOOD, LITTER, SOIL_ORG,
             plotIn_TREE, plotIn_AREA,
             nh, p2eu, a, w) %>%
      pivot_longer(cols = AG_LIVE:SOIL_ORG, names_to = 'POOL', values_to = 'CARB_ACRE') %>%
      mutate(CARB_ACRE = CARB_ACRE * 0.90718474)

  } else if (byComponent){
    ## Get both pool and component here
    t <- t %>%
      select(all_of(grpBy[!(grpBy %in% c('POOL', 'COMPONENT'))]),
             ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN,
             fa, AG_OVER_LIVE, AG_UNDER_LIVE, BG_OVER_LIVE, BG_UNDER_LIVE,
             DOWN_DEAD, STAND_DEAD, LITTER, SOIL_ORG,
             plotIn_TREE, plotIn_AREA,
             nh, p2eu, a, w) %>%
      pivot_longer(cols = AG_OVER_LIVE:SOIL_ORG,
                   names_to = 'COMPONENT', values_to = 'CARB_ACRE') %>%
      mutate(CARB_ACRE = CARB_ACRE * 0.90718474)

  } else{
    ## Totals
    t <- t %>%
      mutate(CARB_ACRE = AG_OVER_LIVE + AG_UNDER_LIVE +
               BG_OVER_LIVE + BG_UNDER_LIVE +
               DOWN_DEAD + STAND_DEAD + LITTER + SOIL_ORG) %>%
      mutate(CARB_ACRE = CARB_ACRE * 0.90718474) %>%
      select(all_of(grpBy[!(grpBy %in% c('POOL', 'COMPONENT'))]),
             ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN,
             fa, CARB_ACRE, plotIn_TREE, plotIn_AREA,
             nh, p2eu, a, w)
  }

  ######## ------------------ TREE ESTIMATES + CV

  ## Strata level estimates
  tEst <- t %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              aStrat = mean(fa * r_t, na.rm = TRUE),
              cStrat = mean(CARB_ACRE * r_t, na.rm = TRUE),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              ## Strata level variances
              av = stratVar(ESTN_METHOD, fa, aStrat, ndif, a, nh),
              cv = stratVar(ESTN_METHOD, CARB_ACRE, cStrat, ndif, a, nh),
              # Strata level covariances
              cvStrat_c = stratVar(ESTN_METHOD, CARB_ACRE, cStrat, ndif, a, nh, fa, aStrat)) %>%

    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(cEst = unitMean(ESTN_METHOD, a, nh, w, cStrat),
              aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              # Estimation of unit variance
              cVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cv, cStrat, cEst),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, av, aStrat, aEst),
              cvEst_c = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_c, cStrat, cEst, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE),
              plotIn_TREE = sum(plotIn_TREE, na.rm = TRUE))

  out <- list(tEst = tEst)

  return(out)
}



















