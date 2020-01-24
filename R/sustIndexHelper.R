sustIndexHelper1 <- function(x, plts, db, grpBy, byPlot, minLive){

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
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp', 'DESIGNCD')) %>%
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & DESIGNCD == 1) %>%

    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    ## AGENTCD at remeasurement, died during the measurement interval
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'TPA_UNADJ', 'DIA', 'AGENTCD', 'MORTYR')), by = c('PLT_CN', 'CONDID')) %>%
    #left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', DIA_BEGIN, DIA_END)), by = c('TRE_CN')) %>%

    left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p', 'DESIGNCD', 'PLOT_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('2', '1')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('2', '1')) %>%
    left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD', 'TPA_UNADJ', 'DIA')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%
    #left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%

    mutate_if(is.factor,
              as.character)

  ## Comprehensive indicator function -- w/ growth accounting
  data$tDI2 <- data$landD2 * data$aD_p2 * data$aD_c2 * data$tD2 * data$typeD2 * data$sp2 #*
    #if_else(data$PLOT_STATUS_CD1 == 1 & data$PLOT_STATUS_CD2 == 1, 1, 0)
  data$tDI1 <- data$landD1 * data$aD_p1 * data$aD_c1 * data$tD1 * data$typeD1 * data$sp1 #*
    #if_else(data$PLOT_STATUS_CD1 == 1 & data$PLOT_STATUS_CD2 == 1, 1, 0)

  ## PREVIOUS and CURRENT attributes
  data <- data %>%
    mutate(TPA_UNADJ1 = TPA_UNADJ1,
           TPA_UNADJ2 = TPA_UNADJ2,
           BAA1 = basalArea(DIA1) * TPA_UNADJ1,
           BAA2 = basalArea(DIA2) * TPA_UNADJ2,
           BUG = if_else(AGENTCD == 10, 1, 0),
           DISEASE = if_else(AGENTCD == 20, 1, 0),
           FIRE = if_else(AGENTCD == 30, 1, 0),
           ANIMAL = if_else(AGENTCD == 40, 1, 0),
           WEATHER = if_else(AGENTCD == 50, 1, 0),
           VEG = if_else(AGENTCD == 60, 1, 0),
           UNKNOWN = if_else(AGENTCD == 70, 1, 0),
           SILV = if_else(AGENTCD == 80, 1, 0),
           ) #%>%

  ## Just what we need
  data <- data %>%
    select(PLT_CN, TRE_CN, SUBP, CONDID, TREE,
          MEASYEAR, MACRO_BREAKPOINT_DIA,
          BUG, DISEASE, FIRE, ANIMAL, WEATHER, VEG, UNKNOWN, SILV, MORTYR,
          REMPER, PLOT_STATUS_CD1, PLOT_STATUS_CD2,
           one_of(str_c(grpP,1), str_c(grpC,1), str_c(grpT,1),
           str_c(grpP,2), str_c(grpC,2), str_c(grpT,2)),
          tDI1, tDI2, #SUBPTYP_GRM1, SUBPTYP_GRM2,
           DIA1, DIA2, BAA1, BAA2, TPA_UNADJ1, TPA_UNADJ2) %>%
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
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR'))


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    # tAll <- data %>%
    #   mutate(YEAR = MEASYEAR) %>%
    #   distinct(PLT_CN, PREV_PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
    #   # Compute estimates at plot level
    #   group_by(.dots = grpBy, PLT_CN, PREV_PLT_CN, REMPER) %>%
    #   summarize(TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE),
    #             BAA = sum(DIA * tDI, na.rm = TRUE),
    #             nStems = length(which(tDI == 1)))
    #
    # t <- tAll %>%
    #   filter(!is.na(PREV_PLT_CN), !is.na(REMPER)) #%>%
    # t <- t %>%
    #   ## FULL JOIN MAKES SURE MISSING VALUES BECOME ZERO
    #   full_join(filter(tAll, PLT_CN %in% t$PLT_CN == FALSE),
    #             by = c('PREV_PLT_CN' = 'PLT_CN', grpBy[grpBy %in% c('YEAR', 'pltID', 'PLOT_STATUS_CD') == FALSE]),
    #             suffix = c('_CURR', '_PREV')) #%>%
    #   mutate(YEAR = YEAR_CURR) %>%
    #   group_by(.dots = grpBy, PLT_CN) %>%
    #   summarize(CURR_TPA = if_else(nStems_PREV > minLive, TPA_CURR, 0),
    #             CURR_BAA = if_else(nStems_PREV > minLive, BAA_CURR, 0),
    #             PREV_TPA = if_else(nStems_PREV > minLive, TPA_PREV, 0),
    #             PREV_BAA = if_else(nStems_PREV > minLive, BAA_PREV, 0),
    #             CHNG_TPA = (CURR_TPA - PREV_TPA)  / REMPER_CURR,
    #             CHNG_BAA = (CURR_BAA - PREV_BAA)  / REMPER_CURR,
    #             TPA_RATE = CHNG_TPA / PREV_TPA,
    #             BAA_RATE = CHNG_BAA / PREV_BAA,
    #             x = rFIA:::projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
    #             y = rFIA:::projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
    #             M = sqrt(x^2 + y^2),
    #             SUST_INDEX = if_else(x < 0, -M, M)) #%>%
    #   select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE, PREV_TPA, CURR_TPA, PREV_BAA, CURR_BAA, nStems_PREV, nStems_CURR)
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
                REMPER = first(REMPER),
                PREV_TPA = if_else(nLive >= minLive, sum(-TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
                PREV_BAA = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
                CHNG_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ * tDI, na.rm = TRUE), 0) / REMPER,
                CHNG_BAA = if_else(nLive >= minLive, sum(BAA * tDI, na.rm = TRUE), 0) / REMPER,
                TPA_RATE = CHNG_TPA / PREV_TPA,
                BAA_RATE = CHNG_BAA / PREV_BAA,
                ## Disturbances
                BUG_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & BUG == 1] * tDI[ONEORTWO == 1 & BUG == 1], na.rm = TRUE), 0) / PREV_BAA,
                DISEASE_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & DISEASE == 1] * tDI[ONEORTWO == 1 & DISEASE == 1], na.rm = TRUE), 0) / PREV_BAA,
                FIRE_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & FIRE == 1] * tDI[ONEORTWO == 1 & FIRE == 1], na.rm = TRUE), 0) / PREV_BAA,
                ANIMAL_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & ANIMAL == 1] * tDI[ONEORTWO == 1 & ANIMAL == 1], na.rm = TRUE), 0) / PREV_BAA,
                WEATHER_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & WEATHER == 1] * tDI[ONEORTWO == 1 & WEATHER == 1], na.rm = TRUE), 0) / PREV_BAA,
                VEG_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & VEG == 1] * tDI[ONEORTWO == 1 & VEG == 1], na.rm = TRUE), 0) / PREV_BAA,
                UNKNOWN_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & UNKNOWN == 1] * tDI[ONEORTWO == 1 & UNKNOWN == 1], na.rm = TRUE), 0) / PREV_BAA,
                SILV_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & SILV == 1] * tDI[ONEORTWO == 1 & SILV == 1], na.rm = TRUE), 0) / PREV_BAA,
                x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
                y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
                M = sqrt(x^2 + y^2),
                SUST_INDEX = if_else(x < 0, -M, M),
                nStems = length(which(tDI == 1))) %>%
      ungroup() %>%
      select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE, PREV_TPA, PREV_BAA, BUG_RATE,
             DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE,
             UNKNOWN_RATE, SILV_RATE, nStems, nLive)


      #           PREV_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ1 * tDI1, na.rm = TRUE), 0),
      #           CURR_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ2 * tDI2, na.rm = TRUE), 0),
      #           CHNG_TPA = (CURR_TPA - PREV_TPA) / REMPER,
      #           PREV_BAA = if_else(nLive >= minLive, sum(BAA1 * tDI1, na.rm = TRUE), 0),
      #           CURR_BAA = if_else(nLive >= minLive, sum(BAA2 * tDI2, na.rm = TRUE), 0),
      #           CHNG_BAA = (CURR_BAA - PREV_BAA) / REMPER,
      #           TPA_RATE = CHNG_TPA / PREV_TPA,
      #           BAA_RATE = CHNG_BAA / PREV_BAA,
      #           x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
      #           y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
      #           M = sqrt(x^2 + y^2),
      #           SUST_INDEX = if_else(x < 0, -M, M),
      #           nStems = length(which(tDI == 1))) %>%
      # ungroup() %>%
      # select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE, PREV_TPA, CURR_TPA, PREV_BAA, CURR_BAA, nStems, nLive)

    # t <- data %>%
    #   mutate(YEAR = MEASYEAR) %>%
    #   distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
    #   # Compute estimates at plot level
    #   group_by(.dots = grpBy, PLT_CN) %>%
    #   summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
    #             REMPER = first(REMPER),
    #             PREV_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
    #             CURR_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ[ONEORTWO == 2] * tDI[ONEORTWO == 2], na.rm = TRUE), 0),
    #             CHNG_TPA = (CURR_TPA - PREV_TPA) / REMPER,
    #             PREV_BAA = if_else(nLive >= minLive, sum(BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
    #             CURR_BAA = if_else(nLive >= minLive, sum(BAA[ONEORTWO == 2] * tDI[ONEORTWO == 2], na.rm = TRUE), 0),
    #             CHNG_BAA = (CURR_BAA - PREV_BAA) / REMPER,
    #             TPA_RATE = CHNG_TPA / PREV_TPA,
    #             BAA_RATE = CHNG_BAA / PREV_BAA,
    #             x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
    #             y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
    #             M = sqrt(x^2 + y^2),
    #             SUST_INDEX = if_else(x < 0, -M, M),
    #             nStems = length(which(tDI == 1))) %>%
    #   ungroup() %>%
    #   select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE, PREV_TPA, CURR_TPA, PREV_BAA, CURR_BAA, nStems, nLive)

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
    # ### Plot-level estimates
    # a <- data %>%
    #   ## Will be lots of trees here, so CONDPROP listed multiple times
    #   ## Adding PROP_BASIS so we can handle adjustment factors at strata level
    #   distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
    #   group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
    #   summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
    #             plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
    #   left_join(select(a_ga, PLT_CN, PROP_BASIS, aGrpBy, fa_ga, plotIn_ga), by = c('PLT_CN', 'PROP_BASIS', aGrpBy))


    ### Compute total TREES in domain of interest
    t <- data %>%
      distinct(PLT_CN, TRE_CN, ONEORTWO, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
      summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
                REMPER = first(REMPER),
                MEASYEAR = first(MEASYEAR),
                PREV_TPA = if_else(nLive >= minLive, sum(-TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
                PREV_BAA = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
                CHNG_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ * tDI, na.rm = TRUE), 0),
                CHNG_BAA = if_else(nLive >= minLive, sum(BAA * tDI, na.rm = TRUE), 0),
                ## Disturbances
                BUG_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & BUG == 1] * tDI[ONEORTWO == 1 & BUG == 1], na.rm = TRUE), 0),
                DISEASE_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & DISEASE == 1] * tDI[ONEORTWO == 1 & DISEASE == 1], na.rm = TRUE), 0),
                FIRE_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & FIRE == 1] * tDI[ONEORTWO == 1 & FIRE == 1], na.rm = TRUE), 0),
                ANIMAL_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & ANIMAL == 1] * tDI[ONEORTWO == 1 & ANIMAL == 1], na.rm = TRUE), 0),
                WEATHER_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & WEATHER == 1] * tDI[ONEORTWO == 1 & WEATHER == 1], na.rm = TRUE), 0),
                VEG_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & VEG == 1] * tDI[ONEORTWO == 1 & VEG == 1], na.rm = TRUE), 0),
                UNKNOWN_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & UNKNOWN == 1] * tDI[ONEORTWO == 1 & UNKNOWN == 1], na.rm = TRUE), 0),
                SILV_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & SILV == 1] * tDI[ONEORTWO == 1 & SILV == 1], na.rm = TRUE), 0),
                BUG_YEAR = if_else(nLive >= minLive, mean(-MORTYR[BUG == 1] * tDI[BUG == 1], na.rm = TRUE), 0),
                DISEASE_YEAR = if_else(nLive >= minLive, mean(-MORTYR[DISEASE == 1] * tDI[DISEASE == 1], na.rm = TRUE), 0),
                FIRE_YEAR = if_else(nLive >= minLive, mean(-MORTYR[FIRE == 1] * tDI[FIRE == 1], na.rm = TRUE), 0),
                ANIMAL_YEAR = if_else(nLive >= minLive, mean(-MORTYR[ANIMAL == 1] * tDI[ANIMAL == 1], na.rm = TRUE), 0),
                WEATHER_YEAR = if_else(nLive >= minLive, mean(-MORTYR[WEATHER == 1] * tDI[WEATHER == 1], na.rm = TRUE), 0),
                VEG_YEAR = if_else(nLive >= minLive, mean(-MORTYR[VEG == 1] * tDI[VEG == 1], na.rm = TRUE), 0),
                UNKNOWN_YEAR = if_else(nLive >= minLive, mean(-MORTYR[UNKNOWN == 1] * tDI[UNKNOWN == 1], na.rm = TRUE), 0),
                SILV_YEAR = if_else(nLive >= minLive, mean(-MORTYR[SILV == 1] * tDI[SILV == 1], na.rm = TRUE), 0),

                plotIn = if_else(sum(tDI) > 0, 1, 0))
  }

  pltOut <- list(a = NULL, t = t)
  return(pltOut)

}

sustIndexHelper2 <- function(x, popState, t, grpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR

  }

  ######## ------------------ TREE ESTIMATES + CV

  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD, REMPER)), by = 'PLT_CN') %>%
    ## Need this for covariance later on
    #Add adjustment factors
    #Add adjustment factors
    # mutate(tAdj = case_when(
    #   ## When NA, stay NA
    #   is.na(SUBPTYP_GRM) ~ NA_real_,
    #   ## If the proportion was measured for a macroplot,
    #   ## use the macroplot value
    #   SUBPTYP_GRM == 0 ~ 0,
    #   SUBPTYP_GRM == 1 ~ as.numeric(ADJ_FACTOR_SUBP),
    #   SUBPTYP_GRM == 2 ~ as.numeric(ADJ_FACTOR_MICR),
    #   SUBPTYP_GRM == 3 ~ as.numeric(ADJ_FACTOR_MACR)),
    #   ct = CHNG_TPA * tAdj,
    #   cb = CHNG_BAA * tAdj,
    #   pt = PREV_TPA * tAdj,
    #   pb = PREV_BAA * tAdj) %>%
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
      pb = PREV_BAA * tAdj,
    ## Adjusted for time since disturbance
      # bug = BUG_RATE * tAdj / (MEASYEAR - BUG_YEAR),
      # disease = DISEASE_RATE * tAdj / (MEASYEAR - DISEASE_YEAR),
      # fire = FIRE_RATE * tAdj / (MEASYEAR - FIRE_YEAR),
      # animal = ANIMAL_RATE * tAdj / (MEASYEAR - ANIMAL_YEAR),
      # weather = WEATHER_RATE * tAdj / (MEASYEAR - WEATHER_YEAR),
      # veg = VEG_RATE * tAdj / (MEASYEAR - VEG_YEAR),
      # unknown = UNKNOWN_RATE * tAdj / (MEASYEAR - UNKNOWN_YEAR),
      # silv = SILV_RATE * tAdj / (MEASYEAR - SILV_YEAR)
    # bug = BUG_RATE * tAdj / (REMPER - (MEASYEAR - BUG_YEAR)),
    # disease = DISEASE_RATE * tAdj / (REMPER - (MEASYEAR - DISEASE_YEAR)),
    # fire = FIRE_RATE * tAdj / (REMPER - (MEASYEAR - FIRE_YEAR)),
    # animal = ANIMAL_RATE * tAdj / (REMPER - (MEASYEAR - ANIMAL_YEAR)),
    # weather = WEATHER_RATE * tAdj / (REMPER - (MEASYEAR - WEATHER_YEAR)),
    # veg = VEG_RATE * tAdj / (REMPER - (MEASYEAR - VEG_YEAR)),
    # unknown = UNKNOWN_RATE * tAdj / (REMPER - (MEASYEAR - UNKNOWN_YEAR)),
    # silv = SILV_RATE * tAdj / (REMPER - (MEASYEAR - SILV_YEAR))
    bug = BUG_RATE * tAdj / REMPER,
    disease = DISEASE_RATE * tAdj / REMPER,
    fire = FIRE_RATE * tAdj / REMPER,
    animal = ANIMAL_RATE * tAdj / REMPER,
    weather = WEATHER_RATE * tAdj / REMPER,
    veg = VEG_RATE * tAdj / REMPER,
    unknown = UNKNOWN_RATE * tAdj / REMPER,
    silv = SILV_RATE * tAdj / REMPER
    ) %>%
    ## Computing change
    mutate(ct = (ct) / REMPER,
           cb = (cb) / REMPER) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(ctPlot = sum(ct, na.rm = TRUE),
              cbPlot = sum(cb, na.rm = TRUE),
              ptPlot = sum(pt, na.rm = TRUE),
              pbPlot = sum(pb, na.rm = TRUE),
              bugPlot = sum(bug, na.rm = TRUE),
              diseasePlot = sum(disease, na.rm = TRUE),
              firePlot = sum(fire, na.rm = TRUE),
              animalPlot = sum(animal, na.rm = TRUE),
              weatherPlot = sum(weather, na.rm = TRUE),
              vegPlot = sum(veg, na.rm = TRUE),
              unPlot = sum(unknown, na.rm = TRUE),
              silvPlot = sum(silv, na.rm = TRUE),
              plotIn_t = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              ctStrat = mean(ctPlot * r_t, na.rm = TRUE),
              cbStrat = mean(cbPlot * r_t, na.rm = TRUE),
              ptStrat = mean(ptPlot * r_t, na.rm = TRUE),
              pbStrat = mean(pbPlot * r_t, na.rm = TRUE),
              bugStrat = mean(bugPlot * r_t, na.rm = TRUE),
              diseaseStrat = mean(diseasePlot * r_t, na.rm = TRUE),
              fireStrat = mean(firePlot * r_t, na.rm = TRUE),
              animalStrat = mean(animalPlot * r_t, na.rm = TRUE),
              weatherStrat = mean(weatherPlot * r_t, na.rm = TRUE),
              vegStrat = mean(vegPlot * r_t, na.rm = TRUE),
              unStrat = mean(unPlot * r_t, na.rm = TRUE),
              silvStrat = mean(silvPlot * r_t, na.rm = TRUE),

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
              bugv = stratVar(ESTN_METHOD, bugPlot, bugStrat, ndif, a, nh),
              diseasev = stratVar(ESTN_METHOD, diseasePlot, diseaseStrat, ndif, a, nh),
              firev = stratVar(ESTN_METHOD, firePlot, fireStrat, ndif, a, nh),
              animalv = stratVar(ESTN_METHOD, animalPlot, animalStrat, ndif, a, nh),
              weatherv = stratVar(ESTN_METHOD, weatherPlot, weatherStrat, ndif, a, nh),
              vegv = stratVar(ESTN_METHOD, vegPlot, vegStrat, ndif, a, nh),
              unv = stratVar(ESTN_METHOD, unPlot, unStrat, ndif, a, nh),
              silvv = stratVar(ESTN_METHOD, silvPlot, silvStrat, ndif, a, nh),

              # Strata level covariances
              cvStrat_ct = stratVar(ESTN_METHOD, ctPlot, ctStrat, ndif, a, nh, ptPlot, ptStrat),
              cvStrat_cb = stratVar(ESTN_METHOD, cbPlot, cbStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_bug = stratVar(ESTN_METHOD, bugPlot, bugStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_disease = stratVar(ESTN_METHOD, diseasePlot, diseaseStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_fire = stratVar(ESTN_METHOD, firePlot, fireStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_animal = stratVar(ESTN_METHOD, animalPlot, animalStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_weather = stratVar(ESTN_METHOD, weatherPlot, weatherStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_veg = stratVar(ESTN_METHOD, vegPlot, vegStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_un = stratVar(ESTN_METHOD, unPlot, unStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_silv = stratVar(ESTN_METHOD, silvPlot, silvStrat, ndif, a, nh, pbPlot, pbStrat),

    ) %>%

    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
              cbEst = unitMean(ESTN_METHOD, a, nh, w, cbStrat),
              ptEst = unitMean(ESTN_METHOD, a, nh, w, ptStrat),
              pbEst = unitMean(ESTN_METHOD, a, nh, w, pbStrat),
              bugEst = unitMean(ESTN_METHOD, a, nh, w, bugStrat),
              diseaseEst = unitMean(ESTN_METHOD, a, nh, w, diseaseStrat),
              fireEst = unitMean(ESTN_METHOD, a, nh, w, fireStrat),
              animalEst = unitMean(ESTN_METHOD, a, nh, w, animalStrat),
              weatherEst = unitMean(ESTN_METHOD, a, nh, w, weatherStrat),
              vegEst = unitMean(ESTN_METHOD, a, nh, w, vegStrat),
              unEst = unitMean(ESTN_METHOD, a, nh, w, unStrat),
              silvEst = unitMean(ESTN_METHOD, a, nh, w, silvStrat),
              nh = first(nh),
              # Estimation of unit variance
              ctVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ctv, ctStrat, ctEst),
              cbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cbv, cbStrat, cbEst),
              ptVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ptv, ptStrat, ptEst),
              pbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pbv, pbStrat, pbEst),
              bugVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bugv, bugStrat, bugEst),
              diseaseVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, diseasev, diseaseStrat, diseaseEst),
              fireVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, firev, fireStrat, fireEst),
              animalVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, animalv, animalStrat, animalEst),
              weatherVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, weatherv, weatherStrat, weatherEst),
              vegVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, vegv, vegStrat, vegEst),
              unVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, unv, unStrat, unEst),
              silvVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, silvv, silvStrat, silvEst),
              ## Covariances
              cvEst_ct = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ct, ctStrat, ctEst, ptStrat, ptEst),
              cvEst_cb = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_cb, cbStrat, cbEst, pbStrat, pbEst),
              cvEst_bug = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_bug, bugStrat, bugEst, pbStrat, pbEst),
              cvEst_disease = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_disease, diseaseStrat, diseaseEst, pbStrat, pbEst),
              cvEst_fire = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_fire, fireStrat, fireEst, pbStrat, pbEst),
              cvEst_animal = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_animal, animalStrat, animalEst, pbStrat, pbEst),
              cvEst_weather = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_weather, weatherStrat, weatherEst, pbStrat, pbEst),
              cvEst_veg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_veg, vegStrat, vegEst, pbStrat, pbEst),
              cvEst_un = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_un, unStrat, unEst, pbStrat, pbEst),
              cvEst_silv = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_silv, silvStrat, silvEst, pbStrat, pbEst),
              plotIn_t = sum(plotIn_t, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = NULL)

  return(out)
}

siHelper1 <- function(x, plts, db, grpBy, byPlot, minLive){

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
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR',
                            'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp', 'DESIGNCD',
                            'drought_sev', 'wet_sev', 'all_sev', 'grow_drought_sev', 'grow_wet_sev', 'grow_all_sev')) %>%
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & DESIGNCD == 1) %>%

    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    ## AGENTCD at remeasurement, died during the measurement interval
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'TPA_UNADJ', 'DIA', 'AGENTCD', 'MORTYR')), by = c('PLT_CN', 'CONDID')) %>%
    #left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', DIA_BEGIN, DIA_END)), by = c('TRE_CN')) %>%

    left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p', 'DESIGNCD', 'PLOT_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('2', '1')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('2', '1')) %>%
    left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD', 'TPA_UNADJ', 'DIA')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%
    #left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%

    mutate_if(is.factor,
              as.character)

  ## Comprehensive indicator function -- w/ growth accounting
  data$tDI2 <- data$landD2 * data$aD_p2 * data$aD_c2 * data$tD2 * data$typeD2 * data$sp2 #*
  #if_else(data$PLOT_STATUS_CD1 == 1 & data$PLOT_STATUS_CD2 == 1, 1, 0)
  data$tDI1 <- data$landD1 * data$aD_p1 * data$aD_c1 * data$tD1 * data$typeD1 * data$sp1 #*
  #if_else(data$PLOT_STATUS_CD1 == 1 & data$PLOT_STATUS_CD2 == 1, 1, 0)

  ## PREVIOUS and CURRENT attributes
  data <- data %>%
    mutate(TPA_UNADJ1 = TPA_UNADJ1,
           TPA_UNADJ2 = TPA_UNADJ2,
           BAA1 = basalArea(DIA1) * TPA_UNADJ1,
           BAA2 = basalArea(DIA2) * TPA_UNADJ2,
           BUG = if_else(AGENTCD == 10, 1, 0),
           DISEASE = if_else(AGENTCD == 20, 1, 0),
           FIRE = if_else(AGENTCD == 30, 1, 0),
           ANIMAL = if_else(AGENTCD == 40, 1, 0),
           WEATHER = if_else(AGENTCD == 50, 1, 0),
           VEG = if_else(AGENTCD == 60, 1, 0),
           UNKNOWN = if_else(AGENTCD == 70, 1, 0),
           SILV = if_else(AGENTCD == 80, 1, 0),
    ) #%>%

  ## Just what we need
  data <- data %>%
    select(PLT_CN, TRE_CN, SUBP, CONDID, TREE, CONDPROP_UNADJ,
           MEASYEAR, MACRO_BREAKPOINT_DIA, PROP_BASIS,
           BUG, DISEASE, FIRE, ANIMAL, WEATHER, VEG, UNKNOWN, SILV, MORTYR,
           drought_sev, wet_sev, all_sev, grow_drought_sev, grow_wet_sev, grow_all_sev,
           REMPER, PLOT_STATUS_CD1, PLOT_STATUS_CD2,
           one_of(str_c(grpP,1), str_c(grpC,1), str_c(grpT,1),
                  str_c(grpP,2), str_c(grpC,2), str_c(grpT,2)),
           tDI1, tDI2, #SUBPTYP_GRM1, SUBPTYP_GRM2,
           DIA1, DIA2, BAA1, BAA2, TPA_UNADJ1, TPA_UNADJ2) %>%
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
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR'))


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
                REMPER = first(REMPER),
                PREV_TPA = if_else(nLive >= minLive, sum(-TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
                PREV_BAA = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
                CHNG_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ * tDI, na.rm = TRUE), 0) / REMPER,
                CHNG_BAA = if_else(nLive >= minLive, sum(BAA * tDI, na.rm = TRUE), 0) / REMPER,
                TPA_RATE = CHNG_TPA / PREV_TPA,
                BAA_RATE = CHNG_BAA / PREV_BAA,
                ## Disturbances
                BUG_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & BUG == 1] * tDI[ONEORTWO == 1 & BUG == 1], na.rm = TRUE), 0) / PREV_BAA,
                DISEASE_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & DISEASE == 1] * tDI[ONEORTWO == 1 & DISEASE == 1], na.rm = TRUE), 0) / PREV_BAA,
                FIRE_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & FIRE == 1] * tDI[ONEORTWO == 1 & FIRE == 1], na.rm = TRUE), 0) / PREV_BAA,
                ANIMAL_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & ANIMAL == 1] * tDI[ONEORTWO == 1 & ANIMAL == 1], na.rm = TRUE), 0) / PREV_BAA,
                WEATHER_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & WEATHER == 1] * tDI[ONEORTWO == 1 & WEATHER == 1], na.rm = TRUE), 0) / PREV_BAA,
                VEG_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & VEG == 1] * tDI[ONEORTWO == 1 & VEG == 1], na.rm = TRUE), 0) / PREV_BAA,
                UNKNOWN_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & UNKNOWN == 1] * tDI[ONEORTWO == 1 & UNKNOWN == 1], na.rm = TRUE), 0) / PREV_BAA,
                SILV_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & SILV == 1] * tDI[ONEORTWO == 1 & SILV == 1], na.rm = TRUE), 0) / PREV_BAA,
                DROUGHT_SEV = first(drought_sev),
                WET_SEV = first(wet_sev),
                ALL_SEV = first(all_sev),
                GROW_DROUGHT_SEV = first(grow_drought_sev),
                GROW_WET_SEV = first(grow_wet_sev),
                GROW_ALL_SEV = first(grow_all_sev),
                x = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$x,
                y = projectPnts(TPA_RATE, BAA_RATE, 1, 0)$y,
                M = sqrt(x^2 + y^2),
                SUST_INDEX = if_else(x < 0, -M, M),
                nStems = length(which(tDI == 1))) %>%
      ungroup() %>%
      select(grpBy, SUST_INDEX, TPA_RATE, BAA_RATE, BUG_RATE,
             DISEASE_RATE, FIRE_RATE, ANIMAL_RATE, WEATHER_RATE, VEG_RATE,
             UNKNOWN_RATE, SILV_RATE, DROUGHT_SEV, WET_SEV, ALL_SEV,
             GROW_DROUGHT_SEV, GROW_WET_SEV, GROW_ALL_SEV,
             PREV_TPA, PREV_BAA, nStems, nLive)

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
                DROUGHT_SEV = first(drought_sev),
                WET_SEV = first(wet_sev),
                ALL_SEV = first(all_sev),
                GROW_DROUGHT_SEV = first(grow_drought_sev),
                GROW_WET_SEV = first(grow_wet_sev),
                GROW_ALL_SEV = first(grow_all_sev),
                CONDPROP_UNADJ = first(CONDPROP_UNADJ)) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = grpBy) %>%
      summarize(fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                DROUGHT_SEV = sum(CONDPROP_UNADJ * DROUGHT_SEV * aDI, na.rm = TRUE),
                WET_SEV = sum(CONDPROP_UNADJ * WET_SEV * aDI, na.rm = TRUE),
                ALL_SEV = sum(CONDPROP_UNADJ * ALL_SEV * aDI, na.rm = TRUE),
                GROW_DROUGHT_SEV = sum(CONDPROP_UNADJ * GROW_DROUGHT_SEV * aDI, na.rm = TRUE),
                GROW_WET_SEV = sum(CONDPROP_UNADJ * GROW_WET_SEV * aDI, na.rm = TRUE),
                GROW_ALL_SEV = sum(CONDPROP_UNADJ * GROW_ALL_SEV * aDI, na.rm = TRUE))

    ### Compute total TREES in domain of interest
    t <- data %>%
      distinct(PLT_CN, TRE_CN, ONEORTWO, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
      summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
                REMPER = first(REMPER),
                MEASYEAR = first(MEASYEAR),
                PREV_TPA = if_else(nLive >= minLive, sum(-TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
                PREV_BAA = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE), 0),
                CHNG_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ * tDI, na.rm = TRUE), 0),
                CHNG_BAA = if_else(nLive >= minLive, sum(BAA * tDI, na.rm = TRUE), 0),
                ## Disturbances
                BUG_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & BUG == 1] * tDI[ONEORTWO == 1 & BUG == 1], na.rm = TRUE), 0),
                DISEASE_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & DISEASE == 1] * tDI[ONEORTWO == 1 & DISEASE == 1], na.rm = TRUE), 0),
                FIRE_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & FIRE == 1] * tDI[ONEORTWO == 1 & FIRE == 1], na.rm = TRUE), 0),
                ANIMAL_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & ANIMAL == 1] * tDI[ONEORTWO == 1 & ANIMAL == 1], na.rm = TRUE), 0),
                WEATHER_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & WEATHER == 1] * tDI[ONEORTWO == 1 & WEATHER == 1], na.rm = TRUE), 0),
                VEG_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & VEG == 1] * tDI[ONEORTWO == 1 & VEG == 1], na.rm = TRUE), 0),
                UNKNOWN_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & UNKNOWN == 1] * tDI[ONEORTWO == 1 & UNKNOWN == 1], na.rm = TRUE), 0),
                SILV_RATE = if_else(nLive >= minLive, sum(-BAA[ONEORTWO == 1 & SILV == 1] * tDI[ONEORTWO == 1 & SILV == 1], na.rm = TRUE), 0),
                BUG_YEAR = if_else(nLive >= minLive, mean(-MORTYR[BUG == 1] * tDI[BUG == 1], na.rm = TRUE), 0),
                DISEASE_YEAR = if_else(nLive >= minLive, mean(-MORTYR[DISEASE == 1] * tDI[DISEASE == 1], na.rm = TRUE), 0),
                FIRE_YEAR = if_else(nLive >= minLive, mean(-MORTYR[FIRE == 1] * tDI[FIRE == 1], na.rm = TRUE), 0),
                ANIMAL_YEAR = if_else(nLive >= minLive, mean(-MORTYR[ANIMAL == 1] * tDI[ANIMAL == 1], na.rm = TRUE), 0),
                WEATHER_YEAR = if_else(nLive >= minLive, mean(-MORTYR[WEATHER == 1] * tDI[WEATHER == 1], na.rm = TRUE), 0),
                VEG_YEAR = if_else(nLive >= minLive, mean(-MORTYR[VEG == 1] * tDI[VEG == 1], na.rm = TRUE), 0),
                UNKNOWN_YEAR = if_else(nLive >= minLive, mean(-MORTYR[UNKNOWN == 1] * tDI[UNKNOWN == 1], na.rm = TRUE), 0),
                SILV_YEAR = if_else(nLive >= minLive, mean(-MORTYR[SILV == 1] * tDI[SILV == 1], na.rm = TRUE), 0),
                plotIn = if_else(sum(tDI) > 0, 1, 0)) #%>%
      #left_join()
  }

  pltOut <- list(a = a, t = t)
  return(pltOut)

}

siHelper2 <- function(x, popState, a, t, grpBy, method){

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
      DROUGHT_SEV =  DROUGHT_SEV * aAdj,
      WET_SEV = WET_SEV * aAdj,
      ALL_SEV = ALL_SEV * aAdj,
      GROW_DROUGHT_SEV = GROW_DROUGHT_SEV * aAdj,
      GROW_WET_SEV = GROW_WET_SEV * aAdj,
      GROW_ALL_SEV = GROW_ALL_SEV * aAdj,
      fa = fa * aAdj) %>%
    ungroup()

  ## Strata level estimates
  tEst <- t %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD, REMPER)), by = 'PLT_CN') %>%
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
    pb = PREV_BAA * tAdj,
    ## Adjusted for time since disturbance
    # bug = BUG_RATE * tAdj / (MEASYEAR - BUG_YEAR),
    # disease = DISEASE_RATE * tAdj / (MEASYEAR - DISEASE_YEAR),
    # fire = FIRE_RATE * tAdj / (MEASYEAR - FIRE_YEAR),
    # animal = ANIMAL_RATE * tAdj / (MEASYEAR - ANIMAL_YEAR),
    # weather = WEATHER_RATE * tAdj / (MEASYEAR - WEATHER_YEAR),
    # veg = VEG_RATE * tAdj / (MEASYEAR - VEG_YEAR),
    # unknown = UNKNOWN_RATE * tAdj / (MEASYEAR - UNKNOWN_YEAR),
    # silv = SILV_RATE * tAdj / (MEASYEAR - SILV_YEAR)
    # bug = BUG_RATE * tAdj / (REMPER - (MEASYEAR - BUG_YEAR)),
    # disease = DISEASE_RATE * tAdj / (REMPER - (MEASYEAR - DISEASE_YEAR)),
    # fire = FIRE_RATE * tAdj / (REMPER - (MEASYEAR - FIRE_YEAR)),
    # animal = ANIMAL_RATE * tAdj / (REMPER - (MEASYEAR - ANIMAL_YEAR)),
    # weather = WEATHER_RATE * tAdj / (REMPER - (MEASYEAR - WEATHER_YEAR)),
    # veg = VEG_RATE * tAdj / (REMPER - (MEASYEAR - VEG_YEAR)),
    # unknown = UNKNOWN_RATE * tAdj / (REMPER - (MEASYEAR - UNKNOWN_YEAR)),
    # silv = SILV_RATE * tAdj / (REMPER - (MEASYEAR - SILV_YEAR))
    bug = BUG_RATE * tAdj / REMPER,
    disease = DISEASE_RATE * tAdj / REMPER,
    fire = FIRE_RATE * tAdj / REMPER,
    animal = ANIMAL_RATE * tAdj / REMPER,
    weather = WEATHER_RATE * tAdj / REMPER,
    veg = VEG_RATE * tAdj / REMPER,
    unknown = UNKNOWN_RATE * tAdj / REMPER,
    silv = SILV_RATE * tAdj / REMPER
  ) %>%
    ## Computing change
    mutate(ct = (ct) / REMPER,
           cb = (cb) / REMPER) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(ctPlot = sum(ct, na.rm = TRUE),
              cbPlot = sum(cb, na.rm = TRUE),
              ptPlot = sum(pt, na.rm = TRUE),
              pbPlot = sum(pb, na.rm = TRUE),
              bugPlot = sum(bug, na.rm = TRUE),
              diseasePlot = sum(disease, na.rm = TRUE),
              firePlot = sum(fire, na.rm = TRUE),
              animalPlot = sum(animal, na.rm = TRUE),
              weatherPlot = sum(weather, na.rm = TRUE),
              vegPlot = sum(veg, na.rm = TRUE),
              unPlot = sum(unknown, na.rm = TRUE),
              silvPlot = sum(silv, na.rm = TRUE),
              plotIn_t = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    left_join(select(aAdj, PLT_CN, grpBy, fa, DROUGHT_SEV, WET_SEV, ALL_SEV, GROW_DROUGHT_SEV, GROW_WET_SEV, GROW_ALL_SEV), by = c('PLT_CN', grpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              ctStrat = mean(ctPlot * r_t, na.rm = TRUE),
              cbStrat = mean(cbPlot * r_t, na.rm = TRUE),
              ptStrat = mean(ptPlot * r_t, na.rm = TRUE),
              pbStrat = mean(pbPlot * r_t, na.rm = TRUE),
              bugStrat = mean(bugPlot * r_t, na.rm = TRUE),
              diseaseStrat = mean(diseasePlot * r_t, na.rm = TRUE),
              fireStrat = mean(firePlot * r_t, na.rm = TRUE),
              animalStrat = mean(animalPlot * r_t, na.rm = TRUE),
              weatherStrat = mean(weatherPlot * r_t, na.rm = TRUE),
              vegStrat = mean(vegPlot * r_t, na.rm = TRUE),
              unStrat = mean(unPlot * r_t, na.rm = TRUE),
              silvStrat = mean(silvPlot * r_t, na.rm = TRUE),
              faStrat = mean(fa * r_t, na.rm = TRUE),
              dStrat = mean(DROUGHT_SEV * r_t, na.rm = TRUE),
              wStrat = mean(WET_SEV * r_t, na.rm = TRUE),
              aStrat = mean(ALL_SEV * r_t, na.rm = TRUE),
              gdStrat = mean(GROW_DROUGHT_SEV * r_t, na.rm = TRUE),
              gwStrat = mean(GROW_WET_SEV * r_t, na.rm = TRUE),
              gaStrat = mean(GROW_ALL_SEV * r_t, na.rm = TRUE),
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
              bugv = stratVar(ESTN_METHOD, bugPlot, bugStrat, ndif, a, nh),
              diseasev = stratVar(ESTN_METHOD, diseasePlot, diseaseStrat, ndif, a, nh),
              firev = stratVar(ESTN_METHOD, firePlot, fireStrat, ndif, a, nh),
              animalv = stratVar(ESTN_METHOD, animalPlot, animalStrat, ndif, a, nh),
              weatherv = stratVar(ESTN_METHOD, weatherPlot, weatherStrat, ndif, a, nh),
              vegv = stratVar(ESTN_METHOD, vegPlot, vegStrat, ndif, a, nh),
              unv = stratVar(ESTN_METHOD, unPlot, unStrat, ndif, a, nh),
              silvv = stratVar(ESTN_METHOD, silvPlot, silvStrat, ndif, a, nh),
              fav = stratVar(ESTN_METHOD, fa, faStrat, ndif, a, nh),
              dv = stratVar(ESTN_METHOD, DROUGHT_SEV, dStrat, ndif, a, nh),
              wv = stratVar(ESTN_METHOD, WET_SEV, wStrat, ndif, a, nh),
              av = stratVar(ESTN_METHOD, ALL_SEV, aStrat, ndif, a, nh),
              gdv = stratVar(ESTN_METHOD, GROW_DROUGHT_SEV, gdStrat, ndif, a, nh),
              gwv = stratVar(ESTN_METHOD, GROW_WET_SEV, gwStrat, ndif, a, nh),
              gav = stratVar(ESTN_METHOD, GROW_ALL_SEV, gaStrat, ndif, a, nh),

              # Strata level covariances
              cvStrat_ct = stratVar(ESTN_METHOD, ctPlot, ctStrat, ndif, a, nh, ptPlot, ptStrat),
              cvStrat_cb = stratVar(ESTN_METHOD, cbPlot, cbStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_bug = stratVar(ESTN_METHOD, bugPlot, bugStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_disease = stratVar(ESTN_METHOD, diseasePlot, diseaseStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_fire = stratVar(ESTN_METHOD, firePlot, fireStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_animal = stratVar(ESTN_METHOD, animalPlot, animalStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_weather = stratVar(ESTN_METHOD, weatherPlot, weatherStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_veg = stratVar(ESTN_METHOD, vegPlot, vegStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_un = stratVar(ESTN_METHOD, unPlot, unStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_silv = stratVar(ESTN_METHOD, silvPlot, silvStrat, ndif, a, nh, pbPlot, pbStrat),
              cvStrat_d = stratVar(ESTN_METHOD, DROUGHT_SEV, dStrat, ndif, a, nh, fa, faStrat),
              cvStrat_w = stratVar(ESTN_METHOD, WET_SEV, wStrat, ndif, a, nh, fa, faStrat),
              cvStrat_a = stratVar(ESTN_METHOD, ALL_SEV, aStrat, ndif, a, nh, fa, faStrat),
              cvStrat_gd = stratVar(ESTN_METHOD, GROW_DROUGHT_SEV, gdStrat, ndif, a, nh, fa, faStrat),
              cvStrat_gw = stratVar(ESTN_METHOD, GROW_WET_SEV, gwStrat, ndif, a, nh, fa, faStrat),
              cvStrat_ga = stratVar(ESTN_METHOD, GROW_ALL_SEV, gaStrat, ndif, a, nh, fa, faStrat),

    ) %>%

    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
              cbEst = unitMean(ESTN_METHOD, a, nh, w, cbStrat),
              ptEst = unitMean(ESTN_METHOD, a, nh, w, ptStrat),
              pbEst = unitMean(ESTN_METHOD, a, nh, w, pbStrat),
              bugEst = unitMean(ESTN_METHOD, a, nh, w, bugStrat),
              diseaseEst = unitMean(ESTN_METHOD, a, nh, w, diseaseStrat),
              fireEst = unitMean(ESTN_METHOD, a, nh, w, fireStrat),
              animalEst = unitMean(ESTN_METHOD, a, nh, w, animalStrat),
              weatherEst = unitMean(ESTN_METHOD, a, nh, w, weatherStrat),
              vegEst = unitMean(ESTN_METHOD, a, nh, w, vegStrat),
              unEst = unitMean(ESTN_METHOD, a, nh, w, unStrat),
              silvEst = unitMean(ESTN_METHOD, a, nh, w, silvStrat),
              faEst = unitMean(ESTN_METHOD, a, nh, w, faStrat),
              dEst = unitMean(ESTN_METHOD, a, nh, w, dStrat),
              wEst = unitMean(ESTN_METHOD, a, nh, w, wStrat),
              aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
              gdEst = unitMean(ESTN_METHOD, a, nh, w, gdStrat),
              gwEst = unitMean(ESTN_METHOD, a, nh, w, gwStrat),
              gaEst = unitMean(ESTN_METHOD, a, nh, w, gaStrat),

              nh = first(nh),
              # Estimation of unit variance
              ctVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ctv, ctStrat, ctEst),
              cbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cbv, cbStrat, cbEst),
              ptVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ptv, ptStrat, ptEst),
              pbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pbv, pbStrat, pbEst),
              bugVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bugv, bugStrat, bugEst),
              diseaseVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, diseasev, diseaseStrat, diseaseEst),
              fireVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, firev, fireStrat, fireEst),
              animalVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, animalv, animalStrat, animalEst),
              weatherVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, weatherv, weatherStrat, weatherEst),
              vegVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, vegv, vegStrat, vegEst),
              unVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, unv, unStrat, unEst),
              silvVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, silvv, silvStrat, silvEst),
              faVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, fav, faStrat, faEst),
              dVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, dv, dStrat, dEst),
              wVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, wv, wStrat, wEst),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, av, aStrat, aEst),
              gdVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, gdv, gdStrat, gdEst),
              gwVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, gwv, gwStrat, gwEst),
              gaVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, gav, gaStrat, gaEst),

              ## Covariances
              cvEst_ct = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ct, ctStrat, ctEst, ptStrat, ptEst),
              cvEst_cb = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_cb, cbStrat, cbEst, pbStrat, pbEst),
              cvEst_bug = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_bug, bugStrat, bugEst, pbStrat, pbEst),
              cvEst_disease = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_disease, diseaseStrat, diseaseEst, pbStrat, pbEst),
              cvEst_fire = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_fire, fireStrat, fireEst, pbStrat, pbEst),
              cvEst_animal = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_animal, animalStrat, animalEst, pbStrat, pbEst),
              cvEst_weather = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_weather, weatherStrat, weatherEst, pbStrat, pbEst),
              cvEst_veg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_veg, vegStrat, vegEst, pbStrat, pbEst),
              cvEst_un = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_un, unStrat, unEst, pbStrat, pbEst),
              cvEst_silv = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_silv, silvStrat, silvEst, pbStrat, pbEst),

              cvEst_d = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_d, dStrat, dEst, faStrat, faEst),
              cvEst_w = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_w, wStrat, wEst, faStrat, faEst),
              cvEst_a = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_a, aStrat, aEst, faStrat, faEst),
              cvEst_gd = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_gd, gdStrat, gdEst, faStrat, faEst),
              cvEst_gw = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_gw, gwStrat, gwEst, faStrat, faEst),
              cvEst_ga = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ga, gaStrat, gaEst, faStrat, faEst),


              plotIn_t = sum(plotIn_t, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = NULL)

  return(out)
}

diMortHelper1_old1 <- function(x, plts, db, grpBy, aGrpBy, byPlot, minLive){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)

  # Only subplots from cond change matrix
  db$SUBP_COND_CHNG_MTRX <- filter(db$SUBP_COND_CHNG_MTRX, SUBPTYP == 1)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', 'PREV_PLT_CN', 'REMPER', grpP, 'aD_p', 'sp')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'PREVCOND', 'TRE_CN', 'PREV_TRE_CN', 'SUBP', 'TREE', grpT, 'tD', 'typeD')), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$TREE_GRM_COMPONENT, c('TRE_CN', 'SUBPTYP_GRM', 'TPAGROW_UNADJ', 'TPARECR_UNADJ', 'TPAREMV_UNADJ', 'TPAMORT_UNADJ', 'COMPONENT', DIA_BEGIN, DIA_END)), by = c('TRE_CN')) %>%
    left_join(select(db$TREE_GRM_MIDPT, c('TRE_CN', 'DIA')), by = c('TRE_CN'), suffix = c('', '.mid')) %>%
    left_join(select(db$SUBP_COND_CHNG_MTRX, SUBP:SUBPTYP_PROP_CHNG), by = c('PLT_CN', 'CONDID'), suffix = c('', '.subp')) %>%
    left_join(select(db$COND, PLT_CN, CONDID, COND_STATUS_CD), by = c('PREV_PLT_CN.subp' = 'PLT_CN', 'PREVCOND.subp' = 'CONDID'), suffix = c('', '.chng')) %>%
    left_join(select(db$PLOT, c('PLT_CN', grpP, 'sp', 'aD_p')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', grpC, 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('', '.prev')) %>%
    left_join(select(db$TREE, c('TRE_CN', grpT, 'typeD', 'tD')), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('', '.prev')) %>%
    mutate_if(is.factor,
              as.character) %>%
    mutate(SUBPTYP_PROP_CHNG = SUBPTYP_PROP_CHNG * .25,
           TPAGROW_UNADJ = TPAGROW_UNADJ, ## NEEDS TO BE PREVIOUS
           TPAREMV_UNADJ = TPAREMV_UNADJ * REMPER,
           TPAMORT_UNADJ = TPAMORT_UNADJ * REMPER,
           TPARECR_UNADJ = TPARECR_UNADJ,
           aChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.chng == 1 & !is.null(CONDPROP_UNADJ), 1, 0),
           tChng = ifelse(COND_STATUS_CD == 1 & COND_STATUS_CD.prev == 1, 1, 0))

  # If previous attributes are unavailable for trees, default to current (otherwise we get NAs for early inventories)
  data$tD.prev <- ifelse(is.na(data$tD.prev), data$tD, data$tD.prev)
  data$typeD.prev <- ifelse(is.na(data$typeD.prev), data$typeD, data$typeD.prev)
  data$landD.prev <- ifelse(is.na(data$landD.prev), data$landD, data$landD.prev)
  data$aD_p.prev <- ifelse(is.na(data$aD_p.prev), data$aD_p, data$aD_p.prev)
  data$aD_c.prev <- ifelse(is.na(data$aD_c.prev), data$aD_c, data$aD_c.prev)
  data$sp.prev <- ifelse(is.na(data$sp.prev), data$sp, data$sp.prev)

  ## Comprehensive indicator function -- w/ growth accounting
  data$aDI_ga <- data$landD * data$aD_p * data$aD_c * data$sp * data$aChng
  data$tDI_ga <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev * data$tChng
  data$tDI_ga_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp #* data$tChng

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD.prev * data$aD_p.prev * data$aD_c.prev * data$tD.prev * data$typeD.prev * data$sp.prev
  data$tDI_r <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp



  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, .keep_all = TRUE) %>%
      # Compute estimates at plot level
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(nLive = length(which(TPAGROW_UNADJ > 0 | TPAMORT_UNADJ > 0 | TPAREMV_UNADJ > 0)) - length(which(TPARECR_UNADJ > 0)),
                RECR_TPA = if_else(nLive >= minLive, sum(TPARECR_UNADJ * tDI_r, na.rm = TRUE), 0),
                MORT_TPA = if_else(nLive >= minLive, sum(TPAMORT_UNADJ * tDI, na.rm = TRUE), 0),
                REMV_TPA = if_else(nLive >= minLive, sum(TPAREMV_UNADJ * tDI, na.rm = TRUE), 0),
                PREV_TPA = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI, na.rm = TRUE) - RECR_TPA + MORT_TPA + REMV_TPA, 0),
                REMPER = first(REMPER),
                MORT_RATE = if_else(nLive >= minLive, 1 - ((1 - (MORT_TPA / PREV_TPA))^(1/REMPER)), 0),
                HARV_RATE = if_else(nLive >= minLive, 1 - ((1 - (REMV_TPA / PREV_TPA))^(1/REMPER)), 0),
                RECR_RATE = if_else(nLive >= minLive, ((1 + (RECR_TPA / PREV_TPA))^(1/REMPER)) - 1, 0),
                LAMBDA = if_else(nLive >= minLive, RECR_RATE - MORT_RATE - HARV_RATE, 0),
                BAA = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI * basalArea(DIA_END), na.rm = TRUE), 0),
                ## For previous BAA, we use DIA_BEGIN to eliminate recruitment
                ## Have to add HARV and MORT back though becuase not included in GROW
                HARV_BAA = if_else(nLive >= minLive, sum(TPAREMV_UNADJ * tDI * basalArea(DIA_BEGIN), na.rm = TRUE), 0),
                MORT_BAA = if_else(nLive >= minLive, sum(TPAMORT_UNADJ * tDI * basalArea(DIA_BEGIN), na.rm = TRUE), 0),
                PREV_BAA = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI * basalArea(DIA_BEGIN), na.rm = TRUE)+ HARV_BAA + MORT_BAA, 0),
                BAA_RATE =  if_else(nLive >= minLive, ((1 + ((BAA - PREV_BAA) / PREV_BAA))^(1/REMPER)) - 1, 0),
                x = projectPnts(LAMBDA, BAA_RATE, 1, 0)$x,
                y = projectPnts(LAMBDA, BAA_RATE, 1, 0)$y,
                M = sqrt(x^2 + y^2),
                SUST_INDEX = if_else(x < 0, -M, M),
                nStems = length(which(tDI == 1))) %>%
      ungroup() %>%
      select(grpBy, SUST_INDEX, LAMBDA, BAA_RATE, MORT_RATE, HARV_RATE, RECR_RATE, nLive, nStems)

    a = NULL

  } else {
    ### Plot-level estimates -- growth accounting
    a_ga <- data %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, SUBP, CONDID, .keep_all = TRUE) %>%
      group_by(PLT_CN, PROP_BASIS, .dots = aGrpBy) %>%
      summarize(fa_ga = sum(SUBPTYP_PROP_CHNG * aDI_ga, na.rm = TRUE),
                plotIn_ga = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0))
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
      summarize(nLive = length(which(TPAGROW_UNADJ > 0 | TPAMORT_UNADJ > 0 | TPAREMV_UNADJ > 0)) - length(which(TPARECR_UNADJ > 0)),
                REMPER = first(REMPER),
                ############ No growth accounting
                RECR_TPA = if_else(nLive >= minLive, sum(TPARECR_UNADJ * tDI_r, na.rm = TRUE), 0),
                MORT_TPA = if_else(nLive >= minLive, sum(TPAMORT_UNADJ * tDI, na.rm = TRUE), 0),
                REMV_TPA = if_else(nLive >= minLive, sum(TPAREMV_UNADJ * tDI, na.rm = TRUE), 0),
                PREV_TPA = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI, na.rm = TRUE) - RECR_TPA + MORT_TPA + REMV_TPA, 0),
                MORT_RATE = if_else(nLive >= minLive, 1 - ((1 - (MORT_TPA / PREV_TPA))^(1/REMPER)), 0),
                HARV_RATE = if_else(nLive >= minLive, 1 - ((1 - (REMV_TPA / PREV_TPA))^(1/REMPER)), 0),
                RECR_RATE = if_else(nLive >= minLive, ((1 + (RECR_TPA / PREV_TPA))^(1/REMPER)) - 1, 0),
                LAMBDA = if_else(nLive >= minLive, RECR_RATE - MORT_RATE - HARV_RATE, 0),
                BAA = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI * basalArea(DIA_END), na.rm = TRUE), 0),
                ## For previous BAA, we use DIA_BEGIN to eliminate recruitment
                ## Have to add HARV and MORT back though becuase not included in GROW
                HARV_BAA = if_else(nLive >= minLive, sum(TPAREMV_UNADJ * tDI * basalArea(DIA_BEGIN), na.rm = TRUE), 0),
                MORT_BAA = if_else(nLive >= minLive, sum(TPAMORT_UNADJ * tDI * basalArea(DIA_BEGIN), na.rm = TRUE), 0),
                PREV_BAA = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI * basalArea(DIA_BEGIN), na.rm = TRUE)+ HARV_BAA + MORT_BAA, 0),
                BAA_RATE =  if_else(nLive >= minLive, ((1 + ((BAA - PREV_BAA) / PREV_BAA))^(1/REMPER)) - 1, 0),
                x = projectPnts(LAMBDA, BAA_RATE, 1, 0)$x,
                y = projectPnts(LAMBDA, BAA_RATE, 1, 0)$y,
                M = sqrt(x^2 + y^2),
                SUST_INDEX = if_else(x < 0, -M, M),
                ############# Growth accouting
                RECR_TPA_ga = if_else(nLive >= minLive, sum(TPARECR_UNADJ * tDI_ga_r, na.rm = TRUE), 0),
                MORT_TPA_ga = if_else(nLive >= minLive, sum(TPAMORT_UNADJ * tDI_ga, na.rm = TRUE), 0),
                REMV_TPA_ga = if_else(nLive >= minLive, sum(TPAREMV_UNADJ * tDI_ga, na.rm = TRUE), 0),
                PREV_TPA_ga = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI_ga, na.rm = TRUE) - RECR_TPA + MORT_TPA + REMV_TPA, 0),
                MORT_RATE_ga = if_else(nLive >= minLive, 1 - ((1 - (MORT_TPA_ga / PREV_TPA_ga))^(1/REMPER)), 0),
                HARV_RATE_ga = if_else(nLive >= minLive, 1 - ((1 - (REMV_TPA_ga / PREV_TPA_ga))^(1/REMPER)), 0),
                RECR_RATE_ga = if_else(nLive >= minLive, ((1 + (RECR_TPA_ga / PREV_TPA_ga))^(1/REMPER)) - 1, 0),
                LAMBDA_ga = if_else(nLive >= minLive, RECR_RATE_ga - MORT_RATE_ga - HARV_RATE_ga, 0),
                BAA_ga = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI_ga * basalArea(DIA_END), na.rm = TRUE), 0),
                ## For previous BAA, we use DIA_BEGIN to eliminate recruitment
                ## Have to add HARV and MORT back though becuase not included in GROW
                HARV_BAA_ga = if_else(nLive >= minLive, sum(TPAREMV_UNADJ * tDI_ga * basalArea(DIA_BEGIN), na.rm = TRUE), 0),
                MORT_BAA_ga = if_else(nLive >= minLive, sum(TPAMORT_UNADJ * tDI_ga * basalArea(DIA_BEGIN), na.rm = TRUE), 0),
                PREV_BAA_ga = if_else(nLive >= minLive, sum(TPAGROW_UNADJ * tDI_ga * basalArea(DIA_BEGIN), na.rm = TRUE)+ HARV_BAA + MORT_BAA, 0),
                BAA_RATE_ga =  if_else(nLive >= minLive, ((1 + ((BAA_ga - PREV_BAA_ga) / PREV_BAA_ga))^(1/REMPER)) - 1, 0),
                x_ga = projectPnts(LAMBDA_ga, BAA_RATE_ga, 1, 0)$x,
                y_ga = projectPnts(LAMBDA_ga, BAA_RATE_ga, 1, 0)$y,
                M_ga = sqrt(x_ga^2 + y_ga^2),
                SUST_INDEX_ga = if_else(x_ga < 0, -M_ga, M_ga),
                plotIn = if_else(sum(tDI) > 0 | sum(tDI_ga) > 0, 1, 0))
  }

  pltOut <- list(a = a, t = t)
  return(pltOut)

}

diMortHelper2_old1 <- function(x, popState, a, t, grpBy, aGrpBy, method){

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
      sPlot = case_when(
        GROWTH_ACCT == 'Y' ~ SUST_INDEX_ga * tAdj,
        TRUE ~ SUST_INDEX * tAdj),
      lPlot = case_when(
        GROWTH_ACCT == 'Y' ~ LAMBDA_ga * tAdj,
        TRUE ~ LAMBDA * tAdj),
      mPlot = case_when(
        GROWTH_ACCT == 'Y' ~ MORT_RATE_ga * tAdj,
        TRUE ~ MORT_RATE * tAdj),
      hPlot = case_when(
        GROWTH_ACCT == 'Y' ~ HARV_RATE_ga * tAdj,
        TRUE ~ HARV_RATE * tAdj),
      rPlot = case_when(
        GROWTH_ACCT == 'Y' ~ RECR_RATE_ga * tAdj,
        TRUE ~ RECR_RATE * tAdj),
      bPlot = case_when(
        GROWTH_ACCT == 'Y' ~ BAA_RATE_ga * tAdj,
        TRUE ~ BAA_RATE * tAdj),

      plotIn = case_when(
        GROWTH_ACCT == 'Y' ~ plotIn,
        TRUE ~ plotIn)) %>%
    ## Extra step for variance issues
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(sPlot = sum(sPlot, na.rm = TRUE),
              lPlot = sum(lPlot, na.rm = TRUE),
              mPlot = sum(mPlot, na.rm = TRUE),
              hPlot = sum(hPlot, na.rm = TRUE),
              rPlot = sum(rPlot, na.rm = TRUE),
              bPlot = sum(bPlot, na.rm = TRUE),
              fa = first(fa),
              plotIn_t = ifelse(sum(plotIn >  0, na.rm = TRUE), 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    ## Joining area data so we can compute ratio variances
    left_join(select(aStrat, aStrat, av, ESTN_UNIT_CN, STRATUM_CN, ESTN_METHOD, aGrpBy), by = c('ESTN_UNIT_CN', 'ESTN_METHOD', 'STRATUM_CN', aGrpBy)) %>%
    ## Strata level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              sStrat = mean(sPlot * r_t, na.rm = TRUE),
              lStrat = mean(lPlot * r_t, na.rm = TRUE),
              mStrat = mean(mPlot * r_t, na.rm = TRUE),
              hStrat = mean(hPlot * r_t, na.rm = TRUE),
              rStrat = mean(rPlot * r_t, na.rm = TRUE),
              bStrat = mean(bPlot * r_t, na.rm = TRUE),
              aStrat = first(aStrat),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              sv = stratVar(ESTN_METHOD, sPlot, sStrat, ndif, a, nh),
              lv = stratVar(ESTN_METHOD, lPlot, lStrat, ndif, a, nh),
              mv = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh),
              hv = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh),
              rv = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh),
              bv = stratVar(ESTN_METHOD, bPlot, bStrat, ndif, a, nh),

              # Strata level covariances
              cvStrat_s = stratVar(ESTN_METHOD, sPlot, sStrat, ndif, a, nh, fa, aStrat),
              cvStrat_l = stratVar(ESTN_METHOD, lPlot, lStrat, ndif, a, nh, fa, aStrat),
              cvStrat_m = stratVar(ESTN_METHOD, mPlot, mStrat, ndif, a, nh, fa, aStrat),
              cvStrat_h = stratVar(ESTN_METHOD, hPlot, hStrat, ndif, a, nh, fa, aStrat),
              cvStrat_r = stratVar(ESTN_METHOD, rPlot, rStrat, ndif, a, nh, fa, aStrat),
              cvStrat_b = stratVar(ESTN_METHOD, bPlot, bStrat, ndif, a, nh, fa, aStrat)
    ) %>%

    ## Estimation unit
    left_join(select(aEst, ESTN_UNIT_CN, aEst, aVar, aGrpBy), by = c('ESTN_UNIT_CN', aGrpBy)) %>%
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(sEst = unitMean(ESTN_METHOD, a, nh, w, sStrat),
              lEst = unitMean(ESTN_METHOD, a, nh, w, lStrat),
              mEst = unitMean(ESTN_METHOD, a, nh, w, mStrat),
              hEst = unitMean(ESTN_METHOD, a, nh, w, hStrat),
              rEst = unitMean(ESTN_METHOD, a, nh, w, rStrat),
              bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
              #aEst = first(aEst),
              # Estimation of unit variance
              sVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, sv, sStrat, sEst),
              lVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, lv, lStrat, lEst),
              mVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mv, mStrat, mEst),
              hVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, hv, hStrat, hEst),
              rVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, rv, rStrat, rEst),
              bVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bv, bStrat, bEst),
              ## Covariances
              cvEst_s = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_s, sStrat, sEst, aStrat, aEst),
              cvEst_l = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_l, lStrat, lEst, aStrat, aEst),
              cvEst_m = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_m, mStrat, mEst, aStrat, aEst),
              cvEst_h = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_h, hStrat, hEst, aStrat, aEst),
              cvEst_r = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_r, rStrat, rEst, aStrat, aEst),
              cvEst_b = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_b, bStrat, bEst, aStrat, aEst),
              plotIn_t = sum(plotIn_t, na.rm = TRUE))

  out <- list(tEst = tEst, aEst = aEst)

  return(out)
}

diMortHelper1_old <- function(x, plts, db, grpBy, byPlot, minLive){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  ## Carrying out filter across all tables
  #db <- clipFIA(db, mostRecent = FALSE)


  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% grpBy]
  grpC <- names(db$COND)[names(db$COND) %in% grpBy & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% grpBy & names(db$TREE) %in% c(grpP, grpC) == FALSE]

  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- select(db$PLOT, c('PLT_CN', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD', grpP, 'aD_p', 'sp', PREV_PLT_CN, REMPER)) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    left_join(select(db$TREE, c('PLT_CN', 'CONDID', 'DIA', 'SUBP', 'TREE', grpT, 'tD', 'typeD', 'PREV_TRE_CN', 'PREVDIA', 'STATUSCD', 'TPA_UNADJ')), by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$TREE, CN, STATUSCD, TPA_UNADJ), by = c('PREV_TRE_CN' = 'CN'), suffix = c('', '.prev')) %>%
    ## FOR A TREE TO BE CONSIDERED INGROWTH
    mutate(TPARECR_UNADJ = case_when(!is.na(PREV_PLT_CN) &
                                       is.na(PREV_TRE_CN) &
                                       is.na(PREVDIA) &
                                       STATUSCD == 1 &
                                       DIA >= 5 ~ TPA_UNADJ,
                                     TRUE ~ 0),
           TPAMORT_UNADJ = case_when(!is.na(PREV_PLT_CN) &
                                       !is.na(PREV_TRE_CN) &
                                       !is.na(PREVDIA) &
                                       STATUSCD == 2 &
                                       STATUSCD.prev == 1 &
                                       DIA >= 5 ~ TPA_UNADJ,
                                     TRUE ~ 0),
           BAA = case_when(!is.na(PREV_PLT_CN) &
                             !is.na(PREV_TRE_CN) &
                             !is.na(PREVDIA) &
                             STATUSCD == 1 &
                             STATUSCD.prev == 1 &
                             DIA >= 5 ~ basalArea(DIA) * TPA_UNADJ,
                           TRUE ~ 0),
           PREV_BAA = case_when(!is.na(PREV_PLT_CN) &
                                  !is.na(PREV_TRE_CN) &
                                  !is.na(PREVDIA) &
                                  STATUSCD.prev == 1 &
                                  DIA >= 5 ~ basalArea(PREVDIA) * TPA_UNADJ.prev,
                                TRUE ~ 0)) %>%
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
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR')) #%>%
  ## Trying this, not sure
  #filter(!is.na(PREV_TRE_CN))


  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD_p * data$aD_c * data$sp
  data$tDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$typeD * data$sp
  data$pDI <- data$landD * data$aD_p * data$aD_c * data$tD * data$sp


  if (byPlot){
    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    # t <- data %>%
    #   mutate(YEAR = MEASYEAR) %>%
    #   distinct(PLT_CN, SUBP, CONDID, TREE, .keep_all = TRUE) %>%
    #   group_by(.dots = grpBy, PLT_CN) %>%
    #   summarize(H = divIndex(grp, state  * tDI, index = 'H'),
    #             S = divIndex(grp, state * tDI, index = 'S'),
    #             Eh = divIndex(grp, state * tDI, index = 'Eh'),
    #             nStems = length(which(tDI == 1)))

    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, CONDID, TREE, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(nLive = length(which(tDI == 1 & STATUSCD.prev == 1)),
                MORT_TPA = if_else(nLive >= minLive, sum(TPAMORT_UNADJ  * tDI, na.rm = TRUE), 0),
                RECR_TPA = if_else(nLive >= minLive, sum(TPARECR_UNADJ  * tDI, na.rm = TRUE), 0),
                PREV_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ.prev[STATUSCD.prev == 1]  * tDI, na.rm = TRUE), 0),
                REMPER = first(REMPER),
                Q1 = if_else(nLive >= minLive, 1 - ((1 - (MORT_TPA / PREV_TPA))^(1/REMPER)), 0),
                P1 = if_else(nLive >= minLive, ((1 + (RECR_TPA / PREV_TPA))^(1/REMPER)) - 1, 0),
                S1 = if_else(nLive >= minLive, P1 - Q1, 0),
                BAA = if_else(nLive >= minLive, sum(BAA  * tDI, na.rm = TRUE), 0),
                PREV_BAA = if_else(nLive >= minLive, sum(PREV_BAA  * tDI, na.rm = TRUE), 0),
                BA1 =  if_else(nLive >= minLive, ((1 + ((BAA - PREV_BAA) / PREV_BAA))^(1/REMPER)) - 1, 0),
                x = projectPnts(S1, BA1, 1, 0)$x,
                y = projectPnts(S1, BA1, 1, 0)$y,
                M = sqrt(x^2 + y^2),
                M = if_else(x < 0, -M, M),
                nStems = length(which(tDI == 1))) %>%
      ungroup() %>%
      select(grpBy, M, Q1, P1, S1, BA1, nLive, nStems)
    a = NULL

  } else {
    # variable is computed at the stand (condition level), and we continue to use the ratio of means estimator to get at
    #  average of the attribute of interest weighted by the area in which it occurs.
    t <- data %>%
      distinct(PLT_CN, CONDID, TREE, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN, PROP_BASIS, CONDID) %>%
      # filter(tDI > 0) %>%
      summarize(condArea = first(CONDPROP_UNADJ),
                nLive = length(which(tDI == 1 & STATUSCD.prev == 1)),
                MORT_TPA = if_else(nLive >= minLive, sum(TPAMORT_UNADJ  * tDI, na.rm = TRUE), 0),
                RECR_TPA = if_else(nLive >= minLive, sum(TPARECR_UNADJ  * tDI, na.rm = TRUE), 0),
                PREV_TPA = if_else(nLive >= minLive, sum(TPA_UNADJ.prev[STATUSCD.prev == 1]  * tDI, na.rm = TRUE), 0),
                REMPER = first(REMPER),
                Q1 = if_else(nLive >= minLive, 1 - ((1 - (MORT_TPA / PREV_TPA))^(1/REMPER)), 0) * condArea,
                P1 = if_else(nLive >= minLive, ((1 + (RECR_TPA / PREV_TPA))^(1/REMPER)) - 1, 0) * condArea,
                S1 = if_else(nLive >= minLive, P1 - Q1, 0),
                BAA = if_else(nLive >= minLive, sum(BAA  * tDI, na.rm = TRUE), 0),
                PREV_BAA = if_else(nLive >= minLive, sum(PREV_BAA  * tDI, na.rm = TRUE), 0),
                BA1 =  if_else(nLive >= minLive, ((1 + ((BAA - PREV_BAA) / PREV_BAA))^(1/REMPER)) - 1, 0) * condArea,
                x = projectPnts(S1, BA1, 1, 0)$x,
                y = projectPnts(S1, BA1, 1, 0)$y,
                M = sqrt(x^2 + y^2),
                M = if_else(x < 0, -M, M),
                plotIn = ifelse(sum(tDI >  0, na.rm = TRUE), 1,0),
                aDI = ifelse(sum(tDI > 0, na.rm = TRUE), 1, 0)) %>%
      group_by(.dots = grpBy, PROP_BASIS, PLT_CN) %>%
      summarize(Q1 = sum(Q1, na.rm = TRUE),
                P1 = sum(P1, na.rm = TRUE),
                S1 = sum(S1, na.rm = TRUE),
                BA1 = sum(BA1, na.rm = TRUE),
                M = sum(M, na.rm = TRUE),
                fa = sum(condArea * aDI, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE))
  }

  pltOut <- list(t = t)
  return(pltOut)
}



diMortHelper2_old <- function(x, popState, t, grpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    #aGrpBy <- c(aGrpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR

  }

  ## Strata level estimates
  tEst <- t %>%
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
      fa = fa * aAdj,
      Q1 = Q1 * aAdj,
      P1 = P1 * aAdj,
      S1 = S1 * aAdj,
      BA1 = BA1 * aAdj,
      M = M * aAdj) %>%
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(a_t = length(unique(PLT_CN)) / first(P2POINTCNT),
              aStrat = mean(fa * a_t, na.rm = TRUE),
              qStrat = mean(Q1 * a_t, na.rm = TRUE),
              pStrat = mean(P1 * a_t, na.rm = TRUE),
              sStrat = mean(S1 * a_t, na.rm = TRUE),
              bStrat = mean(BA1 * a_t, na.rm = TRUE),
              mStrat = mean(M * a_t, na.rm = TRUE),
              plotIn_AREA = sum(plotIn, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(P2POINTCNT),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU),
              p2eu = first(p2eu),
              ndif = nh - n,
              ## Strata level variances
              av = stratVar(ESTN_METHOD, fa, aStrat, ndif, a, nh),
              qv = stratVar(ESTN_METHOD, Q1, qStrat, ndif, a, nh),
              pv = stratVar(ESTN_METHOD, P1, pStrat, ndif, a, nh),
              sv = stratVar(ESTN_METHOD, S1, sStrat, ndif, a, nh),
              bv = stratVar(ESTN_METHOD, BA1, bStrat, ndif, a, nh),
              mv = stratVar(ESTN_METHOD, M, mStrat, ndif, a, nh),
              # Strata level covariances
              cvStrat_q = stratVar(ESTN_METHOD, Q1, qStrat, ndif, a, nh, fa, aStrat),
              cvStrat_p = stratVar(ESTN_METHOD, P1, pStrat, ndif, a, nh, fa, aStrat),
              cvStrat_s = stratVar(ESTN_METHOD, S1, sStrat, ndif, a, nh, fa, aStrat),
              cvStrat_b = stratVar(ESTN_METHOD, BA1, bStrat, ndif, a, nh, fa, aStrat),
              cvStrat_m = stratVar(ESTN_METHOD, M, mStrat, ndif, a, nh, fa, aStrat)) %>%
    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh,  w, aStrat),
              qEst = unitMean(ESTN_METHOD, a, nh,  w, qStrat),
              pEst = unitMean(ESTN_METHOD, a, nh,  w, pStrat),
              sEst = unitMean(ESTN_METHOD, a, nh,  w, sStrat),
              bEst = unitMean(ESTN_METHOD, a, nh,  w, bStrat),
              mEst = unitMean(ESTN_METHOD, a, nh,  w, mStrat),
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, av, aStrat, aEst),
              qVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, qv, qStrat, qEst),
              pVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pv, pStrat, pEst),
              sVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, sv, sStrat, sEst),
              bVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, bv, bStrat, bEst),
              mVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, mv, mStrat, mEst),
              cvEst_q = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_q, qStrat, qEst, aStrat, aEst),
              cvEst_p = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_p, pStrat, pEst, aStrat, aEst),
              cvEst_s = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_s, sStrat, sEst, aStrat, aEst),
              cvEst_b = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_b, bStrat, bEst, aStrat, aEst),
              cvEst_m = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_m, mStrat, mEst, aStrat, aEst),
              plotIn_AREA = sum(plotIn_AREA, na.rm = TRUE))

  out <- list(tEst = tEst)

  return(out)
}






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
                H_a = sum(hPlot, na.rm = TRUE) / AREA_TOTAL,
                Eh_a = sum(EhPlot, na.rm = TRUE) / AREA_TOTAL,
                S_a = sum(sPlot, na.rm = TRUE) / AREA_TOTAL,
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
  #gc()

  #Return a dataframe
  d

}
