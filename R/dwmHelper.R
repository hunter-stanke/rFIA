dwmHelper1 <- function(x, plts, db, grpBy, byPlot){

  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]


  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$COND_DWM_CALC, by = c('PLT_CN', 'CONDID'))

  ## Comprehensive indicator function
  data$aDI <- data$landD * data$aD * data$sp



  if (byPlot){

    grpBy <- c('YEAR', grpBy, 'PLOT_STATUS_CD')
    grpSyms <- syms(grpBy)

    t <- data %>%
      mutate(YEAR = INVYR) %>%
      distinct(PLT_CN, CONDID, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(!!!grpSyms, PLT_CN) %>%
      summarize(VOL_1HR = sum(FWD_SM_VOLCF_ADJ * aDI, na.rm = TRUE),
        VOL_10HR = sum(FWD_MD_VOLCF_ADJ * aDI, na.rm = TRUE),
        VOL_100HR = sum(FWD_LG_VOLCF_ADJ * aDI, na.rm = TRUE),
        VOL_1000HR = sum(CWD_VOLCF_ADJ * aDI, na.rm = TRUE),
        VOL_PILE = sum(PILE_VOLCF_ADJ * aDI, na.rm = TRUE),
        BIO_DUFF = sum(DUFF_BIOMASS* aDI / 2000, na.rm = TRUE),
        BIO_LITTER = sum(LITTER_BIOMASS * aDI / 2000, na.rm = TRUE),
        BIO_1HR = sum(FWD_SM_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
        BIO_10HR = sum(FWD_MD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
        BIO_100HR = sum(FWD_LG_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
        BIO_1000HR = sum(CWD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
        BIO_PILE = sum(PILE_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
        CARB_DUFF = sum(DUFF_CARBON* aDI / 2000, na.rm = TRUE),
        CARB_LITTER = sum(LITTER_CARBON * aDI / 2000, na.rm = TRUE),
        CARB_1HR = sum(FWD_SM_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
        CARB_10HR = sum(FWD_MD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
        CARB_100HR = sum(FWD_LG_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
        CARB_1000HR = sum(CWD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
        CARB_PILE = sum(PILE_CARBON_ADJ * aDI / 2000, na.rm = TRUE)) %>%
      mutate(BIO = BIO_LITTER + BIO_DUFF + BIO_1HR + BIO_10HR + BIO_100HR + BIO_1000HR + BIO_PILE,
             VOL = VOL_1HR + VOL_10HR+VOL_100HR+VOL_1000HR+VOL_PILE,
             CARB = CARB_LITTER + CARB_DUFF + CARB_1HR + CARB_10HR + CARB_100HR + CARB_1000HR + CARB_PILE) %>%
      as.data.frame()

  } else {

    grpSyms <- syms(grpBy)


    # Compute estimates
    t <- data %>%
      distinct(STRATUM_CN, PLT_CN, CONDID, COND_STATUS_CD, .keep_all = TRUE) %>%
      lazy_dt() %>%
      group_by(STRATUM_CN, PROP_BASIS, PLT_CN, !!!grpSyms) %>%
      summarize(vsmPlot = sum(FWD_SM_VOLCF_ADJ * aDI, na.rm = TRUE),
                vmdPlot = sum(FWD_MD_VOLCF_ADJ * aDI, na.rm = TRUE),
                vlgPlot = sum(FWD_LG_VOLCF_ADJ * aDI, na.rm = TRUE),
                vcPlot = sum(CWD_VOLCF_ADJ * aDI, na.rm = TRUE),
                vpPlot = sum(PILE_VOLCF_ADJ * aDI, na.rm = TRUE),
                bdPlot = sum(DUFF_BIOMASS* aDI / 2000, na.rm = TRUE),
                blPlot = sum(LITTER_BIOMASS * aDI / 2000, na.rm = TRUE),
                bsmPlot = sum(FWD_SM_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bmdPlot = sum(FWD_MD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                blgPlot = sum(FWD_LG_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bcPlot = sum(CWD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bpPlot = sum(PILE_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                cdPlot = sum(DUFF_CARBON* aDI / 2000, na.rm = TRUE),
                clPlot = sum(LITTER_CARBON * aDI / 2000, na.rm = TRUE),
                csmPlot = sum(FWD_SM_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                cmdPlot = sum(FWD_MD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                clgPlot = sum(FWD_LG_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                ccPlot = sum(CWD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                cpPlot = sum(PILE_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                fa = sum(CONDPROP_UNADJ * aDI, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      mutate(vPlot = vsmPlot + vmdPlot+ vlgPlot+ vcPlot+ vpPlot,
             bPlot = bdPlot+ blPlot+ bsmPlot+ bmdPlot+ blgPlot+ bcPlot+ bpPlot,
             cPlot = cdPlot+ clPlot+ csmPlot+ cmdPlot + clgPlot + ccPlot+ cpPlot) %>%
      as.data.frame()
  }

  pltOut <- list(t = t)
  return(pltOut)

}



dwmHelper2 <- function(x, popState, t, grpBy, method){

  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$P1POINTCNT <- popState[[x]]$P1POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }

  grpSyms = syms(grpBy)

  ## Strata level estimates
  tEst <- t %>%
    lazy_dt() %>%
    ## Rejoin with population tables
    inner_join(select(popState[[x]], -c(STATECD)), by = c('STRATUM_CN', 'PLT_CN')) %>%
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
              a = dplyr::first(AREA_USED),
              w = dplyr::first(P1POINTCNT) / dplyr::first(P1PNTCNT_EU),
              p2eu = dplyr::first(p2eu),

              ## dtplyr is fast, but requires a few extra steps, so we'll finish
              ## means and variances in subseqent mutate step

              aStrat = sum(fa , na.rm = TRUE),
              vsmStrat = sum(vsmPlot , na.rm = TRUE),
              vmdStrat = sum(vmdPlot , na.rm = TRUE),
              vlgStrat = sum(vlgPlot , na.rm = TRUE),
              vcStrat = sum(vcPlot , na.rm = TRUE),
              vpStrat = sum(vpPlot , na.rm = TRUE),
              vStrat = sum(vPlot , na.rm = TRUE),
              bdStrat = sum(bdPlot , na.rm = TRUE),
              blStrat = sum(blPlot , na.rm = TRUE),
              bsmStrat = sum(bsmPlot , na.rm = TRUE),
              bmdStrat = sum(bmdPlot , na.rm = TRUE),
              blgStrat = sum(blgPlot , na.rm = TRUE),
              bcStrat = sum(bcPlot , na.rm = TRUE),
              bpStrat = sum(bpPlot , na.rm = TRUE),
              bStrat = sum(bPlot , na.rm = TRUE),
              cdStrat = sum(cdPlot , na.rm = TRUE),
              clStrat = sum(clPlot , na.rm = TRUE),
              csmStrat = sum(csmPlot , na.rm = TRUE),
              cmdStrat = sum(cmdPlot , na.rm = TRUE),
              clgStrat = sum(clgPlot , na.rm = TRUE),
              ccStrat = sum(ccPlot , na.rm = TRUE),
              cpStrat = sum(cpPlot , na.rm = TRUE),
              cStrat = sum(cPlot , na.rm = TRUE),
              plotIn = sum(plotIn, na.rm = TRUE),

              # ## Strata level variances
              av = sum(fa^2, na.rm = TRUE),
              vsmVar = sum(vsmPlot^2, na.rm = TRUE),
              vmdVar = sum(vmdPlot^2, na.rm = TRUE),
              vlgVar = sum(vlgPlot^2, na.rm = TRUE),
              vcVar = sum(vcPlot^2, na.rm = TRUE),
              vpVar = sum(vpPlot^2, na.rm = TRUE),
              vVar = sum(vPlot^2, na.rm = TRUE),
              bdVar = sum(bdPlot^2, na.rm = TRUE),
              blVar = sum(blPlot^2, na.rm = TRUE),
              bsmVar = sum(bsmPlot^2, na.rm = TRUE),
              bmdVar = sum(bmdPlot^2, na.rm = TRUE),
              blgVar = sum(blgPlot^2, na.rm = TRUE),
              bcVar = sum(bcPlot^2, na.rm = TRUE),
              bpVar = sum(bpPlot^2, na.rm = TRUE),
              bVar = sum(bPlot^2, na.rm = TRUE),
              cdVar = sum(cdPlot^2, na.rm = TRUE),
              clVar = sum(clPlot^2, na.rm = TRUE),
              csmVar = sum(csmPlot^2, na.rm = TRUE),
              cmdVar = sum(cmdPlot^2, na.rm = TRUE),
              clgVar = sum(clgPlot^2, na.rm = TRUE),
              ccVar = sum(ccPlot^2, na.rm = TRUE),
              cpVar = sum(cpPlot^2, na.rm = TRUE),
              cVar = sum(cPlot^2, na.rm = TRUE),
              # ## Strata level co-variances
              cvStrat_vsm = sum(vsmPlot*fa, na.rm = TRUE),
              cvStrat_vmd = sum(vmdPlot*fa, na.rm = TRUE),
              cvStrat_vlg = sum(vlgPlot*fa, na.rm = TRUE),
              cvStrat_vc = sum(vcPlot*fa, na.rm = TRUE),
              cvStrat_vp = sum(vpPlot*fa, na.rm = TRUE),
              cvStrat_v = sum(vPlot*fa, na.rm = TRUE),
              cvStrat_bd = sum(bdPlot*fa, na.rm = TRUE),
              cvStrat_bl = sum(blPlot*fa, na.rm = TRUE),
              cvStrat_bsm = sum(bsmPlot*fa, na.rm = TRUE),
              cvStrat_bmd = sum(bmdPlot*fa, na.rm = TRUE),
              cvStrat_blg = sum(blgPlot*fa, na.rm = TRUE),
              cvStrat_bc = sum(bcPlot*fa, na.rm = TRUE),
              cvStrat_bp = sum(bpPlot*fa, na.rm = TRUE),
              cvStrat_b = sum(bPlot*fa, na.rm = TRUE),
              cvStrat_cd = sum(cdPlot*fa, na.rm = TRUE),
              cvStrat_cl = sum(clPlot*fa, na.rm = TRUE),
              cvStrat_csm = sum(csmPlot*fa, na.rm = TRUE),
              cvStrat_cmd = sum(cmdPlot*fa, na.rm = TRUE),
              cvStrat_clg = sum(clgPlot*fa, na.rm = TRUE),
              cvStrat_cc = sum(ccPlot*fa, na.rm = TRUE),
              cvStrat_cp = sum(cpPlot*fa, na.rm = TRUE),
              cvStrat_c = sum(cPlot*fa, na.rm = TRUE)) %>%
    mutate(aStrat = aStrat / nh,
           vsmStrat = vsmStrat / nh,
           vmdStrat = vmdStrat / nh,
           vlgStrat = vlgStrat / nh,
           vcStrat = vcStrat / nh,
           vpStrat = vpStrat / nh,
           vStrat = vStrat / nh,
           bdStrat = bdStrat / nh,
           blStrat = blStrat / nh,
           bsmStrat = bsmStrat / nh,
           bmdStrat = bmdStrat / nh,
           blgStrat = blgStrat / nh,
           bcStrat = bcStrat / nh,
           bpStrat = bpStrat / nh,
           bStrat = bStrat / nh,
           cdStrat = cdStrat / nh,
           clStrat = clStrat / nh,
           csmStrat = csmStrat / nh,
           cmdStrat = cmdStrat / nh,
           clgStrat = clgStrat / nh,
           ccStrat = ccStrat / nh,
           cpStrat = cpStrat / nh,
           cStrat = cStrat / nh,

           adj = nh * (nh-1),

           # ## Strata level variances
           av = (av - (nh*aStrat^2)) / adj,
           vsmVar = (vsmVar - (nh*vsmStrat^2)) / adj,
           vmdVar =(vmdVar - (nh*vmdStrat^2)) / adj,
           vlgVar = (vlgVar - (nh*vlgStrat^2)) / adj,
           vcVar = (vcVar - (nh*vcStrat^2)) / adj,
           vpVar = (vpVar - (nh*vpStrat^2)) / adj,
           vVar = (vVar - (nh*vStrat^2)) / adj,
           bdVar = (bdVar - (nh*bdStrat^2)) / adj,
           blVar = (blVar - (nh*blStrat^2)) / adj,
           bsmVar = (bsmVar - (nh*bsmStrat^2)) / adj,
           bmdVar = (bmdVar - (nh*bmdStrat^2)) / adj,
           blgVar = (blgVar - (nh*blgStrat^2)) / adj,
           bcVar = (bcVar - (nh*bcStrat^2)) / adj,
           bpVar = (bpVar - (nh*bpStrat^2)) / adj,
           bVar = (bVar - (nh*bStrat^2)) / adj,
           cdVar = (cdVar - (nh*cdStrat^2)) / adj,
           clVar = (clVar - (nh*clStrat^2)) / adj,
           csmVar = (csmVar - (nh*csmStrat^2)) / adj,
           cmdVar = (cmdVar - (nh*cmdStrat^2)) / adj,
           clgVar = (clgVar - (nh*clgStrat^2)) / adj,
           ccVar = (ccVar - (nh*ccStrat^2)) / adj,
           cpVar = (cpVar - (nh*cpStrat^2)) / adj,
           cVar = (cVar - (nh*cStrat^2)) / adj,

           # ## Strata level co-variances
           cvEst_vsm = (cvStrat_vsm - (nh * vsmStrat * aStrat)) / adj,
           cvEst_vmd = (cvStrat_vmd - (nh * vmdStrat * aStrat)) / adj,
           cvEst_vlg = (cvStrat_vlg - (nh * vlgStrat * aStrat)) / adj,
           cvEst_vc = (cvStrat_vc - (nh * vcStrat * aStrat)) / adj,
           cvEst_vp = (cvStrat_vp - (nh * vpStrat * aStrat)) / adj,
           cvEst_v = (cvStrat_v - (nh * vStrat * aStrat)) / adj,
           cvEst_bd = (cvStrat_bd - (nh * bdStrat * aStrat)) / adj,
           cvEst_bl = (cvStrat_bl - (nh * blStrat * aStrat)) / adj,
           cvEst_bsm = (cvStrat_bsm - (nh * bsmStrat * aStrat)) / adj,
           cvEst_bmd = (cvStrat_bmd - (nh * bmdStrat * aStrat)) / adj,
           cvEst_blg = (cvStrat_blg - (nh * blgStrat * aStrat)) / adj,
           cvEst_bc = (cvStrat_bc - (nh * bcStrat * aStrat)) / adj,
           cvEst_bp = (cvStrat_bp - (nh * bpStrat * aStrat)) / adj,
           cvEst_b = (cvStrat_b - (nh * bStrat * aStrat)) / adj,
           cvEst_cd = (cvStrat_cd - (nh * cdStrat * aStrat)) / adj,
           cvEst_cl = (cvStrat_cl - (nh * clStrat * aStrat)) / adj,
           cvEst_csm = (cvStrat_csm - (nh * csmStrat * aStrat)) / adj,
           cvEst_cmd = (cvStrat_cmd - (nh * cmdStrat * aStrat)) / adj,
           cvEst_clg = (cvStrat_clg - (nh * clgStrat * aStrat)) / adj,
           cvEst_cc = (cvStrat_cc - (nh * ccStrat * aStrat)) / adj,
           cvEst_cp = (cvStrat_cp - (nh * cpStrat * aStrat)) / adj,
           cvEst_c = (cvStrat_c - (nh * cStrat * aStrat)) / adj) %>%
    as.data.frame() %>%
    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
              vsmEst = unitMean(ESTN_METHOD, a, nh, w, vsmStrat),
              vmdEst = unitMean(ESTN_METHOD, a, nh, w, vmdStrat),
              vlgEst = unitMean(ESTN_METHOD, a, nh, w, vlgStrat),
              vcEst = unitMean(ESTN_METHOD, a, nh, w, vcStrat),
              vpEst = unitMean(ESTN_METHOD, a, nh, w, vpStrat),
              vEst = unitMean(ESTN_METHOD, a, nh, w, vStrat),
              bdEst = unitMean(ESTN_METHOD, a, nh, w, bdStrat),
              blEst = unitMean(ESTN_METHOD, a, nh, w, blStrat),
              bsmEst = unitMean(ESTN_METHOD, a, nh, w, bsmStrat),
              bmdEst = unitMean(ESTN_METHOD, a, nh, w, bmdStrat),
              blgEst = unitMean(ESTN_METHOD, a, nh, w, blgStrat),
              bcEst = unitMean(ESTN_METHOD, a, nh, w, bcStrat),
              bpEst = unitMean(ESTN_METHOD, a, nh, w, bpStrat),
              bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
              cdEst = unitMean(ESTN_METHOD, a, nh, w, cdStrat),
              clEst = unitMean(ESTN_METHOD, a, nh, w, clStrat),
              csmEst = unitMean(ESTN_METHOD, a, nh, w, csmStrat),
              cmdEst = unitMean(ESTN_METHOD, a, nh, w, cmdStrat),
              clgEst = unitMean(ESTN_METHOD, a, nh, w, clgStrat),
              ccEst = unitMean(ESTN_METHOD, a, nh, w, ccStrat),
              cpEst = unitMean(ESTN_METHOD, a, nh, w, cpStrat),
              cEst = unitMean(ESTN_METHOD, a, nh, w, cStrat),
              plotIn = sum(plotIn, na.rm = TRUE),
              N = dplyr::first(p2eu),
              A = dplyr::first(a),
              # Estimation of unit variance
              aVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, av, aStrat, aEst),
              vsmVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, vsmVar, vsmStrat, vsmEst),
              vmdVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, vmdVar, vmdStrat, vmdEst),
              vlgVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, vlgVar, vlgStrat, vlgEst),
              vcVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, vcVar, vcStrat, vcEst),
              vpVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, vpVar, vpStrat, vpEst),
              vVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, vVar, vStrat, vEst),
              bdVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bdVar, bdStrat, bdEst),
              blVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, blVar, blStrat, blEst),
              bsmVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bsmVar, bsmStrat, bsmEst),
              bmdVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bmdVar, bmdStrat, bmdEst),
              blgVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, blgVar, blgStrat, blgEst),
              bcVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bcVar, bcStrat, bcEst),
              bpVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bpVar, bpStrat, bpEst),
              bVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, bVar, bStrat, cbEst),
              cdVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cdVar, cdStrat, cdEst),
              clVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, clVar, clStrat, clEst),
              csmVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, csmVar, csmStrat, csmEst),
              cmdVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cmdVar, cmdStrat, cmdEst),
              clgVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, clgVar, clgStrat, clgEst),
              ccVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, ccVar, ccStrat, ccEst),
              cpVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cpVar, cpStrat, cpEst),
              cVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cVar, cStrat, cEst),
              # Estimation of unit covariance
              cvEst_vsm = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_vsm, vsmStrat, vsmEst, aStrat, aEst),
              cvEst_vmd = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_vmd, vmdStrat, vmdEst, aStrat, aEst),
              cvEst_vlg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_vlg, vlgStrat, vlgEst, aStrat, aEst),
              cvEst_vc = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_vc, vcStrat, vcEst, aStrat, aEst),
              cvEst_vp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_vp, vpStrat, vpEst, aStrat, aEst),
              cvEst_v = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_v, vStrat, vEst, aStrat, aEst),
              cvEst_bd = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_bd, bdStrat, bdEst, aStrat, aEst),
              cvEst_bl = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_bl, blStrat, blEst, aStrat, aEst),
              cvEst_bsm = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_bsm, bsmStrat, bsmEst, aStrat, aEst),
              cvEst_bmd = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_bmd, bmdStrat, bmdEst, aStrat, aEst),
              cvEst_blg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_blg, blgStrat, blgEst, aStrat, aEst),
              cvEst_bc = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_bc, bcStrat, bcEst, aStrat, aEst),
              cvEst_bp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_bp, bpStrat, bpEst, aStrat, aEst),
              cvEst_b = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_b, bStrat, cbEst, aStrat, aEst),
              cvEst_cd = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_cd, cdStrat, cdEst, aStrat, aEst),
              cvEst_cl = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_cl, clStrat, clEst, aStrat, aEst),
              cvEst_csm = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_csm, csmStrat, csmEst, aStrat, aEst),
              cvEst_cmd = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_cmd, cmdStrat, cmdEst, aStrat, aEst),
              cvEst_clg = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_clg, clgStrat, clgEst, aStrat, aEst),
              cvEst_cc = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_cc, ccStrat, ccEst, aStrat, aEst),
              cvEst_cp = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_cp, cpStrat, cpEst, aStrat, aEst),
              cvEst_c = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, dplyr::first(p2eu), w, cvEst_c, cStrat, cEst, aStrat, aEst))

  out <- list(tEst = tEst)

  return(out)
}













dwmHelper <- function(x, combos, data, grpBy, totals, SE){
  # Update domain indicator for each each column speficed in grpBy
  ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Area domain indicator for each column in
    aObs <- as.character(combos[[x]][[grpBy[n]]]) == as.character(data[[grpBy[n]]])
    if (length(which(is.na(aObs))) == length(aObs)) aObs <- 1
    ad <- data$aDI * aObs * ad
  }

  if(SE){
    data$aDI <- ad
    data$aDI[is.na(data$aDI)] <- 0
    cwd <- data %>%
      distinct(PLT_CN, CONDID, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(vsmPlot = sum(FWD_SM_VOLCF_ADJ * aDI, na.rm = TRUE),
                vmdPlot = sum(FWD_MD_VOLCF_ADJ * aDI, na.rm = TRUE),
                vlgPlot = sum(FWD_LG_VOLCF_ADJ * aDI, na.rm = TRUE),
                vcPlot = sum(CWD_VOLCF_ADJ * aDI, na.rm = TRUE),
                vpPlot = sum(PILE_VOLCF_ADJ * aDI, na.rm = TRUE),
                vPlot = sum(vsmPlot, vmdPlot, vlgPlot, vcPlot, vpPlot, na.rm = TRUE),
                bdPlot = sum(DUFF_BIOMASS* aDI / 2000, na.rm = TRUE),
                blPlot = sum(LITTER_BIOMASS * aDI / 2000, na.rm = TRUE),
                bsmPlot = sum(FWD_SM_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bmdPlot = sum(FWD_MD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                blgPlot = sum(FWD_LG_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bcPlot = sum(CWD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bpPlot = sum(PILE_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bPlot = sum(bdPlot, blPlot, bsmPlot, bmdPlot, blgPlot, bcPlot, bpPlot, na.rm = TRUE),
                cdPlot = sum(DUFF_CARBON* aDI / 2000, na.rm = TRUE),
                clPlot = sum(LITTER_CARBON * aDI / 2000, na.rm = TRUE),
                csmPlot = sum(FWD_SM_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                cmdPlot = sum(FWD_MD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                clgPlot = sum(FWD_LG_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                ccPlot = sum(CWD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                cpPlot = sum(PILE_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                cPlot = sum(cdPlot, clPlot, csmPlot, cmdPlot, clgPlot, ccPlot, cpPlot, na.rm = TRUE),
                fa = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
                a = dplyr::first(AREA_USED),
                p1EU = dplyr::first(P1PNTCNT_EU),
                p1 = dplyr::first(P1POINTCNT),
                p2 = dplyr::first(P2POINTCNT)) %>%
      # Stratum level
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                vsmStrat = mean(vsmPlot, na.rm = TRUE),
                vmdStrat = mean(vmdPlot, na.rm = TRUE),
                vlgStrat = mean(vlgPlot, na.rm = TRUE),
                vcStrat = mean(vcPlot, na.rm = TRUE),
                vpStrat = mean(vpPlot, na.rm = TRUE),
                vStrat = mean(vPlot, na.rm = TRUE),
                bdStrat = mean(bdPlot, na.rm = TRUE),
                blStrat = mean(blPlot, na.rm = TRUE),
                bsmStrat = mean(bsmPlot, na.rm = TRUE),
                bmdStrat = mean(bmdPlot, na.rm = TRUE),
                blgStrat = mean(blgPlot, na.rm = TRUE),
                bcStrat = mean(bcPlot, na.rm = TRUE),
                bpStrat = mean(bpPlot, na.rm = TRUE),
                bStrat = mean(bPlot, na.rm = TRUE),
                cdStrat = mean(cdPlot, na.rm = TRUE),
                clStrat = mean(clPlot, na.rm = TRUE),
                csmStrat = mean(csmPlot, na.rm = TRUE),
                cmdStrat = mean(cmdPlot, na.rm = TRUE),
                clgStrat = mean(clgPlot, na.rm = TRUE),
                ccStrat = mean(ccPlot, na.rm = TRUE),
                cpStrat = mean(cpPlot, na.rm = TRUE),
                cStrat = mean(cPlot, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = dplyr::first(a),
                w = dplyr::first(p1) / dplyr::first(p1EU), # Stratum weight
                nh = dplyr::first(p2), # Number plots in stratum
                # Strata level variances
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * dplyr::first(a) / nh),
                            (sum(fa^2) - sum(nh * aStrat^2)) / (nh * (nh-1))),
                vsmV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(vsmPlot * dplyr::first(a) / nh),
                              (sum(vsmPlot^2) - sum(nh * vsmStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vmdV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(vmdPlot * dplyr::first(a) / nh),
                              (sum(vmdPlot^2) - sum(nh * vmdStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vlgV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(vlgPlot * dplyr::first(a) / nh),
                              (sum(vlgPlot^2) - sum(nh * vlgStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vcV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(vcPlot * dplyr::first(a) / nh),
                             (sum(vcPlot^2) - sum(nh * vcStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vpV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(vpPlot * dplyr::first(a) / nh),
                             (sum(vpPlot^2) - sum(nh * vpStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vV = ifelse(first(ESTN_METHOD == 'simple'),
                            var(vPlot * dplyr::first(a) / nh),
                            (sum(vPlot^2) - sum(nh * vStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bdV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(bdPlot * dplyr::first(a) / nh),
                              (sum(bdPlot^2) - sum(nh * bdStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                blV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(blPlot * dplyr::first(a) / nh),
                              (sum(blPlot^2) - sum(nh * blStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bsmV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(bsmPlot * dplyr::first(a) / nh),
                              (sum(bsmPlot^2) - sum(nh * bsmStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bmdV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(bmdPlot * dplyr::first(a) / nh),
                              (sum(bmdPlot^2) - sum(nh * bmdStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                blgV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(blgPlot * dplyr::first(a) / nh),
                              (sum(blgPlot^2) - sum(nh * blgStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bcV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(bcPlot * dplyr::first(a) / nh),
                             (sum(bcPlot^2) - sum(nh * bcStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bpV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(bpPlot * dplyr::first(a) / nh),
                             (sum(bpPlot^2) - sum(nh * bpStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bV = ifelse(first(ESTN_METHOD == 'simple'),
                            var(bPlot * dplyr::first(a) / nh),
                            (sum(bPlot^2) - sum(nh * bStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cdV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(cdPlot * dplyr::first(a) / nh),
                             (sum(cdPlot^2) - sum(nh * cdStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                clV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(clPlot * dplyr::first(a) / nh),
                             (sum(clPlot^2) - sum(nh * clStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                csmV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(csmPlot * dplyr::first(a) / nh),
                              (sum(csmPlot^2) - sum(nh * csmStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cmdV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(cmdPlot * dplyr::first(a) / nh),
                              (sum(cmdPlot^2) - sum(nh * cmdStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                clgV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(clgPlot * dplyr::first(a) / nh),
                              (sum(clgPlot^2) - sum(nh * clgStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                ccV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(ccPlot * dplyr::first(a) / nh),
                             (sum(ccPlot^2) - sum(nh * ccStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cpV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(cpPlot * dplyr::first(a) / nh),
                             (sum(cpPlot^2) - sum(nh * cpStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cV = ifelse(first(ESTN_METHOD == 'simple'),
                            var(cPlot * dplyr::first(a) / nh),
                            (sum(cPlot^2) - sum(nh * cStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                # Strata level covariances
                vsmCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,vsmPlot),
                               (sum(fa*vsmPlot) - sum(nh * aStrat *vsmStrat)) / (nh * (nh-1))), # Stratified and double cases
                vmdCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,vmdPlot),
                               (sum(fa*vmdPlot) - sum(nh * aStrat *vmdStrat)) / (nh * (nh-1))), # Stratified and double cases
                vlgCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,vlgPlot),
                               (sum(fa*vlgPlot) - sum(nh * aStrat *vlgStrat)) / (nh * (nh-1))), # Stratified and double cases
                vcCV = ifelse(first(ESTN_METHOD == 'simple'),
                              cov(fa,vcPlot),
                              (sum(fa*vcPlot) - sum(nh * aStrat *vcStrat)) / (nh * (nh-1))), # Stratified and double cases
                vpCV = ifelse(first(ESTN_METHOD == 'simple'),
                              cov(fa,vpPlot),
                              (sum(fa*vpPlot) - sum(nh * aStrat *vpStrat)) / (nh * (nh-1))), # Stratified and double cases
                vCV = ifelse(first(ESTN_METHOD == 'simple'),
                             cov(fa,vPlot),
                             (sum(fa*vPlot) - sum(nh * aStrat *vStrat)) / (nh * (nh-1))), # Stratified and double cases
                bdCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,bdPlot),
                               (sum(fa*bdPlot) - sum(nh * aStrat *bdStrat)) / (nh * (nh-1))), # Stratified and double cases
                blCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,blPlot),
                               (sum(fa*blPlot) - sum(nh * aStrat *blStrat)) / (nh * (nh-1))), # Stratified and double cases
                bsmCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,bsmPlot),
                               (sum(fa*bsmPlot) - sum(nh * aStrat *bsmStrat)) / (nh * (nh-1))), # Stratified and double cases
                bmdCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,bmdPlot),
                               (sum(fa*bmdPlot) - sum(nh * aStrat *bmdStrat)) / (nh * (nh-1))), # Stratified and double cases
                blgCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,blgPlot),
                               (sum(fa*blgPlot) - sum(nh * aStrat *blgStrat)) / (nh * (nh-1))), # Stratified and double cases
                bcCV = ifelse(first(ESTN_METHOD == 'simple'),
                              cov(fa,bcPlot),
                              (sum(fa*bcPlot) - sum(nh * aStrat *bcStrat)) / (nh * (nh-1))), # Stratified and double cases
                bpCV = ifelse(first(ESTN_METHOD == 'simple'),
                              cov(fa,bpPlot),
                              (sum(fa*bpPlot) - sum(nh * aStrat *bpStrat)) / (nh * (nh-1))), # Stratified and double cases
                bCV = ifelse(first(ESTN_METHOD == 'simple'),
                             cov(fa,bPlot),
                             (sum(fa*bPlot) - sum(nh * aStrat *bStrat)) / (nh * (nh-1))), # Stratified and double cases
                cdCV = ifelse(first(ESTN_METHOD == 'simple'),
                              cov(fa,cdPlot),
                              (sum(fa*cdPlot) - sum(nh * aStrat *cdStrat)) / (nh * (nh-1))), # Stratified and double cases
                clCV = ifelse(first(ESTN_METHOD == 'simple'),
                              cov(fa,blPlot),
                              (sum(fa*clPlot) - sum(nh * aStrat *clStrat)) / (nh * (nh-1))), # Stratified and double cases
                csmCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,csmPlot),
                               (sum(fa*csmPlot) - sum(nh * aStrat *csmStrat)) / (nh * (nh-1))), # Stratified and double cases
                cmdCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,cmdPlot),
                               (sum(fa*cmdPlot) - sum(nh * aStrat *cmdStrat)) / (nh * (nh-1))), # Stratified and double cases
                clgCV = ifelse(first(ESTN_METHOD == 'simple'),
                               cov(fa,clgPlot),
                               (sum(fa*clgPlot) - sum(nh * aStrat *clgStrat)) / (nh * (nh-1))), # Stratified and double cases
                ccCV = ifelse(first(ESTN_METHOD == 'simple'),
                              cov(fa,ccPlot),
                              (sum(fa*ccPlot) - sum(nh * aStrat *ccStrat)) / (nh * (nh-1))), # Stratified and double cases
                cpCV = ifelse(first(ESTN_METHOD == 'simple'),
                              cov(fa,cpPlot),
                              (sum(fa*cpPlot) - sum(nh * aStrat *cpStrat)) / (nh * (nh-1))), # Stratified and double cases
                cCV = ifelse(first(ESTN_METHOD == 'simple'),
                             cov(fa,cPlot),
                             (sum(fa*cPlot) - sum(nh * aStrat *cStrat)) / (nh * (nh-1)))) %>% # Stratified and double cases
      ## Estimation unit level
      group_by(ESTN_UNIT_CN) %>%
      summarize(aEst = unitMean(ESTN_METHOD, a, nh, w, aStrat),
                vsmEst = unitMean(ESTN_METHOD, a, nh, w, vsmStrat),
                vmdEst = unitMean(ESTN_METHOD, a, nh, w, vmdStrat),
                vlgEst = unitMean(ESTN_METHOD, a, nh, w, vlgStrat),
                vcEst = unitMean(ESTN_METHOD, a, nh, w, vcStrat),
                vpEst = unitMean(ESTN_METHOD, a, nh, w, vpStrat),
                vEst = unitMean(ESTN_METHOD, a, nh, w, vStrat),
                bdEst = unitMean(ESTN_METHOD, a, nh, w, bdStrat),
                blEst = unitMean(ESTN_METHOD, a, nh, w, blStrat),
                bsmEst = unitMean(ESTN_METHOD, a, nh, w, bsmStrat),
                bmdEst = unitMean(ESTN_METHOD, a, nh, w, bmdStrat),
                blgEst = unitMean(ESTN_METHOD, a, nh, w, blgStrat),
                bcEst = unitMean(ESTN_METHOD, a, nh, w, bcStrat),
                bpEst = unitMean(ESTN_METHOD, a, nh, w, bpStrat),
                bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
                cdEst = unitMean(ESTN_METHOD, a, nh, w, cdStrat),
                clEst = unitMean(ESTN_METHOD, a, nh, w, clStrat),
                csmEst = unitMean(ESTN_METHOD, a, nh, w, csmStrat),
                cmdEst = unitMean(ESTN_METHOD, a, nh, w, cmdStrat),
                clgEst = unitMean(ESTN_METHOD, a, nh, w, clgStrat),
                ccEst = unitMean(ESTN_METHOD, a, nh, w, ccStrat),
                cpEst = unitMean(ESTN_METHOD, a, nh, w, cpStrat),
                cEst = unitMean(ESTN_METHOD, a, nh, w, cStrat),
                plotIn = sum(plotIn, na.rm = TRUE),
                # Estimation of unit variance
                aVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, av, aStrat, aEst),
                vsmVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, vsmV, vsmStrat, vsmEst),
                vmdVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, vmdV, vmdStrat, vmdEst),
                vlgVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, vlgV, vlgStrat, vlgEst),
                vcVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, vcV, vcStrat, vcEst),
                vpVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, vpV, vpStrat, vpEst),
                vVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, vV, vStrat, vEst),
                bdVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bdV, bdStrat, bdEst),
                blVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, blV, blStrat, blEst),
                bsmVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bsmV, bsmStrat, bsmEst),
                bmdVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bmdV, bmdStrat, bmdEst),
                blgVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, blgV, blgStrat, blgEst),
                bcVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bcV, bcStrat, bcEst),
                bpVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bpV, bpStrat, bpEst),
                bVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bV, bStrat, cbEst),
                cdVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, cdV, cdStrat, cdEst),
                clVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, clV, clStrat, clEst),
                csmVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, csmV, csmStrat, csmEst),
                cmdVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, cmdV, cmdStrat, cmdEst),
                clgVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, clgV, clgStrat, clgEst),
                ccVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, ccV, ccStrat, ccEst),
                cpVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, cpV, cpStrat, cpEst),
                cVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, cV, cStrat, cEst),
                # Estimation of unit covariance
                vsmCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, vsmCV, vsmStrat, vsmEst, aStrat, aEst),
                vmdCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, vmdCV, vmdStrat, vmdEst, aStrat, aEst),
                vlgCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, vlgCV, vlgStrat, vlgEst, aStrat, aEst),
                vcCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, vcCV, vcStrat, vcEst, aStrat, aEst),
                vpCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, vpCV, vpStrat, vpEst, aStrat, aEst),
                vCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, vCV, vStrat, vEst, aStrat, aEst),
                bdCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bdCV, bdStrat, bdEst, aStrat, aEst),
                blCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, blCV, blStrat, blEst, aStrat, aEst),
                bsmCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bsmCV, bsmStrat, bsmEst, aStrat, aEst),
                bmdCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bmdCV, bmdStrat, bmdEst, aStrat, aEst),
                blgCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, blgCV, blgStrat, blgEst, aStrat, aEst),
                bcCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bcCV, bcStrat, bcEst, aStrat, aEst),
                bpCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bpCV, bpStrat, bpEst, aStrat, aEst),
                bCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bCV, bStrat, cbEst, aStrat, aEst),
                cdCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cdCV, cdStrat, cdEst, aStrat, aEst),
                clCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, clCV, clStrat, clEst, aStrat, aEst),
                csmCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, csmCV, csmStrat, csmEst, aStrat, aEst),
                cmdCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cmdCV, cmdStrat, cmdEst, aStrat, aEst),
                clgCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, clgCV, clgStrat, clgEst, aStrat, aEst),
                ccCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, ccCV, ccStrat, ccEst, aStrat, aEst),
                cpCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cpCV, cpStrat, cpEst, aStrat, aEst),
                cCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cCV, cStrat, cEst, aStrat, aEst)) %>%
      # Compute totals
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                VOL_DUFF = NA,
                VOL_LITTER = NA,
                VOL_1HR = sum(vsmEst, na.rm = TRUE),
                VOL_10HR = sum(vmdEst, na.rm = TRUE),
                VOL_100HR = sum(vlgEst, na.rm = TRUE),
                VOL_1000HR = sum(vcEst, na.rm = TRUE),
                VOL_PILE = sum(vpEst, na.rm = TRUE),
                VOL = sum(vEst, na.rm = TRUE),
                BIO_DUFF = sum(bdEst, na.rm = TRUE),
                BIO_LITTER = sum(blEst, na.rm = TRUE),
                BIO_1HR = sum(bsmEst, na.rm = TRUE),
                BIO_10HR = sum(bmdEst, na.rm = TRUE),
                BIO_100HR = sum(blgEst, na.rm = TRUE),
                BIO_1000HR = sum(bcEst, na.rm = TRUE),
                BIO_PILE = sum(bpEst, na.rm = TRUE),
                BIO = sum(bEst, na.rm = TRUE),
                CARB_DUFF = sum(cdEst, na.rm = TRUE),
                CARB_LITTER = sum(clEst, na.rm = TRUE),
                CARB_1HR = sum(csmEst, na.rm = TRUE),
                CARB_10HR = sum(cmdEst, na.rm = TRUE),
                CARB_100HR = sum(clgEst, na.rm = TRUE),
                CARB_1000HR = sum(ccEst, na.rm = TRUE),
                CARB_PILE = sum(cpEst, na.rm = TRUE),
                CARB = sum(cEst, na.rm = TRUE),
                nPlots = sum(plotIn, na.rm = TRUE),
                # Per Acre
                VOL_DUFF_ACRE = NA,
                VOL_LITTER_ACRE = NA,
                VOL_1HR_ACRE = VOL_1HR / AREA_TOTAL,
                VOL_10HR_ACRE = VOL_10HR / AREA_TOTAL,
                VOL_100HR_ACRE = VOL_100HR / AREA_TOTAL,
                VOL_1000HR_ACRE = VOL_1000HR / AREA_TOTAL,
                VOL_PILE_ACRE = VOL_PILE / AREA_TOTAL,
                VOL_ACRE = VOL / AREA_TOTAL,
                BIO_DUFF_ACRE = BIO_DUFF / AREA_TOTAL,
                BIO_LITTER_ACRE = BIO_LITTER / AREA_TOTAL,
                BIO_1HR_ACRE = BIO_1HR / AREA_TOTAL,
                BIO_10HR_ACRE = BIO_10HR / AREA_TOTAL,
                BIO_100HR_ACRE = BIO_100HR / AREA_TOTAL,
                BIO_1000HR_ACRE = BIO_1000HR / AREA_TOTAL,
                BIO_PILE_ACRE = BIO_PILE / AREA_TOTAL,
                BIO_ACRE = BIO / AREA_TOTAL,
                CARB_DUFF_ACRE = CARB_DUFF / AREA_TOTAL,
                CARB_LITTER_ACRE = CARB_LITTER / AREA_TOTAL,
                CARB_1HR_ACRE = CARB_1HR / AREA_TOTAL,
                CARB_10HR_ACRE = CARB_10HR / AREA_TOTAL,
                CARB_100HR_ACRE = CARB_100HR / AREA_TOTAL,
                CARB_1000HR_ACRE = CARB_1000HR / AREA_TOTAL,
                CARB_PILE_ACRE = CARB_PILE / AREA_TOTAL,
                CARB_ACRE = CARB / AREA_TOTAL,
                # Sampling Errors totals
                AREA_TOTAL_SE = sqrt(sum(aVar, na.rm = TRUE)) / AREA_TOTAL * 100,
                VOL_DUFF_SE = NA,
                VOL_LITTER_SE = NA,
                VOL_1HR_SE = sqrt(sum(vsmVar, na.rm = TRUE)) / VOL_1HR * 100,
                VOL_10HR_SE = sqrt(sum(vmdVar, na.rm = TRUE)) / VOL_10HR * 100,
                VOL_100HR_SE = sqrt(sum(vlgVar, na.rm = TRUE)) / VOL_100HR * 100,
                VOL_1000HR_SE = sqrt(sum(vcVar, na.rm = TRUE)) / VOL_1000HR * 100,
                VOL_PILE_SE = sqrt(sum(vpVar, na.rm = TRUE)) / VOL_PILE * 100,
                VOL_SE = sqrt(sum(vVar, na.rm = TRUE)) / VOL * 100,
                BIO_DUFF_SE = sqrt(sum(bdVar, na.rm = TRUE)) / BIO_DUFF * 100,
                BIO_LITTER_SE = sqrt(sum(blVar, na.rm = TRUE)) / BIO_LITTER * 100,
                BIO_1HR_SE = sqrt(sum(bsmVar, na.rm = TRUE)) / BIO_1HR * 100,
                BIO_10HR_SE = sqrt(sum(bmdVar, na.rm = TRUE)) / BIO_10HR * 100,
                BIO_100HR_SE = sqrt(sum(blgVar, na.rm = TRUE)) / BIO_100HR * 100,
                BIO_1000HR_SE = sqrt(sum(bcVar, na.rm = TRUE)) / BIO_1000HR * 100,
                BIO_PILE_SE = sqrt(sum(bpVar, na.rm = TRUE)) / BIO_PILE * 100,
                BIO_SE = sqrt(sum(bVar, na.rm = TRUE)) / BIO * 100,
                CARB_DUFF_SE = sqrt(sum(cdVar, na.rm = TRUE)) / CARB_DUFF * 100,
                CARB_LITTER_SE = sqrt(sum(clVar, na.rm = TRUE)) / CARB_LITTER * 100,
                CARB_1HR_SE = sqrt(sum(csmVar, na.rm = TRUE)) / CARB_1HR * 100,
                CARB_10HR_SE = sqrt(sum(cmdVar, na.rm = TRUE)) / CARB_10HR * 100,
                CARB_100HR_SE = sqrt(sum(clgVar, na.rm = TRUE)) / CARB_100HR * 100,
                CARB_1000HR_SE = sqrt(sum(ccVar, na.rm = TRUE)) / CARB_1000HR * 100,
                CARB_PILE_SE = sqrt(sum(cpVar, na.rm = TRUE)) / CARB_PILE * 100,
                CARB_SE = sqrt(sum(cVar, na.rm = TRUE)) / CARB * 100,
                # Per Acre variances
                vsmVar = (1/AREA_TOTAL^2) * (sum(vsmVar, na.rm = TRUE) + (VOL_1HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_1HR_ACRE * sum(vsmCV, na.rm = TRUE))),
                vmdVar = (1/AREA_TOTAL^2) * (sum(vmdVar, na.rm = TRUE) + (VOL_10HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_10HR_ACRE * sum(vmdCV, na.rm = TRUE))),
                vlgVar = (1/AREA_TOTAL^2) * (sum(vlgVar, na.rm = TRUE) + (VOL_100HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_100HR_ACRE * sum(vlgCV, na.rm = TRUE))),
                vcVar = (1/AREA_TOTAL^2) * (sum(vcVar, na.rm = TRUE) + (VOL_1000HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_1000HR_ACRE *sum(vcCV, na.rm = TRUE))),
                vpVar = (1/AREA_TOTAL^2) * (sum(vpVar, na.rm = TRUE) + (VOL_PILE_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_PILE_ACRE * sum(vpCV, na.rm = TRUE))),
                vVar = (1/AREA_TOTAL^2) * (sum(vVar, na.rm = TRUE) + (VOL_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_ACRE * sum(vCV, na.rm = TRUE))),
                bdVar = (1/AREA_TOTAL^2) * (sum(bdVar, na.rm = TRUE) + (BIO_DUFF_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_DUFF_ACRE * sum(bdCV, na.rm = TRUE))),
                blVar = (1/AREA_TOTAL^2) * (sum(blVar, na.rm = TRUE) + (BIO_LITTER_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_LITTER_ACRE * sum(blCV, na.rm = TRUE))),
                bsmVar = (1/AREA_TOTAL^2) * (sum(bsmVar, na.rm = TRUE) + (BIO_1HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_1HR_ACRE * sum(bsmCV, na.rm = TRUE))),
                bmdVar = (1/AREA_TOTAL^2) * (sum(bmdVar, na.rm = TRUE) + (BIO_10HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_10HR_ACRE * sum(bmdCV, na.rm = TRUE))),
                blgVar = (1/AREA_TOTAL^2) * (sum(blgVar, na.rm = TRUE) + (BIO_100HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_100HR_ACRE * sum(blgCV, na.rm = TRUE))),
                bcVar = (1/AREA_TOTAL^2) * (sum(bcVar, na.rm = TRUE) + (BIO_1000HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_1000HR_ACRE * sum(bcCV, na.rm = TRUE))),
                bpVar = (1/AREA_TOTAL^2) * (sum(bpVar, na.rm = TRUE) + (BIO_PILE_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_PILE_ACRE * sum(bpCV, na.rm = TRUE))),
                bVar = (1/AREA_TOTAL^2) * (sum(bVar, na.rm = TRUE) + (BIO_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_ACRE * sum(bCV, na.rm = TRUE))),
                cdVar = (1/AREA_TOTAL^2) * (sum(cdVar, na.rm = TRUE) + (CARB_DUFF_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_DUFF_ACRE * sum(cdCV, na.rm = TRUE))),
                clVar = (1/AREA_TOTAL^2) * (sum(clVar, na.rm = TRUE) + (CARB_LITTER_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_LITTER_ACRE * sum(clCV, na.rm = TRUE))),
                csmVar = (1/AREA_TOTAL^2) * (sum(csmVar, na.rm = TRUE) + (CARB_1HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_1HR_ACRE * sum(csmCV, na.rm = TRUE))),
                cmdVar = (1/AREA_TOTAL^2) * (sum(cmdVar, na.rm = TRUE) + (CARB_10HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_10HR_ACRE * sum(cmdCV, na.rm = TRUE))),
                clgVar = (1/AREA_TOTAL^2) * (sum(clgVar, na.rm = TRUE) + (CARB_100HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_100HR_ACRE * sum(clgCV, na.rm = TRUE))),
                ccVar = (1/AREA_TOTAL^2) * (sum(ccVar, na.rm = TRUE) + (CARB_1000HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_1000HR_ACRE * sum(ccCV, na.rm = TRUE))),
                cpVar = (1/AREA_TOTAL^2) * (sum(cpVar, na.rm = TRUE) + (CARB_PILE_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_PILE_ACRE * sum(cpCV, na.rm = TRUE))),
                cVar = (1/AREA_TOTAL^2) * (sum(cVar, na.rm = TRUE) + (CARB_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_ACRE * sum(cCV, na.rm = TRUE))),
                # Per acre sampling errors
                VOL_DUFF_ACRE_SE = NA,
                VOL_LITTER_ACRE_SE = NA,
                VOL_1HR_ACRE_SE = sqrt(sum(vsmVar, na.rm = TRUE)) / VOL_1HR_ACRE * 100,
                VOL_10HR_ACRE_SE = sqrt(sum(vmdVar, na.rm = TRUE)) / VOL_10HR_ACRE * 100,
                VOL_100HR_ACRE_SE = sqrt(sum(vlgVar, na.rm = TRUE)) / VOL_100HR_ACRE * 100,
                VOL_1000HR_ACRE_SE = sqrt(sum(vcVar, na.rm = TRUE)) / VOL_1000HR_ACRE * 100,
                VOL_PILE_ACRE_SE = sqrt(sum(vpVar, na.rm = TRUE)) / VOL_PILE_ACRE * 100,
                VOL_ACRE_SE = sqrt(sum(vVar, na.rm = TRUE)) / VOL_ACRE * 100,
                BIO_DUFF_ACRE_SE = sqrt(sum(bdVar, na.rm = TRUE)) / BIO_DUFF_ACRE * 100,
                BIO_LITTER_ACRE_SE = sqrt(sum(blVar, na.rm = TRUE)) / BIO_LITTER_ACRE * 100,
                BIO_1HR_ACRE_SE = sqrt(sum(bsmVar, na.rm = TRUE)) / BIO_1HR_ACRE * 100,
                BIO_10HR_ACRE_SE = sqrt(sum(bmdVar, na.rm = TRUE)) / BIO_10HR_ACRE * 100,
                BIO_100HR_ACRE_SE = sqrt(sum(blgVar, na.rm = TRUE)) / BIO_100HR_ACRE * 100,
                BIO_1000HR_ACRE_SE = sqrt(sum(bcVar, na.rm = TRUE)) / BIO_1000HR_ACRE * 100,
                BIO_PILE_ACRE_SE = sqrt(sum(bpVar, na.rm = TRUE)) / BIO_PILE_ACRE * 100,
                BIO_ACRE_SE = sqrt(sum(bVar, na.rm = TRUE)) / BIO_ACRE * 100,
                CARB_DUFF_ACRE_SE = sqrt(sum(cdVar, na.rm = TRUE)) / CARB_DUFF_ACRE * 100,
                CARB_LITTER_ACRE_SE = sqrt(sum(clVar, na.rm = TRUE)) / CARB_LITTER_ACRE * 100,
                CARB_1HR_ACRE_SE = sqrt(sum(csmVar, na.rm = TRUE)) / CARB_1HR_ACRE * 100,
                CARB_10HR_ACRE_SE = sqrt(sum(cmdVar, na.rm = TRUE)) / CARB_10HR_ACRE * 100,
                CARB_100HR_ACRE_SE = sqrt(sum(clgVar, na.rm = TRUE)) / CARB_100HR_ACRE * 100,
                CARB_1000HR_ACRE_SE = sqrt(sum(ccVar, na.rm = TRUE)) / CARB_1000HR_ACRE * 100,
                CARB_PILE_ACRE_SE = sqrt(sum(cpVar, na.rm = TRUE)) / CARB_PILE_ACRE * 100,
                CARB_ACRE_SE = sqrt(sum(cVar, na.rm = TRUE)) / CARB_ACRE * 100)

    if (totals) {
      cwd <- cwd %>%
        select(names(cwd)[str_detect(names(cwd), 'Var', negate = TRUE)], nPlots)
    } else {
      cwd <- cwd %>%
        select(names(cwd)[str_detect(names(cwd), 'Var', negate = TRUE) & str_detect(names(cwd), 'ACRE')], nPlots)
    }

    # Rejoin w/ groupby names
    cwd <- data.frame(combos[[x]], cwd)

  } else {
    ### BELOW DOES NOT PRODUCE SAMPLING ERRORS, use EXPNS instead (much quicker)
    cwd <- data %>%
      filter(EVAL_TYP == 'EXPDWM') %>%
      distinct(ESTN_UNIT_CN, STRATUM_CN, PLT_CN, CND_CN, .keep_all = TRUE) %>%
      group_by(.dots = grpBy, PLT_CN) %>%
      summarize(VOL_1HR = sum(FWD_SM_VOLCF_ADJ * aDI * EXPNS, na.rm = TRUE),
                VOL_10HR = sum(FWD_MD_VOLCF_ADJ * aDI * EXPNS, na.rm = TRUE),
                VOL_100HR = sum(FWD_LG_VOLCF_ADJ * aDI * EXPNS, na.rm = TRUE),
                VOL_1000HR = sum(CWD_VOLCF_ADJ * aDI * EXPNS, na.rm = TRUE),
                VOL_PILE = sum(PILE_VOLCF_ADJ * aDI * EXPNS, na.rm = TRUE),
                VOL = sum(VOL_1HR, VOL_10HR, VOL_100HR, VOL_1000HR, VOL_PILE, na.rm = TRUE),
                BIO_DUFF = sum(DUFF_BIOMASS* aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_LITTER = sum(LITTER_BIOMASS * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_1HR = sum(FWD_SM_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_10HR = sum(FWD_MD_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_100HR = sum(FWD_LG_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_1000HR = sum(CWD_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_PILE = sum(PILE_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO = sum(BIO_DUFF, BIO_LITTER, BIO_1HR, BIO_10HR, BIO_100HR, BIO_1000HR, BIO_PILE, na.rm = TRUE),
                CARB_DUFF = sum(DUFF_CARBON* aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_LITTER = sum(LITTER_CARBON * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_1HR = sum(FWD_SM_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_10HR = sum(FWD_MD_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_100HR = sum(FWD_LG_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_1000HR = sum(CWD_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_PILE = sum(PILE_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB = sum(CARB_DUFF, CARB_LITTER, CARB_1HR, CARB_10HR, CARB_100HR, CARB_1000HR, CARB_PILE, na.rm = TRUE),
                fa = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      group_by(.dots = grpBy) %>%
      summarize(AREA_TOTAL = sum(fa, na.rm = TRUE),
                VOL_DUFF = NA,
                VOL_LITTER = NA,
                VOL_1HR = sum(VOL_1HR, na.rm = TRUE),
                VOL_10HR = sum(VOL_10HR, na.rm = TRUE),
                VOL_100HR = sum(VOL_100HR, na.rm = TRUE),
                VOL_1000HR = sum(VOL_1000HR, na.rm = TRUE),
                VOL_PILE = sum(VOL_PILE, na.rm = TRUE),
                VOL = sum(VOL, na.rm = TRUE),
                BIO_DUFF = sum(BIO_DUFF, na.rm = TRUE),
                BIO_LITTER = sum(BIO_LITTER, na.rm = TRUE),
                BIO_1HR = sum(BIO_1HR, na.rm = TRUE),
                BIO_10HR = sum(BIO_10HR, na.rm = TRUE),
                BIO_100HR = sum(BIO_100HR, na.rm = TRUE),
                BIO_1000HR = sum(BIO_1000HR, na.rm = TRUE),
                BIO_PILE = sum(BIO_PILE, na.rm = TRUE),
                BIO = sum(BIO, na.rm = TRUE),
                CARB_DUFF = sum(CARB_DUFF, na.rm = TRUE),
                CARB_LITTER = sum(CARB_LITTER, na.rm = TRUE),
                CARB_1HR = sum(CARB_1HR, na.rm = TRUE),
                CARB_10HR = sum(CARB_10HR, na.rm = TRUE),
                CARB_100HR = sum(CARB_100HR, na.rm = TRUE),
                CARB_1000HR = sum(CARB_1000HR, na.rm = TRUE),
                CARB_PILE = sum(CARB_PILE, na.rm = TRUE),
                CARB = sum(CARB, na.rm = TRUE),
                VOL_DUFF_ACRE = NA,
                VOL_LITTER_ACRE = NA,
                VOL_1HR_ACRE = VOL_1HR / AREA_TOTAL,
                VOL_10HR_ACRE = VOL_10HR / AREA_TOTAL,
                VOL_100HR_ACRE = VOL_100HR / AREA_TOTAL,
                VOL_1000HR_ACRE = VOL_1000HR / AREA_TOTAL,
                VOL_PILE_ACRE = VOL_PILE / AREA_TOTAL,
                VOL_ACRE = VOL / AREA_TOTAL,
                BIO_DUFF_ACRE = BIO_DUFF / AREA_TOTAL,
                BIO_LITTER_ACRE = BIO_LITTER / AREA_TOTAL,
                BIO_1HR_ACRE = BIO_1HR / AREA_TOTAL,
                BIO_10HR_ACRE = BIO_10HR / AREA_TOTAL,
                BIO_100HR_ACRE = BIO_100HR / AREA_TOTAL,
                BIO_1000HR_ACRE = BIO_1000HR / AREA_TOTAL,
                BIO_PILE_ACRE = BIO_PILE / AREA_TOTAL,
                BIO_ACRE = BIO / AREA_TOTAL,
                CARB_DUFF_ACRE = CARB_DUFF / AREA_TOTAL,
                CARB_LITTER_ACRE = CARB_LITTER / AREA_TOTAL,
                CARB_1HR_ACRE = CARB_1HR / AREA_TOTAL,
                CARB_10HR_ACRE = CARB_10HR / AREA_TOTAL,
                CARB_100HR_ACRE = CARB_100HR / AREA_TOTAL,
                CARB_1000HR_ACRE = CARB_1000HR / AREA_TOTAL,
                CARB_PILE_ACRE = CARB_PILE / AREA_TOTAL,
                CARB_ACRE = CARB / AREA_TOTAL,
                nPlots = sum(plotIn, na.rm = TRUE))
    # Remove the total values if told to do so
    if (totals) {
      cwd <- cwd %>%
        select(grpBy, names(cwd)[str_detect(names(cwd), 'Var', negate = TRUE)], nPlots)
    } else {
      cwd <- cwd %>%
        select(grpBy, names(cwd)[str_detect(names(cwd), 'Var', negate = TRUE) & str_detect(names(cwd), 'ACRE')], nPlots)
    }
  } # End SE Conditional

  # Do some cleanup
  #gc()

  #Return a dataframe
  cwd

}
