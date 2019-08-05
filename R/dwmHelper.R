dwmHelper <- function(x, combos, data, grpBy, totals, SE){
  # Update domain indicator for each each column speficed in grpBy
  ad = 1
  for (n in 1:ncol(combos[[x]])){
    # Area domain indicator for each column in
    aObs <- combos[[x]][[grpBy[n]]] == data[[grpBy[n]]]
    ad <- data$aDI * aObs * ad
  }

  if(SE){
    data$aDI <- ad
    data$aDI[is.na(data$aDI)] <- 0
    cwd <- data %>%
      filter(EVAL_TYP == 'EXPDWM') %>%
      distinct(PLT_CN, CONDID, EVALID, COND_STATUS_CD, .keep_all = TRUE) %>%
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN) %>%
      summarize(vsmPlot = sum(FWD_SM_VOLCF_ADJ * aDI, na.rm = TRUE),
                vmdPlot = sum(FWD_MD_VOLCF_ADJ * aDI, na.rm = TRUE),
                vlgPlot = sum(FWD_LG_VOLCF_ADJ * aDI, na.rm = TRUE),
                vcPlot = sum(CWD_VOLCF_ADJ * aDI, na.rm = TRUE),
                vpPlot = sum(PILE_VOLCF_ADJ * aDI, na.rm = TRUE),
                vPlot = sum(vsmPlot, vmdPlot, vlgPlot, vcPlot, vpPlot, na.rm = TRUE),
                bsmPlot = sum(FWD_SM_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bmdPlot = sum(FWD_MD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                blgPlot = sum(FWD_LG_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bcPlot = sum(CWD_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bpPlot = sum(PILE_DRYBIO_ADJ * aDI / 2000, na.rm = TRUE),
                bPlot = sum(bsmPlot, bmdPlot, blgPlot, bcPlot, bpPlot, na.rm = TRUE),
                csmPlot = sum(FWD_SM_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                cmdPlot = sum(FWD_MD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                clgPlot = sum(FWD_LG_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                ccPlot = sum(CWD_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                cpPlot = sum(PILE_CARBON_ADJ * aDI / 2000, na.rm = TRUE),
                cPlot = sum(csmPlot, cmdPlot, clgPlot, ccPlot, cpPlot, na.rm = TRUE),
                fa = sum(CONDPROP_UNADJ * aDI * aAdj, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0),
                a = first(AREA_USED),
                p1EU = first(P1PNTCNT_EU),
                p1 = first(P1POINTCNT),
                p2 = first(P2POINTCNT)) %>%
      # Stratum level
      group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN) %>%
      summarize(aStrat = mean(fa, na.rm = TRUE),
                vsmStrat = mean(vsmPlot, na.rm = TRUE),
                vmdStrat = mean(vmdPlot, na.rm = TRUE),
                vlgStrat = mean(vlgPlot, na.rm = TRUE),
                vcStrat = mean(vcPlot, na.rm = TRUE),
                vpStrat = mean(vpPlot, na.rm = TRUE),
                vStrat = mean(vPlot, na.rm = TRUE),
                bsmStrat = mean(bsmPlot, na.rm = TRUE),
                bmdStrat = mean(bmdPlot, na.rm = TRUE),
                blgStrat = mean(blgPlot, na.rm = TRUE),
                bcStrat = mean(bcPlot, na.rm = TRUE),
                bpStrat = mean(bpPlot, na.rm = TRUE),
                bStrat = mean(bPlot, na.rm = TRUE),
                csmStrat = mean(csmPlot, na.rm = TRUE),
                cmdStrat = mean(cmdPlot, na.rm = TRUE),
                clgStrat = mean(clgPlot, na.rm = TRUE),
                ccStrat = mean(ccPlot, na.rm = TRUE),
                cpStrat = mean(cpPlot, na.rm = TRUE),
                cStrat = mean(cPlot, na.rm = TRUE),
                plotIn = sum(plotIn, na.rm = TRUE),
                a = first(a),
                w = first(p1) / first(p1EU), # Stratum weight
                nh = first(p2), # Number plots in stratum
                # Strata level variances
                av = ifelse(first(ESTN_METHOD == 'simple'),
                            var(fa * first(a) / nh),
                            (sum(fa^2) - sum(nh * aStrat^2)) / (nh * (nh-1))),
                vsmV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(vsmPlot * first(a) / nh),
                              (sum(vsmPlot^2) - sum(nh * vsmStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vmdV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(vmdPlot * first(a) / nh),
                              (sum(vmdPlot^2) - sum(nh * vmdStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vlgV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(vlgPlot * first(a) / nh),
                              (sum(vlgPlot^2) - sum(nh * vlgStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vcV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(vcPlot * first(a) / nh),
                             (sum(vcPlot^2) - sum(nh * vcStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vpV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(vpPlot * first(a) / nh),
                             (sum(vpPlot^2) - sum(nh * vpStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                vV = ifelse(first(ESTN_METHOD == 'simple'),
                            var(vPlot * first(a) / nh),
                            (sum(vPlot^2) - sum(nh * vStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bsmV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(bsmPlot * first(a) / nh),
                              (sum(bsmPlot^2) - sum(nh * bsmStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bmdV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(bmdPlot * first(a) / nh),
                              (sum(bmdPlot^2) - sum(nh * bmdStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                blgV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(blgPlot * first(a) / nh),
                              (sum(blgPlot^2) - sum(nh * blgStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bcV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(bcPlot * first(a) / nh),
                             (sum(bcPlot^2) - sum(nh * bcStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bpV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(bpPlot * first(a) / nh),
                             (sum(bpPlot^2) - sum(nh * bpStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                bV = ifelse(first(ESTN_METHOD == 'simple'),
                            var(bPlot * first(a) / nh),
                            (sum(bPlot^2) - sum(nh * bStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                csmV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(csmPlot * first(a) / nh),
                              (sum(csmPlot^2) - sum(nh * csmStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cmdV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(cmdPlot * first(a) / nh),
                              (sum(cmdPlot^2) - sum(nh * cmdStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                clgV = ifelse(first(ESTN_METHOD == 'simple'),
                              var(clgPlot * first(a) / nh),
                              (sum(clgPlot^2) - sum(nh * clgStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                ccV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(ccPlot * first(a) / nh),
                             (sum(ccPlot^2) - sum(nh * ccStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cpV = ifelse(first(ESTN_METHOD == 'simple'),
                             var(cpPlot * first(a) / nh),
                             (sum(cpPlot^2) - sum(nh * cpStrat^2)) / (nh * (nh-1))), # Stratified and double cases
                cV = ifelse(first(ESTN_METHOD == 'simple'),
                            var(cPlot * first(a) / nh),
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
                bsmEst = unitMean(ESTN_METHOD, a, nh, w, bsmStrat),
                bmdEst = unitMean(ESTN_METHOD, a, nh, w, bmdStrat),
                blgEst = unitMean(ESTN_METHOD, a, nh, w, blgStrat),
                bcEst = unitMean(ESTN_METHOD, a, nh, w, bcStrat),
                bpEst = unitMean(ESTN_METHOD, a, nh, w, bpStrat),
                bEst = unitMean(ESTN_METHOD, a, nh, w, bStrat),
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
                bsmVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bsmV, bsmStrat, bsmEst),
                bmdVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bmdV, bmdStrat, bmdEst),
                blgVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, blgV, blgStrat, blgEst),
                bcVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bcV, bcStrat, bcEst),
                bpVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bpV, bpStrat, bpEst),
                bVar = unitVar(method = 'var', ESTN_METHOD, a, nh, w, bV, bStrat, cbEst),
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
                bsmCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bsmCV, bsmStrat, bsmEst, aStrat, aEst),
                bmdCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bmdCV, bmdStrat, bmdEst, aStrat, aEst),
                blgCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, blgCV, blgStrat, blgEst, aStrat, aEst),
                bcCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bcCV, bcStrat, bcEst, aStrat, aEst),
                bpCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bpCV, bpStrat, bpEst, aStrat, aEst),
                bCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, bCV, bStrat, cbEst, aStrat, aEst),
                csmCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, csmCV, csmStrat, csmEst, aStrat, aEst),
                cmdCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cmdCV, cmdStrat, cmdEst, aStrat, aEst),
                clgCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, clgCV, clgStrat, clgEst, aStrat, aEst),
                ccCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, ccCV, ccStrat, ccEst, aStrat, aEst),
                cpCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cpCV, cpStrat, cpEst, aStrat, aEst),
                cCV = unitVar(method = 'cov', ESTN_METHOD, a, nh, w, cCV, cStrat, cEst, aStrat, aEst)) %>%
      # Compute totals
      summarize(AREA_TOTAL = sum(aEst, na.rm = TRUE),
                VOL_1HR = sum(vsmEst, na.rm = TRUE),
                VOL_10HR = sum(vmdEst, na.rm = TRUE),
                VOL_100HR = sum(vlgEst, na.rm = TRUE),
                VOL_1000HR = sum(vcEst, na.rm = TRUE),
                VOL_PILE = sum(vpEst, na.rm = TRUE),
                VOL = sum(vEst, na.rm = TRUE),
                BIO_1HR = sum(bsmEst, na.rm = TRUE),
                BIO_10HR = sum(bmdEst, na.rm = TRUE),
                BIO_100HR = sum(blgEst, na.rm = TRUE),
                BIO_1000HR = sum(bcEst, na.rm = TRUE),
                BIO_PILE = sum(bpEst, na.rm = TRUE),
                BIO = sum(bEst, na.rm = TRUE),
                CARB_1HR = sum(csmEst, na.rm = TRUE),
                CARB_10HR = sum(cmdEst, na.rm = TRUE),
                CARB_100HR = sum(clgEst, na.rm = TRUE),
                CARB_1000HR = sum(ccEst, na.rm = TRUE),
                CARB_PILE = sum(cpEst, na.rm = TRUE),
                CARB = sum(cEst, na.rm = TRUE),
                nPlots = sum(plotIn, na.rm = TRUE),
                # Per Acre
                VOL_1HR_ACRE = VOL_1HR / AREA_TOTAL,
                VOL_10HR_ACRE = VOL_10HR / AREA_TOTAL,
                VOL_100HR_ACRE = VOL_100HR / AREA_TOTAL,
                VOL_1000HR_ACRE = VOL_1000HR / AREA_TOTAL,
                VOL_PILE_ACRE = VOL_PILE / AREA_TOTAL,
                VOL_ACRE = VOL / AREA_TOTAL,
                BIO_1HR_ACRE = BIO_1HR / AREA_TOTAL,
                BIO_10HR_ACRE = BIO_10HR / AREA_TOTAL,
                BIO_100HR_ACRE = BIO_100HR / AREA_TOTAL,
                BIO_1000HR_ACRE = BIO_1000HR / AREA_TOTAL,
                BIO_PILE_ACRE = BIO_PILE / AREA_TOTAL,
                BIO_ACRE = BIO / AREA_TOTAL,
                CARB_1HR_ACRE = CARB_1HR / AREA_TOTAL,
                CARB_10HR_ACRE = CARB_10HR / AREA_TOTAL,
                CARB_100HR_ACRE = CARB_100HR / AREA_TOTAL,
                CARB_1000HR_ACRE = CARB_1000HR / AREA_TOTAL,
                CARB_PILE_ACRE = CARB_PILE / AREA_TOTAL,
                CARB_ACRE = CARB / AREA_TOTAL,
                # Sampling Errors totals
                AREA_TOTAL_SE = sqrt(sum(aVar, na.rm = TRUE)) / AREA_TOTAL * 100,
                VOL_1HR_SE = ifelse(nPlots > 1, sqrt(sum(vsmVar, na.rm = TRUE)) / VOL_1HR * 100,0),
                VOL_10HR_SE = ifelse(nPlots > 1, sqrt(sum(vmdVar, na.rm = TRUE)) / VOL_10HR * 100,0),
                VOL_100HR_SE = ifelse(nPlots > 1, sqrt(sum(vlgVar, na.rm = TRUE)) / VOL_100HR * 100,0),
                VOL_1000HR_SE = ifelse(nPlots > 1, sqrt(sum(vcVar, na.rm = TRUE)) / VOL_1000HR * 100,0),
                VOL_PILE_SE = ifelse(nPlots > 1,sqrt(sum(vpVar, na.rm = TRUE)) / VOL_PILE * 100,0),
                VOL_SE = ifelse(nPlots > 1, sqrt(sum(vVar, na.rm = TRUE)) / VOL * 100,0),
                BIO_1HR_SE = ifelse(nPlots > 1, sqrt(sum(bsmVar, na.rm = TRUE)) / BIO_1HR * 100,0),
                BIO_10HR_SE = ifelse(nPlots > 1, sqrt(sum(bmdVar, na.rm = TRUE)) / BIO_10HR * 100,0),
                BIO_100HR_SE = ifelse(nPlots > 1, sqrt(sum(blgVar, na.rm = TRUE)) / BIO_100HR * 100,0),
                BIO_1000HR_SE = ifelse(nPlots > 1, sqrt(sum(bcVar, na.rm = TRUE)) / BIO_1000HR * 100,0),
                BIO_PILE_SE = ifelse(nPlots > 1, sqrt(sum(bpVar, na.rm = TRUE)) / BIO_PILE * 100,0),
                BIO_SE = ifelse(nPlots > 1, sqrt(sum(bVar, na.rm = TRUE)) / BIO * 100,0),
                CARB_1HR_SE = ifelse(nPlots > 1, sqrt(sum(csmVar, na.rm = TRUE)) / CARB_1HR * 100,0),
                CARB_10HR_SE = ifelse(nPlots > 1, sqrt(sum(cmdVar, na.rm = TRUE)) / CARB_10HR * 100,0),
                CARB_100HR_SE = ifelse(nPlots > 1, sqrt(sum(clgVar, na.rm = TRUE)) / CARB_100HR * 100,0),
                CARB_1000HR_SE = ifelse(nPlots > 1, sqrt(sum(ccVar, na.rm = TRUE)) / CARB_1000HR * 100,0),
                CARB_PILE_SE = ifelse(nPlots > 1, sqrt(sum(cpVar, na.rm = TRUE)) / CARB_PILE * 100,0),
                CARB_SE = ifelse(nPlots > 1, sqrt(sum(cVar, na.rm = TRUE)) / CARB * 100,0),
                # Per Acre variances
                vsmVar = (1/AREA_TOTAL^2) * (sum(vsmVar, na.rm = TRUE) + (VOL_1HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_1HR_ACRE * sum(vsmCV, na.rm = TRUE))),
                vmdVar = (1/AREA_TOTAL^2) * (sum(vmdVar, na.rm = TRUE) + (VOL_10HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_10HR_ACRE * sum(vmdCV, na.rm = TRUE))),
                vlgVar = (1/AREA_TOTAL^2) * (sum(vlgVar, na.rm = TRUE) + (VOL_100HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_100HR_ACRE * sum(vlgCV, na.rm = TRUE))),
                vcVar = (1/AREA_TOTAL^2) * (sum(vcVar, na.rm = TRUE) + (VOL_1000HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_1000HR_ACRE *sum(vcCV, na.rm = TRUE))),
                vpVar = (1/AREA_TOTAL^2) * (sum(vpVar, na.rm = TRUE) + (VOL_PILE_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_PILE_ACRE * sum(vpCV, na.rm = TRUE))),
                vVar = (1/AREA_TOTAL^2) * (sum(vVar, na.rm = TRUE) + (VOL_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * VOL_ACRE * sum(vCV, na.rm = TRUE))),
                bsmVar = (1/AREA_TOTAL^2) * (sum(bsmVar, na.rm = TRUE) + (BIO_1HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_1HR_ACRE * sum(bsmCV, na.rm = TRUE))),
                bmdVar = (1/AREA_TOTAL^2) * (sum(bmdVar, na.rm = TRUE) + (BIO_10HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_10HR_ACRE * sum(bmdCV, na.rm = TRUE))),
                blgVar = (1/AREA_TOTAL^2) * (sum(blgVar, na.rm = TRUE) + (BIO_100HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_100HR_ACRE * sum(blgCV, na.rm = TRUE))),
                bcVar = (1/AREA_TOTAL^2) * (sum(bcVar, na.rm = TRUE) + (BIO_1000HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_1000HR_ACRE * sum(bcCV, na.rm = TRUE))),
                bpVar = (1/AREA_TOTAL^2) * (sum(bpVar, na.rm = TRUE) + (BIO_PILE_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_PILE_ACRE * sum(bpCV, na.rm = TRUE))),
                bVar = (1/AREA_TOTAL^2) * (sum(bVar, na.rm = TRUE) + (BIO_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * BIO_ACRE * sum(bCV, na.rm = TRUE))),
                csmVar = (1/AREA_TOTAL^2) * (sum(csmVar, na.rm = TRUE) + (CARB_1HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_1HR_ACRE * sum(csmCV, na.rm = TRUE))),
                cmdVar = (1/AREA_TOTAL^2) * (sum(cmdVar, na.rm = TRUE) + (CARB_10HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_10HR_ACRE * sum(cmdCV, na.rm = TRUE))),
                clgVar = (1/AREA_TOTAL^2) * (sum(clgVar, na.rm = TRUE) + (CARB_100HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_100HR_ACRE * sum(clgCV, na.rm = TRUE))),
                ccVar = (1/AREA_TOTAL^2) * (sum(ccVar, na.rm = TRUE) + (CARB_1000HR_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_1000HR_ACRE * sum(ccCV, na.rm = TRUE))),
                cpVar = (1/AREA_TOTAL^2) * (sum(cpVar, na.rm = TRUE) + (CARB_PILE_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_PILE_ACRE * sum(cpCV, na.rm = TRUE))),
                cVar = (1/AREA_TOTAL^2) * (sum(cVar, na.rm = TRUE) + (CARB_ACRE^2 * sum(aVar, na.rm = TRUE) - 2 * CARB_ACRE * sum(cCV, na.rm = TRUE))),
                # Per acre sampling errors
                VOL_1HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(vsmVar, na.rm = TRUE)) / VOL_1HR_ACRE * 100,0),
                VOL_10HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(vmdVar, na.rm = TRUE)) / VOL_10HR_ACRE * 100,0),
                VOL_100HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(vlgVar, na.rm = TRUE)) / VOL_100HR_ACRE * 100,0),
                VOL_1000HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(vcVar, na.rm = TRUE)) / VOL_1000HR_ACRE * 100,0),
                VOL_PILE_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(vpVar, na.rm = TRUE)) / VOL_PILE_ACRE * 100,0),
                VOL_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(vVar, na.rm = TRUE)) / VOL_ACRE * 100,0),
                BIO_1HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(bsmVar, na.rm = TRUE)) / BIO_1HR_ACRE * 100,0),
                BIO_10HR_ACRE_SE = ifelse(nPlots > 1,sqrt(sum(bmdVar, na.rm = TRUE)) / BIO_10HR_ACRE * 100,0),
                BIO_100HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(blgVar, na.rm = TRUE)) / BIO_100HR_ACRE * 100,0),
                BIO_1000HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(bcVar, na.rm = TRUE)) / BIO_1000HR_ACRE * 100,0),
                BIO_PILE_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(bpVar, na.rm = TRUE)) / BIO_PILE_ACRE * 100,0),
                BIO_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(bVar, na.rm = TRUE)) / BIO_ACRE * 100,0),
                CARB_1HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(csmVar, na.rm = TRUE)) / CARB_1HR_ACRE * 100,0),
                CARB_10HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(cmdVar, na.rm = TRUE)) / CARB_10HR_ACRE * 100,0),
                CARB_100HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(clgVar, na.rm = TRUE)) / CARB_100HR_ACRE * 100,0),
                CARB_1000HR_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(ccVar, na.rm = TRUE)) / CARB_1000HR_ACRE * 100,0),
                CARB_PILE_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(cpVar, na.rm = TRUE)) / CARB_PILE_ACRE * 100,0),
                CARB_ACRE_SE = ifelse(nPlots > 1, sqrt(sum(cVar, na.rm = TRUE)) / CARB_ACRE * 100,0))
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
                BIO_1HR = sum(FWD_SM_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_10HR = sum(FWD_MD_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_100HR = sum(FWD_LG_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_1000HR = sum(CWD_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO_PILE = sum(PILE_DRYBIO_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                BIO = sum(BIO_1HR, BIO_10HR, BIO_100HR, BIO_1000HR, BIO_PILE, na.rm = TRUE),
                CARB_1HR = sum(FWD_SM_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_10HR = sum(FWD_MD_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_100HR = sum(FWD_LG_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_1000HR = sum(CWD_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB_PILE = sum(PILE_CARBON_ADJ * aDI * EXPNS / 2000, na.rm = TRUE),
                CARB = sum(CARB_1HR, CARB_10HR, CARB_100HR, CARB_1000HR, CARB_PILE, na.rm = TRUE),
                fa = sum(CONDPROP_UNADJ * aDI * aAdj * EXPNS, na.rm = TRUE),
                plotIn = ifelse(sum(aDI >  0, na.rm = TRUE), 1,0)) %>%
      group_by(.dots = grpBy) %>%
      summarize(AREA_TOTAL = sum(fa, na.rm = TRUE),
                VOL_1HR = sum(VOL_1HR, na.rm = TRUE),
                VOL_10HR = sum(VOL_10HR, na.rm = TRUE),
                VOL_100HR = sum(VOL_100HR, na.rm = TRUE),
                VOL_1000HR = sum(VOL_1000HR, na.rm = TRUE),
                VOL_PILE = sum(VOL_PILE, na.rm = TRUE),
                VOL = sum(VOL, na.rm = TRUE),
                BIO_1HR = sum(BIO_1HR, na.rm = TRUE),
                BIO_10HR = sum(BIO_10HR, na.rm = TRUE),
                BIO_100HR = sum(BIO_100HR, na.rm = TRUE),
                BIO_1000HR = sum(BIO_1000HR, na.rm = TRUE),
                BIO_PILE = sum(BIO_PILE, na.rm = TRUE),
                BIO = sum(BIO, na.rm = TRUE),
                CARB_1HR = sum(CARB_1HR, na.rm = TRUE),
                CARB_10HR = sum(CARB_10HR, na.rm = TRUE),
                CARB_100HR = sum(CARB_100HR, na.rm = TRUE),
                CARB_1000HR = sum(CARB_1000HR, na.rm = TRUE),
                CARB_PILE = sum(CARB_PILE, na.rm = TRUE),
                CARB = sum(CARB, na.rm = TRUE),
                VOL_1HR_ACRE = VOL_1HR / AREA_TOTAL,
                VOL_10HR_ACRE = VOL_10HR / AREA_TOTAL,
                VOL_100HR_ACRE = VOL_100HR / AREA_TOTAL,
                VOL_1000HR_ACRE = VOL_1000HR / AREA_TOTAL,
                VOL_PILE_ACRE = VOL_PILE / AREA_TOTAL,
                VOL_ACRE = VOL / AREA_TOTAL,
                BIO_1HR_ACRE = BIO_1HR / AREA_TOTAL,
                BIO_10HR_ACRE = BIO_10HR / AREA_TOTAL,
                BIO_100HR_ACRE = BIO_100HR / AREA_TOTAL,
                BIO_1000HR_ACRE = BIO_1000HR / AREA_TOTAL,
                BIO_PILE_ACRE = BIO_PILE / AREA_TOTAL,
                BIO_ACRE = BIO / AREA_TOTAL,
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
  gc()

  #Return a dataframe
  cwd


}
