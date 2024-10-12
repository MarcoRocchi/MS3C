library("mstate")

source("./src/utils/features.r")

build_m6 <- function(pre, post, preop, relapse_status, dead_status, post_times, relapse_times, dead_times) {
    cat("Building model 6: Pre chemo -> Post chemo -> Relapse -> Dead with competing risk\n")

    in_state <- 1

    data <- cbind(pre, preop, relapse_status, dead_status, post_times, relapse_times, dead_times, in_state)

    tmat <- transMat(x = list(c(2), c(3, 4), c(4), c()), names = c("Pre", "Post", "Relapse", "Dead"))

    print(tmat)

    data_long <- msprep(time = c(NA, names.post_time, names.relapse_time, names.dead_time),
                        status = c(NA, "in_state", names.relapse_indicator, names.dead_indicator),
                        data = data,
                        trans = tmat,
                        keep = c(names.radiomics, names.pre_operative))

    data_long[data_long$trans == 2, names.radiomics] <- post
    data_long[data_long$trans == 3, names.radiomics] <- post


    data_long <- expand.covs(data_long, 
                            c(names.radiomics, names.pre_operative), 
                            append = TRUE,
                            longnames = FALSE)

    print(events(data_long))

    model <- coxph(Surv(Tstart, Tstop, status) ~
        AA.1 + AB.1 + AC.1 + AD.1 + AE.1 + AF.1 + AG.1 + AH.1 + AI.1 + AJ.1 + AK.1 + AL.1 + AM.1 + AN.1 +
        AO.1 + AP.1 + AQ.1 + AR.1 + AS.1 + AT.1 + AU.1 + AV.1 + AW.1 + AX.1 + AY.1 + AZ.1 + BA.1 +
        AA.2 + AB.2 + AC.2 + AD.2 + AE.2 + AF.2 + AG.2 + AH.2 + AI.2 + AJ.2 + AK.2 + AL.2 + AM.2 + AN.2 +
        AO.2 + AP.2 + AQ.2 + AR.2 + AS.2 + AT.2 + AU.2 + AV.2 + AW.2 + AX.2 + AY.2 + AZ.2 + BA.2 +
        AA.3 + AB.3 + AC.3 + AD.3 + AE.3 + AF.3 + AG.3 + AH.3 + AI.3 + AJ.3 + AK.3 + AL.3 + AM.3 + AN.3 +
        AO.3 + AP.3 + AQ.3 + AR.3 + AS.3 + AT.3 + AU.3 + AV.3 + AW.3 + AX.3 + AY.3 + AZ.3 + BA.3 +
        CA.1 + CB.1 + CC.1 + CD.1 + CE.1 + CF.1 + CG.1 + CH.1 +
        CA.2 + CB.2 + CC.2 + CD.2 + CE.2 + CF.2 + CG.2 + CH.2 +
        CA.3 + CB.3 + CC.3 + CD.3 + CE.3 + CF.3 + CG.3 + CH.3 +
        CA.4 + CB.4 + CC.4 + CD.4 + CE.4 + CF.4 + CG.4 + CH.4 +
        strata(trans),
        data = data_long)

    print(cox.zph(model))

    print(summary(model))

    return(model)
}