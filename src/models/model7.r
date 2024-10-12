library("mstate")

source("./src/utils/features.r")

build_m7 <- function(pre, 
                    post, 
                    preop, 
                    surgery,
                    relapse_status, 
                    dead_status, 
                    post_times, 
                    surgery_times, 
                    relapse_times, 
                    dead_times) {
    
    cat("Building model 7: Pre chemo -> Post chemo -> Surgery -> Relapse -> Dead with competing risk\n")

    in_state <- 1

    data <- cbind(
                pre,
                preop, 
                surgery, 
                relapse_status, 
                dead_status, 
                post_times, 
                surgery_times,
                relapse_times, 
                dead_times,
                in_state)

    tmat <- transMat(x = list(c(2), c(3), c(4, 5), c(5), c()), names = c("Pre", "Post", "Surgery", "Relapse", "Dead"))

    print(tmat)

    data_long <- msprep(time = c(NA, names.post_time, names.surgery_time, names.relapse_time, names.dead_time),
                    status = c(NA, "in_state", "in_state", names.relapse_indicator, names.dead_indicator),
                    data = data,
                    trans = tmat,
                    keep = c(names.radiomics, names.pre_operative, names.surgery)
                )

    data_long[data_long$trans == 2, names.radiomics] <- post

    data_long <- expand.covs(data_long, 
                            c(names.radiomics, names.pre_operative, names.surgery), 
                            append = TRUE,
                            longnames = FALSE)

    print(events(data_long))

    model <- coxph(Surv(Tstart, Tstop, status) ~
        AA.1 + AB.1 + AC.1 + AD.1 + AE.1 + AF.1 + AG.1 + AH.1 + AI.1 + AJ.1 + AK.1 + AL.1 + AM.1 + AN.1 +
        AO.1 + AP.1 + AQ.1 + AR.1 + AS.1 + AT.1 + AU.1 + AV.1 + AW.1 + AX.1 + AY.1 + AZ.1 + BA.1 +
        AA.2 + AB.2 + AC.2 + AD.2 + AE.2 + AF.2 + AG.2 + AH.2 + AI.2 + AJ.2 + AK.2 + AL.2 + AM.2 + AN.2 +
        AO.2 + AP.2 + AQ.2 + AR.2 + AS.2 + AT.2 + AU.2 + AV.2 + AW.2 + AX.2 + AY.2 + AZ.2 + BA.2 +
        CA.1 + CB.1 + CC.1 + CD.1 + CE.1 + CF.1 + CG.1 + CH.1 +
        CA.2 + CB.2 + CC.2 + CD.2 + CE.2 + CF.2 + CG.2 + CH.2 +
        CA.3 + CB.3 + CC.3 + CD.3 + CE.3 + CF.3 + CG.3 + CH.3 +
        CA.4 + CB.4 + CC.4 + CD.4 + CE.4 + CF.4 + CG.4 + CH.4 +
        CA.5 + CB.5 + CC.5 + CD.5 + CE.5 + CF.5 + CG.5 + CH.5 +
        DA.3 + DB.3 + DC.3 + DD.3 + DE.3 + DF.3 + 
        DA.4 + DB.4 + DC.4 + DD.4 + DE.4 + DF.4 +
        DA.5 + DB.5 + DC.5 + DD.5 + DE.5 + DF.5 +
        strata(trans),
        data = data_long)
    
    print(cox.zph(model))

    print(summary(model))

    return(model)
}