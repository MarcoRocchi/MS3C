library("mstate")

source("./src/utils/features.r")

build_m4 <- function(pre, post, preop, surgery, relapse_status, post_times, surgery_times, relapse_times) {
    cat("Building model 4: Pre chemo -> Post chemo -> Surgery -> Relapse\n")

    in_state <- 1

    data <- cbind(pre, post, preop, surgery, relapse_status, post_times, surgery_times, relapse_times, in_state)

    tmat <- transMat(x = list(c(2), c(3), c(4), c()), names = c("Pre", "Post", "Surgery", "Relapse"))

    print(tmat)

    data_long <- msprep(time = c(NA, names.post_time, names.surgery_time, names.relapse_time),
                        status = c(NA, "in_state", "in_state", names.relapse_indicator),
                        data = data,
                        trans = tmat,
                        keep = c(names.radiomics, names.pre_operative, names.surgery)
                    )

    data_long[data_long$trans == 2, names.radiomics] <- post

    data_long <- expand.covs(data_long, 
                            c(names.radiomics, names.surgery), 
                            append = TRUE,
                            longnames = FALSE)

    print(events(data_long))

    model <- coxph(Surv(Tstart, Tstop, status) ~
        AA.1 + AB.1 + AC.1 + AD.1 + AE.1 + AF.1 + AG.1 + AH.1 + AI.1 + AJ.1 + AK.1 + AL.1 + AM.1 + AN.1 +
        AO.1 + AP.1 + AQ.1 + AR.1 + AS.1 + AT.1 + AU.1 + AV.1 + AW.1 + AX.1 + AY.1 + AZ.1 + BA.1 +
        AA.2 + AB.2 + AC.2 + AD.2 + AE.2 + AF.2 + AG.2 + AH.2 + AI.2 + AJ.2 + AK.2 + AL.2 + AM.2 + AN.2 +
        AO.2 + AP.2 + AQ.2 + AR.2 + AS.2 + AT.2 + AU.2 + AV.2 + AW.2 + AX.2 + AY.2 + AZ.2 + BA.2 +
        CA + CB + CC + CD + CE + 
        DA.3 + DB.3 + DC.3,
        data = data_long)

    print(summary(model))

    return(model)
}