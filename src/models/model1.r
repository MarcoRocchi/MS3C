library("mstate")

source("./src/utils/features.r")

build_m1 <- function(pre, post, preop, dead_status, post_times, dead_times) {
    cat("Building model 1: Pre chemo -> Post chemo -> Dead\n")

    in_state <- 1

    data <- cbind(pre, post, preop, dead_status, post_times, dead_times, in_state)

    tmat <- transMat(x = list(c(2), c(3), c()), names = c("Pre", "Post", "Dead"))

    print(tmat)

    data_long <- msprep(time = c(NA, names.post_time, names.dead_time),
                        status = c(NA, "in_state", names.dead_indicator),
                        data = data,
                        trans = tmat,
                        keep = c(paste0("A", names.radiomics), paste0("B", names.radiomics), names.pre_operative))

    data_long <- expand.covs(data_long, 
                            c(paste0("A", names.radiomics), paste0("B", names.radiomics)), 
                            append = TRUE,
                            longnames = FALSE)

    print(events(data_long))

    model <- coxph(Surv(Tstart, Tstop, status) ~
        AAA.1 + AAB.1 + AAC.1 + AAD.1 + AAE.1 + AAF.1 + AAG.1 + AAH.1 + AAI.1 + AAJ.1 + AAK.1 + AAL.1 + AAM.1 + AAN.1 +
        AAO.1 + AAP.1 + AAQ.1 + AAR.1 + AAS.1 + AAT.1 + AAU.1 + AAV.1 + AAW.1 + AAX.1 + AAY.1 + AAZ.1 + ABA.1 +
        AAA.2 + AAB.2 + AAC.2 + AAD.2 + AAE.2 + AAF.2 + AAG.2 + AAH.2 + AAI.2 + AAJ.2 + AAK.2 + AAL.2 + AAM.2 + AAN.2 +
        AAO.2 + AAP.2 + AAQ.2 + AAR.2 + AAS.2 + AAT.2 + AAU.2 + AAV.2 + AAW.2 + AAX.2 + AAY.2 + AAZ.2 + ABA.2 + 
        BAA.1 + BAB.1 + BAC.1 + BAD.1 + BAE.1 + BAF.1 + BAG.1 + BAH.1 + BAI.1 + BAJ.1 + BAK.1 + BAL.1 + BAM.1 + BAN.1 +
        BAO.1 + BAP.1 + BAQ.1 + BAR.1 + BAS.1 + BAT.1 + BAU.1 + BAV.1 + BAW.1 + BAX.1 + BAY.1 + BAZ.1 + BBA.1 +
        BAA.2 + BAB.2 + BAC.2 + BAD.2 + BAE.2 + BAF.2 + BAG.2 + BAH.2 + BAI.2 + BAJ.2 + BAK.2 + BAL.2 + BAM.2 + BAN.2 +
        BAO.2 + BAP.2 + BAQ.2 + BAR.2 + BAS.2 + BAT.2 + BAU.2 + BAV.2 + BAW.2 + BAX.2 + BAY.2 + BAZ.2 + BBA.2 +
        CA + CB + CC + CD + CE + 
        strata(trans), 
        data = data_long)

    print(summary(model))

    return(model)
}