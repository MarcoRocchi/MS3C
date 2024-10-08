library("mstate")

source("./src/utils/features.r")

build_m5 <- function(pre, 
                    post, 
                    preop, 
                    surgery,
                    relapse_status, 
                    dead_status, 
                    post_times, 
                    surgery_times, 
                    relapse_times, 
                    dead_times) {
    
    cat("Building model 5: Pre chemo -> Post chemo -> Surgery -> Relapse -> Dead\n")

    in_state <- 1

    data <- cbind(
                pre,
                post,
                preop, 
                surgery, 
                relapse_status, 
                dead_status, 
                post_times, 
                surgery_times,
                relapse_times, 
                dead_times,
                in_state)

    tmat <- transMat(x = list(c(2), c(3), c(4), c(5), c()), names = c("Pre", "Post", "Surgery", "Relapse", "Dead"))

    print(tmat)

    data_long <- msprep(time = c(NA, names.post_time, names.surgery_time, names.relapse_time, names.dead_time),
                    status = c(NA, "in_state", "in_state", names.relapse_indicator, names.dead_indicator),
                    data = data,
                    trans = tmat,
                    keep = c(
                            paste0("A", names.radiomics), 
                            paste0("B", names.radiomics), 
                            names.pre_operative, 
                            names.surgery)
                        )

    data_long <- expand.covs(data_long, 
                            c(paste0("A", names.radiomics), paste0("B", names.radiomics), names.surgery), 
                            append = TRUE,
                            longnames = FALSE)

    print(events(data_long))

    model <- coxph(Surv(Tstart, Tstop, status) ~
        AAA.1 + AAB.1 + AAC.1 + AAD.1 + AAE.1 + AAF.1 + AAG.1 + AAH.1 + AAI.1 + AAJ.1 + AAK.1 + AAL.1 + AAM.1 + AAN.1 +
        AAO.1 + AAP.1 + AAQ.1 + AAR.1 + AAS.1 + AAT.1 + AAU.1 + AAV.1 + AAW.1 + AAX.1 + AAY.1 + AAZ.1 + ABA.1 +
        AAA.2 + AAB.2 + AAC.2 + AAD.2 + AAE.2 + AAF.2 + AAG.2 + AAH.2 + AAI.2 + AAJ.2 + AAK.2 + AAL.2 + AAM.2 + AAN.2 +
        AAO.2 + AAP.2 + AAQ.2 + AAR.2 + AAS.2 + AAT.2 + AAU.2 + AAV.2 + AAW.2 + AAX.2 + AAY.2 + AAZ.2 + ABA.2 +
        AAA.3 + AAB.3 + AAC.3 + AAD.3 + AAE.3 + AAF.3 + AAG.3 + AAH.3 + AAI.3 + AAJ.3 + AAK.3 + AAL.3 + AAM.3 + AAN.3 +
        AAO.3 + AAP.3 + AAQ.3 + AAR.3 + AAS.3 + AAT.3 + AAU.3 + AAV.3 + AAW.3 + AAX.3 + AAY.3 + AAZ.3 + ABA.3 +
        AAA.4 + AAB.4 + AAC.4 + AAD.4 + AAE.4 + AAF.4 + AAG.4 + AAH.4 + AAI.4 + AAJ.4 + AAK.4 + AAL.4 + AAM.4 + AAN.4 +
        AAO.4 + AAP.4 + AAQ.4 + AAR.4 + AAS.4 + AAT.4 + AAU.4 + AAV.4 + AAW.4 + AAX.4 + AAY.4 + AAZ.4 + ABA.4 +
        BAA.1 + BAB.1 + BAC.1 + BAD.1 + BAE.1 + BAF.1 + BAG.1 + BAH.1 + BAI.1 + BAJ.1 + BAK.1 + BAL.1 + BAM.1 + BAN.1 +
        BAO.1 + BAP.1 + BAQ.1 + BAR.1 + BAS.1 + BAT.1 + BAU.1 + BAV.1 + BAW.1 + BAX.1 + BAY.1 + BAZ.1 + BBA.1 +
        BAA.2 + BAB.2 + BAC.2 + BAD.2 + BAE.2 + BAF.2 + BAG.2 + BAH.2 + BAI.2 + BAJ.2 + BAK.2 + BAL.2 + BAM.2 + BAN.2 +
        BAO.2 + BAP.2 + BAQ.2 + BAR.2 + BAS.2 + BAT.2 + BAU.2 + BAV.2 + BAW.2 + BAX.2 + BAY.2 + BAZ.2 + BBA.2 +
        BAA.3 + BAB.3 + BAC.3 + BAD.3 + BAE.3 + BAF.3 + BAG.3 + BAH.3 + BAI.3 + BAJ.3 + BAK.3 + BAL.3 + BAM.3 + BAN.3 +
        BAO.3 + BAP.3 + BAQ.3 + BAR.3 + BAS.3 + BAT.3 + BAU.3 + BAV.3 + BAW.3 + BAX.3 + BAY.3 + BAZ.3 + BBA.3 +
        BAA.4 + BAB.4 + BAC.4 + BAD.4 + BAE.4 + BAF.4 + BAG.4 + BAH.4 + BAI.4 + BAJ.4 + BAK.4 + BAL.4 + BAM.4 + BAN.4 +
        BAO.4 + BAP.4 + BAQ.4 + BAR.4 + BAS.4 + BAT.4 + BAU.4 + BAV.4 + BAW.4 + BAX.4 + BAY.4 + BAZ.4 + BBA.4 +
        CA + CB + CC + CD + CE + 
        DA.1 + DB.1 + DC.1 + DA.2 + DB.2 + DC.2 + DA.3 + DB.3 + DC.3 + DA.4 + DB.4 + DC.4 +
        strata(trans), 
        data = data_long)

    print(summary(model))

    return(model)
}