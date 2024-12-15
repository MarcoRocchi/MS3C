library("mstate")

source("./src/data/features.r")
source("./src/data/fill_dataset.r")

get_tmat <- function() {
    tmat <- transMat(
        x = list(c(2), c()),
        names = c("Pre chemo", "Relapse")
    )

    return(tmat)
}

expand_dataset <- function(dataset) {
    cat("\nBuilding long dataset for model 0")

    tmat <- get_tmat()

    data <- cbind(
        dataset$radiomics_pre,
        dataset$relapse_times,
        dataset$relapse_status 
    )

    data_long <- msprep(
        time = c(NA, names.relapse_time),
        status = c(NA, names.relapse_indicator),
        data = data,
        trans = tmat,
        keep = c(names.radiomics_pre)
    )

    return(data_long)
}

split_by_transition <- function(dataset, patients_count) {
    dataset <- group_split(dataset, dataset$trans)

    #Pre chemo -> Post chemo
    d <- insert_missing_patients(dataset[[1]], patients_count)
    t1_data <- list(
        features = as.matrix(d[names.radiomics_pre]),
        new_features = 1:(length(names.radiomics_pre)),
        times = as.matrix(d["time"]),
        status = as.matrix(d["status"]),
        #TODO
        transition_coefficient = 1
    )

    return(list(t1_data))
}

build_m0 <- function(dataset) {
    cat("Building model 0: Pre chemo -> Relapse\n")

    data <- expand_dataset(dataset)

    data_long <- expand.covs(
        data, 
        c(names.radiomics_pre), 
        append = TRUE,
        longnames = FALSE
    )

    print(events(data_long))

    #Rimosso BAH
    model <- coxph(Surv(Tstart, Tstop, status) ~
        AAA.1 + AAB.1 + AAC.1 + AAD.1 + AAE.1 + AAF.1 + AAG.1 + AAH.1 + AAI.1 + AAJ.1 + AAK.1 + AAL.1 + AAM.1 + AAN.1 +
        AAO.1 + AAP.1 + AAQ.1 + AAR.1 + AAS.1 + AAT.1 + AAU.1 + AAV.1 + AAW.1 + AAX.1 + AAY.1 + AAZ.1 + ABA.1 +
        BAA.2 + BAB.2 + BAC.2 + BAD.2 + BAE.2 + BAF.2 + BAG.2 + BAI.2 + BAJ.2 + BAK.2 + BAL.2 + BAM.2 + BAN.2 +
        BAO.2 + BAP.2 + BAQ.2 + BAR.2 + BAS.2 + BAT.2 + BAU.2 + BAV.2 + BAW.2 + BAX.2 + BAY.2 + BAZ.2 + BBA.2 +
        BAA.3 + BAB.3 + BAC.3 + BAD.3 + BAE.3 + BAF.3 + BAG.3 + BAI.3 + BAJ.3 + BAK.3 + BAL.3 + BAM.3 + BAN.3 +
        BAO.3 + BAP.3 + BAQ.3 + BAR.3 + BAS.3 + BAT.3 + BAU.3 + BAV.3 + BAW.3 + BAX.3 + BAY.3 + BAZ.3 + BBA.3 +
        BAA.4 + BAB.4 + BAC.4 + BAD.4 + BAE.4 + BAF.4 + BAG.4 + BAI.4 + BAJ.4 + BAK.4 + BAL.4 + BAM.4 + BAN.4 +
        BAO.4 + BAP.4 + BAQ.4 + BAR.4 + BAS.4 + BAT.4 + BAU.4 + BAV.4 + BAW.4 + BAX.4 + BAY.4 + BAZ.4 + BBA.4 +
        CA.1 + CB.1 + CC.1 + CD.1 + CE.1 + CF.1 + CG.1 + CH.1 +
        CA.2 + CB.2 + CC.2 + CD.2 + CE.2 + CF.2 + CG.2 + CH.2 +
        CA.3 + CB.3 + CC.3 + CD.3 + CE.3 + CF.3 + CG.3 + CH.3 +
        CA.4 + CB.4 + CC.4 + CD.4 + CE.4 + CF.4 + CG.4 + CH.4 +
        strata(trans),
        data = data_long)

    #print(cox.zph(model))

    print(summary(model))

    return(model)
}