library("mstate")

source("./src/data/features.r")
source("./src/data/fill_dataset.r")

get_optimal_parameters <- function() {
    return(list(eta = 0.1, gamma = 500, mu = 1e-4, k = 6))
}

get_tmat <- function() {
    tmat <- transMat(
        x = list(c(2), c(3), c()), 
        names = c("Pre", "Post", "Relapse")
    )

    return(tmat)
}

get_classification_features <- function(dataset) {
    data <- cbind(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative
    )

    return(data)
}

expand_dataset <- function(dataset) {
    cat("\nBuilding long dataset for model 2")

    in_state <- 1

    tmat <- get_tmat()

    data <- cbind(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative, 
        dataset$post_times, 
        dataset$relapse_times,
        dataset$relapse_status,
        in_state
    )

    data_long <- msprep(
        time = c(NA, names.post_time, names.relapse_time),
        status = c(NA, "in_state", names.relapse_indicator),
        data = data,
        trans = tmat,
        keep = c(names.radiomics_pre, names.radiomics_post, names.pre_operative)
    )

    return(data_long)
}

split_by_transition <- function(dataset, patients_count) {
    dataset <- group_split(dataset, dataset$trans)

    #Pre chemo -> Post chemo
    d <- insert_missing_patients(dataset[[1]], patients_count)
    t1_data <- list(
        features = as.matrix(d[c(names.radiomics_pre, names.pre_operative)]),
        new_features = 1:(length(names.radiomics_pre) + length(names.pre_operative)),
        times = as.matrix(d["time"]),
        status = as.matrix(d["status"]),
        transition_weight = 1
    )

    #Post chemo -> Relapse
    d <- insert_missing_patients(dataset[[2]], patients_count)
    t2_data <- list(
        features = as.matrix(d[c(names.radiomics_post, names.pre_operative)]),
        new_features = 1:length(names.radiomics_post),        
        times = as.matrix(d["time"]),
        status = as.matrix(d["status"]),
        transition_weight = 1
    )

    return(list(t1_data, t2_data))
}

build_model <- function(dataset) {
    cat("Building model 2: Pre chemo -> Post chemo -> Relapse\n")

    tmat <- get_tmat()

    print(tmat)

    print(events(dataset))

    data_long <- expand.covs(dataset, 
        c(names.radiomics_pre, names.radiomics_post, names.pre_operative), 
            append = TRUE,
            longnames = FALSE)

    model <- coxph(Surv(time, status) ~
        AAA.1 + AAB.1 + AAC.1 + AAD.1 + AAE.1 + AAF.1 + AAG.1 + AAH.1 + AAI.1 + AAJ.1 + AAK.1 + AAL.1 + AAM.1 + AAN.1 +
        AAO.1 + AAP.1 + AAQ.1 + AAR.1 + AAS.1 + AAT.1 + AAU.1 + AAV.1 + AAW.1 + AAX.1 + AAY.1 + AAZ.1 + ABA.1 +
        BAA.2 + BAB.2 + BAC.2 + BAD.2 + BAE.2 + BAF.2 + BAG.2 + BAH.2 + BAI.2 + BAJ.2 + BAK.2 + BAL.2 + BAM.2 + BAN.2 +
        BAO.2 + BAP.2 + BAQ.2 + BAR.2 + BAS.2 + BAT.2 + BAU.2 + BAV.2 + BAW.2 + BAX.2 + BAY.2 + BAZ.2 + BBA.2 +
        CA.1 + CB.1 + CC.1 + CD.1 + CE.1 + CF.1 + CG.1 + CH.1 +
        CA.2 + CB.2 + CC.2 + CD.2 + CE.2 + CF.2 + CG.2 + CH.2 +
        strata(trans),
        data = data_long)

    print(cox.zph(model))

    print(summary(model))

    return(model)
}