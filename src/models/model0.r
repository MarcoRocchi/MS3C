library("mstate")

source("./src/data/features.r")
source("./src/data/fill_dataset.r")

get_optimal_parameters <- function() {
    return(list(eta = 0.0005, tau = 2, k = 3))
}

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
        transition_weight = 1
    )

    return(list(t1_data))
}

build_model <- function(dataset) {
    cat("Building model 0: Pre chemo -> Relapse\n")

    print(events(dataset))

    model <- coxph(Surv(Tstart, Tstop, status) ~
        AAA + AAB + AAC + AAD + AAE + AAF + AAG + AAH + AAI + AAJ + AAK + AAL + AAM + 
        AAN + AAO + AAP + AAQ + AAR + AAS + AAT + AAU + AAV + AAW + AAX + AAY + AAZ + ABA,
        data = dataset)

    print(cox.zph(model))

    print(summary(model))

    return(model)
}