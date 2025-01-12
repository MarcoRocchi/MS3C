library("mstate")

source("./src/data/features.r")
source("./src/data/fill_dataset.r")

get_optimal_parameters <- function() {
    return(list(eta = 0.1, gamma = 500, mu = 0.001, k = 4))
}

get_tmat <- function() {
    tmat <- transMat(
        x = list(c(2), c(3), c()), 
        names = c("Pre", "Post", "Dead")
    )

    return(tmat)
}

expand_dataset <- function(dataset) {
    cat("\nBuilding long dataset for model 1")

    in_state <- 1

    tmat <- get_tmat()

    data <- cbind(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative, 
        dataset$post_times, 
        dataset$dead_times,
        dataset$dead_status,
        in_state
    )

    data_long <- msprep(
        time = c(NA, names.post_time, names.dead_time),
        status = c(NA, "in_state", names.dead_indicator),
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

    #Post chemo -> Dead
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
    cat("Building model 1: Pre chemo -> Post chemo -> Dead\n")

    tmat <- get_tmat()

    print(tmat)

    print(events(dataset))

    model <- coxph(Surv(time, status) ~
        AAA + AAB + AAC + AAD + AAE + AAF + AAG + AAH + AAI + AAJ + AAK + AAL + AAM + AAN +
        AAO + AAP + AAQ + AAR + AAS + AAT + AAU + AAV + AAW + AAX + AAY + AAZ + ABA +
        BAA + BAB + BAC + BAD + BAE + BAF + BAG + BAH + BAI + BAJ + BAK + BAL + BAM + BAN +
        BAO + BAP + BAQ + BAR + BAS + BAT + BAU + BAV + BAW + BAX + BAY + BAZ + BBA +
        CA + CB + CC + CD + CE + CF + CG + CH +
        CA + CB + CC + CD + CE + CF + CG + CH +
        strata(trans),
        data = dataset
    )
    
    print(cox.zph(model))

    print(summary(model))

    return(model)
}