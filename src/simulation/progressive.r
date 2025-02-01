library("mstate")
source("./src/data/fill_dataset.r")

get_tmat <- function() {
    tmat <- transMat(
        x = list(c(2), c(3), c()), 
        names = c("State 1", "State 2", "State 3")
    )

    return(tmat)
}

get_classification_features <- function(dataset) {
    data <- cbind(
        c("t1_1", "t1_2", "t1_3", "t1_4", "t1_5"),
        c("t2_1", "t2_2", "t2_3", "t2_4", "t2_5")
    )

    return(data)
}

split_by_transition <- function(dataset, patients_count) {
    dataset <- group_split(dataset, dataset$trans)

    d <- insert_missing_patients(dataset[[1]], patients_count)
    t1_data <- list(
        features = as.matrix(d[c("t1_1", "t1_2", "t1_3", "t1_4", "t1_5")]),
        new_features = 1:5,
        times = as.matrix(d["time"]),
        status = as.matrix(d["status"]),
        transition_weight = 1
    )

    d <- insert_missing_patients(dataset[[2]], patients_count)
    t2_data <- list(
        features = as.matrix(d[c("t2_1", "t2_2", "t2_3", "t2_4", "t2_5")]),
        new_features = 1:5,
        times = as.matrix(d["time"]),
        status = as.matrix(d["status"]),
        transition_weight = 1
    )

    return(list(t1_data, t2_data))
}

build_model <- function(dataset) {
    cat("Building model simulation progressive\n")

    tmat <- get_tmat()

    print(tmat)

    print(events(dataset))

    data_long <- expand.covs(dataset, 
        c("t1_1", "t1_2", "t1_3", "t1_4", "t1_5", "t1_6", "t1_7", "t1_8", "t1_9", "t1_10",
         "t2_1", "t2_2", "t2_3", "t2_4", "t2_5", "t2_6", "t2_7", "t2_8", "t2_9", "t2_10"), 
            append = TRUE,
            longnames = FALSE)

    model <- coxph(Surv(time, status) ~
        t1_1.1 + t1_2.1 + t1_3.1 + t1_4.1 + t1_5.1 + t1_6.1 + t1_7.1 + t1_8.1 + t1_9.1 + t1_10.1 +
        t2_1.2 + t2_2.2 + t2_3.2 + t2_4.2 + t2_5.2 + t2_6.2 + t2_7.2 + t2_8.2 + t2_9.2 + t2_10.2 +
        strata(trans),
        data = data_long)
    
    print(cox.zph(model))

    print(summary(model))

    return(model)
}

generate_sample_data <- function(n, 
                                time_params = list(mean1, sd1, mean2, sd2), 
                                haz_params = list(mean1, sd1, mean2, sd2), 
                                censor_params = list(p1, p2),
                                covariates_params = list(mean1, sd1, mean2, sd2)) {
    set.seed(123)
    tmat <- get_tmat()
    
    time1 <- sort(abs(rnorm(n, mean = time_params$mean1, sd = time_params$sd1)))
    time2 <- sort(abs(rnorm(n, mean = time_params$mean2, sd = time_params$sd2)))
    Haz <- data.frame(
        time = c(time1, time2),
        Haz = c(
            cumsum(abs(rnorm(n, mean = haz_params$mean1, sd = haz_params$sd1))), 
            cumsum(abs(rnorm(n, mean = haz_params$mean2, sd = haz_params$sd2)))),
        trans = rep(1:2, each = length(time1))
    )
    
    sample_data <- mssample(
        Haz = Haz,
        trans = tmat,
        tvec = time,
        M = n,
        clock = "reset",
        output = "data"
    )
    
    colnames(sample_data)[colnames(sample_data) == "duration"] <- "time"
    
    t1 <- matrix(rnorm(n * 10, mean = covariates_params$mean1, sd = covariates_params$sd1), n, 10)
    t2 <- matrix(rnorm(n * 10, mean = covariates_params$mean2, sd = covariates_params$sd2), n, 10)
    
    t1 <- scale(t1)
    t2 <- scale(t2)

    for (i in 1:n) {
        sample_data[sample_data$id == i, paste0("t1_", 1:10)] <- 
            matrix(rep(t(t1[i, ]), each = sum(sample_data$id == i)), ncol = 10, byrow = FALSE)
        sample_data[sample_data$id == i, paste0("t2_", 1:10)] <- 
            matrix(rep(t(t2[i, ]), each = sum(sample_data$id == i)), ncol = 10, byrow = FALSE)
    }
    
    cens1 <- runif(n, 0, 1) < censor_params$p1
    cens2 <- runif(n, 0, 1) < censor_params$p2
    
    for (i in 1:n) {
        if (cens1[i]) {
            sample_data <- sample_data[!(sample_data$id == i & sample_data$trans == 2), ]
            sample_data[sample_data$id == i, "status"] <- 0
        }
    }

    for (i in 1:n) {
        if (cens2[i]) {
            sample_data[(sample_data$id == i & sample_data$trans == 2), "status"] <- 0
        }
    }

    return(sample_data)
}