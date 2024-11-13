source("./src/utils/features.r")

library("mstate")

expand_transitions <- function(dataset) {
    message("Building long dataset")

    tmat <- transMat(
        x = list(c(2), c(3), c()), 
        names = c("Pre chemo", "Post chemo", "Relapse")
    )

    in_state <- 1

    data <- cbind(
                dataset$radiomics_pre,
                dataset$radiomics_post,
                dataset$pre_operative, 
                dataset$relapse_status, 
                dataset$post_times, 
                dataset$relapse_times, 
                in_state)

    data_long <- msprep(
        time = c(NA, names.post_time, names.relapse_time),
        status = c(NA, "in_state", names.relapse_indicator),
        data = data,
        trans = tmat,
        keep = c(names.radiomics_pre, names.radiomics_post, names.pre_operative)
    )

    return(data_long)
}