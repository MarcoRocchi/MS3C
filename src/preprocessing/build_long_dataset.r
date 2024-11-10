source("./src/utils/features.r")

library("mstate")

expand_transitions <- function(dataset) {
    message("Building long dataset")

    tmat <- transMat(
        x = list(c(2), c(3), c(4, 5), c(5), c()), 
        names = c("Pre chemo", "Post chemo", "Post surgery", "Relapse", "Dead")
    )

    in_state <- 1

    data <- cbind(
                dataset$radiomics_pre,
                dataset$radiomics_post,
                dataset$pre_operative, 
                dataset$surgery, 
                dataset$relapse_status, 
                dataset$dead_status, 
                dataset$post_times, 
                dataset$surgery_times,
                dataset$relapse_times, 
                dataset$dead_times,
                in_state)

    data_long <- msprep(
        time = c(NA, names.post_time, names.surgery_time, names.relapse_time, names.dead_time),
        status = c(NA, "in_state", "in_state", names.relapse_indicator, names.dead_indicator),
        data = data,
        trans = tmat,
        keep = c(names.radiomics_pre, names.radiomics_post, names.pre_operative, names.surgery)
    )

    return(data_long)
}