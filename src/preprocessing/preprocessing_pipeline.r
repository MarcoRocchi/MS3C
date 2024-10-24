source("./src/preprocessing/standardize.r")
source("./src/preprocessing/sort.r")
source("./src/preprocessing/compute_frequencies.r")
source("./src/preprocessing/compute_atrisk.r")
source("./src/preprocessing/compute_tied.r")

do_preprocessing <- function(features, times, responses) {
    features <- standardize_columns(as.matrix(features))
    list[features, times, responses] <- sort_times(features, times, responses)
    frequencies <- compute_frequencies(features)
    atrisk <- compute_atrisk(times)
    tied <- compute_tied(times)

    return(list(
            features = features,
            times = times,
            responses = responses,
            frequencies = frequencies,
            atrisk = atrisk,
            tied = tied
        )
    )
}

preprocess <- function(features) {
    features <- standardize_columns(features)

    return(features)
}

preprocess_pre_operative <- function(features) {
    age <- features[, 1]
    age <- standardize_columns(age)
    features[, 1] <- age

    mtx <- features[, 2]
    mtx <- standardize_columns(mtx)
    features[, 2] <- mtx

    return(features)
}