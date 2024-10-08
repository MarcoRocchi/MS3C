source("./src/preprocessing/standardize.r")
source("./src/preprocessing/sort.r")
source("./src/preprocessing/compute_frequencies.r")
source("./src/preprocessing/compute_atrisk.r")
source("./src/preprocessing/compute_tied.r")

#TODO remove
do_preprocessing <- function(features, times, responses) {
    #features <- standardize_matrix(features)
    #features <- standardize_columns(features)
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
    features <- standardize_matrix(features)
    features <- standardize_columns(features)

    return(features)
}

preprocess_pre_operative <- function(features) {
    age <- features[, 2]
    age <- standardize_columns(age)
    features[, 2] <- age

    mtx <- features[, 4]
    mtx <- standardize_columns(mtx)
    features[, 4] <- mtx
    
    return(features)
}