#responses: a boolean array of the same size as Y that is 1 for
# observations that are right-censored and 0 for
# observations that are observed exactly
sort_times <- function(features, times, responses) {
    sorted_indices <- order(times)
    features <- features[sorted_indices, ]
    times <- times[sorted_indices]
    responses <- ifelse(responses[sorted_indices] == 1, 0, 1)

    return(list(features = features, times = times, censoring = responses, sorted_indices = sorted_indices))
}

compute_tied <- function(times) {
    tied <- c(FALSE, diff(times) == 0)
    tied <- c(FALSE, tied) | c(tied, FALSE)

    return(tied)
}

#An array of the same size as Y containing non-negative integer counts.
#The jth element of this vector gives the number of times the jth element
#of Y and the jth row of X were observed.
#Default is 1 observation per row of X and Y
compute_frequencies <- function(features) {
    freq <- rep(1, nrow(features))

    return(freq)
}

compute_atrisk <- function(times) {
    x <- times
    y <- rev(times)

    atrisk <- rep(0, length(x))

    for (i in 1:length(x)) {
        last <- 1
        for (j in 1:length(y)) {
            if (y[j] == x[i]) {
                last <- j
            }
        }
        atrisk[i] <- last
    }

    atrisk <- length(times) + 1 - atrisk

    return(atrisk)
}

prepare_data <- function(features, times, responses) {
    list[features, times, censoring, sorted_indices] <- sort_times(features, times, responses)
    frequencies <- compute_frequencies(features)
    atrisk <- compute_atrisk(times)
    tied <- compute_tied(times)

    return(
        list(
            features = features,
            patients = sorted_indices,
            times = times,
            censoring = censoring,
            frequencies = frequencies,
            atrisk = atrisk,
            tied = tied
        )
    )    
}