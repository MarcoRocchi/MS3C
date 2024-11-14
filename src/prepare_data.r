sort_times <- function(features, times, responses) {
  sorted_indices <- order(times)
  features <- features[sorted_indices, ]
  times <- times[sorted_indices]
  responses <- !responses[sorted_indices]

  return(list(features = features, times = times, censoring = responses))
}

compute_tied <- function(times) {
  tied <- c(FALSE, diff(times) == 0)
  tied <- c(FALSE, tied) | c(tied, FALSE)

  return(tied)
}

compute_frequencies <- function(features) {
  #TODO
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
    list[features, times, censoring] <- sort_times(features, times, responses)
    frequencies <- compute_frequencies(features)
    atrisk <- compute_atrisk(times)
    tied <- compute_tied(times)

    return(
        list(
            features = features,
            times = times,
            censoring = censoring,
            frequencies = frequencies,
            atrisk = atrisk,
            tied = tied
        )
    )    
}