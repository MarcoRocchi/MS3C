sort_times <- function(features, times, responses) {
  sorted_indices <- order(times)
  features <- features[sorted_indices, ]
  times <- times[sorted_indices]
  responses <- !responses[sorted_indices]

  return(list(features = features, times = times, censoring = responses))
}