standardize_matrix <- function(features) {
  m <- mean(features)
  stdev <- sd(features)
  features <- (features - m) / stdev

  return(features)
}

standardize_columns <- function(features) {
  features <- scale(features)

  return(features)
}