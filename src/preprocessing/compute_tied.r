compute_tied <- function(times) {
  tied <- c(FALSE, diff(times) == 0)
  tied <- c(FALSE, tied) | c(tied, FALSE)

  return(tied)
}