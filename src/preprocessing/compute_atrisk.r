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