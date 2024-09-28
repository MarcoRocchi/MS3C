neglogparlike <- function(w, features, freq, cens, atrisk) {
    obsfreq <- as.matrix(freq * (!cens), nrow=length(freq), ncol=1)
    xb <- features %*% w
    r <- exp(xb)
    risksum <- rev(cumsum(rev(freq * r)))
    risksum <- risksum[atrisk]
    l <- -(t(obsfreq) %*% (xb - log(risksum)))

    return(l)
}