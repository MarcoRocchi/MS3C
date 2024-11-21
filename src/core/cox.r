#Computes the negative partial log-likelihood of the multiview Cox
neglogparlike <- function(w, data) {
    l <- 0

    start <- 1
    stop <- 0

    for (d in data) {
        stop <- stop + ncol(d$features)
        obsfreq <- as.matrix(d$frequencies * (!d$censoring), nrow = length(d$frequencies), ncol = 1)
        xb <- d$features %*% w[start:stop]
        r <- exp(xb)
        risksum <- rev(cumsum(rev(d$frequencies * r)))
        risksum <- risksum[d$atrisk]
        l <- l - (t(obsfreq) %*% (xb - log(risksum)))
        start <- start + ncol(d$features)
    }
    
    return(l)
}