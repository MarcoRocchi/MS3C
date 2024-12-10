#Computes the negative partial log-likelihood of the stratified Cox
neglogparlike <- function(w, data) {
    l <- 0

    for (i in 1:length(data)) {
        d <- data[[i]]
        obsfreq <- as.matrix(d$frequencies * (!d$censoring), nrow = length(d$frequencies), ncol = 1)
        xb <- d$features %*% w[[i]]
        r <- exp(xb)
        risksum <- rev(cumsum(rev(d$frequencies * r)))
        risksum <- risksum[d$atrisk]
        l <- l - as.double(t(obsfreq) %*% (xb - log(risksum)))
    }
    
    return(l)
}