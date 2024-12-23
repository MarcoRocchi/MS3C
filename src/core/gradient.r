source("./src/core/cox.r")

#TODO nome
compute_log_likelihood <- function(w, data) {
    dl <- c()

    for (i in 1:length(data)) {
        d <- data[[i]]
        #Element i of risksum contains the sum of the expected survival times 
        #of all patients in the risk set of patient i
        r <- exp(d$features %*% w[[i]])
        risksum <- rev(cumsum(rev(d$frequencies * r)))
        risksum <- risksum[d$atrisk]

        p <- ncol(d$features)
        xr <- d$features * kronecker(matrix(1, 1, p), r * d$frequencies)
        revxr <- xr[nrow(xr):1, ]
        xrsum <- apply(revxr, 2, cumsum)
        xrsum <- xrsum[nrow(xrsum):1, ]
        xrsum <- xrsum[d$atrisk, ]
        a <- xrsum / kronecker(matrix(1, 1, p), risksum)
        obsfreq <- d$frequencies * as.numeric((!d$censoring))
        dl[[i]] <- -t(t(obsfreq) %*% (d$features - a))
    }

    return(dl)
}

evaluate_gradient <- function(data, w, eta) {
    gradient <- compute_log_likelihood(w, data)
    likelihood <- neglogparlike(w, data)

    return(list(gradient = gradient, likelihood = likelihood))
}