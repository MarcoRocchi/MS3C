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

#Computes the gradient of the stratified cox
compute_gradient <- function(w, data) {
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

evaluate_gradient <- function(data, w) {
    gradient <- compute_gradient(w, data)
    likelihood <- neglogparlike(w, data)

    return(list(gradient = gradient, likelihood = likelihood))
}

compute_variance <- function(wzp, X, time, status) {
    if (all(wzp == 0)) {
        return(1)
    }

    n <- nrow(X)
    p <- ncol(X)
    eta <- as.vector(X %*% wzp)
    eta <- eta - max(eta)
    exp_eta <- exp(eta)
    info <- matrix(0, p, p)
    
    for (i in 1:n) {
        if (status[i] == 1) {
            risk_set <- which(time >= time[i])
            S0 <- sum(exp_eta[risk_set])
            S1 <- colSums(X[risk_set, , drop = FALSE] * exp_eta[risk_set])
            S2 <- t(X[risk_set, , drop = FALSE]) %*% (X[risk_set, , drop = FALSE] * exp_eta[risk_set])
            info <- info + (S2 / S0) - (S1 %*% t(S1)) / (S0^2)
        }
    }
    info <- info + diag(ncol(info)) 
    inv_info <- solve(info)    
    variances <- diag(inv_info)    
    avg_variance <- mean(variances)
    return(avg_variance)
}