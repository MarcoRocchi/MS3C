source("./src/core/likelihood.r")

#TODO nome
compute_log_likelihood <- function(w, data) {
    start <- 1
    stop <- 0
    dl <- c()

    for (d in data) {
        stop <- stop + ncol(d$features)
        #TODO controllare flusso variabile censoring
        obsfreq <- d$frequencies * (!d$censoring)
        r <- exp(d$features %*% w[start:stop])
        risksum <- rev(cumsum(rev(d$frequencies * r)))
        risksum <- risksum[d$atrisk]
        p <- ncol(d$features)
        xr <- d$features * kronecker(matrix(1, 1, p), r * d$frequencies)
        xrsum <- matrix(rev(cumsum(rev(xr))), nrow = nrow(d$features))
        xrsum <- xrsum[d$atrisk, ]
        a <- xrsum / kronecker(matrix(1, 1, p), risksum)
        dl <- append(dl, -t(t(obsfreq) %*% (d$features - a)))
        start <- start + ncol(d$features)
    }

    return(dl)
}

evaluate_gradient <- function(data, w, eta, features_count) {
    grad_w <- compute_log_likelihood(w, data)
    grad_cs <- numeric(features_count)

    start_i <- 1
    stop_i <- 0

    for (d1 in data) {
        stop_i <- stop_i + ncol(d1$features)
        start_j <- 1
        stop_j <- 0
        
        for (d2 in data) {
            stop_j <- stop_j + ncol(d2$features)
            grad_cs[start_i:stop_i] <- 
                grad_cs[start_i:stop_i] + 
                t(d1$features) %*% (d1$features %*% w[start_i:stop_i] - d2$features %*% w[start_j:stop_j])
            start_j <- start_j + ncol(d2$features)
        }

        start_i <- start_i + ncol(d1$features)
    }

    grad_w <- grad_w + eta * grad_cs

    func_val <- 0
    func_val <- func_val + neglogparlike(w, data)

    return(list(grad_w = grad_w, funcVal = func_val))
}