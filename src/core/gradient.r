source("./src/core/likelihood.r")

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

evaluate_gradient <- function(data, w, eta) {
    grad_w <- compute_log_likelihood(w, data)
    #QUI
    grad_cs <- numeric(ncol(radiomics))

    for (i in seq_along(ncol(radiomics))) {
        for (j in seq_along(ncol(radiomics))) {
            grad_cs <- grad_cs + t(radiomics) %*% (radiomics %*% w - radiomics %*% w)
        }
    }

    grad_w <- grad_w + eta * grad_cs

    func_val <- 0
    func_val <- func_val + neglogparlike(w, radiomics, freq, cens, atrisk)

    return(list(grad_w = grad_w, funcVal = func_val))
}