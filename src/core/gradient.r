source("./src/core/likelihood.r")

compute_gradient <- function(w, radiomics, freq, cens, atrisk) {
    obsfreq <- freq * (!cens)
    r <- exp(radiomics %*% w)
    risksum <- rev(cumsum(rev(freq * r)))
    risksum <- risksum[atrisk]
    p <- ncol(radiomics)
    xr <- radiomics * kronecker(matrix(1, 1, p), r * freq)
    xrsum <- matrix(rev(cumsum(rev(xr))), nrow = nrow(radiomics))
    xrsum <- xrsum[atrisk, ]
    a <- xrsum / kronecker(matrix(1, 1, p), risksum)
    dl <- -t(t(obsfreq) %*% (radiomics - a))

    return(dl)
}

evaluate_gradient <- function(radiomics, freq, cens, atrisk, w, eta) {
    grad_w <- compute_gradient(w, radiomics, freq, cens, atrisk)
    
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