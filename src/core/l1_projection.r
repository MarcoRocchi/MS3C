l1_projection <- function(v, beta) {
    z <- matrix(0, nrow = nrow(v), ncol = ncol(v))
    vp <- v - beta / 2
    z[v > beta / 2] <- vp[v > beta / 2]
    vn <- v + beta / 2
    z[v < -beta / 2] <- vn[v < -beta / 2]
    l1_comp_val <- sum(abs(z))
    return(list(z = z, l1_comp_val = l1_comp_val))
}