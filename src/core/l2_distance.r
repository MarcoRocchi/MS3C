l2_distance <- function(a, b) {
    if (nrow(a) == 1) {
        a <- rbind(a, rep(0, ncol(a)))
        b <- rbind(b, rep(0, ncol(b)))
    }

    aa <- colSums(a * a)
    bb <- colSums(b * b)
    ab <- crossprod(a, b)
    
    d <- kronecker(matrix(1, 1, length(bb)), aa) + kronecker(matrix(1, length(aa), 1), t(bb)) - 2 * ab

    d <- Re(d)
    d <- pmax(d, 0)

    return(d)
}