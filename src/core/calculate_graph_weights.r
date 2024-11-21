#Estimate S by fixing w
CalGraphWeight <- function(features, non_repeated_features, w, k) {
    #TODO ogni transizone ha il suo w
    mu <- 10e-4
    alpha <- 1
    n <- nrow(features)

    distX <- l2_distance(t(non_repeated_features), t(non_repeated_features))

    y <- features %*% w
    disty <- l2_distance(t(y), t(y))
    S <- matrix(0, n, n)
    
    idx <- t(apply(mu * distX + disty, 2, order))
        
    for (i in 1:n) {
        idxa0 <- idx[i, 2:(k + 1)]
        dfi <- disty[i, idxa0]
        dxi <- distX[i, idxa0]
        distK <- sum((mu * dxi + dfi) / alpha) / k
        distFull <- mu * distX[i, ] + disty[i, ]
        d <- (distK - distFull)
        d[d <= 0] <- 0
        d[i] <- 0
        S[i, ] <- d
        S[i, ] <- d / sum(d)
    }

    SS <- (S + t(S)) / 2
    D <- diag(rowSums(SS))
    L <- D - SS

    return(list(L = L, SS = SS))
}