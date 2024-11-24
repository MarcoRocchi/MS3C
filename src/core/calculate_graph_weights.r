#Estimate S by fixing w
CalGraphWeight <- function(features, non_repeated_features, w, k) {
    #TODO ogni transizone ha il suo w
    #TODO rendere mu parametrizzabile
    mu <- 1e-3
    alpha <- 1
    n <- nrow(features)

    #TODO restore non_repeated features but preprocessed
    distX <- as.matrix(dist(features, diag = TRUE, upper = TRUE))
    distX <- distX * distX

    y <- features %*% w
    disty <- as.matrix(dist(y, diag = TRUE, upper = TRUE))
    disty <- disty * disty

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
    D <- diag(colSums(SS))
    L <- D - SS

    return(list(L = L, SS = SS))
}