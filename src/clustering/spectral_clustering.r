source("./src/clustering/discretisation.r")

library(Matrix)
library(Rcpp)

spectral_clustering <- function(w, num_clusters) {
    degs <- rowSums(w)
    D <- as.matrix(Diagonal(x = degs))

    L <- D - w

    degs[degs == 0] <- .Machine$double.eps
    D <- as.matrix(Diagonal(x = 1 / sqrt(degs)))
    L <- D %*% L %*% D

    eig_result <- eigen(L)
    eigenvalue <- diag(x = rev(tail(eig_result$values, num_clusters)))
    U <- t(as.matrix(tail(eig_result$vectors, num_clusters)))
    a <- order(diag(eigenvalue))
    eigenvalue <- eigenvalue[, a]
    U <- U[, a]
    eigengap <- abs(diff(eigenvalue))
    U <- U[, 1:num_clusters]

    flag <- 0
    Cluster <- vector("list", length(num_clusters))

    for (ck in num_clusters) {
        Cindex <- which(num_clusters == ck)
        UU <- U[, 1:ck]
        UU <- UU / sqrt(rowSums(UU^2))
        EigenvectorsDiscrete <- discretisation(UU)$eigenvectorsDiscrete
        temp <- max.col(EigenvectorsDiscrete)

        Cluster[[Cindex]] <- temp
    }

    if (length(num_clusters) == 1) {
        group <- Cluster[[1]]
    } else {
        group <- Cluster
    }

    return(list(group = group, eigengap = eigengap))
}