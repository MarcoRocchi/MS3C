discretisation <- function(EigenVectors) {
    #TODO
    EigenVectors[is.na(EigenVectors)] <- 0
    n <- nrow(EigenVectors)
    k <- ncol(EigenVectors)
    
    vm <- sqrt(rowSums(EigenVectors^2))
    EigenVectors <- EigenVectors / (vm + .Machine$double.eps)
    
    R <- matrix(0, nrow = k, ncol = k)
    R[, 1] <- EigenVectors[round(n / 2), ]
    c <- numeric(n)
    
    for (j in 2:k) {
        c <- c + abs(EigenVectors %*% R[, j - 1])
        minimum <- min(c)
        i <- which.min(c)
        R[, j] <- EigenVectors[i, ]
    }
    
    lastObjectiveValue <- 0
    exitLoop <- 0
    nbIterationsDiscretisation <- 0
    nbIterationsDiscretisationMax <- 20
    
    while (exitLoop == 0) {
        nbIterationsDiscretisation <- nbIterationsDiscretisation + 1
        EigenvectorsDiscrete <- discretisationEigenVectorData(EigenVectors %*% R)
        svd_result <- svd(t(EigenvectorsDiscrete) %*% EigenVectors + .Machine$double.eps)
        NcutValue <- 2 * (n - sum(svd_result$d))
        
        if (abs(NcutValue - lastObjectiveValue) < .Machine$double.eps || nbIterationsDiscretisation > nbIterationsDiscretisationMax) {
            exitLoop <- 1
        } else {
            lastObjectiveValue <- NcutValue
            R <- svd_result$v %*% t(svd_result$u)
        }
    }

    return(list(eigenvectorsDiscrete = EigenvectorsDiscrete, eigenVectors = EigenVectors))
}

discretisationEigenVectorData <- function(EigenVector) {
    n <- nrow(EigenVector)
    k <- ncol(EigenVector)
    
    Maximum <- apply(EigenVector, 1, max)
    J <- apply(EigenVector, 1, which.max)
    
    Y <- Matrix::sparseMatrix(i = 1:n, j = J, x = 1, dims = c(n, k))
    return(Y)
}