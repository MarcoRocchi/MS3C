compute_distance <- function(data) {
    features <- list()

    for (i in 1:length(data)) {
        d <- data[[i]]
        order <- data[[i]]$patients
        new_features <- data[[i]]$new_features
        features <- cbind(features, d$features[order, new_features])
    }

    return(as.matrix(dist(features, diag = TRUE, upper = TRUE)))
}

compute_risk_distance <- function(data, w, n) {
    y <- matrix(0, n, 1)

    for (i in 1:length(data)) {
        order <- data[[i]]$patients
        y <- y + (data[[i]]$transition_weight * data[[i]]$features[order, ] %*% w[[i]])
    }

    return(as.matrix(dist(y, diag = TRUE, upper = TRUE)))
}

#Estimate S by fixing w
estimate_similarity <- function(data, n, w, mu, k) {
    patient_distances <- compute_distance(data)
    patient_distances <- patient_distances * patient_distances

    risk_distances <- compute_risk_distance(data, w, n)
    risk_distances <- risk_distances * risk_distances

    S <- matrix(0, n, n)
    
    idx <- t(apply(mu * patient_distances + risk_distances, 2, order))
        
    for (i in 1:n) {
        idxa0 <- t(idx[i, 2:(k + 1)])
        dfi <- t(risk_distances[i, idxa0])
        dxi <- t(patient_distances[i, idxa0])
        distk <- sum(mu * dxi + dfi) / k
        distance <- t(mu * patient_distances[i, ] + risk_distances[i, ])
        d <- (distk - distance)
        d[d <= 0] <- 0
        d[i] <- 0
        S[i, ] <- d / sum(d)
    }

    SS <- (S + t(S)) / 2
    D <- diag(colSums(SS))
    L <- D - SS

    return(list(L = L, S = SS))
}