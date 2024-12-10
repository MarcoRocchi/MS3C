compute_distance <- function(data) {
    features <- list()

    for (d in data) {
        features <- cbind(features, d$features[d$patients, d$new_features])
    }

    return(as.matrix(dist(features, diag = TRUE, upper = TRUE)))
}

compute_survival_distance <- function(data, w, n) {
    y <- matrix(0, n, 1)

    for (i in length(data)) {
        #TODO sort by same order
        #TODO multiply by weight
        y <- y + (data[[i]]$features %*% w[[i]])[data[[i]]$patients]
    }

    return(as.matrix(dist(y, diag = TRUE, upper = TRUE)))
}

#Estimate S by fixing w
estimate_similarity <- function(data, n, w, k) {
    #TODO rendere mu parametrizzabile
    mu <- 1e-3
    alpha <- 1

    patient_distances <- compute_distance(data)
    patient_distances <- patient_distances * patient_distances

    survival_distances <- compute_survival_distance(data, w, n)
    survival_distances <- survival_distances * survival_distances

    S <- matrix(0, n, n)
    
    #TODO verificare
    idx <- t(apply(mu * patient_distances + survival_distances, 2, order))
        
    for (i in 1:n) {
        idxa0 <- idx[i, 2:(k + 1)]
        dfi <- survival_distances[i, idxa0]
        dxi <- patient_distances[i, idxa0]
        distk <- sum((mu * dxi + dfi) / alpha) / k
        distance <- mu * patient_distances[i, ] + survival_distances[i, ]
        d <- (distk - distance)
        d[d <= 0] <- 0
        d[i] <- 0
        S[i, ] <- t(d)
        S[i, ] <- t(d) / sum(d)
    }

    SS <- (S + t(S)) / 2
    D <- diag(colSums(SS))
    L <- D - SS

    return(list(L = L, S = SS))
}