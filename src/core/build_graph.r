source("./src/core/gradient.r")
source("./src/core/calculate_graph_weights.r")
source("./src/core/l1_projection.r")
source("./src/core/l2_distance.r")

library(fastmatrix)

build_graph <- function(data, non_repeated_features, lambda, eta, tau, w_init) {
    features_count <- 0
    for (d in data) {
        features_count <- features_count + ncol(d$features)
    }

    if (missing(w_init)) {
        wz <- matrix(0, features_count, 1)
    } else {
        wz <- unname(w_init)
    }
    
    wz_old <- wz

    t <- 1
    t_old <- 0

    iterations <- 0

    gamma <- 1
    gamma_inc <- 2

    max_iterations <- 100

    #TODO
    funcVal <- as.numeric(list(0, 0))

    while (!((iterations > max_iterations) || (iterations > (max_iterations / 2) &&
            abs(funcVal[length(funcVal)] - funcVal[length(funcVal) - 1]) <= 1e-6))) {

        cat(sprintf("\nIteration: %d", iterations))

        alpha <- (t_old - 1) / t
        ws <- (1 + alpha) * wz - alpha * wz_old

        list[gws, Fs] <- evaluate_gradient(data, ws, eta, features_count)
                
        k <- 10

        all_features <- data[[1]]$features

        for (i in 2:length(data)) {
            all_features <- cbind(all_features, data[[i]]$features)
        }

        graph <- CalGraphWeight(all_features, non_repeated_features, ws, k)
        
        #TODO non ha senso. Devo usare i pazienti sempre nello stesso ordine
        tmp <- tau * crossprod(all_features, graph$L)
        tmp2 <- all_features %*% ws
        res <- tmp %*% tmp2

        gws <- gws + res
        
        innerIter <- 0
        maxInnerIter <- 1000

        r_sum <- 1
        Fzp_gamma <- 0
        Fzp <- 1

        wzp <- NULL

        while (innerIter < maxInnerIter && r_sum > 1e-20 && (is.nan(Fzp) || Fzp > Fzp_gamma)) {
            wzp <- l1_projection(ws - gws / gamma, lambda / gamma)
            Fzp <- neglogparlike(wzp$z, data)
            delta_wzp <- wzp$z - ws
            r_sum <- matrix.norm(delta_wzp, type = "Frobenius") ^ 2
            Fzp_gamma <- Fs + sum(delta_wzp * gws) + gamma / 2 * matrix.norm(delta_wzp, type = "Frobenius") ^ 2
            gamma <- gamma * gamma_inc
            innerIter <- innerIter + 1
        }

        wz_old <- wz
        wz <- wzp$z

        funcVal <- c(funcVal, Fzp)

        iterations <- iterations + 1
        t_old <- t
        t <- 0.5 * (1 + sqrt(1 + 4 * t^2))
    }

    return(list(wzp = wzp, funcval = funcVal, s = graph$SS, l = graph$L))
}