source("./src/core/gradient.r")
source("./src/core/calculate_graph_weights.r")
source("./src/core/l1_projection.r")

library(fastmatrix)

optimize <- function(data, lambda, eta, tau, w_init) {
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

    gradient_history <- as.numeric(c())

    while (iterations < max_iterations && (!((iterations >= (max_iterations / 2) &&
        abs(gradient_history[length(gradient_history)] - gradient_history[length(gradient_history) - 1]) <= 1e-6)))) {

        cat(sprintf("\nIteration: %d", iterations))

        alpha <- (t_old - 1) / t
        ws <- (1 + alpha) * wz - alpha * wz_old

        list[gws, Fs] <- evaluate_gradient(data, ws, eta, features_count)
                
        k <- 10

        all_features <- all_features[data[[1]]$patients, ]

        for (i in 2:length(data)) {
            new_data <- data[[i]]$features
            new_data <- new_data[data[[i]]$patients, ]
            all_features <- cbind(all_features, new_data)
        }

        #TODO Naming
        graph <- CalGraphWeight(all_features, non_repeated_features, ws, k)
        
        tmp <- tau * crossprod(all_features, graph$L)
        tmp2 <- all_features %*% ws
        res <- tmp %*% tmp2

        gws <- gws + res
        
        inner_iter <- 0
        max_inner_iter <- 1000

        r_sum <- 1
        Fzp_gamma <- 0
        Fzp <- 1

        wzp <- NULL

        while (inner_iter < max_inner_iter && (r_sum > 1e-20 && (is.nan(Fzp) || Fzp > Fzp_gamma))) {
            wzp <- l1_projection(ws - gws / gamma, lambda / gamma)
            Fzp <- neglogparlike(wzp$z, data)
            delta_wzp <- wzp$z - ws
            r_sum <- matrix.norm(delta_wzp, type = "Frobenius") ^ 2
            Fzp_gamma <- Fs + sum(delta_wzp * gws) + gamma / 2 * matrix.norm(delta_wzp, type = "Frobenius") ^ 2
            gamma <- gamma * gamma_inc
            inner_iter <- inner_iter + 1
        }

        wz_old <- wz
        wz <- wzp$z

        gradient_history <- c(gradient_history, Fzp)

        iterations <- iterations + 1
        t_old <- t
        t <- 0.5 * (1 + sqrt(1 + 4 * t^2))
    }

    return(list(wzp = wzp, gradient_history = gradient_history, s = graph$SS, l = graph$L))
}