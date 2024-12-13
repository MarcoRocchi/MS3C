source("./src/core/gradient.r")
source("./src/core/similarity.r")
source("./src/core/l1_projection.r")

library(fastmatrix)

initialize_null <- function(data) {
    result <- list()    

    for (i in 1:length(data)) {
        result[[i]] <- matrix(0, ncol(d$features), 1)
    }

    return(result)
}

inner_loop <- function(data, ws, gws, lambda, gamma, gamma_inc, likelihood_cox) {
    inner_iter <- 0
    max_inner_iter <- 1000

    r_sum <- 1
    likelihood_gamma <- 0
    likelihood <- 1

    wzp <- initialize_null(data)

    while (inner_iter < max_inner_iter && r_sum > 1e-20 && likelihood > likelihood_gamma) {
        total_penalty <- 0

        for (i in 1:length(data)) {
            #TODO verificare ok l1_projection parziale
            wzp[[i]] <- l1_projection(ws[[i]] - gws[[i]] / gamma, lambda / gamma)
        }    

        likelihood <- neglogparlike(lapply(wzp, function(x) lapply(x, return(x$z))), data)

        for (i in 1:length(data))  {
            delta_wzp <- wzp[[i]]$z - ws[[i]]
            r_sum <- matrix.norm(delta_wzp, type = "Frobenius") ^ 2
            penalty <- sum(delta_wzp * gws[[i]]) + gamma / 2 * matrix.norm(delta_wzp, type = "Frobenius") ^ 2
            total_penalty <- total_penalty + penalty
        }

        likelihood_gamma <- likelihood_cox + total_penalty
        gamma <- gamma * gamma_inc
        inner_iter <- inner_iter + 1
    }

    return(list(wzp, gamma, likelihood))
}

optimize <- function(data, n, lambda, eta, tau, w_init) {
    wz <- initialize_null(data)
    wz_old <- initialize_null(data)
    ws <- initialize_null(data)

    t <- 1
    t_old <- 0

    gamma <- 1
    gamma_inc <- 2

    iterations <- 0
    max_iterations <- 100

    likelihood_difference <- 0
    likelihood_history <- as.numeric(c())

    while (iterations < max_iterations && (!((iterations >= (max_iterations / 2) && likelihood_difference <= 1e-6)))) {
        cat(sprintf("\nIteration: %d", iterations))

        alpha <- (t_old - 1) / t

        for (i in 1:length(data)) {
            ws[[i]] <- (1 + alpha) * wz[[i]] - alpha * wz_old[[i]]
        }

        list[gws, likelihood_cox] <- evaluate_gradient(data, ws, eta)
                
        k <- 10

        graph <- estimate_similarity(data, n, ws, k)

        for (i in 1:length(data)) {
            features <- data[[i]]$features[data[[i]]$patients, ]
            transition_weight <- data[[i]]$transition_weight
            gws[[i]] <- gws[[i]] + tau * transition_weight * crossprod(features, graph$L) %*% (features %*% ws[[i]])
        }

        list[wzp, gamma, likelihood] <- inner_loop(data, ws, gws, lambda, gamma, gamma_inc, likelihood_cox)        

        for (i in 1:length(data)) {
            wz_old[[i]] <- wz[[i]]
            wz[[i]] <- wzp[[i]]$z
        }

        likelihood_history <- c(likelihood_history, likelihood)
        likelihood_difference <- abs(
            likelihood_history[length(likelihood_history)] - 
            likelihood_history[length(likelihood_history) - 1])
        iterations <- iterations + 1
        t_old <- t
        t <- 0.5 * (1 + sqrt(1 + 4 * t^2))
    }

    return(list(weights = wzp, likelihood_history = likelihood_history, S = graph$S, L = graph$L))
}