source("./src/core/cox.r")
source("./src/core/similarity.r")
source("./src/core/l1_projection.r")

library(fastmatrix)
library(ggplot2)

plot_training_process <- function(likelihood_history) {
    # Create a data frame for plotting
    training_data <- data.frame(
        Epoch = 1:length(likelihood_history),
        Loss = likelihood_history
    )
    
    # Plot the training process
    plot <- ggplot(training_data, aes(x = Epoch, y = Loss)) +
        geom_line(color = "blue", size = 1) +
        ggtitle("Loss Function Over Epochs") +
        xlab("Epoch") +
        ylab("Loss") +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12)
        )

        print(plot)
}

initialize_null <- function(data) {
    result <- list()    

    for (i in 1:length(data)) {
        d <- data[[i]]
        result[[i]] <- matrix(0, ncol(d$features), 1)
    }

    return(result)
}

inner_loop <- function(data, ws, gws, eta, tau, tau_inc, likelihood_cox) {
    inner_iter <- 0
    max_inner_iter <- 1000

    r_sum <- 1
    likelihood_gamma <- 0
    likelihood <- 1

    wzp <- initialize_null(data)

    while (inner_iter < max_inner_iter && r_sum > 1e-20 && likelihood > likelihood_gamma) {
        total_penalty <- 0

        for (i in 1:length(data)) {
            wzp[[i]] <- l1_projection(ws[[i]] - gws[[i]] / tau, eta / tau)
        }    

        likelihood <- neglogparlike(lapply(wzp, function(x) lapply(x, return(x$z))), data)

        for (i in 1:length(data))  {
            delta_wzp <- wzp[[i]]$z - ws[[i]]
            r_sum <- matrix.norm(delta_wzp, type = "Frobenius") ^ 2
            penalty <- sum(delta_wzp * gws[[i]]) + tau / 2 * matrix.norm(delta_wzp, type = "Frobenius") ^ 2
            total_penalty <- total_penalty + penalty
        }

        likelihood_gamma <- likelihood_cox + total_penalty
        tau <- tau * tau_inc
        inner_iter <- inner_iter + 1
    }

    return(list(wzp, tau, likelihood))
}

optimize <- function(data, n, eta, gamma, mu, k) {
    wz <- initialize_null(data)
    wz_old <- initialize_null(data)
    ws <- initialize_null(data)

    t <- 1
    t_old <- 0

    tau <- 1
    tau_inc <- 2

    epochs <- 0
    max_epochs <- 1000

    likelihood_history <- as.numeric(c())

    best_wzp <- initialize_null(data)
    best_likelihood <- NULL
    best_graph <- NULL
    not_improving <- 0

    while (epochs < max_epochs && not_improving < 10) {
        alpha <- (t_old - 1) / t
        if (alpha > 0) {
            alpha <- -log(alpha)
        }

        for (i in 1:length(data)) {
            ws[[i]] <- (1 + alpha) * wz[[i]] - alpha * wz_old[[i]]
        }

        list[gws, likelihood_cox] <- evaluate_gradient(data, ws)

        graph <- estimate_similarity(data, n, ws, mu, k)

        hy <- gamma * mu

        for (i in 1:length(data)) {
            features <- data[[i]]$features[data[[i]]$patients, ]
            variance <- compute_variance(wz[[i]], features, data[[i]]$time, 1 - data[[i]]$censoring)
            transition_weight <- 1 / variance
            data[[i]]$transition_weight <- transition_weight
            gws[[i]] <- gws[[i]] + hy * transition_weight * crossprod(features, graph$L) %*% (features %*% ws[[i]])
        }

        list[wzp, tau, likelihood] <- inner_loop(data, ws, gws, eta, tau, tau_inc, likelihood_cox)        
        
        for (i in 1:length(data)) {
            wz_old[[i]] <- wz[[i]]
            wz[[i]] <- wzp[[i]]$z
        }

        if (is.null(best_likelihood) || best_likelihood > likelihood) {
            best_likelihood <- likelihood
            best_wzp <- wzp
            best_graph <- graph
            not_improving <- 0

        } else {
            not_improving <- not_improving + 1
        }

        likelihood_history <- c(likelihood_history, likelihood)

        epochs <- epochs + 1
        
        t_old <- t
        t <- 1 + sqrt(1 + 4 * t^2)
    }

    #plot_training_process(likelihood_history)
    return(list(weights = best_wzp, likelihood_history = likelihood_history, S = best_graph$S, L = best_graph$L))
}