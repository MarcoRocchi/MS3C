source("./src/core/optimize.r")
source("./src/clustering/spectral_clustering.r")
source("./src/validation/logrank.r")
source("./src/validation/concordance.r")

grid_search <- function(dataset, patients_count, verbose = FALSE) {
    cat("\nSearching optimal parameters")

    e_opt <- 0
    g_opt <- 0
    k_opt <- 0
    m_opt <- 0
    maximum_value <- 0

    candidate_eta <- list(1e-4, 1e-3, 1e-2, 0.1, 0.5)
    candidate_mu <- list(1e-5, 1e-4, 1e-3, 1e-2)
    candidate_k <- list(3, 4, 5, 6)
    candidate_gamma <- list(5e2, 2e3, 3e3, 4e3)

    for (e in candidate_eta) {
        for (g in candidate_gamma) {
            for (k in candidate_k) {
                for (m in candidate_mu) {
                    result <- optimize(dataset, patients_count, e, g, m, k)
                    optimal_clusters_number <- sum(eigen(result$L, only.values = TRUE)$values < 1e-10)
                    if (optimal_clusters_number > 1 && optimal_clusters_number <= 5) {
                        clusters <- spectral_clustering(result$S, optimal_clusters_number)
                        current_lr <- logrank(mstate_dataset, clusters$group)[["logtest"]]["pvalue"]
                        if (current_lr <= 0.05) {
                            current_c <- concordance_index(dataset, result$weights, patients_count)
                            if (current_c > maximum_value) {
                                maximum_value <- current_c
                                e_opt <- e
                                g_opt <- g
                                m_opt <- m
                                k_opt <- k
                            }
                        }
                    }
                }
            }
        }
    }

    if (verbose) {
        cat(sprintf("\nMaximum concordance: %f", maximum_value))
        cat("\n eta ", e_opt)
        if (e_opt == candidate_eta[[1]] || e_opt == candidate_eta[[length(candidate_eta)]]) {
            cat(" !!!!!")
        }

        cat("\n gamma ", g_opt)
        if (g_opt == candidate_gamma[[1]] || g_opt == candidate_gamma[[length(candidate_gamma)]]) {
            cat(" !!!!!")
        }

        cat("\n mu ", m_opt)
        if (m_opt == candidate_mu[[1]] || m_opt == candidate_mu[[length(candidate_mu)]]) {
            cat(" !!!!!")
        }

        cat("\n k ", k_opt)
        if (k_opt == candidate_k[[1]] || k_opt == candidate_k[[length(candidate_k)]]) {
            cat(" !!!!!")
        }
    }

    return(list(
        eta = e_opt,
        gamma = g_opt,
        mu = m_opt,
        k = k_opt,
        concordance = maximum_value
    ))
}