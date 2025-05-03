Sys.setenv(LANG = "en")

source("./src/simulation/progressive.r")
source("./src/data/preprocessing.r")
source("./src/data/prepare_data.r")
source("./src/core/optimize.r")
source("./src/clustering/spectral_clustering.r")
source("./src/validation/concordance.r")
source("./src/validation/logrank.r")
source("./src/validation/classificator.r")
source("./src/plot/plot_km.r")

library(dplyr)
library(gsubfn)
library(clue)

set.seed(3)

n_g1 <- 50
n_g2 <- 30

g1 <- generate_sample_data(
    n_g1, 
    time_params = list(mean1 = 10, sd1 = 2, mean2 = 40, sd2 = 10), 
    haz_params = list(mean1 = 0.1, sd1 = 0.04, mean2 = 0.4, sd2 = 0.05), 
    censor_params = list(p1 = 0.2, p2 = 0.4),
    covariates_params = list(mean1 = 7, sd1 = 2, mean2 = 50, sd2 = 8),
    categorical_params = list(p1 = list(0.3, 0.5, 0.1, 0.1, 0.2), p2 = list(0.1, 0.1, 0.1, 0.2, 0.5))
)
g2 <- generate_sample_data(
    n_g2, 
    time_params = list(mean1 = 45, sd1 = 5, mean2 = 120, sd2 = 10), 
    haz_params = list(mean1 = 0.4, sd1 = 0.04, mean2 = 0.4, sd2 = 0.04), 
    censor_params = list(p1 = 0, p2 = 0.5),
    covariates_params = list(mean1 = 45, sd1 = 4, mean2 = 150, sd2 = 10),
    categorical_params = list(p1 = list(0.05, 0.35, 0.15, 0.4, 0.15), p2 = list(0.4, 0.6, 0.1, 0.005, 0.05))
)

g1$id <- g1$id
g2$id <- g2$id + max(g1$id)

dataset <- rbind(g1, g2)

classification_dataset <- get_classification_features(dataset)

true_label <- rep(1, n_g1)
true_label <- c(true_label, rep(2, n_g2))

patients_count <- length(unique(dataset$id))

build_model(dataset)

dataset <- split_by_transition(dataset, patients_count)

i <- 1

for (d in dataset) {
    dataset[[i]] <- prepare_data(d)
    i <- i + 1
}

cat("\nSearching optimal parameters")
e_opt <- 1e-04
g_opt <- 500
k_opt <- 6
m_opt <- 0.001

if (FALSE) {
    maximum_value <- 0

    candidate_eta <- list(1e-4, 1e-3, 1e-2, 0.1, 0.5)
    candidate_mu <- list(1e-5, 1e-4, 1e-3, 1e-2)
    candidate_k <- list(4, 5, 6, 7, 8)
    candidate_gamma <- list(5e2, 2e3, 3e3, 4e3)

    for (e in candidate_eta) {
        for (g in candidate_gamma) {
            for (k in candidate_k) {
                for (m in candidate_mu) {
                    result <- optimize(dataset, patients_count, e, g, m, k)
                    optimal_clusters_number <- sum(eigen(result$L, only.values = TRUE)$values < 1e-10)
                    if (optimal_clusters_number > 1) {
                        clusters <- spectral_clustering(result$S, optimal_clusters_number)
                        current_concordance <- concordance_index(dataset, result$weights, patients_count)
                        if (current_concordance > maximum_value) {
                            maximum_value <- current_concordance
                            e_opt <- e
                            g_opt <- g
                            m_opt <- m
                            k_opt <- k

                            cat("\nConcordance: ", current_concordance)
                        }
                    }
                }
            }
        }
    }

    cat("\n ", maximum_value)
}

cat("\nBuilding similarity graph")
result <- optimize(dataset, patients_count, e_opt, g_opt, m_opt, k_opt)
cat("\nConcordance:", concordance_index(dataset, result$weights, patients_count))

cat("\nPerforming spectral clustering")
optimal_clusters_number <- sum(eigen(result$L, only.values = TRUE)$values < 1e-20)
cat(sprintf("\nOptimal clusters number: %d", optimal_clusters_number))
clusters <- spectral_clustering(result$S, optimal_clusters_number)

true_label <- as.numeric(true_label)
clusters$group <- as.numeric(clusters$group)
assignment <- solve_LSAP(table(true_label, clusters$group), maximum = TRUE)
optimal_clusters <- clusters$group

plot_km_curve(dataset, clusters$group)

for (i in 1:2) {
    correct <- sum(true_label == i & clusters$group == assignment[i])
    total <- sum(true_label == i)
    percentage <- correct / total * 100
    cat(sprintf("\nLabel %d: %.2f%% correctly classified", i, percentage))
}

cat("End of simulation study")