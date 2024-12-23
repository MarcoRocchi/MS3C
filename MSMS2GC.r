#TODO nomenclatura omogenea
#TODO frequencies, tied non serve
#TODO verificare meccanismo ordinamento individui
#TODO Test logrank (anche tempo all'evento)
Sys.setenv(LANG = "en")

source("./src/data/load_data.r")
source("./src/data/preprocessing.r")
source("./src/data/prepare_data.r")
source("./src/models/model0.r")
source("./src/core/optimize.r")
source("./src/clustering/spectral_clustering.r")
source("./src/plot/km_plot.r")
source("./src/validation/concordance.r")

library(mstate)
library(dplyr)
library(gsubfn)

dataset <- load_data()
patients_count <- dataset$patients_count
dataset <- preprocess(dataset)

mstate_dataset <- expand_dataset(dataset)
dataset <- split_by_transition(mstate_dataset, patients_count)

i <- 1

for (d in dataset) {
    dataset[[i]] <- prepare_data(d)
    i <- i + 1
}

cat("\nBuilding similarity graph")

l_opt <- 0
e_opt <- 0
t_opt <- 0
maximum_value <- 0

candidate_lambda <- list(0.05, 0.5, 1, 1.5, 2.5, 3.5)
candidate_eta <- list(0.001, 0.01, 0.1, 0.5)
candidate_tau <- list(0.5, 1, 3, 4, 5, 6, 7, 10)

search_optimal_parameters <- FALSE

if (search_optimal_parameters) {
    for (a in candidate_lambda) {
        for (b in candidate_eta) {
            for (c in candidate_tau) {
                result <- optimize(dataset, patients_count, a, b, c)
                current_concordance <- concordance_index(dataset, result$weights, patients_count)
                if (current_concordance > maximum_value) {
                    maximum_value <- current_concordance
                    l_opt <- a
                    e_opt <- b
                    t_opt <- c
                }
            }
        }
    }
}

lambda <- 0.05
eta <- 0.001
tau <- 3

result <- optimize(dataset, patients_count, lambda, eta, tau)

cat("\nPerforming spectral clustering")
optimal_clusters_number <- sum(eigen(result$L, only.values = TRUE)$values < 1e-10)
cat(sprintf("\nOptimal clusters number: %d", optimal_clusters_number))
clusters <- spectral_clustering(result$S, 2)

#TODO
#plot_clusters(dataset$non_repeated_features, times, responses, clusters$group)
cat("\nEnd")