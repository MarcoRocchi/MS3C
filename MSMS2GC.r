#TODO nomenclatura omogenea
#TODO verificare correttezza con Matlab S2GC
#TODO frequencies, tied non serve

Sys.setenv(LANG = "en")

source("./src/data/load_data.r")
source("./src/data/preprocessing.r")
source("./src/data/prepare_data.r")
source("./src/models/model0.r")
source("./src/core/optimize.r")
source("./src/clustering/spectral_clustering.r")
source("./src/plot/km_plot.r")

library(mstate)
library(dplyr)
library(gsubfn)

dataset <- load_data()
patients_count <- dataset$patients_count
dataset <- preprocess(dataset)

dataset <- expand_dataset(dataset)
dataset <- split_by_transition(dataset, patients_count)

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

candidate_lambda <- list(0.01, 0.1, 0.5, 1, 1.5, 5, 10)
candidate_eta <- list(0.01, 0.1, 0.5, 1, 1.5, 5, 10)
candidate_tau <- list(0.01, 0.1, 0.5, 1, 1.5, 5, 10)

for (a in candidate_lambda){
    for (b in candidate_eta) {
        for (c in candidate_tau) {
            result <- build_graph(processed_dataset, dataset$non_repeated_features, a, b, c)
            if (sum(eigen(result$l, only.values = TRUE)$values < 1e-10) > maximum_value) {
                maximum_value <- sum(eigen(result$l, only.values = TRUE)$values < 1e-10)
                l_opt <- a
                e_opt <- b
                t_opt <- c
            }
        }
    }
}

lambda <- 1
eta <- 0.1
tau <- 10

result <- optimize(dataset, patients_count, lambda, eta, tau)

cat("\nPerforming spectral clustering")
optimal_clusters_number <- sum(eigen(result$L, only.values = TRUE)$values < 1e-10)
cat(sprintf("\nOptimal clusters number: %d", optimal_clusters_number))
clusters <- spectral_clustering(result$S, 2)

#TODO
#plot_clusters(dataset$non_repeated_features, times, responses, clusters$group)
cat("\nEnd")