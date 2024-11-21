#TODO nomenclatura omogenea
#TODO output anche codice paziente (sort)
#TODO verificare correttezza con Matlab S2GC
#TODO frequencies non serve

Sys.setenv(LANG = "en")

source("./src/load_data.r")
source("./src/preprocessing.r")
source("./src/prepare_data.r")
source("./src/models/model6.r")
source("./src/core/build_graph.r")
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

processed_dataset <- list()

for (d in dataset$data) {
    processed_dataset[[length(processed_dataset) + 1]] <- prepare_data(d$features, d$times, d$status)
}

cat("\nBuilding similarity graph")
lambda <- 1
eta <- 0.1
tau <- 10
result <- build_graph(processed_dataset, dataset$non_repeated_features, lambda, eta, tau)

cat("\nPerforming spectral clustering")
optimal_clusters_number <- sum(eigen(result$l, only.values = TRUE)$values == 0)
#TODO
optimal_clusters_number <- 2
cat(sprintf("\nOptimal clusters number: %d", optimal_clusters_number))
clusters <- spectral_clustering(result$s, optimal_clusters_number)

#TODO
#plot_clusters(dataset$non_repeated_features, times, responses, clusters$group)

cat("\nEnd")