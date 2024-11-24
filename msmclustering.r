#TODO nomenclatura omogenea
#TODO verificare correttezza con Matlab S2GC
#TODO frequencies non serve

Sys.setenv(LANG = "en")

source("./src/data/load_data.r")
source("./src/data/preprocessing.r")
source("./src/data/prepare_data.r")
source("./src/models/model0.r")
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

l_opt <- 0
e_opt <- 0
t_opt <- 0
eig <- 0

lambda <- 1
eta <- 0.1
tau <- 10

#for (a in list(0.01, 0.1, 0.5, 1, 1.5, 5, 10)){
#    for (b in list(0.01, 0.1, 0.5, 1, 1.5, 5, 10)) {
#        for (c in list(0.01, 0.1, 0.5, 1, 1.5, 5, 10)) {
#            result <- build_graph(processed_dataset, dataset$non_repeated_features, a, b, c)
#            if (sum(eigen(result$l, only.values = TRUE)$values < 1e-10) > eig) {
#                eig <- sum(eigen(result$l, only.values = TRUE)$values < 1e-10)
#                l_opt <- a
#                e_opt <- b
#                t_opt <- c
#            }
#        }
#    }
#}

result <- build_graph(processed_dataset, dataset$non_repeated_features, lambda, eta, tau)

cat("\nPerforming spectral clustering")
optimal_clusters_number <- sum(eigen(result$l, only.values = TRUE)$values < 1e-10)
cat(sprintf("\nOptimal clusters number: %d", optimal_clusters_number))
clusters <- spectral_clustering(result$s, 4)

#TODO
#plot_clusters(dataset$non_repeated_features, times, responses, clusters$group)
cat("\nEnd")