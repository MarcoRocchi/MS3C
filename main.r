#TODO nomenclatura omogenea
#TODO output anche codice paziente (sort)

Sys.setenv(LANG = "en")

source("./src/load_data.r")
source("./src/preprocessing/preprocessing_pipeline.r")
source("./src/core/build_graph.r")
source("./src/clustering/spectral_clustering.r")
source("./src/plot/km_plot.r")

library(gsubfn)

message("Loading data")
list[features, times, responses] <- load_data_post()

message("Preprocessing data")
list[features, times, responses, frequencies, atrisk, tied] <- do_preprocessing(features, times, responses)

#TODO grid search
lambda <- 1
eta <- 0.1
tau <- 10

message("Building similarity graph")
result <- build_graph(features, frequencies, responses, atrisk, lambda, eta, tau)

#TODO grid search
optimal_clusters_number <- 4
message("Performing spectral clustering")
clusters <- spectral_clustering(result$s, optimal_clusters_number)

plot_clusters(features, times, responses, clusters$group)

message("End")