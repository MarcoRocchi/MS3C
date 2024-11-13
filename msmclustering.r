#TODO nomenclatura omogenea
#TODO output anche codice paziente (sort)

Sys.setenv(LANG = "en")

source("./src/load_data.r")
source("./src/preprocessing/preprocessing_pipeline.r")
source("./src/preprocessing/build_long_dataset.r")
source("./src/core/build_graph.r")
source("./src/clustering/spectral_clustering.r")
source("./src/plot/km_plot.r")

library(mstate)
library(dplyr)
library(gsubfn)

dataset <- load_data()
patients_count <- dataset$patients_count
dataset <- preprocess(dataset)
dataset <- expand_transitions(dataset)
#TODO add one level of abstraction?
dataset <- group_split(dataset, dataset$trans)

#TODO
for (d in dataset) {
    for (i in 1:patients_count) {
        patients <- as.list(d["id"])$id
        if (!i %in% patients) {
            print("Missing patient")
        }  
    }
}

processed_dataset <- list()

for (d in dataset) {
    #TODO transition specific
    #TODO portare dentro al codice dei modelli
    features <- as.matrix(d[c(names.radiomics_pre, names.pre_operative)])
    times <- as.matrix(d["time"])
    status <- as.matrix(d["status"])

    processed_dataset[[length(processed_dataset) + 1]] <- do_preprocessing(features, times, status)
}

#TODO gestire ordinamento diverso dei pazienti nei vari stati

#------------------------------------------------------------------------------------------------------------

#TODO grid search
lambda <- 1
eta <- 0.1
tau <- 10

#GRID
if (FALSE) {
    for (l in 1:20) {
        for (e in 1:20) {
            for (t in 1:100) {
                result <- build_graph(features, frequencies, responses, atrisk, l * 0.2, e * 0.1, t * 0.5)
                optimal_clusters_number <- sum(eigen(result$l, only.values = TRUE)$values == 0)

                if (optimal_clusters_number != 0) {
                    print("AAAA")
                }

            }
        }
    }
}
#------------------------------------------------------------------------------------------------------------

message("Building similarity graph")
result <- build_graph(processed_dataset, lambda, eta, tau)

message("Performing spectral clustering")
optimal_clusters_number <- sum(eigen(result$l, only.values = TRUE)$values == 0)
#TODOs
optimal_clusters_number <- 2
clusters <- spectral_clustering(result$s, optimal_clusters_number)

plot_clusters(features, times, responses, clusters$group)

message("End")