Sys.setenv(LANG = "en")

source("./src/data/load_data.r")
source("./src/data/preprocessing.r")
source("./src/data/prepare_data.r")
source("./src/models/model3.r")
source("./src/core/optimize.r")
source("./src/clustering/spectral_clustering.r")
source("./src/validation/concordance.r")
source("./src/validation/agreement.r")
source("./src/validation/logrank.r")
source("./src/plot/plot_km.r")
source("./src/validation/classifier.r")

library(mstate)
library(dplyr)
library(gsubfn)

dataset <- load_data()
patients_count <- dataset$patients_count
dataset <- preprocess(dataset)
classification_dataset <- get_classification_features(dataset)

mstate_dataset <- expand_dataset(dataset)
dataset <- split_by_transition(mstate_dataset, patients_count)

i <- 1

for (d in dataset) {
    dataset[[i]] <- prepare_data(d)
    i <- i + 1
}

cat("\nBuilding similarity graph")
list[eta, gamma, mu, k] <- get_optimal_parameters()
result <- optimize(dataset, patients_count, eta, gamma, mu, k)
cat("\nConcordance:", concordance_index(dataset, result$weights, patients_count))

cat("\nPerforming spectral clustering")
optimal_clusters_number <- sum(eigen(result$L, only.values = TRUE)$values < 1e-10)
cat(sprintf("\nOptimal clusters number: %d", optimal_clusters_number))
clusters <- spectral_clustering(result$S, optimal_clusters_number)

cat("\nComputing clustering logrank\n")
print(logrank(mstate_dataset, clusters$group))

plot_km_curve(dataset, clusters$group)

compute_classifier(classification_dataset, clusters$group)

#TODO
#cat("\nClustering agreement:\n")
#print(compute_survival_cindex(dataset[[4]]$times, dataset[[4]]$censoring, clusters$group))

cat("\nEnd")