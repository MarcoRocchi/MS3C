#TODO nomenclatura omogenea
#TODO output anche codice paziente (sort)

Sys.setenv(LANG = "en")

source("./src/load_data.r")
source("./src/preprocessing/preprocessing_pipeline.r")
source("./src/core/build_graph.r")
source("./src/clustering/spectral_clustering.r")
source("./src/plot/km_plot.r")
source("./src/models/model6.r")

library(gsubfn)

message("Loading data")
dataset <- load_data()

message("Preprocessing data")
preprocessed_pre <- preprocess(dataset$radiomics_pre)
preprocessed_post <- preprocess(dataset$radiomics_post)
preprocessed_pre_operative <- preprocess_pre_operative(dataset$pre_operative)

msm <- build_m6(
    preprocessed_pre,
    preprocessed_post,
    preprocessed_pre_operative,
    dataset$relapse_status,
    dataset$dead_status,
    dataset$post_times,
    dataset$relapse_times,
    dataset$dead_times
)

dataset$radiomics_pre <- subset(dataset$radiomics_pre, select = -c(AAM))
coefs <- msm$coefficients[c(1:12, 14:27)]
#coefs <- msm$coefficients[1:26]

list[features, times, responses, frequencies, atrisk, tied] <- 
    do_preprocessing(dataset$radiomics_pre, dataset$relapse_times, dataset$relapse_status)

#TODO grid search
lambda <- 1
eta <- 0.1
tau <- 10

#GRID
if(FALSE) {
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
#------------

message("Building similarity graph")
result <- build_graph(features, frequencies, responses, atrisk, lambda, eta, tau, coefs)

message("Performing spectral clustering")
#optimal_clusters_number <- sum(eigen(result$l, only.values = TRUE)$values == 0)
optimal_clusters_number <- 2
clusters <- spectral_clustering(result$s, optimal_clusters_number)

plot_clusters(features, times, responses, clusters$group)

message("End")