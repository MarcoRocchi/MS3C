#TODO nomenclatura omogenea
#TODO output anche codice paziente (sort)

Sys.setenv(LANG = "en")

source("./src/load_data.r")
source("./src/preprocessing/preprocessing_pipeline.r")
source("./src/core/build_graph.r")
source("./src/core/build_graph_cox.r")
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

dataset$radiomics_pre <- subset(dataset$radiomics_pre, select = -c(AH, AI))
dataset$radiomics_post <- subset(dataset$radiomics_post, select = -c(AH, AI))
coefs <- msm$coefficients[1:25]

list[features, times, responses, frequencies, atrisk, tied] <- 
    do_preprocessing_s2gc(dataset$radiomics_pre, dataset$relapse_times, dataset$relapse_status)

#TODO grid search
lambda <- 1
eta <- 0.1
tau <- 10

message("Building similarity graph")
result <- build_graph(features, frequencies, responses, atrisk, lambda, eta, tau)

#TODO grid search
optimal_clusters_number <- 2
message("Performing spectral clustering")
clusters <- spectral_clustering(result$s, optimal_clusters_number)

plot_clusters(features, times, responses, clusters$group)

message("End")