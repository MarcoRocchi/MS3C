Sys.setenv(LANG = "en")

source("./src/data/load_data.r")
source("./src/data/preprocessing.r")
source("./src/data/prepare_data.r")
source("./src/models/model0.r")
source("./src/core/optimize.r")
source("./src/validation/concordance.r")
source("./src/validation/logrank.r")

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

cat("\nSearching optimal parameters")

e_opt <- 0
t_opt <- 0
k_opt <- 0
m_opt <- 0
maximum_value <- 0

candidate_eta <- list(0.001, 0.01, 0.1, 0.5)
candidate_mu <- list(1e-6, 1e-4, 1e-3, 1e-2)
candidate_k <- list(3, 4, 5, 6)
candidate_tau <- list(1, 2, 3, 4)

for (e in candidate_eta) {
    for (t in candidate_tau) {
        for (k in candidate_k) {
            for (m in candidate_mu) {
                result <- optimize(dataset, patients_count, e, t, m, k)
                current_concordance <- concordance_index(dataset, result$weights, patients_count)
                if (current_concordance > maximum_value) {
                    maximum_value <- current_concordance
                    e_opt <- e
                    t_opt <- t
                    m_opt <- m
                    k_opt <- k

                    data <- c("Temporary", e_opt, t_opt, m_opt, k_opt, current_concordance)
                    write(data, file = "optimal_parameters.txt", append = FALSE, sep = " // ")
                }
            }
        }
    }
}

result <- optimize(dataset, patients_count, e_opt, t_opt, m_opt, k_opt)
optimal_clusters_number <- sum(eigen(result$L, only.values = TRUE)$values < 1e-10)
cat("\nConcordance:", maximum_value)
cat(sprintf("\nOptimal clusters number: %d", optimal_clusters_number))

cat("\n eta ", e_opt)
cat("\n tau ", t_opt)
cat("\n mu ", m_opt)
cat("\n k ", k_opt)

data <- c(e_opt, t_opt, m_opt, k_opt, current_concordance)
write(data, file = "optimal_parameters.txt", append = FALSE, sep = " // ")