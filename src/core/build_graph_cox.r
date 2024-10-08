source("./src/core/calculate_graph_weights.r")

build_graph_cox <- function(features, responses, times) {
    responses <- as.data.frame(responses)
    times <- as.data.frame(times)

    data <- cbind(features, responses, times)

    features_names <- paste(colnames(features), collapse = " + ")

    f <- as.formula(paste("Surv(times, responses) ~ ", features_names))

    model <- coxph(f, data = data)
    
    summary(model)

    k <- 10
    #graph <- CalGraphWeight(features, ws, k)

    #return(list(s = graph$SS))
}