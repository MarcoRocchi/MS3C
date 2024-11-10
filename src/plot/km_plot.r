library(survival)

#TODO generalizzare
plot_clusters <- function(features, times, responses, clusters, xlabel, ylabel) {
    data <- as.data.frame(cbind(features, times, responses, clusters))

    #TODO verificare formula
    model <- survfit(Surv(times, responses) ~ clusters, data = data)
    plot(model, conf.int = FALSE, xlab = xlabel, ylab = ylabel,
        col = c("blue", "green"))
    
    legend(
        x = "topright", 
        fill = c("blue","green"),
        legend = c("Cluster 1", "Cluster 2")
    )
}