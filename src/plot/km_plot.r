library(survival)

#TODO generalizzare
plot_clusters <- function(features, times, responses, clusters) {
    data <- as.data.frame(cbind(features, times, responses, clusters))

    #TODO verificare formula
    model <- survfit(Surv(times, responses) ~ clusters, data = data)
    plot(model, conf.int = FALSE, xlab = "Time (weeks)", ylab = "p",
        col = c("blue", "red", "green", "orange"))
    
    legend(
        x = "topright", 
        fill = c("blue", "red", "green", "orange"),
        legend = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
    )
}