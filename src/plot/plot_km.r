library(survival)
library(survminer)
library(ggsurvfit)

plot_km_curve <- function(dataset, group) {
    d <- dataset[[5]]

    order <- d$patients
    features <- d$features[order, ]
    times <- d$times[order]
    status <- d$censoring[order]

    model <- survfit(Surv(times, status) ~ group, data = as.data.frame(cbind(features, group)))
    survp <- ggsurvplot(
        model, 
        risk.table = FALSE, 
        data = as.data.frame(cbind(features, times, status, group)),
        xlab = "Time [years]",
        ylab = "Survival probability",
        conf.int.style = "step",
        conf.int = FALSE,
        censor = FALSE,
        title = "Relapse -> Dead",
        legend = "bottom",
        legend.title = "",
        legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5",
                        "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10",
                        "Cluster 11", "Cluster 12", "Cluster 13", "Cluster 14", "Cluster 15"
        )
    )

    print(survp)
}
