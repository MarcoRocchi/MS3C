library(survival)
library(survminer)
library(ggsurvfit)

plot_km_curve <- function(dataset, group) {
    d <- dataset[[2]]

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
        conf.int.style = "step",
        conf.int = FALSE,
        censor = FALSE,
        title = "Post chemo -> Relapse",
        legend = "bottom",
        legend.title = "",
        legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3")
    )

    print(survp)
}
