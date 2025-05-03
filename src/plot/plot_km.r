library(survival)
library(survminer)
library(ggsurvfit)

plot_km_curve <- function(dataset, group) {
    d <- dataset[[2]]

    order <- d$patients
    features <- d$features[order, ]
    times <- d$times[order]
    status <- d$censoring[order]

    legend_labels <- paste("Cluster", seq_len(length(unique(group))))

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
        legend.labs = legend_labels
    )

    print(survp)
}
