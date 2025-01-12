logrank <- function(dataset, groups) {
    
    dataset$group <- groups[dataset$id]
    
    model <- coxph(Surv(time, status) ~ group + strata(trans),
        data = dataset
    )

    return(summary(model))
}
