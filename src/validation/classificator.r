library(caret)
library(pROC)

compute_classificator <- function(data, target, classes) {
    auc_values <- c()

    for (cluster in 1:classes) {
        labels <- factor(ifelse(target == cluster, "Yes", "No"), levels = c("No", "Yes"))

        data <- data.frame(data, label = labels)

        train_control <- trainControl(
            method = "cv", 
            number = 5, 
            savePredictions = TRUE, 
            classProbs = TRUE, 
            summaryFunction = twoClassSummary
        )

        model <- train(
            label ~ ., 
            data = data, 
            method = "glm", 
            family = "binomial", 
            trControl = train_control, 
            metric = "ROC"
        )

        roc_curve <- roc(model$pred$obs, model$pred$Yes)
        auc <- auc(roc_curve)
        auc_values <- c(auc_values, auc)
    }

    average_auc <- mean(auc_values)
    minimum_auc <- min(auc_values)
    cat(sprintf("\nAverage AUC: %f", average_auc))
    cat(sprintf("\nMinimum AUC: %f", minimum_auc))
}