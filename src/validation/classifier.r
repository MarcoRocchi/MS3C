library(caret)
library(pROC)

compute_classifier <- function(data, target) {
    labels <- factor(make.names(target))

    data <- data.frame(data, label = labels)

    train_control <- trainControl(
        method = "cv", 
        number = 5, 
        savePredictions = TRUE, 
        classProbs = TRUE, 
        summaryFunction = multiClassSummary
    )

    model <- train(
        label ~ ., 
        data = data, 
        method = "multinom", 
        trControl = train_control, 
        metric = "ROC"
    )

    predictions <- predict(model, data, type = "prob")
    roc_curves <- multiclass.roc(labels, predictions)
    auc <- auc(roc_curves)
    
    cat(sprintf("\nAUC: %f", auc))

    return(auc)
}