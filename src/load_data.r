source("./src/preprocessing/standardize.r")

#TODO rimuovere
load_data <- function() {
    data1 <- as.matrix(read.csv("./data/BREAST_Gene_Expression.csv"))
    data2 <- as.matrix(read.csv("./data/BREAST_Methy_Expression.csv"))
    data3 <- as.matrix(read.csv("./data/BREAST_Mirna_Expression.csv"))
    
    #TODO non mischiare funzionalità logiche
    data1 <- standardize_matrix(matrix(as.double(t(data1[, -1])), nrow = 105))
    data1 <- standardize_columns(data1)
    data2 <- standardize_matrix(matrix(as.double(t(data2[, -1])), nrow = 105))
    data2 <- standardize_columns(data2)
    data3 <- standardize_matrix(matrix(as.double(t(data3[, -1])), nrow = 105))
    data3 <- standardize_columns(data3)

    features <- cbind(data1, data2, data3)

    y_data <- read.csv("./data/BREAST_Survival.csv")
    times <- as.matrix(y_data["Survival"])
    responses <- as.matrix(y_data["Death"])

    return(list(features = features, times = times, responses = responses))
}