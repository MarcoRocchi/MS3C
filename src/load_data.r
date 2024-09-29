source("./src/preprocessing/standardize.r")

library("readxl")

#TODO non mischiare funzionalità logiche

#TODO rimuovere
load_data_s2gc <- function() {
    data1 <- as.matrix(read.csv("./data/BREAST_Gene_Expression.csv"))
    data2 <- as.matrix(read.csv("./data/BREAST_Methy_Expression.csv"))
    data3 <- as.matrix(read.csv("./data/BREAST_Mirna_Expression.csv"))
    
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

load_data_pre <- function() {
    data <- read_excel("../../Data/V2/2 - Selected features/1 - Pre-chemo.xlsx")
    radiomics <- as.matrix(data[, -1])

    radiomics <- standardize_matrix(matrix(as.double(radiomics), nrow = 99))
    radiomics <- standardize_columns(radiomics)

    times_data <- read_excel("../../Data/V2/1 - Cleaned/5 - Times.xlsx")
    survival_data <- read_excel("../../Data/V2/1 - Cleaned/6 - DFS.xlsx")

    responses <- as.matrix(survival_data["Recidiva (0/1)"])
    times <- as.matrix(times_data["Delta pre post"] + times_data["Delta post surgery"] + survival_data["DFS"])

    return(list(features = radiomics, times = times, responses = responses))
}

load_data_post <- function() {
    data <- read_excel("../../Data/V2/2 - Selected features/2 - Post-chemo.xlsx")
    radiomics <- as.matrix(data[, -1])

    radiomics <- standardize_matrix(matrix(as.double(radiomics), nrow = 99))
    radiomics <- standardize_columns(radiomics)

    times_data <- read_excel("../../Data/V2/1 - Cleaned/5 - Times.xlsx")
    survival_data <- read_excel("../../Data/V2/1 - Cleaned/6 - DFS.xlsx")

    responses <- as.matrix(survival_data["Recidiva (0/1)"])
    times <- as.matrix(times_data["Delta post surgery"] + survival_data["DFS"])

    return(list(features = radiomics, times = times, responses = responses))
}