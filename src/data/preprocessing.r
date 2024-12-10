preprocess <- function(dataset) {
    cat("\nPreprocessing data")

    #TODO verificare come funziona s2gc
    dataset$radiomics_pre <- preprocess_radiomics(dataset$radiomics_pre)
    dataset$radiomics_post <- preprocess_radiomics(dataset$radiomics_post)
    dataset$pre_operative <- preprocess_pre_operative(dataset$pre_operative)

    return(dataset)
}

preprocess_radiomics <- function(features) {  
    features <- scale(features)

    return(features)
}

preprocess_pre_operative <- function(features) {
    age <- features[, 1]
    age <- scale(age)
    features[, 1] <- age

    mtx <- features[, 2]
    mtx <- scale(mtx)
    features[, 2] <- mtx

    return(features)
}