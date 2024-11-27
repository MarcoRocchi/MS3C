preprocess <- function(dataset) {
    cat("\nPreprocessing data")

    #TODO verificare come funziona s2gc
    dataset$radiomics_pre <- t(preprocess_radiomics(t(dataset$radiomics_pre)))
    dataset$radiomics_post <- t(preprocess_radiomics(t(dataset$radiomics_post)))
    dataset$pre_operative <- t(preprocess_pre_operative(t(dataset$pre_operative)))

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