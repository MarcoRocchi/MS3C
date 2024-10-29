Sys.setenv(LANG = "en")

source("./src/load_data.r")
source("./src/preprocessing/preprocessing_pipeline.r")
source("./src/models/model1.r")
source("./src/models/model2.r")
source("./src/models/model3.r")
source("./src/models/model4.r")
source("./src/models/model5.r")
source("./src/models/model6.r")
source("./src/models/model7.r")

dataset <- load_data()

dataset$radiomics_pre <- preprocess(dataset$radiomics_pre)
dataset$radiomics_post <- preprocess(dataset$radiomics_post)
dataset$pre_operative <- preprocess_pre_operative(dataset$pre_operative)

build <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)

if (build[1]) {
    model1 <- build_m1(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative,
        dataset$dead_status,
        dataset$post_times,
        dataset$dead_times
    )

    cat("AIC", AIC(model1), "\n")
    cat("BIC", BIC(model1), "\n")
}

if (build[2]) {
    model2 <- build_m2(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative,
        dataset$relapse_status,
        dataset$post_times,
        dataset$relapse_times
    )

    cat("AIC", AIC(model2), "\n")
    cat("BIC", BIC(model2), "\n")
}

if (build[3]) {
    model3 <- build_m3(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative,
        dataset$relapse_status,
        dataset$dead_status,
        dataset$post_times,
        dataset$relapse_times,
        dataset$dead_times
    )

    cat("AIC", AIC(model3), "\n")
    cat("BIC", BIC(model3), "\n")
}

if (build[4]) {
    model4 <- build_m4(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative,
        dataset$surgery,
        dataset$relapse_status,
        dataset$post_times,
        dataset$surgery_times,
        dataset$relapse_times
    )

    cat("AIC", AIC(model4), "\n")
    cat("BIC", BIC(model4), "\n")
}

if (build[5]) {
    model5 <- build_m5(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative,
        dataset$surgery,
        dataset$relapse_status,
        dataset$dead_status,
        dataset$post_times,
        dataset$surgery_times,
        dataset$relapse_times,
        dataset$dead_times
    )

    cat("AIC", AIC(model5), "\n")
    cat("BIC", BIC(model5), "\n")
}

if (build[6]) {
    model6 <- build_m6(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative,
        dataset$relapse_status,
        dataset$dead_status,
        dataset$post_times,
        dataset$relapse_times,
        dataset$dead_times
    )

    cat("AIC", AIC(model6), "\n")
    cat("BIC", BIC(model6), "\n")
}

if (build[7]) {
    model7 <- build_m7(
        dataset$radiomics_pre,
        dataset$radiomics_post,
        dataset$pre_operative,
        dataset$surgery,
        dataset$relapse_status,
        dataset$dead_status,
        dataset$post_times,
        dataset$surgery_times,
        dataset$relapse_times,
        dataset$dead_times
    )

    cat("AIC", AIC(model7), "\n")
    cat("BIC", BIC(model7), "\n")
}

cat("End")