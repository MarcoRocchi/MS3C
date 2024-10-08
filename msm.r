Sys.setenv(LANG = "en")

source("./src/load_data.r")
source("./src/preprocessing/preprocessing_pipeline.r")
#source("./src/models/model_h0.r")
source("./src/models/model1.r")
source("./src/models/model2.r")
source("./src/models/model3.r")
source("./src/models/model4.r")
source("./src/models/model5.r")

dataset <- load_data()

dataset$radiomics_pre <- preprocess(dataset$radiomics_pre)
dataset$radiomics_post <- preprocess(dataset$radiomics_post)
dataset$pre_operative <- preprocess_pre_operative(dataset$pre_operative)

#model_h0 <- build_h0(
#    dataset$radiomics_pre,
#    dataset$relapse_status,
#    dataset$relapse_times
#)

model1 <- build_m1(
    dataset$radiomics_pre,
    dataset$radiomics_post,
    dataset$pre_operative,
    dataset$dead_status,
    dataset$post_times,
    dataset$dead_times
)

model2 <- build_m2(
    dataset$radiomics_pre,
    dataset$radiomics_post,
    dataset$pre_operative,
    dataset$relapse_status,
    dataset$post_times,
    dataset$relapse_times
)

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

cat("End")