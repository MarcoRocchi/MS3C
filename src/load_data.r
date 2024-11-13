library("readxl")

source("./src/utils/features.r")

rename <- function(df, names, prefix = "") {
    for (i in 1:dim(df)[2]) {
        name <- colnames(df)[i]
        new_name <- paste0(prefix, names[name])
        colnames(df)[i] <- new_name
    }

    return(df)
}

load_data <- function() {
    #TODO embedding
    message("Loading data")

    radiomics_pre <- read_excel("../Data/V2/2 - Selected features/1 - Pre-chemo.xlsx")
    radiomics_pre <- rename(radiomics_pre[, -1], names.radiomics_pre, "")
    radiomics_pre <- as.matrix(radiomics_pre)

    patients_count <- nrow(radiomics_pre)

    radiomics_post <- read_excel("../Data/V2/2 - Selected features/2 - Post-chemo.xlsx")
    radiomics_post <- rename(radiomics_post[, -1], names.radiomics_post, "")
    radiomics_post <- as.matrix(radiomics_post)

    pre_operative <- read_excel("../Data/V2/1 - Cleaned/3 - Pre operative.xlsx")
    pre_operative["synchronous"] <- ifelse(pre_operative["Sinc 1/Meta 2"] == 0, 1, 0)
    pre_operative["metachronous"] <- ifelse(pre_operative["Sinc 1/Meta 2"] == 1, 1, 0)
    pre_operative["M"] <- ifelse(pre_operative["SEX"] == 1, 1, 0)
    pre_operative["F"] <- ifelse(pre_operative["SEX"] == 0, 1, 0)
    pre_operative["no"] <- ifelse(pre_operative["Malattia extrahep sinc fegato (0/1)"] == 0, 1, 0)
    pre_operative["sì"] <- ifelse(pre_operative["Malattia extrahep sinc fegato (0/1)"] == 1, 1, 0)
    pre_operative <- rename(pre_operative[, -1], names.pre_operative)
    pre_operative <- pre_operative[names.pre_operative]
    pre_operative <- as.matrix(pre_operative)

    surgery <- read_excel("../Data/V2/1 - Cleaned/4 - Surgery.xlsx")
    surgery["morb severa sì"] <- ifelse(surgery["Morb severa"] == 1, 1, 0)
    surgery["morb severa no"] <- ifelse(surgery["Morb severa"] == 0, 1, 0)
    surgery["Infective morbidity sì"] <- ifelse(surgery["Infective morbidity (0/1)"] == 1, 1, 0)
    surgery["Infective morbidity no"] <- ifelse(surgery["Infective morbidity (0/1)"] == 0, 1, 0)
    surgery["r0 = Margine almeno 1 mm sì"] <- ifelse(surgery["r0 = Margine almeno 1 mm"] == 1, 1, 0)
    surgery["r0 = Margine almeno 1 mm no"] <- ifelse(surgery["r0 = Margine almeno 1 mm"] == 0, 1, 0)
    surgery <- rename(surgery[, -1], names.surgery)
    surgery <- surgery[names.surgery]
    surgery <- as.matrix(surgery)

    times <- read_excel("../Data/V2/1 - Cleaned/5 - Times.xlsx")
    dfs <- read_excel("../Data/V2/1 - Cleaned/6 - DFS.xlsx")
    os <- read_excel("../Data/V2/1 - Cleaned/7 - OS.xlsx")

    relapse_status <- as.matrix(dfs["Recidiva (0/1)"])
    dead_status <- as.matrix(os["Morto =1"])

    relapse_times <- as.matrix((times["Delta pre post"] + times["Delta post surgery"] + dfs["DFS"]) / 12)
    dimnames(relapse_times) <- list(NULL, names.relapse_time)
    dead_times <- as.matrix((times["Delta pre post"] + times["Delta post surgery"] + os["OS"]) / 12)
    dimnames(dead_times) <- list(NULL, names.dead_time)

    post_times <- as.matrix((times["Delta pre post"]) / 12)
    dimnames(post_times) <- list(NULL, names.post_time)
    surgery_times <- as.matrix((times["Delta pre post"] + times["Delta post surgery"]) / 12)
    dimnames(surgery_times) <- list(NULL, names.surgery_time)

    return(list(
        patients_count = patients_count,
        radiomics_pre = radiomics_pre, 
        radiomics_post = radiomics_post,
        pre_operative = pre_operative,
        surgery = surgery,
        relapse_status = relapse_status,
        dead_status = dead_status,
        relapse_times = relapse_times,
        dead_times = dead_times,
        post_times = post_times,
        surgery_times = surgery_times)
    )
}