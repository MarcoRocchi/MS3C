library("readxl")

source("./src/data/features.r")

rename <- function(df, names, prefix = "") {
    for (i in 1:dim(df)[2]) {
        name <- colnames(df)[i]
        new_name <- paste0(prefix, names[name])
        colnames(df)[i] <- new_name
    }

    return(df)
}

load_data <- function() {
    cat("\nLoading data")

    radiomics_pre <- read_excel("../Data/V2/2 - Selected features/1 - Pre-chemo.xlsx")
    radiomics_pre <- rename(radiomics_pre[, -1], names.radiomics_pre, "")
    radiomics_pre <- as.matrix(radiomics_pre)

    patients_count <- nrow(radiomics_pre)

    radiomics_post <- read_excel("../Data/V2/2 - Selected features/2 - Post-chemo.xlsx")
    radiomics_post <- rename(radiomics_post[, -1], names.radiomics_post, "")
    radiomics_post <- as.matrix(radiomics_post)

    pre_operative <- read_excel("../Data/V2/1 - Cleaned/3 - Pre operative.xlsx")
    pre_operative["Sinc 1/Meta 2_1"] <- ifelse(pre_operative["Sinc 1/Meta 2"] == 0, 0.05255485, -0.39586723)
    pre_operative["Sinc 1/Meta 2_2"] <- ifelse(pre_operative["Sinc 1/Meta 2"] == 0, 0.17627974, 0.28472582)
    pre_operative["SEX_1"] <- ifelse(pre_operative["SEX"] == 1, -0.25222108, -0.21486229)
    pre_operative["SEX_2"] <- ifelse(pre_operative["SEX"] == 1, -0.96437500, -0.77958780)
    pre_operative["extrahep_1"] <- ifelse(pre_operative["Malattia extrahep sinc fegato (0/1)"] == 0, 1.43862840, 1.86265170)
    pre_operative["extrahep_2"] <- ifelse(pre_operative["Malattia extrahep sinc fegato (0/1)"] == 0, 0.20406507, 0.22706775)
    pre_operative <- rename(pre_operative[, -1], names.pre_operative)
    pre_operative <- pre_operative[names.pre_operative]
    pre_operative <- as.matrix(pre_operative)

    surgery <- read_excel("../Data/V2/1 - Cleaned/4 - Surgery.xlsx")
    surgery["Morb severa_1"] <- ifelse(surgery["Morb severa"] == 1, -0.45834910, -0.36790693)
    surgery["Morb severa_2"] <- ifelse(surgery["Morb severa"] == 1, -0.77975430, -0.65975250)
    surgery["Infective morbidity_1"] <- ifelse(surgery["Infective morbidity (0/1)"] == 1, -0.43539820, 0.09909006)
    surgery["Infective morbidity_2"] <- ifelse(surgery["Infective morbidity (0/1)"] == 1, 0.84917235, 0.65826590)
    surgery["r0 = Margine almeno 1 mm_1"] <- ifelse(surgery["r0 = Margine almeno 1 mm"] == 1, 0.61577220, 0.75887620)
    surgery["r0 = Margine almeno 1 mm_2"] <- ifelse(surgery["r0 = Margine almeno 1 mm"] == 1, -0.18973102, -0.18494628)
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
    surgery_times <- as.matrix((times["Delta pre post"] + times["Delta post surgery"]) / 12)
    dimnames(surgery_times) <- list(NULL, names.surgery_time)
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
        surgery_times = surgery_times,
        dead_times = dead_times,
        post_times = post_times,
        surgery_times = surgery_times)
    )
}