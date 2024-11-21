insert_missing_patients <- function(d, patients_count) {
    patients <- as.list(d["id"])$id

    for (j in 1:patients_count) {
        if (!j %in% patients) {
            new_patient <- setNames(numeric(ncol(d)), names(d))
            new_patient["id"] <- j
            d <- rbind(d, new_patient)
        }
    }

    return(arrange(d, id))
}