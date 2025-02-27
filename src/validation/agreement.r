# Computes the concordance index for clustering agreement with survival data
# T: Numeric vector of survival times
# E: Binary event indicator (0 = event observed, 1 = censored)
# C: Integer vector of cluster assignments
compute_survival_cindex <- function(T, E, C) {
    n <- length(T)
    concordant <- 0
    comparable <- 0

    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            if (E[i] == 0 && E[j] == 0) {
                comparable <- comparable + 1
                if ((T[i] < T[j] && C[i] < C[j]) || (T[i] > T[j] && C[i] > C[j])) {
                    concordant <- concordant + 1
                } else if (T[i] == T[j]) {
                    concordant <- concordant + 0.5
                }
            } else {
                if (E[i] == 0) {
                    if (T[j] > T[i]) {
                        comparable <- comparable + 1
                        if (C[j] > C[i]) {
                            concordant <- concordant + 1
                        } else {
                            if (C[i] == C[j]) {
                                concordant <- concordant + 0.5
                            }
                        }   
                    }
                } else if (E[j] == 0) {
                    if (T[i] > T[j]) {
                        comparable <- comparable + 1
                        if (C[i] > C[j]) {
                            concordant <- concordant + 1
                        } else {
                            if (C[i] == C[j]) {
                                concordant <- concordant + 0.5
                            }
                        }   
                    }
                }
            }
        }
    }

    c_index <- concordant / comparable
    return(c_index)
}
