concordance_index <- function(data, weights, patients_count) {
    concordance <- 0

    for (i in 1:length(data)) {
        d <- data[[i]]
        concordant <- 0
        discordant <- 0
        
        risks <- d$features %*% weights[[i]]$z

        for (a in 1:patients_count) {
            for (b in a:patients_count) {
                if (a != b) {
                    if (d$censoring[a] == 0 && d$censoring[b] == 0) {
                        #Not censored
                        if ((risks[a] > risks[b] && d$times[a] < d$times[b]) || 
                            (risks[a] < risks[b] && d$times[a] > d$times[b])) {
                            concordant <- concordant + 1
                        } else {
                            discordant <- discordant + 1
                        }
                    } else {
                        if (d$censoring[a] == 0) {
                            #Individual in position a not censored
                            if (d$times[b] > d$times[a]) {
                                if (risks[a] > risks[b]) {
                                    concordant <- concordant + 1
                                } else {
                                    discordant <- discordant + 1
                                }   
                            }
                        } else if (d$censoring[b] == 0) {
                            #Individual in position b not censored
                            if (d$times[a] > d$times[b]) {
                                if (risks[b] > risks[a]) {
                                    concordant <- concordant + 1
                                } else {
                                    discordant <- discordant + 1
                                }   
                            }
                        }
                    }
                }
            }
        }

        c <- concordant / (concordant + discordant)
        concordance <- concordance + c
    }

    return(concordance)
}