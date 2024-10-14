library("readxl")

gt <- as.matrix(read_excel("breast-gt.xlsx", col_names = FALSE))
output <- as.matrix(read_excel("breast.xlsx"))

diff <- gt - output

cat("Maximum error: ", max(unlist(diff)))