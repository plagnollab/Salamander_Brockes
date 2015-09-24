## D0  baseline
## D15 15 days
## D150 150 days

all.results <- list.files("results_std_deseq/", pattern = "*_fullresults.csv", full.names = TRUE)

for (file in all.results) {
  message("File ", file)
  data <- read.csv(file, stringsAsFactors = FALSE)
  data <- data[ order(data$pval), ]
  print(head(data))
}
