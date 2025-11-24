## this is an example for survival analysis @Zimeng Yu

library(dplyr)
library(openxlsx)
library(survival)
library(survminer)

getwd()
setwd("/Users/folder")

filtered_df <- read.csv("data.csv")
all_features <- read.csv("all_features.csv")

source("km_analysis_functions.R")

all_features <- all_features$x


### for one feature use
result <- km_one_feature("X3.gene_classifier_subtype")
print(result)

### for all features use
all_results <- lapply(all_features, function(feat) {
  tryCatch({
    km_one_feature(feat)
  }, error = function(e) {
    cat("Error with feature:", feat, "-", e$message, "\n")
    return(NULL)
  })
})


### Combine results into one dataframe
results_df <- do.call(rbind, all_results[!sapply(all_results, is.null)])
print(results_df)

### Save results
write.csv(results_df, "KM_analysis_results.csv", row.names = FALSE)
