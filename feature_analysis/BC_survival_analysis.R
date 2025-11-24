# install.packages("openxlsx")
# install.packages("survminer")

library(dplyr)
library(openxlsx)
library(survival)
library(survminer)

getwd()
setwd("/Users/yuzimeng/Desktop/Yale/courses_fall_2025/IHI/final")

filtered_df <- read.csv("filtered_data.csv")
all_features <- read.csv("all_features.csv")$x

# feature_df <- read.xlsx("selected_features.xlsx")
# 
# time <- filtered_df$overall_survival_months
# event <- filtered_df$overall_survival
# 
# clin_feature <- feature_df$Clinical
# exp_feature <- feature_df$GeneExpression
# mut_feature <- feature_df$Mutation
# comb_feature <- feature_df$Combined
# 
# 
# all_features <- unique(c(clin_feature, exp_feature, mut_feature, comb_feature))

# replace_dict <- c(
#   "type_of_breast_surgery_MASTECTOMY" = "type_of_breast_surgery",
#   "pam50_+_claudin-low_subtype_LumA" = "pam50_+_claudin-low_subtype",
#   "primary_tumor_laterality_Left" = "primary_tumor_laterality",
#   "primary_tumor_laterality_Right" = "primary_tumor_laterality",
#   "neoplasm_histologic_grade_3.0" = "neoplasm_histologic_grade",
#   "3-gene_classifier_subtype_ER+/HER2- Low Prolif" = "3-gene_classifier_subtype"
# )
# 
# all_features <- dplyr::recode(all_features, !!!replace_dict)
# all_features <- unique(all_features)
# # all_features <- all_features %>% drop_na()
# # write.csv(all_features, "all_features.csv", row.names = FALSE)
# cat("Total unique features to test:", length(all_features), "\n\n")

feature_name <- "X3.gene_classifier_subtype"
x <- filtered_df[[feature_name]]

cat("Data type:", class(x), "\n")
cat("Unique values:", length(unique(na.omit(x))), "\n\n")

SurvOS_full <- Surv(time = filtered_df$overall_survival_months,
                    event = filtered_df$overall_survival == 0)

###### ---------- is numerical ----------- #######

if (is.numeric(x)) {
  quartiles <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  group <- cut(x, breaks = quartiles, 
               labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)"), 
               include.lowest = TRUE)
  plot_title <- paste0(feature_name, " (Quartile groups)")
  
} else {
  # Categorical / Binary
  if (is.character(x)) {
    x <- trimws(x)                  # remove leading/trailing whitespace
    x[x == ""] <- NA                # convert empty string "" → real NA
  }
  group <- as.factor(x)
  
  # Now NAs are properly marked
  cat("   After cleaning empty strings:\n")
  print(table(group, useNA = "always"))
  
  # Optional: drop NA only for plotting (recommended for clean figures)
  plot_data <- filtered_df[!is.na(group), ]   # only complete cases for plot
  group_plot <- droplevels(group[!is.na(group)])
  
  if (nlevels(group_plot) > 12) {
    cat("   Too many levels → skipping\n")
    return(invisible(NULL))
  }
  
  plot_title <- paste0(feature_name, " (n=", nrow(plot_data), " complete)")
}

# Show distribution (including NA)
print(table(group, useNA = "always"))
# Show group distribution by counts
print(table(group, useNA = "ifany"))


####### --------- Survival Analysis Model ------------- #######
Surv_clean <- Surv(time = plot_data$overall_survival_months,
                   event = plot_data$overall_survival == 0)
fit <- survfit(Surv_clean ~ group_plot)

# # Log-rank test
logrank <- survdiff(Surv_clean ~ group_plot)
p_value <- 1 - pchisq(logrank$chisq, df = nlevels(group_plot) - 1)
p_logrank <- 1 - pchisq(survdiff(Surv_clean ~ group_plot)$chisq, df = nlevels(group_plot)-1)

###### ---------- If only 2 groups → also compute HR ----------- #######

hr_text <- "Multi-group (HR not shown)"

if (nlevels(group_plot) == 2) {
  cox_model <- coxph(Surv_clean ~ group_plot, data = plot_data)
  
  # Extract coefficients safely
  beta   <- coef(cox_model)           # named vector
  hr     <- round(exp(beta), 3)       # HR
  ci     <- round(exp(confint(cox_model)), 3)  # 95% CI matrix
  p_val  <- round(summary(cox_model)$coefficients[5], 4)  # p-value column 5
  
  # Build clean text
  ref_level  <- levels(group_plot)[1]
  risk_level <- levels(group_plot)[2]
  
  hr_text <- paste0("HR (", risk_level, " vs ", ref_level, ") = ", hr,
                    "\n95% CI: ", ci[1], "–", ci[2], 
                    ", p = ", ifelse(p_val < 0.001, "<0.001", p_val))
} else {
  hr_text <- "Multi-group (HR not shown)"
}

###### ---------- Plot  ----------- #######

p <- ggsurvplot(
  fit,
  data = plot_data,
  title = plot_title,
  pval = TRUE,
  pval.coord = c(5, 0.1),
  risk.table = TRUE,
  risk.table.height = 0.3,
  legend.title = feature_name,
  legend.labs = levels(group_plot),
  palette = "jco",
  xlab = "Overall Survival (months)",
  ylab = "Survival Probability",
  ggtheme = theme_bw(base_size = 14)
)

print(p)

png_file <- paste0("KM_plots/KM_", feature_name, "_p", signif(p_logrank, 3), ".png")
ggsave(png_file, plot = p$plot, width = 9, height = 7, dpi = 300, bg = "white")

cat("\nLog-rank p-value:", signif(p_logrank, 3),
    ifelse(p_logrank < 0.05, " → SIGNIFICANT!\n", "\n"))
cat(hr_text, "\n")

