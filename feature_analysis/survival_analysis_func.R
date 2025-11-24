library(tidyr)
library(survival)
library(survminer)
library(dplyr)

getwd()
setwd("/Users/yuzimeng/Desktop/Yale/courses_fall_2025/IHI/final")

filtered_df <- read.csv("filtered_data.csv")
all_features <- read.csv("all_features.csv")
all_features <- all_features$x


SurvOS_full <- Surv(time = filtered_df$overall_survival_months,
                    event = filtered_df$overall_survival == 0)

km_one_feature <- function(feature_name, data = filtered_df) {
  
  cat("\n=== Analyzing:", feature_name, "===\n")
  x <- data[[feature_name]]
  
  if (is.null(x) || all(is.na(x)) || length(unique(na.omit(x))) < 2) {
    cat("   Skipped: no variation or all NA\n")
    return(invisible(NULL))
  }
  
  # ——— Force known ordinal/grade variables to categorical ———
  grade_vars <- c("neoplasm_histologic_grade", "NEOPLASM_HISTOLOGIC_GRADE",
                  "cellularity", "mitotic_index", "tubule_formation",
                  "nuclear_pleomorphism", "histological_grade", "GRADE")
  if (feature_name %in% grade_vars) {
    cat("   → Treating as categorical (clinical grade)\n")
    x <- as.factor(x)
  }
  
  # ——— Initialize variables that will be used later ———
  plot_data <- data
  group_plot <- NULL
  plot_title <- feature_name
  
  # ——— Grouping logic ———
  if (is.numeric(x) && !is.factor(x)) {
    # Truly continuous → quartiles
    q <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    if (length(unique(q)) < 4) {
      cat("   Low spread → median split\n")
      group_plot <- ifelse(x <= median(x, na.rm = TRUE), "Low", "High")
    } else {
      group_plot <- cut(x, breaks = q,
                        labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)"),
                        include.lowest = TRUE)
    }
    plot_title <- paste0(feature_name, " (Quartiles)")
    
    # Remove NA from numeric grouping
    valid <- !is.na(group_plot)
    plot_data <- data[valid, ]
    group_plot <- droplevels(group_plot[valid])
    
  } else {
    # Categorical (including forced grades)
    if (is.character(x)) {
      x <- trimws(x)
      x[x == ""] <- NA
    }
    group <- as.factor(x)
    cat("   After cleaning:\n")
    print(table(group, useNA = "always"))
    
    # Remove NA for plotting
    valid <- !is.na(group)
    plot_data <- data[valid, ]
    group_plot <- droplevels(group[valid])
    
    if (nlevels(group_plot) == 0 || nrow(plot_data) < 10) {
      cat("   Too few complete cases → skipping\n")
      return(invisible(NULL))
    }
    if (nlevels(group_plot) > 12) {
      cat("   Too many levels → skipping\n")
      return(invisible(NULL))
    }
    
    # Nice labels for grade
    if (feature_name %in% grade_vars) {
      levels(group_plot) <- paste("Grade", levels(group_plot))
      plot_title <- "Histologic Grade"
    } else {
      plot_title <- paste0(feature_name, " (n=", nrow(plot_data), ")")
    }
  }
  
  # Final safety check
  if (nlevels(group_plot) < 1) {
    cat("   No valid groups after cleaning → skipping\n")
    return(invisible(NULL))
  }
  
  # ——— Survival object on clean data ———
  Surv_clean <- Surv(time = plot_data$overall_survival_months,
                     event = plot_data$overall_survival == 0)
  
  # ——— KM fit ———
  fit <- survfit(Surv_clean ~ group_plot)
  
  # ——— Log-rank p-value ———
  p_logrank <- try({
    test <- survdiff(Surv_clean ~ group_plot)
    1 - pchisq(test$chisq, df = nlevels(group_plot) - 1)
  }, silent = TRUE)
  if (inherits(p_logrank, "try-error")) p_logrank <- NA
  
  
  ###### ——— HR calculation: binary → single HR, multi-group → pairwise HRs ——— ######
  hr_text <- "Multi-group (HR not shown)"   # default for plot annotation
  
  if (nlevels(group_plot) == 2) {
    # Binary: single HR
    cox_fit <- try(coxph(Surv_clean ~ group_plot, data = plot_data), silent = TRUE)
    if (!inherits(cox_fit, "try-error") && nrow(summary(cox_fit)$coefficients) > 0) {
      hr_val <- round(exp(coef(cox_fit)), 2)
      ci_val <- round(exp(confint(cox_fit)), 2)
      p_cox  <- summary(cox_fit)$coefficients[5]
      p_text <- ifelse(p_cox < 0.001, "<0.001", round(p_cox, 3))
      
      ref  <- levels(group_plot)[1]
      risk <- levels(group_plot)[2]
      
      hr_text <- paste0("HR (", risk, " vs ", ref, ") = ", hr_val,
                        "\n95% CI: ", ci_val[1], "–", ci_val[2], ", p = ", p_text)
    } else {
      hr_text <- "Cox model failed"
    }
  } 
  else if (nlevels(group_plot) > 2) {
    # Multi-group: pairwise vs first level
    cox_multi <- try(coxph(Surv_clean ~ group_plot, data = plot_data), silent = TRUE)
    if (!inherits(cox_multi, "try-error")) {
      summ <- summary(cox_multi)
      ref_level <- levels(group_plot)[1]
      
      hr_lines <- sapply(2:nlevels(group_plot), function(i) {
        hr  <- round(exp(coef(cox_multi))[i-1], 2)
        ci  <- round(exp(confint(cox_multi))[i-1, ], 2)
        p   <- summ$coefficients[i-1, 5]
        p_txt <- ifelse(p < 0.001, "<0.001", round(p, 3))
        paste0(levels(group_plot)[i], " vs ", ref_level, 
               ": HR = ", hr, " (", ci[1], "–", ci[2], ", p = ", p_txt, ")")
      })
      hr_text <- paste("Pairwise HRs (vs", ref_level, "):", 
                       paste(hr_lines, collapse = "\n"), sep = "\n")
    } else {
      hr_text <- "Cox model failed (multi-group)"
    }
  }
  
  
  # ——— Plot ———
  # p <- ggsurvplot(
  #   fit, data = plot_data,
  #   title = plot_title,
  #   pval = TRUE, pval.coord = c(5, 0.1),
  #   risk.table = TRUE, risk.table.height = 0.3,
  #   legend.title = feature_name,
  #   legend.labs = levels(group_plot),
  #   palette = "jco",
  #   xlab = "Overall Survival (months)",
  #   ggtheme = theme_bw(base_size = 14)
  # )
  
  # print(p)
  
  # safe_name <- gsub("[^A-Za-z0-9_]", "_", feature_name)
  # png_file <- paste0("KM_plots/KM_", safe_name, "_p", signif(p_logrank, 3), ".png")
  # ggsave(png_file, plot = p$plot, width = 9.5, height = 7.5, dpi = 300, bg = "white")
  # cat("   Saved: ", png_file, "\n")
  
  cat("Log-rank p =", ifelse(is.na(p_logrank), "NA", signif(p_logrank, 3)),
      ifelse(!is.na(p_logrank) && p_logrank < 0.05, " → SIGNIFICANT!\n", "\n"))
  
  return(data.frame(
    Feature = feature_name,
    Type = ifelse(is.numeric(data[[feature_name]]), "Numeric", "Categorical"),
    Groups = paste(levels(group_plot), collapse = " | "),
    N = nrow(plot_data),
    Logrank_p = ifelse(is.na(p_logrank), NA, p_logrank),
    HR_info = hr_text,                              # Now contains full pairwise info!
    Significant = !is.na(p_logrank) && p_logrank < 0.05,
    stringsAsFactors = FALSE
  ))
}

results_list <- list()
for (feat in all_features) {
  result <- km_one_feature(feat)
  if (!is.null(result)) results_list[[feat]] <- result
  Sys.sleep(0.5)
}

final_results <- do.call(rbind, results_list) %>% arrange(Logrank_p)

write.csv(final_results, "KM_Survival_Results_All_Features.csv", row.names = FALSE)
cat("\n=== DONE ===\n")
cat("Significant features (p < 0.05):\n")
print(final_results[final_results$Significant, c("Feature", "Logrank_p", "HR_info")])
