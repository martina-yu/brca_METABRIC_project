library(tidyr)
library(survival)
library(survminer)
library(dplyr)

km_one_feature <- function(feature_name, data = filtered_df) {
  
  cat("\n=== Analyzing:", feature_name, "===\n")
  x <- data[[feature_name]]
  
  # Check if feature exists and has variation
  if (is.null(x) || all(is.na(x)) || length(unique(na.omit(x))) < 2) {
    cat("   Skipped: no variation or all NA\n")
    return(invisible(NULL))
  }
  
  # ——— Initialize variables ———
  plot_data <- data
  group_plot <- NULL
  plot_title <- feature_name
  
  # ——— Grouping logic ———
  if (is.numeric(x) && !is.factor(x)) {
    # Numeric → quartiles
    quartiles <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    
    if (length(unique(quartiles)) < 4) {
      cat("   Low spread → median split\n")
      group <- ifelse(x <= median(x, na.rm = TRUE), "Low", "High")
    } else {
      group <- cut(x, breaks = quartiles, 
                   labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)"), 
                   include.lowest = TRUE)
    }
    
    plot_title <- paste0(feature_name, " (Quartiles)")
    
    # Remove NA from numeric grouping
    valid <- !is.na(group)
    plot_data <- data[valid, ]
    group_plot <- droplevels(as.factor(group[valid]))
    
  } else {
    # Categorical / Binary
    if (is.character(x)) {
      x <- trimws(x)                  # remove leading/trailing whitespace
      x[x == ""] <- NA                # convert empty string "" → real NA
    }
    group <- as.factor(x)
    
    cat("   After cleaning empty strings:\n")
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
    
    plot_title <- paste0(feature_name, " (n=", nrow(plot_data), " complete)")
  }
  
  # Final safety check
  if (nlevels(group_plot) < 2) {
    cat("   Less than 2 groups → skipping\n")
    return(invisible(NULL))
  }
  
  # Show group distribution
  cat("   Group distribution:\n")
  print(table(group_plot))
  
  
  ####### --------- Survival Analysis Model ------------- #######
  Surv_clean <- Surv(time = plot_data$overall_survival_months,
                     event = plot_data$overall_survival == 0)
  
  fit <- survfit(Surv_clean ~ group_plot)
  
  # Log-rank test
  p_logrank <- try({
    test <- survdiff(Surv_clean ~ group_plot)
    1 - pchisq(test$chisq, df = nlevels(group_plot) - 1)
  }, silent = TRUE)
  
  if (inherits(p_logrank, "try-error")) {
    p_logrank <- NA
  }
  
  
  ###### ---------- HR Calculation: Binary or Multi-group ----------- #######
  hr_text <- "Multi-group (HR not shown)"
  
  if (nlevels(group_plot) == 2) {
    # Binary: single HR
    cox_model <- try(coxph(Surv_clean ~ group_plot), silent = TRUE)
    
    if (!inherits(cox_model, "try-error") && nrow(summary(cox_model)$coefficients) > 0) {
      beta   <- coef(cox_model)
      hr     <- round(exp(beta), 3)
      ci     <- round(exp(confint(cox_model)), 3)
      p_val  <- summary(cox_model)$coefficients[5]
      
      ref_level  <- levels(group_plot)[1]
      risk_level <- levels(group_plot)[2]
      
      hr_text <- paste0("HR (", risk_level, " vs ", ref_level, ") = ", hr,
                        "\n95% CI: ", ci[1], "–", ci[2], 
                        ", p = ", ifelse(p_val < 0.001, "<0.001", round(p_val, 4)))
      
      cat("\n", hr_text, "\n")
    } else {
      hr_text <- "Cox model failed"
    }
    
  } else if (nlevels(group_plot) > 2) {
    # Multi-group: pairwise HRs vs first level (reference)
    cox_multi <- try(coxph(Surv_clean ~ group_plot), silent = TRUE)
    
    if (!inherits(cox_multi, "try-error")) {
      summ <- summary(cox_multi)
      ref_level <- levels(group_plot)[1]
      
      # Pre-calculate confidence intervals to avoid environment issues
      all_ci <- round(exp(confint(cox_multi)), 3)
      
      cat("\n=== Pairwise Hazard Ratios (vs", ref_level, ") ===\n")
      
      hr_lines <- sapply(2:nlevels(group_plot), function(i) {
        hr    <- round(exp(coef(cox_multi))[i-1], 3)
        ci    <- all_ci[i-1, ]
        p     <- summ$coefficients[i-1, 5]
        p_txt <- ifelse(p < 0.001, "<0.001", round(p, 4))
        
        line <- paste0(levels(group_plot)[i], " vs ", ref_level, 
                       ": HR = ", hr, " (95% CI: ", ci[1], "–", ci[2], 
                       ", p = ", p_txt, ")")
        cat(line, "\n")
        return(line)
      })
      
      hr_text <- paste("Pairwise HRs (vs", ref_level, "):", 
                       paste(hr_lines, collapse = "\n"), sep = "\n")
      
    } else {
      hr_text <- "Cox model failed (multi-group)"
    }
  }
  
  
  ###### ---------- Plot ----------- #######
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
  
  # Save plot
  safe_name <- gsub("[^A-Za-z0-9_]", "_", feature_name)
  png_file <- paste0("KM_plots/KM_", safe_name, "_p", signif(p_logrank, 3), ".png")
  ggsave(png_file, plot = p$plot, width = 9.5, height = 7.5, dpi = 300, bg = "white")
  cat("   Saved: ", png_file, "\n")
  
  # Print summary
  cat("\nLog-rank p =", ifelse(is.na(p_logrank), "NA", signif(p_logrank, 3)),
      ifelse(!is.na(p_logrank) && p_logrank < 0.05, " → SIGNIFICANT!\n", "\n"))
  
  
  ###### ---------- Return dataframe ----------- #######
  return(data.frame(
    Feature = feature_name,
    Type = ifelse(is.numeric(data[[feature_name]]), "Numeric", "Categorical"),
    Groups = paste(levels(group_plot), collapse = " | "),
    N = nrow(plot_data),
    Logrank_p = ifelse(is.na(p_logrank), NA, p_logrank),
    HR_info = hr_text,
    Significant = !is.na(p_logrank) && p_logrank < 0.05,
    stringsAsFactors = FALSE
  ))
}

