# brca_METABRIC_project


# README – Kaplan–Meier Survival Analysis  

1. **Loads** the processed METABRIC dataset.
2. **Automatically** performs Kaplan–Meier survival analysis for **every feature**.
3. Smart handling of variable types:  
   - Continuous variables → divided into quartiles (or median split if low spread)  
   - Known ordinal variables (e.g. histologic grade 1/2/3) → forced to categorical  
   - Categorical / binary variables → used as-is (empty strings and whitespace cleaned to NA)  
4. **Correct survival coding**: `overall_survival == 0` = death (event), `1` = alive (censored)  
5. For each feature it automatically:  
   - Builds the KM curve  
   - Performs log-rank test  
   - Calculates HR + 95% CI (single HR for binary variables; pairwise HRs vs reference for multi-level variables)  
   - Shows a publication-quality plot with p-value and HR directly on the figure  
   - all results saved as high-resolution PNG in the folder `KM_plots/`  
6. Produces one clean CSV file:  
   `KM_Survival_Results_{name_you_like}.csv`  
   Contains for every feature:  
   - Feature name  
   - Type (Numeric / Categorical)  
   - Groups used  
   - Number of patients (complete cases)  
   - Log-rank p-value  
   - Full HR information (single or pairwise)  
   - Significance flag

### Key files generated
- `KM_plots/` → all KM plots (PNG)  
- `data/` → containing filtered_data, selected features by different models, and unique features
- `KM_analysis_results.csv` → final results 