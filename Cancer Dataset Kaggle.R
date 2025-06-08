pacman::p_load(
  tidyverse,      # for data management and viz
  slider,         # for calculating moving averages
  tidyquant,       # for calculating moving averages within ggplot
  janitor,
  survminer,
  lubridate,
  survival
)
df <- rio::import("china_cancer_patients_synthetic.csv")
# Convert DiagnosisDate to Date format and add months
df <- df %>%
  mutate(
    DiagnosisDate = as.Date(DiagnosisDate),  # Convert to Date format
    FollowUpEndDate = DiagnosisDate %m+% months(FollowUpMonths)  # Add months using lubridate
  )

# Use this function for clean age ranges:
create_age_ranges_clean <- function(data, age_column, n_groups) {
  breaks <- quantile(data[[age_column]],
                    probs = seq(0, 1, 1/n_groups),
                    na.rm = TRUE)

  breaks_rounded <- round(breaks)

  labels <- character(n_groups)
  for(i in 1:n_groups) {
    if(i == 1) {
      labels[i] <- paste0("age_", breaks_rounded[i], "_", breaks_rounded[i+1])
    } else {
      labels[i] <- paste0("age_", breaks_rounded[i]+1, "_", breaks_rounded[i+1])
    }
  }

  age_categories <- cut(data[[age_column]],
                       breaks = breaks,
                       labels = labels,
                       include.lowest = TRUE)

  return(age_categories)
}

# Apply it:
df$age_4groups_ranges <- create_age_ranges_clean(df, "Age", 4)
mean_size <- mean(df$TumorSize, na.rm = TRUE)
sd_size <- sd(df$TumorSize, na.rm = TRUE)

# Using 1 standard deviation cuts (68% of data in middle)
# Calculate the actual cutoff points
small_cutoff <- mean_size - sd_size
large_cutoff <- mean_size + sd_size

# Alternative version showing actual size ranges
df$TumorSize_WithRanges <- cut(df$TumorSize,
                              breaks = c(-Inf, small_cutoff, large_cutoff, Inf),
                              labels = c(
                                paste0("Small (<", round(small_cutoff, 1), " cm)"),
                                paste0("Average (", round(small_cutoff, 1), "-", round(large_cutoff, 1), " cm)"),
                                paste0("Large (>", round(large_cutoff, 1), " cm)")
                              ))

# Chemotherapy Sessions - Clinical Categories
df$chemo_clinical_groups <- cut(df$ChemotherapySessions,
                               breaks = c(-1, 0, 3, 6, 12, Inf),
                               labels = c("No Chemo", "Short Course (1-3)",
                                         "Standard Course (4-6)", "Extended Course (7-12)",
                                         "Intensive (>12)"),
                               include.lowest = TRUE)

# Radiation Sessions - Clinical Categories
df$radiation_clinical_groups <- cut(df$RadiationSessions,
                                   breaks = c(-1, 0, 5, 15, 25, 35, Inf),
                                   labels = c("No Radiation", "Palliative (1-5)",
                                             "Short Course (6-15)", "Standard (16-25)",
                                             "Extended (26-35)", "Intensive (>35)"),
                                   include.lowest = TRUE)

grouping_rules <- list(
  # Tumor types - group by organ system
  TumorType = list(
    "Respiratory" = c("Lung", "Bronchial", "Pleural"),
    "Digestive" = c("Stomach", "Colon", "Liver", "Pancreatic"),
    "Reproductive" = c("Breast", "Ovarian", "Prostate", "Cervical"),
    "Hematologic" = c("Leukemia", "Lymphoma", "Myeloma"),
    "Other" = "OTHER"  # Catch-all
  ),

  # Provinces - group by region
  Province = list(
    "Eastern" = c("Beijing", "Shanghai", "Tianjin", "Jiangsu", "Zhejiang", "Fujian", "Shandong", "Guangdong", "Hainan", "Hebei", "Liaoning"),
    "Central" = c("Anhui", "Jiangxi", "Henan", "Hubei", "Hunan", "Shanxi", "Heilongjiang", "Jilin"),
    "Western" = c("Sichuan", "Yunnan", "Guizhou", "Shaanxi", "Gansu", "Qinghai", "Tibet", "Xinjiang", "Inner Mongolia", "Ningxia"),
    "Other" = "OTHER"
  ),


  # Cancer stages - simplify
  CancerStage = list(
    "Early" = c("I"),
    "Intermediate" = c("II", "III"),
    "Advanced" = c("IV")
  )
)

# Function to apply domain-specific grouping
apply_domain_grouping <- function(data, var_name, grouping_rules) {
  if(!var_name %in% names(grouping_rules)) {
    return(data[[var_name]])  # Return original if no rules
  }

  rules <- grouping_rules[[var_name]]
  new_var <- rep("Other", nrow(data))

  for(group_name in names(rules)) {
    if(group_name == "Other") next
    mask <- data[[var_name]] %in% rules[[group_name]]
    new_var[mask] <- group_name
  }

  # Handle NAs
  new_var[is.na(data[[var_name]])] <- "Missing"

  return(factor(new_var))
}

cat("\n\n=== DOMAIN-SPECIFIC GROUPING ===\n")
for(var in names(grouping_rules)) {
  new_var_name <- paste0(var, "_domain")
  df[[new_var_name]] <- apply_domain_grouping(df, var, grouping_rules)

  cat(paste0(var, " (domain grouped):\n"))
  print(table(df[[new_var_name]]))
  cat("\n")
}

df$ethnicity_binary <- ifelse(df$Ethnicity == "Han", "Han", "Ethnic Minority")
create_comorbidity_binary_quick <- function(data, comorbidity_col) {

  # Get the comorbidity column
  comorbid <- data[[comorbidity_col]]

  # Replace blank with "Missing"
  comorbid[comorbid == "" | is.na(comorbid) | trimws(comorbid) == ""] <- "Missing"

  # Create binary indicators for each condition
  data$Comorbid_Diabetes <- ifelse(
    grepl("Diabetes", comorbid, ignore.case = TRUE) & comorbid != "Missing", 1,
    ifelse(comorbid == "Missing", NA, 0)
  )

  data$Comorbid_Hypertension <- ifelse(
    grepl("Hypertension", comorbid, ignore.case = TRUE) & comorbid != "Missing", 1,
    ifelse(comorbid == "Missing", NA, 0)
  )

  data$Comorbid_Hepatitis_B <- ifelse(
    grepl("Hepatitis B", comorbid, ignore.case = TRUE) & comorbid != "Missing", 1,
    ifelse(comorbid == "Missing", NA, 0)
  )

  # Create summary variables
  data$Comorbid_Any <- ifelse(
    comorbid == "Missing", NA,
    ifelse(comorbid == "None", 0, 1)
  )

  # Count total comorbidities per patient
  data$Comorbid_Count <- ifelse(
    comorbid == "Missing", NA,
    ifelse(comorbid == "None", 0,
           data$Comorbid_Diabetes + data$Comorbid_Hypertension + data$Comorbid_Hepatitis_B)
  )

  # Clean the original comorbidity column
  data$Comorbidities_Cleaned <- comorbid

  return(data)
}

df <- create_comorbidity_binary_quick(df, "Comorbidities")

df$GeneticMutation[is.na(df$GeneticMutation)] <- "Missing"
df$genetic_binary <- ifelse(df$GeneticMutation == "Missing", "Missing",
                           ifelse(grepl("none|negative|no|normal",
                                       tolower(df$GeneticMutation)),
                                 "No Mutation", "Mutation Present"))
df$event <- ifelse(is.na(df$SurvivalStatus) | df$SurvivalStatus == 'Alive', 0, 1)

survobj <-Surv(time = df$FollowUpMonths,
               event = df$event)



# Overall survival
km_fit <- survfit(survobj ~ 1, data = df)
plot(km_fit, main = "Overall Survival Curve")

# =============================================================================
# LOG-RANK TEST COMPARISON FUNCTION
# =============================================================================

# Function to run log-rank tests across multiple variables
run_logrank_comparison <- function(data, survival_object, variables,
                                  alpha = 0.05, sort_by_pvalue = TRUE) {

  # Initialize results data frame
  results <- data.frame(
    Variable = character(),
    Groups = character(),
    N_Total = numeric(),
    N_Events = numeric(),
    Chi_Square = numeric(),
    DF = numeric(),
    P_Value = numeric(),
    Significant = character(),
    Effect_Size = character(),
    stringsAsFactors = FALSE
  )

  cat("=== LOG-RANK TESTS ===\n")
  cat("Testing variables:", paste(variables, collapse = ", "), "\n")
  cat("Significance level:", alpha, "\n\n")

  # Loop through each variable
  for(var in variables) {

    # Check if variable exists in data
    if(!var %in% names(data)) {
      cat("Warning: Variable", var, "not found in data. Skipping...\n")
      next
    }

    # Check if variable has sufficient groups
    var_table <- table(data[[var]], useNA = "ifany")
    if(length(var_table) < 2) {
      cat("Warning: Variable", var, "has less than 2 groups. Skipping...\n")
      next
    }

    # Remove missing values for this analysis
    complete_cases <- !is.na(data[[var]]) & !is.na(survival_object[,1]) & !is.na(survival_object[,2])

    if(sum(complete_cases) < 10) {
      cat("Warning: Variable", var, "has insufficient complete cases. Skipping...\n")
      next
    }

    # Run log-rank test
    tryCatch({
      formula_str <- paste("survival_object[complete_cases] ~", var)
      logrank_result <- survdiff(as.formula(formula_str), data = data[complete_cases, ])

      # Extract results
      chi_square <- logrank_result$chisq
      df <- length(logrank_result$n) - 1
      p_value <- 1 - pchisq(chi_square, df)

      # Determine significance
      significant <- ifelse(p_value < alpha, "Yes", "No")

      # Effect size interpretation (rough guide based on chi-square)
      effect_size <- ifelse(chi_square < 3.84, "Small",
                           ifelse(chi_square < 10.83, "Medium", "Large"))

      # Group information
      groups <- paste(names(var_table), collapse = ", ")
      n_total <- sum(logrank_result$n)
      n_events <- sum(logrank_result$obs)

      # Add to results
      results <- rbind(results, data.frame(
        Variable = var,
        Groups = groups,
        N_Total = n_total,
        N_Events = n_events,
        Chi_Square = round(chi_square, 3),
        DF = df,
        P_Value = round(p_value, 6),
        Significant = significant,
        Effect_Size = effect_size,
        stringsAsFactors = FALSE
      ))

      cat("✓ Completed:", var, "(p =", round(p_value, 4), ")\n")

    }, error = function(e) {
      cat("✗ Error with variable", var, ":", e$message, "\n")
    })
  }

  # Sort results by p-value if requested
  if(sort_by_pvalue && nrow(results) > 0) {
    results <- results[order(results$P_Value), ]
  }

  cat("\n=== LOG-RANK TEST RESULTS SUMMARY ===\n")
  cat("Total variables tested:", nrow(results), "\n")
  cat("Significant results (p <", alpha, "):", sum(results$Significant == "Yes"), "\n")

  return(results)
}


# =============================================================================
# ENHANCED VERSION WITH MEDIAN SURVIVAL TIMES
# =============================================================================

run_logrank_detailed <- function(data, survival_object, variables, alpha = 0.05) {

  # Initialize detailed results
  detailed_results <- data.frame(
    Variable = character(),
    Group = character(),
    N = numeric(),
    Events = numeric(),
    Median_Survival = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    stringsAsFactors = FALSE
  )

  # Initialize summary results
  summary_results <- data.frame(
    Variable = character(),
    N_Groups = numeric(),
    Chi_Square = numeric(),
    P_Value = numeric(),
    Significant = character(),
    stringsAsFactors = FALSE
  )

  for(var in variables) {

    if(!var %in% names(data)) next

    # Remove missing values
    complete_cases <- !is.na(data[[var]]) & !is.na(survival_object[,1]) & !is.na(survival_object[,2])

    if(sum(complete_cases) < 10) next

    tryCatch({
      # Fit Kaplan-Meier by groups
      formula_str <- paste("survival_object[complete_cases] ~", var)
      km_fit <- survfit(as.formula(formula_str), data = data[complete_cases, ])

      # Log-rank test
      logrank_result <- survdiff(as.formula(formula_str), data = data[complete_cases, ])

      # Extract median survival times and confidence intervals
      medians <- quantile(km_fit, probs = 0.5)

      # Get group-specific results
      groups <- levels(as.factor(data[[var]][complete_cases]))

      for(i in 1:length(groups)) {
        group_name <- groups[i]

        detailed_results <- rbind(detailed_results, data.frame(
          Variable = var,
          Group = group_name,
          N = logrank_result$n[i],
          Events = logrank_result$obs[i],
          Median_Survival = ifelse(is.na(medians$quantile[i]), Inf, round(medians$quantile[i], 1)),
          CI_Lower = ifelse(is.na(medians$lower[i]), NA, round(medians$lower[i], 1)),
          CI_Upper = ifelse(is.na(medians$upper[i]), NA, round(medians$upper[i], 1)),
          stringsAsFactors = FALSE
        ))
      }

      # Summary results
      chi_square <- logrank_result$chisq
      df <- length(logrank_result$n) - 1
      p_value <- 1 - pchisq(chi_square, df)

      summary_results <- rbind(summary_results, data.frame(
        Variable = var,
        N_Groups = length(groups),
        Chi_Square = round(chi_square, 3),
        P_Value = round(p_value, 6),
        Significant = ifelse(p_value < alpha, "Yes", "No"),
        stringsAsFactors = FALSE
      ))

    }, error = function(e) {
      cat("Error with", var, ":", e$message, "\n")
    })
  }

  return(list(
    summary = summary_results[order(summary_results$P_Value), ],
    detailed = detailed_results
  ))
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Define variables of interest for survival analysis
survival_variables <- c(
  "Gender",
  "Metastasis",
  "TreatmentType",
  "age_4groups_ranges",
  "TumorSize_WithRanges",
  "chemo_clinical_groups",
  "radiation_clinical_groups",
  "TumorType_domain",
  "Province_domain",
  "CancerStage_domain",
  "ethnicity_binary",
  "Comorbid_Diabetes",
  "Comorbid_Hypertension",
  "Comorbid_Hepatitis_B",
  "Comorbid_Any",
  "Comorbid_Count",
  "Comorbidities_Cleaned",
  "genetic_binary"
)

# Filter to only include variables that exist in your dataset
existing_variables <- survival_variables[survival_variables %in% names(df)]

cat("=== AVAILABLE VARIABLES FOR TESTING ===\n")
cat("Variables found in dataset:\n")
for(var in existing_variables) {
  n_groups <- length(unique(df[[var]]))
  cat(sprintf("  %-25s: %d groups\n", var, n_groups))
}

# Run basic log-rank comparison
cat("\n" , rep("=", 60), "\n")
logrank_results <- run_logrank_comparison(df, survobj, existing_variables)

# Display results table
cat("\n=== LOG-RANK TEST COMPARISON TABLE ===\n")
print(logrank_results)

# Run detailed analysis
cat("\n" , rep("=", 60), "\n")
detailed_results <- run_logrank_detailed(df, survobj, existing_variables)

cat("\n=== SUMMARY TABLE ===\n")
print(detailed_results$summary)

cat("\n=== DETAILED RESULTS (First 20 rows) ===\n")
print(head(detailed_results$detailed, 20))

# =============================================================================
# EXPORT RESULTS
# =============================================================================

# Save results to CSV files
write.csv(logrank_results, "logrank_comparison.csv", row.names = FALSE)
write.csv(detailed_results$summary, "logrank_summary.csv", row.names = FALSE)
write.csv(detailed_results$detailed, "logrank_detailed.csv", row.names = FALSE)

cat("\n=== RESULTS EXPORTED ===\n")
cat("✓ logrank_comparison.csv - Basic comparison table\n")
cat("✓ logrank_summary.csv - Summary of all tests\n")
cat("✓ logrank_detailed.csv - Detailed group-specific results\n")

# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================

cat("\n\n=== INTERPRETATION GUIDE ===\n")
cat("Chi-Square Values:\n")
cat("  < 3.84   : Small effect (not significant at p=0.05)\n")
cat("  3.84-10.83: Medium effect\n")
cat("  > 10.83  : Large effect\n\n")

cat("P-Value Interpretation:\n")
cat("  < 0.001  : Highly significant (***)\n")
cat("  < 0.01   : Very significant (**)\n")
cat("  < 0.05   : Significant (*)\n")
cat("  ≥ 0.05   : Not significant\n\n")

cat("Next Steps:\n")
cat("1. Focus on variables with p < 0.05\n")
cat("2. Include significant variables in Cox model\n")
cat("3. Check for interactions between significant variables\n")
cat("4. Validate results with larger sample if needed\n")

# Show most significant results
if(nrow(logrank_results) > 0) {
  significant_vars <- logrank_results[logrank_results$Significant == "Yes", ]
  if(nrow(significant_vars) > 0) {
    cat("\n MOST SIGNIFICANT VARIABLES:\n")
    print(significant_vars[1:min(5, nrow(significant_vars)), c("Variable", "P_Value", "Chi_Square")])
  }
}



# =============================================================================
# 1. ENHANCED METASTASIS PLOT
# =============================================================================

if("Metastasis" %in% names(df)) {
  cat("Creating enhanced Metastasis plot...\n")

  fit_metastasis <- survfit(survobj ~ Metastasis, data = df)

  enhanced_metastasis <- ggsurvplot(
    fit_metastasis,
    data = df,
    pval = TRUE,
    conf.int = T,
    pval.size = 5,
    risk.table = TRUE,
    risk.table.height = 0.3,
    xlab = "Follow-up time (months)",
    ylab = "Survival Probability (%)",
    title = "Kaplan-Meier Survival Curves by Metastasis Status",
    surv.scale = "percent",
    break.time.by = 12,
    # BLACK AND WHITE STYLING
    palette = c("black", "grey50"),  # Black and grey
    linetype = c("solid", "dashed"), # Different line types
    size = 1.2,                      # Thicker lines
    legend.title = "Metastasis",
    font.main = 14,
    font.x = 12,
    font.y = 12,
    surv.median.line = "hv",
    ggtheme = theme_classic() +      # Clean white background
      theme(
        panel.grid = element_blank(),
        legend.position = "top",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = "black")
      )
  )

  print(enhanced_metastasis)

  ggsave("survival_metastasis.png", plot = enhanced_metastasis$plot,
         width = 12, height = 8, dpi = 300)

  cat("✓ Enhanced Metastasis plot created\n")
} else {
  cat("⚠ Metastasis variable not found\n")
}

# =============================================================================
# 2. ENHANCED TUMOR SIZE PLOT
# =============================================================================

fit_tumorsize_ranges <- survfit(survobj ~ TumorSize_WithRanges, data = df)

# Plot with size ranges
ranges_tumorsize_plot <- ggsurvplot(
  fit_tumorsize_ranges,
  data = df,
  pval = TRUE,
  pval.size = 5,
  conf.int = T,
  risk.table = TRUE,
  risk.table.height = 0.3,
  xlab = "Follow-up time (months)",
  ylab = "Survival Probability (%)",
  title = "Kaplan-Meier Survival Curves by Tumor Size",
  surv.scale = "percent",
  break.time.by = 12,
  # BLACK AND WHITE STYLING
  palette = c("black", "grey40", "grey70"),     # Different greys
  linetype = c("solid", "dashed", "dotted"),    # Different line types
  legend.title = "Tumor Size",
  legend.labs = c("Small (<4.1 cm)", "Average (4.1-8.6 cm)", "Large (>8.6 cm)"),
  size = 1.2,                      # Thicker lines
  font.main = 14,
  font.x = 12,
  font.y = 12,
  font.legend = 10,  # Smaller font for longer labels
  surv.median.line = "hv",
  ggtheme = theme_classic() +      # Clean white background
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black")
    )
)

print(ranges_tumorsize_plot)

# Save ranges plot
ggsave("survival_tumorsize.png", plot = ranges_tumorsize_plot$plot,
       width = 12, height = 8, dpi = 300)

# =============================================================================
# 3. ENHANCED CANCER STAGE PLOT
# =============================================================================

# Use CancerStage_domain and fix the factor ordering
df$Stage <- df$CancerStage_domain

# CRITICAL: Set the correct factor order (Early → Intermediate → Advanced)
df$Stage <- factor(df$Stage, levels = c("Early", "Intermediate", "Advanced"))

# Verify the ordering is correct
cat("Stage distribution (should be in correct order):\n")
print(table(df$Stage, useNA = "ifany"))

# =============================================================================
# FIX 2: CREATE PLOT WITH CUSTOM LEGEND LABELS
# =============================================================================

if(sum(!is.na(df$Stage)) > 0) {

  # Create survival fit
  fit_stage_final <- survfit(survobj ~ Stage, data = df)

  # Create plot with custom legend labels
  perfect_stage_plot <- ggsurvplot(
    fit_stage_final,
    data = df,
    pval = TRUE,
    pval.size = 5,
    conf.int = T,
    risk.table = TRUE,
    risk.table.height = 0.3,
    xlab = "Follow-up time (months)",
    ylab = "Survival Probability (%)",
    title = "Kaplan-Meier Survival Curves by Cancer Stage",
    surv.scale = "percent",
    break.time.by = 12,
    # BLACK AND WHITE STYLING
    palette = c("black", "grey40", "grey70"),     # Different greys
    linetype = c("solid", "dashed", "dotted"),    # Different line types
    legend.title = "Cancer Stage",
    # CUSTOM LEGEND LABELS with stage numbers
    size = 1.2,                      # Thicker lines
    legend.labs = c("Early (I)", "Intermediate (II-III)", "Advanced (IV)"),
    font.main = 14,
    font.x = 12,
    font.y = 12,
    surv.median.line = "hv",
    ggtheme = theme_classic() +      # Clean white background
      theme(
        panel.grid = element_blank(),
        legend.position = "top",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = "black")
      )
  )

  print(perfect_stage_plot)

  # Save the corrected plot
  ggsave("survival_cancer_stage.png", plot = perfect_stage_plot$plot,
         width = 12, height = 8, dpi = 300)

} else {
  cat("No valid Stage data found\n")
}
# =============================================================================
# 4. ENHANCED ETHNICITY PLOT
# =============================================================================

if("ethnicity_binary" %in% names(df)) {
  cat("Creating B&W ethnicity plot...\n")

  fit_ethnicity_bw <- survfit(survobj ~ ethnicity_binary, data = df)

  # Black and white version
  bw_ethnicity_clean <- ggsurvplot(
    fit_ethnicity_bw,
    data = df,
    pval = TRUE,
    pval.size = 5,
    conf.int = T,
    risk.table = TRUE,
    risk.table.height = 0.3,
    xlab = "Follow-up time (months)",
    ylab = "Survival Probability (%)",
    title = "Kaplan-Meier Survival Curves by Ethnicity",
    surv.scale = "percent",
    break.time.by = 12,
    # BLACK AND WHITE STYLING
    palette = c("black", "grey50"),
    linetype = c("solid", "dashed"),
    size = 1.2,
    legend.title = "Ethnicity",
    # Clean legend labels
    legend.labs = c("Ethnic Minority", "Han Chinese"),
    font.main = 14,
    font.x = 12,
    font.y = 12,
    surv.median.line = "hv",
    ggtheme = theme_classic() +
      theme(
        panel.grid = element_blank(),
        legend.position = "top",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = "black")
      )
  )

  print(bw_ethnicity_clean)

  # Save B&W version
  ggsave("survival_ethnicity.png", plot = bw_ethnicity_clean$plot,
         width = 12, height = 8, dpi = 300)

}

#fitting the cox model
df_cox <-  survival::coxph(
              Surv(FollowUpMonths, event) ~ Metastasis + TumorSize_WithRanges + ethnicity_binary,
              data = df
              )
# Stage dropped, Hazard Ratios problematic.
#StageIntermediate: HR = 76,260,110 (!)
#StageAdvanced: HR = 109,005,800 (!)
#The extreme HR values suggest perfect or near-perfect separation - meaning perhaps:

#All or most Early stage patients survived
#All or most Advanced stage patients died
#This creates numerical instability in the Cox model.

summary(df_cox)

test_df_cox <- survival::cox.zph(df_cox)
survminer::ggcoxzph(test_df_cox)
ggforest(df_cox, data = df)
