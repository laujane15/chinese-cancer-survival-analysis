# chinese-cancer-survival-analysis
Comprehensive survival analysis of 10,000 Chinese cancer patients using Kaplan-Meier curves, log-rank tests, and Cox proportional hazards modeling in R. Includes data preprocessing, categorical variable grouping, and visualizations.

# Key Findings:
Metastasis was the strongest predictor of mortality risk. Patients with metastatic disease had a 2.8-fold increased hazard of death (HR 3.8, 95% CI: 3.5-4.2, p<0.001) compared to the baseline of those without metastasis.
Tumor size demonstrated a dose-response relationship with survival. Compared to small tumors (<4.1 cm), patients with average-sized tumors (4.1-8.6 cm) had a 90% increased hazard (HR 1.9, 95% CI: 1.6-2.3, p<0.001), while those with large tumors (>8.6 cm) had a 150% increased hazard of death (HR 2.5, 95% CI: 2.1-3.1, p<0.001).
Ethnicity showed a modest but statistically significant association with survival. Han Chinese patients had a 30% increased hazard of death compared to ethnic minorities (HR 1.3, 95% CI: 1.1-1.5, p=0.002).
# Model Performance:
Concordance index: 0.71 (indicating good discriminatory ability)
Global p-value: <0.001 (highly significant model)
Number of events: 2,210 deaths out of 10,000 patients

# Clinical Implications:
These results confirm that metastasis status and tumor size are critical prognostic factors for survival in this Chinese cancer patient cohort. The findings support current staging systems that emphasize metastatic spread and primary tumor characteristics. The observed ethnic differences in survival warrant further investigation and may reflect disparities in healthcare access, genetic factors, or unmeasured confounders.
# Model limitations: 
Cancer stage which was a strong predictor from Kaplan-Meier estimates was excluded from the final model due to numerical instability, likely reflecting the impact of staging.
