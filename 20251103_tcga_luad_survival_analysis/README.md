# TCGA LUAD Survival Analysis

This analysis was performed using expression, clincal and survival data obtained from the TCGA Lung Adenocarcinoma (LUAD) cohort via [Xena Browser](https://xenabrowser.net/datapages/?cohort=TCGA%20Lung%20Adenocarcinoma%20(LUAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).

This remains work-in-progress: so far I have fitted regularised Cox proportional hazards models using CoxnetSurvivalAnalysis (sksurv) using PCA-transformed expression features and engineered clinical variables via cross-validation (scikit-learn).

Since the dataset is very small I have not used an independent holdout set for evaluation, but instead used the average valiation set C-index from cross-valiation as an estimate of model performance. My best model so far has a C-index of ~75%.

Each step in this analysis is saved in a separate notebook.

## Data summary

* ~500 samples

* ~20000 genes

* ~100 clinical features

Expression:

* gene-level

* RSEM-normalized counts (RNA-Seq by Expectation-Maximization)

* log2(x+1) transformed

Survival:

* event_col: 'OS'

* time_col: 'OS.time'

* censoring rate: ~64%

## Next steps

* Expression data EDA
  * Invetigte the outlier in the plot of PCA1 vs PCA2 and consider removing
  * Filter genes using a threshold for the proportion of the modal value 
  * Try WGCNA + eigen genes as features

* Try alternative models
  * SurvivalTree
  * GradientBoostingSurvivalAnalysis

* Re-train using a different seed to check CI estimate stability - 75% feels a little too high for this data set

* Look for an alternative LUAD data set to use for external validation 


