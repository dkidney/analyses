# TCGA LUAD Survival Analysis

This analysis was performed using expression, clincal and survival data obtained from the TCGA Lung Adenocarcinoma (LUAD) cohort via [Xena Browser](https://xenabrowser.net/datapages/?cohort=TCGA%20Lung%20Adenocarcinoma%20(LUAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).

So far I have fitted regularised Cox proportional hazards models using CoxnetSurvivalAnalysis (sksurv) using feature-engineered variables via cross-validation (scikit-learn).

Since the dataset is relatively small (see below) I have not used an independent holdout set for evaluation, but instead used the average valiation set C-index from cross-valiation as an estimate of model performance. 

I have tried two dimension reduction techniques so far, both of which gave C-index scores of ~75:

1. Using PCA transformed expression features within cross-valiation - C-index = 77% (bootstrap: mean 0.758, 95% CI (0.704, 0.807))

1. Using a Weighted Gene Co-expression Network analysis approach where eigengenes were derived prior to cross-validation - C-index 75% (bootstrap: mean 0.751, 95% CI (0.701, 0.801))

The WGCNA method was ~50x faster to fit, had fewer non-zero coefficients and had similar bootstrap CI results, so this would be my preferred model.

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

* Try alternative models - e.g. SurvivalTree, GradientBoostingSurvivalAnalysis

* Look for an alternative LUAD data set to use for external validation 


