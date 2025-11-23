# Survival analysis of LUAD data from TGCA

## Raw data

* ~500 samples

* ~20000 genes

* ~100 clinical features

## Raw expression features:

* gene-level

* RSEM-normalized counts (RNA-Seq by Expectation-Maximization)

* log2(x+1) transformed

## Raw survival features:

* event_col: 'OS'

* time_col: 'OS.time'

* censoring rate: ~64%

## Analysis notes

So far I've fitting a CoxnetSurvivalAnalysis model using engineerd clinical features and PCA transformed expression features. Since the dataset is so small I have not used a holdout set for evaluation but instead used the average of the validatation results from cross-valiation as an estimate of model performance. 

Next steps:

* Try SurvivalTree, GradientBoostingSurvivalAnalysis, ComponentwiseGradientBoostingSurvivalAnalysis

* look for an alternative TCGA LUAD data set to use for external validation 


