# Survival analysis of lung adenocarcinoma (LUAD) data from TCGA

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

So far I've fitted a CoxnetSurvivalAnalysis model using engineerd clinical features and PCA-transformed expression features. Since the dataset is very small I have not used an independent holdout set for evaluation, but instead used the average valiation set C-index from cross-valiation as an estimate of model performance. 

Next steps:

* try SurvivalTree, GradientBoostingSurvivalAnalysis

* look for an alternative TCGA LUAD data set to use for external validation 

* try graph-based extraction of co-expressed genes using kNN + Leiden/Louvain


