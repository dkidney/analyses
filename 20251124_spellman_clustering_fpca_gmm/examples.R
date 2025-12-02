library(curveclust)

spellman |> head()
spellman |> plot_curves(x='time', y='expression', subject='gene', colour='peak', facet='peak')



# Required packages
# library(dplyr)
library(conflicted)
library(tidyverse)
library(fda)        # for spline smoothing
library(fdapace)    # for FPCA / PACE
library(mclust)     # model-based clustering + ARI
# library(cluster)    # silhouette
# library(ggplot2)    # plotting
# library(purrr)

conflicted::conflicts_prefer(
	dplyr::select,
	dplyr::filter
)

# ---------------------------------------------------------------------
# Assume df is your input long dataframe with columns: gene, time, expression
# ---------------------------------------------------------------------

# Optional: if you have true labels for genes (e.g., cell-cycle phases)
# supply as a data frame 'labels_df' with columns gene and true_label.
# If not available, ARI step will be skipped.
# labels_df <- data.frame(gene = c(...), true_label = c(...))

# check data: column classes, row order
df = spellman |>
	mutate(
		time = as.numeric(time),
		gene = as.character(gene),
		expression = as.numeric(expression),
		peak = as.character(peak),
	) |>
	arrange(gene, time)

df

df$gene |> unique() |> length() # 791

# ---------------------------------------------------------------------
# 1) Prepare a gene-by-(time,expr) list for both pipelines
# ---------------------------------------------------------------------
# make a grouped data object
df_grouped <- df %>%
	group_by(gene, peak) %>%
	summarise(
		Lt = list(time),
		Ly = list(expression),
		.groups = "drop"
	)
df_grouped

genes <- grouped_data$gene

# ---------------------------------------------------------------------
# 2) Spline-based smoothing per gene -> extract basis coefficients
#    (Using fda::smooth.basis with a B-spline basis)
# ---------------------------------------------------------------------

# TUNE these parameters:
nbasis <- 8        # number of spline basis functions (try 6-12)
penalty_order <- 2 # second derivative penalty
lambda  <- 1e-2    # smoothing parameter (tune or choose via CV/GCV)

# create a basis object
basis <- fda::create.bspline.basis(
	rangeval = range(df$time, na.rm = TRUE),
	nbasis = 8,
	norder = 4 # basis degree + 1 = order
)
basis$type
basis$params # internal knots?
basis$rangeval # x-axis range
basis$nbasis # n basis functions = degrees of freedom
basis$names #
summary(basis)
plot(basis)

# create an fdPar object
fdParobj <- fda::fdPar(
	fdobj = basis,
	Lfdobj = penalty_order,
	lambda = lambda
)
fdParobj

# Function to smooth a single gene (time vector t, obs y)
smooth_gene <- function(t, y, fdParobj) {
	# check NAs and enough points
	if(sum(!is.na(y)) < 4) {
		# if too few points, return NA coefficients
		return(rep(NA, fdParobj$basis$nbasis))
	}
	res <- try(fda::smooth.basis(argvals = t, y = y, fdParobj = fdParobj), silent = TRUE)
	if(inherits(res, "try-error")) return(rep(NA, fdParobj$basis$nbasis))
	as.numeric(res$fd$coefs) # vector of length nbasis
}

# build coefficient matrix: genes x nbasis
coef_mat <- purrr::map2(df_grouped$Lt, df_grouped$Ly, ~ smooth_gene(.x, .y, fdParobj)) %>%
	do.call(rbind, .)
rownames(coef_mat) <- genes
colnames(coef_mat) <- paste0("b", 1:ncol(coef_mat))
head(coef_mat)

# Optionally: remove genes with NA coefficients (too sparse)
ok_genes <- which(rowSums(is.na(coef_mat)) == 0)
coef_mat <- coef_mat[ok_genes, , drop = FALSE]
genes_coef <- rownames(coef_mat)

# Standardize features (important before clustering)
coef_mat_scaled <- scale(coef_mat)

# ---------------------------------------------------------------------
# 3) FPCA (PACE) using fdapace
#    Supply Lt and Ly lists (only those genes we can use)
# ---------------------------------------------------------------------
# Build Lt, Ly lists for fdapace (use same subset of genes as for coef_mat, or use all)
Lt <- df_grouped$Lt
Ly <- df_grouped$Ly
names(Lt) <- names(Ly) <- df_grouped$gene

# If you want to restrict to genes with good smoothing, use genes_coef
Lt_sub <- Lt[genes_coef]
Ly_sub <- Ly[genes_coef]

# Run FPCA (tune options as needed)
# dataType = "Sparse" is appropriate for sparse / irregular sampling
fpca_res <- FPCA(Ly = Ly_sub, Lt = Lt_sub,
				 optns = list(dataType = "Sparse",
				 			 methodSelectK = "FVE", # select K by fraction of variance explained
				 			 FVEthreshold = 0.99))   # aim to explain e.g. 99% variance

scores <- fpca_res$xiEst         # genes x K (rows correspond to names(Lt_sub))
rownames(scores) <- names(Lt_sub)
head(scores)

# standardize FPC scores before clustering
scores_scaled <- scale(scores)

# ---------------------------------------------------------------------
# 4) Clustering both feature sets with mclust
#    We'll let Mclust pick best model and number of clusters by BIC,
#    but you can force G (n clusters) if desired.
# ---------------------------------------------------------------------

# 4A: cluster on spline coefficients
mc_coef <- Mclust(coef_mat_scaled)   # returns best G by BIC by default
summary(mc_coef)
clusters_coef <- mc_coef$classification
# make sure cluster vector aligns to genes_coef
names(clusters_coef) <- genes_coef

# 4B: cluster on FPCA scores
mc_fpca <- Mclust(scores_scaled)
summary(mc_fpca)
clusters_fpca <- mc_fpca$classification
names(clusters_fpca) <- rownames(scores)

# ---------------------------------------------------------------------
# 5) Compare clusters
#    - BIC: use mc_coef$bic and mc_fpca$bic (mclust stores BIC per model)
#    - Silhouette: compute average silhouette width
#    - ARI: if you have true labels, compute Adjusted Rand Index
# ---------------------------------------------------------------------

# BIC comparison: higher is better
bic_coef_best <- mc_coef$bic
bic_fpca_best <- mc_fpca$bic
cat("Mclust BIC (coef):", max(bic_coef_best, na.rm=TRUE), "\n")
cat("Mclust BIC (fpca):", max(bic_fpca_best, na.rm=TRUE), "\n")

# # Silhouette: need a distance matrix and cluster labels
# # For spline-coefs
# dist_coef <- dist(coef_mat_scaled)
# sil_coef <- silhouette(clusters_coef, dist_coef)
# avg_sil_coef <- mean(sil_coef[, 'sil_width'])
# cat("Average silhouette (spline-coefs):", round(avg_sil_coef, 3), "\n")

# # For fpca scores
# dist_fpca <- dist(scores_scaled)
# sil_fpca <- silhouette(clusters_fpca, dist_fpca)
# avg_sil_fpca <- mean(sil_fpca[, 'sil_width'])
# cat("Average silhouette (FPCA):", round(avg_sil_fpca, 3), "\n")

# ARI if true labels available
if(exists("labels_df")) {
	labvec <- labels_df$true_label[match(names(clusters_coef), labels_df$gene)]
	# remove NAs if any
	ok <- !is.na(labvec)
	if(sum(ok) > 0) {
		ari_coef <- adjustedRandIndex(clusters_coef[ok], labvec[ok])
		cat("ARI (spline-coefs):", ari_coef, "\n")
	}
	labvec2 <- labels_df$true_label[match(names(clusters_fpca), labels_df$gene)]
	ok2 <- !is.na(labvec2)
	if(sum(ok2) > 0) {
		ari_fpca <- adjustedRandIndex(clusters_fpca[ok2], labvec2[ok2])
		cat("ARI (FPCA):", ari_fpca, "\n")
	}
}

# ---------------------------------------------------------------------
# 6) Visualize cluster mean trajectories (reconstruct smoothed curves)
#    - For spline coefs: evaluate fd object on a dense grid
#    - For FPCA: use fpca_res to reconstruct mean + score*phi
# ---------------------------------------------------------------------

# Dense time grid
tt <- seq(0, 120, length.out = 101)

# 6A: spline-coefs mean curves per cluster
# Recreate fd objects for each gene from coef_mat and basis
fd_coefs <- t(coef_mat) # nbasis x nGenes
fdobj_genes <- fd(fd_coefs, basis) # columns = genes
# evaluate
mat_eval <- eval.fd(tt, fdobj_genes) # matrix 101 x nGenes

plot_df_coef <- data.frame(
	time = rep(tt, times = ncol(mat_eval)),
	expr = as.numeric(mat_eval),
	gene = rep(colnames(mat_eval), each = nrow(mat_eval))
)
plot_df_coef$cluster <- clusters_coef[match(plot_df_coef$gene, names(clusters_coef))]

# compute cluster means
cluster_mean_coef <- plot_df_coef %>%
	group_by(cluster, time) %>%
	summarize(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")

ggplot(cluster_mean_coef, aes(x = time, y = mean_expr, color = factor(cluster))) +
	geom_line(size = 1) +
	labs(title = "Cluster mean curves (spline coefficients)", color = "cluster")

# 6B: FPCA-based reconstructions per cluster
# fpca_res has phi (matrix of eigenfunctions evaluated at fpca_res$workGrid)
# xiEst are scores
# reconstruct: mu + phi %*% xi
workGrid <- fpca_res$workGrid   # grid where phi is evaluated
phi <- fpca_res$phi             # matrix length(workGrid) x K
mu  <- fpca_res$mu              # mean curve length(workGrid)
recon <- t( t(mu) + phi %*% t(scores) ) # rows = genes, columns = workGrid
# recon currently genes x grid; transpose to grid x gene
recon_t <- t(recon) # grid x gene
# build df for cluster mean plotting
recon_df <- as.data.frame(recon_t)
colnames(recon_df) <- colnames(recon_t) # may be numeric
recon_df$time <- workGrid
recon_long <- pivot_longer(recon_df, cols = -time, names_to = "gene", values_to = "expr")
recon_long$cluster <- clusters_fpca[match(recon_long$gene, names(clusters_fpca))]

cluster_mean_fpca <- recon_long %>%
	group_by(cluster, time) %>%
	summarize(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")

ggplot(cluster_mean_fpca, aes(x = time, y = mean_expr, color = factor(cluster))) +
	geom_line(size = 1) +
	labs(title = "Cluster mean curves (FPCA)", color = "cluster")

# ---------------------------------------------------------------------
# 7) Additional diagnostics and extensions
# ---------------------------------------------------------------------
# - Try different nbasis or lambda for spline approach.
# - Try different FPCA options, e.g. methodSelectK, FVEthreshold.
# - Try clustering in reduced space: for spline coefficients, you could do PCA on coefficients first,
#   then cluster on first few PCs (to reduce collinearity).
# - Try funHDDC / funFEM or clustering directly on FPCA scores with model constraints.
# - Evaluate cluster stability via bootstrapping (fpc::clusterboot or manual subsampling).






results = spellman |> cluster_curves(x='time', y='expression', subject='gene', dfs=10)
results |> head()
results |> plot_curves(x='time', y='y_hat', subject='gene', colour='peak', facet='peak')
results |> plot_curves(x='time', y='y_hat', subject='gene', colour='peak', facet='cluster')


fit_bspline_regression(spellman_subset, x='time', y='expression', subject='gene', df=5)



# fit5 = fit_bspline_regression(spellman, x='time', y='expression', subject='gene', df=5)
# fit6 = fit_bspline_regression(spellman, x='time', y='expression', subject='gene', df=6)
# aic(fit5)
# aic(fit6)

# spellman_clusters |> with(plot_curves(x='expression', y='time', subject='gene', group='peak', facet='cluster'))


basis = bspline_basis(spellman, x='time', df=5)
basis$degree
basis$intercept
basis$df
basis$knots
basis$matrix[1:5, ]
plot_basis(basis)
