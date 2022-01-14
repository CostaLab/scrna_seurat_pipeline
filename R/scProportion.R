### adjusted from https://github.com/rpolicastro/scProportionTest.git
#' Single-Cell Utilities
#'
#' @slot meta_data Seurat meta-data
#' @slot results Results for various analysis
#'
#' @rdname sc_utils-class
#' @export

setClass(
	"sc_utils",
	representation(
		meta_data = "data.table",
		results = "list"
	),
	prototype(
		meta_data = data.frame(),
		results = list()
	)
)

#' Single-Cell Utilities Constructor
#'
#' @import methods
#' @import Seurat
#' @import data.table
#'
#' @param seurat_object Seurat object with meta-data
#'
#' @rdname sc_utils-class
#' @export

sc_utils <- function(seurat_object) {

	## Pull meta-data out of Seurat object.
	metadata <- as.data.table(
		seurat_object[[]],
		keep.rownames = "cell_id"
	)

	## Create new sc_utils object.
	sc_utils_obj <- new(
		"sc_utils",
		meta_data = metadata
	)

	return(sc_utils_obj)
}

#' Permutation Test For Proportions
#'
#' @param sc_utils_obj sc_utils object
#' @param cluster_identity Column that has cluster names
#' @param sample_1 First sample to compare (ie. control)
#' @param sample_2 Sample to compare to first sample (ie. treatment)
#' @param sample_identity Column that has sample names
#' @param n_permutations Number of permutations
#'
#' @rdname permutation_test-function
#'
#' @export

permutation_test <- function(
	sc_utils_obj,
	cluster_identity = NA,
	sample_1 = NA,
	sample_2 = NA,
	sample_identity = "orig.ident",
	n_permutations = 1000
) {

	## Prepare data.
	meta_data <- copy(sc_utils_obj@meta_data)

	meta_data <- meta_data[
		get(sample_identity) %in% c(sample_1, sample_2),
		c(..sample_identity, ..cluster_identity)
	]

	setnames(
		meta_data,
		old = c(sample_identity, cluster_identity),
		new = c("samples", "clusters")
	)

	meta_data[, clusters := as.character(clusters)]
	cluster_cases <- unique(meta_data[["clusters"]])

	## Get observed differences in fraction.
	obs_diff <- meta_data[, .(count = .N), by = .(samples, clusters)]
	obs_diff <- obs_diff[
		CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
		on = .(samples, clusters)
	][
		is.na(count), count := 0
	][]
	obs_diff[, fraction := count / sum(count), by = samples]
	obs_diff <- dcast(obs_diff, clusters ~ samples, value.var = "fraction")
	obs_diff[, obs_log2FD := log2(get(sample_2)) - log2(get(sample_1))]

	## Permutation test.
	perm_results <- matrix(NA, nrow(obs_diff), n_permutations)
	rownames(perm_results) <- sort(cluster_cases)

	for (i in seq_len(n_permutations)) {
		permuted <- copy(meta_data)
		permuted[["samples"]] <- sample(permuted[["samples"]])
		permuted <- permuted[, .(count = .N), by = .(samples, clusters)]
		permuted <- permuted[
			CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
			on = .(samples, clusters)
		][
			is.na(count), count := 0
		][]
		permuted[, fraction := count / sum(count), by = samples]
		permuted <- dcast(permuted, clusters ~ samples, value.var = "fraction")
		permuted[, perm_log2FD := log2(get(sample_2)) - log2(get(sample_1))]

		perm_results[, i] <- permuted[["perm_log2FD"]]
	}

	increased <- rowSums(apply(perm_results, 2, function(x) obs_diff[["obs_log2FD"]] <= x))
	increased <- (increased + 1) / (n_permutations + 1)

	decreased <- rowSums(apply(perm_results, 2, function(x) obs_diff[["obs_log2FD"]] >= x))
	decreased <- (decreased + 1) / (n_permutations + 1)

	obs_diff[, pval := ifelse(obs_log2FD > 0, increased[.I], decreased[.I])]
	obs_diff[, FDR := p.adjust(pval, "fdr")]

	## Boostrap log2FD CI.
	boot_results <- matrix(NA, nrow(obs_diff), n_permutations)
	rownames(boot_results) <- sort(cluster_cases)

	for (i in seq_len(n_permutations)) {
		booted <- copy(meta_data)
		booted[, clusters := sample(clusters, replace = TRUE), by = samples]
		booted <- booted[, .(count = .N), by = .(samples, clusters)]
		booted <- booted[
			CJ(samples = samples, clusters = cluster_cases, unique = TRUE),
			on = .(samples, clusters)
		][
			is.na(count), count := 0
		][]
		booted[, fraction := count / sum(count), by = samples]
		booted <- dcast(booted, clusters ~ samples, value.var = "fraction")
		booted[, boot_log2FD := log2(get(sample_2)) - log2(get(sample_1))]

		boot_results[, i] <- booted[["boot_log2FD"]]
	}

	boot_results[!is.finite(boot_results)] <- NA
	boot_mean <- rowMeans(boot_results, na.rm = TRUE)
	boot_ci <- t(apply(boot_results, 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))
	boot_ci <- as.data.table(boot_ci)
	setnames(boot_ci, old = c(1, 2), new = c("boot_CI_2.5", "boot_CI_97.5"))

	obs_diff[, boot_mean_log2FD := boot_mean]
	obs_diff <- cbind(obs_diff, boot_ci)

	## Store results and return object.
	sc_utils_obj@results$permutation <- obs_diff
	return(sc_utils_obj)
}

#' Plot Permutation Results
#'
#' @importFrom forcats fct_reorder
#' @importFrom dplyr desc case_when
#' @importFrom ggplot2 ggplot aes geom_pointrange theme_bw geom_hline
#' coord_flip scale_color_manual
#'
#' @param sc_utils_obj sc_utils object
#' @param FDR_threshold FDR value cutoff for significance
#' @param log2FD_threshold Absolute value of log2FD cutoff for significance
#' @param order_clusters Whether to order the clusters by observed log2FD
#'
#' @rdname permutation_plot-function
#'
#' @export

permutation_plot <- function(
	plot_data,
	FDR_threshold = 0.05,
	log2FD_threshold = log2(1.5),
	order_clusters = TRUE,
  scale_color = c("salmon", "grey")
) {

	## Retrieve results.
	#plot_data <- copy(sc_utils_obj@results$permutation)
	plot_data <- data.table::as.data.table(plot_data)

	## Mark the significant results.
	plot_data[, significance := ifelse(
		FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold,
		paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
		"n.s."
	)]

	plot_data[, significance := factor(significance, levels = c(
		paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
		"n.s."
	))]

	## Order the clusters by observed log2FD if requested.
	if (order_clusters) {
		plot_data[, clusters := forcats::fct_reorder(factor(clusters), desc(obs_log2FD))]
	}

	## Plot the results.
	p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
		geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) +
		theme_bw() +
		geom_hline(yintercept = log2FD_threshold, lty = 2) +
		geom_hline(yintercept = -log2FD_threshold, lty = 2) +
		geom_hline(yintercept = 0) +
		scale_color_manual(values = scale_color) +
		coord_flip()

	return(p)
}
