#!/usr/bin/env Rscript

# Check pathway/enrichment outputs produced by data_factory functions:
# - generate_scrna_go
# - generate_scrna_kegg
# - generate_scrna_reactome
# - generate_scrna_hallmark
#
# Usage:
#   Rscript tools/check_pathway_outputs.R <scrna_rds> [preview_n]
#
# Example:
#   Rscript tools/check_pathway_outputs.R ./save/scrna_phase_comparing.Rds 5

suppressPackageStartupMessages({
  library(archive)
})

load_object <- function(file_name) {
  con <- archive::file_read(file = file_name)
  on.exit(close(con), add = TRUE)
  readRDS(file = con)
}

load_store_list <- function(value) {
  if (is.list(value)) {
    return(value)
  }
  if (is.character(value) && length(value) == 1) {
    if (!file.exists(value)) {
      stop(sprintf("Cannot find output file for value: %s", value))
    }
    return(load_object(value))
  }
  stop("Unsupported storage type in scrna@tools.")
}

count_terms <- function(x) {
  if (is.null(x)) {
    return(0L)
  }
  out <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (is.null(out)) {
    return(NA_integer_)
  }
  nrow(out)
}

summarize_one_direction <- function(lst) {
  if (is.null(lst) || !is.list(lst)) {
    return(list(cluster_n = 0L, non_null_n = 0L, total_terms = 0L, unknown_n = 0L))
  }
  term_counts <- vapply(lst, count_terms, integer(1))
  list(
    cluster_n = length(lst),
    non_null_n = sum(!is.na(term_counts) & term_counts > 0),
    total_terms = sum(term_counts[!is.na(term_counts)]),
    unknown_n = sum(is.na(term_counts))
  )
}

get_first_non_empty_df <- function(lst) {
  if (is.null(lst) || !is.list(lst)) {
    return(NULL)
  }
  for (nm in names(lst)) {
    df <- tryCatch(as.data.frame(lst[[nm]]), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      return(list(cluster = nm, df = df))
    }
  }
  NULL
}

print_preview <- function(tag, lst, preview_n) {
  first <- get_first_non_empty_df(lst)
  if (is.null(first)) {
    cat(sprintf("  %s preview: no non-empty results\n", tag))
    return(invisible(NULL))
  }
  keep_cols <- intersect(
    c("ID", "Description", "pvalue", "p.adjust", "qvalue", "geneID", "Count"),
    colnames(first$df)
  )
  if (length(keep_cols) == 0) {
    keep_cols <- colnames(first$df)[1:min(6, ncol(first$df))]
  }
  df_show <- first$df[, keep_cols, drop = FALSE]
  n_show <- min(preview_n, nrow(df_show))
  cat(sprintf("  %s preview (cluster=%s, top %d rows):\n", tag, first$cluster, n_show))
  print(utils::head(df_show, n_show), row.names = FALSE)
}

print_summary <- function(method_name, store_list, preview_n = 5L) {
  if (!is.list(store_list) || length(store_list) < 2) {
    cat(sprintf("[WARN] %s: unexpected object structure\n", method_name))
    return(invisible(NULL))
  }

  up <- store_list[[1]]
  down <- store_list[[2]]

  up_s <- summarize_one_direction(up)
  down_s <- summarize_one_direction(down)

  cat(sprintf("\n[%s]\n", method_name))
  cat(sprintf(
    "  up:   clusters=%d, with_terms=%d, total_terms=%d, unknown=%d\n",
    up_s$cluster_n, up_s$non_null_n, up_s$total_terms, up_s$unknown_n
  ))
  cat(sprintf(
    "  down: clusters=%d, with_terms=%d, total_terms=%d, unknown=%d\n",
    down_s$cluster_n, down_s$non_null_n, down_s$total_terms, down_s$unknown_n
  ))
  print_preview("up", up, preview_n = preview_n)
  print_preview("down", down, preview_n = preview_n)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript tools/check_pathway_outputs.R <scrna_rds> [preview_n]")
}

scrna_rds <- args[[1]]
preview_n <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
if (is.na(preview_n) || preview_n < 1) {
  preview_n <- 5L
}

if (!file.exists(scrna_rds)) {
  stop("scrna_rds not found: ", scrna_rds)
}

scrna <- load_object(scrna_rds)
if (!("tools" %in% slotNames(scrna))) {
  stop("Input object does not have @tools slot.")
}

tool_names <- names(scrna@tools)
targets <- tool_names[grepl("^(go|kegg|reactome|hallmark)_", tool_names)]

if (length(targets) == 0) {
  cat("No go_/kegg_/reactome_/hallmark_ outputs found in scrna@tools.\n")
  quit(save = "no", status = 0)
}

cat("Found pathway outputs:\n")
for (nm in targets) {
  cat(" - ", nm, "\n", sep = "")
}

for (nm in targets) {
  val <- scrna@tools[[nm]]
  store_list <- tryCatch(
    load_store_list(val),
    error = function(e) e
  )

  if (inherits(store_list, "error")) {
    cat(sprintf("\n[ERROR] %s: %s\n", nm, store_list$message))
    next
  }

  print_summary(nm, store_list, preview_n = preview_n)
}

cat("\nDone.\n")
