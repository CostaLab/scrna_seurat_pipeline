#!/usr/bin/env Rscript

# Download KEGG tables and build offline clusterProfiler caches.
# Usage:
#   Rscript tools/prepare_kegg_offline_cache.R [species] [base_dir]
# species:
#   - all (default): build both hsa + mmu
#   - hsa
#   - mmu
# Example:
#   Rscript tools/prepare_kegg_offline_cache.R all ../KEGG_DB
#   # files are saved to ../KEGG_DB

args <- commandArgs(trailingOnly = TRUE)
species_arg <- if (length(args) >= 1) tolower(args[[1]]) else "all"
base_dir <- if (length(args) >= 2) args[[2]] else "../"
output_dir <- file.path(base_dir, "KEGG_DB")

species_codes <- switch(
  species_arg,
  all = c("hsa", "mmu"),
  hsa = c("hsa"),
  mmu = c("mmu"),
  stop("species must be one of: all, hsa, mmu")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

download_if_missing <- function(url, destfile) {
  if (file.exists(destfile) && file.size(destfile) > 0) {
    message("Using existing file: ", normalizePath(destfile, mustWork = FALSE))
    return(invisible(NULL))
  }
  message("Downloading: ", url)
  download.file(url, destfile = destfile, mode = "wb")
}

build_one_species <- function(species_code, out_dir) {
  pathway_list_url <- sprintf("https://rest.kegg.jp/list/pathway/%s", species_code)
  pathway_link_url <- sprintf("https://rest.kegg.jp/link/%s/pathway", species_code)
  gene_info_url <- sprintf("https://rest.kegg.jp/list/%s", species_code)

  pathway_list_file <- file.path(out_dir, sprintf("pathway_%s.list.tsv", species_code))
  pathway_link_file <- file.path(out_dir, sprintf("pathway2gene_%s.tsv", species_code))
  gene_info_file <- file.path(out_dir, sprintf("%s_gene_info.tsv", species_code))

  message("=== Processing ", species_code, " ===")
  download_if_missing(pathway_list_url, pathway_list_file)
  download_if_missing(pathway_link_url, pathway_link_file)
  download_if_missing(gene_info_url, gene_info_file)

  term2name <- read.delim(
    pathway_list_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  colnames(term2name) <- c("term", "name")

  species_label <- switch(
    species_code,
    hsa = "Homo sapiens (human)",
    mmu = "Mus musculus (mouse)"
  )
  term2name$name <- sub(paste0(" - ", species_label, "$"), "", term2name$name)
  term2name <- unique(term2name)

  term2gene <- read.delim(
    pathway_link_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  colnames(term2gene) <- c("term", "gene")
  term2gene$term <- sub("^path:", "", term2gene$term)
  term2gene$gene <- sub(paste0("^", species_code, ":"), "", term2gene$gene)
  term2gene <- unique(term2gene)

  term2gene_rds <- file.path(out_dir, sprintf("kegg_term2gene_%s.rds", species_code))
  term2name_rds <- file.path(out_dir, sprintf("kegg_term2name_%s.rds", species_code))

  saveRDS(term2gene, term2gene_rds)
  saveRDS(term2name, term2name_rds)

  message("Saved files for ", species_code, ":")
  message(" - ", normalizePath(pathway_list_file, mustWork = FALSE))
  message(" - ", normalizePath(pathway_link_file, mustWork = FALSE))
  message(" - ", normalizePath(gene_info_file, mustWork = FALSE))
  message(" - ", normalizePath(term2gene_rds, mustWork = FALSE))
  message(" - ", normalizePath(term2name_rds, mustWork = FALSE))
}

for (sp in species_codes) {
  build_one_species(sp, output_dir)
}

message("Done.")
