#!/usr/bin/env Rscript

# Preview a file (commonly under test/save/meta) with a short summary.
#
# Usage:
#   Rscript tools/preview_meta_file.R <file_path> [n]
#
# Examples:
#   Rscript tools/preview_meta_file.R test/save/meta/scrna_phase_meta.Rds
#   Rscript tools/preview_meta_file.R test/save/meta/scrna_phase_meta.Rds 5

suppressPackageStartupMessages({
  library(archive)
})

load_object <- function(file_name) {
  con <- archive::file_read(file = file_name)
  on.exit(close(con), add = TRUE)
  readRDS(file = con)
}

fmt_bytes <- function(x) {
  if (is.na(x)) return("NA")
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- 1
  val <- as.numeric(x)
  while (val >= 1024 && i < length(units)) {
    val <- val / 1024
    i <- i + 1
  }
  sprintf("%.2f %s", val, units[i])
}

print_df_preview <- function(df, n) {
  cat(sprintf("Rows: %d, Cols: %d\n", nrow(df), ncol(df)))
  show_n <- min(n, nrow(df))
  if (show_n > 0) {
    print(utils::head(df, show_n), row.names = FALSE)
  } else {
    cat("[empty data.frame]\n")
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript tools/preview_meta_file.R <file_path> [n]")
}

file_path <- args[[1]]
n <- if (length(args) >= 2) as.integer(args[[2]]) else 6L
if (is.na(n) || n < 1) n <- 6L

if (!file.exists(file_path)) {
  stop("File not found: ", file_path)
}

ext <- tolower(tools::file_ext(file_path))
info <- file.info(file_path)
preview_mode <- if (ext %in% c("rds", "rda", "rdata")) {
  "R object"
} else if (ext %in% c("csv", "tsv", "txt")) {
  "table/text"
} else {
  "generic text"
}

cat("===== File Summary =====\n")
cat("Path:      ", normalizePath(file_path, mustWork = FALSE), "\n", sep = "")
cat("Extension: ", ifelse(nchar(ext) == 0, "[none]", ext), "\n", sep = "")
cat("Size:      ", fmt_bytes(info$size), "\n", sep = "")
cat("Modified:  ", format(info$mtime, usetz = TRUE), "\n", sep = "")
cat("Preview:   ", preview_mode, " (top ", n, ")\n", sep = "")
cat("========================\n\n")

if (ext %in% c("rds", "rda", "rdata")) {
  obj <- tryCatch(load_object(file_path), error = function(e) NULL)
  if (is.null(obj)) {
    # Fallback for plain, non-archive compressed/non-compressed RDS
    obj <- readRDS(file_path)
  }

  cat("Type:", class(obj)[1], "\n")

  if (is.data.frame(obj)) {
    print_df_preview(obj, n)
  } else if (is.list(obj)) {
    cat("List length:", length(obj), "\n")
    cat("First names:\n")
    print(utils::head(names(obj), n))
    first_non_null <- NULL
    for (i in seq_along(obj)) {
      if (!is.null(obj[[i]])) {
        first_non_null <- obj[[i]]
        break
      }
    }
    if (is.data.frame(first_non_null)) {
      cat("\nFirst non-null element preview:\n")
      print_df_preview(first_non_null, n)
    }
  } else {
    utils::str(obj, max.level = 1)
  }

  quit(save = "no", status = 0)
}

if (ext %in% c("csv", "tsv", "txt")) {
  sep <- if (ext == "tsv") "\t" else ","
  if (ext == "txt") {
    lines <- readLines(file_path, warn = FALSE)
    cat(sprintf("Total lines: %d\n", length(lines)))
    print(utils::head(lines, n))
  } else {
    df <- utils::read.table(file_path, sep = sep, header = TRUE, check.names = FALSE)
    print_df_preview(df, n)
  }
  quit(save = "no", status = 0)
}

# Generic fallback: print first lines as text
lines <- readLines(file_path, warn = FALSE)
cat(sprintf("Total lines: %d\n", length(lines)))
print(utils::head(lines, n))
