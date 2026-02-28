#!/usr/bin/env Rscript

# Preview all files under a partition folder (default: test/save/partition).
#
# Usage:
#   Rscript tools/preview_partition_folder.R [options]
#
# Examples:
#   Rscript tools/preview_partition_folder.R
#   Rscript tools/preview_partition_folder.R -d test/save/partition -n 5
#   Rscript tools/preview_partition_folder.R --help

suppressPackageStartupMessages({
  library(archive)
  library(optparse)
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

print_rds_preview <- function(file_path, n) {
  obj <- tryCatch(load_object(file_path), error = function(e) NULL)
  if (is.null(obj)) {
    obj <- tryCatch(readRDS(file_path), error = function(e) NULL)
  }

  if (is.null(obj)) {
    cat("Type: [unreadable R object]\n")
    return(invisible(NULL))
  }

  cat("Type:", class(obj)[1], "\n")

  if (is.data.frame(obj)) {
    print_df_preview(obj, n)
    return(invisible(NULL))
  }

  if (is.list(obj)) {
    cat("List length:", length(obj), "\n")
    nm <- names(obj)
    if (!is.null(nm)) {
      cat("First names:\n")
      print(utils::head(nm, n))
    }
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
    return(invisible(NULL))
  }

  utils::str(obj, max.level = 1)
}

print_text_preview <- function(file_path, n) {
  lines <- readLines(file_path, warn = FALSE)
  cat(sprintf("Total lines: %d\n", length(lines)))
  print(utils::head(lines, n))
}

preview_one_file <- function(file_path, n) {
  info <- file.info(file_path)
  ext <- tolower(tools::file_ext(file_path))
  mode <- if (ext %in% c("rds", "rda", "rdata")) "R object" else "text"

  cat("\n===== File Summary =====\n")
  cat("Path:      ", normalizePath(file_path, mustWork = FALSE), "\n", sep = "")
  cat("Extension: ", ifelse(nchar(ext) == 0, "[none]", ext), "\n", sep = "")
  cat("Size:      ", fmt_bytes(info$size), "\n", sep = "")
  cat("Modified:  ", format(info$mtime, usetz = TRUE), "\n", sep = "")
  cat("Preview:   ", mode, " (top ", n, ")\n", sep = "")
  cat("========================\n")

  if (ext %in% c("rds", "rda", "rdata")) {
    print_rds_preview(file_path, n)
  } else {
    print_text_preview(file_path, n)
  }
}

option_list <- list(
  make_option(
    c("-d", "--partition_dir"),
    type = "character",
    default = "test/save/partition",
    help = "Partition folder to preview [default %default]"
  ),
  make_option(
    c("-n", "--n"),
    type = "integer",
    default = 6L,
    help = "Top rows/lines to preview per file [default %default]"
  )
)

parser <- OptionParser(
  usage = "Rscript %prog [options]",
  option_list = option_list,
  description = "Preview all files under a partition folder."
)
opt <- parse_args(parser)

partition_dir <- opt$partition_dir
n <- as.integer(opt$n)
if (is.na(n) || n < 1) n <- 6L

if (!dir.exists(partition_dir)) {
  stop("partition_dir does not exist: ", partition_dir)
}

files <- list.files(partition_dir, full.names = TRUE, recursive = FALSE, all.files = FALSE)
files <- files[file.info(files)$isdir == FALSE]

cat("Partition dir:", normalizePath(partition_dir, mustWork = FALSE), "\n")
cat("Total files:", length(files), "\n")

if (length(files) == 0) {
  cat("No files found.\n")
  quit(save = "no", status = 0)
}

for (fp in files) {
  preview_one_file(fp, n)
}

cat("\nDone.\n")
