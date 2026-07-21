# Utility script to subsample posterior draws for the Shiny app
# Adjust the three parameters below as needed before running.

original_draws_dir <- "/Thesis/3rdChapter/PEP_QC/results/undisturbed/preds"
subsampled_draws_dir <- "/Thesis/3rdChapter/PEP_QC/posterior_draws_200"
K <- 400L

# ---- Helpers ---------------------------------------------------------------
ensure_directory <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path)) {
    stop(sprintf("Unable to create output directory: %s", path))
  }
}

load_preds <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("preds", envir = env, inherits = FALSE)) {
    stop("Object 'preds' not found in Rda file")
  }
  preds <- get("preds", envir = env, inherits = FALSE)
  if (!is.list(preds) || length(preds) == 0L) {
    stop("'preds' must be a non-empty list")
  }
  if (!all(vapply(preds, is.matrix, logical(1L)))) {
    stop("All elements in 'preds' must be matrices")
  }
  first_mat <- preds[[1L]]
  target_dim <- dim(first_mat)
  target_dimnames <- dimnames(first_mat)
  dims_ok <- vapply(preds, function(mat) {
    identical(dim(mat), target_dim) && identical(dimnames(mat), target_dimnames)
  }, logical(1L))
  if (!all(dims_ok)) {
    stop("All matrices in 'preds' must share identical dimensions and dimnames")
  }
  preds
}

subsample_preds <- function(preds, K) {
  n_draws <- length(preds)
  keep_n <- min(K, n_draws)
  if (keep_n < n_draws) {
    # Sample without replacement and sort to keep chronological order minimal
    idx <- sort(sample.int(n_draws, keep_n, replace = FALSE))
  } else {
    idx <- seq_len(n_draws)
  }
  preds[idx]
}

# ---- Main -----------------------------------------------------------------
if (!dir.exists(original_draws_dir)) {
  stop(sprintf("Input directory not found: %s", original_draws_dir))
}

ensure_directory(subsampled_draws_dir)

rda_files <- list.files(original_draws_dir, pattern = "_preds\\.Rda$", full.names = TRUE)
if (length(rda_files) == 0L) {
  stop(sprintf("No '*_preds.Rda' files found in %s", original_draws_dir))
}

for (rda_path in rda_files) {
  scenario_label <- sub("_preds\\.Rda$", "", basename(rda_path))
  message(sprintf("Processing scenario '%s'", scenario_label))
  result <- tryCatch({
    preds_full <- load_preds(rda_path)
    preds_sub <- subsample_preds(preds_full, K)
    preds <- preds_sub
    out_name <- sprintf("%s_preds_subsampled.Rda", scenario_label)
    out_path <- file.path(subsampled_draws_dir, out_name)
    save(preds, file = out_path)
    message(sprintf("  draws found: %d | kept: %d | saved to: %s",
                    length(preds_full), length(preds_sub), out_name))
    NULL
  }, error = function(e) {
    message(sprintf("  skipped due to error: %s", e$message))
    NULL
  })
}

message("Subsampling complete.")
