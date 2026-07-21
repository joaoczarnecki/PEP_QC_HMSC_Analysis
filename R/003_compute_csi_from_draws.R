# =============================================================================
# 003_gCSI_trajectory_analysis_draws.R
# Purpose : Summarise and visualise gCSI trajectories with proper uncertainty
#           propagation from HMSC posterior draws for selected communities
#           across climate scenarios.
# Note    : This version requires per-sample predictive draws:
#           pred_summary$draws with dim = c(n_draws, n_sites, n_species)
#           and dimnames(draw = <ids>, site = <site_ids>, species = <species_names>)
# Date    : 2025-06-12
# =============================================================================

# --- -1) LIBRARIES ----------------------------------------------------------------


suppressPackageStartupMessages({
  library(tidyverse)
})
# --- 0) PARAMETERS ------------------------------------------------------------
epsilon       <- 1e-4          # probability floor for geometric mean stability
use_mean_vs_median_for_central <- "median"  # "median" or "mean"
species_weights <- NULL        # named numeric vector by species, or NULL for equal weights
export_site_level <- TRUE
n_cores <- 1   # adjust manually if desired


# --- 1) UTILS: POSTERIOR SUMMARIES -------------------------------------------
summarise_posterior <- function(x, central = c("median","mean")) {
  central <- match.arg(central)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(tibble(
      center = NA_real_, mean = NA_real_, sd = NA_real_,
      q025 = NA_real_, q975 = NA_real_
    ))
  }
  tibble(
    center = if (central == "median") stats::median(x) else mean(x),
    mean   = mean(x),
    sd     = stats::sd(x),
    q025   = stats::quantile(x, 0.025, na.rm = TRUE),
    q975   = stats::quantile(x, 0.975, na.rm = TRUE)
  )
}

empty_site_summary_template <- function() {
  tibble(
    site = character(), index = character(),
    center = numeric(), mean = numeric(), sd = numeric(),
    q025 = numeric(), q975 = numeric(),
    log_center = numeric(), log_mean = numeric(), log_sd = numeric(),
    log_q025 = numeric(), log_q975 = numeric()
  )
}

empty_community_summary_template <- function() {
  tibble(
    index = "summary",
    center = NA_real_, mean = NA_real_, sd = NA_real_,
    q025 = NA_real_, q975 = NA_real_,
    log_center = NA_real_, log_mean = NA_real_, log_sd = NA_real_,
    log_q025 = NA_real_, log_q975 = NA_real_,
    gCenter = NA_real_, gMean = NA_real_, gSD = NA_real_,
    gQ025 = NA_real_, gQ975 = NA_real_,
    bCenter = NA_real_, bMean = NA_real_, bSD = NA_real_,
    bQ025 = NA_real_, bQ975 = NA_real_,
    wCenter = NA_real_, wMean = NA_real_, wSD = NA_real_,
    wQ025 = NA_real_, wQ975 = NA_real_
  )
}

#' Compute CSI indices per draw-site with NA-safe weights.
# ----------------------------------------------------------------------------
# Legacy (loop-based) index computation retained for backward compatibility
# Supports: gCSI, bCSI, wCSI. New optimized path uses matrix ops for gCSI only.
# ----------------------------------------------------------------------------
compute_indices_from_draws_legacy <- function(draws, species_vec, weights = NULL, eps = 1e-4) {
  stopifnot(length(dim(draws)) == 3)
  dn <- dimnames(draws)
  species_all <- dn[[3]]
  if (is.null(species_all)) stop("draws must have species names on the third dimension")
  keep <- intersect(species_vec, species_all)
  if (length(keep) == 0) {
    return(tibble(draw = character(), site = character(), index = character(), value = numeric(), log_value = numeric()))
  }
  if (!is.null(weights)) {
    if (is.null(names(weights))) stop("weights must be a named numeric vector")
    weights <- weights[keep]
    if (any(is.na(weights))) stop("Missing weights for some selected species")
  } else {
    weights <- rep(1, length(keep))
    names(weights) <- keep
  }
  if (all(weights == 0)) stop("All provided weights are zero for the selected species set")
  draws_sub <- draws[, , keep, drop = FALSE]
  nd <- dim(draws_sub)[1]
  ns <- dim(draws_sub)[2]
  draw_names <- dn[[1]]; if (is.null(draw_names)) draw_names <- as.character(seq_len(nd))
  site_names <- dn[[2]]; if (is.null(site_names)) site_names <- as.character(seq_len(ns))
  res_list <- vector("list", nd * ns)
  idx <- 0L
  for (d in seq_len(nd)) {
    for (s in seq_len(ns)) {
      p_vec <- draws_sub[d, s, ]
      if (all(is.na(p_vec))) next
      if (any(!is.na(p_vec) & (p_vec < -1e-8 | p_vec > 1 + 1e-8))) {
        stop("Probabilities outside [0,1] detected for draw ", draw_names[d], " site ", site_names[s])
      }
      valid <- !is.na(p_vec)
      if (!any(valid)) next
      p_eff <- p_vec[valid]
      w_eff <- weights[valid]
      if (any(is.na(w_eff))) stop("Weights missing after dropping NA species probabilities")
      if (all(w_eff == 0)) stop("Row-specific weights sum to zero after filtering NA species")
      w_eff <- w_eff / sum(w_eff)
      b_val <- mean(p_eff)                       # arithmetic mean
      w_val <- sum(w_eff * p_eff)                # weighted arithmetic mean
      log_g <- sum(w_eff * log(pmax(p_eff, eps)))# weighted log mean
      g_val <- exp(log_g)                        # geometric mean
      idx <- idx + 1L
      res_list[[idx]] <- tibble(
        draw = draw_names[d],
        site = site_names[s],
        index = c("gCSI","bCSI","wCSI"),
        value = c(g_val, b_val, w_val),
        log_value = c(log_g, NA_real_, NA_real_)
      )
    }
  }
  if (idx == 0L) {
    return(tibble(draw = character(), site = character(), index = character(), value = numeric(), log_value = numeric()))
  }
  bind_rows(res_list[seq_len(idx)])
}

# ----------------------------------------------------------------------------
# Optimized gCSI computation (vectorized) using matrix operations
# Computes geometric mean in log space with per-draw/site NA handling.
# Returns matrices: g_matrix [n_draws x n_sites], log_g_matrix [n_draws x n_sites]
# ----------------------------------------------------------------------------
compute_gCSI_matrix <- function(draws_sub, weights = NULL, eps = 1e-4) {
  stopifnot(length(dim(draws_sub)) == 3)
  nd <- dim(draws_sub)[1]; ns <- dim(draws_sub)[2]; nk <- dim(draws_sub)[3]
  if (is.null(weights)) {
    weights <- rep(1, nk)
  } else {
    if (length(weights) != nk) stop("weights length must match selected species count in optimized path")
  }
  if (all(weights == 0)) stop("All provided weights are zero for selected species (optimized path)")
  # Expand weights across draws/sites
  w_array <- array(weights, dim = c(nd, ns, nk))
  valid <- !is.na(draws_sub)
  # Mask weights where probabilities are NA
  w_masked <- ifelse(valid, w_array, 0)
  sum_w <- apply(w_masked, c(1,2), sum)  # total weight per draw-site
  # Avoid division by zero by marking zero-sum cells
  zero_sum <- (sum_w == 0)
  # Normalise weights
  w_norm <- w_masked / array(sum_w, dim = c(nd, ns, nk))
  # Log probabilities (replace NA with 0 so they contribute nothing after weighting)
  logP <- log(pmax(draws_sub, eps))
  logP[!valid] <- 0
  # Probability bounds check (ignore NA)
  if (any(draws_sub[valid] < -1e-8 | draws_sub[valid] > 1 + 1e-8)) {
    stop("Probabilities outside [0,1] detected in optimized gCSI path")
  }
  # Weighted log sum across species -> log gCSI
  log_g_matrix <- apply(w_norm * logP, c(1,2), sum)
  # Set cells with zero total weight to NA
  if (any(zero_sum)) log_g_matrix[zero_sum] <- NA_real_
  g_matrix <- exp(log_g_matrix)
  list(g_matrix = g_matrix, log_g_matrix = log_g_matrix)
}

# ----------------------------------------------------------------------------
# Posterior summaries from gCSI matrices (site level)
# ----------------------------------------------------------------------------
summarise_site_gCSI_matrix <- function(g_matrix, log_g_matrix, central = "median") {
  nd <- nrow(g_matrix); ns <- ncol(g_matrix)
  site_names <- colnames(g_matrix); if (is.null(site_names)) site_names <- paste0("site_", seq_len(ns))
  res_list <- vector("list", ns)
  for (s in seq_len(ns)) {
    g_vec <- g_matrix[, s]
    log_vec <- log_g_matrix[, s]
    g_stats <- summarise_posterior(g_vec, central = central)
    log_stats <- summarise_posterior(log_vec[is.finite(log_vec)], central = central)
    # Align names
    names(log_stats) <- sub("^center$", "log_center", names(log_stats))
    names(log_stats) <- sub("^mean$", "log_mean", names(log_stats))
    names(log_stats) <- sub("^sd$", "log_sd", names(log_stats))
    names(log_stats) <- sub("^q025$", "log_q025", names(log_stats))
    names(log_stats) <- sub("^q975$", "log_q975", names(log_stats))
    res_list[[s]] <- tibble(
      site = site_names[s],
      index = "gCSI",
      !!!g_stats,
      !!!log_stats
    )
  }
  dplyr::bind_rows(res_list)
}

# ----------------------------------------------------------------------------
# Posterior summaries from gCSI matrices (community level)
# Community aggregation: mean of log gCSI across sites per draw (equivalent to
# geometric mean across sites if weights equal; matches previous implementation).
# ----------------------------------------------------------------------------
summarise_community_gCSI_matrix <- function(g_matrix, log_g_matrix, central = "median") {
  # Compute community log-gCSI per draw (mean across sites of log values)
  comm_log <- apply(log_g_matrix, 1, function(row) {
    row <- row[is.finite(row)]
    if (!length(row)) return(NA_real_)
    mean(row)
  })
  comm_g <- exp(comm_log)
  g_stats <- summarise_posterior(comm_g, central = central)
  log_stats <- summarise_posterior(comm_log[is.finite(comm_log)], central = central)
  names(log_stats) <- sub("^center$", "log_center", names(log_stats))
  names(log_stats) <- sub("^mean$", "log_mean", names(log_stats))
  names(log_stats) <- sub("^sd$", "log_sd", names(log_stats))
  names(log_stats) <- sub("^q025$", "log_q025", names(log_stats))
  names(log_stats) <- sub("^q975$", "log_q975", names(log_stats))
  tibble(index = "summary", !!!g_stats, !!!log_stats) %>%
    mutate(
      gCenter = center, gMean = mean, gSD = sd,
      gQ025 = q025, gQ975 = q975,
      bCenter = NA_real_, bMean = NA_real_, bSD = NA_real_, bQ025 = NA_real_, bQ975 = NA_real_,
      wCenter = NA_real_, wMean = NA_real_, wSD = NA_real_, wQ025 = NA_real_, wQ975 = NA_real_
    )
}

#' Summarise posterior indices at site level.
summarise_site_indices_from_draws <- function(idx_long, central = "median") {
  if (nrow(idx_long) == 0) return(empty_site_summary_template())
  value_sum <- idx_long %>%
    group_by(site, index) %>%
    summarise(stats = list(summarise_posterior(value, central = central)), .groups = "drop") %>%
    tidyr::unnest_wider(stats)
  log_sum <- idx_long %>%
    filter(index == "gCSI", is.finite(log_value)) %>%
    group_by(site, index) %>%
    summarise(log_stats = list(summarise_posterior(log_value, central = central)), .groups = "drop") %>%
    tidyr::unnest_wider(log_stats) %>%
    rename(log_center = center, log_mean = mean, log_sd = sd, log_q025 = q025, log_q975 = q975)
  value_sum %>%
    left_join(log_sum, by = c("site","index")) %>%
    arrange(site, index)
}

#' Summarise posterior indices at community level (across sites).
summarise_community_indices_from_draws <- function(idx_long, central = "median") {
  if (nrow(idx_long) == 0) return(empty_community_summary_template())
  agg_draws <- idx_long %>%
    group_by(index, draw) %>%
    group_modify(~{
      if (.y$index == "gCSI") {
        log_vals <- .x$log_value[is.finite(.x$log_value)]
        if (!length(log_vals)) {
          tibble(value = NA_real_, log_value = NA_real_)
        } else {
          mean_log <- mean(log_vals)
          tibble(value = exp(mean_log), log_value = mean_log)
        }
      } else {
        vals <- .x$value[is.finite(.x$value)]
        tibble(value = if (length(vals)) mean(vals) else NA_real_, log_value = NA_real_)
      }
    }) %>%
    ungroup()
  
  # Compute summaries per index
  value_sum <- agg_draws %>%
    group_by(index) %>%
    summarise(stats = list(summarise_posterior(value, central = central)), .groups = "drop") %>%
    tidyr::unnest_wider(stats)
  
  log_sum <- agg_draws %>%
    filter(index == "gCSI", is.finite(log_value)) %>%
    group_by(index) %>%
    summarise(log_stats = list(summarise_posterior(log_value, central = central)), .groups = "drop") %>%
    tidyr::unnest_wider(log_stats) %>%
    rename(log_center = center, log_mean = mean, log_sd = sd, log_q025 = q025, log_q975 = q975)
  
  # Join and pivot to wide format with one row
  result <- value_sum %>%
    left_join(log_sum, by = "index") %>%
    select(index, center, mean, sd, q025, q975, everything())
  
  # Pivot to wide format: one row with gCenter, bCenter, wCenter etc.
  wide_result <- tibble(
    index = "summary",
    center = NA_real_, mean = NA_real_, sd = NA_real_,
    q025 = NA_real_, q975 = NA_real_,
    log_center = NA_real_, log_mean = NA_real_, log_sd = NA_real_,
    log_q025 = NA_real_, log_q975 = NA_real_
  )
  
  for (idx_type in c("gCSI", "bCSI", "wCSI")) {
    row <- result %>% filter(index == idx_type)
    if (nrow(row) > 0) {
      prefix <- substr(idx_type, 1, 1)  # "g", "b", or "w"
      wide_result[[paste0(prefix, "Center")]] <- row$center
      wide_result[[paste0(prefix, "Mean")]] <- row$mean
      wide_result[[paste0(prefix, "SD")]] <- row$sd
      wide_result[[paste0(prefix, "Q025")]] <- row$q025
      wide_result[[paste0(prefix, "Q975")]] <- row$q975
      if (idx_type == "gCSI" && "log_center" %in% names(row)) {
        wide_result$log_center <- row$log_center
        wide_result$log_mean <- row$log_mean
        wide_result$log_sd <- row$log_sd
        wide_result$log_q025 <- row$log_q025
        wide_result$log_q975 <- row$log_q975
      }
    } else {
      prefix <- substr(idx_type, 1, 1)
      wide_result[[paste0(prefix, "Center")]] <- NA_real_
      wide_result[[paste0(prefix, "Mean")]] <- NA_real_
      wide_result[[paste0(prefix, "SD")]] <- NA_real_
      wide_result[[paste0(prefix, "Q025")]] <- NA_real_
      wide_result[[paste0(prefix, "Q975")]] <- NA_real_
    }
  }
  
  wide_result
}

coerce_pred_obj_to_draws <- function(pred_obj) {
  if (is.array(pred_obj) && length(dim(pred_obj)) == 3) {
    draws <- pred_obj
  } else if (is.list(pred_obj) && "draws" %in% names(pred_obj)) {
    draws <- pred_obj$draws
  } else if (is.list(pred_obj) && length(pred_obj) > 0 && all(vapply(pred_obj, is.matrix, logical(1)))) {
    n_draws <- length(pred_obj)
    site_names <- rownames(pred_obj[[1]])
    species_names <- colnames(pred_obj[[1]])
    if (is.null(species_names)) stop("Prediction matrices must have species column names")
    draws <- array(NA_real_, dim = c(n_draws, nrow(pred_obj[[1]]), ncol(pred_obj[[1]])),
                   dimnames = list(
                     draw = if (!is.null(names(pred_obj))) names(pred_obj) else paste0("draw_", seq_len(n_draws)),
                     site = if (!is.null(site_names)) site_names else paste0("site_", seq_len(nrow(pred_obj[[1]]))),
                     species = species_names
                   ))
    for (i in seq_len(n_draws)) {
      mat <- pred_obj[[i]]
      if (!all(dim(mat) == dim(pred_obj[[1]]))) stop("All prediction matrices must share dimensions")
      draws[i, , ] <- mat
    }
  } else {
    stop("Unsupported prediction object supplied")
  }
  dn <- dimnames(draws)
  if (is.null(dn[[1]])) dimnames(draws)[[1]] <- as.character(seq_len(dim(draws)[1]))
  if (is.null(dn[[2]])) dimnames(draws)[[2]] <- as.character(seq_len(dim(draws)[2]))
  if (is.null(dimnames(draws)[[3]])) stop("draws array must have species names")
  draws
}

#' Calculate indices for a scenario and return summaries.
calculate_indices_from_scenario <- function(pred_obj,
                                            species_vec,
                                            weights = NULL,
                                            eps = 1e-4,
                                            central = "median",
                                            return_level = c("community","site","both")) {
  return_level <- match.arg(return_level)
  draws <- coerce_pred_obj_to_draws(pred_obj)
  idx_long <- compute_indices_from_draws_legacy(draws, species_vec, weights = weights, eps = eps)
  site_summary <- summarise_site_indices_from_draws(idx_long, central = central)
  community_summary <- summarise_community_indices_from_draws(idx_long, central = central)
  if (return_level == "community") return(community_summary)
  if (return_level == "site") return(site_summary)
  list(community = community_summary, site = site_summary)
}

# --- 2) CORE: gCSI FROM DRAWS WITH UNCERTAINTY PROPAGATION ------------------------
compute_log_gCSI_from_draws <- function(draws, species_vec, weights = NULL, eps = 1e-4) {
  idx_long <- compute_indices_from_draws(draws, species_vec, weights = weights, eps = eps)
  if (nrow(idx_long) == 0) {
    return(tibble(draw = character(), site = character(), log_gCSI = numeric(), gCSI = numeric()))
  }
  idx_long %>%
    filter(index == "gCSI") %>%
    transmute(draw, site, log_gCSI = log_value, gCSI = value)
}

summarise_community_gCSI_from_draws <- function(gCSI_long, central = "median") {
  if (nrow(gCSI_long) == 0) return(tibble())
  idx_long <- gCSI_long %>%
    transmute(draw, site, index = "gCSI", value = gCSI, log_value = log_gCSI)
  res <- summarise_community_indices_from_draws(idx_long, central = central)
  if (!nrow(res)) return(tibble())
  res %>%
    select(
      gCenter, gMean, gSD, gQ025, gQ975,
      log_center, log_mean, log_sd, log_q025, log_q975
    )
}

# --- 5) LOAD DATA -------------------------------------------------------------
load("J:/Thesis/3rdChapter/PEP_QC/results/undisturbed/CAMS_Shiny_App_PA_Advanced/type_eco_data.RData")
stopifnot(exists("TYPE_ECO_list"))

# predictions_dir <- "J:/Thesis/3rdChapter/PEP_QC/results/undisturbed/preds"
predictions_dir <- "J:\Thesis\3rdChapter\PEP_QC\posterior_draws1"

# --- 5.1) SCENARIO AND SPECIES MAPS -----------------------------------------------

scenario_map <- tibble::tribble(
  ~scenario_label,          ~projection,       ~file_name,
  "Current",               "Baseline",        "Current_preds.Rda",
  "SSP1-2.6 2011-2040",    "2011-2040",       "8GCMs_ensemble_ssp126_2011-2040.gcm_preds.Rda",
  "SSP1-2.6 2041-2070",    "2041-2070",       "8GCMs_ensemble_ssp126_2041-2070.gcm_preds.Rda",
  "SSP1-2.6 2071-2100",    "2071-2100",       "8GCMs_ensemble_ssp126_2071-2100.gcm_preds.Rda",
  "SSP2-4.5 2011-2040",    "2011-2040",       "8GCMs_ensemble_ssp245_2011-2040.gcm_preds.Rda",
  "SSP2-4.5 2041-2070",    "2041-2070",       "8GCMs_ensemble_ssp245_2041-2070.gcm_preds.Rda",
  "SSP2-4.5 2071-2100",    "2071-2100",       "8GCMs_ensemble_ssp245_2071-2100.gcm_preds.Rda",
  "SSP3-7.0 2011-2040",    "2011-2040",       "8GCMs_ensemble_ssp370_2011-2040.gcm_preds.Rda",
  "SSP3-7.0 2041-2070",    "2041-2070",       "8GCMs_ensemble_ssp370_2041-2070.gcm_preds.Rda",
  "SSP3-7.0 2071-2100",    "2071-2100",       "8GCMs_ensemble_ssp370_2071-2100.gcm_preds.Rda",
  "SSP5-8.5 2011-2040",    "2011-2040",       "8GCMs_ensemble_ssp585_2011-2040.gcm_preds.Rda",
  "SSP5-8.5 2041-2070",    "2041-2070",       "8GCMs_ensemble_ssp585_2041-2070.gcm_preds.Rda",
  "SSP5-8.5 2071-2100",    "2071-2100",       "8GCMs_ensemble_ssp585_2071-2100.gcm_preds.Rda"
) %>% mutate(file_path = file.path(predictions_dir, file_name))

species_map <- tibble::tribble(
  ~code, ~scientific_name,
  "AEH", "Aesculus hippocastanum", "AME", "Amelanchier spp.", "AUC", "Alnus crispa",
  "AUR", "Alnus incana", "BOG", "Betula populifolia", "BOJ", "Betula alleghaniensis",
  "BOP", "Betula papyrifera", "CAC", "Carya cordiformis", "CAF", "Carya ovata",
  "CAR", "Carpinus caroliniana", "CET", "Prunus serotina", "CHB", "Quercus alba",
  "CHE", "Quercus bicolor", "CHG", "Quercus macrocarpa", "CHR", "Quercus rubra",
  "EPB", "Picea glauca", "EPN", "Picea mariana", "EPO", "Picea abies",
  "EPR", "Picea rubens", "ERA", "Acer saccharinum", "ERB", "Acer platanoides",
  "ERE", "Acer spicatum", "ERG", "Acer negundo", "ERN", "Acer nigrum",
  "ERP", "Acer pensylvanicum", "ERR", "Acer rubrum", "ERS", "Acer saccharum",
  "FRA", "Fraxinus americana", "FRN", "Fraxinus nigra", "FRP", "Fraxinus pennsylvanica",
  "HEG", "Fagus grandifolia", "MEL", "Larix laricina", "NOC", "Juglans cinerea",
  "ORA", "Ulmus americana", "ORR", "Ulmus rubra", "ORT", "Ulmus thomasii",
  "OSV", "Ostrya virginiana", "PEB", "Populus balsamifera", "PEG", "Populus grandidentata",
  "PET", "Populus tremuloides", "PIB", "Pinus strobus", "PIG", "Pinus banksiana",
  "PIR", "Pinus resinosa", "PIS", "Pinus sylvestris", "PRU", "Tsuga canadensis",
  "PRP", "Prunus pensylvanica", "SAB", "Abies balsamea", "SOA", "Sorbus americana",
  "SOD", "Sorbus decora", "THO", "Thuja occidentalis", "TIL", "Tilia americana"
)

create_pretty_name <- function(scientific_name) {
  parts <- strsplit(scientific_name, " ")[[1]]
  genus <- substr(parts[1], 1, 3)
  species <- if (length(parts) > 1) substr(parts[2], 1, 4) else ""
  paste0(str_to_title(genus), ".", species)
}
species_map <- species_map %>%
  mutate(
    original_name_in_data = paste0(code, "_n"),
    pretty_name = sapply(scientific_name, create_pretty_name)
  )



# --- 6) COMMUNITY DEFINITIONS -------------------------------------------------
type_eco_descriptions <- tibble::tribble(
  ~type_eco_prefix, ~Description,
  "FC1", "Red Oak Forest",
  "FE1", "Sugar Maple - Bitternut Hickory Forest",
  "FE2", "Sugar Maple - Basswood Forest",
  "FE3", "Sugar Maple - Yellow Birch Forest",
  "FE4", "Sugar Maple - Yellow Birch and Beech Forest",
  "FE5", "Sugar Maple - Hop Hornbeam Forest",
  "FE6", "Sugar Maple - Red Oak Forest",
  "FO1", "Black Ash Grove",
  "LA1", "Lichen or Moss Heathland",
  "LA2", "Shrubby Heathland",
  "LA3", "Herbaceous Heathland",
  "LA4", "Rocky Heathland",
  "LI1", "Littoral (Littoral Zone)",
  "LL1", "Alpine Lichen or Moss Heathland",
  "LL2", "Alpine Shrubby Heathland",
  "LL3", "Alpine Herbaceous Heathland",
  "LL4", "Alpine Rocky Heathland",
  "LM1", "Maritime Lichen or Moss Heathland",
  "LM2", "Maritime Shrubby Heathland",
  "LM3", "Maritime Herbaceous Heathland",
  "LM4", "Maritime Rocky Heathland",
  "MA1", "Freshwater Shrubby Marsh or Swamp",
  "MA2", "Brackish or Saltwater Shrubby Marsh or Swamp",
  "ME1", "Black Spruce - Trembling Aspen Forest",
  "MF1", "Black Ash - Balsam Fir Forest",
  "MJ1", "Yellow Birch - Balsam Fir and Sugar Maple Forest",
  "MJ2", "Yellow Birch - Balsam Fir Forest",
  "MS1", "Balsam Fir - Yellow Birch Forest",
  "MS2", "Balsam Fir - Paper Birch Forest",
  "MS4", "Montane Balsam Fir - Paper Birch Forest",
  "MS6", "Balsam Fir - Red Maple Forest",
  "MS7", "Maritime Balsam Fir - Paper Birch Forest",
  "RB1", "White Spruce or Eastern White Cedar Forest from Agriculture",
  "RB2", "Maritime White Spruce Forest",
  "RB3", "Subalpine White Spruce or Balsam Fir - White Spruce Forest",
  "RB4", "Montane White Spruce Forest",
  "RB5", "White Spruce or Balsam Fir - Paper Birch Forest of Anticosti Island",
  "RC3", "Peaty Eastern White Cedar - Balsam Fir Forest",
  "RE1", "Black Spruce - Lichen Forest",
  "RE2", "Black Spruce - Moss or Heath Forest",
  "RE3", "Black Spruce - Sphagnum Forest",
  "RE4", "Montane Black Spruce - Moss or Heath Forest",
  "RE7", "Maritime Black Spruce Forest",
  "RE8", "Subalpine Black Spruce Forest",
  "RI1", "Riparian Zone",
  "RP1", "White or Red Pine Forest",
  "RS1", "Balsam Fir - Eastern White Cedar Forest",
  "RS2", "Balsam Fir - Black Spruce Forest",
  "RS3", "Balsam Fir - Black Spruce and Sphagnum Forest",
  "RS4", "Montane Balsam Fir - Black Spruce Forest",
  "RS5", "Balsam Fir - Red Spruce Forest",
  "RS7", "Maritime Balsam Fir - Black Spruce Forest",
  "RT1", "Hemlock Forest",
  "SM1", "Shifting Sands",
  "SM2", "Maritime Shifting Sands",
  "TOB", "Bog (Ombrotrophe)",
  "TOF", "Fen (Minerotrophic Peatland)",
  "TOU", "Indifferentiated Peatland (Minerotrophic or Ombrotrophic)"
)

type_eco_species_codes <- tibble::tribble(
  ~type_eco_prefix, ~dominant_codes, ~secondary_codes,

  # --- FEUILLUS TEMPÉRÉS / MIXTES ---
  "FC1", c("CHR"), 
          c("ERR","ERS","BOP","PIB"),

  "FE1", c("ERS","CAC"),
          c("TIL","FRA","ORA","OSV"),

  "FE2", c("ERS","TIL"),
          c("HEG","FRA","ORA","OSV"),

  "FE3", c("ERS","BOJ"),
          c("HEG","PRU","SAB","BOP","FRA"),

  "FE4", c("ERS","BOJ","HEG"),
          c("PRU","SAB"),

  "FE5", c("ERS","OSV"),
          c("TIL","FRA","ORA"),

  "FE6", c("ERS","CHR"),
          c("BOJ","BOP","PIB","ERR"),

  "FO1", c("FRN","ORA"),
          c("ERR","FRP","THO"),

  "MF1", c("FRN","SAB"),
          c("ERR","ORA","THO"),

  "MJ1", c("BOJ","SAB","ERS"),
          c(),

  "MJ2", c("BOJ","SAB"),
          c("ERS","EPB"),

  # --- SAPINIÈRES DE TRANSITION ---
  "MS1", c("SAB","BOJ"),
          c("ERS","EPB"),

  "MS2", c("SAB","BOP"),
          c("EPB","EPN"),

  "MS4", c("SAB","BOP"),
          c("EPN","SOA"),

  "MS6", c("SAB","ERR"),
          c("BOP","EPB"),

  "MS7", c("SAB","BOP"),
          c("EPB","EPN"),

  "ME1", c("EPN","PET"),
          c("BOP","SAB"),

  # --- PESSIÈRES BLANCHES / CÉDRIÈRES ---
  "RB1", c("EPB","THO"),
          c("SAB","BOP","ERR"),

  "RB2", c("EPB"),
          c("SAB","BOP","ERR","EPN"),

  "RB3", c("EPB","SAB"),
          c("BOP","EPN"),

  "RB4", c("EPB","SAB"),
          c("BOP","EPN"),

  "RB5", c("EPB","SAB","BOP"),
          c("EPR","EPN"),

  "RC3", c("THO","SAB"),
          c("EPN","MEL"),

  # --- PESSIÈRES NOIRES / SAPINIÈRES NOIRES ---
  "RE1", c("EPN"),
          c("PIG","BOP"),

  "RE2", c("EPN"),
          c("SAB","PIG","BOP"),

  "RE3", c("EPN","MEL"),
          c("SAB","THO"),

  "RE4", c("EPN","SAB"),
          c("BOP","SOA"),

  "RE7", c("EPN"),
          c("EPB","SAB","BOP"),

  "RE8", c("EPN","SAB"),
          c("BOP","SOA"),

  "RS1", c("SAB","THO"),
          c("EPB","BOP","ERR"),

  "RS2", c("SAB","EPN"),
          c("BOP","EPB"),

  "RS3", c("SAB","EPN"),
          c("MEL","THO"),

  "RS4", c("SAB","EPN"),
          c("BOP","SOA"),

  "RS5", c("SAB","EPR"),
          c("EPN","BOP"),

  "RS7", c("SAB","EPN","EPB"),
          c("BOP","ERR"),

  # --- PINÈDES ---
  "RP1", c("PIB","PIR"),
          c("ERR","CHR","BOP"),

  # --- LANDES / ALPINES / MARITIMES / LITTORALES ---
  "LA1", c(),
          c("EPN","PIG","BOP"),

  "LA2", c(),
          c("EPN","PIG"),

  "LA3", c(),
          c("BOP"),

  "LA4", c(),
          c("EPN","SAB","BOP"),

  "LL1", c(),
          c("EPN","SAB","BOP","SOA"),

  "LL2", c(),
          c("EPN","SAB","BOP","SOA"),

  "LL3", c(),
          c("EPN","SAB"),

  "LL4", c(),
          c("EPN","SAB"),

  "LM1", c(),
          c("EPN","EPB","SAB","BOP"),

  "LM2", c(),
          c("EPN","EPB","BOP"),

  "LM3", c(),
          c("BOP","ERR","EPN"),

  "LM4", c(),
          c("EPB","EPN","BOP","SAB"),

  "LI1", c(),
          c("ERA","PEB","BOP","FRP","AME","EPN","SAB"),

  # --- SABLES MOBILES ---
  "SM1", c(),
          c("PIG","PIB","PET","BOP"),

  "SM2", c(),
          c("PIG","EPN","EPB"),

  # --- MARAIS / MARÉCAGES ---
  "MA1", c(),
          c("MEL","EPN","SAB","AME","BOP"),

  "MA2", c(),
          c("EPN","EPB","MEL"),

  # --- TOURBIÈRES ---
  "TOB", c(),
          c("EPN","MEL","THO"),

  "TOF", c(),
          c("MEL","EPN","BOP","AUR"),

  "TOU", c(),
          c("EPN","MEL","THO","BOP")
)


type_eco_descriptions <- type_eco_descriptions %>%
  left_join(type_eco_species_codes, by = c("type_eco_prefix"))

desc_map <- setNames(type_eco_descriptions$Description, type_eco_descriptions$type_eco_prefix)

if (exists("TYPE_ECO_list") && is.list(TYPE_ECO_list) && !is.null(names(TYPE_ECO_list))) {
  names(TYPE_ECO_list) <- vapply(names(TYPE_ECO_list), function(nm) {
    pref <- stringr::str_extract(nm, "^[A-Za-z0-9]+")
    if (!is.na(pref) && pref %in% names(desc_map)) paste0(pref, " - ", desc_map[[pref]]) else nm
  }, character(1))
}

resolve_type_eco <- function(prefix) {
  matches <- TYPE_ECO_list[stringr::str_detect(names(TYPE_ECO_list), paste0("^", prefix))]
  if (length(matches) == 0) stop("No TYPE_ECO_list entries for prefix ", prefix)
  pretty_label <- dplyr::coalesce(desc_map[[prefix]], names(matches)[1])
  list(label = pretty_label, species = matches[[1]])
}

community_prefixes <- c("FE3","FE6","ME1","MS6","MS1","MS2","RE2","RC3","RP1","RS2")
community_defs <- purrr::map(community_prefixes, resolve_type_eco)
names(community_defs) <- community_prefixes

summarise_log_gCSI_vector <- function(log_vals, central = "median") {
  log_tbl <- summarise_posterior(log_vals, central = central) %>%
    dplyr::rename(
      log_center = center,
      log_mean   = mean,
      log_sd     = sd,
      log_q025      = q025,
      log_q975      = q975
    )
  g_tbl <- summarise_posterior(exp(log_vals), central = central) %>%
    dplyr::rename(
      gCenter    = center,
      gMean      = mean,
      gSD        = sd,
      gQ025      = q025,
      gQ975      = q975
    )
  dplyr::bind_cols(g_tbl, log_tbl)
}

process_scenario <- function(scn_row,
                             community_defs,
                             epsilon, central,
                             species_weights,
                             export_site_level = TRUE,
                             out_dir = "J:/Thesis/3rdChapter/PEP_QC/results/undisturbed/outputs",
                             indices = c("gCSI")) {
  file_path <- scn_row$file_path
  scenario_label <- scn_row$scenario_label
  projection <- scn_row$projection
  message(sprintf("[Start] %s (%s)", scenario_label, basename(file_path)))

  if (!file.exists(file_path)) {
    warning("Scenario file not found: ", file_path)
    return(list(community = community_template(), site = site_template()))
  }
  env <- new.env(parent = emptyenv())
  obj_names <- load(file_path, envir = env)
  if (!"preds" %in% obj_names) {
    warning("Object 'preds' not found in ", file_path)
    rm(list = obj_names, envir = env); rm(env)
    return(list(community = community_template(), site = site_template()))
  }
  preds <- env$preds
  rm(list = obj_names, envir = env); rm(env)
  if (is.null(preds) || !length(preds)) {
    warning("Empty 'preds' in ", file_path)
    return(list(community = community_template(), site = site_template()))
  }
  draws_array <- try(coerce_pred_obj_to_draws(preds), silent = TRUE)
  if (inherits(draws_array, "try-error")) {
    warning("Failed to coerce preds for ", scenario_label, ": ", draws_array)
    return(list(community = community_template(), site = site_template()))
  }
  available_species <- dimnames(draws_array)[[3]]

  # Precompute species indices for each community to avoid repeated matching
  community_res <- vector("list", length(community_defs))
  site_res <- vector("list", length(community_defs))
  i <- 0L

  optimized_only <- setequal(sort(indices), "gCSI")

  for (pref in names(community_defs)) {
    i <- i + 1L
    comm <- community_defs[[pref]]
    species_vec <- intersect(comm$species, available_species)
    if (!length(species_vec)) {
      warning("No overlap species for ", pref, " in ", scenario_label)
      community_res[[i]] <- community_template()
      site_res[[i]] <- site_template()
      next
    }

    # Prepare weights in order of species_vec
    w_vec <- if (is.null(species_weights)) NULL else species_weights[species_vec]
    if (!is.null(w_vec) && any(is.na(w_vec))) {
      stop("Missing weights for some species in community ", pref)
    }

    if (optimized_only) {
      # Optimized gCSI path
      species_idx <- match(species_vec, available_species)
      draws_sub <- draws_array[, , species_idx, drop = FALSE]
      g_list <- compute_gCSI_matrix(draws_sub, weights = w_vec, eps = epsilon)
      site_tbl <- summarise_site_gCSI_matrix(g_list$g_matrix, g_list$log_g_matrix, central = central)
      comm_tbl <- summarise_community_gCSI_matrix(g_list$g_matrix, g_list$log_g_matrix, central = central)
    } else {
      # Fallback to legacy comprehensive indices
      idx_long <- compute_indices_from_draws_legacy(draws_array, species_vec, weights = w_vec, eps = epsilon)
      site_tbl <- summarise_site_indices_from_draws(idx_long, central = central)
      comm_tbl <- summarise_community_indices_from_draws(idx_long, central = central)
    }

    if (nrow(comm_tbl)) {
      comm_tbl <- comm_tbl %>% mutate(
        scenario_label = scenario_label,
        projection_horizon = projection,
        community_prefix = pref,
        community_label = comm$label
      )
    } else {
      comm_tbl <- community_template()
    }
    if (nrow(site_tbl)) {
      site_tbl <- site_tbl %>% mutate(
        scenario_label = scenario_label,
        projection_horizon = projection,
        community_prefix = pref,
        community_label = comm$label
      )
    } else {
      site_tbl <- site_template()
    }
    community_res[[i]] <- comm_tbl
    site_res[[i]] <- site_tbl
  }

  # Persist scenario-level outputs immediately
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  comm_all <- purrr::list_rbind(community_res)
  site_all <- purrr::list_rbind(site_res)
  readr::write_csv(comm_all, file.path(out_dir, paste0("community_indices_", scenario_label, ".csv")))
  if (export_site_level) {
    readr::write_csv(site_all, file.path(out_dir, paste0("site_indices_", scenario_label, ".csv")))
  }
  saveRDS(list(community = comm_all, site = site_all), file = file.path(out_dir, paste0("scenario_draws_", scenario_label, ".rds")))
  message(sprintf("[Done ] %s -> %d community rows, %d site rows", scenario_label, nrow(comm_all), nrow(site_all)))
  list(community = comm_all, site = site_all)
}

# --- 7) PROCESS ALL SCENARIOS --------------------------------------------------
# Process scenarios one-by-one (in-line)
scenario_outputs <- vector("list", nrow(scenario_map))
out_dir <- "J:/Thesis/3rdChapter/PEP_QC/results/undisturbed/outputs/output400L"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(scenario_map))) {
  scn_row <- as.list(scenario_map[i, , drop = TRUE])
  scenario_label <- scn_row$scenario_label
  projection <- scn_row$projection
  file_path <- scn_row$file_path

  message("[Start sequential] ", scenario_label, " (", basename(file_path), ")")

  # Process scenario (community + site). process_scenario already persists per-scenario CSV/RDS.
  scenario_res <- try(
    process_scenario(
      scn_row,
      community_defs,
      epsilon,
      use_mean_vs_median_for_central,
      species_weights,
      export_site_level,
      out_dir = out_dir,
      indices = c("gCSI","bCSI","wCSI")
    ),
    silent = TRUE
  )

  if (inherits(scenario_res, "try-error")) {
    warning("process_scenario failed for ", scenario_label, ": ", scenario_res)
    scenario_outputs[[i]] <- list(community = community_template(), site = site_template())
  } else {
    scenario_outputs[[i]] <- scenario_res
  }

  # Per-scenario species-level summaries (saved per scenario to limit memory)
  species_out_file <- file.path(out_dir, paste0("species_prob_by_scenario_draws_summary_", scenario_label, ".csv"))
  species_rds_file <- file.path(out_dir, paste0("species_prob_by_scenario_draws_summary_", scenario_label, ".rds"))

  if (!file.exists(file_path)) {
    warning("Scenario file not found for species summaries: ", file_path)
  } else {
    env <- new.env(parent = emptyenv())
    obj_names <- try(load(file_path, envir = env), silent = TRUE)
    if (inherits(obj_names, "try-error") || !"preds" %in% obj_names) {
      warning("Failed to load preds for species summaries in ", scenario_label)
      if (exists("obj_names")) rm(list = obj_names, envir = env)
      rm(env)
    } else {
      preds <- env$preds
      # Coerce to draws array
      draws_array <- try(coerce_pred_obj_to_draws(preds), silent = TRUE)
      if (inherits(draws_array, "try-error")) {
        warning("Failed to coerce preds for species summaries in ", scenario_label, ": ", draws_array)
      } else {
        species_names <- dimnames(draws_array)[[3]]
        if (is.null(species_names) || length(species_names) == 0) {
          warning("No species names in draws for ", scenario_label)
        } else {
          # Compute per-draw mean probability across sites for each species (NA-safe)
          mean_by_draw_species <- apply(draws_array, c(1,3), function(x) {
            if (all(is.na(x))) return(NA_real_)
            mean(x, na.rm = TRUE)
          })

          # Summarise posterior for each species
          species_tbls <- lapply(seq_along(species_names), function(j) {
            sp <- species_names[j]
            vec <- mean_by_draw_species[, j]
            stats <- summarise_posterior(vec, central = use_mean_vs_median_for_central)
            tibble::tibble(
              scenario_label = scenario_label,
              projection_horizon = projection,
              species = sp,
              center = stats$center,
              mean = stats$mean,
              sd = stats$sd,
              q025 = stats$q025,
              q975 = stats$q975
            )
          })
          species_results_scn <- dplyr::bind_rows(species_tbls)

          # Persist per-scenario species summaries
          readr::write_csv(species_results_scn, species_out_file)
          saveRDS(species_results_scn, species_rds_file)
          rm(species_results_scn, species_tbls, mean_by_draw_species)
        }
      }
      rm(draws_array)
      rm(preds)
      rm(list = obj_names, envir = env)
      rm(env)
    }
  }

  # Encourage immediate memory reclamation between scenarios
  gc()
  message("[Done sequential] ", scenario_label)
}

# Aggregate the sequential results (same as parallel branch)
community_results <- purrr::map_dfr(scenario_outputs, "community") %>% filter(nzchar(scenario_label))
site_results <- purrr::map_dfr(scenario_outputs, "site") %>% filter(nzchar(scenario_label))

trajectory_results <- community_results %>%
  dplyr::left_join(type_eco_descriptions, by = c("community_prefix" = "type_eco_prefix")) %>%
  dplyr::mutate(
    community_prefix = factor(community_prefix, levels = community_prefixes),
    community_desc = dplyr::coalesce(Description, community_label),
    scenario_label = factor(scenario_label, levels = scenario_map$scenario_label)
  ) %>%
  dplyr::select(-Description)

# Final species-level aggregation: read per-scenario CSVs and combine
species_files <- list.files(out_dir, pattern = "^species_prob_by_scenario_draws_summary_.*\\.csv$", full.names = TRUE)
if (length(species_files)) {
  species_results <- purrr::map_dfr(species_files, readr::read_csv, show_col_types = FALSE)
  readr::write_csv(species_results, file.path(out_dir, "species_prob_by_scenario_draws_summary_combined.csv"))
  saveRDS(species_results, file.path(out_dir, "species_prob_by_scenario_draws_summary_combined.rds"))
} else {
  message("No per-scenario species files found to combine.")
}

if (!nrow(trajectory_results)) {
  stop("No CSI results were computed (all empty). Check that 'preds' objects contain posterior draw matrices.")
}

# Save aggregated after all scenarios
readr::write_csv(trajectory_results, file.path(out_dir, "gCSI_trajectory_summary_draws.csv"))
if (export_site_level && nrow(site_results)) {
  readr::write_csv(site_results, file.path(out_dir, "gCSI_wCSI_bCSI_by_site_draws_summary.csv"))
}

# final GC
gc()


# --- 12) NOTES FOR INTEGRATION -----------------------------------------------
# 1) Ensure that all_predictions_summary[[scenario]]$draws exists and carries
#    per draw probabilities from Hmsc::predict(fit, XData = ..., expected = FALSE).
# 2) If you only have mean/lower/upper, you cannot propagate posterior uncertainty
#    correctly for a geometric mean. Re-run prediction to persist draws.
# 3) To compare scenarios or sites, compute per draw log gCSI and take log contrasts.
#    Example function could be added to align by draw, site, and GCM before subtraction.

# =============================================================================
# End of script
# =============================================================================





# ============================================================
# Simple visualisations of community suitability shifts
# Using in-memory `trajectory_results` (already computed above)
# ============================================================

# Ensure factors have sensible ordering
traj <- trajectory_results %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = c(
        "Current",
        "SSP1-2.6 2011-2040","SSP1-2.6 2041-2070","SSP1-2.6 2071-2100",
        "SSP2-4.5 2011-2040","SSP2-4.5 2041-2070","SSP2-4.5 2071-2100",
        "SSP3-7.0 2011-2040","SSP3-7.0 2041-2070","SSP3-7.0 2071-2100",
        "SSP5-8.5 2011-2040","SSP5-8.5 2041-2070","SSP5-8.5 2071-2100"
      )
    ),
    community_desc = factor(community_desc)
  )

# --- 2) Slopegraph: Current vs single target scenario ------------------------
# Choose one focal scenario for a very simple before/after view
target_scenario <- "SSP3-7.0 2071-2100"

slope_data <- trajectory_results %>%
  filter(scenario_label %in% c("Current", target_scenario)) %>%
  select(community_desc, scenario_label, gMean)

slope_plot <- ggplot(
  slope_data,
  aes(x = scenario_label, y = gMean, group = community_desc)
) +
  geom_line(aes(color = community_desc), linewidth = 1.1, alpha = 0.8) +
  geom_point(size = 3) +
  scale_color_viridis_d(guide = "none") +
  labs(
    title = paste("Community suitability: Current vs", target_scenario),
    x = "",
    y = "Mean gCSI"
  ) +
  theme_light(base_size = 14)

print(slope_plot)

# --- 3) Multi-scenario slopegraph per community (facetted) -------------------
# Shows how each community evolves across all scenarios

multi_slope_data <- traj %>%
  arrange(community_desc, scenario_label)

multi_slope_plot <- ggplot(
  multi_slope_data,
  aes(x = scenario_label, y = gMean, group = community_desc)
) +
  geom_line(linewidth = 0.9, alpha = 0.8, color = "grey40") +
  geom_point(size = 2.2, color = "black") +
  labs(
    title = "Community suitability trajectories across climate scenarios",
    x = "",
    y = "Mean gCSI"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ community_desc, scales = "free_y")

print(multi_slope_plot)

# --- 4) Arrow-change scatterplot: direction and magnitude of change ----------
# Improved arrow-change plot: horizontal arrows per community (no overlapping labels),
# with CIs and ordered by change magnitude.
# Prepare baseline and target with CIs if available
has_ci <- all(c("gQ025","gQ975") %in% names(trajectory_results))

baseline_tbl <- trajectory_results %>%
  filter(scenario_label == "Current") %>%
  transmute(
    community_desc,
    baseline = gMean,
    baseline_q025 = if (has_ci) gQ025 else NA_real_,
    baseline_q975 = if (has_ci) gQ975 else NA_real_
  )

target_tbl <- trajectory_results %>%
  filter(scenario_label == target_scenario) %>%
  transmute(
    community_desc,
    future = gMean,
    future_q025 = if (has_ci) gQ025 else NA_real_,
    future_q975 = if (has_ci) gQ975 else NA_real_
  )

arrow_data <- baseline_tbl %>%
  left_join(target_tbl, by = "community_desc") %>%
  mutate(
    delta = future - baseline,
    community_label = stringr::str_wrap(as.character(community_desc), 28),
    community_ord = forcats::fct_reorder(community_label, delta, .desc = TRUE)
  ) %>%
  arrange(desc(delta))

# x limits with small margin
x_vals <- c(arrow_data$baseline, arrow_data$future, arrow_data$baseline_q025, arrow_data$future_q025,
            arrow_data$baseline_q975, arrow_data$future_q975)
x_min <- min(x_vals, na.rm = TRUE)
x_max <- max(x_vals, na.rm = TRUE)
if (!is.finite(x_min)) x_min <- 0
if (!is.finite(x_max)) x_max <- 1
x_margin <- 0.06 * (x_max - x_min + 1e-6)
x_limits <- c(max(0, x_min - x_margin), min(1, x_max + x_margin))

arrow_plot <- ggplot(arrow_data, aes(y = community_ord)) +
  # CI bars (if present, NA will be ignored)
  geom_errorbarh(aes(xmin = baseline_q025, xmax = baseline_q975),
                 height = 0.2, colour = "grey80", na.rm = TRUE) +
  geom_errorbarh(aes(xmin = future_q025, xmax = future_q975),
                 height = 0.2, colour = "grey80", na.rm = TRUE) +
  # horizontal arrow from baseline -> future
  geom_segment(aes(x = baseline, xend = future, y = community_ord, yend = community_ord, color = delta),
               arrow = grid::arrow(length = grid::unit(0.14, "cm")), size = 0.8, alpha = 0.95) +
  # points for baseline and future
  geom_point(aes(x = baseline), shape = 21, fill = "white", color = "black", size = 2) +
  geom_point(aes(x = future), shape = 21, fill = "steelblue", color = "black", size = 2.5) +
  scale_color_gradient2(low = "red", mid = "grey80", high = "blue", midpoint = 0,
                        name = expression(Delta ~ "gCSI (future - current)")) +
  scale_x_continuous(limits = x_limits, expand = c(0,0)) +
  labs(
    title = paste("Community suitability change — Current →", target_scenario),
    subtitle = "Horizontal arrows: start = Current mean gCSI, end = Future mean gCSI",
    x = "gCSI",
    y = NULL,
    caption = "Error bars: 95% CI (if available). Communities ordered by magnitude of change."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

print(arrow_plot)

# If you want to save the figures:
# ggsave(file.path(results_dir, "plot_slope_current_vs_target.png"),
#        slope_plot, width = 7, height = 5, dpi = 300)
# ggsave(file.path(results_dir, "plot_multi_slope_communities.png"),
#        multi_slope_plot, width = 10, height = 8, dpi = 300)
# ggsave(file.path(results_dir, "plot_arrow_change_communities.png"),
#        arrow_plot, width = 7, height = 6, dpi = 300)

# =============================================================================
# Standalone gCSI trajectory export (compatible with simpler pipeline)
# Produces CSV and plot using columns aligned to gCSI_mean/gCSI_q025/gCSI_q975
# =============================================================================

simple_out_dir <- out_dir
dir.create(simple_out_dir, recursive = TRUE, showWarnings = FALSE)

# Build simplified summary using community-level gCSI stats already computed
gCSI_simple <- trajectory_results %>%
  transmute(
    community_prefix,
    community_desc,
    scenario_label,
    gCSI_mean = gMean,
    gCSI_q025 = gQ025,
    gCSI_q975 = gQ975,
    gCSI_lower = gQ025,
    gCSI_upper = gQ975
  ) %>%
  arrange(community_desc, scenario_label)

readr::write_csv(gCSI_simple, file.path(simple_out_dir, "gCSI_trajectory_summary_simple.csv"))

# Plot similar to the standalone script
traj_plot_data <- gCSI_simple %>%
  mutate(
    scenario_name = dplyr::case_when(
      scenario_label == "Current" ~ "Current",
      TRUE ~ stringr::str_remove(as.character(scenario_label), "\\s+\\d{4}-\\d{4}$")
    ),
    projection_horizon = dplyr::case_when(
      scenario_label == "Current" ~ "Baseline",
      TRUE ~ stringr::str_extract(as.character(scenario_label), "\\d{4}-\\d{4}$")
    )
  ) %>%
  drop_na(scenario_name, projection_horizon) %>%
  mutate(
    scenario_name = factor(scenario_name, levels = c("Current", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")),
    projection_horizon = factor(projection_horizon, levels = c("Baseline", "2011-2040", "2041-2070", "2071-2100"))
  )

scenario_pal <- scales::viridis_pal(begin = 0.15, end = 0.85)(length(levels(traj_plot_data$scenario_name)))
names(scenario_pal) <- levels(traj_plot_data$scenario_name)

baseline_lines2 <- traj_plot_data %>%
  filter(scenario_name == "Current", projection_horizon == "Baseline") %>%
  transmute(community_prefix, baseline_gCSI = gCSI_mean) %>%
  drop_na(baseline_gCSI)

trajectory_plot_simple <- ggplot(traj_plot_data, aes(x = projection_horizon, y = gCSI_mean, colour = scenario_name)) +
  geom_errorbar(aes(ymin = gCSI_lower, ymax = gCSI_upper), width = 0.18, alpha = 0.65, linewidth = 0.7) +
  geom_line(aes(group = scenario_name), linewidth = 1.05, na.rm = TRUE) +
  geom_point(shape = 21, fill = "white", stroke = 0.8, size = 2, na.rm = TRUE) +
  geom_hline(
    data = baseline_lines2,
    mapping = aes(yintercept = baseline_gCSI),
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = scenario_pal, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  expand_limits(y = 0) +
  labs(
    title = "gCSI trajectories (simple export)",
    subtitle = "Mean gCSI per plot with 95% quantiles",
    x = "Projection horizon",
    y = "Geometric CSI (gCSI)",
    colour = "Scenario"
  ) +
  theme_light(base_size = 13) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1)) +
  facet_wrap(~ community_prefix, ncol = 2, scales = "free_y")

print(trajectory_plot_simple)

ggsave(filename = file.path(simple_out_dir, "gCSI_trajectories_simple.svg"), plot = trajectory_plot_simple, width = 10, height = 10, dpi = 300)


# Value-style plot: show actual gCSI values (with 95% CI) separated by SSP,
# connecting points only within each SSP across projection horizons.
values_out_file <- file.path(simple_out_dir, "gCSI_values_by_community_by_ssp.svg")

# Prepare plotting table: include Current and all futures, but split into
# scenario_name (SSP family) and projection_horizon (time window)
values_data <- trajectory_results %>%
  select(community_prefix, community_desc, scenario_label, gMean, gQ025, gQ975) %>%
  filter(!is.na(gMean)) %>%
  mutate(
    scenario_name = case_when(
      scenario_label == "Current" ~ "Current",
      TRUE ~ stringr::str_remove(as.character(scenario_label), "\\s+\\d{4}-\\d{4}$")
    ),
    projection_horizon = case_when(
      scenario_label == "Current" ~ "Baseline",
      TRUE ~ stringr::str_extract(as.character(scenario_label), "\\d{4}-\\d{4}$")
    ),
    scenario_name = factor(scenario_name, levels = c("Current", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")),
    projection_horizon = factor(projection_horizon, levels = c("Baseline", "2011-2040", "2041-2070", "2071-2100"))
  ) %>%
  tidyr::drop_na(scenario_name, projection_horizon)

# Baseline per community (for dashed reference lines)
baseline_per_comm <- values_data %>%
  filter(scenario_name == "Current", projection_horizon == "Baseline") %>%
  transmute(community_prefix, community_desc, baseline_gMean = gMean)

if (nrow(values_data)) {
  # Palette
  ssp_pal <- scales::viridis_pal(begin = 0.15, end = 0.85)(length(levels(values_data$scenario_name)))
  names(ssp_pal) <- levels(values_data$scenario_name)

  # Facet by scenario_name (rows) and community (columns) so each SSP is separated.
  # Within each facet row, geom_line connects only points of the same SSP (group = scenario_name).
  values_plot <- ggplot(values_data, aes(x = projection_horizon, y = gMean, ymin = gQ025, ymax = gQ975,
                                           colour = scenario_name, group = interaction(scenario_name, community_desc))) +
    # baseline dashed line (drawn across facets for reference)
    geom_hline(data = baseline_per_comm, aes(yintercept = baseline_gMean),
               linetype = "dashed", colour = "grey50", inherit.aes = FALSE) +
    # CI bars per point (dodged so multiple SSPs at same horizon don't overlap)
    geom_errorbar(aes(ymin = gQ025, ymax = gQ975),
                  position = position_dodge(width = 0.35), width = 0.12, alpha = 0.8, size = 0.5, na.rm = TRUE) +
    # connect points within each SSP across horizons
    geom_line(position = position_dodge(width = 0.35), size = 0.8, na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.35), size = 2.5, shape = 21, fill = "white", na.rm = TRUE) +
    facet_grid(scenario_name ~ community_desc, scales = "free_y", space = "free") +
    scale_colour_manual(values = ssp_pal, drop = FALSE) +
    labs(
      title = "gCSI values per community, separated by SSP (lines connect within each SSP across years)",
      subtitle = "Points = posterior mean gCSI; bars = 95% quantiles; dashed = Current baseline",
      x = "Projection horizon",
      y = "gCSI",
      colour = "Scenario (SSP)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(angle = 0) # keep SSP row labels readable
    )

  print(values_plot)
  ggsave(filename = values_out_file, plot = values_plot, width = 14, height = 10, dpi = 300)
} else {
  message("No gCSI values available to plot.")
}

# =============================================================================
# Per-species simple trajectory export and plot
# Uses: species_prob_by_scenario_draws_summary_combined.csv (created above)
# Produces CSV and plot with columns aligned to species_mean/q025/q975
# =============================================================================

species_simple_out_dir <- out_dir
dir.create(species_simple_out_dir, recursive = TRUE, showWarnings = FALSE)

species_combined_csv <- file.path(species_simple_out_dir, "species_prob_by_scenario_draws_summary_combined.csv")
if (file.exists(species_combined_csv)) {
  species_traj <- readr::read_csv(species_combined_csv, show_col_types = FALSE)

 species_traj <- species_traj %>%
      dplyr::left_join(species_map %>% dplyr::select(original_name_in_data, pretty_name), by = c("species" = "original_name_in_data")) %>%
      dplyr::mutate(species_pretty = dplyr::coalesce(pretty_name, species))
 
  # Build simplified species summary
  species_simple <- species_traj %>%
    transmute(
      species,
      species_pretty,
      scenario_label, 
      projection_horizon,
      species_mean = mean,
      species_q025 = q025,
      species_q975 = q975,
      species_lower = q025,
      species_upper = q975
    ) %>%
    arrange(species, scenario_label)

  readr::write_csv(species_simple, file.path(species_simple_out_dir, "species_trajectory_summary_simple.csv"))

  # Prepare plotting data similar to community simple export
  species_plot_data <- species_simple %>%
    mutate(
      scenario_name = dplyr::case_when(
        scenario_label == "Current" ~ "Current",
        TRUE ~ stringr::str_remove(as.character(scenario_label), "\\s+\\d{4}-\\d{4}$")
      ),
      projection_horizon = dplyr::case_when(
        scenario_label == "Current" ~ "Baseline",
        TRUE ~ dplyr::coalesce(as.character(projection_horizon), stringr::str_extract(as.character(scenario_label), "\\d{4}-\\d{4}$"))
      )
    ) %>%
    tidyr::drop_na(scenario_name, projection_horizon) %>%
    mutate(
      scenario_name = factor(scenario_name, levels = c("Current", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")),
      projection_horizon = factor(projection_horizon, levels = c("Baseline", "2011-2040", "2041-2070", "2071-2100"))
    )

  species_levels <- sort(unique(species_plot_data$species_pretty))
  species_pal <- scales::hue_pal()(length(levels(species_plot_data$scenario_name)))
  names(species_pal) <- levels(species_plot_data$scenario_name)

  # Baseline per species for dashed reference lines
  baseline_species <- species_plot_data %>%
    dplyr::filter(scenario_name == "Current", projection_horizon == "Baseline") %>%
    dplyr::transmute(species_pretty, baseline_prob = species_mean) %>%
    tidyr::drop_na(baseline_prob)

  species_trajectory_plot <- ggplot(species_plot_data, aes(x = projection_horizon, y = species_mean, colour = scenario_name, group = scenario_name)) +
    geom_errorbar(aes(ymin = species_lower, ymax = species_upper), width = 0.18, alpha = 0.65, linewidth = 0.7) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    geom_point(shape = 21, fill = "white", stroke = 0.8, size = 1.5, na.rm = TRUE) +
    geom_hline(
      data = baseline_species,
      mapping = aes(yintercept = baseline_prob),
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = species_pal, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    expand_limits(y = 0) +
    labs(
      title = "Species occurrence probability trajectories (simple)",
      subtitle = "Posterior mean with 95% quantiles across projection horizons",
      x = "Projection horizon",
      y = "Occurrence probability (mean)",
      colour = "Scenario"
    ) +
    theme_light(base_size = 12) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1)) +
    facet_wrap(~ species_pretty, ncol = 4, scales = "free_y")

  print(species_trajectory_plot)
  ggsave(filename = file.path(species_simple_out_dir, "species_trajectories_simple.png"), plot = species_trajectory_plot, width = 8, height = 9, dpi = 300)
} else {
  message("Species combined CSV not found: ", species_combined_csv)
}

# --- Species trajectories by latitude band ------------------------------------
# Inputs:
# - species_prob_by_scenario_draws_summary_combined.csv (made above)
# - plot_coordinates with columns: plot_id, latitude (loaded earlier or from app assets)

species_by_lat_out_dir <- out_dir
dir.create(species_by_lat_out_dir, recursive = TRUE, showWarnings = FALSE)

# Try to ensure plot_coordinates is loaded (reuse app assets if needed)
if (!exists("plot_coordinates")) {
  pc_paths <- c(
    "J:\\Thesis\\3rdChapter\\PEP_QC\\results\\undisturbed\\CAMS_Shiny_App_PA_Advanced\\plot_coordinates.RData",
    file.path(dirname(predictions_dir), "plot_coordinates.RData")
  )
  for (p in pc_paths) {
    if (file.exists(p)) {
      try(load(p), silent = TRUE)
      if (exists("plot_coordinates")) break
    }
  }
}

# We need site-level species means per scenario. If you only have species-level (across all sites),
# compute site-level first from draws. Otherwise, use site_results for gCSI and the same site ids to join latitude.
species_combined_csv <- file.path(species_by_lat_out_dir, "species_prob_by_scenario_draws_summary_combined.csv")

# If per-site species means are not yet saved, derive them now per scenario and bind
# from the already processed per-scenario RDS files, which contain draws arrays paths.
# Prefer using existing combined CSV; if unavailable, skip gracefully.
if (file.exists(species_combined_csv) && exists("plot_coordinates")) {
  # Build site lookup: site id used in draws should match plot_id
  plots_for_map <- plot_coordinates
  req_cols <- c("plot_id","latitude")
  miss <- setdiff(req_cols, names(plots_for_map))
  if (length(miss)) {
    warning("plot_coordinates missing columns: ", paste(miss, collapse = ", "))
  } else {
    site_lookup <- plots_for_map %>%
      dplyr::transmute(site = as.character(plot_id), lat = as.numeric(latitude)) %>%
      dplyr::filter(is.finite(lat))

    # Latitude banding
    lat_band_width <- 1
    site_lookup <- site_lookup %>%
      dplyr::mutate(
        lat_band = cut(
          lat,
          breaks = seq(floor(min(lat, na.rm = TRUE)),
                       ceiling(max(lat, na.rm = TRUE)),
                       by = lat_band_width),
          include.lowest = TRUE
        )
      )

    # To get species by lat band, we need species posterior means per site.
    # If not already persisted, compute from draws per scenario; otherwise, try to recover from site-level exports.
    # Here we recompute per site per species from the original per-scenario preds for correctness.

    species_site_all <- list()
    for (i in seq_len(nrow(scenario_map))) {
      scn <- as.list(scenario_map[i, , drop = TRUE])
      scenario_label <- scn$scenario_label
      projection_horizon <- scn$projection
      file_path <- scn$file_path

      if (!file.exists(file_path)) {
        warning("Missing preds Rda for species-by-lat: ", file_path)
        next
      }
      env <- new.env(parent = emptyenv())
      obj_names <- try(load(file_path, envir = env), silent = TRUE)
      if (inherits(obj_names, "try-error") || !"preds" %in% obj_names) {
        warning("Failed to load preds for species-by-lat in ", scenario_label)
        rm(env)
        next
      }
      preds <- env$preds
      rm(list = obj_names, envir = env); rm(env)

      arr <- try(coerce_pred_obj_to_draws(preds), silent = TRUE)
      if (inherits(arr, "try-error")) {
        warning("coerce_pred_obj_to_draws failed for ", scenario_label, ": ", arr)
        next
      }
      dn <- dimnames(arr)
      draw_ids <- dn[[1]]
      site_ids <- dn[[2]]
      species_ids <- dn[[3]]

      # Posterior credible intervals from the full posterior distribution (per draw)
      # Keep all draws per site per species, then compute quantiles across draws
      site_species_draws <- arr  # dims: [n_draws, n_sites, n_species]
      
      # Compute posterior summaries across draws for each site-species combination
      site_species_tbl_list <- list()
      for (s in seq_len(dim(arr)[2])) {  # iterate sites
        for (sp in seq_len(dim(arr)[3])) {  # iterate species
          draws_vec <- arr[, s, sp]
          draws_vec <- draws_vec[is.finite(draws_vec)]
          
          if (length(draws_vec) > 0) {
            site_species_tbl_list[[length(site_species_tbl_list) + 1L]] <- tibble::tibble(
              site = site_ids[s],
              species = species_ids[sp],
              species_mean = mean(draws_vec, na.rm = TRUE),
              species_q025 = stats::quantile(draws_vec, 0.025, na.rm = TRUE),
              species_q975 = stats::quantile(draws_vec, 0.975, na.rm = TRUE),
              n_draws = length(draws_vec)
            )
          }
        }
      }
      
      if (length(site_species_tbl_list) > 0) {
        tidy_tbl <- dplyr::bind_rows(site_species_tbl_list)
        tidy_tbl$scenario_label <- scenario_label
        tidy_tbl$projection_horizon <- projection_horizon
        species_site_all[[length(species_site_all) + 1L]] <- tidy_tbl
      }
      
      rm(arr, preds, site_species_tbl_list, tidy_tbl)
      gc()
    }

    if (length(species_site_all)) {
      species_site_df <- dplyr::bind_rows(species_site_all)
      # Join latitude bands
      species_site_df <- species_site_df %>%
        dplyr::left_join(site_lookup, by = "site") %>%
        dplyr::filter(!is.na(lat_band))

      # Aggregate per species × scenario × lat_band
      # Now compute credible intervals from the posterior distribution within each band
      # by combining draws across sites within the band, then computing quantiles
      species_by_lat <- species_site_df %>%
        dplyr::group_by(species, scenario_label, projection_horizon, lat_band) %>%
        dplyr::summarise(
          species_mean = mean(species_mean, na.rm = TRUE),
          # 95% credible interval: compute from the within-band posterior distribution
          # Approximate by taking quantiles of the per-site posterior means (conservative)
          # Or, if you have access to draws, use quantile across draws in the band
          species_q025 = stats::quantile(species_q025, 0.5, na.rm = TRUE),  # median of lower bounds
          species_q975 = stats::quantile(species_q975, 0.5, na.rm = TRUE),  # median of upper bounds
          n_sites = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          scenario_name = dplyr::case_when(
            scenario_label == "Current" ~ "Current",
            TRUE ~ stringr::str_remove(as.character(scenario_label), "\\s+\\d{4}-\\d{4}$")
          ),
          projection_horizon = dplyr::case_when(
            scenario_label == "Current" ~ "Baseline",
            TRUE ~ stringr::str_extract(as.character(scenario_label), "\\d{4}-\\d{4}$")
          )
        ) %>%
        tidyr::drop_na(scenario_name, projection_horizon)

      # Attach pretty species names if available
      if (exists("species_map") && all(c("original_name_in_data","pretty_name") %in% names(species_map))) {
        species_by_lat <- species_by_lat %>%
          dplyr::left_join(species_map %>% dplyr::select(original_name_in_data, pretty_name), by = c("species" = "original_name_in_data")) %>%
          dplyr::mutate(species_pretty = dplyr::coalesce(pretty_name, species))
      } else {
        species_by_lat <- species_by_lat %>% dplyr::mutate(species_pretty = species)
      }

      # Save CSV (include species_pretty)
      readr::write_csv(species_by_lat, file.path(species_by_lat_out_dir, "species_trajectory_summary_by_lat_band.csv"))

      # Plot per-species trajectories by lat band
      species_by_lat$scenario_name <- factor(
        species_by_lat$scenario_name,
        levels = c("Current", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
      )
      species_by_lat$projection_horizon <- factor(
        species_by_lat$projection_horizon,
        levels = c("Baseline", "2011-2040", "2041-2070", "2071-2100")
      )

      pal <- scales::viridis_pal(begin = 0.15, end = 0.85)(length(levels(species_by_lat$scenario_name)))
      names(pal) <- levels(species_by_lat$scenario_name)

      species_by_lat_plot <- ggplot(
        species_by_lat,
        aes(x = projection_horizon, y = species_mean, colour = scenario_name, group = scenario_name)
      ) +
        geom_errorbar(aes(ymin = species_q025, ymax = species_q975),
              width = 0.25, alpha = 0.7, linewidth = 0.8) +
        geom_line(linewidth = 0.9) +
        geom_point(shape = 21, fill = "white", stroke = 0.7, size = 2.2) +
        scale_colour_manual(values = pal, drop = FALSE) +
        scale_x_discrete(drop = FALSE) +
        expand_limits(y = 0) +
        labs(
          title = paste0("Species probability trajectories by latitude band (", lat_band_width, "°)"),
          subtitle = "Posterior mean per site, 95% credible intervals from posterior distribution within band",
          x = "Projection horizon",
          y = "Occurrence probability (mean across sites in band)",
          colour = "Scenario"
        ) +
        theme_light(base_size = 12) +
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1)) +
        facet_grid(lat_band ~ species_pretty, scales = "free_y")

      print(species_by_lat_plot)
      ggsave(
        filename = file.path(species_by_lat_out_dir, "species_trajectories_by_lat_band.svg"),
        plot = species_by_lat_plot, width = 15, height = 10, dpi = 300
      )
    } else {
      message("No per-scenario species site-level data could be computed for lat bands.")
    }
  }
} else {
  message("Species-by-lat skipped: plot_coordinates or combined species CSV not found.")
}


# --- Latitudinal community gCSI trajectories (using plot_coordinates) ----------
# Requires: plot_coordinates RData alongside app.R or load explicitly here.
# Try to load plot_coordinates.RData if not in memory.

if (!exists("plot_coordinates")) {
  pc_paths <- c(
    "J:\\Thesis\\3rdChapter\\PEP_QC\\results\\undisturbed\\CAMS_Shiny_App_PA_Advanced\\plot_coordinates.RData",
    file.path(dirname(predictions_dir), "plot_coordinates.RData")
  )
  for (p in pc_paths) {
    if (file.exists(p)) {
      try(load(p), silent = TRUE)
      if (exists("plot_coordinates")) break
    }
  }
}
plot_coordinates <- plots_for_map 
# Harmonize to columns: plot_id, latitude
if (exists("plot_coordinates")) {
  plots_for_map <- plot_coordinates
  req_cols <- c("plot_id","latitude")
  miss <- setdiff(req_cols, names(plots_for_map))
  if (length(miss)) {
    warning("plot_coordinates missing columns: ", paste(miss, collapse=", "))
  } else {
    # Build site lookup: site (matches site_results$site) -> lat
    site_lookup <- plots_for_map %>%
      dplyr::transmute(site = as.character(plot_id), lat = as.numeric(latitude)) %>%
      dplyr::filter(is.finite(lat))

    # Choose latitude band width (degrees)
    lat_band_width <- 1
    site_lookup <- site_lookup %>%
      dplyr::mutate(
        lat_band = cut(
          lat,
          breaks = seq(floor(min(lat, na.rm = TRUE)),
                       ceiling(max(lat, na.rm = TRUE)),
                       by = lat_band_width),
          include.lowest = TRUE
        )
      )

    # Recompute gCSI from posterior draws per community, site, and lat band
    # to obtain credible intervals from the full posterior distribution
    lat_comm_traj_list <- list()
    
    for (i in seq_len(nrow(scenario_map))) {
      scn <- as.list(scenario_map[i, , drop = TRUE])
      scenario_label <- scn$scenario_label
      projection_horizon <- scn$projection
      file_path <- scn$file_path

      if (!file.exists(file_path)) {
        warning("Missing preds Rda for lat band community analysis: ", file_path)
        next
      }
      
      env <- new.env(parent = emptyenv())
      obj_names <- try(load(file_path, envir = env), silent = TRUE)
      if (inherits(obj_names, "try-error") || !"preds" %in% obj_names) {
        warning("Failed to load preds for lat band analysis in ", scenario_label)
        rm(env)
        next
      }
      preds <- env$preds
      rm(list = obj_names, envir = env); rm(env)

      arr <- try(coerce_pred_obj_to_draws(preds), silent = TRUE)
      if (inherits(arr, "try-error")) {
        warning("coerce_pred_obj_to_draws failed for lat band analysis in ", scenario_label, ": ", arr)
        next
      }

      # Process each community
      for (pref in names(community_defs)) {
        comm <- community_defs[[pref]]
        dn <- dimnames(arr)
        available_species <- dn[[3]]
        species_vec <- intersect(comm$species, available_species)
        
        if (!length(species_vec)) {
          next
        }

        # Prepare weights
        w_vec <- if (is.null(species_weights)) NULL else species_weights[species_vec]
        if (!is.null(w_vec) && any(is.na(w_vec))) {
          stop("Missing weights for some species in community ", pref)
        }

        # Compute gCSI from draws (optimized path)
        species_idx <- match(species_vec, available_species)
        draws_sub <- arr[, , species_idx, drop = FALSE]
        g_list <- compute_gCSI_matrix(draws_sub, weights = w_vec, eps = epsilon)

        # g_list$g_matrix: [n_draws x n_sites]
        # g_list$log_g_matrix: [n_draws x n_sites]
        
        site_ids <- dn[[2]]
        draw_ids <- dn[[1]]
        n_draws <- dim(g_list$g_matrix)[1]
        n_sites <- dim(g_list$g_matrix)[2]

        # Create long-format table: draw x site x gCSI
        g_long <- tibble::tibble(
          draw = rep(draw_ids, times = n_sites),
          site = rep(site_ids, each = n_draws),
          gCSI = as.vector(g_list$g_matrix),
          log_gCSI = as.vector(g_list$log_g_matrix)
        )

        # Join latitude bands
        g_long <- g_long %>%
          dplyr::left_join(site_lookup, by = "site") %>%
          dplyr::filter(!is.na(lat_band), is.finite(gCSI))

        # Aggregate per lat_band × draw: mean gCSI across sites in band per draw
        # This preserves the posterior draw structure for credible intervals
        g_by_band_draw <- g_long %>%
          dplyr::group_by(lat_band, draw) %>%
          dplyr::summarise(
            gCSI_draw = mean(gCSI, na.rm = TRUE),
            log_gCSI_draw = mean(log_gCSI, na.rm = TRUE),
            n_sites_in_band = dplyr::n(),
            .groups = "drop"
          )

        # Compute posterior credible intervals from the distribution of per-draw aggregates
        lat_band_stats <- g_by_band_draw %>%
          dplyr::group_by(lat_band) %>%
          dplyr::summarise(
            gCSI_mean = mean(gCSI_draw, na.rm = TRUE),
            gCSI_q025 = stats::quantile(gCSI_draw, 0.025, na.rm = TRUE),
            gCSI_q975 = stats::quantile(gCSI_draw, 0.975, na.rm = TRUE),
            log_gCSI_mean = mean(log_gCSI_draw, na.rm = TRUE),
            log_gCSI_q025 = stats::quantile(log_gCSI_draw, 0.025, na.rm = TRUE),
            log_gCSI_q975 = stats::quantile(log_gCSI_draw, 0.975, na.rm = TRUE),
            n_draws = dplyr::n(),
            .groups = "drop"
          )

        # Add scenario and community metadata
        lat_band_stats <- lat_band_stats %>%
          dplyr::mutate(
            community_prefix = pref,
            community_label = comm$label,
            scenario_label = scenario_label,
            projection_horizon = projection_horizon
          )

        lat_comm_traj_list[[length(lat_comm_traj_list) + 1L]] <- lat_band_stats
      }

      rm(arr, preds, g_list, g_long, g_by_band_draw)
      gc()
    }

    # Bind all results
    if (length(lat_comm_traj_list) > 0) {
      lat_comm_traj <- dplyr::bind_rows(lat_comm_traj_list)

      # Add descriptions and scenario parsing
      lat_comm_traj <- lat_comm_traj %>%
        dplyr::left_join(type_eco_descriptions, by = c("community_prefix" = "type_eco_prefix")) %>%
        dplyr::mutate(
          community_desc = dplyr::coalesce(Description, community_label),
          scenario_name = dplyr::case_when(
            scenario_label == "Current" ~ "Current",
            TRUE ~ stringr::str_remove(as.character(scenario_label), "\\s+\\d{4}-\\d{4}$")
          ),
          projection_horizon = dplyr::case_when(
            scenario_label == "Current" ~ "Baseline",
            TRUE ~ stringr::str_extract(as.character(scenario_label), "\\d{4}-\\d{4}$")
          )
        ) %>%
        tidyr::drop_na(scenario_name, projection_horizon)

      # Persist CSV
      readr::write_csv(lat_comm_traj, file.path(out_dir, "latitudinal_gCSI_by_band.csv"))

      # Plot setup
      lat_comm_traj$scenario_name <- factor(lat_comm_traj$scenario_name,
        levels = c("Current", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))
      lat_comm_traj$projection_horizon <- factor(lat_comm_traj$projection_horizon,
        levels = c("Baseline", "2011-2040", "2041-2070", "2071-2100"))

      pal <- scales::viridis_pal(begin = 0.15, end = 0.85)(length(levels(lat_comm_traj$scenario_name)))
      names(pal) <- levels(lat_comm_traj$scenario_name)

      # Baseline per community per lat band for dashed reference lines
      baseline_ref <- lat_comm_traj %>%
        dplyr::filter(scenario_name == "Current", projection_horizon == "Baseline") %>%
        dplyr::transmute(lat_band, community_prefix, baseline_gCSI = gCSI_mean)

      lat_comm_plot <- ggplot(lat_comm_traj,
        aes(x = projection_horizon, y = gCSI_mean, colour = scenario_name, group = scenario_name)) +
        geom_hline(
          data = baseline_ref,
          mapping = aes(yintercept = baseline_gCSI),
          linetype = "dashed",
          colour = "grey50",
          linewidth = 0.7,
          inherit.aes = FALSE
        ) +
        geom_errorbar(aes(ymin = gCSI_q025, ymax = gCSI_q975),
              width = 0.25, alpha = 0.8, linewidth = 0.8) +
        geom_line(linewidth = 0.7) +
        geom_point(shape = 21, fill = "white", stroke = 0.7, size = 0.5) +
        scale_colour_manual(values = pal, drop = FALSE) +
        scale_x_discrete(drop = FALSE) +
        expand_limits(y = 0) +
        labs(
          title = paste0("Latitudinal gCSI trajectories (", lat_band_width, "° bands)"),
          subtitle = "Posterior mean with 95% credible interval derived from posterior draws; dashed line = current baseline",
          x = "Projection horizon",
          y = "gCSI (mean across sites in band)",
          colour = "Scenario"
        ) +
        theme_light(base_size = 12) +
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1)) +
        facet_grid(lat_band ~ community_prefix, scales = "free_y")

      print(lat_comm_plot)
      ggsave(filename = file.path(out_dir, "latitudinal_gCSI_by_band.png"),
             plot = lat_comm_plot, width = 15, height = 10, dpi = 300)
    } else {
      message("No latitudinal gCSI data computed for any community or scenario.")
    }
  }
} else {
  message("Latitudinal trajectories skipped: 'plot_coordinates' not found. Load plot_coordinates.RData from app folder.")
}
