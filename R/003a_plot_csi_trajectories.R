# =============================================================================
# 004_gCSI_trajectory_visuals_from_trajectory_results.R
# Purpose : Create interpretable visualisations using the exported file
#           trajectory_results.csv with columns:
#           gCSI_mean, gCSI_lower, gCSI_upper, gCSI_q025, gCSI_q975,
#           scenario_label, scenario_key, community_prefix,
#           community_label, Description, community_desc
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# -------------------------------------------------------------------------
# 0) LOAD EXPORTED RESULTS
# -------------------------------------------------------------------------

results_dir <- "J:/Thesis/3rdChapter/PEP_QC/results/undisturbed/outputs"

trajectory_results <- readr::read_csv(
  file.path(results_dir, "gCSI_trajectory_summary_draws.csv"),#trajectory_results.csv
  show_col_types = FALSE
)

# Defensive check
required_cols <- c(
  "gCSI_mean", "gCSI_lower", "gCSI_upper",
  "gCSI_q025", "gCSI_q975",
  "scenario_label",
  "community_desc"
)

missing_cols <- setdiff(required_cols, names(trajectory_results))
if (length(missing_cols) > 0) {
  stop(
    "trajectory_results.csv is missing the following columns:\n",
    paste(missing_cols, collapse = ", "), "\n",
    "Adapt the script to your actual structure."
  )
}
# index	center	mean	sd	q025	q975	log_center	log_mean	log_sd	log_q025	log_q975	gCenter	gMean	gSD	gQ025	gQ975	bCenter	bMean	bSD	bQ025	bQ975	wCenter	wMean	wSD	wQ025	wQ975	scenario_label	projection_horizon	community_prefix	community_label	community_desc
# summary	NA	NA	NA	NA	NA	-3.022383143	-3.002351471	0.035082285	-3.022806248	-2.964869773	0.048685057	0.04969065	0.001760546	0.048664462	0.051571593	0.3235046139181269	0.3235164437678573	3.13699856350303e-4	0.323223616	0.3238193273476415	0.3235046139181269	0.3235164437678573	3.13699856350303e-4	0.323223616	0.3238193273476415	Current	Baseline	FE3	Sugar Maple - Yellow Birch Forest	Sugar Maple - Yellow Birch Forest
# summary	NA	NA	NA	NA	NA	-3.357546063	-3.369732357	0.032532243	-3.404145727	-3.345677337	0.034820602	0.034410916	0.00111049	0.033237127	0.03523647	0.25727163843800116	0.2578414271990559	0.001435122	0.25680334786054365	0.2593638269844647	0.25727163843800116	0.2578414271990559	0.001435122	0.25680334786054365	0.2593638269844647	Current	Baseline	FE6	Sugar Maple - Red Oak Forest	Sugar Maple - Red Oak Forest

# -------------------------------------------------------------------------
# 1) RENAME TO GENERIC NAMES USED BY THE PLOTTING LOGIC
# -------------------------------------------------------------------------

trajectory_results <- trajectory_results %>%
  dplyr::rename(
    center    = gCSI_mean,
    hdi_lower = gCSI_lower,
    hdi_upper = gCSI_upper,
    q025      = gCSI_q025,
    q975      = gCSI_q975
  )

# If projection_horizon is not present, reconstruct it from scenario_label
if (!"projection_horizon" %in% names(trajectory_results)) {
  trajectory_results <- trajectory_results %>%
    dplyr::mutate(
      projection_horizon = dplyr::case_when(
        scenario_label == "Current" ~ "Baseline",
        TRUE ~ stringr::str_extract(scenario_label, "\\d{4}-\\d{4}")
      )
    )
}

# -------------------------------------------------------------------------
# 2) PREP DATA FOR PLOTTING
# -------------------------------------------------------------------------

trajectory_plot_data <- trajectory_results %>%
  dplyr::mutate(
    scenario_name = dplyr::case_when(
      scenario_label == "Current" ~ "Current",
      TRUE ~ stringr::str_extract(as.character(scenario_label),
                                  "^SSP[0-9]-[0-9]\\.[0-9]")
    ),
    projection_horizon = dplyr::coalesce(projection_horizon, "Baseline")
  ) %>%
  tidyr::drop_na(scenario_name, projection_horizon) %>%
  dplyr::mutate(
    scenario_name = factor(
      scenario_name,
      levels = c("Current", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    ),
    projection_horizon = factor(
      projection_horizon,
      levels = c("Baseline", "2011-2040", "2041-2070", "2071-2100")
    )
  )

# Palette for scenarios
scenario_palette <- scales::viridis_pal(begin = 0.15, end = 0.85)(
  length(levels(trajectory_plot_data$scenario_name))
)
names(scenario_palette) <- levels(trajectory_plot_data$scenario_name)

hdi_prob <- 0.89  # used when you computed gCSI_lower/gCSI_upper

# -------------------------------------------------------------------------
# 3) RELATIVE TRAJECTORIES: % CHANGE VS CURRENT BASELINE
# -------------------------------------------------------------------------

baseline_tbl <- trajectory_plot_data %>%
  dplyr::filter(
    scenario_name == "Current",
    projection_horizon == "Baseline"
  ) %>%
  dplyr::select(
    community_desc,
    baseline_center = center
  )

delta_plot_data <- trajectory_plot_data %>%
  dplyr::left_join(baseline_tbl, by = "community_desc") %>%
  dplyr::mutate(
    delta_pct       = 100 * (center    - baseline_center) / baseline_center,
    delta_hdi_lower = 100 * (hdi_lower - baseline_center) / baseline_center,
    delta_hdi_upper = 100 * (hdi_upper - baseline_center) / baseline_center
  )

delta_trajectory_plot <- ggplot(
  delta_plot_data,
  aes(x = projection_horizon, y = delta_pct, colour = scenario_name)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(
    aes(ymin = delta_hdi_lower, ymax = delta_hdi_upper),
    width = 0.18, alpha = 0.7, linewidth = 0.7
  ) +
  geom_line(aes(group = scenario_name), linewidth = 1.05, na.rm = TRUE) +
  geom_point(shape = 21, fill = "white", stroke = 0.8, size = 3.0, na.rm = TRUE) +
  scale_colour_manual(values = scenario_palette, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title    = "Relative gCSI trajectories",
    subtitle = "Percent change vs current baseline (mean with HDI)",
    x        = "Projection horizon",
    y        = "Δ gCSI (%) relative to current baseline",
    colour   = "Scenario"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 25, hjust = 1)
  ) +
  facet_wrap(~ community_desc, ncol = 2, scales = "free_y")

print(delta_trajectory_plot)

ggplot2::ggsave(
  filename = file.path(results_dir, "gCSI_delta_trajectories.png"),
  plot   = delta_trajectory_plot,
  width  = 10,
  height = 10,
  dpi    = 300
)

# -------------------------------------------------------------------------
# 4) END-OF-CENTURY COMPARISON (2071–2100) – FOREST PLOT STYLE
# -------------------------------------------------------------------------

end_century_data <- trajectory_plot_data %>%
  dplyr::filter(projection_horizon == "2071-2100")

end_century_plot <- ggplot(
  end_century_data,
  aes(x = scenario_name, y = center, ymin = hdi_lower, ymax = hdi_upper)
) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_pointrange(position = position_dodge(width = 0.4), size = 0.4) +
  scale_colour_manual(values = scenario_palette, drop = FALSE) +
  labs(
    title    = "End-of-century gCSI by scenario (2071–2100)",
    subtitle = "Mean gCSI with credible intervals per tree community",
    x        = "Scenario",
    y        = "gCSI (mean with interval)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  ) +
  facet_wrap(~ community_desc, ncol = 2, scales = "free_y")

print(end_century_plot)

ggplot2::ggsave(
  filename = file.path(results_dir, "gCSI_end_century_forest_plot.png"),
  plot   = end_century_plot,
  width  = 10,
  height = 10,
  dpi    = 300
)

# -------------------------------------------------------------------------
# 5) HEATMAP OF MEAN gCSI (DESCRIPTIVE OVERVIEW)
# -------------------------------------------------------------------------

heatmap_data <- trajectory_plot_data %>%
  dplyr::mutate(
    scenario_compact = dplyr::case_when(
      scenario_label == "Current" ~ "Current – Baseline",
      TRUE ~ paste0(
        stringr::str_extract(as.character(scenario_label),
                             "^SSP[0-9]-[0-9]\\.[0-9]"),
        " ",
        as.character(projection_horizon)
      )
    ),
    scenario_compact = factor(
      scenario_compact,
      levels = unique(scenario_compact)
    )
  )

gCSI_heatmap <- ggplot(
  heatmap_data,
  aes(x = community_desc, y = scenario_compact, fill = center)
) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", center)), size = 2.8) +
  scale_fill_viridis_c(option = "D", name = "Mean gCSI") +
  labs(
    title = "Overview of mean gCSI across scenarios and communities",
    x     = "Tree community",
    y     = "Scenario and projection horizon"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1)
  )

print(gCSI_heatmap)

ggplot2::ggsave(
  filename = file.path(results_dir, "gCSI_heatmap_mean.png"),
  plot   = gCSI_heatmap,
  width  = 11,
  height = 7,
  dpi    = 300
)

# =============================================================================
# End – visualisations based on trajectory_results.csv
# =============================================================================


# =============================================================================
# Visualising gCSI Trajectories and Climate Impacts (summary-based)
# Uses: trajectory_results.csv
# Produces:
#   1) Heatmap of proportional change in gCSI
#   2) Pareto ("resilience") scatter plot
#   3) Alluvial (“Sankey-style”) plot of gCSI categories
#   4) Bump chart of ranking shifts
#   5) Volcano-like plot: effect size vs uncertainty
# =============================================================================

# --- 0) Libraries -------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggridges)     # not strictly required here, but handy if you extend
  library(ggalluvial)
  library(viridis)      # for viridis palettes
})

# --- 1) Load trajectory summary file -----------------------------------------
# Adjust path if needed
traj_path <- "F:/Thesis/3rdChapter/PEP_QC/results/undisturbed/outputs/trajectory_results.csv"

traj <- readr::read_csv(traj_path, show_col_types = FALSE)

# Expected columns (from your head()):
# gCSI_mean, gCSI_lower, gCSI_upper, gCSI_q025, gCSI_q975,
# scenario_label, scenario_key, community_prefix, community_label,
# Description, community_desc

# --- 2) Basic transforms: scenario ordering, log scale, deltas ----------------
eps <- 1e-9

# Order scenarios: "Current" first, others in their existing order
scenario_levels <- traj %>%
  distinct(scenario_label) %>%
  pull()
scenario_levels <- c("Current", setdiff(scenario_levels, "Current"))

traj <- traj %>%
  mutate(
    scenario_label = factor(scenario_label, levels = scenario_levels),
    log_center = log(gCSI_mean + eps)
  )

# Baseline per community (Current)
baseline <- traj %>%
  filter(scenario_label == "Current") %>%
  select(community_desc,
         gCSI_mean_baseline = gCSI_mean,
         log_baseline       = log_center)

# Attach baseline and compute deltas for other scenarios
traj_delta <- traj %>%
  left_join(baseline, by = "community_desc") %>%
  mutate(
    delta_raw   = gCSI_mean - gCSI_mean_baseline,
    delta_log   = log_center - log_baseline,
    pct_change  = exp(delta_log) - 1,               # multiplicative change
    ci_width    = gCSI_upper - gCSI_lower           # crude uncertainty proxy
  )

# -----------------------------------------------------------------------------#
# 1) Heatmap of proportional change (ΔgCSI on log scale → exp(Δ) − 1)
# -----------------------------------------------------------------------------#
heat_data <- traj_delta %>%
  filter(scenario_label != "Current")

p_heat <- ggplot(heat_data, aes(
  x     = scenario_label,
  y     = community_desc,
  fill  = pct_change
)) +
  geom_tile(color = "grey30") +
  scale_fill_viridis(option = "B", direction = 1,
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Proportional change in gCSI relative to current climate",
    x        = "Scenario",
    y        = "Community",
    fill     = "ΔgCSI (%)\n(exp(Δlog) − 1)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid  = element_blank()
  )

print(p_heat)

# -----------------------------------------------------------------------------#
# 2) Pareto "resilience" plot:
#    x = current gCSI, y = proportional change under each scenario
# -----------------------------------------------------------------------------#
pareto_data <- traj_delta %>%
  filter(scenario_label != "Current")

p_pareto <- ggplot(pareto_data, aes(
  x     = gCSI_mean_baseline,
  y     = pct_change,
  colour = scenario_label,
  label = community_desc
)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_colour_viridis_d(option = "D") +
  labs(
    title = "Community resilience: baseline suitability vs climate-induced change",
    x     = "Baseline gCSI (Current climate)",
    y     = "Proportional change in gCSI\n(exp(Δlog) − 1)",
    colour = "Scenario"
  ) +
  theme_minimal(base_size = 13)

print(p_pareto)

# -----------------------------------------------------------------------------#
# 3) Alluvial (“Sankey-style”) plot of gCSI categories across scenarios
#    Categories based on gCSI_mean:
#      Low      : gCSI < 0.33
#      Moderate : 0.33–0.66
#      High     : > 0.66
# -----------------------------------------------------------------------------#
# install.packages("ggalluvial")  # if not already installed

traj_cat <- traj %>%
  mutate(
    scenario_label = factor(scenario_label, levels = scenario_levels),
    gCSI_cat = cut(
      gCSI_mean,
      breaks = c(-Inf, 0.33, 0.66, Inf),
      labels = c("Low", "Moderate", "High")
    )
  )

stopifnot(ggalluvial::is_alluvia_form(
  traj_cat,
  axes = 1,
  silent = TRUE
))

p_alluvial <- ggplot(
  traj_cat,
  aes(
    x           = scenario_label,
    stratum     = gCSI_cat,
    alluvium    = community_desc,
    fill        = gCSI_cat,
    label       = gCSI_cat
  )
) +
  ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "frontback",
                        alpha = 0.6) +
  ggalluvial::geom_stratum(color = "grey20") +
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Transitions in gCSI categories across climate scenarios",
    x     = "Scenario",
    y     = "Number of communities",
    fill  = "gCSI category"
  ) +
  theme_minimal(base_size = 13)

print(p_alluvial)

# -----------------------------------------------------------------------------#
# 4) Bump chart: ranking shifts of communities under each scenario
# -----------------------------------------------------------------------------#
# Bump chart: panel per SSP (exclude Current) and legend for communities
bump_data <- traj %>%
  dplyr::mutate(
    scenario_name = dplyr::case_when(
      scenario_label == "Current" ~ "Current",
      TRUE ~ stringr::str_extract(as.character(scenario_label),
                                  "^SSP[0-9]-[0-9]\\.[0-9]")
    ),
    projection_horizon = dplyr::coalesce(
      stringr::str_extract(as.character(scenario_label), "\\d{4}-\\d{4}"),
      "Baseline"
    ),
    projection_horizon = factor(
      projection_horizon,
      levels = c("Baseline", "2011-2040", "2041-2070", "2071-2100")
    ),
    scenario_name = factor(
      scenario_name,
      levels = c("Current", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    )
  ) %>%
  # keep only SSP panels (drop Current) — remove filter() if you want Current too
  dplyr::filter(scenario_name != "Current") %>%
  dplyr::group_by(scenario_name, projection_horizon) %>%
  dplyr::mutate(rank = rank(-gCSI_mean, ties.method = "first")) %>%
  dplyr::ungroup()

max_rank <- max(bump_data$rank, na.rm = TRUE)

# Widescreen bump chart with colorblind-friendly viridis palette and stronger lines
p_bump <- ggplot(bump_data, aes(
  x      = projection_horizon,
  y      = rank,
  group  = community_desc,
  colour = community_desc
)) +
  geom_line(size = 1.1, alpha = 0.95, na.rm = TRUE) +
  geom_point(size = 2.8, shape = 21, fill = "white", stroke = 0.6, na.rm = TRUE) +
  facet_wrap(~ scenario_name, ncol = 2, scales = "free_x") +
  scale_y_reverse(breaks = seq_len(max_rank)) +
  # viridis is perceptually uniform and colorblind-friendly; widen the gamut for distinct colors
  scale_colour_viridis_d(option = "D", begin = 0.05, end = 0.95,
                         guide = guide_legend(ncol = 3, byrow = TRUE,
                                              override.aes = list(size = 3))) +
  labs(
    title = "Ranking shifts in community suitability across projection horizons",
    subtitle = "Panels per SSP; lines coloured by community",
    x = "Projection horizon",
    y = "Rank (1 = highest gCSI)",
    colour = "Community"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

print(p_bump)

# Save widescreen output (adjust path/filenames as desired)
if (exists("results_dir")) {
  ggplot2::ggsave(
    filename = file.path(results_dir, "gCSI_bumpchart_widescreen.png"),
    plot   = p_bump,
    width  = 14,   # widescreen
    height = 10,
    dpi    = 300
  )
}

# -----------------------------------------------------------------------------#
# 5) Volcano-like plot: effect size (Δlog gCSI) vs uncertainty (CI width)
# -----------------------------------------------------------------------------#
volc_data <- traj_delta %>%
  filter(scenario_label != "Current")

p_volcano <- ggplot(volc_data, aes(
  x      = delta_log,
  y      = ci_width,
  colour = community_desc,
  label  = scenario_label
)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 100) +
  scale_colour_viridis_d(option = "E", guide = "none") +
  labs(
    title = "Climate effect vs uncertainty on gCSI",
    x     = "Effect size (Δ log gCSI)",
    y     = "Uncertainty (gCSI_upper − gCSI_lower)"
  ) +
  theme_minimal(base_size = 13)

print(p_volcano)

# =============================================================================
# End of script
# =============================================================================
