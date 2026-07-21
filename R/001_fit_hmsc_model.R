# ---
# title: "Joint Species Distribution Models in R for Species Abundance"
# author: "Liz, JPC."
# date: "2025-08-12"
# ---

#===============================================================================
# 0. SETUP AND PACKAGE LOADING
#===============================================================================

# for the actual analysis
# install.packages("devtools")
# devtools::install_github("hmsc-r/hmsc", force = TRUE)
# install.packages("sf")
# install.packages("webshot2")
# install.packages("ggtext")
# install.packages("ggpubr")
# install.packages("kableExtra")
# install.packages("knitr")
# for parallel processing
.libPaths(c("H:/R_lib", .libPaths()))
library(tidyverse)
library(Hmsc)
library(coda)
library(readr)
library(sf)
library(doParallel)
library(foreach)
# for aesthetics
library(ggtext)
library(knitr)
library(kableExtra)
library(ggpubr)
# library(conflicted)
library(BayesLogit)
# ===============================================================================
# 1. TABLE TO MAP SPECIES (DICTIONARY)
# ===============================================================================
print("--- 0. Criando a tabela de mapeamento de espécies ---")

# Create a dataframe 
species_map <- tibble::tribble(
  ~code, ~french_name, ~scientific_name,
  "AEH", "Marronnier d'Inde", "Aesculus hippocastanum",
  "AME", "Amélanchier", "Amelanchier spp.",
  "AUC", "Aulne crispé", "Alnus crispa",
  "AUR", "Aulne rugueux", "Alnus incana",
  "BOG", "Bouleau à feuilles de peuplier (gris)", "Betula populifolia",
  "BOJ", "Bouleau jaune", "Betula alleghaniensis",
  "BOP", "Bouleau à papier (blanc)", "Betula papyrifera",
  "CAC", "Caryer cordiforme", "Carya cordiformis",
  "CAF", "Caryer ovale (à noix douce)", "Carya ovata",
  "CAR", "Charme de Caroline", "Carpinus caroliniana",
  "CET", "Cerisier tardif", "Prunus serotina",
  "CHB", "Chêne blanc", "Quercus alba",
  "CHE", "Chêne bicolore", "Quercus bicolor",
  "CHG", "Chêne à gros fruits", "Quercus macrocarpa",
  "CHR", "Chêne rouge", "Quercus rubra",
  "EPB", "Épinette blanche", "Picea glauca",
  "EPN", "Épinette noire", "Picea mariana",
  "EPO", "Épinette de Norvège", "Picea abies",
  "EPR", "Épinette rouge", "Picea rubens",
  "ERA", "Érable argenté", "Acer saccharinum",
  "ERB", "Érable de Norvège", "Acer platanoides",
  "ERE", "Érable à épis", "Acer spicatum",
  "ERG", "Érable à Giguère", "Acer negundo",
  "ERN", "Érable noir", "Acer nigrum",
  "ERP", "Érable de Pennsylvanie", "Acer pensylvanicum",
  "ERR", "Érable rouge", "Acer rubrum",
  "ERS", "Érable à sucre", "Acer saccharum",
  "FRA", "Frêne d'Amérique (blanc)", "Fraxinus americana",
  "FRN", "Frêne noir", "Fraxinus nigra",
  "FRP", "Frêne de Pennsylvanie (rouge)", "Fraxinus pennsylvanica",
  "HEG", "Hêtre à grandes feuilles", "Fagus grandifolia",
  "MEL", "Mélèze laricin", "Larix laricina",
  "NOC", "Noyer cendré", "Juglans cinerea",
  "ORA", "Orme d'Amérique", "Ulmus americana",
  "ORR", "Orme rouge", "Ulmus rubra",
  "ORT", "Orme de Thomas", "Ulmus thomasii",
  "OSV", "Ostryer de Virginie", "Ostrya virginiana",
  "PEB", "Peuplier baumier", "Populus balsamifera",
  "PEG", "Peuplier à grandes dents", "Populus grandidentata",
  "PET", "Peuplier faux-tremble", "Populus tremuloides",
  "PIB", "Pin blanc", "Pinus strobus",
  "PIG", "Pin gris", "Pinus banksiana",
  "PIR", "Pin rouge", "Pinus resinosa",
  "PIS", "Pin sylvestre", "Pinus sylvestris",
  "PRU", "Pruche de l'Est", "Tsuga canadensis",
  "PRP", "Cerisier de Pennsylvanie", "Prunus pensylvanica",
  "SAB", "Sapin baumier", "Abies balsamea",
  "SOA", "Sorbier d'Amérique", "Sorbus americana",
  "SOD", "Sorbier des montagnes", "Sorbus decora",
  "THO", "Thuya occidental", "Thuja occidentalis",
  "TIL", "Tilleul d'Amérique", "Tilia americana"
)

# Create acronym (ex: Abi.bals)
create_pretty_name <- function(scientific_name) {
  parts <- strsplit(scientific_name, " ")[[1]]
  genus <- substr(parts[1], 1, 3)
  species <- if (length(parts) > 1) substr(parts[2], 1, 4) else ""
  # Capitaliza o gênero e junta com a espécie
  paste0(str_to_title(genus), ".", species)
}

# apply 
species_map <- species_map %>%
  mutate(
    original_name_in_data = paste0(code, "_n"),
    pretty_name = sapply(scientific_name, create_pretty_name)
  )

print("Table species_map created:")
print(head(species_map))


#===============================================================================
# 2. DATA LOADING AND PREPARATION
#===============================================================================

# --- Define Paths ---
output_dir <- "/Thesis/3rdChapter/PEP_QC/results/undisturbed"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# plots_sf_env_clim_tot <- read_csv("/Thesis/3rdChapter/PEP_QC/databases/plots_sf_env_clim_4th_inv.csv")
plots_sf_env_clim <- read_csv("/Thesis/3rdChapter/PEP_QC/databases/undisturbed_plots_sf_env_clim_4th_inv.csv")
# --- Define Predictor and Species Variables ---
# Environmental predictors to be considered
# plots_sf_env_clim %>% group_by() %>%
#   summarise(CEC = sum(is.na(CEC)), 
#             pH = sum(is.na(pH)), 
#             Clay = sum(is.na(Clay)), 
#             OM = sum(is.na(OM)))
colnames(plots_sf_env_clim)
# # number total of plots
# plots_sf_env_clim_tot %>% group_by() %>%
#   summarise(total_plots = n())

# env_predictors <- c("MAP", "MAT", "OM", "slope","drainage_class") 
env_predictors <-  c( "MAP", "MAT", "pH",
                    "Clay", "OM", "rtp", "twi",
                    "slope", "drainage_class")

# ensure drainage_class is a factor
plots_sf_env_clim$drainage_class <- as.factor(plots_sf_env_clim$drainage_class)

#"pH","Clay","rtp","twi","surface_deposit" retire aspect_class

# Plot histograms of species abundance

# plots_sf_env_clim %>%
#   pivot_longer(cols = ends_with("_n"), names_to = "species", values_to = "abundance") %>%
#   dplyr::filter(abundance > 0) %>%
#   ggplot(aes(x = abundance)) +
#   geom_histogram(bins = 30) +
#   facet_wrap(~ species, scales = "free_x") +
#   theme_minimal()


# Species abundance columns (assuming they end with "_n")
species_abun_vars <- plots_sf_env_clim %>%
  select(ends_with("_n") & !starts_with("NA")) %>%
  colnames()

# --- Create the Initial Species Abundance Matrix (Y_full) ---
Y_full <- plots_sf_env_clim %>%
  dplyr::select(plot_id, all_of(species_abun_vars)) %>%
  mutate(across(all_of(species_abun_vars), ~ as.integer(replace_na(., 0)))) %>%
  as.data.frame() %>%
  column_to_rownames(var = "plot_id") %>%
  as.matrix()

# --- ** REMOVE RARE SPECIES** ---
# 1. Calculate prevalence (number of occurrences) for each species
prevalence <- colSums(Y_full > 0, na.rm = TRUE)

# 2. Define the threshold (e.g., 10% of the total number of plots)
n_plots <- nrow(Y_full)
prevalence_threshold <- n_plots * 0.05  # Change 0.1 to 0.05 for 5%

# 3. Identify the species to keep
species_to_keep <- names(prevalence[prevalence >= prevalence_threshold])

# 4. Filter the abundance matrix to keep only the common species
Y <- Y_full[, species_to_keep]

print(paste("Removed", ncol(Y_full) - ncol(Y), "rare species (present in < 5% of plots)."))
print(paste("Proceeding with", ncol(Y), "species."))

# --- Create the Environmental Covariate Data Frame (XData) ---
XData <- plots_sf_env_clim %>%
  select(plot_id, all_of(env_predictors)) %>%
  mutate_if(is.character, as.factor) %>%
  as.data.frame() %>%  # Convert to data.frame to allow row names
  column_to_rownames(var = "plot_id")

str(XData)
#===============================================================================
# 3. DATA ALIGNMENT (CRITICAL STEP)
#===============================================================================
# Now, align all data based on the plots that have complete environmental data
# AND are present in the filtered Y matrix.

# Find common plot IDs between the filtered Y matrix and XData with complete cases
common_plots <- intersect(rownames(Y), rownames(na.omit(XData)))

# Filter all data objects to this final, clean set of plots
Y <- Y[common_plots, ]

XData <- XData[common_plots, ]

# Remove species that are no longer present after filtering plots
Y <- Y[, colSums(Y) > 0]



print(paste("Data aligned. Final dimensions:"))
print(paste("Species Matrix (Y):", nrow(Y), "plots,", ncol(Y), "species"))
print(paste("Covariate Matrix (XData):", nrow(XData), "plots,", ncol(XData), "covariates"))

# Adapt to presence absence if needed

Y[Y > 0] <- 1

#===============================================================================
# 4. HMSC MODEL SETUP
#===============================================================================

# --- Define the model formula ---
XFormula <- as.formula(paste("~", paste(env_predictors, collapse = " + ")))

# --- Define the study design and random effects ---
# The study design links plots to the random effect levels.
studyDesign <- data.frame(site = as.factor(rownames(Y)))
rownames(studyDesign) <- rownames(Y)

# # (A) Pick the lon/lat columns from the central DB (robust to name variants)
# lon_col <- if ("longitude" %in% names(plots_sf_env_clim)) "longitude" else "Longitude"
# lat_col <- if ("latitude"  %in% names(plots_sf_env_clim)) "latitude"  else "Latitude"

# if (!all(c(lon_col, lat_col) %in% names(plots_sf_env_clim))) {
#   stop("Could not find longitude/latitude columns in plots_sf_env_clim.")
# }

# # (B) Build projected coordinates ONCE (meters), then keep them in the same table
# #     Pick a metric CRS suitable for your region. Example for Québec:
# target_crs <- 32198  # NAD83 / Quebec Lambert (change if needed)

# plots_sf_env_clim <- plots_sf_env_clim %>%
#   # make an sf using the lon/lat picked above
#   sf::st_as_sf(coords = c(lon_col, lat_col), crs = 4326, remove = FALSE) %>%
#   # project to meters
#   sf::st_transform(crs = target_crs) %>%
#   # extract projected x,y and drop geometry
#   dplyr::mutate(
#     x = sf::st_coordinates(geometry)[,1],
#     y = sf::st_coordinates(geometry)[,2]
#   ) %>%
#   sf::st_drop_geometry()

# # (C) Quick sanity checks
# stopifnot(is.numeric(plots_sf_env_clim$x), is.numeric(plots_sf_env_clim$y))
# if (anyNA(plots_sf_env_clim$x) || anyNA(plots_sf_env_clim$y)) {
#   stop("Missing projected coordinates (x/y) for some plots.")
# }

# -------------------------------------------------------------
# 2) Build coords aligned to studyDesign (row order must match)
# -------------------------------------------------------------
# coords <- plots_sf_env_clim %>%
#   dplyr::filter(plot_id %in% rownames(studyDesign)) %>%
#   # enforce the exact same order as studyDesign rows
#   dplyr::slice(match(rownames(studyDesign), plot_id)) %>%
#   dplyr::select(x, y) %>%
#   as.matrix()

# rownames(coords) <- rownames(studyDesign)  # keep rownames identical

# # final safety check
# stopifnot(identical(rownames(coords), rownames(studyDesign)))

# -------------------------------------------------------------
# 3) Spatial random level for Hmsc (with a reasonable nfMax)
# -------------------------------------------------------------
# rL_site <- HmscRandomLevel(sData = coords)
# rL_site$nfMax <- 15  # good starting point; tune 10–20 as needed

# Define the random level itself
# rL.site <- HmscRandomLevel(units = unique(studyDesign$site)) # run with no spatial structure

rL_site <- HmscRandomLevel(units = studyDesign$site)

# --- Construct the Hmsc model object ---
# Using "lognormal poisson" for abundance data, as it handles overdispersion.
# or "probit" for presence-absence data.
mod <- Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  distr = "probit", # distribution for presence-absence data 
   #"lognormal poisson", # distribution for abundance if needed
  studyDesign = studyDesign,
  ranLevels = list(site = rL_site)
)

print("Hmsc model object created successfully.")

#===============================================================================
# 5. MCMC SAMPLING SETUP
#===============================================================================

nChains <- 4 # A robust number for checking convergence
thin = 100
samples = 1000

#===============================================================================
# 6. MCMC SAMPLING (THE BUSINESS) - SIMPLIFIED WORKFLOW
#===============================================================================

set.seed(1304)
  ptm = proc.time()

    m1 = sampleMcmc(mod, samples = samples, thin = thin,
    adaptNf = rep(ceiling(0.4*samples*thin),1),
    transient = ceiling(0.5*samples*thin),
    nChains = nChains, nParallel = nChains,
    initPar = "fixed effects")

  computational.time = proc.time() - ptm
print(paste("Computational time:", computational.time))
# Save the model object
save(m1, file = file.path(output_dir, "hmsc", "hmsc_presence_model.Rdata"))

print("Model object saved.")

#===============================================================================
# 7. MODEL CONVERGENCE DIAGNOSTICS
#===============================================================================
print("--- 7. Checking Model Convergence ---")

# Convert the Hmsc object to a coda object for diagnostics
mpost <- convertToCodaObject(m1)

# Calculate Effective Sample Size (ESS) and Potential Scale Reduction Factor (PSRF)
ess.beta <- effectiveSize(mpost$Beta) %>%
  as_tibble() %>% dplyr::rename(ess_beta = value)

psrf.beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf %>%
  as_tibble() %>% dplyr::rename(psrf_beta = `Point est.`)

# Create diagnostic plots
diag_ess_plot <- ggplot(ess.beta, aes(x = ess_beta)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Effective Sample Size (ESS)", x = "ESS", y = "Count") +
  theme_bw()

diag_psrf_plot <- ggplot(psrf.beta, aes(x = psrf_beta)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = 1.1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Gelman-Rubin Diagnostic (PSRF)", x = "PSRF (Point Estimate)", y = "Count") +
  theme_bw()

diag_all <- ggarrange(diag_ess_plot, diag_psrf_plot, ncol = 2)

# Save and print the plot
ggsave(file.path(output_dir, "hmsc_figures", "convergence_diagnostics.png"), diag_all, width = 7, height = 3.5, bg = "white")
print(diag_all)
print("Convergence diagnostics saved. Check for ESS > 1000 and PSRF < 1.1 for most parameters.")

#===============================================================================
# 8. MODEL FIT AND EXPLANATORY POWER
#===============================================================================
print("--- 8. Evaluating Model Fit ---")

# For abundance models (lognormal poisson), explanatory R2 is a good metric.
preds <- computePredictedValues(m1)
MF <- evaluateModelFit(hM = m1, predY = preds)


# Create a prevalence table to join with R2 results
prevalence <- colSums(m1$Y > 0, na.rm = TRUE) %>%
  as_tibble(rownames = "Species") %>%
  dplyr::rename(prevalence = value)

# prevalence <- colSums(Y > 0, na.rm = TRUE) %>%
#   as_tibble(rownames = "Species") %>%
#   dplyr::rename(prevalence = value) %>%
#   arrange(desc(prevalence))

# Combine R2 and prevalence
model_fit_df <- data.frame(
  Species = colnames(m1$Y),
  ExplanatoryTJUR = MF$TjurR2
) %>%
  left_join(prevalence, by = "Species") %>%
  left_join(species_map, by = c("Species" = "original_name_in_data"))


# After MF is calculated and model_fit_df is created:
species_with_zero_sd <- model_fit_df %>%
  dplyr::filter(is.na(ExplanatoryTJUR))

print("Species with zero standard deviation (TJUR is NA):")
print(species_with_zero_sd)

# Check prevalence for these species
prevalence_of_problem_species <- prevalence %>%
  dplyr::filter(Species %in% species_with_zero_sd$Species)
print(prevalence_of_problem_species)

# Plot TJUR vs. Prevalence
tjur_vs_prevalence_plot <- ggplot(model_fit_df, aes(x = prevalence, y = ExplanatoryTJUR)) +
  geom_point(color = "darkblue", size = 3, alpha = 0.7) +
  labs(
    title = "Model Explanatory Power vs. Species Prevalence",
    x = "Prevalence (Number of Plots Occupied)",
    y = expression(paste("Explanatory ", TjurR^2))
  ) +
  geom_text(aes(label = pretty_name), vjust = 1.5, size = 3) +
  theme_bw()

ggsave(file.path(output_dir,"hmsc_figures", "TJUR_vs_prevalence.svg"), tjur_vs_prevalence_plot, 
       width = 6, height = 4, bg = "white")
print(tjur_vs_prevalence_plot)

# Display TJUR table
print("Explanatory TjurR2 per species:")
model_fit_df %>%
  arrange(desc(ExplanatoryTJUR)) %>%
  knitr::kable() %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE)
library("webshot2")
save_kable(
  model_fit_df %>%
    arrange(desc(ExplanatoryTJUR)) %>%
    knitr::kable() %>%
    kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE),
  file.path(output_dir, "hmsc_figures" , "Explanatory_TjurR2_per_species.pdf")
)
# Plot AUC and RMSE vs. Prevalence if needed for presence-absence models
# str(MF)
# List of 3
#  $ RMSE  : num [1:13] 0.314 0.186 0.217 0.34 0.327 ...
#  $ AUC   : num [1:13] 0.918 0.986 0.989 0.929 0.939 ...
#  $ TjurR2: num [1:13] 0.429 0.617 0.519 0.445 0.43 ...
# Create figures folder if missing
dir.create(file.path(output_dir, "hmsc_figures"), recursive = TRUE, showWarnings = FALSE)

# Append AUC and RMSE from MF to the model fit table
model_fit_df <- model_fit_df %>%
  mutate(
    AUC = MF$AUC,
    RMSE = MF$RMSE,
    TjurR2 = MF$TjurR2
  )

# Export modelfit table pdfs of AUC and RMSE vs Prevalence table 
save_kable(
  model_fit_df %>%
    select(Species, pretty_name, prevalence, AUC, RMSE) %>%
    arrange(desc(AUC)) %>%
    knitr::kable() %>%
    kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE),
  file.path("J:/Thesis/3rdChapter/PEP_QC/results/undisturbed/hmsc_figures", "AUC_RMSE_per_Species.pdf")
)



# AUC vs Prevalence plot
auc_vs_prevalence_plot <- ggplot(model_fit_df %>% filter(!is.na(AUC)), aes(x = prevalence, y = AUC)) +
  geom_point(color = "black", size = 2, alpha = 0.7) +
  # geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(label = pretty_name), vjust = 1.5, size = 3) +
  labs(
    title = "AUC vs Species Prevalence",
    x = "Prevalence (Number of Plots Occupied)",
    y = "AUC"
  ) +
  theme_bw()

# RMSE vs Prevalence plot
rmse_vs_prevalence_plot <- ggplot(model_fit_df %>% filter(!is.na(RMSE)), aes(x = prevalence, y = RMSE)) +
  geom_point(color = "black", size = 2, alpha = 0.7) +
  # geom_hline(yintercept = median(model_fit_df$RMSE, na.rm = TRUE), color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(label = pretty_name), vjust = 1.5, size = 3) +
  labs(
    title = "RMSE vs Species Prevalence",
    x = "Prevalence (Number of Plots Occupied)",
    y = "RMSE"
  ) +
  theme_bw()

# Combine and save plots
combined_auc_rmse <- ggarrange(auc_vs_prevalence_plot, rmse_vs_prevalence_plot, ncol = 2, nrow = 1)
ggsave(file.path(output_dir, "hmsc_figures", "AUC_RMSE_vs_prevalence.svg"), combined_auc_rmse, width = 12, height = 4, bg = "white")
ggsave(file.path(output_dir, "hmsc_figures", "AUC_vs_prevalence.png"), auc_vs_prevalence_plot, width = 6, height = 4, bg = "white")
ggsave(file.path(output_dir, "hmsc_figures", "RMSE_vs_prevalence.png"), rmse_vs_prevalence_plot, width = 6, height = 4, bg = "white")

# Save a table of metrics per species
save_kable(
  model_fit_df %>%
    select(Species, pretty_name, prevalence, AUC, RMSE) %>%
    arrange(desc(AUC)),
  file.path(output_dir, "hmsc_figures", "Model_Metrics_per_Species.pdf")
)

# Print plot objects for interactive sessions / logs
print(auc_vs_prevalence_plot)
print(rmse_vs_prevalence_plot)

# Summary table for AUC RMSE and TjurR2 across all species

kable_summary <- model_fit_df %>%
  knitr::kable() %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE)
save_kable(
  kable_summary,
  file.path(output_dir, "hmsc_figures", "Model_Fit_Summary.pdf")
)

#===============================================================================
# 9. VARIANCE PARTITIONING
#===============================================================================
print("--- 9. Partitioning Variance ---")

VP <- computeVariancePartitioning(m1)

# Prepare data for plotting
vp_df <- VP$vals %>%
  as_tibble(rownames = "variable") %>%
  pivot_longer(
    cols = -variable,
    names_to = "Species",
    values_to = "value"
  )

# save vp_df
save_kable(
  vp_df %>%
    knitr::kable() %>%
    kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE),
  file.path(output_dir, "hmsc_figures" , "Variance_Partitioning_per_Species.pdf")
)

# Order species by prevalence for a more organized plot
species_order <- prevalence %>% arrange(prevalence) %>% pull(Species)
vp_df$Species <- factor(vp_df$Species, levels = species_order)

# --- 1. Calculate the Mean Importance of Each Variable ---
# We group by 'variable' and calculate its average contribution across all species.
variable_importance <- vp_df %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  arrange(desc(mean_value)) # Arrange from most to least important

print("Average importance of each variable:")
print(variable_importance)
# save kable
save_kable(
  variable_importance %>%
    knitr::kable() %>%
    kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE),
  file.path(output_dir, "hmsc_figures" , "Variable_Importance_Variance_Partitioning.pdf")
)
# --- 2. Create the Ordered Factor Levels ---
# Get the names of your variables in their ordered sequence
ordered_vars <- (variable_importance$variable)
fixed_effects <- setdiff(ordered_vars, "Random: site")

# Extended Okabe-Ito-like palette to 10 distinct, colorblind-friendly colors.
okabe_ito_palette <- c(
  "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733",
  "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"
)

# Create a named vector to map each variable to a specific color
# The random effect will always be grey.
color_mapping <- c(
  # Assign Okabe-Ito colors to the fixed effects
  setNames(rep_len(okabe_ito_palette, length.out = length(fixed_effects)), fixed_effects),
  # Assign grey to the random effect
  "Random: site" = "grey50"
)

# --- 3. Create the Plot with the Custom Colorblind-Safe Palette ---
# Join vp_df with species_map to get pretty names (select only needed columns to avoid .x/.y suffixes)
vp_df <- vp_df %>%
  left_join(
    species_map %>% select(original_name_in_data, pretty_name),
    by = c("Species" = "original_name_in_data")
  ) %>%
  # coalesce any duplicated pretty_name columns created by the join (.x/.y)
  mutate(
    pretty_name = dplyr::coalesce(
      if ("pretty_name" %in% names(.)) .data$pretty_name else NA_character_,
      if ("pretty_name.x" %in% names(.)) .data$pretty_name.x else NA_character_,
      if ("pretty_name.y" %in% names(.)) .data$pretty_name.y else NA_character_
    )
  ) %>%
  # drop any residual ".x" / ".y" suffixed columns
  select(-matches("\\.x$|\\.y$"))
  
# Order species by prevalence using pretty names
species_order_df <- model_fit_df  %>%
  arrange(prevalence)

vp_df$pretty_name <- factor(vp_df$pretty_name, levels = species_order_df$pretty_name)

vp_plot <- ggplot(vp_df, aes(x = value, y = pretty_name, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Variance Partitioning of Species Abundance",
    #subtitle = "Variables are stacked from most to least important (left to right)",
    x = "Proportion of Explained Variance",
    y = "Species",
    fill = "Predictor Group"
  ) +
  scale_x_continuous(labels = scales::percent) +
  
  # Apply the custom, named color palette
  scale_fill_manual(values = color_mapping) +
  
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8)
  )

# --- 4. Save and Print the New Plot ---
ggsave(
  file.path(file.path(output_dir, "hmsc_figures", "variance_partitioning_ordered_colorblind.svg")),
  vp_plot,
  height = 9,
  width = 7,
  bg = "white"
)

print("Ordered variance partitioning plot has been saved.")
print(vp_plot)

#===============================================================================
# 10. SPECIES NICHES (BETA COEFFICIENTS)
#===============================================================================
print("--- 10. Visualizing Species Niches (Beta Coefficients) using ggplot ---")

# --- 5a. Get Posterior Estimates for Beta Coefficients ---
postBeta <- getPostEstimate(m1, parName = "Beta")

# --- 5b. Prepare the Data for Plotting in a Tidy Format ---

# Automatically get the correct predictor names from the model object
predictor_names <- m1$covNames

# Extract the mean of the posterior distribution for Beta
means_df <- postBeta$mean %>%
  as_tibble() %>%
  # Add predictor names as a column
  mutate(variable = predictor_names, .before = 1) %>%
  # Pivot to a long format
  pivot_longer(
    cols = -variable,
    names_to = "Species",
    values_to = "Mean"
  )

# Extract the support levels (probability of being positive or negative)
supported_df <- postBeta$support %>%
  as_tibble() %>%
  mutate(variable = predictor_names, .before = 1) %>%
  pivot_longer(
    cols = -variable,
    names_to = "Species",
    values_to = "Support"
  ) %>%
  # Join with the mean values
  left_join(means_df, by = c("variable", "Species")) %>%
  # Filter for only the strongly supported effects
  dplyr::filter(
    (Support > 0.95 | Support < 0.05),
    variable != "(Intercept)" # Still exclude the intercept
  ) %>%
  # Create a 'sign' column for coloring the tile borders
  mutate(sign = ifelse(Mean > 0, "+", "-")) %>%
  # Join with species_map to get pretty names
  left_join(species_map, by = c("Species" = "original_name_in_data")) %>%
  # Filter out any without pretty names
  dplyr::filter(!is.na(pretty_name))

print(supported_df, n = 123)

# --- 5c. Order Species and Variables for a Clean Plot ---

# Order species by prevalence using pretty names
species_order <- model_fit_df %>% arrange(desc(prevalence)) %>% pull(pretty_name)
supported_df$pretty_name <- factor(supported_df$pretty_name, levels = species_order)

# Order environmental variables for a logical layout
variable_order <- c(
  "(Intercept)", "MAP", "MAT", "pH",
  "Clay", "OM", "rtp", "twi",
  "slope",
  "drainage_class10", "drainage_class11", "drainage_class16",
  "drainage_class20", "drainage_class21", "drainage_class22",
  "drainage_class30", "drainage_class31", "drainage_class33", "drainage_class34",
  "drainage_class40", "drainage_class41", "drainage_class43",
  "drainage_class50", "drainage_class51", "drainage_class53",
  "drainage_class60"
)

supported_df$variable <- factor(supported_df$variable, levels = variable_order)

# Create meaningful labels for drainage classes based on the system
drainage_labels <- c(
  "(Intercept)" = "(Intercept)",
  "MAP" = "MAP",
  "MAT" = "MAT",
  "pH" = "pH",
  "Clay" = "Clay",
  "OM" = "OM",
  "rtp" = "rtp",
  "twi" = "twi",
  "slope" = "slope",
  "drainage_class10" = "Rapid",
  "drainage_class11" = "Rapid (Lateral drainage)",
  "drainage_class16" = "Complex drainage",
  "drainage_class20" = "Good",
  "drainage_class21" = "Good (Lateral drainage)",
  "drainage_class22" = "Good (Frozen horizon)",
  "drainage_class30" = "Moderate",
  "drainage_class31" = "Moderate (Lateral drainage)",
  "drainage_class33" = "Moderate (Anthropogenic improvement)",
  "drainage_class34" = "Moderate (Anthropogenic slowing)",
  "drainage_class40" = "Imperfect",
  "drainage_class41" = "Imperfect (Lateral drainage)",
  "drainage_class43" = "Imperfect (Anthropogenic improvement)",
  "drainage_class50" = "Bad",
  "drainage_class51" = "Bad (Lateral drainage)",
  "drainage_class53" = "Bad (Anthropogenic improvement)",
  "drainage_class60" = "Very bad"
)

# --- 5d. Create the Final Heatmap Plot ---
p_beta <- ggplot(supported_df, aes(x = variable, y = pretty_name, fill = Mean, color = sign)) +
  geom_tile(linewidth = 0.5) + # Use linewidth instead of lwd
  theme_classic() +
  
  # Use a diverging color scale for the fill (magnitude of the effect)
  scale_fill_gradient2(
    low = "darkblue", mid = "grey90", high = "darkred",
    midpoint = 0, name = "Effect Size (Mean Beta)"
  ) +
  
  # Use a manual color scale for the border (sign of the effect)
  scale_color_manual(values = c("+" = "red", "-" = "blue")) +
  
  # Hide the color guide for the border
  guides(color = "none") +
  
  # Apply the meaningful labels to the x-axis
  scale_x_discrete(labels = drainage_labels) +
  
  labs(
    title = "Significant Species Responses to Environmental Predictors",
    subtitle = "Showing effects with >95% posterior support. Border color indicates sign.",
    x = "Environmental Variable",
    y = "Species"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# --- 5e. Save and Display the Plot ---
ggsave(
  file.path(output_dir,"hmsc_figures", "beta_coefficients_heatmap.svg"),
  p_beta,
  width = 7,
  height = 9,
  bg = "white"
)

print("Custom species niches heatmap has been saved.")
print(p_beta)


#===============================================================================
# 11. GRADIENT VISUALIZATION
#===============================================================================
print("--- 11. Plotting Predictions Along an Environmental Gradient ---")

# Example: Mean Annual Temperature (MAT)
# You can change "MAT" to any other continuous variable in your XFormula
gradient_var <- "slope" # e.g., "MAT", "MAP", "slope", "Clay", "OM", "pH", "rtp", "twi", "drainage_class"

# Construct the gradient
gradient <- constructGradient(m1, focalVariable = gradient_var)

# Predict along the gradient
predY_gradient <- predict(m1, XData = gradient$XDataNew, studyDesign = gradient$studyDesignNew,
              ranLevels = gradient$rLNew, expected = TRUE)

# Select the 13 most prevalent species to plot
species_to_plot <- prevalence %>% arrange(desc(prevalence)) %>% head(13) %>% pull(Species)

# Get corresponding pretty names
species_to_plot_pretty <- species_map %>%
  filter(original_name_in_data %in% species_to_plot) %>%
  arrange(match(original_name_in_data, species_to_plot)) %>%
  pull(pretty_name)

# Plot the gradient for the selected species
# Define a vector of gradient variables to loop over
gradient_vars <- c("MAP", "MAT", "slope", "Clay", "OM", "pH", "rtp", "twi")

# Loop over each gradient variable
for (gradient_var in gradient_vars) {
  # Construct the gradient
  gradient <- constructGradient(m1, focalVariable = gradient_var)
  
  # Predict along the gradient
  predY_gradient <- predict(m1, XData = gradient$XDataNew, studyDesign = gradient$studyDesignNew,
                            ranLevels = gradient$rLNew, expected = TRUE)
  
  # Select the 13 most prevalent species to plot
  species_to_plot <- prevalence %>% arrange(desc(prevalence)) %>% head(13) %>% pull(Species)
  
  # Get corresponding pretty names
  species_to_plot_pretty <- species_map %>%
    filter(original_name_in_data %in% species_to_plot) %>%
    arrange(match(original_name_in_data, species_to_plot)) %>%
    pull(pretty_name)
  
  # Plot the gradient for the selected species
  svg(file.path(output_dir, "hmsc_figures", paste0("gradient_plot_", gradient_var, ".svg")), width = 12, height = 8)
  par(mfrow = c(3, 5)) # Create a grid for plots
  for (i in 1:length(species_to_plot)) {
    sp <- species_to_plot[i]
    pretty_name <- species_to_plot_pretty[i]
    plotGradient(m1, gradient, predY_gradient, yshow = c(0, max(m1$Y[, sp])),
                 measure = "Y", index = which(colnames(m1$Y) == sp), showData = TRUE, q = c(0.025, 0.5, 0.975),
                 xlabel = gradient_var, ylabel = "Predicted Occurrence Probability",
                 main = pretty_name, showPosteriorSupport = FALSE)
  }
  par(mfrow = c(1, 1)) # Reset plotting device
  dev.off()
  print(paste("Gradient plots for", gradient_var, "saved to SVG file with pretty names."))
}



#===============================================================================
#  12. SPECIES CO-OCCURRENCE (using ggcorrplot)
#===============================================================================
print("--- Analyzing Species Co-occurrence Patterns ---")

# Ensure ggcorrplot is installed
if (!requireNamespace("ggcorrplot", quietly = TRUE)) install.packages("ggcorrplot")
library(ggcorrplot)

# --- 12a. Compute the Association Matrix ---
OmegaCor <- computeAssociations(m1)

# --- 12b. Prepare Data for Plotting ---
# Extract the matrix of mean posterior estimates for the correlations (Omega)
# These are the effect sizes we want to visualize.
cor_matrix_mean <- OmegaCor[[1]]$mean

# as tibble for inspection
cor_matrix_mean_df <- as_tibble(cor_matrix_mean, rownames = "Species1") %>%
  pivot_longer(
    cols = -Species1,
    names_to = "Species2",
    values_to = "Correlation"
  ) %>% save_kable(
    file.path(output_dir, "hmsc_figures" , "Species_Correlation_Mean_Values.pdf")
  ) %>% save_kable(
    file.path(output_dir, "hmsc_figures" , "Species_Correlation_Mean_Values.pdf")
  )

# Create a matrix of p-values (or support levels) to identify significant associations.
# Here, we calculate the probability of the association being positive or negative.
# A support value close to 1 means strong positive support.
# A support value close to 0 means strong negative support.
support_matrix <- OmegaCor[[1]]$support

# as tibble for inspection
support_matrix_df <- as_tibble(support_matrix, rownames = "Species1") %>%
  pivot_longer(
    cols = -Species1,
    names_to = "Species2",
    values_to = "Support"
  ) %>% save_kable(
    file.path(output_dir, "hmsc_figures" , "Species_Correlation_Support_Levels.pdf")
  )

# Define the significance level. 0.95 is a common choice.
supportLevel <- 0.95

# Create a "p-value" matrix for masking. We will consider an association
# non-significant if its support is between (1 - supportLevel) and supportLevel.
# For supportLevel = 0.95, this is the interval [0.05, 0.95].
p_matrix <- ifelse(support_matrix > supportLevel | support_matrix < (1 - supportLevel), 0.01, 1)

# Create a mapping from original names to pretty names
name_map <- species_map %>%
  dplyr::filter(original_name_in_data %in% colnames(cor_matrix_mean)) %>%
  select(original_name_in_data, pretty_name) %>%
  deframe()

# Get the new pretty names in the correct order of the matrix
new_labels <- name_map[colnames(cor_matrix_mean)]

# Rename the rows and columns of both matrices to use pretty names
colnames(cor_matrix_mean) <- new_labels
rownames(cor_matrix_mean) <- new_labels
colnames(support_matrix) <- new_labels
rownames(support_matrix) <- new_labels


# --- 12c. Create the Correlation Plot ---
# hc.order = TRUE automatically reorders the species to group similar ones together.
# library(RColorBrewer) # For color palettes
# display.brewer.pal(n = 3, name = "BrBG") # Check the color palette
# pal <- brewer.pal(n = 3, name = "BrBG") # Use this palette for the plot
omega_plot <- ggcorrplot(
  cor_matrix_mean,
  hc.order = TRUE,          # Hierarchically cluster the species
  type = "lower",           # Show only the lower triangle of the matrix
  p.mat = p_matrix,         # Use the support matrix to mask non-significant cells
  insig = "blank",          # Hide non-significant associations completely
  lab = FALSE,              # Do not show the correlation values as text
  colors = c("#5AB4AC", "#F5F5F5", "#D8B365"), # Color scheme for neg -> zero -> pos
  title = "Significant Residual Species Co-occurrences (Support > 95%)",
  ggtheme = theme_minimal() # Use a clean theme
) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(face = "bold")
  ) +
  labs(fill = "Correlation")

# --- 12d. Save and Display the Plot ---
ggsave(
  file.path(output_dir,"hmsc_figures",  "species_associations_ggcorrplot.svg"),
  omega_plot,
  width = 9,
  height = 7,
  bg = "white"
)

print("Species co-occurrence plot (ggcorrplot) has been saved.")
print(omega_plot)

# print pdf with corr values
save_kable(
  as.data.frame(cor_matrix_mean) %>%
    rownames_to_column(var = "Species") %>%
    knitr::kable() %>%
    kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE),
  file.path(output_dir, "hmsc_figures" , "Species_Correlation_Matrix.pdf")
)

#===============================================================================
# 13. SPECIES CO-OCCURRENCE (NETWORK)
#===============================================================================
print("--- 13. Visualizing Species Co-occurrence as a Network ---")

# --- 7a. Compute the Association Matrix ---
OmegaCor <- computeAssociations(m1)

# --- 7b. Prepare Data for the Network Graph ---

# Get the matrix of mean posterior correlations
cor_matrix_mean <- OmegaCor[[1]]$mean

# Get the support matrix
support_matrix <- OmegaCor[[1]]$support

# Create an adjacency matrix, setting non-significant correlations to 0
adj_matrix <- cor_matrix_mean
adj_matrix[support_matrix < 0.95 & support_matrix > 0.05] <- 0
diag(adj_matrix) <- 0 # No self-loops

# --- 7c. Create and Refine the igraph Object ---

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(igraph)

species_network <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# Add edge attributes
E(species_network)$layout_weight <- abs(E(species_network)$weight)
E(species_network)$association_type <- ifelse(E(species_network)$weight > 0, "Positive", "Negative")

# --- Add pretty_name as vertex attribute ---
# Map original species names to pretty_name using species_map
vertex_names <- V(species_network)$name
pretty_labels <- species_map$pretty_name[match(vertex_names, species_map$original_name_in_data)]
V(species_network)$pretty_label <- pretty_labels

# Remove isolated nodes
species_network <- delete_vertices(species_network, degree(species_network) == 0)

# --- 7d. Create the Network Plot using ggraph ---
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
library(ggraph)
set.seed(1304)  # For reproducible layout
network_plot <- ggraph(species_network, layout = "fr", weights = layout_weight) +
  geom_edge_link(aes(color = association_type, width = layout_weight), alpha = 0.8) +
  geom_node_point(size = 5, color = "black") +
  geom_node_text(aes(label = pretty_label), repel = TRUE, size = 3.5) +
  scale_edge_color_manual(values = c("Negative" = "#E41A1C", "Positive" = "#4DAF4A"), name = "Association Type") +
  scale_edge_width(range = c(0.5, 2.5), name = "Strength (Abs. Correlation)") +
  labs(
    title = "Significant Residual Species Associations",
    subtitle = "Nodes are species (pretty names), edges represent positive or negative co-occurrence"
  ) +
  theme_graph(base_family = 'sans')

# --- 7e. Save and Display the Plot ---
ggsave(
  file.path(output_dir,"hmsc_figures", "species_network_plot.svg"),
  network_plot,
  width = 12,
  height = 8,
  bg = "white"
)

print("Species association network plot has been saved.")
print(network_plot)

ggsave(file.path(output_dir, "hmsc_figures", "species_association_network.svg"), network_plot, width = 12, height = 12)
print(network_plot)


#===============================================================================