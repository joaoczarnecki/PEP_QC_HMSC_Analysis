# ---
# title: "A Complete Workflow for Community-Assisted Migration Analysis using JSDMs"
# author: "Liz, JPC."
# date: "2025-08-12"
# ---

#===============================================================================
# 0. SETUP AND PACKAGE LOADING
#===============================================================================
# Ensure all necessary packages are installed and loaded
packages_to_install <- c("tidyverse", "Hmsc", "coda", "doParallel", "sf", 
                         "ggraph", "igraph", "knitr", "kableExtra", "ggcorrplot", 
                         "shiny", "shinythemes", "leaflet", "DT", "terra", "abind")
new_packages <- packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cran.rstudio.com/")

# Load libraries
library(tidyverse)
library(Hmsc)
library(coda)
library(doParallel)
library(sf)
library(ggraph)
library(igraph)
library(knitr)
library(kableExtra)
library(ggcorrplot)
library(shiny)
library(shinythemes)
library(leaflet)
library(DT)
library(terra)
library(abind)

#===============================================================================
# 1. DEFINE PATHS AND CREATE DICTIONARIES
#===============================================================================
print("--- 1. Defining paths and creating dictionaries ---")

# --- Define Paths ---
output_dir <- "/Thesis/3rdChapter/PEP_QC/results/undisturbed/"
shiny_app_dir <- file.path(output_dir, "CAMS_Shiny_App_PA_Advanced")
dir.create(shiny_app_dir, showWarnings = FALSE, recursive = TRUE)

# --- Species Mapping Dictionary ---
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

# --- Potential Vegetation (Végétation Potentielle) Dictionary ---
type_eco_descriptions <- tibble::tribble(
  ~type_eco_prefix, ~Description,
  "FC1", "Chênaie rouge (Red Oak Forest)",
  "FE1", "Érablière à caryer cordiforme (Sugar Maple - Bitternut Hickory Forest)",
  "FE2", "Érablière à tilleul (Sugar Maple - Basswood Forest)",
  "FE3", "Érablière à bouleau jaune (Sugar Maple - Yellow Birch Forest)",
  "FE4", "Érablière à bouleau jaune et hêtre (Sugar Maple - Yellow Birch and Beech Forest)",
  "FE5", "Érablière à ostryer (Sugar Maple - Hop Hornbeam Forest)",
  "FE6", "Érablière à chêne rouge (Sugar Maple - Red Oak Forest)",
  "FO1", "Ormaie à frêne noir (Black Ash Grove)",
  "LA1", "Lande à lichens (ou à mousses) (Lichen or Moss Heathland)",
  "LA2", "Lande arbustive (Shrubby Heathland)",
  "LA3", "Lande herbacée (Herbaceous Heathland)",
  "LA4", "Lande rocheuse (Rocky Heathland)",
  "LI1", "Littoral (Littoral Zone)",
  "LL1", "Lande alpine à lichens (ou à mousses) (Alpine Lichen or Moss Heathland)",
  "LL2", "Lande alpine arbustive (Alpine Shrubby Heathland)",
  "LL3", "Lande alpine herbacée (Alpine Herbaceous Heathland)",
  "LL4", "Lande alpine rocheuse (Alpine Rocky Heathland)",
  "LM1", "Lande maritime à lichens (ou à mousses) (Maritime Lichen or Moss Heathland)",
  "LM2", "Lande maritime arbustive (Maritime Shrubby Heathland)",
  "LM3", "Lande maritime herbacée (Maritime Herbaceous Heathland)",
  "LM4", "Lande maritime rocheuse (Maritime Rocky Heathland)",
  "MA1", "Marais ou marécage arbustif d’eau douce (Freshwater Shrubby Marsh or Swamp)",
  "MA2", "Marais ou marécage arbustif d’eau saumâtre ou salée (Brackish or Saltwater Shrubby Marsh or Swamp)",
  "ME1", "Pessière noire à peuplier faux-tremble (Black Spruce - Trembling Aspen Forest)",
  "MF1", "Frênaie noire à sapin (Black Ash - Balsam Fir Forest)",
  "MJ1", "Bétulaie jaune à sapin et érable à sucre (Yellow Birch - Balsam Fir and Sugar Maple Forest)",
  "MJ2", "Bétulaie jaune à sapin (Yellow Birch - Balsam Fir Forest)",
  "MS1", "Sapinière à bouleau jaune (Balsam Fir - Yellow Birch Forest)",
  "MS2", "Sapinière à bouleau blanc (Balsam Fir - Paper Birch Forest)",
  "MS4", "Sapinière à bouleau blanc montagnarde (Montane Balsam Fir - Paper Birch Forest)",
  "MS6", "Sapinière à érable rouge (Balsam Fir - Red Maple Forest)",
  "MS7", "Sapinière à bouleau blanc maritime (Maritime Balsam Fir - Paper Birch Forest)",
  "RB1", "Pessière blanche ou cédrière issue d’agriculture (White Spruce or Eastern White Cedar Forest from Agriculture)",
  "RB2", "Pessière blanche maritime (Maritime White Spruce Forest)",
  "RB3", "Pessière blanche subalpine ou sapinière à épinette blanche subalpine (Subalpine White Spruce or Balsam Fir - White Spruce Forest)",
  "RB4", "Pessière blanche montagnarde (Montane White Spruce Forest)",
  "RB5", "Pessière blanche ou sapinière à bouleau blanc de l’île d’Anticosti (White Spruce or Balsam Fir - Paper Birch Forest of Anticosti Island)",
  "RC3", "Cédrière tourbeuse à sapin (Peaty Eastern White Cedar - Balsam Fir Forest)",
  "RE1", "Pessière noire à lichens (Black Spruce - Lichen Forest)",
  "RE2", "Pessière noire à mousses ou à éricacées (Black Spruce - Moss or Heath Forest)",
  "RE3", "Pessière noire à sphaignes (Black Spruce - Sphagnum Forest)",
  "RE4", "Pessière noire à mousses ou à éricacées montagnarde (Montane Black Spruce - Moss or Heath Forest)",
  "RE7", "Pessière noire maritime (Maritime Black Spruce Forest)",
  "RE8", "Pessière noire subalpine (Subalpine Black Spruce Forest)",
  "RI1", "Rive (Riparian Zone)",
  "RP1", "Pinède blanche ou pinède rouge (White or Red Pine Forest)",
  "RS1", "Sapinière à thuya (Balsam Fir - Eastern White Cedar Forest)",
  "RS2", "Sapinière à épinette noire (Balsam Fir - Black Spruce Forest)",
  "RS3", "Sapinière à épinette noire et sphaignes (Balsam Fir - Black Spruce and Sphagnum Forest)",
  "RS4", "Sapinière à épinette noire montagnarde (Montane Balsam Fir - Black Spruce Forest)",
  "RS5", "Sapinière à épinette rouge (Balsam Fir - Red Spruce Forest)",
  "RS7", "Sapinière à épinette noire maritime (Maritime Balsam Fir - Black Spruce Forest)",
  "RT1", "Prucheraie (Hemlock Forest)",
  "SM1", "Sables mobiles (Shifting Sands)",
  "SM2", "Sables mobiles maritimes (Maritime Shifting Sands)",
  "TOB", "Tourbière ombrotrophe (Bog)",
  "TOF", "Tourbière minérotrophe (Fen)",
  "TOU", "Tourbière indifférenciée (minérotrophe ou ombrotrophe) (Indifferentiated Peatland (Minerotrophic or Ombrotrophic))"
)



# type_eco_species <- tibble::tribble(
#   ~type_eco_prefix, ~dominant_species, ~secondary_species,

#   "FC1", c("Quercus rubra"), c("Acer rubrum","Acer saccharum","Betula papyrifera","Pinus strobus"),

#   "FE1", c("Acer saccharum","Carya cordiformis"), 
#           c("Tilia americana","Fraxinus americana","Ulmus americana","Ostrya virginiana"),

#   "FE2", c("Acer saccharum","Tilia americana"),
#           c("Fagus grandifolia","Fraxinus americana","Ulmus americana","Ostrya virginiana"),

#   "FE3", c("Acer saccharum","Betula alleghaniensis"),
#           c("Fagus grandifolia","Tsuga canadensis","Abies balsamea","Betula papyrifera","Fraxinus americana"),

#   "FE4", c("Acer saccharum","Betula alleghaniensis","Fagus grandifolia"),
#           c("Tsuga canadensis","Abies balsamea"),

#   "FE5", c("Acer saccharum","Ostrya virginiana"),
#           c("Tilia americana","Fraxinus americana","Ulmus americana"),

#   "FE6", c("Acer saccharum","Quercus rubra"),
#           c("Betula alleghaniensis","Betula papyrifera","Pinus strobus","Acer rubrum"),

#   "FO1", c("Fraxinus nigra","Ulmus americana"),
#           c("Acer rubrum","Fraxinus pennsylvanica","Thuja occidentalis"),

#   "MF1", c("Fraxinus nigra","Abies balsamea"),
#           c("Acer rubrum","Ulmus americana","Thuja occidentalis"),

#   "MJ1", c("Betula alleghaniensis","Abies balsamea","Acer saccharum"),
#           c(),

#   "MJ2", c("Betula alleghaniensis","Abies balsamea"),
#           c("Acer saccharum","Picea glauca"),

#   "MS1", c("Abies balsamea","Betula alleghaniensis"),
#           c("Acer saccharum","Picea glauca"),

#   "MS2", c("Abies balsamea","Betula papyrifera"),
#           c("Picea glauca","Picea mariana"),

#   "MS4", c("Abies balsamea","Betula papyrifera"),
#           c("Picea mariana","Sorbus americana"),

#   "MS6", c("Abies balsamea","Acer rubrum"),
#           c("Betula papyrifera","Picea glauca"),

#   "MS7", c("Abies balsamea","Betula papyrifera"),
#           c("Picea glauca","Picea mariana"),

#   "ME1", c("Picea mariana","Populus tremuloides"),
#           c("Betula papyrifera","Abies balsamea"),

#   "RB1", c("Picea glauca","Thuja occidentalis"),
#           c("Abies balsamea","Betula papyrifera","Acer rubrum"),

#   "RB2", c("Picea glauca"),
#           c("Abies balsamea","Betula papyrifera","Acer rubrum","Picea mariana"),

#   "RB3", c("Picea glauca","Abies balsamea"),
#           c("Betula papyrifera","Picea mariana"),

#   "RB4", c("Picea glauca","Abies balsamea"),
#           c("Betula papyrifera","Picea mariana"),

#   "RB5", c("Picea glauca","Abies balsamea","Betula papyrifera"),
#           c("Picea rubens","Picea mariana"),

#   "RC3", c("Thuja occidentalis","Abies balsamea"),
#           c("Picea mariana","Larix laricina"),

#   "RE1", c("Picea mariana"),
#           c("Pinus banksiana","Betula papyrifera"),

#   "RE2", c("Picea mariana"),
#           c("Abies balsamea","Pinus banksiana","Betula papyrifera"),

#   "RE3", c("Picea mariana","Larix laricina"),
#           c("Abies balsamea","Thuja occidentalis"),

#   "RE4", c("Picea mariana","Abies balsamea"),
#           c("Betula papyrifera","Sorbus americana"),

#   "RE7", c("Picea mariana"),
#           c("Picea glauca","Abies balsamea","Betula papyrifera"),

#   "RE8", c("Picea mariana","Abies balsamea"),
#           c("Betula papyrifera","Sorbus americana"),

#   "RS1", c("Abies balsamea","Thuja occidentalis"),
#           c("Picea glauca","Betula papyrifera","Acer rubrum"),

#   "RS2", c("Abies balsamea","Picea mariana"),
#           c("Betula papyrifera","Picea glauca"),

#   "RS3", c("Abies balsamea","Picea mariana"),
#           c("Larix laricina","Thuja occidentalis"),

#   "RS4", c("Abies balsamea","Picea mariana"),
#           c("Betula papyrifera","Sorbus americana"),

#   "RS5", c("Abies balsamea","Picea rubens"),
#           c("Picea mariana","Betula papyrifera"),

#   "RS7", c("Abies balsamea","Picea mariana","Picea glauca"),
#           c("Betula papyrifera","Acer rubrum"),

#   "RP1", c("Pinus strobus","Pinus resinosa"),
#           c("Acer rubrum","Quercus rubra","Betula papyrifera"),

#   # Heathlands, littoral, dunes (tree-scarce)
#   "LA1", c(), c("Picea mariana","Pinus banksiana","Betula papyrifera"),
#   "LA2", c(), c("Picea mariana","Pinus banksiana"),
#   "LA3", c(), c("Betula papyrifera"),
#   "LA4", c(), c("Picea mariana","Abies balsamea","Betula papyrifera"),

#   "LL1", c(), c("Picea mariana","Abies balsamea","Betula papyrifera","Sorbus americana"),
#   "LL2", c(), c("Picea mariana","Abies balsamea","Betula papyrifera","Sorbus americana"),
#   "LL3", c(), c("Picea mariana","Abies balsamea"),
#   "LL4", c(), c("Picea mariana","Abies balsamea"),

#   "LM1", c(), c("Picea mariana","Picea glauca","Abies balsamea","Betula papyrifera"),
#   "LM2", c(), c("Picea mariana","Picea glauca","Betula papyrifera"),
#   "LM3", c(), c("Betula papyrifera","Acer rubrum","Picea mariana"),
#   "LM4", c(), c("Picea glauca","Picea mariana","Betula papyrifera","Abies balsamea"),

#   "LI1", c(), c("Acer saccharinum","Populus balsamifera","Betula papyrifera",
#                 "Fraxinus pennsylvanica","Salix spp.","Picea mariana","Abies balsamea"),

#   "SM1", c(), c("Pinus banksiana","Pinus strobus","Populus tremuloides","Betula papyrifera"),
#   "SM2", c(), c("Pinus banksiana","Picea mariana","Picea glauca"),

#   "MA1", c(), c("Larix laricina","Picea mariana","Abies balsamea","Salix spp.","Betula papyrifera"),
#   "MA2", c(), c("Picea mariana","Picea glauca","Larix laricina"),

#   "TOB", c(), c("Picea mariana","Larix laricina","Thuja occidentalis"),
#   "TOF", c(), c("Larix laricina","Picea mariana","Betula papyrifera","Alnus incana"),
#   "TOU", c(), c("Picea mariana","Larix laricina","Thuja occidentalis","Betula papyrifera")
# )



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


#===============================================================================
# 2. DATA PREPARATION FOR PRESENCE-ABSENCE MODEL
#===============================================================================
print("--- 2. Loading and preparing data for Presence-Absence model ---")

# --- Load Source Data ---

plots_sf_env_clim <- read_csv("/Thesis/3rdChapter/PEP_QC/databases/undisturbed_plots_sf_env_clim_4th_inv.csv")
station_pe_4th_inv <- st_read("/Thesis/3rdChapter/PEP_QC/PEP.gpkg", layer = "station_pe") %>% st_drop_geometry()
# pep_gpkg <- st_read("/Thesis/3rdChapter/PEP_QC/PEP.gpkg")

# write_sf(station_pe_4th_inv, "/Thesis/3rdChapter/PEP_QC/databases/station_pe_4th_inv.gpkg", layer = "station_pe", delete_layer = TRUE)
# --- Define Predictors and Species ---

# env_predictors <- c("MAP", "MAT", "OM", "slope", "drainage_class")
env_predictors <- c( "MAP", "MAT", "pH",
                    "Clay", "OM", "rtp", "twi",
                    "slope", "drainage_class")

species_abun_vars <- plots_sf_env_clim %>%
  select(ends_with("_n") & !starts_with("NA")) %>%
  colnames()

# --- Create and Filter Species Abundance Matrix (Y) ---
Y_full <- plots_sf_env_clim %>%
  select(plot_id, all_of(species_abun_vars)) %>%
  mutate(across(all_of(species_abun_vars), ~ as.integer(replace_na(., 0)))) %>%
  column_to_rownames(var = "plot_id") %>%
  as.matrix()

prevalence <- colSums(Y_full > 0, na.rm = TRUE)
prevalence_threshold <- nrow(Y_full) * 0.05
species_to_keep <- names(prevalence[prevalence >= prevalence_threshold])
Y_abundance <- Y_full[, species_to_keep]

# --- Convert to Presence-Absence ---
Y_pa <- (Y_abundance > 0) * 1

# --- Create Environmental Covariate Data Frame (XData) ---
XData <- plots_sf_env_clim %>%
  select(plot_id, all_of(env_predictors)) %>%
  column_to_rownames(var = "plot_id")
# as factor drainage_class
# XData$drainage_class <- as.factor(XData$drainage_class)

# --- Create Potential Vegetation Mapping ---
type_eco_mapping <- station_pe_4th_inv %>%
  mutate(type_eco_prefix = str_sub(type_eco, 1, 3)) %>%
  group_by(id_pe) %>%
  summarise(type_eco = first(na.omit(type_eco)),
            type_eco_prefix = first(na.omit(type_eco_prefix)),
            .groups = "drop") %>%
  dplyr::filter(!is.na(type_eco))

# --- Align All Data Objects ---
common_plots <- intersect(rownames(Y_pa), rownames(na.omit(XData)))
common_plots <- intersect(common_plots, type_eco_mapping$id_pe)

Y_pa <- Y_pa[common_plots, ]
XData <- XData[common_plots, ]
type_eco_aligned <- type_eco_mapping %>% dplyr::filter(id_pe %in% common_plots)
Y_pa <- Y_pa[, colSums(Y_pa) > 0]

print(paste("Data aligned. Final dimensions: Y_pa =", nrow(Y_pa), "x", ncol(Y_pa)))

#===============================================================================
# 3. HMSC MODEL (PROBIT) LOAD EXISTING MODEL
#===============================================================================

#--- Load the fitted model ---
# Load the first model (assuming it's already loaded as m1)
# Hmsc object with 5572 sampling units, 13 species, 26 covariates, 1 traits and 1 random levels
# Posterior MCMC sampling with 4 chains each with 1000 samples, thin 100 and transient 50000
load(file.path(output_dir, "hmsc", "hmsc_presence_model.Rdata"))  # Loads as m1
mod1 <- m1  # Just to clarify


# Check objects 

print(m1)

print("--- Comparing models using WAIC ---")

# Compare models with waic
waic_mod1 <- computeWAIC(mod1, byColumn = FALSE)

print(paste("WAIC for Model 1:", waic_mod1))

print("--- 3. Setting up and fitting the Hmsc PROBIT model ---")


#===============================================================================
# 4. PREPARE DATA FOR SHINY APP
#===============================================================================
print("--- 4. Preparing and saving all objects for the Shiny app ---")

shiny_app_dir <- file.path(output_dir, "CAMS_Shiny_App_PA_Advanced")
dir.create(shiny_app_dir, showWarnings = FALSE, recursive = TRUE)

# --- 4.1 Prepare Future Climate Scenarios ---
fut_clim_raw <- read_csv("/Thesis/3rdChapter/PEP_QC/climate/plots_coord_12_GCMsY.csv")
gcm_scenarios <- unique(fut_clim_raw$GCM)
# gcm_scenarios <- gcm_scenarios[gcm_scenarios %in% c("8GCMs_ensemble_ssp370_2071-2100.gcm")] # fast test_include only ssp3 7.0 2071-2100 scenarios
XData_future_list <- list()
for (scenario in gcm_scenarios) {
  fut_clim_scenario <- fut_clim_raw %>%
    dplyr::filter(GCM == scenario) %>%
    select(id1, MAT, MAP) %>%
    rename(plot_id = id1) %>%
    mutate(plot_id = as.character(plot_id))
  
  XData_future_list[[scenario]] <- XData %>%
    as_tibble(rownames = "plot_id") %>%
    left_join(fut_clim_scenario, by = "plot_id") %>%
    mutate(MAT = ifelse(is.na(MAT.y), MAT.x, MAT.y), MAP = ifelse(is.na(MAP.y), MAP.x, MAP.y)) %>%
    select(-ends_with(".x"), -ends_with(".y")) %>%
    column_to_rownames("plot_id")
}
XData_future_list[["Current"]] <- XData

# --- 4.2 Pre-calculate Predictions with Uncertainty ---
predictions_file <- file.path(shiny_app_dir, "all_scenario_predictions.RData")

# Sequential (no parallelization) prediction + summarise
if (!file.exists(predictions_file)) {
  suppressPackageStartupMessages({
    library(Hmsc)
    library(abind)
  })

  predictions_dir <- dirname(predictions_file)
  dir.create(predictions_dir, showWarnings = FALSE, recursive = TRUE)

  # 1) Produce and save raw posterior predictions per scenario (separate files)
  for (scenario_name in names(XData_future_list)) {
    message(sprintf("Predicting for scenario: %s", scenario_name))

    X_new <- XData_future_list[[scenario_name]]

    # Align factor levels and column order to the fitted model
    fac_cols <- names(Filter(is.factor, m1$XData))
    for (fc in fac_cols) {
      if (fc %in% colnames(X_new)) {
        X_new[[fc]] <- factor(X_new[[fc]], levels = levels(m1$XData[[fc]]))
      }
    }
    stopifnot(all(colnames(m1$XData) %in% colnames(X_new)))
    X_new <- X_new[, colnames(m1$XData), drop = FALSE]

    # Predict fixed-effects-only (marginal) probabilities
    preds <- predict(
      m1,
      XData           = X_new,
      expected        = TRUE,
      predictEtaMean  = TRUE
    )

    # Save per-scenario raw posterior predictions
    scenario_file <- file.path(predictions_dir, paste0(scenario_name, "_preds.Rda"))
    save(preds, file = scenario_file)

    rm(preds, X_new); gc()
  }

  # 2) Load per-scenario prediction files and summarise to mean / 2.5% / 97.5%
  all_predictions_summary <- list()
  for (scenario_name in names(XData_future_list)) {
    scenario_file <- file.path(predictions_dir, paste0(scenario_name, "_preds.Rda"))
    if (!file.exists(scenario_file)) {
      warning("Missing prediction file for scenario: ", scenario_name); next
    }
    load(scenario_file) # loads 'preds'

    # preds is expected to be a list (one element per posterior draw) of [n_sites x n_species] matrices
    pred_array <- abind::abind(preds, along = 3)           # [n_sites x n_species x n_draws]
    pred_mean  <- apply(pred_array, c(1, 2), mean)
    pred_lower <- apply(pred_array, c(1, 2), quantile, probs = 0.025)
    pred_upper <- apply(pred_array, c(1, 2), quantile, probs = 0.975)

    all_predictions_summary[[scenario_name]] <- list(
      mean  = as.data.frame(pred_mean),
      lower = as.data.frame(pred_lower),
      upper = as.data.frame(pred_upper)
    )

    rm(preds, pred_array, pred_mean, pred_lower, pred_upper); gc()
  }

  # 3) Save consolidated summary
  save(all_predictions_summary, file = predictions_file)

} else {
  message(sprintf("Loading existing predictions summary from: %s", predictions_file))
  load(predictions_file)  # loads all_predictions_summary
}

save(all_predictions_summary, file = file.path(output_dir, "all_predictions_summary.RData"))
#Read prediction summary
load(file.path(output_dir, "all_predictions_summary.RData"))  # loads all_predictions_summary

# # COMPUTE SPECIES RICHNESS (S), COMMUNITY WEIGHTED MEANS (predT),
# # REGIONS OF COMMON PROFILE (RCP)
# S=rowSums(current_preds)
# str(current_preds)
# predT = (current_preds[[1]]%*%m$Tr)/matrix(rep(S,m$nt),ncol=m$nt)
# RCP = kmeans(current_preds[[1]], 7)
# RCP$cluster = as.factor(RCP$cluster)
# # EXTRACT THE OCCURRENCE PROBABILITIES OF ONE EXAMPLE SPECIES
# pred_Cm = predYR[,50]
# # MAKE A DATAFRAME OF THE DATA TO BE PLOTTED
# mapData=data.frame(xy,S,predT,pred_Cm,RCP$cluster)


#read all prediction summary files
# load("H:/Thesis/3rdChapter/PEP_QC/results/undisturbed/CAMS_Shiny_App_PA_Advanced/all_predictions_summary.RData")

# --- 4.3 Prepare Supporting Data Objects ---
library(terra)
install.packages("terra")

# Build TYPE_ECO_list from type_eco_descriptions (uses dominant_codes + secondary_codes)
library(purrr)
TYPE_ECO_list <- type_eco_descriptions %>%
  mutate(
    all_codes = map2(dominant_codes, secondary_codes, ~ unique(c(.x, .y))),
    all_codes = map(all_codes, ~ .x[!is.na(.x) & nzchar(.x)]),
    # convert species codes to the original names used in the data (e.g. "AEH" -> "AEH_n")
    species_list = map(all_codes, ~ {
      orig_names <- paste0(.x, "_n")
      # keep only species actually present in the predictions / data (safe filtering)
      orig_names[orig_names %in% colnames(Y_pa)]
    })
  ) %>%
  select(type_eco_prefix, species_list) %>%
  deframe()

plots_for_map <- plots_sf_env_clim %>%
  dplyr::filter(plot_id %in% common_plots) %>%
  select(plot_id, longitude, latitude)

study_area_sf <- st_as_sf(plots_for_map, coords = c("longitude", "latitude"), crs = 4326)
study_area_raster <- rast(ext(study_area_sf) + 0.1, resolution = 0.05)
crs(study_area_raster) <- "EPSG:4326"
# plot(study_area_raster)
# --- 4.4 Save All Objects for the App ---
save(all_predictions_summary, file = file.path(shiny_app_dir, "all_predictions_summary.RData"))
save(TYPE_ECO_list, type_eco_descriptions, file = file.path(shiny_app_dir, "type_eco_data.RData"))
save(species_map, file = file.path(shiny_app_dir, "species_map.RData"))
save(plots_for_map, file = file.path(shiny_app_dir, "plot_coordinates.RData"))
save(study_area_raster, file = file.path(shiny_app_dir, "basemap_raster.RData"))
save(XData_future_list, file = file.path(shiny_app_dir, "XData_future_list.RData"))
print(paste("All necessary data objects for the advanced Shiny app have been saved to:", shiny_app_dir))


#===============================================================================
print("Script completed successfully.") # END OF SCRIPT
#===============================================================================

# Plot densities of MAP and MAT for all scenarios in two panels
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare data in long format
all_data <- list()
for (scenario_name in names(XData_future_list)) {
  scenario_data <- XData_future_list[[scenario_name]] %>%
    mutate(scenario = scenario_name) %>%
    select(MAP, MAT, scenario)
  all_data <- bind_rows(all_data, scenario_data)
}

all_data_long <- all_data %>%
  pivot_longer(cols = c(MAP, MAT), names_to = "variable", values_to = "value")

# Create the plot with two panels
p <- ggplot(all_data_long, aes(x = value, fill = scenario)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Density of MAP and MAT across scenarios",
       x = "Value",
       y = "Density") +
  scale_fill_viridis_d() +
  theme_minimal()

print(p)
# save in output dir hmsc_figures
ggsave::ggsave("MAP_densities_all_scenarios.png", path = file.path(output_dir, "hmsc_figures"))

# create a plot of interval ranges for MAP and MAT across scenarios
interval_data <- all_data_long %>%
  group_by(scenario, variable) %>%
  summarise(
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    mean_val = mean(value, na.rm = TRUE),
    median_val = median(value, na.rm = TRUE),
    .groups = "drop"
  )
p_intervals <- ggplot(interval_data, aes(x = scenario, color = scenario)) +
  geom_linerange(aes(ymin = min, ymax = max), size = 1) +
  geom_point(aes(y = mean_val), shape = 16, size = 3) +  # Point for mean
  geom_point(aes(y = median_val), shape = 17, size = 3) +  # Point for median (triangle)
  # Label scenario name beside the line, rotated vertically
  geom_text(aes(y = max, label = scenario), angle = 90, vjust = -0.5, hjust = 0.5, size = 3, show.legend = FALSE) +  
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Interval Ranges of MAP and MAT across scenarios",
       x = "Scenario",
       y = "Value") +
  scale_color_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

print(p_intervals)
ggsave("MAP_MAT_interval_ranges_all_scenarios.png", path = file.path(output_dir, "hmsc_figures"), width = 10, height = 7)


# Exploratory Analysis to see difference in predicted probabilities
# create a table with mean predicted probabilities for each species, scenario every 0.5 latitude degree band
library(dplyr)

mean_probabilities <- data.frame()
for (scenario in names(all_predictions_summary)) {
  preds_mean <- all_predictions_summary[[scenario]]$mean
  preds_mean$plot_id <- rownames(preds_mean)
  
  preds_long <- preds_mean %>%
    pivot_longer(cols = -plot_id, names_to = "species", values_to = "mean_probability") %>%
    left_join(plots_for_map, by = c("plot_id" = "plot_id")) %>%
    mutate(lat_band = floor(latitude * 2) / 2,
           scenario = scenario)
  
  mean_probabilities <- rbind(mean_probabilities, preds_long)
}

# Create a table summarizing mean predicted probabilities per species scenario every 0.5 latitude degree band
mean_prob_summary <- mean_probabilities %>%
  group_by(species, scenario, lat_band) %>%
  summarise(mean_probability = mean(mean_probability, na.rm = TRUE), .groups = "drop")

# Compute differences: for each future scenario, diff = future - current
current_scenario <- "Current"
future_scenarios <- setdiff(names(all_predictions_summary), current_scenario)

diff_prob_summary <- mean_prob_summary %>%
  pivot_wider(names_from = scenario, values_from = mean_probability) %>%
  mutate(across(all_of(future_scenarios), ~ . - .data[[current_scenario]], .names = "diff_{.col}")) %>%
  select(species, lat_band, starts_with("diff_")) %>%
  pivot_longer(cols = starts_with("diff_"), names_to = "scenario", values_to = "diff_probability") %>%
  mutate(scenario = str_remove(scenario, "^diff_")) %>%
  filter(!is.na(diff_probability))

# Plot only probability differences across latitude bands (no climate)
library(ggplot2)
ggplot(diff_prob_summary, aes(x = lat_band, y = diff_probability, color = scenario)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ species, scales = "fixed") +
  labs(
    title = "Differences in Mean Predicted Probabilities across Latitude Bands",
    x = "Latitude Band (°)",
    y = "Difference in Mean Predicted Probability"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")




# plot betas to see which variables have the most influence
# Extract beta coefficients

beta_coeffs <- getPostEstimate(m1, parName = "Beta")
beta_means <- beta_coeffs$mean
beta_lower <- beta_coeffs$supportNeg  # 2.5% quantile
beta_upper <- beta_coeffs$support # 97.5% quantile
beta_df <- as.data.frame(beta_means)
beta_df$predictor <- m1$covNames
beta_long <- beta_df %>%
  pivot_longer(cols = -predictor, names_to = "species", values_to = "beta_mean")

# Add lower and upper bounds
beta_lower_df <- as.data.frame(beta_lower)
beta_lower_df$predictor <- m1$covNames
beta_lower_long <- beta_lower_df %>%
  pivot_longer(cols = -predictor, names_to = "species", values_to = "beta_lower")

beta_upper_df <- as.data.frame(beta_upper)
beta_upper_df$predictor <- m1$covNames
beta_upper_long <- beta_upper_df %>%
  pivot_longer(cols = -predictor, names_to = "species", values_to = "beta_upper")

# Combine into one data frame
beta_combined <- beta_long %>%
  left_join(beta_lower_long, by = c("predictor", "species")) %>%
  left_join(beta_upper_long, by = c("predictor", "species"))

# Plot beta coefficients with error bars
ggplot(beta_combined, aes(x = predictor, y = beta_mean, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = beta_lower, ymax = beta_upper), width = 0.2, position = position_dodge(0.9)) +
  labs(title = "Beta Coefficients for Environmental Predictors",
       x = "Environmental Predictor",
       y = "Beta Coefficient") +
  theme_minimal() +
  theme(legend.position = "bottom")

# scatter plot comparing mean probabilities between Current and 8GCMs_ensemble_ssp370_2071-2100 scenarios
library(ggplot2)
# Filter to Current and 8GCMs_ensemble_ssp370_2071-2100 scenarios
current_scenario <- "Current"
future_scenario <- "8GCMs_ensemble_ssp370_2071-2100.gcm"  # Adjust if the exact name differs

scatter_data <- mean_probabilities %>%
  filter(scenario %in% c(current_scenario, future_scenario)) %>%
  pivot_wider(names_from = scenario, values_from = mean_probability) %>%
  rename(Current = !!current_scenario, Future = !!future_scenario)

ggplot(scatter_data, aes(x = Current, y = Future, color = type_eco_prefix)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Scatter Plot: Current vs Future Mean Predicted Probabilities",
       x = "Current Mean Probability",
       y = "Future Mean Probability") +
  theme_minimal() +
  theme(legend.position = "bottom")




#===============================================================================
# Analysis of pobability changes across communities and climate scenarios
#===============================================================================

# Load necessary libraries
library(tidyverse)
library(Hmsc)
library(coda)
library(sf)
library(terra)
library(abind)
# Load the saved predictions summary
load("H:/Thesis/3rdChapter/PEP_QC/results/undisturbed/CAMS_Shiny_App_PA_Advanced/all_predictions_summary.RData")
# Load species mapping
load("H:/Thesis/3rdChapter/PEP_QC/results/undisturbed/CAMS_Shiny_App_PA_Advanced/species_map.RData")
# Load plot coordinates
load("H:/Thesis/3rdChapter/PEP_QC/results/undisturbed/CAMS_Shiny_App_PA_Advanced/plot_coordinates.RData")
# Load ecological type data
load("H:/Thesis/3rdChapter/PEP_QC/results/undisturbed/CAMS_Shiny_App_PA_Advanced/type_eco_data.RData")  



# ===========================
# Analysis: compare Current vs SSP3-7.0 for a selected potential vegetation
# ============================
library(tidyverse)
library(ggplot2)

# User-selectable: set the target potential vegetation prefix (e.g. "FE3", "FE2", "FE1", "FE4", "FE2", "FE23", etc.)
# Bioclimatic domains (Potential Vegetation codes) of southern and central Quebec:
# FE34 – Sugar maple–Bitternut hickory forest → warmest (thermophilous) domain
# FE33 – Sugar maple–Basswood forest → temperate, nutrient-rich domain
# FE32 – Sugar maple–Yellow birch forest → temperate mesothermal domain
# FE23 – Balsam fir–Yellow birch forest → mixed boreal domain
# FE22 – Balsam fir–Paper birch forest → boreal domain
# FE21 – Black spruce–Feather moss forest → subarctic boreal domain

# Publishable density plots of predicted probabilities by community and scenario
# - One density per panel (community x scenario)
# - Dashed vertical lines show 95% interval; solid line shows median
# - Scenarios ordered: Current, 2041-2070, 2071-2100

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)

# 1) Choose communities (potential vegetation prefix) to plot
# Edit this vector as needed (only those present will be used)
target_prefixes <- c("FE1","FE2","FE3","FE4","FE5","FE6")

# Keep only prefixes that exist in type_eco_aligned
available_prefixes <- intersect(target_prefixes, unique(type_eco_aligned$type_eco_prefix))
if (length(available_prefixes) == 0) stop("None of the chosen prefixes exist in 'type_eco_aligned'.")

# 2) Find scenario names
scenario_names <- names(all_predictions_summary)
if (is.null(scenario_names) || length(scenario_names) == 0) {
  stop("all_predictions_summary not found or has no names.")
}

# Find "Current"
current_name <- scenario_names[grepl("^Current$|^Baseline$|present|current", scenario_names, ignore.case = TRUE)]
if (length(current_name) == 0) {
  # fallback: first scenario
  current_name <- scenario_names[1]
} else {
  current_name <- current_name[1]
}

# Helper to find a scenario containing both ssp370 (or ssp3) and a year window
find_scn <- function(year_window_regex) {
  hit <- scenario_names[grepl("(ssp370|ssp3|SSP3)", scenario_names, ignore.case = TRUE) &
                          grepl(year_window_regex, scenario_names, ignore.case = TRUE)]
  if (length(hit)) hit[1] else character(0)
}

scn_2041_2070 <- find_scn("2041[-_]?2070")
scn_2071_2100 <- find_scn("2071[-_]?2100")

# If not found with strict pattern, relax a bit
if (!nzchar(scn_2041_2070)) {
  hit <- scenario_names[grepl("(2041|2070)", scenario_names) & grepl("(ssp370|ssp3|SSP3)", scenario_names, ignore.case = TRUE)]
  if (length(hit)) scn_2041_2070 <- hit[1]
}
if (!nzchar(scn_2071_2100)) {
  hit <- scenario_names[grepl("(2071|2100)", scenario_names) & grepl("(ssp370|ssp3|SSP3)", scenario_names, ignore.case = TRUE)]
  if (length(hit)) scn_2071_2100 <- hit[1]
}

scenario_order <- c(current_name, scn_2041_2070, scn_2071_2100) %>% .[nzchar(.)]
scenario_order <- unique(scenario_order) # remove duplicates if any
if (length(scenario_order) < 2) {
  stop("Could not identify the requested future scenarios (2041-2070 and 2071-2100). Check names(all_predictions_summary).")
}

# 3) Helper: subset a summary matrix (rows=plots, cols=species) to plots/species of a community
get_summary_matrix <- function(scn_name, part = "mean") {
  if (!scn_name %in% names(all_predictions_summary)) stop("Scenario not found: ", scn_name)
  mat_df <- all_predictions_summary[[scn_name]][[part]]
  if (is.null(rownames(mat_df))) rownames(mat_df) <- seq_len(nrow(mat_df))
  as.matrix(mat_df)
}

# 4) Build data for densities: vectorize mean probabilities for each (community, scenario)
dens_df <- list()

for (prefix in available_prefixes) {
  # plots in community
  plots_in_type <- type_eco_aligned %>%
    dplyr::filter(type_eco_prefix == prefix) %>%
    pull(id_pe) %>%
    as.character()
  if (length(plots_in_type) == 0) next

  # species present in those plots (based on Y_pa)
  species_in_type <- colnames(Y_pa)[colSums(Y_pa[rownames(Y_pa) %in% plots_in_type, , drop = FALSE] > 0, na.rm = TRUE) > 0]
  if (length(species_in_type) == 0) next

  for (scn in scenario_order) {
    mat <- get_summary_matrix(scn, "mean")
    rows_keep <- intersect(rownames(mat), plots_in_type)
    cols_keep <- intersect(colnames(mat), species_in_type)
    if (length(rows_keep) == 0 || length(cols_keep) == 0) next

    vals <- as.vector(mat[rows_keep, cols_keep, drop = FALSE])
    vals <- vals[is.finite(vals) & vals >= 0 & vals <= 1]
    if (!length(vals)) next

    dens_df[[length(dens_df) + 1]] <- tibble(
      value = vals,
      community = prefix,
      scenario = scn
    )
  }
}

dens_df <- bind_rows(dens_df)
if (nrow(dens_df) == 0) stop("No overlapping data to plot for the selected communities and scenarios.")

# Order scenarios for plotting
dens_df <- dens_df %>%
  mutate(
    scenario = factor(scenario, levels = scenario_order, labels = scenario_order),
    community = factor(community, levels = available_prefixes)
  )

# 5) Compute per-panel quantiles (median + 95% interval)
panel_q <- dens_df %>%
  group_by(community, scenario) %>%
  summarise(
    q025 = quantile(value, 0.025, na.rm = TRUE),
    med  = median(value, na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# 6) Build nicer facet labels for communities (prefix + description if available)
prefix_labels <- type_eco_descriptions %>%
  dplyr::filter(type_eco_prefix %in% levels(dens_df$community)) %>%
  mutate(label = paste0(type_eco_prefix, " — ", Description)) %>%
  select(type_eco_prefix, label)

lab_map <- setNames(prefix_labels$label, prefix_labels$type_eco_prefix)
# Default to prefix if missing in dictionary
missing <- setdiff(levels(dens_df$community), names(lab_map))
if (length(missing)) lab_map[missing] <- missing

# 7) Plot
p_dens <- ggplot(dens_df, aes(x = value)) +
  geom_density(color = "black", fill = "grey75", linewidth = 0.6, adjust = 1) +
  geom_vline(data = panel_q, aes(xintercept = q025), linetype = "dashed", linewidth = 0.5) +
  geom_vline(data = panel_q, aes(xintercept = q975), linetype = "dashed", linewidth = 0.5) +
  geom_vline(data = panel_q, aes(xintercept = med), linetype = "solid", linewidth = 0.6) +
  facet_grid(community ~ scenario, scales = "free_y", labeller = labeller(community = lab_map)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Predicted presence probability distributions by community and climate scenario",
    x = "Predicted presence probability",
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    strip.text.y = element_text(face = "bold", size = 10),
    strip.text.x = element_text(size = 10),
    panel.spacing.x = unit(8, "pt"),
    panel.spacing.y = unit(8, "pt"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9)
  )

print(p_dens)

# 8) Save figure
out_dir <- file.path(shiny_app_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = file.path(out_dir, paste0("density_by_community_scenarios_", paste(available_prefixes, collapse = "_"), ".png")),
  plot = p_dens,
  width = max(8, 3 + 2.5 * length(scenario_order)),
  height = max(6, 2 + 1.5 * length(available_prefixes)),
  dpi = 300
)

# ===========================
# Species-wise probability distributions across scenarios
# ===========================

# Choose a manageable set of species (top by prevalence and available in predictions)
first_scn_name <- names(all_predictions_summary)[1]
pred_species <- colnames(all_predictions_summary[[first_scn_name]]$mean)

sp_in_data <- intersect(colnames(Y_pa), pred_species)
sp_rank <- sort(colSums(Y_pa[, sp_in_data, drop = FALSE], na.rm = TRUE), decreasing = TRUE)
target_species <- names(sp_rank)[seq_len(min(12, length(sp_rank)))]

# Labels for facets and colors
sp_labels_tbl <- species_map %>%
  dplyr::filter(original_name_in_data %in% target_species) %>%
  transmute(col = original_name_in_data,
            label = paste0(pretty_name, " (", code, ")"))
label_map <- setNames(sp_labels_tbl$label, sp_labels_tbl$col)

# Scenario order (reuse if defined earlier; otherwise keep current ordering)
if (!exists("scenario_order") || !length(scenario_order)) {
  scenario_order <- names(all_predictions_summary)
}

# Build a long data.frame of probabilities by species and scenario
dens_sp <- list()
for (sp in target_species) {
  for (scn in scenario_order) {
    m <- all_predictions_summary[[scn]]$mean
    if (!sp %in% colnames(m)) next
    v <- as.numeric(m[, sp])
    v <- v[is.finite(v) & v >= 0 & v <= 1]
    if (!length(v)) next
    dens_sp[[length(dens_sp) + 1]] <- tibble::tibble(
      value = v,
      species = sp,
      scenario = scn
    )
  }
}
dens_sp <- dplyr::bind_rows(dens_sp)
if (nrow(dens_sp) == 0) stop("No data available for the selected species and scenarios.")

dens_sp <- dens_sp %>%
  dplyr::mutate(
    species_label = factor(label_map[species], levels = unique(label_map[target_species])),
    scenario = factor(scenario, levels = scenario_order)
  )

# Per-panel summary stats
panel_stats <- dens_sp %>%
  dplyr::group_by(species_label, scenario) %>%
  dplyr::summarise(
    lo = unname(stats::quantile(value, 0.025, na.rm = TRUE)),
    md = stats::median(value, na.rm = TRUE),
    hi = unname(stats::quantile(value, 0.975, na.rm = TRUE)),
    .groups = "drop"
  )

# Distinct colors per species (consistent across all its panels)
nsp <- nlevels(dens_sp$species_label)
col_pal <- grDevices::hcl.colors(nsp, palette = "Dark 3")
names(col_pal) <- levels(dens_sp$species_label)
fill_pal <- grDevices::adjustcolor(col_pal, alpha.f = 0.35)
names(fill_pal) <- names(col_pal)
# Split target_species into 3 groups for separate plots
n_species <- length(target_species)
group_size <- ceiling(n_species / 3)
species_groups <- split(target_species, ceiling(seq_along(target_species) / group_size))

# Function to create plot for a group
create_species_plot <- function(sp_group) {
  dens_sp_sub <- dens_sp %>% dplyr::filter(species %in% sp_group)
  panel_stats_sub <- panel_stats %>% dplyr::filter(species_label %in% label_map[sp_group])
  
  nsp_sub <- length(unique(dens_sp_sub$species_label))
  col_pal_sub <- grDevices::hcl.colors(nsp_sub, palette = "Dark 3")
  names(col_pal_sub) <- unique(dens_sp_sub$species_label)
  fill_pal_sub <- grDevices::adjustcolor(col_pal_sub, alpha.f = 0.35)
  names(fill_pal_sub) <- names(col_pal_sub)
  
  p <- ggplot2::ggplot(dens_sp_sub, ggplot2::aes(x = value, color = species_label, fill = species_label)) +
    ggplot2::geom_density(linewidth = 0.6, adjust = 1) +
    ggplot2::geom_vline(data = panel_stats_sub, ggplot2::aes(xintercept = lo, color = species_label),
                        linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
    ggplot2::geom_vline(data = panel_stats_sub, ggplot2::aes(xintercept = md, color = species_label),
                        linetype = "solid", linewidth = 0.6, inherit.aes = FALSE) +
    ggplot2::geom_vline(data = panel_stats_sub, ggplot2::aes(xintercept = hi, color = species_label),
                        linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
    ggplot2::facet_grid(rows = ggplot2::vars(scenario), cols = ggplot2::vars(species_label),
                        scales = "free_y") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::scale_color_manual(values = col_pal_sub, guide = "none") +
    ggplot2::scale_fill_manual(values = fill_pal_sub, guide = "none") +
    ggplot2::labs(
      title = "Predicted presence probability distributions by species and scenario",
      x = "Predicted presence probability",
      y = "Density"
    ) +
    ggplot2::theme_light(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      strip.text.y = ggplot2::element_text(face = "bold", size = 10),
      strip.text.x = ggplot2::element_text(size = 10),
      panel.spacing.x = grid::unit(8, "pt"),
      panel.spacing.y = grid::unit(8, "pt"),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 9)
    )
  return(p)
}

# Create the 3 plots
p1 <- create_species_plot(species_groups[[1]])
p2 <- create_species_plot(species_groups[[2]])
p3 <- create_species_plot(species_groups[[3]])

# Arrange them vertically (requires patchwork or gridExtra; assuming patchwork is installed)
# install.packages("patchwork") # Uncomment if patchwork is not installed
library(patchwork)
p_combined <- p1 / p2 / p3
print(p_combined)

# Save the combined plot
ggplot2::ggsave(
  filename = file.path(output_dir, "hmsc_figures", "density_by_species_scenarios_split.png"),
  plot = p_combined,
  width = 12,
  height = 18,
  dpi = 300
)

# Save
ggplot2::ggsave(
  filename = file.path(output_dir, "hmsc_figures", paste0(
    "density_by_species_scenarios_",
    paste0(unique(target_species), collapse = "_"), ".png"
  )),
  plot = p_species,
  width = max(8, 3 + 2.5 * length(levels(dens_sp$scenario))),
  height = max(6, 1.5 + 1.2 * nlevels(dens_sp$species_label)),
  dpi = 300
)







