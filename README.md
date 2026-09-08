# PEP_QC_HMSC_Analysis

Reproduction pipeline for *"Beyond Single Species: Using Community Suitability to
Inform Climate Adaptation in the Boreal-Temperate Ecotone"* (Czarnecki de Liz et al.).

Everything needed to reproduce the paper is in **`PEP_QC_HMSC_Analysis.Rmd`**, using
the plot database **`databases/undisturbed_plots_sf_env_clim_4th_inv.csv`**. Set
`BASE_DIR` at the top of the Rmd to this repository's local path and knit.

The companion Shiny application (TCAMS) built on top of these outputs lives in the
sibling repository **PEP_QC_TCAMS**.

## Future climate data (not included — regenerate before knitting)

Section 10 of the Rmd reads `climate/plots_coord_12_GCMsY.csv` (8-GCM-ensemble MAT/MAP
under 4 SSPs × 3 horizons, per plot), not included here. To regenerate it:

1. From `databases/undisturbed_plots_sf_env_clim_4th_inv.csv`, export one row per plot
   with columns `id1, id2, lat, long, elev` (`id1`/`id2` = `plot_id`).
2. Run this file through [ClimateNA](https://climatena.ca/) in batch mode, requesting
   the 8-GCM ensemble mean for SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5 × 2011-2040,
   2041-2070, 2071-2100 (13 runs, or one multi-scenario batch job).
3. Concatenate ClimateNA's outputs into a single CSV with columns `id1, GCM, MAT, MAP`
   (one row per plot × scenario; `GCM` = the scenario name, e.g.
   `8GCMs_ensemble_ssp370_2071-2100.gcm`) and place it at `climate/plots_coord_12_GCMsY.csv`.
