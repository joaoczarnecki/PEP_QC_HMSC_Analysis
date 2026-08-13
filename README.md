# PEP_QC_HMSC_Analysis

Analysis/reproduction pipeline for the manuscript *"Beyond Single Species:
Using Community Suitability to Inform Climate Adaptation in the
Boreal-Temperate Ecotone"* (Czarnecki de Liz et al.). Fits a Bayesian Joint
Species Distribution Model (JSDM) with the
[Hmsc](https://github.com/hmsc-r/HMSC) R package to presence/absence data
for 13 dominant tree species from Quebec's permanent forest inventory
(5,572 undisturbed plots), then derives Community Suitability Indices (CSI)
and projects them under four SSP climate scenarios.

The companion Shiny application built on top of these outputs lives in the
sibling repository **PEP_QC_TCAMS**.

## Pipeline (`R/`)

Run in numeric order. Each script was copied here, under a clean name, from
the working (unversioned) `PEP_QC/script/` folder — see "Provenance &
assumptions" below for how the canonical version of each step was chosen out
of several parallel/duplicate copies found during this reorganization.

| Script | Purpose |
|---|---|
| `000_study_area_map.R` | Study-area map figure |
| `001_fit_hmsc_model.R` | Builds `Y`/`XData`/`studyDesign`, fits the probit Hmsc model (4 chains, thin=100, samples=1000, transient=50000, adaptNf=40000), runs convergence diagnostics (ESS/PSRF on `Beta`), model fit (Tjur R²/AUC/RMSE, in-sample), variance partitioning, species associations, and environmental gradient plots |
| `002_predict_future_scenarios.R` | Loads the fitted model, builds per-scenario `XData` (Current + 12 GCM-ensemble × SSP × horizon combinations, varying only MAT/MAP — all other covariates held at their current values), predicts occurrence probability per posterior draw, and packages everything the Shiny app needs (`species_map`, `TYPE_ECO_list`, prediction summaries) |
| `003_compute_csi_from_draws.R` | Recomputes bCSI/wCSI/gCSI **from the full posterior draws array** (not from already-aggregated summaries — see caveat below), per site and per predefined community, for all 13 scenarios, with 95% credible intervals |
| `003a_plot_csi_trajectories.R` | Figures for CSI trajectories across scenarios/horizons |
| `004_prepare_app_draws_subsample.R` | Subsamples the (very large) full posterior draw arrays down to a size the Shiny app can ship/load quickly |
| `utils_spatial_study_design_UNUSED.R` | Utility to build a **spatially explicit** random level (`HmscRandomLevel(sData = coords)`), mirroring the approach used for the Route random effect in Tikhonov et al.'s bird case study. **Not currently called anywhere in the pipeline** — the fitted model uses a plain, non-spatial site-level random effect (`HmscRandomLevel(units = ...)`). Kept for reference / future work; see caveat below. |

## Known methodological caveats found during this review

Documented here so they can be addressed in the manuscript text and/or a
future re-fit, rather than silently living only in code:

1. **Convergence diagnostics (ESS/PSRF) are computed only for `Beta`**
   (species environmental responses), not for `Omega` (residual
   associations) or the random-effect hyperparameters. If the manuscript
   states convergence was checked "overall", either narrow that claim or
   extend `001_fit_hmsc_model.R` to also check `Omega`/`V`/`Alpha`.
2. **`evaluateModelFit`/`computePredictedValues` are in-sample** (no
   `createPartition()`/cross-validation is used anywhere in this pipeline).
   The reported AUC/Tjur R² are explanatory, not predictive/validated —
   make sure the manuscript's wording matches.
3. **Future predictions use `predictEtaMean = TRUE`**, i.e. the random
   effect's posterior *mean* is used rather than propagating the full
   per-draw random-effect realization. This is a legitimate simplification
   but narrows/simplifies the uncertainty propagated into the CSI credible
   intervals compared to a "full" predictive draw — worth an explicit
   sentence in the methods if not already there.
4. **Two non-equivalent gCSI pipelines existed** in the original `script/`
   folder: one applies the geometric mean to already-aggregated
   mean/lower/upper summaries (`003_gCSI_trajectory_analysis.R`, *not*
   copied into this repo), the other (kept here as
   `003_compute_csi_from_draws.R`) recomputes it from the raw posterior
   draws, which is the mathematically correct way to propagate uncertainty
   through a non-linear (geometric-mean) statistic. **This repo assumes the
   draws-based version is the one whose numbers were used in the
   manuscript** — please confirm.
5. **The site-level random effect is non-spatial** (`HmscRandomLevel(units =
   ...)`), even though a ready-made utility for a spatial random level
   (`utils_spatial_study_design_UNUSED.R`) and a fully commented-out block
   inside `001_fit_hmsc_model.R` show a spatial version was drafted and
   abandoned. If spatial autocorrelation among plots is a plausible concern,
   this is worth flagging as a limitation (or revisiting) rather than
   leaving as silent dead code.
6. **Two harmless debug residues were removed** from
   `001_fit_hmsc_model.R` during this reorganization: a bare `m1` statement
   left over mid-script (before `m1` existed — would raise an "object not
   found" error on a clean re-run) and a bare comparison expression with no
   effect (`unique(studyDesign$site) == studyDesign$site`). Neither affected
   the already-fitted model saved to disk; they only blocked a from-scratch
   re-execution of the script.
7. **Community/vegetation-type trajectories are a province-wide average of
   the community's defining species set, not an average restricted to plots
   actually classified as that ecological type in the forest inventory.**
   Traced in `003_compute_csi_from_draws.R`: `process_scenario()` loads the
   full posterior prediction array for **all** ~5,572 modelled plots
   (`draws_array <- coerce_pred_obj_to_draws(preds)`) and, per community,
   only subsets the **species** dimension to that community's member species
   (`species_vec <- intersect(comm$species, available_species)`;
   `draws_sub <- draws_array[, , species_idx, drop = FALSE]`) — the **site**
   dimension is never filtered by `plot_id`/inventory-classified `type_eco`
   anywhere in the pipeline (confirmed in `compute_indices_from_draws_legacy()`
   and `compute_gCSI_matrix()`, both of which iterate over every site in the
   array). `summarise_community_indices_from_draws()` /
   `summarise_community_gCSI_matrix()` then average across *all* of those
   sites. In other words, a community's gCSI/bCSI/wCSI trajectory answers
   "how suitable is the *whole study area*, on average, for this species
   assemblage" — not "how suitable are the areas *currently occupied* by
   this vegetation type." This is consistent with how the manuscript defines
   CSI (a suitability index computable at any location for a target
   community), but is worth stating explicitly in the methods/discussion, since
   a reader could otherwise assume "FE3 trajectory" means plots inventoried
   as FE3 specifically.

## Provenance & assumptions made during this reorganization

The original `PEP_QC/script/` folder contained many parallel/dated copies of
each pipeline step (no git history to disambiguate them — see below). The
files copied into `R/` here were chosen as follows:

- **Fitting script**: `001_JSDM_PEP_QC_Modeling_20250814.R` (most recently
  modified, numbered/pipeline-style name) over ~8 older/alternative fitting
  scripts (`HMSC_PEP_QC_20251016.r`, `JSDM_PEP_QC_Refactored.R`,
  `JSDM_PEP_QC_DESIGNECO_20250814.R`, abundance-model variants, etc.), which
  were left untouched in the original workspace rather than duplicated here.
- **Prediction script**: `002_JSDM_PEP_TCAMS_predictions.r` (the only one
  matching the numbered pipeline and referenced by the app).
- **CSI script**: `003_gCSI_trajectory_analysis_draws.R`, over its
  near-duplicate `003_gCSI_trajectory_analysis copy.R` (byte-identical to
  `003_gCSI_trajectory_analysis.R`), the intermediate
  `003_JSDM_PEP_communities_trajectories_lat_shifts.r`, and the experimental
  `CSI_test.r` (a `doParallel` variant) — see caveat #4 above for *why* the
  draws-based version specifically was chosen, not just which file is newest.
- **Fitted model object**: the deployed Shiny app snapshot loads predictions
  already computed from a model, but there are **five** candidate fitted
  model files in `PEP_QC/results/undisturbed/hmsc/`
  (`hmsc_presence_model.Rdata`, two intermediate-MCMC-length versions dated
  14 Oct, one dated 19 Oct, one dated 22 Oct 2025). Based on file
  modification dates, **`hmsc_presence_model22oct25.Rda` is assumed to be
  the final one** — this is an assumption made during this review, not
  something confirmed against the manuscript's reported AUC=0.93/Tjur
  R²=0.37. Please verify before treating this as ground truth.

### Large data & model artefacts (not in this repo)

None of the following are tracked in this git repository — they are far too
large for a normal GitHub remote (plain git, no LFS configured) and remain
in the original workspace:

| Artefact | Approx. size | Where it lives |
|---|---|---|
| Fitted Hmsc model objects (5 candidates) | 0.8–3.9 GB each | `PEP_QC/results/undisturbed/hmsc/` |
| Full-resolution per-scenario posterior prediction arrays (13 scenarios) | ~2.2–2.3 GB each (~29 GB total) | `PEP_QC/results/undisturbed/` (root) and `.../preds/` (appear duplicated 1:1 — consolidate to one location before any further work) |
| Subsampled posterior draws, K≈200 per scenario | ~215 MB each (~2.9 GB total) | `PEP_QC/posterior_draws_200/` |
| Raw forest-inventory geopackage | ~16.6 GB | `PEP_QC/PEP.gpkg` |
| Per-site CSI draws (fine-grained) | up to 460 MB per file | `PEP_QC/results/undisturbed/outputs/output400L/site_indices_*.csv`, `gCSI_wCSI_bCSI_by_site_draws_summary.csv` |

Only the lightweight, already-aggregated summary CSVs and figures needed to
reproduce the manuscript's tables/plots were copied into `results/` and
`figures/` here.

### Other historical/duplicate content intentionally left out

`PEP_QC/script/presence_approach/` (a second, fully independent
implementation of fit → analysis → Shiny app) and `PEP_QC/old/` (six-plus
archived app/model versions, including the only on-disk evidence of a real
shinyapps.io deployment history) were **not** migrated here. They still
exist, untouched, in the original workspace if anything needs to be
recovered from them.

## Setup

Git was not available on the machine used to prepare this folder. See
`SETUP_GIT.md` for the commands to run once Git is installed.
