# create_spatial_studyDesign_for_Hmsc.R
# Utility to build a studyDesign data.frame and HmscRandomLevel objects
# for all levels listed in ranLevels, using available coordinates.
#
# Usage:
#   res <- create_spatial_studyDesign_spatial(plots_sf_env_clim, ranLevels = c("Site","Plot"))
#   studyDesign <- res$studyDesign            # data.frame with columns named as ranLevels, rownames = plot_id
#   ranLevels_hmsc <- res$ranLevels_hmsc      # named list of HmscRandomLevel objects
#
# Requirements: plots_sf_env_clim must contain a unique plot identifier column (default "plot_id")
# and two coordinate columns (one of longitude/latitude or x/y).

create_spatial_studyDesign_spatial <- function(data,
                                              ranLevels,
                                              plot_id_col = "plot_id",
                                              coord_cols = c("longitude", "latitude", "x", "y")) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame or tibble.")
  if (!plot_id_col %in% names(data)) stop("plot_id column not found in data.")
  # find coordinate pair
  present_coords <- coord_cols[coord_cols %in% names(data)]
  if (length(present_coords) < 2) stop("Need at least two coordinate columns (e.g. longitude+latitude or x+y).")
  coord_x <- present_coords[1]; coord_y <- present_coords[2]
  # unique plots with coordinates
  coords_df <- unique(data[, c(plot_id_col, coord_x, coord_y), drop = FALSE])
  # ensure one row per plot_id
  if (any(duplicated(coords_df[[plot_id_col]]))) {
    # collapse duplicates by taking mean of coordinates
    coords_df <- aggregate(coords_df[, c(coord_x, coord_y)], by = list(coords_df[[plot_id_col]]), FUN = function(z) mean(as.numeric(z), na.rm = TRUE))
    names(coords_df)[1] <- plot_id_col
  }
  plot_ids <- as.character(coords_df[[plot_id_col]])
  # build studyDesign: for each level create a factor column with length = number of plots
  sd <- data.frame(row_id = plot_ids, stringsAsFactors = FALSE)
  for (lvl in ranLevels) {
    if (lvl %in% names(data)) {
      # use value at plot level (match by plot_id)
      lvl_vals <- data[[lvl]][match(plot_ids, as.character(data[[plot_id_col]]))]
      sd[[lvl]] <- factor(lvl_vals)
    } else {
      # fallback: use unique plot_id (each plot is its own unit for that level)
      sd[[lvl]] <- factor(plot_ids)
    }
  }
  rownames(sd) <- sd$row_id
  sd$row_id <- NULL
  # create HmscRandomLevel objects with sData aggregated per level using mean coordinates
  if (!requireNamespace("Hmsc", quietly = TRUE)) stop("Package 'Hmsc' is required.")
  ranLevels_hmsc <- list()
  for (lvl in ranLevels) {
    # map each plot to its level value
    if (lvl %in% names(data)) {
      level_map <- data[[lvl]][match(plot_ids, as.character(data[[plot_id_col]]))]
    } else {
      level_map <- plot_ids
    }
    # aggregate coordinates by level (mean)
    agg <- aggregate(coords_df[, c(coord_x, coord_y)], by = list(level = level_map), FUN = function(z) mean(as.numeric(z), na.rm = TRUE))
    rownames(agg) <- as.character(agg$level)
    sMat <- as.matrix(agg[, c(coord_x, coord_y), drop = FALSE])
    colnames(sMat) <- c("x", "y")
    # create HmscRandomLevel with sData; HmscRandomLevel will accept a matrix with rownames = level names
    ranLevels_hmsc[[lvl]] <- Hmsc::HmscRandomLevel(sData = sMat)
  }
  # Return studyDesign and ranLevels list suitable for passing to Hmsc()
  return(list(
    studyDesign = as.data.frame(sd),
    ranLevels_hmsc = ranLevels_hmsc
  ))
}

# ---- Example ----
# If you have plots_sf_env_clim loaded, and want levels Site and Plot:
res <- create_spatial_studyDesign_spatial(plots_sf_env_clim, ranLevels = c("Site","Plot"))
studyDesign <- res$studyDesign
ranLevels <- res$ranLevels_hmsc
