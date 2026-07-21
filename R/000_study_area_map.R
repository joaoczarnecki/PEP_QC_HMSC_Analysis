# ---
# title: "Study Area Map - Permanent Sample Plots (PEP) - Quebec 4th Inventory"
# author: "Liz, JPC."
# date: "2025-08-14"
# ---

#===============================================================================
# 0. PACKAGES
#===============================================================================
.libPaths("H:/R_lib")

pkgs <- c("tidyverse", "sf", "ggplot2", "rnaturalearth", "rnaturalearthdata",
          "rnaturalearthhires", "ggspatial", "patchwork", "scales", "viridis")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(tidyverse); library(sf); library(ggplot2)
library(rnaturalearth); library(rnaturalearthdata); library(rnaturalearthhires)
library(ggspatial); library(patchwork); library(scales); library(viridis)

#===============================================================================
# 1. DATA LOADING
#===============================================================================
plots_df <- read_csv(
  "g:/Thesis/3rdChapter/PEP_QC/databases/undisturbed_plots_sf_env_clim_4th_inv.csv",
  show_col_types = FALSE
)

plots_sf <- plots_df %>%
  dplyr::select(plot_id, latitude, longitude, MAT, MAP) %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

cat("Plots:", nrow(plots_sf), "\n")

#===============================================================================
# 2. SPATIAL REFERENCE LAYERS
#===============================================================================
LAMBERT_QC <- 32198   # NAD83 / Quebec Lambert

# Raw boundaries (WGS84)
provinces <- ne_states(country = "canada",
                        returnclass = "sf", scale = "large") %>% st_make_valid()
us_states  <- ne_states(country = "united states of america",
                         returnclass = "sf", scale = "large") %>% st_make_valid()
world      <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid()

quebec_sf  <- provinces %>% filter(name == "Quebec") %>% st_make_valid()

context_sf <- bind_rows(
  provinces %>% filter(name %in% c("Ontario", "New Brunswick",
                                    "Nova Scotia", "Newfoundland and Labrador",
                                    "Prince Edward Island")),
  us_states  %>% filter(name %in% c("Maine", "Vermont", "New Hampshire",
                                     "New York", "Massachusetts"))
) %>% st_make_valid()

# Study extent (WGS84 → Lambert)
bbox_wgs <- st_bbox(c(xmin = -80, xmax = -57, ymin = 44.5, ymax = 53),
                    crs = st_crs(4326))
bbox_lam  <- st_bbox(st_transform(st_as_sfc(bbox_wgs), LAMBERT_QC))

# Clip polygon with generous buffer so polygon fills don't get truncated
clip_lam  <- st_as_sfc(bbox_lam) %>% st_buffer(500000) %>% st_make_valid()

# Helper: safe pre-project + intersect (more robust than st_crop for complex polys)
safe_project <- function(x, crs, clip) {
  x %>%
    st_transform(crs) %>%
    st_make_valid() %>%
    st_intersection(clip) %>%
    st_make_valid()
}

world_p   <- safe_project(world,      LAMBERT_QC, clip_lam)
context_p <- safe_project(context_sf, LAMBERT_QC, clip_lam)
quebec_p  <- safe_project(quebec_sf,  LAMBERT_QC, clip_lam)
plots_p   <- st_transform(plots_sf,   LAMBERT_QC)

# 52nd parallel in Lambert
parallel_52 <- st_linestring(matrix(c(-82, 52, -56, 52), ncol = 2, byrow = TRUE)) %>%
  st_sfc(crs = 4326) %>% st_sf() %>% st_transform(LAMBERT_QC)

label_52N <- st_coordinates(
  st_transform(st_sfc(st_point(c(-58, 52.5)), crs = 4326), LAMBERT_QC))

#===============================================================================
# 3. COLOURS & THEME
#===============================================================================
clr_ocean      <- "#B8D4E8"
clr_land       <- "#EAE8E2"
clr_quebec     <- "#C8E6C9"
clr_neighbours <- "#D8D3CB"
clr_border     <- "#888888"
clr_parallel   <- "#CC3333"

theme_map <- theme_void(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = clr_ocean, colour = NA),
    panel.border     = element_rect(fill = NA, colour = "grey35", linewidth = 0.7),
    panel.grid.major = element_line(colour = "white", linewidth = 0.25),
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5,
                                     margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 9,  hjust = 0.5, colour = "grey30",
                                     margin = margin(b = 6)),
    plot.margin      = margin(6, 6, 6, 6)
  )

#===============================================================================
# 4. MAIN MAP
#===============================================================================
main_map <- ggplot() +

  # Layer order: land → neighbours → Quebec → parallel → plots
  geom_sf(data = world_p,   fill = clr_land,       colour = "grey65",  linewidth = 0.10) +
  geom_sf(data = context_p, fill = clr_neighbours,  colour = clr_border, linewidth = 0.35) +
  geom_sf(data = quebec_p,  fill = clr_quebec,      colour = "grey20",  linewidth = 0.70) +

  geom_sf(data = parallel_52,
          colour = clr_parallel, linewidth = 0.5, linetype = "dashed") +
  annotate("text",
           x = label_52N[1], y = label_52N[2],
           label = "52°N", colour = clr_parallel,
           size = 3, fontface = "italic", hjust = 0) +

  geom_sf(data = plots_p,
          aes(colour = MAT), size = 0.70, alpha = 0.85, shape = 16) +

  scale_colour_viridis_c(
    name   = "MAT (°C)",
    option = "plasma",
    breaks = pretty_breaks(n = 5),
    guide  = guide_colourbar(
      barwidth = 0.8, barheight = 6,
      title.position = "top", title.hjust = 0.5
    )
  ) +

  coord_sf(
    xlim   = c(bbox_lam["xmin"], bbox_lam["xmax"]),
    ylim   = c(bbox_lam["ymin"], bbox_lam["ymax"]),
    expand = FALSE,
    crs    = LAMBERT_QC,
    datum  = st_crs(4326)
  ) +

  annotation_scale(
    location = "bl", width_hint = 0.22, text_cex = 0.7,
    bar_cols = c("grey20", "white"), line_col = "grey20", text_col = "grey20"
  ) +

  annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.55, "in"), pad_y = unit(0.15, "in"),
    style = north_arrow_fancy_orienteering(text_size = 8,
                                            fill = c("grey20", "white"))
  ) +

  labs(
    title    = "Permanent Sample Plots – 4th Forest Inventory",
    subtitle = paste0("Quebec, Canada  |  n = ",
                      format(nrow(plots_sf), big.mark = ","),
                      " plots  |  below 52°N latitude")
  ) +
  theme_map

#===============================================================================
# 5. INSET MAP
#===============================================================================
world_inset  <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid()
canada_inset <- ne_countries(country = "canada",
                              scale = "medium", returnclass = "sf") %>% st_make_valid()
usa_inset    <- ne_countries(country = "united states of america",
                              scale = "medium", returnclass = "sf") %>% st_make_valid()
study_box    <- st_as_sfc(bbox_wgs)

inset_map <- ggplot() +
  geom_sf(data = world_inset,  fill = "#E8E6DF", colour = "grey65", linewidth = 0.15) +
  geom_sf(data = canada_inset, fill = "#D8D3CB", colour = "grey50", linewidth = 0.25) +
  geom_sf(data = usa_inset,    fill = "#D8D3CB", colour = "grey50", linewidth = 0.25) +
  geom_sf(data = quebec_sf,    fill = "#A8D5B5", colour = "grey20", linewidth = 0.40) +
  geom_sf(data = study_box,    fill = NA, colour = clr_parallel,
          linewidth = 0.85, linetype = "solid") +
  coord_sf(xlim = c(-145, -50), ylim = c(40, 75), expand = FALSE) +
  theme_void(base_size = 8) +
  theme(
    panel.background = element_rect(fill = clr_ocean,  colour = "grey35", linewidth = 0.6),
    plot.background  = element_rect(fill = "white",    colour = NA),
    plot.margin      = margin(1, 1, 1, 1)
  )

#===============================================================================
# 6. COMPOSE & SAVE
#===============================================================================
final_fig <- main_map +
  inset_element(inset_map,
                left = 0.61, bottom = 0.61, right = 0.99, top = 0.99)

output_dir <- "g:/Thesis/3rdChapter/PEP_QC/results/maps"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(output_dir, "study_area_map.png"),
       final_fig, width = 10, height = 8, dpi = 300, bg = "white")

ggsave(file.path(output_dir, "study_area_map.svg"),
       final_fig, width = 10, height = 8, device = "svg", bg = "white")

cat("Saved to:", output_dir, "\n")
