#######################################################################
# Genetic diversity analysis with wingen (algatr)
# https://thewanglab.github.io/algatr/articles/wingen_vignette.html
#######################################################################

library(algatr)
library(wingen)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

## ------------------------------------------------
## 0. Paths (match earlier pipeline scripts)
## ------------------------------------------------
data_dir <- path.expand("~/Desktop/LandscapeGenomicsPipeline/")
out_dir <- file.path(data_dir, "outputs", "02-genetic-diversity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# All US states
us_states <- rnaturalearth::ne_states(
  country     = "United States of America",
  returnclass = "sf"
)
# Keep only East Coast states (adjust list as you like)
east_coast_abbrev <- c(
  "ME","NH","MA","RI","CT","NY","NJ",
  "DE","MD","VA","NC","SC","GA","FL",
  "AL", "LA", "TN", "MS", "AR", "KY", 
  "OH", "MI", "IN", "IL", "MO", "IA",
  "PA", "TX", "OK", "KS"
)
us_states_ec <- us_states %>%
  filter(postal %in% east_coast_abbrev)

## ------------------------------------------------
## 0. Input checks and CRS helpers
## ------------------------------------------------

required_objs <- c("wclim", "coord", "vcf_pruned")
missing_objs <- required_objs[!vapply(required_objs, exists, logical(1), inherits = TRUE)]
if (length(missing_objs) > 0) {
  stop(
    "Missing required object(s): ", paste(missing_objs, collapse = ", "), "\n",
    "Run `00-Preprocessing.R` first in the same R session, or load saved objects before running this script."
  )
}

if (!all(c("x", "y") %in% names(coord))) {
  stop("`coord` must contain numeric columns named `x` (lon) and `y` (lat).")
}

# Ensure sf points have an explicit CRS (assume WGS84 lon/lat unless you know otherwise)
coord_sf <- st_as_sf(coord, coords = c("x", "y"), crs = 4326, remove = FALSE)

# Pick a sensible projected CRS for window operations (local UTM based on centroid)
pick_utm_epsg <- function(lon, lat) {
  zone <- floor((lon + 180) / 6) + 1
  if (lat >= 0) 32600 + zone else 32700 + zone
}

center_lon <- mean(coord$x, na.rm = TRUE)
center_lat <- mean(coord$y, na.rm = TRUE)
proj_epsg <- pick_utm_epsg(center_lon, center_lat)
proj_crs <- st_crs(proj_epsg)

# Asimina range
asimina_range_path_gpkg <- file.path(data_dir, "outputs", "01-rangemap", "Asimina_triloba_range_ecol.gpkg")
asimina_range_path_shp <- file.path(data_dir, "outputs", "01-rangemap", "Asimina_triloba_range_ecol.shp")
if (file.exists(asimina_range_path_gpkg)) {
  asimina_range <- st_read(asimina_range_path_gpkg, layer = "range", quiet = TRUE)
} else if (file.exists(asimina_range_path_shp)) {
  asimina_range <- st_read(asimina_range_path_shp, quiet = TRUE)
} else {
  stop(
    "Asimina range file not found.\nExpected one of:\n  ",
    asimina_range_path_gpkg, "\n  ", asimina_range_path_shp,
    "\nRun `01-Atriloba_rangemap.R` first."
  )
}
asimina_range_proj <- st_transform(asimina_range, crs = proj_crs)

## ------------------------------------------------
## 1. Environmental distances
## ------------------------------------------------

# Extract raster values at sample locations as a matrix/data.frame
# (wclim is a terra SpatRaster from earlier)
coord_vect_ll <- terra::vect(coord, geom = c("x", "y"), crs = "EPSG:4326")
env_mat <- terra::extract(wclim, coord_vect_ll)[, -1, drop = FALSE]

# Calculate environmental distances
env_dist_obj <- env_dist(env_mat)

# Example: inspect distance distribution for BIO1
png(filename = file.path(out_dir, "env_dist_bio1.png"), width = 1600, height = 1200, res = 200)
plot(env_dist_obj$bio1, main = "Environmental distance for BIO1")
dev.off()

## ------------------------------------------------
## 2. Project coordinates and environmental layer
## ------------------------------------------------

# Project coordinates to a local projected CRS (UTM) so `res` is in meters
coords_proj <- st_transform(coord_sf, crs = proj_crs)
states_proj <- st_transform(us_states_ec, crs = proj_crs)
png(filename = file.path(out_dir, paste0("projected_sample_locations_EPSG", proj_epsg, ".png")),
    width = 1600, height = 1200, res = 200)
plot(coords_proj["geometry"], main = paste0("Projected sample locations (EPSG:", proj_epsg, ")"))
dev.off()

# Pick a base environmental layer for plotting (BIO6) and aggregate for speed
# Convert to terra SpatRaster if needed
bio15 <- wclim[["bio1"]]

# Aggregate to coarser resolution (factor = 5)
envlayer <- terra::aggregate(bio1, fact = 5)


# Reproject envlayer to match projected coordinates
envlayer <- terra::project(envlayer, proj_crs$wkt)
png(filename = file.path(out_dir, paste0("bio1_aggregated_projected_EPSG", proj_epsg, ".png")),
    width = 1600, height = 1200, res = 200)
plot(envlayer, main = "Aggregated & projected BIO1")
plot(st_geometry(states_proj), add = TRUE, border = "black", lwd = 0.5)
dev.off()

## ------------------------------------------------
## 3. Moving-window (wingen) genetic diversity
## ------------------------------------------------

# Generate raster for sliding window (spacing & buffer in map units of coords_proj)
png(
  filename = file.path(out_dir, "coords_to_raster_window_grid.png"),
  width = 1600,
  height = 1200,
  res = 200
)
liz_lyr <- coords_to_raster(
  coords_proj,
  res    = 50000,  # cell size (e.g. 50 km)
  buffer = 5,
  plot   = TRUE
)
dev.off()

# Preview window and sample counts
sample_count_prev <- preview_gd(liz_lyr, coords_proj, wdim = 3, fact = 0)
sample_count_prev0 <- sample_count_prev
sample_count_prev0[sample_count_prev0 == 0] <- NA

p_prev <- ggplot_count(sample_count_prev0) +
  ggtitle("Preview sample count") +
  # 1) Asimina triloba native range (soft grey polygon)
  geom_sf(
    data  = asimina_range_proj,
    fill  = "grey90",
    alpha = 0.4,
    color = "grey60",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(na.value = "grey") +
  geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.3)
print(p_prev)
ggsave(
  filename = file.path(out_dir, "wingen_preview_sample_count.png"),
  plot = p_prev,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Run moving-window genetic diversity using the pruned VCF
# (use vcf_pruned from your earlier LD-pruning step in preprocessing
wgd <- window_gd(
  gen    = vcf_pruned,
  coords = coords_proj,
  lyr    = liz_lyr,
  stat   = "pi",
  wdim   = 3,
  fact   = 0
)

## ------------------------------------------------
## 4. Visualizing moving-window results
## ------------------------------------------------
## ------------------------------------------------
## 4. Visualizing moving-window results (improved π map)
## ------------------------------------------------

# Projected sample points for overlay
pts_proj <- st_transform(coord_sf, crs = proj_crs)

# Get π range and nice breaks
pi_rng <- terra::global(wgd[["pi"]], range, na.rm = TRUE)
pi_min <- as.numeric(pi_rng[1, 1])
pi_max <- as.numeric(pi_rng[1, 2])
pi_brks <- pretty(c(pi_min, pi_max), n = 5)

# Crop view to Asimina range bbox (slightly buffered)
bb <- st_bbox(asimina_range_proj)
x_pad <- (bb["xmax"] - bb["xmin"]) * 0.05
y_pad <- (bb["ymax"] - bb["ymin"]) * 0.05

p_pi <- ggplot_gd(wgd) +
  ggtitle("Moving window π") +
  # 1) Asimina triloba native range (soft grey polygon)
  geom_sf(
    data  = asimina_range_proj,
    fill  = "grey90",
    alpha = 0.4,
    color = "grey60",
    linewidth = 0.4
  ) +
  # 2) State borders (thin outlines)
  geom_sf(
    data  = states_proj,
    fill  = NA,
    color = "grey40",
    linewidth = 0.25
  ) +
  # 3) Sample points
  geom_sf(
    data  = pts_proj,
    color = "black",
    fill  = "white",
    shape = 21,
    size  = 1.2,
    stroke = 0.25
  ) +
  # 4) π color scale
  scale_fill_viridis_c(
    option = "magma",
    na.value = NA,
    limits = c(pi_min, pi_max),
    breaks = pi_brks
  ) +
  # 5) Zoom to Asimina range with a small buffer
  coord_sf(
    xlim = c(bb["xmin"] - x_pad, bb["xmax"] + x_pad),
    ylim = c(bb["ymin"] - y_pad, bb["ymax"] + y_pad)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_pi)

ggsave(
  filename = file.path(out_dir, "wingen_moving_window_pi.png"),
  plot = p_pi,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)

# Sample count map: set unsampled (0) cells to NA so they can be colored grey
# 0-count -> NA so unsampled is grey
wgd_count0 <- wgd
wgd_count0[["sample_count"]][wgd_count0[["sample_count"]] == 0] <- NA

# projected sample points for overlay
pts_proj <- st_transform(coord_sf, crs = proj_crs)

# compute good legend breaks
mx <- terra::global(wgd[["sample_count"]], "max", na.rm = TRUE)[1, 1]
mx <- if (is.na(mx)) 1 else as.numeric(mx)
brks <- pretty(c(1, mx), n = 5)

p_count <- ggplot_count(wgd_count0) +
  ggtitle("Sample count (moving window)") +
  # 1) Asimina triloba native range (soft grey polygon)
  geom_sf(
    data  = asimina_range_proj,
    fill  = "grey90",
    alpha = 0.4,
    color = "grey60",
    linewidth = 0.4
  ) +
  geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.25) +
  geom_sf(data = pts_proj, color = "black", fill = "white", shape = 21, size = 1.2, stroke = 0.25) +
  scale_fill_viridis_c(
    option = "cividis",     # easier for colorblind + print
    na.value = "grey85",
    breaks = brks,
    limits = c(0, mx)
  ) +
  coord_sf(
    xlim = st_bbox(states_proj)[c("xmin", "xmax")],
    ylim = st_bbox(states_proj)[c("ymin", "ymax")]
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "wingen_moving_window_sample_count.png"),
  plot = p_count, width = 8, height = 6, units = "in", dpi = 300
)

## ------------------------------------------------
## 5. Kriging and masking
## ------------------------------------------------

# Kriging of genetic diversity metrics
kgd <- krig_gd(
  wgd,
  index      = 1:2,   # e.g. first two summary metrics in wgd
  disagg_grd = 5
)

summary(kgd)

# Kriged pi map
p_kpi <- ggplot_gd(kgd) + ggtitle("Kriged pi")
print(p_kpi)
ggsave(
  filename = file.path(out_dir, "wingen_kriged_pi.png"),
  plot = p_kpi,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Kriged sample count map with readable background
p_kcount <- ggplot_count(kgd) +
  ggtitle("Kriged sample counts") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )
print(p_kcount)
ggsave(
  filename = file.path(out_dir, "wingen_kriged_sample_counts.png"),
  plot = p_kcount,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Mask kriged results by sample count
mgd_1 <- mask_gd(kgd, kgd[["sample_count"]], minval = 1)
mgd_2 <- mask_gd(kgd, kgd[["sample_count"]], minval = 2)
mgd_3 <- mask_gd(kgd, kgd[["sample_count"]], minval = 3)

p_masked <- ggplot_gd(mgd_2) + ggtitle("Kriged & masked pi (≥2 samples)") +
  # 1) Asimina triloba native range (soft grey polygon)
  geom_sf(
    data  = asimina_range_proj,
    fill  = "grey90",
    alpha = 0.4,
    color = "grey60",
    linewidth = 0.4
  ) +
  # 2) State borders (thin outlines)
  geom_sf(
    data  = states_proj,
    fill  = NA,
    color = "grey40",
    linewidth = 0.25
  ) +
  # 3) Sample points
  geom_sf(
    data  = pts_proj,
    color = "black",
    fill  = "white",
    shape = 21,
    size  = 1.2,
    stroke = 0.25
  ) +
  scale_fill_viridis_c(
    option = "magma",
    na.value = NA
  )
print(p_masked)
ggsave(
  filename = file.path(out_dir, "wingen_kriged_masked_pi_min2samples.png"),
  plot = p_masked,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
