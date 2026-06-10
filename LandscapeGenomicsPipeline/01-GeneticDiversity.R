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

._lgpr <- Sys.getenv("LGP_PROJECT_ROOT", "~/Desktop/LandscapeGenomicsPipeline")
options(lgp.project_root = sub("/+$", "", path.expand(getOption("lgp.project_root", ._lgpr))))
suppressPackageStartupMessages(base::source(
  base::file.path(getOption("lgp.project_root"), "Scripts", "lgp_pipeline_cache.R"),
  encoding = "UTF-8"
))
rm(._lgpr)

## ------------------------------------------------
## 0. Paths (match earlier pipeline scripts)
## ------------------------------------------------
# Expects in InputData/: Pawpaw_range_Little_1977.zip (shapefile in zip),
#   kml_13792.kmz (Eastern Continental Divide). CRS is read from files; if
#   missing, WGS84 is assumed (see ensure_lonlat_crs()).
data_dir <- lgp_project_root()
out_dir <- lgp_outputs_step_dir("01-genetic-diversity")

preproc_dir <- lgp_preprocess_dir()
load_preproc_if_needed <- function(name, rds) {
  if (!exists(name, inherits = TRUE)) {
    p <- file.path(preproc_dir, rds)
    if (!file.exists(p)) {
      stop("Missing `", name, "` — run Scripts/00-Preprocessing.R or load RDS (expected ", p, ")")
    }
    assign(name, lgp_read_rds(p), envir = .GlobalEnv)
  }
}
load_preproc_if_needed("coord", "coord.rds")
load_preproc_if_needed("wclim", "selected_wclim.rds") # script expects name `wclim` historically
load_preproc_if_needed("vcf_pruned", "vcf_pruned.rds")

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
# Name must not be `coord_sf` — that masks ggplot2::coord_sf() and breaks coord/layers.
coord_pts_sf <- st_as_sf(coord, coords = c("x", "y"), crs = 4326, remove = FALSE)

# Pick a sensible projected CRS for window operations (local UTM based on centroid)
pick_utm_epsg <- function(lon, lat) {
  zone <- floor((lon + 180) / 6) + 1
  if (lat >= 0) 32600 + zone else 32700 + zone
}

center_lon <- mean(coord$x, na.rm = TRUE)
center_lat <- mean(coord$y, na.rm = TRUE)
proj_epsg <- pick_utm_epsg(center_lon, center_lat)
proj_crs <- st_crs(proj_epsg)

# Map limits: union bboxes so overlays (e.g. ECD line) are not clipped by state extent alone
st_bbox_union <- function(...) {
  objs <- list(...)
  bb <- sf::st_bbox(objs[[1]])
  for (i in seq_along(objs)[-1]) {
    b2 <- sf::st_bbox(objs[[i]])
    bb["xmin"] <- min(bb["xmin"], b2["xmin"])
    bb["ymin"] <- min(bb["ymin"], b2["ymin"])
    bb["xmax"] <- max(bb["xmax"], b2["xmax"])
    bb["ymax"] <- max(bb["ymax"], b2["ymax"])
  }
  bb
}

# --- Range map: Little (1977) pawpaw range (shapefile in zip)
read_shp_from_zip <- function(zip_path) {
  if (!file.exists(zip_path)) stop("Zip not found: ", zip_path)
  zl <- utils::unzip(zip_path, list = TRUE)
  shp <- zl$Name[grepl("\\.shp$", zl$Name, ignore.case = TRUE)][1]
  if (is.na(shp) || !nzchar(shp)) stop("No .shp file inside ", zip_path)
  zip_norm <- normalizePath(zip_path, winslash = "/", mustWork = TRUE)
  vsipath <- paste0("/vsizip/", zip_norm, "/", shp)
  sf::read_sf(vsipath, quiet = TRUE)
}

ensure_lonlat_crs <- function(x, label) {
  x <- sf::st_zm(x, drop = TRUE, what = "ZM")
  if (is.na(sf::st_crs(x))) {
    warning(
      "`", label, "` has no CRS; assuming WGS84 (EPSG:4326). ",
      "If the map is misaligned, set the correct CRS before reprojection."
    )
    sf::st_crs(x) <- 4326
  }
  x
}

pawpaw_zip <- file.path(data_dir, "InputData", "Pawpaw_range_Little_1977.zip")
if (!file.exists(pawpaw_zip)) {
  stop(
    "Pawpaw range zip not found.\nPlace `Pawpaw_range_Little_1977.zip` in:\n  ",
    dirname(pawpaw_zip)
  )
}
asimina_range <- ensure_lonlat_crs(read_shp_from_zip(pawpaw_zip), "Pawpaw range (Little 1977)")
asimina_range_proj <- st_transform(asimina_range, crs = proj_crs)

# Eastern Continental Divide (KMZ): try GDAL on archive + each KML/layer (doc.kml is not always the data)
read_kmz_as_sf <- function(kmz_path) {
  if (!file.exists(kmz_path)) stop("KMZ not found: ", kmz_path)
  nonempty <- function(x) {
    !is.null(x) && inherits(x, "sf") && nrow(x) > 0L && !all(sf::st_is_empty(x))
  }
  try_read <- function(dsn, layer = NULL) {
    if (is.null(layer)) {
      tryCatch(sf::read_sf(dsn, quiet = TRUE), error = function(e) NULL)
    } else {
      tryCatch(sf::read_sf(dsn, layer = layer, quiet = TRUE), error = function(e) NULL)
    }
  }
  try_layers <- function(dsn) {
    layers <- tryCatch(sf::st_layers(dsn), error = function(e) NULL)
    if (is.null(layers) || !nrow(layers)) return(NULL)
    for (i in seq_len(nrow(layers))) {
      r <- try_read(dsn, layers$name[i])
      if (nonempty(r)) return(r)
    }
    NULL
  }
  r <- try_read(kmz_path)
  if (nonempty(r)) return(r)
  r <- try_layers(kmz_path)
  if (nonempty(r)) return(r)
  zl <- utils::unzip(kmz_path, list = TRUE)
  kmls <- zl$Name[grepl("\\.kml$", zl$Name, ignore.case = TRUE)]
  if (!length(kmls)) stop("No .kml file inside ", kmz_path)
  kmz_norm <- normalizePath(kmz_path, winslash = "/", mustWork = TRUE)
  for (k in kmls) {
    vsipath <- paste0("/vsizip/", kmz_norm, "/", k)
    r <- try_read(vsipath)
    if (nonempty(r)) return(r)
    r <- try_layers(vsipath)
    if (nonempty(r)) return(r)
  }
  stop("Could not read non-empty geometries from KMZ: ", kmz_path)
}

prepare_ecd_geometry <- function(x) {
  x <- sf::st_make_valid(x)
  x <- sf::st_zm(x, drop = TRUE, what = "ZM")
  ok <- !sf::st_is_empty(sf::st_geometry(x))
  x <- x[ok, , drop = FALSE]
  if (nrow(x) == 0L) return(x)
  if (any(sf::st_is(x, "GEOMETRYCOLLECTION"))) {
    x <- tryCatch(
      sf::st_collection_extract(x, "LINESTRING"),
      error = function(e) x
    )
  }
  if (nrow(x) > 0L && all(sf::st_is(x, c("POLYGON", "MULTIPOLYGON")))) {
    x <- sf::st_sf(geometry = sf::st_boundary(sf::st_union(sf::st_geometry(x))))
  }
  if (nrow(x) > 0L && all(sf::st_is(x, c("LINESTRING", "MULTILINESTRING")))) {
    ug <- sf::st_union(sf::st_geometry(x))
    ug <- tryCatch(sf::st_line_merge(ug), error = function(e) ug)
    x <- sf::st_sf(geometry = ug)
  }
  x
}

ecd_kmz <- file.path(data_dir, "InputData", "Eastern_watershed.kml") ####Apalachain divide
if (!file.exists(ecd_kmz)) {
  stop(
    "Eastern Continental Divide KMZ not found.\nPlace `kml_13792.kmz` in:\n  ",
    dirname(ecd_kmz)
  )
}
appalach_divide_ll <- ensure_lonlat_crs(read_kmz_as_sf(ecd_kmz), "Eastern Continental Divide (KMZ)")
appalach_divide_ll <- prepare_ecd_geometry(appalach_divide_ll)
if (nrow(appalach_divide_ll) == 0L || all(sf::st_is_empty(sf::st_geometry(appalach_divide_ll)))) {
  stop("Eastern Continental Divide KMZ read as empty geometry after cleaning; check kml_13792.kmz.")
}
appalach_divide_proj <- sf::st_transform(appalach_divide_ll, crs = proj_crs)
if (all(sf::st_is_empty(sf::st_geometry(appalach_divide_proj)))) {
  stop("ECD is empty after projection to EPSG:", proj_epsg, "; check CRS in the KMZ.")
}

## ------------------------------------------------
## Shared map styling helpers
## ------------------------------------------------

map_theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(color = "grey35"),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(14, "pt"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

sample_point_layer <- function(data, size = 1.8) {
  geom_sf(
    data = data,
    color = "#333333",
    fill = "white",
    shape = 21,
    size = size,
    stroke = 0.35,
    inherit.aes = FALSE
  )
}

state_layer <- function(data) {
  geom_sf(
    data = data,
    fill = NA,
    color = "grey45",
    linewidth = 0.22,
    inherit.aes = FALSE
  )
}

# Little (1977) range polygon — use beneath wingen rasters so π / counts stay visible on top
asimina_range_underlay_layer <- function() {
  ggplot2::geom_sf(
    data  = asimina_range_proj,
    fill  = "grey90",
    alpha = 0.4,
    color = "grey60",
    linewidth = 0.4,
    inherit.aes = FALSE
  )
}

# ggplot_gd / ggplot_count append new layers after geom_tile, which paints the range on top.
# Prepend underlay layer(s) so they draw first (behind the diversity grid).
layer_beneath_wingen_raster <- function(p, ...) {
  bottom <- ggplot2::ggplot()
  for (ly in list(...)) {
    bottom <- bottom + ly
  }
  p$layers <- c(bottom$layers, p$layers)
  p
}

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
coords_proj <- st_transform(coord_pts_sf, crs = proj_crs)
states_proj <- st_transform(us_states_ec, crs = proj_crs)
# Plot limits including ECD (otherwise coord_sf can clip the divide off-map)
map_bb <- st_bbox_union(states_proj, appalach_divide_proj)

# Union outline segments that lie along natural-earth coastline (Atlantic/Gulf): draw lighter on top of state borders
us_coastal_border_sf <- tryCatch(
  {
    outer <- sf::st_boundary(sf::st_union(sf::st_make_valid(sf::st_geometry(states_proj))))
    outer_sf <- sf::st_sf(geometry = outer, crs = sf::st_crs(states_proj))
    coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
    coast <- sf::st_make_valid(coast)
    coast <- sf::st_transform(coast, proj_crs)
    buf <- sf::st_buffer(sf::st_union(coast), dist = 120000)
    hits <- sf::st_intersection(sf::st_make_valid(outer_sf), sf::st_make_valid(buf))
    if (inherits(sf::st_geometry(hits), "sfc_GEOMETRYCOLLECTION")) {
      hits <- tryCatch(
        sf::st_collection_extract(hits, "LINESTRING"),
        error = function(e) hits
      )
    }
    if (nrow(hits) > 0L && !all(sf::st_is_empty(sf::st_geometry(hits)))) hits else outer_sf[FALSE, ]
  },
  error = function(e) {
    warning("Could not build US coastal border overlay: ", conditionMessage(e))
    states_proj[FALSE, ]
  }
)
coastal_border_color <- "grey65"
coastal_border_linewidth <- 0.28
coastal_border_layer <- function() {
  if (nrow(us_coastal_border_sf) < 1L) return(NULL)
  geom_sf(
    data = us_coastal_border_sf,
    fill = NA,
    color = coastal_border_color,
    linewidth = coastal_border_linewidth,
    inherit.aes = FALSE
  )
}
png(filename = file.path(out_dir, paste0("projected_sample_locations_EPSG", proj_epsg, ".png")),
    width = 1600, height = 1200, res = 200)
plot(coords_proj["geometry"], main = paste0("Projected sample locations (EPSG:", proj_epsg, ")"))
dev.off()

# Pick a base environmental layer for plotting (BIO6) and aggregate for speed
# Convert to terra SpatRaster if needed
bio15 <- wclim[["bio15"]]

# Aggregate to coarser resolution (factor = 5)
envlayer <- terra::aggregate(bio15, fact = 5)


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

# π grid + window: `res` (meters in UTM) sets cell size; smaller cells → denser grid and a
# more continuous-looking map (less blocky). `wdim` = odd number of *raster cells* per
# side of the focal cell (e.g. 5 → 5×5). Physical window span ≈ wdim × res per axis.
# Defaults below: ~25 km cells × 5 cells ≈ 125 km window (similar to old 50 km × 3).
pi_res <- 20000 # 40000
pi_wdim <-  5

# Generate raster for sliding window (spacing & buffer in map units of coords_proj)
png(
  filename = file.path(out_dir, "coords_to_raster_window_grid.png"),
  width = 1600,
  height = 1200,
  res = 200
)
liz_lyr <- coords_to_raster(
  coords_proj,
  res    = pi_res,
  buffer = 5,
  plot   = TRUE
)
dev.off()

# Preview window and sample counts
sample_count_prev <- preview_gd(liz_lyr, coords_proj, wdim = pi_wdim, fact = 0)
sample_count_prev0 <- sample_count_prev
sample_count_prev0[sample_count_prev0 == 0] <- NA

p_prev <- ggplot_count(sample_count_prev0) +
  ggtitle("Preview sample count") +
  scale_fill_viridis_c(na.value = "grey") +
  state_layer(states_proj) +
  coastal_border_layer() +
  # ECD drawn last so it is not covered by state/raster layers
  geom_sf(
    data        = appalach_divide_proj,
    color       = "#2F1810",
    linewidth   = 0.8,
    linetype    = "dashed",
    inherit.aes = FALSE
  ) +
  ggplot2::coord_sf(
    crs    = proj_crs,
    xlim   = map_bb[c("xmin", "xmax")],
    ylim   = map_bb[c("ymin", "ymax")],
    expand = FALSE
  )
p_prev <- layer_beneath_wingen_raster(p_prev, asimina_range_underlay_layer())
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
  wdim   = pi_wdim,
  fact   = 0
)

## Manuscript stats: pi across windowed cells with >= 2 sampled populations
pi_min_populations <- 2L
pi_cell_df <- terra::as.data.frame(wgd, xy = TRUE, na.rm = TRUE)
if (!all(c("pi", "sample_count") %in% names(pi_cell_df))) {
  stop(
    "Expected `pi` and `sample_count` columns from window_gd(); got: ",
    paste(names(pi_cell_df), collapse = ", ")
  )
}
pi_cells_min2 <- pi_cell_df[
  !is.na(pi_cell_df$pi) &
    !is.na(pi_cell_df$sample_count) &
    pi_cell_df$sample_count >= pi_min_populations,
  ,
  drop = FALSE
]
if (nrow(pi_cells_min2) < 1L) {
  stop(
    "No windowed cells with sample_count >= ", pi_min_populations,
    "; cannot summarize nucleotide diversity."
  )
}
pi_summary <- data.frame(
  MIN_PI = min(pi_cells_min2$pi),
  MAX_PI = max(pi_cells_min2$pi),
  MEAN_PI = mean(pi_cells_min2$pi),
  MEDIAN_PI = stats::median(pi_cells_min2$pi),
  N_CELLS = nrow(pi_cells_min2),
  MIN_SAMPLE_COUNT = pi_min_populations,
  WINDOW_RES_M = pi_res,
  WINDOW_DIM_CELLS = pi_wdim,
  stringsAsFactors = FALSE
)
pi_summary$MANUSCRIPT_SENTENCE <- paste0(
  "Nucleotide diversity ranged from ",
  signif(pi_summary$MIN_PI, 4),
  " to ",
  signif(pi_summary$MAX_PI, 4),
  ", with a mean of ",
  signif(pi_summary$MEAN_PI, 4),
  " across all windowed cells containing at least two sampled populations."
)
write.table(
  pi_summary,
  file = file.path(out_dir, "nucleotide_diversity_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  pi_cells_min2,
  file = file.path(out_dir, "nucleotide_diversity_cells_min2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message(
  "[01-genetic-diversity] pi summary (n >= ", pi_min_populations, " populations, ",
  pi_summary$N_CELLS, " cells): ",
  pi_summary$MANUSCRIPT_SENTENCE
)

## ------------------------------------------------
## 4. Visualizing moving-window results (improved π map)
## ------------------------------------------------

# Projected sample points for overlay
pts_proj <- st_transform(coord_pts_sf, crs = proj_crs)

# Get π range and nice breaks
pi_rng <- terra::global(wgd[["pi"]], range, na.rm = TRUE)
pi_min <- as.numeric(pi_rng[1, 1])
pi_max <- as.numeric(pi_rng[1, 2])
pi_brks <- pretty(c(pi_min, pi_max), n = 5)

p_pi <- ggplot_gd(wgd) +
  ggtitle("Moving window π") +
  # State borders (match sample-count map: full extent + black outlines)
  geom_sf(
    data  = states_proj,
    fill  = NA,
    color = "black",
    linewidth = 0.25
  ) +
  coastal_border_layer() +
  geom_sf(
    data  = pts_proj,
    color = "black",
    fill  = "white",
    shape = 21,
    size  = 1.2,
    stroke = 0.25
  ) +
  geom_sf(
    data        = appalach_divide_proj,
    color       = "#2F1810",
    linewidth   = 0.8,
    linetype    = "dashed",
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_c(
    option    = "viridis",
    direction = 1,
    na.value  = NA,
    limits    = c(pi_min, pi_max),
    breaks    = pi_brks
  ) +
  ggplot2::coord_sf(
    crs  = proj_crs,
    xlim = map_bb[c("xmin", "xmax")],
    ylim = map_bb[c("ymin", "ymax")]
  ) +
  map_theme_pub()
p_pi <- layer_beneath_wingen_raster(p_pi, asimina_range_underlay_layer())

print(p_pi)

ggsave(
  filename = file.path(out_dir, "wingen_moving_window_pi.png"),
  plot = p_pi,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Sample count map: set unsampled (0) cells to NA so they can be colored grey
# 0-count -> NA so unsampled is grey
wgd_count0 <- wgd
wgd_count0[["sample_count"]][wgd_count0[["sample_count"]] == 0] <- NA

# projected sample points for overlay
pts_proj <- st_transform(coord_pts_sf, crs = proj_crs)

# compute good legend breaks
mx <- terra::global(wgd[["sample_count"]], "max", na.rm = TRUE)[1, 1]
mx <- if (is.na(mx)) 1 else as.numeric(mx)
brks <- pretty(c(1, mx), n = 5)

p_count <- ggplot_count(wgd_count0) +
  ggtitle("Sample count (moving window)") +
  state_layer(states_proj) +
  coastal_border_layer() +
  sample_point_layer(pts_proj, size = 1.7) +
  geom_sf(
    data        = appalach_divide_proj,
    color       = "#2F1810",
    linewidth   = 0.8,
    linetype    = "dashed",
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_c(
    option = "cividis",     # easier for colorblind + print
    na.value = "grey85",
    breaks = brks,
    limits = c(0, mx)
  ) +
  ggplot2::coord_sf(
    crs  = proj_crs,
    xlim = map_bb[c("xmin", "xmax")],
    ylim = map_bb[c("ymin", "ymax")]
  ) +
  map_theme_pub()
p_count <- layer_beneath_wingen_raster(p_count, asimina_range_underlay_layer())

ggsave(
  filename = file.path(out_dir, "wingen_moving_window_sample_count.png"),
  plot = p_count, width = 8, height = 6, units = "in", dpi = 300
)

## ------------------------------------------------
## 5. Kriging and masking
## ------------------------------------------------

# Kriging: prefer wkrig_gd() (gstat). Some R installs break gstat S4 methods
# (spacetime::STFDF vs raster::RasterBrick). We try unloading spacetime, then fall back to
# terra::interpIDW (no gstat) so the pipeline still runs.
try_unload_spacetime <- function() {
  if (!"spacetime" %in% loadedNamespaces()) return(invisible(NULL))
  tryCatch(
    {
      if ("package:spacetime" %in% search()) {
        detach("package:spacetime", character.only = TRUE, unload = TRUE, force = TRUE)
      }
      unloadNamespace("spacetime")
    },
    error = function(e) invisible(NULL)
  )
  invisible(NULL)
}

raster_idw_fallback <- function(r, grd) {
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
  pts <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  if (nrow(pts) < 3L) {
    stop("Interpolation needs at least 3 non-NA raster cells; got ", nrow(pts), ".")
  }
  names(pts)[3L] <- "value"
  v <- terra::vect(pts, geom = c("x", "y"), crs = terra::crs(r))
  ex <- as.vector(terra::ext(grd))
  rad <- sqrt((ex[2] - ex[1])^2 + (ex[4] - ex[3])^2) / 2 + max(terra::res(grd)) * 20
  out <- terra::interpIDW(grd, v, field = "value", radius = rad, power = 2, minPoints = 1L)
  names(out) <- names(r)
  out
}

wgd_krig_layer <- function(r, grd) {
  try_unload_spacetime()
  krig_try <- tryCatch(wkrig_gd(r, grd = grd), error = function(e) e)
  if (!inherits(krig_try, "error")) return(krig_try)
  msg <- conditionMessage(krig_try)
  if (!grepl(
    "coerce|duplicate class|package slot|pkgMethodLabel|STFDF|RasterBrick",
    msg,
    ignore.case = TRUE
  )) {
    stop(msg)
  }
  warning(
    "Variogram kriging (wkrig_gd / gstat) failed (S4 conflict in spatial packages).\n",
    "Using inverse-distance weights (terra::interpIDW) instead. ",
    "Fix long-term with: install.packages(c(\"sp\",\"spacetime\",\"raster\",\"gstat\"), dependencies=TRUE) in a fresh session.\n",
    "Original error: ",
    msg,
    call. = FALSE
  )
  raster_idw_fallback(r, grd)
}

krig_index <- 1:2
krig_disagg_grd <- 2L
grd_krig <- wgd[[krig_index[1]]]
if (krig_disagg_grd > 1L) {
  grd_krig <- terra::disagg(grd_krig, fact = krig_disagg_grd, method = "near")
}
krig_layers <- lapply(krig_index, function(i) wgd_krig_layer(wgd[[i]], grd_krig))
kgd <- terra::rast(krig_layers)
names(kgd) <- names(wgd)[krig_index]

# Kriged π: mask cells outside Asimina (Little 1977) distribution — π → NA outside species range
range_union_sfc <- sf::st_make_valid(sf::st_union(sf::st_geometry(asimina_range_proj)))
clip_range_vect <- if (length(range_union_sfc) > 0L && !all(sf::st_is_empty(range_union_sfc))) {
  terra::vect(range_union_sfc)
} else {
  NULL
}
if (is.null(clip_range_vect) || terra::nrow(clip_range_vect) == 0L) {
  warning("Asimina range union is empty; kriged π will not be masked to species distribution.")
}

mask_pi_to_asimina_range <- function(r_pi, clip_vect = clip_range_vect) {
  if (is.null(clip_vect) || terra::nrow(clip_vect) == 0L) return(r_pi)
  if (!terra::hasValues(r_pi)) return(r_pi)
  clip_v <- terra::project(clip_vect, terra::crs(r_pi))
  mask_r <- terra::rasterize(clip_v, r_pi, touches = TRUE)
  terra::mask(r_pi, mask_r)
}

summary(kgd)

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

# Kriged masked π maps: range + states + points like moving-window π.
# Species range must be drawn *under* π tiles (ggplot_gd puts geom_tile first; adding
# geom_sf after drew grey on top and hid colors). Build: range → tiles → borders → points.
for (spec in list(
  list(mgd = mgd_1, title = "Kriged & masked π", fname = "wingen_kriged_masked_pi.png"),
  list(mgd = mgd_2, title = "Kriged & masked π (n ≥ 2 populations)", fname = "wingen_kriged_masked_pi_min2samples.png")
)) {
  mgd_i <- spec$mgd
  if (!inherits(mgd_i, "SpatRaster")) {
    mgd_i <- terra::rast(mgd_i)
  }
  if (!"pi" %in% names(mgd_i)) {
    stop("Kriged stack has no layer named pi; found: ", paste(names(mgd_i), collapse = ", "))
  }
  # Single layer avoids rast()+subset edge cases that leave an uninitialized SpatRaster
  # (terra::mask then errors: [mask] SpatRaster has no values).
  r_pi <- mgd_i[["pi"]]
  if (!terra::hasValues(r_pi)) {
    stop("pi layer has no cell values (terra::hasValues is FALSE); check kriging / mask_gd output.")
  }
  # Grid-aligned mask: π NA outside Asimina distribution polygon
  r_pi <- mask_pi_to_asimina_range(r_pi)
  x_df <- terra::as.data.frame(r_pi, xy = TRUE)
  pi_col <- setdiff(names(x_df), c("x", "y"))[1]

  pi_rng_i <- terra::global(r_pi, range, na.rm = TRUE)
  pi_min_i <- as.numeric(pi_rng_i[1, 1])
  pi_max_i <- as.numeric(pi_rng_i[1, 2])
  pi_brks_i <- pretty(c(pi_min_i, pi_max_i), n = 5)

  p_masked <- ggplot() +
    geom_sf(
      data      = asimina_range_proj,
      fill      = ggplot2::alpha("grey90", 0.4),
      color     = "grey60",
      linewidth = 0.4
    ) +
    geom_tile(
      data = x_df,
      aes(
        x     = .data[["x"]],
        y     = .data[["y"]],
        fill  = .data[[pi_col]]
      )
    ) +
    geom_sf(
      data      = states_proj,
      fill      = NA,
      color     = "black",
      linewidth = 0.25,
      inherit.aes = FALSE
    ) +
    coastal_border_layer() +
    geom_sf(
      data        = pts_proj,
      color       = "black",
      fill        = "white",
      shape       = 21,
      size        = 1.2,
      stroke      = 0.25,
      inherit.aes = FALSE
    ) +
    geom_sf(
      data        = appalach_divide_proj,
      color       = "#2F1810",
      linewidth   = 0.8,
      linetype    = "dashed",
      inherit.aes = FALSE
    ) +
    ggtitle(spec$title) +
    labs(caption = "White circles show sampled populations; interpolated surface is clipped to the native range.") +
    labs(fill = pi_col) +
    scale_fill_viridis_c(
      option    = "viridis",
      direction = 1,
      na.value  = NA,
      limits    = c(pi_min_i, pi_max_i),
      breaks    = pi_brks_i
    ) +
    ggplot2::coord_sf(
      crs  = st_crs(states_proj),
      xlim = map_bb[c("xmin", "xmax")],
      ylim = map_bb[c("ymin", "ymax")]
    ) +
    map_theme_pub() +
    theme(plot.caption = element_text(color = "grey35", size = 9))

  print(p_masked)
  ggsave(
    filename = file.path(out_dir, spec$fname),
    plot = p_masked,
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
  )
}
