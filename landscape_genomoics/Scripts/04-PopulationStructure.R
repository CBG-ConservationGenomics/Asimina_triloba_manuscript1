#######################################################################
# Population Structure with TESS3 (via algatr)
# https://thewanglab.github.io/algatr/articles/TESS_vignette.html
#######################################################################

library(algatr)
library(tess3r)      # if not already loaded by algatr
library(sf)
library(terra)
library(viridis)
library(ggplot2)

set.seed(1234)

## ------------------------------------------------
## 0. Paths and input loading
## ------------------------------------------------

data_dir <- path.expand("~/Desktop/LandscapeGenomicsPipeline/")
in_dir   <- file.path(data_dir, "outputs", "00-preprocessing")
out_dir  <- file.path(data_dir, "outputs", "04-population-structure")
#plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load from 00-preprocessing if not already in session
if (!exists("coord", inherits = TRUE)) {
  p <- file.path(in_dir, "coord.rds")
  if (!file.exists(p)) stop("Run 00-Preprocessing.R first (expected ", p, ")")
  coord <- readRDS(p)
}
if (!exists("dosage", inherits = TRUE)) {
  p <- file.path(in_dir, "str_dos.rds")
  if (!file.exists(p)) stop("Run 00-Preprocessing.R first (expected ", p, ")")
  dosage <- readRDS(p)
}
if (!exists("selected_wclim", inherits = TRUE)) {
  p <- file.path(in_dir, "selected_wclim.rds")
  if (!file.exists(p)) stop("Run 00-Preprocessing.R first (expected ", p, ")")
  selected_wclim <- readRDS(p)
}

# Build projected coords and kriging raster if not from 02-GeneticDiversity
if (!exists("coords_proj", inherits = TRUE) || !exists("proj_crs", inherits = TRUE)) {
  coord_sf   <- st_as_sf(coord, coords = c("x", "y"), crs = 4326, remove = FALSE)
  pick_utm   <- function(lon, lat) { zone <- floor((lon + 180) / 6) + 1; if (lat >= 0) 32600 + zone else 32700 + zone }
  proj_epsg  <- pick_utm(mean(coord$x, na.rm = TRUE), mean(coord$y, na.rm = TRUE))
  proj_crs   <- st_crs(proj_epsg)
  coords_proj <- st_transform(coord_sf, crs = proj_crs)
}
if (!exists("envlayer", inherits = TRUE) || !inherits(envlayer, "SpatRaster")) {
  envlayer <- terra::aggregate(selected_wclim[[1]], fact = 5)
  envlayer <- terra::project(envlayer, proj_crs$wkt)
}

# Create coarser raster for kriging (1st environmental layer)
krig_raster <- terra::aggregate(envlayer[[1]], fact = 10)

## ------------------------------------------------
## 1. Automatic K selection and Q-matrix
## ------------------------------------------------

tess3_result <- tess_ktest(
  gen        = dosage,
  coords     = as.matrix(coord[, c("x", "y")]),
  Kvals      = 1:10,
  ploidy     = 2,
  K_selection = "auto"
)

tess3_obj <- tess3_result$tess3_obj
bestK     <- tess3_result[["K"]]

qmat <- qmatrix(tess3_obj, K = bestK)

## ------------------------------------------------
## 2. Krige ancestry coefficients
## ------------------------------------------------

# Reproject kriging raster to match projected coordinates
krig_raster_proj <- terra::project(krig_raster, proj_crs$wkt)

krig_admix <- tess_krig(
  qmat = qmat,
  coords = coords_proj,
  grid   = krig_raster_proj
)

## ------------------------------------------------
## 3. Barplots of ancestry (saved)
## ------------------------------------------------

png(file.path(out_dir, "tess_q_barplot_base.png"), width = 1600, height = 800, res = 200)
barplot(qmat, sort.by.Q = TRUE, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry coefficients")
dev.off()

p_qbar <- tess_ggbarplot(qmat, legend = TRUE) +
  ggtitle(paste0("TESS3 ancestry coefficients (K = ", bestK, ")"))

ggsave(
  filename = file.path(out_dir, "tess_q_barplot_ggplot.png"),
  plot     = p_qbar,
  width    = 8,
  height   = 4,
  units    = "in",
  dpi      = 300
)

## ------------------------------------------------
## 4. Admixture map (single, publication-style plot)
## ------------------------------------------------

p_tess <- tess_ggplot(
  krig_admix,
  plot_method = "maxQ",
  ggplot_fill = scale_fill_viridis_d(option = "magma"),  # <‑‑ use _d
  plot_axes   = TRUE,
  coords      = coords_proj
) +
  ggtitle(paste0("TESS3 admixture map (K = ", bestK, ")"))

ggsave(
  filename = file.path(out_dir, "tess_admixture_maxQ.png"),
  plot     = p_tess,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300
)

## ------------------------------------------------
## 5. Optional: base TESS map with custom colormap
## ------------------------------------------------

coords_proj_mat <- st_coordinates(coords_proj)

png(file.path(out_dir, "tess_admixture_map_base.png"),
    width = 1600, height = 1200, res = 200)
plot(
  qmat,
  coords_proj_mat,
  method    = "map.max",
  interpol  = FieldsKrigModel(10),
  main      = paste0("Ancestry coefficients (K = ", bestK, ")"),
  xlab      = "x", ylab = "y",
  col.palette = CreatePalette(),
  resolution  = c(300, 300),
  cex         = 0.4
)
dev.off()
