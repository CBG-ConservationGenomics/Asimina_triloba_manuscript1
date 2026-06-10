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
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

set.seed(1234)

._lgpr <- Sys.getenv("LGP_PROJECT_ROOT", "~/Desktop/LandscapeGenomicsPipeline")
options(lgp.project_root = sub("/+$", "", path.expand(getOption("lgp.project_root", ._lgpr))))
suppressPackageStartupMessages(base::source(
  base::file.path(getOption("lgp.project_root"), "Scripts", "lgp_pipeline_cache.R"),
  encoding = "UTF-8"
))
rm(._lgpr)

## ------------------------------------------------
## 0. Paths and input loading
## ------------------------------------------------

data_dir <- lgp_project_root()
in_dir <- lgp_preprocess_dir()
out_dir <- lgp_outputs_step_dir("02-population-structure")

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
  selected_wclim <- lgp_read_rds(p)
}

# Build projected coords and kriging raster if not from 01-GeneticDiversity
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

str_rds_for_tess <- file.path(in_dir, "str_dos.rds")
coord_rds_path <- file.path(in_dir, "coord.rds")
tess_deps <- c(str_rds_for_tess, coord_rds_path)
tess_cache_path <- file.path(out_dir, "_step_cache", "tess_ktest_qmatrix.rds")
tess_meta <- list(Kmax = 10L)

tess_bundle <- lgp_list_cached(
  tess_cache_path,
  dep_paths = tess_deps,
  meta = tess_meta,
  label = "TESS k-test + qmatrix",
  compute_list = function() {
    tr <- tess_ktest(
      gen         = dosage,
      coords      = as.matrix(coord[, c("x", "y")]),
      Kvals       = seq_len(tess_meta$Kmax),
      ploidy      = 2,
      K_selection = "auto"
    )
    bk <- tr[["K"]]
    to <- tr$tess3_obj
    qx <- qmatrix(to, K = bk)
    list(tess3_result = tr, tess3_obj = to, bestK = bk, qmat = qx)
  }
)

tess3_result <- tess_bundle$tess3_result
tess3_obj <- tess_bundle$tess3_obj
bestK <- tess_bundle$bestK
qmat <- tess_bundle$qmat

# Must match ncol(qmat): selected K and qmatrix column count can disagree in edge cases
K_q <- ncol(as.matrix(qmat))
if (K_q != bestK) {
  warning(
    "ncol(qmat) is ", K_q, " but tess3 selected K = ", bestK,
    "; using K_q for ancestry colors and factor levels."
  )
}

# Same blue→orange ramp as PCA, dendrogram, and ggplot barplot
tess_ac_pal <- grDevices::colorRampPalette(c("#2E5AAC", "#E87B14"))(K_q)
qbar_k_names <- colnames(as.matrix(qmat))
if (is.null(qbar_k_names)) qbar_k_names <- paste0("V", seq_len(K_q))
tess_qbar_fill <- stats::setNames(tess_ac_pal, qbar_k_names)

# tess3r plot(..., method = "map.max"): col.palette is a list of length K; each [[k]] is a
# color ramp for image() (grey → component color), same layout as CreatePalette() but with
# peak colors matching tess_q_barplot_ggplot / tess_ac_pal (column order = qmat columns).
tess_map_col_palette <- lapply(seq_len(K_q), function(k) {
  grDevices::colorRampPalette(c("grey96", tess_ac_pal[k]))(9L)
})

qmat_mx <- as.matrix(qmat)
sample_ids <- if (!is.null(rownames(dosage))) {
  rownames(dosage)
} else if (!is.null(rownames(qmat_mx))) {
  rownames(qmat_mx)
} else {
  as.character(seq_len(nrow(qmat_mx)))
}
rownames(qmat_mx) <- sample_ids

saveRDS(
  list(
    qmat      = qmat_mx,
    bestK     = bestK,
    K_q       = K_q,
    sample_id = sample_ids
  ),
  file.path(out_dir, "tess_qmatrix.rds")
)

# Flat TSV of ancestry coefficients (one row per sample, one column per K component)
q_colnames <- colnames(qmat_mx)
if (is.null(q_colnames) || !length(q_colnames)) {
  q_colnames <- paste0("K", seq_len(K_q))
  colnames(qmat_mx) <- q_colnames
}
q_export <- data.frame(sample_id = sample_ids, qmat_mx, check.names = FALSE)
write.table(
  q_export,
  file = file.path(out_dir, "tess_qmatrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message(
  "[02-population-structure] Wrote TESS Q matrix (K = ", K_q, ", n = ", nrow(q_export),
  ") -> ", file.path(out_dir, "tess_qmatrix.tsv")
)

## ------------------------------------------------
## 2. PCA of genotype matrix (saved)
## ------------------------------------------------

stopifnot(nrow(str_dos) == nrow(qmat))
pca_cache <- file.path(out_dir, "_step_cache", "genotype_pc1_pc2_bundle.rds")
pca_bundle <- lgp_list_cached(
  pca_cache,
  dep_paths = c(str_rds_for_tess, tess_cache_path),
  meta      = tess_meta,
  label     = "PRCOMP(str_dos) + PC1–PC2 variance explained",
  compute_list = function() {
    gg <- stats::prcomp(str_dos, center = TRUE, scale = FALSE)
    list(
      pca_gen = gg,
      pve = as.numeric(100 * gg$sdev^2 / sum(gg$sdev^2))
    )
  }
)
pca_gen <- pca_bundle$pca_gen
pve <- pca_bundle$pve

pca_df <- data.frame(
  PC1     = pca_gen$x[, 1],
  PC2     = pca_gen$x[, 2],
  cluster = factor(max.col(qmat), levels = seq_len(K_q))
)

lab_pc1 <- sprintf("PC1 (%.1f%% variance)", pve[1])
lab_pc2 <- sprintf("PC2 (%.1f%% variance)", pve[2])

p_pca <- ggplot(pca_df, aes(x = .data$PC1, y = .data$PC2, color = .data$cluster)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = tess_ac_pal) +
  coord_equal() +
  labs(
    title    = "PCA of imputed dosage matrix",
    subtitle = paste0("Points colored by dominant TESS3 ancestry component (K = ", K_q, ")"),
    x        = lab_pc1,
    y        = lab_pc2,
    color    = "Dominant\ncomponent"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "genotype_pca_pc1_pc2.png"),
  plot     = p_pca,
  width    = 7,
  height   = 6,
  units    = "in",
  dpi      = 300
)

## ------------------------------------------------
## 2b. Sample map: dominant TESS ancestry (same basemap style as 01-genetic-diversity)
## ------------------------------------------------

# Pawpaw range: Little (1977) shapefile in zip under InputData (same as 01-GeneticDiversity_revised.R)
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

# Same as 01-GeneticDiversity.R: reads .kml or .kmz (nested KML) for Eastern / Appalachian divide.
read_kmz_as_sf <- function(kmz_path) {
  if (!file.exists(kmz_path)) stop("KML/KMZ not found: ", kmz_path)
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
  if (!grepl("\\.kmz$", kmz_path, ignore.case = TRUE)) {
    stop("Could not read non-empty geometries from KML: ", kmz_path)
  }
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
  stop("Could not read non-empty geometries from KML/KMZ: ", kmz_path)
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

us_states <- rnaturalearth::ne_states(
  country     = "United States of America",
  returnclass = "sf"
)
east_coast_abbrev <- c(
  "ME", "NH", "MA", "RI", "CT", "NY", "NJ",
  "DE", "MD", "VA", "NC", "SC", "GA", "FL",
  "AL", "LA", "TN", "MS", "AR", "KY",
  "OH", "MI", "IN", "IL", "MO", "IA",
  "PA", "TX", "OK", "KS"
)
us_states_ec <- dplyr::filter(us_states, .data$postal %in% east_coast_abbrev)
states_proj <- sf::st_transform(us_states_ec, crs = proj_crs)

pawpaw_zip <- file.path(data_dir, "InputData", "Pawpaw_range_Little_1977.zip")
rangemap_root <- file.path(lgp_outputs_base(data_dir), "01-rangemap")
asimina_range_path_gpkg <- file.path(rangemap_root, "Asimina_triloba_range_ecol.gpkg")
asimina_range_path_shp <- file.path(rangemap_root, "Asimina_triloba_range_ecol.shp")

if (file.exists(pawpaw_zip)) {
  asimina_range <- ensure_lonlat_crs(read_shp_from_zip(pawpaw_zip), "Pawpaw range (Little 1977)")
} else if (file.exists(asimina_range_path_gpkg)) {
  asimina_range <- sf::st_read(asimina_range_path_gpkg, layer = "range", quiet = TRUE)
} else if (file.exists(asimina_range_path_shp)) {
  asimina_range <- sf::st_read(asimina_range_path_shp, quiet = TRUE)
} else {
  stop(
    "Asimina / pawpaw range not found (needed for sample map).\nPlace `Pawpaw_range_Little_1977.zip` in:\n  ",
    dirname(pawpaw_zip),
    "\nOr run the rangemap script so one of these exists:\n  ",
    asimina_range_path_gpkg, "\n  ", asimina_range_path_shp
  )
}
asimina_range_proj <- sf::st_transform(asimina_range, crs = proj_crs)

# Eastern / Appalachian Continental Divide: prefer Eastern_watershed.kml (same as 01-GeneticDiversity.R),
# then optional GeoJSON, else a simplified WGS84 polyline for visualization only.
ecd_kml_path <- file.path(data_dir, "InputData", "Eastern_watershed.kml")
ecd_geojson_path <- file.path(data_dir, "InputData", "eastern_continental_divide.geojson")
ecd_source_note <- " (approximate trace)"
if (file.exists(ecd_kml_path)) {
  appalach_divide_ll <- ensure_lonlat_crs(
    read_kmz_as_sf(ecd_kml_path),
    "Eastern watershed (KML)"
  )
  appalach_divide_ll <- prepare_ecd_geometry(appalach_divide_ll)
  if (nrow(appalach_divide_ll) == 0L || all(sf::st_is_empty(sf::st_geometry(appalach_divide_ll)))) {
    stop(
      "Eastern_watershed.kml read as empty geometry after cleaning; check the file under:\n  ",
      dirname(ecd_kml_path)
    )
  }
  appalach_divide_proj <- sf::st_transform(appalach_divide_ll, crs = proj_crs)
  ecd_source_note <- " (from InputData Eastern_watershed.kml)"
} else if (file.exists(ecd_geojson_path)) {
  appalach_divide_proj <- sf::st_transform(sf::read_sf(ecd_geojson_path, quiet = TRUE), crs = proj_crs)
  ecd_source_note <- " (from InputData GeoJSON)"
} else {
  ecd_xy <- matrix(
    c(
      -78.85, 42.05,
      -79.35, 41.15,
      -79.90, 40.25,
      -80.45, 39.45,
      -80.95, 38.45,
      -81.30, 37.45,
      -81.55, 36.45,
      -81.90, 35.45,
      -82.45, 34.75,
      -83.15, 33.85,
      -83.80, 32.85,
      -84.45, 31.85,
      -85.05, 30.95
    ),
    ncol = 2L,
    byrow = TRUE
  )
  appalach_divide_proj <- sf::st_transform(
    sf::st_as_sf(sf::st_sfc(sf::st_linestring(ecd_xy), crs = 4326L)),
    crs = proj_crs
  )
}

stopifnot(nrow(coords_proj) == nrow(qmat))
pts_admix <- coords_proj
pts_admix$dominant_component <- factor(
  max.col(as.matrix(qmat)),
  levels = seq_len(K_q)
)

p_admix_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data  = asimina_range_proj,
    fill  = "grey90",
    alpha = 0.4,
    color = "grey60",
    linewidth = 0.4
  ) +
  ggplot2::geom_sf(data = states_proj, fill = NA, color = "black", linewidth = 0.25) +
  ggplot2::geom_sf(
    data = appalach_divide_proj,
    color = "#2F1810",
    linewidth = 0.55,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  ggplot2::geom_sf(
    data = pts_admix,
    ggplot2::aes(color = .data$dominant_component),
    size = 2.2,
    alpha = 0.9
  ) +
  ggplot2::scale_color_manual(
    name   = "Dominant\nTESS component",
    values = tess_ac_pal
  ) +
  ggplot2::coord_sf(
    xlim = sf::st_bbox(states_proj)[c("xmin", "xmax")],
    ylim = sf::st_bbox(states_proj)[c("ymin", "ymax")]
  ) +
  ggplot2::labs(
    title = "Samples by dominant TESS ancestry",
    subtitle = paste0(
      "K = ", K_q, " (medium blue → orange ramp); dashed line = Eastern Continental Divide",
      ecd_source_note
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "right"
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "tess_samples_admixture_map.png"),
  plot     = p_admix_map,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300
)

## ------------------------------------------------
## 3. Krige ancestry coefficients
## ------------------------------------------------

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

# Detach raster so terra::rast() / list coercions do not hit broken sp/raster S4 methods.
try_unload_raster <- function() {
  tryCatch(
    {
      if ("package:raster" %in% search()) {
        detach("package:raster", character.only = TRUE, unload = TRUE, force = TRUE)
      }
    },
    error = function(e) invisible(NULL)
  )
  if ("raster" %in% loadedNamespaces()) {
    tryCatch(unloadNamespace("raster"), error = function(e) invisible(NULL))
  }
  invisible(NULL)
}

spatial_s4_conflict_error <- function(e) {
  parts <- character(0L)
  p <- e
  for (i in 1:12L) {
    if (!inherits(p, "condition")) break
    parts <- c(parts, conditionMessage(p))
    np <- p[["parent"]]
    if (is.null(np) || identical(np, p)) break
    p <- np
  }
  txt <- paste(parts, collapse = "\n")
  grepl(
    paste0(
      "coerce|duplicate class|package slot|STFDF|RasterBrick|pkgMethodLabel|",
      "selecting a method for function|rast\\b|evaluating the argument .x."
    ),
    txt,
    ignore.case = TRUE
  )
}

# Avoid terra::rast(list_of_layers): on broken raster/sp installs that can re-trigger S4 coerce.
c_stack_spatraster <- function(lst) {
  n <- length(lst)
  if (n < 1L) stop("c_stack_spatraster: empty layer list")
  if (n == 1L) return(lst[[1L]])
  out <- lst[[1L]]
  for (i in 2:n) out <- c(out, lst[[i]])
  out
}

# Same outputs as algatr::tess_krig (multi-layer SpatRaster, K1…Kk) but uses
# terra::interpIDW — avoids automap/gstat S4 conflicts on some installs.
tess_krig_idw <- function(qmat, coords, grid, correct_kriged_Q = TRUE) {
  K <- ncol(qmat)
  if (!inherits(grid, "SpatRaster")) {
    grid <- terra::rast(grid)
  }
  base_grid <- grid
  krig_df <- if (inherits(coords, "sf")) {
    coords
  } else if (inherits(coords, "SpatVector")) {
    sf::st_as_sf(coords)
  } else {
    cm <- as.matrix(coords)
    sf::st_as_sf(
      data.frame(x = cm[, 1], y = cm[, 2]),
      coords = c("x", "y"),
      crs = sf::st_crs(base_grid)
    )
  }
  if (nrow(krig_df) != nrow(qmat)) {
    stop("tess_krig_idw: nrow(coords) must match nrow(qmat).")
  }
  ex <- as.vector(terra::ext(base_grid))
  rad <- sqrt((ex[2] - ex[1])^2 + (ex[4] - ex[3])^2) / 2 + max(terra::res(base_grid)) * 20

  layers <- vector("list", K)
  xy <- sf::st_coordinates(krig_df)
  for (k in seq_len(K)) {
    Qk <- qmat[, k, drop = TRUE]
    if (length(unique(Qk)) == 1L) {
      warning(
        "Only one unique Q value for K = ", k,
        " (same as tess_krig); layer will be constant."
      )
      layers[[k]] <- terra::init(base_grid, unique(Qk)[1])
      next
    }
    df <- data.frame(x = xy[, 1L], y = xy[, 2L], Q = Qk)
    v <- terra::vect(df, geom = c("x", "y"), crs = terra::crs(base_grid))
    layers[[k]] <- terra::interpIDW(
      base_grid, v, field = "Q", radius = rad, power = 2, minPoints = 1L
    )
  }
  krig_admix <- c_stack_spatraster(layers)
  grid_rs <- terra::resample(grid, krig_admix[[1]])
  krig_admix <- terra::mask(krig_admix, grid_rs)
  if (correct_kriged_Q) {
    krig_admix[krig_admix < 0] <- 0
    krig_admix[krig_admix > 1] <- 1
  }
  names(krig_admix) <- paste0("K", seq_len(K))
  krig_admix
}

tess_krig_safe <- function(qmat, coords, grid) {
  try_unload_spacetime()
  # Default FALSE: algatr::tess_krig → automap often breaks when sp/spacetime/raster S4 is inconsistent.
  # After a clean reinstall of those packages: options(lgp.try_tess_automap = TRUE)
  use_automap <- isTRUE(getOption("lgp.try_tess_automap", FALSE))
  run_idw <- function(reason, err) {
    try_unload_spacetime()
    try_unload_raster()
    warning(
      reason,
      "Using tess_krig_idw() (inverse distance) instead. ",
      "For variogram kriging, reinstall sp, spacetime, raster, gstat, automap from CRAN in a clean R session.\n",
      if (!missing(err)) conditionMessage(err) else "",
      call. = FALSE
    )
    tess_krig_idw(qmat, coords, grid)
  }
  if (!use_automap) {
    message(
      "tess_krig_safe: using IDW admixture surfaces (default). ",
      "Set options(lgp.try_tess_automap = TRUE) to try algatr::tess_krig / automap."
    )
    try_unload_raster()
    return(tess_krig_idw(qmat, coords, grid))
  }
  tryCatch(
    tess_krig(qmat = qmat, coords = coords, grid = grid),
    error = function(e) {
      if (!spatial_s4_conflict_error(e)) stop(e)
      run_idw("tess_krig() / automap failed (S4 conflict in spatial packages). ", e)
    }
  )
}

# Reproject kriging raster to match projected coordinates
krig_raster_proj <- terra::project(krig_raster, proj_crs$wkt)

krig_admix <- tess_krig_safe(
  qmat = qmat,
  coords = coords_proj,
  grid   = krig_raster_proj
)

## ------------------------------------------------
## 4. Barplots of ancestry (saved)
## ------------------------------------------------

# tess3r::barplot.tess3Q() calls graphics::barplot(..., col = colpal, ...). Passing col= here
# duplicates "col" and errors, often leaving an empty PNG if the device is closed anyway.
png(file.path(out_dir, "tess_q_barplot_base.png"), width = 1600, height = 800, res = 200)
barplot(
  qmat,
  sort.by.Q     = TRUE,
  col.palette   = tess_map_col_palette,
  palette.length = 9L,
  border        = NA,
  space         = 0,
  xlab          = "Individuals",
  ylab          = "Ancestry coefficients"
)
dev.off()

p_qbar <- tess_ggbarplot(
  qmat,
  legend      = TRUE,
  ggplot_fill = ggplot2::scale_fill_manual(values = tess_qbar_fill, name = "Ancestry")
) +
  ggtitle(paste0("TESS3 ancestry coefficients (K = ", K_q, ")"))

ggsave(
  filename = file.path(out_dir, "tess_q_barplot_ggplot.png"),
  plot     = p_qbar,
  width    = 8,
  height   = 4,
  units    = "in",
  dpi      = 300
)

## ------------------------------------------------
## 5. Admixture map (single, publication-style plot)
## ------------------------------------------------

p_tess <- tess_ggplot(
  krig_admix,
  plot_method = "maxQ",
  ggplot_fill = scale_fill_viridis_d(option = "magma"),  # <‑‑ use _d
  plot_axes   = TRUE,
  coords      = coords_proj
) +
  ggtitle(paste0("TESS3 admixture map (K = ", K_q, ")"))

ggsave(
  filename = file.path(out_dir, "tess_admixture_maxQ.png"),
  plot     = p_tess,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300
)

## ------------------------------------------------
## 6. Optional: base TESS map with custom colormap
## ------------------------------------------------

coords_proj_mat <- st_coordinates(coords_proj)

png(file.path(out_dir, "tess_admixture_map_base.png"),
    width = 1600, height = 1200, res = 200)
plot(
  qmat,
  coords_proj_mat,
  method    = "map.max",
  interpol  = FieldsKrigModel(10),
  main      = paste0("Ancestry coefficients (K = ", K_q, ")"),
  xlab      = "x", ylab = "y",
  col.palette = tess_map_col_palette,
  resolution  = c(300, 300),
  cex         = 0.4
)
dev.off()

#############################
# Hclust
#############################
# Hierarchical clustering from genetic distance (pipeline artifact)
library(ggtree)
library(ape)

## `in_dir` already set via lgp_preprocess_dir() at top of script
hclust_out_dir <- out_dir
#dir.create(hclust_out_dir, recursive = TRUE, showWarnings = FALSE)

save_dpi <- 600
tess_png_w <- 3600
tess_png_h <- 2400
tess_png_res <- 300

required_pkgs <- c("ggplot2", "ggtree", "ape")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install once and re-run.\n",
    "For ggtree (Bioconductor):\n",
    "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
    "  BiocManager::install('ggtree')"
  )
}

as_dist <- function(x) {
  if (inherits(x, "dist")) return(x)
  if (is.matrix(x)) return(stats::as.dist(x))
  if (is.data.frame(x)) return(stats::as.dist(as.matrix(x)))
  if (is.list(x)) {
    for (nm in c("dist", "dist_mat", "distance", "D")) {
      if (!is.null(x[[nm]])) return(as_dist(x[[nm]]))
    }
  }
  stop("Don't know how to convert `euc_gendist` to a `dist` object.")
}

if (!exists("euc_gendist", inherits = TRUE)) {
  p_g <- file.path(in_dir, "euc_gendist.rds")
  if (!file.exists(p_g)) stop("Missing euc_gendist (expected ", p_g, ")")
  euc_gendist <- readRDS(p_g)
}

dist_matrix <- as_dist(euc_gendist)
hc <- hclust(dist_matrix, method = "average")

newick_path <- file.path(hclust_out_dir, "hclust_average.newick")
ape::write.tree(ape::as.phylo(hc), file = newick_path)

p_circ <- ggtree(hc, layout = "circular") +
  geom_tiplab(size = 0.65) +
  ggtitle("Hierarchical clustering (Euclidean genetic distance)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(hclust_out_dir, "hclust_dendrogram.png"),
  plot     = p_circ,
  width    = 10,
  height   = 10,
  units    = "in",
  dpi      = save_dpi
)

## ------------------------------------------------
## TESS ancestry–colored branches 
## ------------------------------------------------

tess_rds <- file.path(out_dir, "tess_qmatrix.rds")
if (file.exists(tess_rds)) {
  tess <- readRDS(tess_rds)
  q_tess <- as.matrix(tess$qmat)
  K_tess <- tess$bestK
  sid <- tess$sample_id
  stopifnot(length(sid) == nrow(q_tess))
  
  phy <- ape::as.phylo(hc)
  tip_match <- match(phy$tip.label, sid)
  if (anyNA(tip_match)) {
    warning(
      "Some dendrogram tips are missing from tess_qmatrix.rds$sample_id; ",
      "skipping TESS-colored dendrogram plots."
    )
  } else {
    ntip <- ape::Ntip(phy)
    tip_comp <- max.col(q_tess)[tip_match]
    pal <- grDevices::colorRampPalette(c("#2E5AAC", "#E87B14"))(K_tess)
    names(pal) <- as.character(seq_len(K_tess))
    
    edge <- phy$edge
    child_lists <- split(edge[, 2], edge[, 1])
    
    tips_under <- function(node) {
      if (node <= ntip) return(node)
      ch <- child_lists[[as.character(node)]]
      unlist(lapply(ch, tips_under), use.names = FALSE)
    }
    
    edge_majority_comp <- apply(edge, 1L, function(e) {
      ch <- e[2L]
      tips <- if (ch <= ntip) ch else tips_under(ch)
      tab <- table(tip_comp[as.integer(tips)])
      as.integer(names(tab)[which.max(tab)])
    })
    
    edge_col <- pal[as.character(edge_majority_comp)]
    tip_col <- pal[as.character(tip_comp)]
    
    png(
      file.path(hclust_out_dir, "hclust_dendrogram_tess_branches.png"),
      width  = tess_png_w,
      height = tess_png_h,
      res    = tess_png_res
    )
    plot(
      phy,
      type       = "fan",
      edge.color = edge_col,
      tip.color  = tip_col,
      cex        = 0.28,
      main       = "Hierarchical clustering (branches by majority TESS ancestry)",
      sub        = ""
    )
    dev.off()
      }
} else {
  message(
    "No tess_qmatrix.rds found (expected after running the TESS section of this script). ",
    "Skipping TESS-colored dendrogram plots."
  )
}
