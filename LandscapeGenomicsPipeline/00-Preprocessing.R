## ============================
## 0. Reproducible setup
## ============================

## This script assumes packages are already installed.
## For long-term reproducibility, consider using `renv` in the project.

required_pkgs <- c(
  "algatr",
  "tidyverse",
  "vcfR",
  "hexbin",
  "sf",
  "terra",
  "viridis",
  "RStoolbox",
  "rnaturalearth",
  "rnaturalearthdata"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them once (e.g., install.packages(...)) and re-run."
  )
}

## Load libraries (explicit; avoids reliance on search path side-effects)
library(algatr)
library(tidyverse)
library(vcfR)
library(hexbin)
library(sf)
library(terra)
library(viridis)
library(RStoolbox)
library(rnaturalearth)
library(rnaturalearthdata)
library(here)

# Convenience loader for algatr’s suggested packages
alazygatr_packages()

## Pipeline cache helpers (RDS + raster skip). Options: Scripts/lgp_pipeline_cache.R
{
  ._lgpr <- Sys.getenv("LGP_PROJECT_ROOT", "~/Desktop/LandscapeGenomicsPipeline")
  options(lgp.project_root = sub("/+$", "", path.expand(getOption("lgp.project_root", ._lgpr))))
  suppressPackageStartupMessages(
    base::source(base::file.path(getOption("lgp.project_root"), "Scripts", "lgp_pipeline_cache.R"), encoding = "UTF-8")
  )
  rm(._lgpr)
}

## ============================
## 1. Paths & input data
## ============================

## Fixed project path (no `setwd()`): override with options(lgp.project_root = ...) or LGP_PROJECT_ROOT
data_dir <- lgp_project_root()
out_dir <- lgp_preprocess_dir(create_default = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(out_dir, "_step_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

coord_path <- file.path(data_dir, "InputData", "sampling_coords.tsv")

# Read data
coord <- read.table(coord_path, sep = "\t", header = TRUE) # must match order of vcf
stopifnot(all(c("x", "y") %in% names(coord)))

wc_path <- file.path(cache_dir, "worldclim_get_worldclim.rds")

## ============================
## 2. Climate data
## ============================

# Pull WorldClim data at 0.5° resolution (cached vs sampling_coords.tsv)
wclim <- lgp_rds_cached(
  path       = wc_path,
  dep_paths  = coord_path,
  meta       = list(res = 0.5),
  label      = "WorldClim (get_worldclim)",
  compute    = function() get_worldclim(coords = coord, res = 0.5)
)

## ============================
## 3. Genetic distance
## ============================

# 3.1 LD pruning (writes to ~/ per algatr; we skip when reheader output is fresher than source VCF)
raw_ld_vcf     <- file.path(data_dir, "InputData", "filtered_snps4landscape_reheader.vcf")
reheader_pruned <- file.path(data_dir, "InputData", "SNPS_LDpruned_r0.9_n10_reheader.vcf")

if (lgp_should_rerun_external(reheader_pruned, raw_ld_vcf)) {
  ld_prune(
    vcf         = raw_ld_vcf,
    out_name    = "SNPS_LDpruned_r0.9_n10_renamed",
    out_format  = "vcf",
    ld.threshold = 0.9,
    slide.max.n = 10
  )
} else {
  message("[lgp-cache] Skipping ld_prune (reuse): ", basename(reheader_pruned))
}

vcf_pruned <- read.vcfR(reheader_pruned)

# 3.2 Convert VCF to dosage (cached vs pruned VCF)
dos_cache <- file.path(cache_dir, "vcf_to_dosage_matrix.rds")
dosage <- lgp_rds_cached(
  dos_cache,
  dep_paths = reheader_pruned,
  label     = "vcf → dosage matrix",
  compute   = function() vcf_to_dosage(vcf_pruned)
)

# Quick NA check
sum(is.na(dosage))  # should be > 0 before imputation

# 3.3 Impute missing values based on population structure
impute_K <- 1:6
str_dos_cache <- file.path(cache_dir, "str_impute_str_dos.rds")
str_meta <- list(Kmin = min(impute_K), Kmax = max(impute_K), rep = 5L, entropy = TRUE)
str_dos <- lgp_rds_cached(
  str_dos_cache,
  dep_paths = c(coord_path, reheader_pruned),
  meta      = str_meta,
  label     = "str_impute (structure imputation)",
  compute   = function() {
    str_impute(
      gen             = dosage,
      K               = impute_K,
      entropy         = TRUE,
      repetitions     = 5,
      quiet           = FALSE,
      save_output     = TRUE,
      output_filename = paste0("str_imputeK", min(impute_K), "-", max(impute_K), "r5")
    )
  }
)

# Confirm NAs removed
stopifnot(sum(is.na(str_dos)) == 0)

# 3.4 Genetic distance (Euclidean) — reuse when str_impute cache unchanged
euc_cache <- file.path(cache_dir, "euc_gen_dist_euclidean.rds")
euc_gendist <- lgp_rds_cached(
  euc_cache,
  dep_paths = str_dos_cache,
  label     = "euclidean gen_dist(str_dos)",
  compute   = function() gen_dist(str_dos, dist_type = "euclidean")
)

# 3.5 Heatmap of Euclidean distance (skip when jpeg newer than genotype distance cache)
heatmap_path <- file.path(out_dir, "euc_dist_heatmap.jpeg")
lgp_begin_graphics_if_stale(
  heatmap_path,
  plot_expr = { gen_dist_hm(euc_gendist) },
  dep_paths = euc_cache,
  width = 3000,
  height = 3000,
  res = NA,
  units = "px",
  pointsize = 512,
  quality = 100
)

## ============================
## 4. Environmental QC & extraction
## ============================

# US states as sf polygons
us_states <- rnaturalearth::ne_states(
  country     = "United States of America",
  returnclass = "sf"
)

plot_bios_file <- function(bio_idx, filename) {
  outp <- file.path(out_dir, filename)
  deps <- c(coord_path, wc_path)
  lgp_png_if_stale(
    outp,
    plot_expr = {
      plot(wclim[[bio_idx]], col = viridis::turbo(100), axes = FALSE)
      plot(st_geometry(us_states), add = TRUE, border = "black", lwd = 0.5)
      points(coord$x, coord$y, pch = 19, col = "black")
      invisible(TRUE)
    },
    dep_paths = deps,
    width = 1600,
    height = 1200,
    res = NA
  )
}

bio_to_plot <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
invisible(lapply(bio_to_plot, function(i) plot_bios_file(i, paste0("bio", i, "_plot.png"))))

# Convert coordinates to sf (WGS84)
coords_longlat <- st_as_sf(coord, coords = c("x", "y"), crs = 4326)

# Collinearity checks
cors_env <- check_env(wclim) #The extracted values for 55 pairs of variables had correlation coefficients > 0.7. algatr recommends reducing collinearity by removing correlated variables or performing a PCA before proceeeding.
check_vals(wclim, coords_longlat) #Warning: The extracted values for 41 pairs of variables had correlation coefficients > 0.7. algatr recommends reducing collinearity by removing correlated variables or performing a PCA before proceeeding.

# Export the correlation matrix as TSV.
# Note: `cor_matrix` has row names; `col.names = NA` keeps headers aligned when row names are written.
write.table(
  cors_env$cor_matrix,
  file = file.path(out_dir, "bio_colinearity.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Mantel tests among environmental distances
png(filename = file.path(out_dir, "bio_colinearity.png"))
check_results <- check_dists(wclim, coord)
dev.off()

head(check_results$mantel_df)
write.table(check_results$mantel_df,
            file = file.path(out_dir, "mantel_test.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Select a subset of (less collinear) variables
keep_vars <- c("bio2",  "bio7", "bio8", "bio9", "bio13", "bio15", "bio18" )

env_extr_cache <- file.path(cache_dir, "env_extract_scaled.rds")
kv_key <- paste(keep_vars, collapse = ";")
env_ob <- lgp_list_cached(
  env_extr_cache,
  dep_paths = c(coord_path, wc_path),
  meta       = list(keep_vars = kv_key),
  label      = "env extract / scale bio predictors",
  compute_list = function() {
    sw <- wclim[[keep_vars]]
    pts <- vect(coord, geom = c("x", "y"), crs = "EPSG:4326")
    ed_id <- terra::extract(sw, pts)
    ed <- ed_id[, -1L, drop = FALSE]
    list(selected_wclim = sw, env_data = ed, env_scaled = scale(ed, center = TRUE, scale = TRUE))
  }
)
selected_wclim <- env_ob$selected_wclim
env_data <- env_ob$env_data
env_scaled <- env_ob$env_scaled

# Re‑check collinearity for subset
check_env(selected_wclim)
check_vals(selected_wclim, coords_longlat)

## ============================
## 5. Save pipeline artifacts
## ============================
saveRDS(coord, file.path(out_dir, "coord.rds"))
saveRDS(vcf_pruned, file.path(out_dir, "vcf_pruned.rds"))
saveRDS(lgp_pack_for_rds(selected_wclim), file.path(out_dir, "selected_wclim.rds"))
saveRDS(env_data, file.path(out_dir, "env_data.rds"))
saveRDS(env_scaled, file.path(out_dir, "env_scaled.rds"))
saveRDS(euc_gendist, file.path(out_dir, "euc_gendist.rds"))
saveRDS(str_dos, file.path(out_dir, "str_dos.rds"))

# SNP chromosome & position (for Manhattan plots in GEA / LFMM); matches colnames(str_dos)
snp_map <- data.frame(
  snp = colnames(str_dos),
  chr = as.character(vcfR::getCHROM(vcf_pruned)),
  pos = as.integer(vcfR::getPOS(vcf_pruned)),
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(snp_map) == ncol(str_dos),
  length(vcfR::getCHROM(vcf_pruned)) == ncol(str_dos)
)
saveRDS(snp_map, file.path(out_dir, "snp_map.rds"))

