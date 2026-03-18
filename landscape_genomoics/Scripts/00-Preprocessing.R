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
  "ecodist",
  "RStoolbox",
  "here",
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
library(ecodist)
library(RStoolbox)
library(here)
library(rnaturalearth)
library(rnaturalearthdata)

# Convenience loader for algatr’s suggested packages
alazygatr_packages()

## ============================
## 1. Paths & input data
## ============================

## Fixed project path (no `setwd()`)
data_dir <- path.expand("~/Desktop/LandscapeGenomicsPipeline/")
out_dir <- file.path(data_dir, "Outputs", "00-preprocessing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read data
vcf <- read.vcfR("~/Desktop/LandscapeGenomicsPipeline/InputData/filtered_snps4landscape.vcf") # filtered upstream (see project commands)
coord <- read.table("~/Desktop/LandscapeGenomicsPipeline/InputData/sampling_coords.tsv", sep = "\t", header = TRUE) # must match order of vcf
stopifnot(all(c("x", "y") %in% names(coord)))

## ============================
## 2. Climate data
## ============================

# Pull WorldClim data at 0.5° resolution
wclim <- get_worldclim(coords = coord, res = 0.5)

## ============================
## 3. Genetic distance
## ============================

# 3.1 LD pruning (writes to ~/ )
vcf_ldpruned <- ld_prune(
  vcf = "~/Desktop/LandscapeGenomicsPipeline/InputData/filtered_snps4landscape.vcf",
  out_name = "SNPS_LDpruned_r0.9_n10_renamed",
  out_format = "vcf",
  ld.threshold = 0.9,
  slide.max.n = 10
)

vcf_pruned <- read.vcfR("~/Desktop/LandscapeGenomicsPipeline/InputData/SNPS_LDpruned_r0.9_n10_renamed.vcf")

# 3.2 Convert VCF to dosage
dosage <- vcf_to_dosage(vcf_pruned)

# Quick NA check
sum(is.na(dosage))  # should be > 0 before imputation

# 3.3 Impute missing values based on population structure
impute_K <- 1:6
str_dos <- str_impute(
  gen             = dosage,
  K               = impute_K,
  entropy         = TRUE,
  repetitions     = 5,
  quiet           = FALSE,
  save_output     = TRUE,
  output_filename = paste0("str_imputeK", min(impute_K), "-", max(impute_K), "r5")
)

# Confirm NAs removed
stopifnot(sum(is.na(str_dos)) == 0)

# 3.4 Genetic distance (Euclidean)
str_dist    <- as.matrix(ecodist::distance(str_dos, method = "euclidean"))
euc_gendist <- gen_dist(str_dos, dist_type = "euclidean")

# 3.5 Heatmap of Euclidean distance
jpeg(
  filename = file.path(out_dir, "euc_dist_heatmap.jpeg"),
  width    = 3000,
  height   = 3000,
  units    = "px",
  pointsize = 512,
  quality   = 100
)
gen_dist_hm(euc_gendist)
dev.off()

## ============================
## 4. Environmental QC & extraction
## ============================

# US states as sf polygons
us_states <- rnaturalearth::ne_states(
  country     = "United States of America",
  returnclass = "sf"
)

# Look at first BIO variable (e.g. BIO1)
wclim[[1]]

plot_bios <- function(bio_idx, filename) {
  png(filename = file.path(out_dir, filename), width = 1600, height = 1200)
  plot(wclim[[bio_idx]], col = viridis::turbo(100), axes = FALSE)
  plot(st_geometry(us_states), add = TRUE, border = "black", lwd = 0.5)
  points(coord$x, coord$y, pch = 19, col = "black")
  dev.off()
}

bio_to_plot <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
invisible(lapply(bio_to_plot, function(i) plot_bios(i, paste0("bio", i, "_plot.png"))))

# Convert coordinates to sf (WGS84)
coords_longlat <- st_as_sf(coord, coords = c("x", "y"), crs = 4326)

# Collinearity checks
cors_env     <- check_env(wclim) #The extracted values for 55 pairs of variables had correlation coefficients > 0.7. algatr recommends reducing collinearity by removing correlated variables or performing a PCA before proceeeding.
check_result <- check_vals(wclim, coords_longlat) #Warning: The extracted values for 41 pairs of variables had correlation coefficients > 0.7. algatr recommends reducing collinearity by removing correlated variables or performing a PCA before proceeeding.

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

# Select a subset of (less collinear) variables
keep_vars <- c("bio2",  "bio7", "bio8", "bio9", "bio13", "bio15", "bio18" )
selected_wclim <- wclim[[keep_vars]]

# Re‑check collinearity for subset
cors_env_sub  <- check_env(selected_wclim)
check_result2 <- check_vals(selected_wclim, coords_longlat)

# Create SpatVector and extract environmental data
sample_points <- vect(coord, geom = c("x", "y"), crs = "EPSG:4326")
env_dataID    <- terra::extract(selected_wclim, sample_points)
env_data      <- env_dataID[, -1, drop = FALSE]   # drop ID col safely

# Scale predictors for downstream models
env_scaled <- scale(env_data, center = TRUE, scale = TRUE)

## ============================
## 5. Save pipeline artifacts
## ============================
saveRDS(coord, file.path(out_dir, "coord.rds"))
saveRDS(vcf_pruned, file.path(out_dir, "vcf_pruned.rds"))
saveRDS(selected_wclim, file.path(out_dir, "selected_wclim.rds"))
saveRDS(env_data, file.path(out_dir, "env_data.rds"))
saveRDS(env_scaled, file.path(out_dir, "env_scaled.rds"))
saveRDS(euc_gendist, file.path(out_dir, "euc_gendist.rds"))
saveRDS(str_dos, file.path(out_dir, "str_dos.rds"))
