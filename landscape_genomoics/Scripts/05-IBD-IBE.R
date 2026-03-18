#######################################################################
# Isolation by distance (IBD) and isolation by environment (IBE)
# with MMRR and GDM
# References:
# - GDM:  https://thewanglab.github.io/algatr/articles/GDM_vignette.html
# - MMRR: https://thewanglab.github.io/algatr/articles/MMRR_vignette.html
#######################################################################

library(algatr)
gdm_packages()  # loads gdm, ecodist, etc.
library(terra)
library(viridis)
library(dplyr)

set.seed(1234)

## ------------------------------------------------
## 0. Checks and shared objects
## ------------------------------------------------

data_dir <- path.expand("~/Desktop/LandscapeGenomicsPipeline/")
in_dir <- file.path(data_dir, "outputs", "00-preprocessing")
out_dir <- file.path(data_dir, "outputs", "05-ibd-ibe")
plot_dir <- file.path(out_dir, "plots")
results_dir <- file.path(out_dir, "results")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

load_if_missing <- function(obj, rds_name) {
  if (!exists(obj, inherits = TRUE)) {
    p <- file.path(in_dir, rds_name)
    if (!file.exists(p)) stop("Missing required input: ", obj, " (expected ", p, ")")
    assign(obj, readRDS(p), envir = .GlobalEnv)
  }
}

load_if_missing("coord", "coord.rds")
load_if_missing("env_scaled", "env_scaled.rds")
load_if_missing("euc_gendist", "euc_gendist.rds")
load_if_missing("selected_wclim", "selected_wclim.rds")
load_if_missing("str_dos", "str_dos.rds")

# ensure numeric coord matrix for distance functions
coord_xy <- as.matrix(coord[, c("x", "y")])

###################################################################
# 01. Multiple matrix regression with randomization (MMRR)
#     (individual-based)
###################################################################

# Genetic distance matrix (already euclidean from earlier)
Y <- as.matrix(euc_gendist)

# Environmental distances (on scaled predictors)
X <- env_dist(env_scaled)

# Add geographic distance
X[["geodist"]] <- geo_dist(coord_xy)

# Full MMRR
set.seed(10)
results_full <- mmrr_run(Y, X, nperm = 999, stdz = TRUE, model = "full")

# Best-model MMRR (stepwise selection)
set.seed(11)
results_best <- mmrr_run(Y, X, nperm = 999, stdz = TRUE, model = "best")

# Plots – save to disk
png(file.path(plot_dir, "mmrr_vars_full.png"), width = 1600, height = 1200, res = 200)
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "vars", stdz = TRUE)
dev.off()

png(file.path(plot_dir, "mmrr_fitted_full.png"), width = 1600, height = 1200, res = 200)
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "fitted", stdz = TRUE)
dev.off()

png(file.path(plot_dir, "mmrr_cov_full.png"), width = 1600, height = 1200, res = 200)
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "cov", stdz = TRUE)
dev.off()

png(file.path(plot_dir, "mmrr_all_best.png"), width = 1600, height = 1200, res = 200)
mmrr_plot(Y, X, mod = results_best$mod, plot_type = "all", stdz = TRUE)
dev.off()

# Write table out
mmrr_out <- mmrr_table(results_full, digits = 2, summary_stats = TRUE)
write.table(mmrr_out, file.path(results_dir, "mmrr_results_full.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


################################################################
# 02. Generalized Dissimilarity Modeling (GDM)
#     (population-based)
################################################################

## 02a. Allele frequencies per population

sample_names <- rownames(str_dos)
pops         <- sub("-.*", "", sample_names)
pops         <- sub("WC_PA_.*", "WC_PA", pops)
popmap       <- data.frame(sample = sample_names, population = pops)

pop_assignments <- popmap$population[match(rownames(str_dos), popmap$sample)]

calc_pop_allele_freq <- function(dosage_matrix, pop_assignments) {
  stopifnot(nrow(dosage_matrix) == length(pop_assignments))
  pops <- unique(pop_assignments)
  
  pop_mm <- model.matrix(~ 0 + factor(pop_assignments, levels = pops))
  obs    <- !is.na(dosage_matrix)
  dos0   <- replace(dosage_matrix, is.na(dosage_matrix), 0)
  
  allele_sums <- t(pop_mm) %*% dos0
  n_obs       <- t(pop_mm) %*% obs
  
  allele_freqs <- allele_sums / (2 * n_obs)
  allele_freqs[is.na(allele_freqs)] <- NA
  rownames(allele_freqs) <- pops
  allele_freqs
}

allele_freqs <- calc_pop_allele_freq(str_dos, pop_assignments)

## 02b. Genetic distance between populations (Euclidean on freqs)

calc_euclid_genetic_matrix_complete <- function(allele_freqs,
                                                sqrt = TRUE,
                                                fill = c("mean","max","min","zero"),
                                                as_dist = FALSE) {
  fill <- match.arg(fill)
  P <- as.matrix(allele_freqs)
  n <- nrow(P)
  mat <- matrix(NA_real_, n, n, dimnames = list(rownames(P), rownames(P)))
  
  for (i in seq_len(n)) {
    for (j in i:n) {
      valid <- !is.na(P[i, ]) & !is.na(P[j, ])
      d <- if (any(valid)) {
        d2 <- sum((P[i, valid] - P[j, valid])^2)
        if (sqrt) sqrt(d2) else d2
      } else {
        NA_real_
      }
      mat[i, j] <- d
      mat[j, i] <- d
    }
  }
  
  if (any(is.na(mat))) {
    vals <- na.omit(mat[upper.tri(mat)])
    if (length(vals) == 0) stop("No pairwise distances could be computed (no overlapping loci).")
    repl <- switch(fill,
                   mean = mean(vals),
                   max  = max(vals),
                   min  = min(vals),
                   zero = 0
    )
    mat[is.na(mat)] <- repl
  }
  
  if (as_dist) return(as.dist(mat))
  mat
}

gendist_euc_mat <- calc_euclid_genetic_matrix_complete(allele_freqs)

## 02c. Population coordinates and environment

coord_with_samples <- coord %>%
  mutate(sample = sample_names, population = pop_assignments)

centroids <- coord_with_samples %>%
  group_by(population) %>%
  summarize(
    centroid_x = mean(x, na.rm = TRUE),
    centroid_y = mean(y, na.rm = TRUE),
    .groups    = "drop"
  )

coord_pop <- centroids %>%
  left_join(coord_with_samples, by = "population") %>%
  group_by(population) %>%
  slice_min(order_by = (x - centroid_x)^2 + (y - centroid_y)^2,
            n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(population, sample, x, y)

centroid_coords <- as.data.frame(coord_pop[, c("x", "y")])

env_selected <- env_scaled[match(coord_pop$sample, sample_names), , drop = FALSE]

## 02d. GDM with all variables

gdm_full <- gdm_run(
  gendist       = gendist_euc_mat,
  coords        = centroid_coords,
  env           = env_selected,
  model         = "full",
  scale_gendist = TRUE
)

summary(gdm_full$model)

# Dissimilarity vs distance plots
png(file.path(plot_dir, "gdm_dissimilarity.png"), width = 1600, height = 1200, res = 200)
gdm_plot_diss(gdm_full$model)
dev.off()

# I-spline plots
png(file.path(plot_dir, "gdm_isplines_free.png"), width = 1600, height = 1200, res = 200)
gdm_plot_isplines(gdm_full$model, scales = "free")
dev.off()

png(file.path(plot_dir, "gdm_isplines_free_x.png"), width = 1600, height = 1200, res = 200)
gdm_plot_isplines(gdm_full$model, scales = "free_x")
dev.off()

# Write table out
tab <- gdm_table(gdm_full)
write.table(tab,
            file = file.path(results_dir, "gdm_table_full.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## 02e. Variable importance

gdmData <- gdm_format(
  gendist       = gendist_euc_mat,
  coords        = centroid_coords,
  env           = env_selected,
  scale_gendist = TRUE
)

varimp <- gdm.varImp(gdmData, geo = TRUE, nPerm = 50)

# Write table out
varimp_out <- gdm_varimp_table(varimp)
write.table(varimp_out, file.path(results_dir, "gdm_varimp_table.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)



## 02f. Spatial GDM map

# Use selected_wclim as the environmental raster stack for mapping
# Visualize GDM results
gdm_map_res <- gdm_map(gdm_full$model, selected_wclim, centroid_coords)

# Transformed rasters (per variable) and PCA RGB
png(file.path(plot_dir, "gdm_rastTrans.png"), width = 1600, height = 1200, res = 200)
plot(gdm_map_res$rastTrans, col = viridis(100))
dev.off()

png(file.path(plot_dir, "gdm_pcaRastRGB.png"), width = 1600, height = 1200, res = 200)
plot(gdm_map_res$pcaRastRGB)
dev.off()

## 02g. Masking out irrelevant areas

maprgb <- gdm_map_res$pcaRastRGB

map_mask <- extrap_mask(
  coords        = centroid_coords,
  envlayers = selected_wclim,
  method        = "buffer",
  buffer_width  = 1.25
)

png(file.path(plot_dir, "gdm_masked_map.png"), width = 1600, height = 1200, res = 200)
terra::plotRGB(maprgb)
terra::plot(map_mask, col = "white", add = TRUE, legend = FALSE, alpha = 0.6)
dev.off()
