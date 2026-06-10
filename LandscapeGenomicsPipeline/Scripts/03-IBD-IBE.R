#######################################################################
# Isolation by distance (IBD) and isolation by environment (IBE)
# with MMRR and GDM
# References:
# - GDM:  https://thewanglab.github.io/algatr/articles/GDM_vignette.html
# - MMRR: https://thewanglab.github.io/algatr/articles/MMRR_vignette.html
#######################################################################

library(algatr)
gdm_packages()  # loads gdm, ecodist, etc.

# algatr::gdm_map -> gdm_plot_vars: background points filter PC vs RGB rasters
# separately, so nrow(rastvals) != length(rastpcacols); ggplot2 3.5+ errors on
# colour length mismatch. Align rows before building colours.
gdm_plot_vars_fixed <- function(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2",
                                scl = 1, display_axes = FALSE) {
  if (terra::nlyr(pcaRastRGB) > 3) {
    stop("Only three PC layers (RGB) can be used for creating the variable plot (too many provided)")
  }
  if (terra::nlyr(pcaRastRGB) < 3) {
    stop("Need exactly three PC layers (RGB) for creating the variable plot (too few provided)")
  }
  xpc <- data.frame(pcaSamp$x[, 1:3])
  varpc <- data.frame(varnames = rownames(pcaSamp$rotation), pcaSamp$rotation)
  pcavals <- data.frame(terra::extract(pcaRast, coords, ID = FALSE))
  colnames(pcavals) <- colnames(xpc)
  scldat <- min(
    (max(pcavals[, y], na.rm = TRUE) -
      min(pcavals[, y], na.rm = TRUE) / (max(varpc[, y], na.rm = TRUE) - min(varpc[, y], na.rm = TRUE))),
    (max(pcavals[, x], na.rm = TRUE) -
      min(pcavals[, x], na.rm = TRUE) / (max(varpc[, x], na.rm = TRUE) - min(varpc[, x], na.rm = TRUE)))
  )
  varpc <- data.frame(varpc, v1 = scl * scldat * varpc[, x], v2 = scl * scldat * varpc[, y])
  pcavalsRGB <- data.frame(terra::extract(pcaRastRGB, coords, ID = FALSE))
  colnames(pcavalsRGB) <- colnames(xpc)
  ok_pts <- stats::complete.cases(pcavals) & stats::complete.cases(pcavalsRGB)
  pcavals <- pcavals[ok_pts, , drop = FALSE]
  pcavalsRGB <- pcavalsRGB[ok_pts, , drop = FALSE]
  pcacols <- apply(pcavalsRGB, 1, algatr:::create_rgb_vec)
  s <- sample(1:terra::ncell(pcaRast), 10000)
  rastvals <- data.frame(terra::values(pcaRast))[s, ]
  colnames(rastvals) <- colnames(xpc)
  rastvalsRGB <- data.frame(terra::values(pcaRastRGB))[s, ]
  colnames(rastvalsRGB) <- colnames(rastvals)
  ok_bg <- stats::complete.cases(rastvals) & stats::complete.cases(rastvalsRGB)
  rastvals <- rastvals[ok_bg, , drop = FALSE]
  rastvalsRGB <- rastvalsRGB[ok_bg, , drop = FALSE]
  rastpcacols <- apply(rastvalsRGB, 1, algatr:::create_rgb_vec)
  plot <- ggplot2::ggplot() + {
    if (display_axes) ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, col = "gray")
  } + {
    if (display_axes) ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, col = "gray")
  } +
    ggplot2::geom_point(
      data = rastvals,
      ggplot2::aes_string(x = x, y = y),
      colour = rastpcacols,
      size = 4,
      alpha = 0.02
    ) +
    ggplot2::geom_point(
      data = pcavals,
      ggplot2::aes_string(x = x, y = y),
      fill = pcacols,
      col = "black",
      pch = 21,
      size = 3
    ) +
    ggplot2::geom_text(
      data = varpc,
      ggplot2::aes(x = v1, y = v2, label = varnames),
      size = 4,
      vjust = 1
    ) +
    ggplot2::geom_segment(
      data = varpc,
      ggplot2::aes(x = 0, y = 0, xend = v1, yend = v2),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      aspect.ratio = 1
    )
  if (display_axes == FALSE) {
    plot <- plot + ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }
  print(plot)
}
utils::assignInNamespace("gdm_plot_vars", gdm_plot_vars_fixed, ns = "algatr")

library(terra)
library(viridis)
library(dplyr)

set.seed(1234)

._lgpr <- Sys.getenv("LGP_PROJECT_ROOT", "~/Desktop/LandscapeGenomicsPipeline")
options(lgp.project_root = sub("/+$", "", path.expand(getOption("lgp.project_root", ._lgpr))))
suppressPackageStartupMessages(base::source(
  base::file.path(getOption("lgp.project_root"), "Scripts", "lgp_pipeline_cache.R"),
  encoding = "UTF-8"
))
rm(._lgpr)

## ------------------------------------------------
## 0. Checks and shared objects
## ------------------------------------------------

data_dir <- lgp_project_root()
in_dir <- lgp_preprocess_dir()
out_dir <- lgp_outputs_step_dir("03-ibd-ibe")
plot_dir <- file.path(out_dir, "plots")
results_dir <- file.path(out_dir, "results")
step_cache_dir <- file.path(out_dir, "_step_cache")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(step_cache_dir, recursive = TRUE, showWarnings = FALSE)

load_if_missing <- function(obj, rds_name) {
  if (!exists(obj, inherits = TRUE)) {
    p <- file.path(in_dir, rds_name)
    if (!file.exists(p)) stop("Missing required input: ", obj, " (expected ", p, ")")
    assign(obj, lgp_read_rds(p), envir = .GlobalEnv)
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
mmrr_deps <- file.path(in_dir, c("euc_gendist.rds", "env_scaled.rds", "coord.rds"))
mmrr_cache <- file.path(step_cache_dir, "mmrr_YX_models.rds")
mmrr_bundle <- lgp_list_cached(
  mmrr_cache,
  dep_paths = mmrr_deps,
  meta      = list(nperm = 999L, stdz = TRUE),
  label     = "MMRR (full + best + predictor matrices)",
  compute_list = function() {
    Ym <- as.matrix(euc_gendist)
    Xm <- env_dist(env_scaled)
    Xm[["geodist"]] <- geo_dist(coord_xy)
    set.seed(10)
    rf <- mmrr_run(Ym, Xm, nperm = 999, stdz = TRUE, model = "full")
    set.seed(11)
    rb <- mmrr_run(Ym, Xm, nperm = 999, stdz = TRUE, model = "best")
    list(Y = Ym, X = Xm, results_full = rf, results_best = rb)
  }
)
Y <- mmrr_bundle$Y
X <- mmrr_bundle$X
results_full <- mmrr_bundle$results_full
results_best <- mmrr_bundle$results_best

# Plots – save to disk
png(file.path(plot_dir, "mmrr_all_full.png"), width = 1600, height = 1200, res = 200)
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "all", stdz = TRUE)
dev.off()

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

png(file.path(plot_dir, "mmrr_vars_best.png"), width = 1600, height = 1200, res = 200)
mmrr_plot(Y, X, mod = results_best$mod, plot_type = "vars", stdz = TRUE)
dev.off()

png(file.path(plot_dir, "mmrr_cov_best.png"), width = 1600, height = 1200, res = 200)
mmrr_plot(Y, X, mod = results_best$mod, plot_type = "cov", stdz = TRUE)
dev.off()

# Fitted vs observed (best model only): 1:1 line = perfect agreement of predicted and observed distance
png(file.path(plot_dir, "mmrr_fitted_best.png"), width = 1600, height = 1200, res = 200)
print(
  algatr::mmrr_plot_fitted(results_best$mod, Y, X, stdz = TRUE) +
    ggplot2::geom_abline(
      intercept = 0,
      slope     = 1,
      linetype  = "dashed",
      colour    = "grey50",
      linewidth = 0.6
    )
)
dev.off()

# Number of population comparisons
n_individuals <- nrow(Y)
n_comparisons <- n_individuals * (n_individuals - 1) / 2
print(n_comparisons)

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
#
# Reduce uninformative / NA pairwise distances before GDM:
# - Keep loci genotyped in at least `min_pops_per_locus` populations.
# - Keep populations with at least `min_loci_per_pop` non-NA frequency estimates.
# Remaining NA pairs (no overlapping loci between two pops) can be handled with
# `pairwise_na_fill`: "error" stops and writes a report; "mean"/etc. imputes (legacy).

min_pops_per_locus <- 2L   # increase (e.g. 3–5) to require broader sharing across pops
min_loci_per_pop   <- 50L  # lower if small panels; set 0L to disable pop-level filter
pairwise_na_fill   <- "error"  # "error" | "mean" | "max" | "min" | "zero"

`%||%` <- function(a, b) if (!is.null(a)) a else b

filter_allele_freqs_for_gdm <- function(P, min_pops_per_locus, min_loci_per_pop) {
  P <- as.matrix(P)
  keep_loci <- colSums(!is.na(P)) >= min_pops_per_locus
  P <- P[, keep_loci, drop = FALSE]
  if (ncol(P) == 0L) {
    stop("No loci left after min_pops_per_locus filter; reduce min_pops_per_locus.")
  }
  if (is.numeric(min_loci_per_pop) && min_loci_per_pop > 0L) {
    keep_pop <- rowSums(!is.na(P)) >= min_loci_per_pop
    if (sum(keep_pop) < 3L) {
      stop(
        "Fewer than 3 populations after min_loci_per_pop filter (GDM needs ≥3 sites). ",
        "Lower min_loci_per_pop or min_pops_per_locus."
      )
    }
    P <- P[keep_pop, , drop = FALSE]
  }
  P
}

calc_euclid_genetic_matrix_complete <- function(allele_freqs,
                                                sqrt = TRUE,
                                                fill = c("error", "mean", "max", "min", "zero"),
                                                as_dist = FALSE,
                                                na_pair_report = NULL) {
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

  na_upper <- is.na(mat) & upper.tri(mat)
  if (any(na_upper)) {
    if (!is.null(na_pair_report)) {
      ij <- which(na_upper, arr.ind = TRUE)
      utils::write.table(
        data.frame(
          pop_i = rownames(mat)[ij[, 1]],
          pop_j = colnames(mat)[ij[, 2]]
        ),
        file = na_pair_report,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
    if (identical(fill, "error")) {
      stop(
        "Pairwise genetic distance undefined (no shared non-NA loci) for some population pairs. ",
        "See ", na_pair_report %||% "(no report path)", " if written. ",
        "Tighten filters (min_pops_per_locus, min_loci_per_pop) or set pairwise_na_fill to ",
        "\"mean\" / \"max\" / \"min\" / \"zero\" (imputation is only a fallback)."
      )
    }
    vals <- mat[upper.tri(mat)]
    vals <- vals[!is.na(vals)]
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

allele_freqs_gdm <- filter_allele_freqs_for_gdm(
  allele_freqs,
  min_pops_per_locus = min_pops_per_locus,
  min_loci_per_pop   = min_loci_per_pop
)

write.table(
  data.frame(
    n_pops          = nrow(allele_freqs_gdm),
    n_loci          = ncol(allele_freqs_gdm),
    min_pops_per_locus = min_pops_per_locus,
    min_loci_per_pop   = min_loci_per_pop,
    pairwise_na_fill   = pairwise_na_fill
  ),
  file = file.path(results_dir, "gdm_allele_freq_filters.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

gendist_euc_mat <- calc_euclid_genetic_matrix_complete(
  allele_freqs_gdm,
  fill = pairwise_na_fill,
  na_pair_report = file.path(results_dir, "gdm_gendist_na_pairs.txt")
)

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

# Align with populations retained in allele_freqs_gdm / gendist_euc_mat
pops_kept <- rownames(gendist_euc_mat)
coord_pop <- coord_pop %>%
  dplyr::filter(.data$population %in% pops_kept) %>%
  dplyr::slice(match(pops_kept, .data$population))

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

set.seed(1234)
varimp <- gdm.varImp(gdmData, geo = TRUE, nPerm = 999)

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
