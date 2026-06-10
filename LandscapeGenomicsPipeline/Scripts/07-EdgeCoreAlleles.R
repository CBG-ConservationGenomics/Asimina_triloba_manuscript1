#######################################################################
# Edge vs core: per-population allele frequencies, private / near-private
# alleles, site-wise Weir-Cockerham Fst, and candidate SNP tables.
#
# Follows the landscape-genomics "part 2" workflow (freqs -> private ->
# edge/core Fst -> outlier threshold -> candidate summary; optional GEA overlap).
#
# Prerequisites: Scripts/00-Preprocessing.R (str_dos, coord, snp_map).
# InputData/Pawpaw_range_Little_1977.zip for default edge/core (species_range).
# Optional: Scripts/04-GEA.R for LFMM / RDA SNP lists.
#######################################################################

library(algatr)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
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

## ---- User-tunable parameters (override via options(lgp.edge_core.*)) ----

lgp_opt <- function(name, default) {
  v <- getOption(paste0("lgp.edge_core.", name), default)
  if (is.null(v) || length(v) != 1L) default else v
}

EDGE_CORE_METHOD   <- lgp_opt("method", "manual")  # species_range | northern_limit | pop_margin | ecd | tess | manual
EDGE_QUANTILE      <- as.numeric(lgp_opt("edge_quantile", 0.75))
## species_range: edge = closest (1 - EDGE_QUANTILE) fraction to Little (1977) range boundary
## northern_limit: edge = highest-latitude fraction (northern range margin)
## pop_margin: edge = farthest EDGE_QUANTILE fraction from sample centroid
LAT_EDGE_QUANTILE  <- as.numeric(lgp_opt("lat_quantile", EDGE_QUANTILE))
LAT_EDGE_SIDE      <- tolower(trimws(as.character(lgp_opt("lat_edge_side", "north"))))
if (!LAT_EDGE_SIDE %in% c("north", "south")) {
  stop("options(lgp.edge_core.lat_edge_side) must be 'north' or 'south', got: ", LAT_EDGE_SIDE)
}
## NULL = test private alleles in every edge population (default).
## Set to "northern" / "southern" or a list file only to restrict the scan.
PRIVATE_EDGE_SUBSET <- lgp_opt("private_edge_subset", NULL)
if (!is.null(PRIVATE_EDGE_SUBSET) && !nzchar(trimws(as.character(PRIVATE_EDGE_SUBSET)))) {
  PRIVATE_EDGE_SUBSET <- NULL
}
ECD_EDGE_SIDE      <- lgp_opt("ecd_edge_side", "east")    # east | west (relative to ECD line)
TESS_EDGE_CLUSTER  <- lgp_opt("tess_edge_cluster", "K2")  # column in tess_qmatrix.tsv

EDGE_FREQ_HIGH     <- as.numeric(lgp_opt("edge_freq_high", 0.10))
CORE_FREQ_LOW      <- as.numeric(lgp_opt("core_freq_low", 0.05))
OTHER_POP_FREQ_MAX <- as.numeric(lgp_opt("other_pop_freq_max", CORE_FREQ_LOW))
STRICT_PRIVATE_MAX <- as.numeric(lgp_opt("strict_other_max", 0.0))  # max freq elsewhere (strict private)
## private_scope: per_edge_pop (default) | pooled_edge | both
PRIVATE_SCOPE      <- tolower(trimws(as.character(lgp_opt("private_scope", "per_edge_pop"))))
if (!PRIVATE_SCOPE %in% c("per_edge_pop", "pooled_edge", "both")) {
  stop(
    "options(lgp.edge_core.private_scope) must be 'per_edge_pop', 'pooled_edge', or 'both', got: ",
    PRIVATE_SCOPE
  )
}
## "near": enriched in focus pop, rare elsewhere; "strict": present in only one population
PRIVATE_MODE       <- tolower(trimws(as.character(lgp_opt("private_mode", "near"))))
if (!PRIVATE_MODE %in% c("near", "strict")) {
  stop("options(lgp.edge_core.private_mode) must be 'near' or 'strict', got: ", PRIVATE_MODE)
}

FST_QUANTILE       <- as.numeric(lgp_opt("fst_quantile", 0.95))
DELTA_P_MIN        <- as.numeric(lgp_opt("delta_p_min", 0.5))
CAND_EDGE_FIX      <- as.numeric(lgp_opt("cand_edge_freq", 0.8))
CAND_CORE_RARE     <- as.numeric(lgp_opt("cand_core_freq", 0.1))

MIN_POPS_PER_LOCUS <- as.integer(lgp_opt("min_pops_per_locus", 2L))
MIN_INDS_PER_GROUP <- as.integer(lgp_opt("min_inds_per_group", 10L))

INTERSECT_GEA      <- isTRUE(lgp_opt("intersect_gea", TRUE))
MAP_GEA_SNPS       <- isTRUE(lgp_opt("map_gea_snps", TRUE))
MAP_GEA_MAX_SNPS   <- as.integer(lgp_opt("map_gea_max_snps", 30L))

## ---- Paths ----------------------------------------------------------------

data_dir    <- lgp_project_root()
in_dir      <- lgp_preprocess_dir()
out_dir     <- lgp_outputs_step_dir("07-edge-core-alleles")
results_dir <- file.path(out_dir, "results")
plot_dir    <- file.path(out_dir, "plots")
cache_dir   <- file.path(out_dir, "_step_cache")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

load_if_missing <- function(obj, rds_name) {
  if (!exists(obj, inherits = TRUE)) {
    p <- file.path(in_dir, rds_name)
    if (!file.exists(p)) {
      stop("Missing `", obj, "` — run Scripts/00-Preprocessing.R (expected ", p, ")")
    }
    assign(obj, readRDS(p), envir = .GlobalEnv)
  }
}

load_if_missing("coord", "coord.rds")
load_if_missing("str_dos", "str_dos.rds")
load_if_missing("snp_map", "snp_map.rds")

sample_names <- rownames(str_dos)
if (is.null(sample_names) || !length(sample_names)) {
  sample_names <- as.character(seq_len(nrow(str_dos)))
  rownames(str_dos) <- sample_names
}
if (nrow(coord) != length(sample_names)) {
  stop("coord rows (", nrow(coord), ") must match str_dos individuals (", length(sample_names), ").")
}

## ---- Population IDs (same convention as 03-IBD-IBE.R) --------------------

pops <- sub("-.*", "", sample_names)
pops <- sub("WC_PA_.*", "WC_PA", pops)
if (!all(c("x", "y") %in% names(coord))) {
  stop("`coord` must contain columns x and y.")
}
popmap <- data.frame(
  sample_id  = sample_names,
  population = pops,
  x          = coord$x,
  y          = coord$y,
  stringsAsFactors = FALSE
)
pop_assignments <- popmap$population[match(rownames(str_dos), popmap$sample_id)]

## ---- Helpers: allele frequencies & Fst -------------------------------------

calc_pop_allele_freq <- function(dosage_matrix, pop_assignments) {
  stopifnot(nrow(dosage_matrix) == length(pop_assignments))
  pops_u <- unique(pop_assignments)
  pop_mm <- model.matrix(~ 0 + factor(pop_assignments, levels = pops_u))
  colnames(pop_mm) <- pops_u
  obs  <- !is.na(dosage_matrix)
  dos0 <- replace(dosage_matrix, is.na(dosage_matrix), 0)
  allele_sums <- t(pop_mm) %*% dos0
  n_obs       <- t(pop_mm) %*% obs
  allele_freqs <- allele_sums / (2 * n_obs)
  allele_freqs[is.na(allele_freqs)] <- NA_real_
  rownames(allele_freqs) <- pops_u
  allele_freqs
}

weir_cockerham_fst_loci <- function(dosage, group) {
  g <- factor(group, levels = c("edge", "core"))
  if (!all(levels(g) %in% levels(g)[table(g) > 0])) {
    stop("Both edge and core groups need at least one individual.")
  }
  e <- g == "edge"
  c <- g == "core"
  G_e <- dosage[e, , drop = FALSE]
  G_c <- dosage[c, , drop = FALSE]
  n1 <- nrow(G_e)
  n2 <- nrow(G_c)
  p1 <- colMeans(G_e, na.rm = TRUE) / 2
  p2 <- colMeans(G_c, na.rm = TRUE) / 2
  h1 <- colSums(!is.na(G_e))
  h2 <- colSums(!is.na(G_c))
  n  <- h1 + h2
  p  <- (h1 * p1 + h2 * p2) / n
  v  <- (h1 * (p1 - p)^2 + h2 * (p2 - p)^2) / n
  denom <- p * (1 - p)
  fst <- v / denom
  fst[!is.finite(fst) | denom <= 0] <- NA_real_
  fst
}

write_sample_list <- function(samples, path) {
  writeLines(samples, path)
  invisible(path)
}

## Per edge population: private / near-private vs all other pops (edge + core).
calc_edge_pop_private_long <- function(
    af_mat,
    edge_pops,
    core_pops,
    private_mode = PRIVATE_MODE,
    edge_freq_high = EDGE_FREQ_HIGH,
    other_pop_freq_max = OTHER_POP_FREQ_MAX,
    strict_other_max = STRICT_PRIVATE_MAX
) {
  edge_pops <- intersect(edge_pops, rownames(af_mat))
  core_pops <- intersect(core_pops, rownames(af_mat))
  if (!length(edge_pops)) {
    return(data.frame())
  }
  out <- vector("list", length(edge_pops))
  for (ep in edge_pops) {
    f_ep <- af_mat[ep, , drop = TRUE]
    other_rows <- setdiff(rownames(af_mat), ep)
    other_edge <- setdiff(edge_pops, ep)
    if (length(other_rows)) {
      max_other <- apply(af_mat[other_rows, , drop = FALSE], 2, max, na.rm = TRUE)
    } else {
      max_other <- rep(NA_real_, ncol(af_mat))
    }
    if (length(other_edge)) {
      max_other_edge <- apply(af_mat[other_edge, , drop = FALSE], 2, max, na.rm = TRUE)
    } else {
      max_other_edge <- rep(0, ncol(af_mat))
    }
    if (length(core_pops)) {
      max_core <- apply(af_mat[core_pops, , drop = FALSE], 2, max, na.rm = TRUE)
    } else {
      max_core <- rep(NA_real_, ncol(af_mat))
    }
    strict_private <- vapply(
      seq_along(f_ep),
      function(j) {
        v <- af_mat[, j]
        nz <- which(!is.na(v) & v > strict_other_max)
        length(nz) == 1L && rownames(af_mat)[nz] == ep
      },
      FUN.VALUE = logical(1)
    )
    near_private <- !is.na(f_ep) & (f_ep >= edge_freq_high) &
      !is.na(max_other) & (max_other <= other_pop_freq_max)
    private_flag <- if (private_mode == "strict") strict_private else near_private
    out[[ep]] <- data.frame(
      snp = colnames(af_mat),
      edge_population = ep,
      freq_edge_pop = as.numeric(f_ep),
      max_freq_other_pops = as.numeric(max_other),
      max_freq_other_edge = as.numeric(max_other_edge),
      max_freq_core = as.numeric(max_core),
      strict_private = strict_private,
      near_private = near_private,
      private_in_edge_pop = private_flag,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(out)
}

summarize_edge_pop_private <- function(edge_pop_long) {
  if (!nrow(edge_pop_long)) {
    return(data.frame(
      snp = character(),
      n_edge_pops_private = integer(),
      edge_pops_private = character(),
      private_any_edge_pop = logical(),
      stringsAsFactors = FALSE
    ))
  }
  base <- data.frame(snp = unique(edge_pop_long$snp), stringsAsFactors = FALSE)
  hits <- edge_pop_long[edge_pop_long$private_in_edge_pop, , drop = FALSE]
  if (!nrow(hits)) {
    base$n_edge_pops_private <- 0L
    base$edge_pops_private <- NA_character_
    base$private_any_edge_pop <- FALSE
    return(base)
  }
  agg <- hits %>%
    dplyr::group_by(.data$snp) %>%
    dplyr::summarise(
      n_edge_pops_private = dplyr::n(),
      edge_pops_private = paste(sort(unique(.data$edge_population)), collapse = ";"),
      private_any_edge_pop = TRUE,
      .groups = "drop"
    )
  dplyr::left_join(base, agg, by = "snp") %>%
    dplyr::mutate(
      n_edge_pops_private = dplyr::coalesce(.data$n_edge_pops_private, 0L),
      private_any_edge_pop = dplyr::coalesce(.data$private_any_edge_pop, FALSE)
    )
}

## Subset edge populations for private-allele tests (e.g. northern limit only).
resolve_private_edge_pops <- function(
    edge_pops,
    centroids_df,
    subset_opt,
    data_dir_root,
    popmap_df = popmap
) {
  if (is.null(subset_opt) || !nzchar(trimws(as.character(subset_opt)))) {
    return(edge_pops)
  }
  subset_opt <- tolower(trimws(as.character(subset_opt)))
  edge_pops <- intersect(edge_pops, centroids_df$population)
  if (!length(edge_pops)) return(character())
  y_edge <- centroids_df$y[match(edge_pops, centroids_df$population)]
  if (subset_opt %in% c("northern", "north")) {
    thr <- stats::quantile(y_edge, LAT_EDGE_QUANTILE, na.rm = TRUE)
    return(edge_pops[!is.na(y_edge) & y_edge >= thr])
  }
  if (subset_opt %in% c("southern", "south")) {
    thr <- stats::quantile(y_edge, 1 - LAT_EDGE_QUANTILE, na.rm = TRUE)
    return(edge_pops[!is.na(y_edge) & y_edge <= thr])
  }
  path <- subset_opt
  if (!file.exists(path)) {
    path <- file.path(data_dir_root, "InputData", subset_opt)
  }
  if (file.exists(path)) {
    listed <- resolve_manual_populations(
      readLines(path, warn = FALSE),
      popmap_df$sample_id,
      popmap_df$population
    )
    return(intersect(edge_pops, listed))
  }
  stop(
    "Unknown lgp.edge_core.private_edge_subset: ", subset_opt,
    " (use 'northern', 'southern', or a population list file under InputData/)"
  )
}

resolve_manual_populations <- function(lines, sample_ids, populations) {
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) return(character())
  pops_u <- unique(populations)
  out <- lines[lines %in% pops_u]
  si <- match(lines, sample_ids)
  out <- c(out, populations[si[!is.na(si)]])
  derived <- sub("-.*", "", lines)
  derived <- sub("WC_PA_.*", "WC_PA", derived)
  out <- c(out, derived[derived %in% pops_u])
  sort(unique(out[!is.na(out) & nzchar(out)]))
}

read_shp_from_zip <- function(zip_path) {
  if (!file.exists(zip_path)) stop("Zip not found: ", zip_path)
  zl <- utils::unzip(zip_path, list = TRUE)
  shp <- zl$Name[grepl("\\.shp$", zl$Name, ignore.case = TRUE)][1]
  if (is.na(shp) || !nzchar(shp)) stop("No .shp file inside ", zip_path)
  zip_norm <- normalizePath(zip_path, winslash = "/", mustWork = TRUE)
  sf::read_sf(paste0("/vsizip/", zip_norm, "/", shp), quiet = TRUE)
}

ensure_lonlat_crs <- function(x, label) {
  x <- sf::st_zm(x, drop = TRUE, what = "ZM")
  if (is.na(sf::st_crs(x))) {
    warning("`", label, "` has no CRS; assuming WGS84 (EPSG:4326).")
    sf::st_crs(x) <- 4326
  }
  x
}

pick_utm_epsg <- function(lon, lat) {
  zone <- floor((lon + 180) / 6) + 1
  if (lat >= 0) 32600 + zone else 32700 + zone
}

load_pawpaw_range_sf <- function(root = data_dir) {
  pawpaw_zip <- file.path(root, "InputData", "Pawpaw_range_Little_1977.zip")
  if (!file.exists(pawpaw_zip)) {
    return(NULL)
  }
  ensure_lonlat_crs(
    read_shp_from_zip(pawpaw_zip),
    "Pawpaw range (Little 1977)"
  )
}

## Base map: species range (grey underlay, same as 01-GeneticDiversity / 02-PopulationStructure).
add_ec_map_basemap <- function(p, pawpaw_range_sf, states_sf) {
  if (!is.null(pawpaw_range_sf)) {
    p <- p + ggplot2::geom_sf(
      data  = pawpaw_range_sf,
      fill  = "grey90",
      alpha = 0.4,
      colour = "grey60",
      linewidth = 0.4,
      inherit.aes = FALSE
    )
  }
  p + ggplot2::geom_sf(
    data  = states_sf,
    fill  = NA,
    colour = "grey45",
    linewidth = 0.28,
    inherit.aes = FALSE
  )
}

dist_pop_centroids_to_range_km <- function(centroids_df, range_sf_ll) {
  range_sf_ll <- ensure_lonlat_crs(range_sf_ll, "Pawpaw range (Little 1977)")
  crs_ll <- sf::st_crs(range_sf_ll)
  range_union <- sf::st_union(sf::st_make_valid(sf::st_geometry(range_sf_ll)))
  sf::st_crs(range_union) <- crs_ll
  range_boundary <- sf::st_boundary(range_union)
  sf::st_crs(range_boundary) <- crs_ll

  pts <- sf::st_as_sf(centroids_df, coords = c("x", "y"), crs = crs_ll, remove = FALSE)
  utm <- pick_utm_epsg(mean(centroids_df$x, na.rm = TRUE), mean(centroids_df$y, na.rm = TRUE))
  pts_p <- sf::st_transform(pts, utm)
  bound_p <- sf::st_transform(range_boundary, utm)
  d_km <- as.numeric(sf::st_distance(pts_p, bound_p)) / 1000
  inside <- sf::st_within(pts, range_union, sparse = FALSE)[, 1]
  data.frame(
    population = centroids_df$population,
    dist_to_range_boundary_km = d_km,
    inside_species_range = inside,
    stringsAsFactors = FALSE
  )
}

plot_gea_candidate_freq_maps <- function(results_dir,
                                         plot_dir,
                                         af_long,
                                         pop_centroids,
                                         pop_assignment_tbl,
                                         max_snps = 30L,
                                         pawpaw_range_sf = NULL) {
  gea_map_path <- file.path(results_dir, "edge_core_candidates_x_gea.tsv")
  gea_map_dir  <- file.path(plot_dir, "gea_candidate_freq_maps")
  dir.create(gea_map_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(gea_map_path)) {
    message("[07-edge-core] No ", basename(gea_map_path), "; skipping GEA frequency maps.")
    return(invisible(NULL))
  }

  gea_map_snps <- utils::read.delim(gea_map_path, check.names = FALSE)
  if (!"snp" %in% names(gea_map_snps)) {
    stop("edge_core_candidates_x_gea.tsv must contain a `snp` column.")
  }
  if ("fst_edge_core" %in% names(gea_map_snps)) {
    gea_map_snps <- gea_map_snps[order(-gea_map_snps$fst_edge_core), , drop = FALSE]
  }
  snp_plot <- unique(gea_map_snps$snp)
  if (length(snp_plot) > max_snps) {
    message(
      "[07-edge-core] Plotting top ", max_snps, " GEA SNPs by Fst (of ",
      length(snp_plot), " total)."
    )
    snp_plot <- snp_plot[seq_len(max_snps)]
  }

  east_coast_abbrev <- c(
    "ME", "NH", "MA", "RI", "CT", "NY", "NJ", "DE", "MD", "VA", "NC", "SC", "GA", "FL",
    "AL", "LA", "TN", "MS", "AR", "KY", "OH", "MI", "IN", "IL", "MO", "IA", "PA", "TX", "OK", "KS"
  )
  us_states_ec <- rnaturalearth::ne_states(
    country = "United States of America",
    returnclass = "sf"
  ) %>%
    dplyr::filter(.data$postal %in% east_coast_abbrev)
  map_bb <- sf::st_bbox(us_states_ec)

  if (is.null(pawpaw_range_sf)) {
    pawpaw_range_sf <- load_pawpaw_range_sf(lgp_project_root())
  }
  if (is.null(pawpaw_range_sf)) {
    warning(
      "[07-edge-core] Pawpaw_range_Little_1977.zip not found; maps omit species range polygon."
    )
  }
  pawpaw_caption <- "Grey fill = Little (1977) Asimina triloba range (Pawpaw_range_Little_1977.zip)."

  freq_map_df <- af_long %>%
    dplyr::filter(.data$snp %in% snp_plot) %>%
    dplyr::left_join(pop_centroids[, c("population", "x", "y", "n_samples")], by = "population")
  if (!"edge_core" %in% names(freq_map_df)) {
    freq_map_df <- dplyr::left_join(
      freq_map_df,
      pop_assignment_tbl[, c("population", "edge_core")],
      by = "population"
    )
  }
  freq_map_df <- freq_map_df %>%
    dplyr::left_join(
      gea_map_snps[, intersect(names(gea_map_snps), c(
        "snp", "fst_edge_core", "mean_edge_freq", "mean_core_freq",
        "delta_p_edge_minus_core", "chr", "pos"
      )), drop = FALSE],
      by = "snp"
    )
  freq_map_sf <- sf::st_as_sf(freq_map_df, coords = c("x", "y"), crs = 4326, remove = FALSE)

  freq_map_theme <- function() {
    ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "grey98", colour = NA),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.title = ggplot2::element_text(face = "bold", size = 10),
        plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35"),
        legend.title = ggplot2::element_text(face = "bold")
      )
  }

  snp_label <- function(snp_id, row) {
    chr <- if ("chr" %in% names(row)) row$chr[1] else NA
    pos <- if ("pos" %in% names(row)) row$pos[1] else NA
    fst <- if ("fst_edge_core" %in% names(row)) row$fst_edge_core[1] else NA
    sub <- paste0(
      if (!is.na(chr)) paste0(chr, ":", pos) else "",
      if (!is.na(fst)) paste0("  Fst=", signif(fst, 3)) else ""
    )
    paste(snp_id, sub)
  }

  ## shape = 21 (filled circle): geom_sf points ignore fill for shapes 16/17
  pop_pt_shape <- 21L

  for (one_snp in snp_plot) {
    d1 <- freq_map_sf[freq_map_sf$snp == one_snp, , drop = FALSE]
    if (!nrow(d1)) next
    p_map <- add_ec_map_basemap(ggplot2::ggplot(), pawpaw_range_sf, us_states_ec) +
      ggplot2::geom_sf(
        data = d1,
        ggplot2::aes(
          fill    = .data$alt_allele_freq,
          colour  = .data$edge_core,
          size    = .data$n_samples
        ),
        shape  = pop_pt_shape,
        stroke = 0.35,
        alpha  = 0.92
      ) +
      ggplot2::scale_fill_viridis_c(
        option = "viridis", limits = c(0, 1), na.value = "grey80",
        name = "Alt allele\nfrequency"
      ) +
      ggplot2::scale_colour_manual(
        values = c(edge = "#C45A00", core = "grey20"),
        name = "Population\ngroup"
      ) +
      ggplot2::scale_size_continuous(range = c(3, 9), name = "Samples\nper site") +
      ggplot2::coord_sf(
        xlim = map_bb[c("xmin", "xmax")],
        ylim = map_bb[c("ymin", "ymax")],
        expand = FALSE
      ) +
      ggplot2::labs(
        title = "Per-population allele frequency",
        subtitle = snp_label(one_snp, d1),
        caption = paste(
          "Filled circles = site populations (fill = alt-allele freq; border = edge/core).",
          pawpaw_caption
        )
      ) +
      freq_map_theme()
    safe_name <- gsub("[^A-Za-z0-9._-]+", "_", one_snp)
    ggplot2::ggsave(
      file.path(gea_map_dir, paste0("freq_map_", safe_name, ".png")),
      plot = p_map, width = 9, height = 7, dpi = 300
    )
  }

  if (length(snp_plot) > 1L) {
    p_fac <- add_ec_map_basemap(ggplot2::ggplot(), pawpaw_range_sf, us_states_ec) +
      ggplot2::geom_sf(
        data = freq_map_sf,
        ggplot2::aes(
          fill   = .data$alt_allele_freq,
          colour = .data$edge_core,
          size   = .data$n_samples
        ),
        shape  = pop_pt_shape,
        stroke = 0.3,
        alpha  = 0.92
      ) +
      ggplot2::scale_fill_viridis_c(option = "viridis", limits = c(0, 1), name = "Alt freq") +
      ggplot2::scale_colour_manual(
        values = c(edge = "#C45A00", core = "grey20"),
        name = "Group"
      ) +
      ggplot2::scale_size_continuous(range = c(2.5, 7), guide = "none") +
      ggplot2::facet_wrap(~snp, ncol = 2, labeller = ggplot2::label_wrap_gen(width = 28)) +
      ggplot2::coord_sf(
        xlim = map_bb[c("xmin", "xmax")],
        ylim = map_bb[c("ymin", "ymax")],
        expand = FALSE
      ) +
      ggplot2::labs(
        title = "GEA x edge-core candidates: population allele frequencies",
        caption = paste("One panel per SNP; site-level populations.", pawpaw_caption)
      ) +
      freq_map_theme() +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 7), legend.position = "right")
    ggplot2::ggsave(
      file.path(gea_map_dir, "freq_maps_gea_candidates_faceted.png"),
      plot = p_fac,
      width = 12,
      height = max(6, 3 * ceiling(length(snp_plot) / 2)),
      dpi = 300,
      limitsize = FALSE
    )
  }

  message("[07-edge-core] Wrote ", length(snp_plot), " frequency map(s) -> ", gea_map_dir)
  invisible(gea_map_dir)
}

## Map-only: LGP_07_MAP_ONLY=1 Rscript Scripts/07-EdgeCoreAlleles.R
if (identical(Sys.getenv("LGP_07_MAP_ONLY", ""), "1")) {
  af_long <- utils::read.delim(
    file.path(results_dir, "pop_allele_frequencies_long.tsv"), check.names = FALSE
  )
  pop_assignment_tbl <- utils::read.delim(
    file.path(results_dir, "population_edge_core_assignment.tsv"), check.names = FALSE
  )
  pop_centroids <- pop_assignment_tbl[, c("population", "x", "y", "n_samples")]
  plot_gea_candidate_freq_maps(
    results_dir, plot_dir, af_long, pop_centroids, pop_assignment_tbl,
    max_snps = MAP_GEA_MAX_SNPS,
    pawpaw_range_sf = load_pawpaw_range_sf(lgp_project_root())
  )
  quit(save = "no", status = 0)
}

## ---- Edge vs core population assignment ------------------------------------

pop_centroids <- popmap %>%
  dplyr::group_by(.data$population) %>%
  dplyr::summarise(
    x = mean(.data$x, na.rm = TRUE),
    y = mean(.data$y, na.rm = TRUE),
    n_samples = dplyr::n(),
    .groups = "drop"
  )

assign_edge_core_pops <- function(method = EDGE_CORE_METHOD, popmap_df = popmap) {
  method <- match.arg(
    method,
    c("species_range", "northern_limit", "pop_margin", "ecd", "tess", "manual")
  )
  pops_all <- pop_centroids$population
  edge_pops <- character()
  core_pops <- character()
  note <- ""

  if (method == "species_range") {
    pawpaw_zip <- file.path(data_dir, "InputData", "Pawpaw_range_Little_1977.zip")
    if (!file.exists(pawpaw_zip)) {
      stop(
        "method = species_range requires InputData/Pawpaw_range_Little_1977.zip\n",
        "  (Little 1977 Asimina triloba range polygon)."
      )
    }
    asimina_range <- ensure_lonlat_crs(read_shp_from_zip(pawpaw_zip), "Pawpaw range (Little 1977)")
    range_dist <- dist_pop_centroids_to_range_km(pop_centroids, asimina_range)
    thr <- stats::quantile(
      range_dist$dist_to_range_boundary_km,
      probs = 1 - EDGE_QUANTILE,
      na.rm = TRUE
    )
    edge_pops <- range_dist$population[range_dist$dist_to_range_boundary_km <= thr]
    core_pops <- range_dist$population[range_dist$dist_to_range_boundary_km > thr]
    note <- paste0(
      "Edge = populations within ", round(as.numeric(thr), 1),
      " km of Little (1977) range boundary (closest ",
      round(100 * (1 - EDGE_QUANTILE), 1), "%; quantile ", EDGE_QUANTILE, ")"
    )
    range_meta <- list(range_dist = range_dist, asimina_range = asimina_range)
  } else if (method == "manual") {
    edge_file <- file.path(data_dir, "InputData", "edge_populations.txt")
    core_file <- file.path(data_dir, "InputData", "core_populations.txt")
    if (!file.exists(edge_file) || !file.exists(core_file)) {
      stop(
        "method = manual requires:\n  ", edge_file, "\n  ", core_file,
        "\n(one population name per line; see InputData/README_edge_core.txt)"
      )
    }
    edge_pops <- resolve_manual_populations(
      readLines(edge_file, warn = FALSE),
      popmap_df$sample_id,
      popmap_df$population
    )
    core_pops <- resolve_manual_populations(
      readLines(core_file, warn = FALSE),
      popmap_df$sample_id,
      popmap_df$population
    )
    if (!length(edge_pops)) {
      stop(
        "No edge populations matched ", edge_file,
        ". Use site-level population names (e.g. VA_HardR, LA_Hoges) or sample_id lines."
      )
    }
    if (!length(core_pops)) {
      stop(
        "No core populations matched ", core_file,
        ". Use site-level population names (e.g. VA_HardR) or sample_id lines."
      )
    }
    overlap <- intersect(edge_pops, core_pops)
    if (length(overlap)) {
      warning(
        "[07-edge-core] Populations listed as both edge and core (edge wins): ",
        paste(overlap, collapse = ", "),
        call. = FALSE
      )
      core_pops <- setdiff(core_pops, overlap)
    }
    note <- paste0(
      "Manual lists in InputData/ (",
      length(edge_pops), " edge pops, ", length(core_pops), " core pops)"
    )
  } else if (method == "northern_limit") {
    thr <- stats::quantile(pop_centroids$y, LAT_EDGE_QUANTILE, na.rm = TRUE)
    if (LAT_EDGE_SIDE == "north") {
      edge_pops <- pop_centroids$population[pop_centroids$y >= thr]
      core_pops <- pop_centroids$population[pop_centroids$y < thr]
      note <- paste0(
        "Northern range limit: latitude (y) >= ", round(as.numeric(thr), 4),
        " deg (top ", round(100 * (1 - LAT_EDGE_QUANTILE), 1),
        "% of populations; lat_quantile = ", LAT_EDGE_QUANTILE, ")"
      )
    } else {
      thr_s <- stats::quantile(pop_centroids$y, 1 - LAT_EDGE_QUANTILE, na.rm = TRUE)
      edge_pops <- pop_centroids$population[pop_centroids$y <= thr_s]
      core_pops <- pop_centroids$population[pop_centroids$y > thr_s]
      note <- paste0(
        "Southern range limit: latitude (y) <= ", round(as.numeric(thr_s), 4),
        " deg (bottom ", round(100 * (1 - LAT_EDGE_QUANTILE), 1),
        "% of populations; lat_quantile = ", LAT_EDGE_QUANTILE, ")"
      )
    }
    if (!length(edge_pops) || !length(core_pops)) {
      stop("northern_limit split left empty edge or core; adjust lat_quantile.")
    }
  } else if (method == "pop_margin") {
    cx <- mean(pop_centroids$x, na.rm = TRUE)
    cy <- mean(pop_centroids$y, na.rm = TRUE)
    pop_centroids$dist_centroid <- sqrt(
      (pop_centroids$x - cx)^2 + (pop_centroids$y - cy)^2
    )
    thr <- stats::quantile(pop_centroids$dist_centroid, EDGE_QUANTILE, na.rm = TRUE)
    edge_pops <- pop_centroids$population[pop_centroids$dist_centroid >= thr]
    core_pops <- setdiff(pops_all, edge_pops)
    note <- paste0(
      "Populations in top ", 100 * (1 - EDGE_QUANTILE),
      "% distance from range-wide centroid (quantile = ", EDGE_QUANTILE, ")"
    )
  } else if (method == "ecd") {
    ecd_path <- file.path(data_dir, "InputData", "Eastern_watershed.kml")
    if (!file.exists(ecd_path)) stop("ECD KML not found: ", ecd_path)
    med_lon <- stats::median(pop_centroids$x, na.rm = TRUE)
    side_chr <- ifelse(pop_centroids$x >= med_lon, "east", "west")
    note <- paste0(
      "Populations east vs west of median population longitude (",
      round(med_lon, 3), " deg); edge side = ", ECD_EDGE_SIDE
    )
    edge_pops <- pop_centroids$population[side_chr == ECD_EDGE_SIDE]
    core_pops <- pop_centroids$population[side_chr != ECD_EDGE_SIDE]
    if (!length(edge_pops) || !length(core_pops)) {
      stop("ECD split left empty edge or core set; try options(lgp.edge_core.ecd_edge_side).")
    }
  } else if (method == "tess") {
    tess_path <- file.path(
      lgp_outputs_base(), "02-population-structure", "tess_qmatrix.tsv"
    )
    if (!file.exists(tess_path)) {
      stop("TESS Q matrix not found. Run Scripts/02-PopulationStructure.R first: ", tess_path)
    }
    tess_df <- utils::read.delim(tess_path, check.names = FALSE)
    if (!"sample_id" %in% names(tess_df)) {
      names(tess_df)[1] <- "sample_id"
    }
    if (!TESS_EDGE_CLUSTER %in% names(tess_df)) {
      stop("Column ", TESS_EDGE_CLUSTER, " not in ", tess_path)
    }
    tess_df$population <- sub("-.*", "", tess_df$sample_id)
    tess_df$population <- sub("WC_PA_.*", "WC_PA", tess_df$population)
    pop_q <- tess_df %>%
      dplyr::group_by(.data$population) %>%
      dplyr::summarise(
        mean_q = mean(.data[[TESS_EDGE_CLUSTER]], na.rm = TRUE),
        .groups = "drop"
      )
    med <- stats::median(pop_q$mean_q, na.rm = TRUE)
    edge_pops <- pop_q$population[pop_q$mean_q >= med]
    core_pops <- pop_q$population[pop_q$mean_q < med]
    note <- paste0("Populations with mean ", TESS_EDGE_CLUSTER, " >= median (", round(med, 3), ")")
  }

  out <- list(
    edge_populations = sort(unique(edge_pops)),
    core_populations = sort(unique(core_pops)),
    method = method,
    note = note,
    range_dist = NULL,
    asimina_range = NULL
  )
  if (exists("range_meta", inherits = FALSE)) {
    out$range_dist <- range_meta$range_dist
    out$asimina_range <- range_meta$asimina_range
  }
  out
}

edge_core_meta <- assign_edge_core_pops(EDGE_CORE_METHOD)
edge_pops <- edge_core_meta$edge_populations
core_pops <- edge_core_meta$core_populations

pop_assignment_tbl <- pop_centroids %>%
  dplyr::mutate(
    edge_core = dplyr::if_else(
      .data$population %in% edge_pops,
      "edge",
      dplyr::if_else(.data$population %in% core_pops, "core", "unassigned")
    ),
    method = edge_core_meta$method,
    method_note = edge_core_meta$note
  )
if (!is.null(edge_core_meta$range_dist)) {
  pop_assignment_tbl <- dplyr::left_join(
    pop_assignment_tbl,
    edge_core_meta$range_dist,
    by = "population"
  )
}
write.table(
  pop_assignment_tbl,
  file.path(results_dir, "population_edge_core_assignment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pawpaw_range_plot <- edge_core_meta$asimina_range
if (is.null(pawpaw_range_plot)) {
  pawpaw_range_plot <- load_pawpaw_range_sf(data_dir)
}
if (!is.null(pawpaw_range_plot)) {
  pts_ec <- sf::st_as_sf(pop_assignment_tbl, coords = c("x", "y"), crs = 4326, remove = FALSE)
  east_coast_abbrev <- c(
    "ME", "NH", "MA", "RI", "CT", "NY", "NJ", "DE", "MD", "VA", "NC", "SC", "GA", "FL",
    "AL", "LA", "TN", "MS", "AR", "KY", "OH", "MI", "IN", "IL", "MO", "IA", "PA", "TX", "OK", "KS"
  )
  states_ec <- rnaturalearth::ne_states(
    country = "United States of America", returnclass = "sf"
  ) %>%
    dplyr::filter(.data$postal %in% east_coast_abbrev)
  has_range_dist <- "dist_to_range_boundary_km" %in% names(pop_assignment_tbl)
  ec_group_cols <- c(edge = "#C45A00", core = "grey55", unassigned = "#9B59B6")
  p_range <- add_ec_map_basemap(ggplot2::ggplot(), pawpaw_range_plot, states_ec)
  if (has_range_dist) {
    p_range <- p_range +
      ggplot2::geom_sf(
        data = pts_ec,
        ggplot2::aes(
          fill = .data$dist_to_range_boundary_km,
          colour = .data$edge_core,
          size = .data$n_samples
        ),
        shape = 21,
        stroke = 0.35
      ) +
      ggplot2::scale_fill_viridis_c(
        name = "Distance to\nrange (km)", option = "magma", direction = -1
      ) +
      ggplot2::scale_colour_manual(values = ec_group_cols[c("edge", "core")], name = "Group") +
      ggplot2::labs(
        title = "Edge vs core by Little (1977) Asimina triloba range",
        caption = paste(
          "Point fill = km to range boundary; border = edge/core group.",
          "Grey polygon = Little (1977) range (same style as scripts 01–02)."
        )
      )
  } else {
    p_range <- p_range +
      ggplot2::geom_sf(
        data = pts_ec,
        ggplot2::aes(
          fill = .data$edge_core,
          size = .data$n_samples
        ),
        shape = 21,
        colour = "grey15",
        stroke = 0.35
      ) +
      ggplot2::scale_fill_manual(values = ec_group_cols, name = "Group", drop = FALSE) +
      ggplot2::labs(
        title = paste0("Edge vs core populations (", EDGE_CORE_METHOD, ")"),
        caption = "Grey polygon = Little (1977) range underlay when available."
      )
  }
  p_range <- p_range +
    ggplot2::scale_size_continuous(range = c(2.5, 7), name = "Samples") +
    ggplot2::coord_sf(
      xlim = sf::st_bbox(states_ec)[c("xmin", "xmax")],
      ylim = sf::st_bbox(states_ec)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    ggplot2::labs(subtitle = edge_core_meta$note) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(
    file.path(plot_dir, paste0("population_edge_core_", EDGE_CORE_METHOD, ".png")),
    p_range, width = 9, height = 7, dpi = 300
  )
}

popmap$edge_core <- ifelse(
  popmap$population %in% edge_pops, "edge",
  ifelse(popmap$population %in% core_pops, "core", NA_character_)
)
edge_samples <- popmap$sample_id[popmap$edge_core == "edge"]
core_samples <- popmap$sample_id[popmap$edge_core == "core"]
if (length(edge_samples) < MIN_INDS_PER_GROUP || length(core_samples) < MIN_INDS_PER_GROUP) {
  stop(
    "Too few individuals in edge (", length(edge_samples), ") or core (",
    length(core_samples), "); lower min_inds_per_group or change method."
  )
}
write_sample_list(edge_samples, file.path(results_dir, "edge_samples.txt"))
write_sample_list(core_samples, file.path(results_dir, "core_samples.txt"))

message(
  "[07-edge-core] method=", EDGE_CORE_METHOD, ": ",
  length(edge_pops), " edge pops (n=", length(edge_samples), " inds), ",
  length(core_pops), " core pops (n=", length(core_samples), " inds)"
)

## ---- 1. Per-population allele frequencies ----------------------------------

str_dep <- file.path(in_dir, "str_dos.rds")
af_cache <- file.path(cache_dir, "pop_allele_freqs.rds")

allele_freqs <- lgp_rds_cached(
  af_cache,
  dep_paths = str_dep,
  label = "per-population allele frequencies",
  compute = function() calc_pop_allele_freq(str_dos, pop_assignments)
)

## Filter loci with enough population coverage
keep_loci <- colSums(!is.na(allele_freqs)) >= MIN_POPS_PER_LOCUS
allele_freqs <- allele_freqs[, keep_loci, drop = FALSE]
snps <- colnames(allele_freqs)

af_long <- as.data.frame(allele_freqs)
af_long$population <- rownames(af_long)
rownames(af_long) <- NULL
af_long <- tidyr::pivot_longer(
  af_long,
  cols = -"population",
  names_to = "snp",
  values_to = "alt_allele_freq"
)
af_long <- dplyr::left_join(af_long, pop_assignment_tbl[, c("population", "edge_core")], by = "population")
write.table(
  af_long,
  file.path(results_dir, "pop_allele_frequencies_long.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## Wide: one row per SNP, one column per population (allele_freqs is pops x SNPs)
af_wide <- data.frame(snp = snps, t(allele_freqs), check.names = FALSE)
write.table(
  af_wide,
  file.path(results_dir, "pop_allele_frequencies_wide.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## ---- 2. Private & near-private alleles (population level) ----------------

af_mat <- as.matrix(allele_freqs)
edge_rows <- rownames(af_mat) %in% edge_pops
core_rows <- rownames(af_mat) %in% core_pops

max_edge_freq <- apply(af_mat[edge_rows, , drop = FALSE], 2, max, na.rm = TRUE)
mean_edge_freq <- apply(af_mat[edge_rows, , drop = FALSE], 2, mean, na.rm = TRUE)
max_core_freq <- apply(af_mat[core_rows, , drop = FALSE], 2, max, na.rm = TRUE)
mean_core_freq <- apply(af_mat[core_rows, , drop = FALSE], 2, mean, na.rm = TRUE)
min_core_freq <- apply(af_mat[core_rows, , drop = FALSE], 2, min, na.rm = TRUE)

n_pops_nonzero <- colSums(af_mat > STRICT_PRIVATE_MAX, na.rm = TRUE)
strict_private_pop <- n_pops_nonzero == 1L & max_edge_freq > STRICT_PRIVATE_MAX
strict_private_edge <- apply(af_mat, 2, function(v) {
  nz <- which(!is.na(v) & v > STRICT_PRIVATE_MAX)
  if (length(nz) != 1L) return(FALSE)
  rownames(af_mat)[nz] %in% edge_pops
})
strict_private_core <- apply(af_mat, 2, function(v) {
  nz <- which(!is.na(v) & v > STRICT_PRIVATE_MAX)
  if (length(nz) != 1L) return(FALSE)
  rownames(af_mat)[nz] %in% core_pops
})

near_private_edge_pooled <- (max_edge_freq >= EDGE_FREQ_HIGH) & (max_core_freq <= CORE_FREQ_LOW)
quasi_private_edge_pooled <- near_private_edge_pooled & !strict_private_edge
private_edge_pooled <- if (PRIVATE_MODE == "strict") {
  strict_private_edge
} else {
  near_private_edge_pooled
}

private_edge_pops <- resolve_private_edge_pops(
  edge_pops, pop_centroids, PRIVATE_EDGE_SUBSET, data_dir
)
if (!length(private_edge_pops)) {
  stop(
    "No edge populations left for private-allele scan after private_edge_subset=",
    PRIVATE_EDGE_SUBSET, ". Check edge assignment and subset."
  )
}
subset_lab <- if (is.null(PRIVATE_EDGE_SUBSET)) {
  "all edge populations"
} else {
  paste0("subset=", PRIVATE_EDGE_SUBSET)
}
message(
  "[07-edge-core] Private-allele scan (per edge pop, ", subset_lab, "): ",
  length(private_edge_pops), " populations"
)
write.table(
  data.frame(
    population = private_edge_pops,
    latitude = pop_centroids$y[match(private_edge_pops, pop_centroids$population)],
    stringsAsFactors = FALSE
  ),
  file.path(results_dir, "populations_used_for_private_scan.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
edge_pop_private_long <- calc_edge_pop_private_long(
  af_mat, private_edge_pops, core_pops
)
edge_pop_private_summary <- summarize_edge_pop_private(edge_pop_private_long)

private_df <- data.frame(
  snp = snps,
  max_edge_freq = max_edge_freq,
  mean_edge_freq = mean_edge_freq,
  max_core_freq = max_core_freq,
  mean_core_freq = mean_core_freq,
  delta_p_edge_minus_core = mean_edge_freq - mean_core_freq,
  n_pops_with_allele = n_pops_nonzero,
  private_scope = PRIVATE_SCOPE,
  private_mode = PRIVATE_MODE,
  strict_private_any_pop = strict_private_pop,
  strict_private_edge_pop = strict_private_edge,
  strict_private_core_pop = strict_private_core,
  near_private_edge_pooled = near_private_edge_pooled,
  private_edge_pooled = private_edge_pooled,
  stringsAsFactors = FALSE
)
private_df <- dplyr::left_join(private_df, edge_pop_private_summary, by = "snp")
private_df <- dplyr::mutate(
  private_df,
  private_any_edge_pop = dplyr::coalesce(.data$private_any_edge_pop, FALSE),
  n_edge_pops_private = dplyr::coalesce(.data$n_edge_pops_private, 0L)
)
if (PRIVATE_SCOPE == "pooled_edge") {
  private_df$private_edge_pop <- private_df$private_edge_pooled
} else if (PRIVATE_SCOPE == "per_edge_pop") {
  private_df$private_edge_pop <- private_df$private_any_edge_pop
} else {
  private_df$private_edge_pop <- private_df$private_edge_pooled | private_df$private_any_edge_pop
}
private_df <- dplyr::left_join(private_df, snp_map, by = "snp")

write.table(
  edge_pop_private_long,
  file.path(results_dir, "private_per_edge_population.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  private_df,
  file.path(results_dir, "private_and_near_private_snps.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## ---- 3. Edge vs core Weir-Cockerham Fst (individuals pooled by group) ----

dos_sub <- str_dos[, snps, drop = FALSE]
grp <- popmap$edge_core[match(rownames(dos_sub), popmap$sample_id)]
fst_vec <- weir_cockerham_fst_loci(dos_sub, grp)

fst_df <- data.frame(
  snp = snps,
  fst_edge_core = fst_vec,
  stringsAsFactors = FALSE
)
fst_df <- dplyr::left_join(fst_df, private_df[, c(
  "snp", "mean_edge_freq", "mean_core_freq", "delta_p_edge_minus_core",
  "private_edge_pop", "private_any_edge_pop", "n_edge_pops_private", "edge_pops_private",
  "strict_private_edge_pop", "near_private_edge_pooled"
)], by = "snp")
fst_df <- dplyr::left_join(fst_df, snp_map, by = "snp")

fst_thr <- stats::quantile(fst_df$fst_edge_core, FST_QUANTILE, na.rm = TRUE)
fst_df$fst_outlier <- !is.na(fst_df$fst_edge_core) & fst_df$fst_edge_core >= fst_thr
fst_df$fst_quantile_threshold <- as.numeric(fst_thr)

write.table(
  fst_df,
  file.path(results_dir, "snp_fst_edge_core.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  data.frame(
    fst_quantile = FST_QUANTILE,
    fst_threshold = as.numeric(fst_thr),
    n_snps = nrow(fst_df),
    n_fst_outliers = sum(fst_df$fst_outlier, na.rm = TRUE),
    edge_freq_high = EDGE_FREQ_HIGH,
    core_freq_low = CORE_FREQ_LOW,
    other_pop_freq_max = OTHER_POP_FREQ_MAX,
    private_scope = PRIVATE_SCOPE,
    private_mode = PRIVATE_MODE,
    edge_core_method = EDGE_CORE_METHOD,
    method_note = edge_core_meta$note,
    stringsAsFactors = FALSE
  ),
  file.path(results_dir, "fst_threshold_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## ---- 4. Candidate SNP table (Fst outliers + private / enriched) ------------

candidates <- fst_df %>%
  dplyr::mutate(
    candidate_high_fst = .data$fst_outlier,
    candidate_private_edge = .data$private_edge_pop,
    candidate_strict_private_edge = .data$strict_private_edge_pop,
    candidate_near_private_edge = .data$near_private_edge_pooled,
    candidate_private_any_edge_pop = .data$private_any_edge_pop,
    candidate_delta_p = .data$delta_p_edge_minus_core >= DELTA_P_MIN,
    candidate_fixation_pattern = (.data$mean_edge_freq >= CAND_EDGE_FIX) &
      (.data$mean_core_freq <= CAND_CORE_RARE),
    candidate_any = .data$candidate_high_fst |
      .data$candidate_private_edge |
      .data$candidate_fixation_pattern
  )

write.table(
  candidates,
  file.path(results_dir, "edge_core_snp_candidates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  dplyr::filter(candidates, .data$candidate_any),
  file.path(results_dir, "edge_core_snp_candidates_filtered.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## ---- 5. Optional: GEA / RDA overlap ----------------------------------------

if (INTERSECT_GEA) {
  gea_dir <- lgp_outputs_step_dir("04-gea", root = data_dir)
  gea_paths <- list(
    lfmm = file.path(gea_dir, "results", "sign_lfmm_ridge_table.csv"),
    rda  = file.path(gea_dir, "results", "rda_p0.05.csv")
  )
  gea_snps <- character()
  if (file.exists(gea_paths$lfmm)) {
    lf <- utils::read.csv(gea_paths$lfmm, stringsAsFactors = FALSE)
    if ("snp" %in% names(lf)) gea_snps <- c(gea_snps, lf$snp)
  }
  if (file.exists(gea_paths$rda)) {
    rd <- utils::read.csv(gea_paths$rda, stringsAsFactors = FALSE)
    if ("snp" %in% names(rd)) gea_snps <- c(gea_snps, rd$snp)
  }
  if (length(gea_snps)) {
    gea_snps <- unique(gea_snps)
    gea_hits <- candidates %>%
      dplyr::filter(.data$snp %in% gea_snps, .data$candidate_any) %>%
      dplyr::mutate(in_gea = TRUE)
    write.table(
      gea_hits,
      file.path(results_dir, "edge_core_candidates_x_gea.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    message("[07-edge-core] GEA overlap: ", nrow(gea_hits), " edge-core candidates also in GEA lists.")
  } else {
    message("[07-edge-core] GEA files not found; skipping intersection.")
  }
}

## ---- 6. Plots: Fst distribution & genome-wide scan -------------------------

p_fst_dens <- ggplot(fst_df, aes(x = .data$fst_edge_core)) +
  geom_histogram(bins = 60, fill = "#2E5AAC", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = fst_thr, linetype = "dashed", colour = "#E87B14", linewidth = 0.9) +
  labs(
    title = "Edge vs core Weir-Cockerham Fst (per SNP)",
    subtitle = paste0(
      "Dashed line: ", 100 * FST_QUANTILE,
      "th percentile (Fst = ", signif(fst_thr, 4), "); ",
      edge_core_meta$note
    ),
    x = expression(F[ST]),
    y = "SNP count"
  ) +
  theme_minimal()
ggsave(
  file.path(plot_dir, "fst_edge_core_distribution.png"),
  p_fst_dens, width = 8, height = 5, dpi = 300
)

if (all(c("chr", "pos") %in% names(candidates)) && sum(!is.na(candidates$chr)) > 0) {
  manh <- candidates %>%
    dplyr::filter(!is.na(.data$chr), !is.na(.data$pos)) %>%
    dplyr::mutate(
      chr = as.character(.data$chr),
      pos = as.numeric(.data$pos),
      highlight = .data$candidate_any
    )
  chr_ord <- unique(manh$chr[order(as.numeric(gsub("[^0-9].*", "", manh$chr)), manh$chr)])
  manh$chr <- factor(manh$chr, levels = chr_ord)
  p_manh <- ggplot(manh, aes(x = .data$pos, y = .data$fst_edge_core, colour = .data$highlight)) +
    geom_point(alpha = 0.35, size = 0.5) +
    scale_colour_manual(
      values = c(`FALSE` = "grey60", `TRUE` = "#E87B14"),
      labels = c(`FALSE` = "Other", `TRUE` = "Edge-core candidate")
    ) +
    facet_wrap(~chr, scales = "free_x", ncol = 3) +
    geom_hline(yintercept = fst_thr, linetype = "dashed", colour = "#4E0707", linewidth = 0.4) +
    labs(
      title = "Genome-wide edge vs core Fst",
      x = "Position",
      y = expression(F[ST]),
      colour = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave(
    file.path(plot_dir, "fst_edge_core_manhattan.png"),
    p_manh, width = 12, height = 8, dpi = 300
  )
}

## ---- 7. Maps: per-population allele frequency for GEA x edge-core SNPs ----

if (MAP_GEA_SNPS) {
  plot_gea_candidate_freq_maps(
    results_dir, plot_dir, af_long, pop_centroids, pop_assignment_tbl,
    max_snps = MAP_GEA_MAX_SNPS,
    pawpaw_range_sf = pawpaw_range_plot
  )
}

message(
  "[07-edge-core] Done. Summary: ",
  sum(candidates$candidate_any, na.rm = TRUE), " candidates; ",
  "Fst threshold = ", signif(fst_thr, 4), " (q", FST_QUANTILE, "). Outputs: ", out_dir
)
