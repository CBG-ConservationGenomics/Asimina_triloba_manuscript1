# LD decay from PLINK output (run plink separately; see example below).

# Cached / consistent paths across the LandscapeGenomicsPipeline (see Scripts/lgp_pipeline_cache.R).
._lgpr <- Sys.getenv("LGP_PROJECT_ROOT", "~/Desktop/LandscapeGenomicsPipeline")
options(lgp.project_root = sub("/+$", "", path.expand(getOption("lgp.project_root", ._lgpr))))
suppressPackageStartupMessages(base::source(
  base::file.path(getOption("lgp.project_root"), "Scripts", "lgp_pipeline_cache.R"),
  encoding = "UTF-8"
))
rm(._lgpr)

data_dir <- lgp_project_root()
out_dir <- lgp_outputs_step_dir("05-linkage-decay")
plot_dir <- file.path(out_dir, "plots")
results_dir <- file.path(out_dir, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# GEA CSVs live under Outputs/…/04-gea/results (see Scripts/lgp_pipeline_cache.R for path resolution).
gea_results_dir <- file.path(lgp_outputs_base(data_dir), "04-gea", "results")
dir.create(gea_results_dir, recursive = TRUE, showWarnings = FALSE)
# Primary VCF location for this pipeline (legacy typo "Inputs/" kept as fallback only).
vcf_candidates <- c(
  file.path(data_dir, "InputData", "filtered_snps4landscape_reheader.vcf"),
  file.path(data_dir, "Inputs", "filtered_snps4landscape_reheader.vcf"),
  file.path(data_dir, "Pawpaw_snps", "remove", "filtered_snps4landscape_reheader.vcf")
)
vcf_for_ld <- vcf_candidates[file.exists(vcf_candidates)][1L]
if (is.na(vcf_for_ld)) {
  stop(
    "Cannot find filtered_snps4landscape_reheader.vcf. Tried:\n",
    paste("  ", vcf_candidates, collapse = "\n"),
    call. = FALSE
  )
}

# VCF ID column often '.'. Pipeline / GEA SNP names match CHROM + '_' + POS (algatr colnames).
# Set synthetic IDs during VCF load so top_snps.txt matches PLINK’s variant names.
# Template: PLINK replaces '@' → chrom code and '#' → position; literal '_' in between → e.g.
# Chr1 + 870969 → "Chr1_870969".
# Alternate: options(lgp.plink_set_missing_var_ids = "@:#") for "Chr1:870969".
# Rare duplicate positions/biallelic splits may need allele tokens in template (see plink docs).
pls_opt <- getOption("lgp.plink_set_missing_var_ids", "@_#")
plink_set_missing_ids <- trimws(as.character(pls_opt)[1L])
if (!nzchar(plink_set_missing_ids)) plink_set_missing_ids <- "@_#"

plink_vcf_intro <- c(
  "--vcf", vcf_for_ld,
  "--double-id", "--allow-extra-chr",
  "--set-missing-var-ids", plink_set_missing_ids
)

plink_ld_prefix <- file.path(out_dir, "pawpaw_LD")
plink_ld_ldgz <- paste0(plink_ld_prefix, ".ld.gz")

if (lgp_should_rerun_external(plink_ld_ldgz, vcf_for_ld)) {
  message("[lgp-cache] Running PLINK genome-wide pairwise LD (this can take a long time).")
  system2(
    "plink",
    args = c(
      plink_vcf_intro,
      "--maf", "0.05", "--geno", "0.1",
      "--r2", "gz",
      "--ld-window", "99999",
      "--ld-window-kb", "300",
      "--ld-window-r2", "0",
      "--out", plink_ld_prefix
    )
  )
} else {
  message("[lgp-cache] Skipping PLINK genome-wide LD — reused existing ", basename(plink_ld_ldgz))
}

library(data.table)
library(dplyr)
library(ggplot2)

# Mirrors Scripts/04-GEA.R — split CHROM_bp marker names → SNP, CHR, BP (bedtools LD windows).
gea_chr_bp_from_snp_markers <- function(snp_ids) {
  u <- unique(trimws(as.character(snp_ids)))
  u <- u[!is.na(u) & nzchar(u)]
  rx <- "^(.+)_([0-9]+)$"
  ok <- grepl(rx, u)
  if (any(!ok)) {
    warning(
      "gea_chr_bp_from_snp_markers: dropping ",
      sum(!ok),
      " marker(s) without trailing _<position> digits.",
      call. = FALSE
    )
  }
  u <- u[ok]
  if (!length(u)) {
    return(data.frame(SNP = character(0), CHR = character(0), BP = integer(0), stringsAsFactors = FALSE))
  }
  data.frame(
    SNP = u,
    CHR = sub(rx, "\\1", u),
    BP  = as.integer(sub(rx, "\\2", u)),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Concatenated unique SNP list for PLINK --ld-snp-list (significant LFMM ridge + pRDA Z outliers)
# Unique SNPs from sign_lfmm_ridge_table.csv (GEA LFMM ridge hits; tab-separated), merged with
# outliers from best_pRDA_sigZ_SNPS.csv (rda_snps or SNP).
# ---------------------------------------------------------------------------
top_snps_path <- file.path(out_dir, "top_snps.txt")
ld_top_prefix <- file.path(out_dir, "LD_top")
ld_top_ldgz <- paste0(ld_top_prefix, ".ld.gz")

pick_gea_csv <- function(name) {
  cands <- c(
    file.path(gea_results_dir, name),
    file.path(data_dir, "Outputs", "04-gea", "results", name),
    file.path(data_dir, "outputs", "04-gea", "results", name)
  )
  hit <- cands[file.exists(cands)][1L]
  if (is.na(hit)) NA_character_ else hit
}

lfmm_sign_csv <- pick_gea_csv("sign_lfmm_ridge_table.csv")
best_sigz_csv <- pick_gea_csv("best_pRDA_sigZ_SNPS.csv")

snps_from_sign_lfmm_ridge_table <- function(path) {
  if (is.na(path) || identical(path, "") || !file.exists(path)) return(character(0))
  tb <- data.table::fread(path, sep = "auto", data.table = FALSE, showProgress = FALSE)
  nm <- intersect(names(tb), c("snp", "SNP"))
  if (length(nm) != 1L) {
    warning("sign_lfmm_ridge_table.csv: expected column snp or SNP — ", path)
    return(character(0))
  }
  v <- unique(trimws(as.character(tb[[nm]])))
  v[!is.na(v) & nzchar(v)]
}

snps_from_best_sigz_csv <- function(path) {
  if (is.na(path) || identical(path, "") || !file.exists(path)) return(character(0))
  tb <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("rda_snps" %in% names(tb)) {
    v <- tb[["rda_snps"]]
  } else if ("SNP" %in% names(tb)) {
    v <- tb[["SNP"]]
  } else {
    warning(
      "best_pRDA_sigZ_SNPS.csv: expected columns rda_snps (algatr tidy) or SNP — ",
      path
    )
    return(character(0))
  }
  v <- unique(trimws(as.character(v)))
  v[!is.na(v) & nzchar(v)]
}

lfmm_snps <- snps_from_sign_lfmm_ridge_table(lfmm_sign_csv)
sigz_snps <- snps_from_best_sigz_csv(best_sigz_csv)

pRDA_gea_positions <- gea_chr_bp_from_snp_markers(sigz_snps)
lfmm_gea_positions <- gea_chr_bp_from_snp_markers(lfmm_snps)
utils::write.table(
  pRDA_gea_positions,
  file      = file.path(gea_results_dir, "pRDA_snp_positions.tsv"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
utils::write.csv(lfmm_gea_positions, file.path(gea_results_dir, "lfmm_snp_positions.csv"), row.names = FALSE)
message(
  "Wrote SNP position tables for LD/bedtools: pRDA (n=", nrow(pRDA_gea_positions),
  "), LFMM sign-hit (n=", nrow(lfmm_gea_positions), ") → ",
  basename(gea_results_dir)
)

top_snps_union <- sort(unique(c(lfmm_snps, sigz_snps)))

if (length(top_snps_union) > 0L) {
  writeLines(top_snps_union, top_snps_path)
  data.table::fwrite(
    data.frame(
      snp           = top_snps_union,
      source        = ifelse(top_snps_union %in% sigz_snps,
        ifelse(top_snps_union %in% lfmm_snps, "both", "best_pRDA_sigZ_SNPS.csv"),
        "sign_lfmm_ridge_table.csv"
      ),
      stringsAsFactors = FALSE
    ),
    file.path(results_dir, "top_snps_sources.csv")
  )
  message(
    "Wrote top_snps.txt: ", length(top_snps_union), " unique SNPs (sign LFMM ridge: ",
    length(lfmm_snps), " unique; best_pRDA_sigZ: ",
    length(sigz_snps), " unique)."
  )
} else {
  message(
    "top_snps.txt not written — could not load SNPs. Expected one of:\n",
    "  ", file.path(gea_results_dir, "sign_lfmm_ridge_table.csv"), "\n",
    "  ", file.path(gea_results_dir, "best_pRDA_sigZ_SNPS.csv"),
    "\n(Run 04-GEA.R first, or adjust paths.)"
  )
}

ld_ldgz <- paste0(plink_ld_prefix, ".ld.gz")
ld <- fread(ld_ldgz)
ld <- ld %>% mutate(dist = abs(BP_B - BP_A))

# Pairwise LD matrix from PLINK is already stored at ld_ldgz (gzip). Below we export
# binned summaries and threshold summaries under results_dir.

decay <- ld %>%
  group_by(bin = cut(dist, breaks = seq(0, 300000, by = 1000))) %>%
  summarise(mean_r2 = mean(R2, na.rm = TRUE),
            mid_dist = mean(dist, na.rm = TRUE)) %>%
  ungroup()

decay_tbl <- decay %>%
  dplyr::mutate(mid_dist_kb = .data$mid_dist / 1000)

data.table::fwrite(decay_tbl, file.path(results_dir, "ld_decay_mean_r2_by_kb_bin.csv"))

p_decay <- ggplot(decay_tbl, aes(x = .data$mid_dist_kb, y = .data$mean_r2)) +
  geom_line() +
  labs(
    title = "Linkage disequilibrium decay (mean pairwise r²)",
    x = "Distance (kb)",
    y = expression(italic(r)^2)
  ) +
  theme_classic()

ggsave(
  filename = file.path(plot_dir, "ld_decay_mean_r2.png"),
  plot = p_decay,
  width = 8,
  height = 5,
  dpi = 200,
  bg = "white"
)

target_r2 <- 0.2
decay_below <- decay %>%
  dplyr::arrange(.data$mid_dist) %>%
  dplyr::filter(.data$mean_r2 <= target_r2)

ld_thresh <- if (nrow(decay_below) > 0L) {
  decay_below %>% dplyr::slice(1L) %>% dplyr::pull(.data$mid_dist)
} else {
  NA_real_
}
ld_thresh_kb <- ld_thresh / 1000

thresh_tbl <- data.frame(
  target_r2                   = target_r2,
  distance_bp_first_mean_le_q = ld_thresh,
  distance_kb_first_mean_le_q = ld_thresh_kb,
  plink_ld_prefix             = plink_ld_prefix,
  pairwise_ld_file            = basename(ld_ldgz),
  stringsAsFactors            = FALSE
)
data.table::fwrite(thresh_tbl, file.path(results_dir, "ld_decay_distance_at_target_r2.csv"))

if (is.na(ld_thresh)) {
  warning(
    "No distance bin had mean r² ≤ ", target_r2,
    "; ld_decay_distance_at_target_r2.csv has NA distances."
  )
}

print(p_decay)

# Optional: inspect console
ld_thresh
ld_thresh_kb

# PLINK focal SNP pairwise LD (--ld-window-kb must align with genomic half-span below).
focal_ld_window_kb <- suppressWarnings(as.integer(getOption(
  "lgp.focal_ld_window_kb",
  200L
)[1]))
if (!is.finite(focal_ld_window_kb) || focal_ld_window_kb < 1L) focal_ld_window_kb <- 200L
focal_ld_halfwidth_bp <- as.integer(focal_ld_window_kb) * 1000L

ld_bed_min_halfwidth_bp <- suppressWarnings(as.integer(getOption(
  "lgp.ld_bed_min_halfwidth_bp",
  5000L
)[1]))
if (!is.finite(ld_bed_min_halfwidth_bp) || ld_bed_min_halfwidth_bp < 100L) {
  ld_bed_min_halfwidth_bp <- 5000L
}

ld_block_r2 <- suppressWarnings(as.numeric(getOption(
  "lgp.ld_block_r2_threshold",
  0.2
)[1]))
if (!is.finite(ld_block_r2) || ld_block_r2 <= 0 || ld_block_r2 > 1) ld_block_r2 <- 0.2

####################################
# Focal-SNPs LD: SNP list is generated above from LFMM ridge + best pRDA Z hits
# (sign_lfmm_ridge_table.csv ∪ best_pRDA_sigZ_SNPS.csv → top_snps.txt).
###################################

if (file.exists(top_snps_path)) {
  tl <- trimws(readLines(top_snps_path, warn = FALSE))
  tl <- unique(tl[nzchar(tl) & !is.na(tl)])
  n_top <- length(tl)

  plink_stderr <- tempfile(pattern = "plink_focal_stderr_", tmpdir = out_dir, fileext = ".log")

  # Pairwise --r2 output needs ≥2 variant IDs understood by PLINK.
  if (n_top < 2L) {
    warning(
      "Skipping focal PLINK LD: top_snps.txt has ", n_top,
      " unique id(s); need ≥2 to produce an LD pairs file.",
      "\n(Add more outliers in GEA, or widen significance so the merged list crosses two SNPs.)"
    )
  } else {
    focal_ld_plain <- paste0(ld_top_prefix, ".ld")
    focal_ld_product <- if (file.exists(ld_top_ldgz)) ld_top_ldgz else focal_ld_plain

    run_focal_plink <- lgp_should_rerun_external(focal_ld_product, c(vcf_for_ld, top_snps_path))

    if (!run_focal_plink) {
      message("[lgp-cache] Skipping PLINK focal-SNP LD — reused ", basename(focal_ld_product))
    } else {
    rc_ld <- suppressWarnings(system2(
      "plink",
      args = c(
        plink_vcf_intro,
        "--ld-snp-list", top_snps_path,
        "--r2", "gz",
        "--ld-window-kb", as.character(focal_ld_window_kb),
        "--ld-window", "99999",
        "--ld-window-r2", "0",
        "--out", ld_top_prefix
      ),
      stderr = plink_stderr
    ))

    if (file.exists(plink_stderr)) {
      slog <- paste(readLines(plink_stderr, warn = FALSE), collapse = "\n")
      if (nzchar(trimws(slog))) {
        message("PLINK focal LD log (from stderr):\n", slog)
      }
      unlink(plink_stderr)
    }

    if (identical(rc_ld, NA_integer_)) {
      warning("`plink` not found on PATH — install PLINK or add it to PATH and re-run.")
    } else if (rc_ld != 0L) {
      warning(
        "PLINK focal LD exited with status ", rc_ld,
        ". If IDs do not resolve, check VCF #CHROM names vs top_snps (PLINK ",
        "`--set-missing-var-ids` ", shQuote(plink_set_missing_ids, type = "sh"), "); see LD_top.log."
      )
    }

    ld_top_plain <- paste0(ld_top_prefix, ".ld")
    if (!file.exists(ld_top_ldgz) && !file.exists(ld_top_plain)) {
      sibs <- list.files(dirname(ld_top_prefix), pattern = "^LD_top[.]")
      hint <- if (length(sibs)) paste0(": ", paste(sort(sibs), collapse = ", "), ".") else "."
      warning(
        "No LD_top.ld.gz or LD_top.ld after focal PLINK", hint,
        "\nUsually this means focal SNPs in top_snps.txt are absent from --vcf, or incompatible IDs."
      )
    }
    }
  }
} else {
  message("Skipping focal-SNP LD — no ", basename(top_snps_path), " (see GEA step above).")
}

########
# LD windows around focal SNPs (after PLINK --ld-snp-list, etc.)
########

# Outputs expected from GEA step (adjust filenames here if you change 04-GEA.R).
pRDA_snp_pos_path <- file.path(gea_results_dir, "pRDA_snp_positions.tsv")
lfmm_snp_pos_path <- file.path(gea_results_dir, "lfmm_snp_positions.csv")

ld_top_file <- if (file.exists(ld_top_ldgz)) {
  ld_top_ldgz
} else if (file.exists(paste0(ld_top_prefix, ".ld"))) {
  paste0(ld_top_prefix, ".ld")
} else {
  NA_character_
}

plink_ld_pick_name <- function(nms, want) {
  k <- match(tolower(want), tolower(nms))[1]
  if (is.na(k)) NA_character_ else nms[k]
}

plink_ld_pick_r_col <- function(nms) {
  intersect(c("R2", "R"), nms)[1]
}

# Symmetric aggregation: focal SNPs may appear only as SNP_B in PLINK output.
gea_aggregate_ld_pairs_symmetric <- function(ld_df, r2_min) {
  nms <- names(ld_df)
  rcol <- plink_ld_pick_r_col(nms)
  snpa <- plink_ld_pick_name(nms, "SNP_A")
  snpb <- plink_ld_pick_name(nms, "SNP_B")
  bpa <- plink_ld_pick_name(nms, "BP_A")
  bpb <- plink_ld_pick_name(nms, "BP_B")
  chra <- plink_ld_pick_name(nms, "CHR_A")
  chrb <- plink_ld_pick_name(nms, "CHR_B")

  if (is.na(rcol)) {
    warning("LD file has no R2 or R column — cannot summarise blocks.")
    return(dplyr::tibble(
      focal_id = character(0),
      BP = integer(0),
      chr_tag_ld = character(0),
      max_highld_dist_bp = integer(0),
      n_pairs_ge_r2 = integer(0)))
  }
  if (any(is.na(c(snpa, snpb, bpa, bpb)))) {
    warning("LD columns need SNP_A, SNP_B, BP_A, BP_B — found:\n  ",
      paste(nms, collapse = ", "))
    return(dplyr::tibble(
      focal_id = character(0),
      BP = integer(0),
      chr_tag_ld = character(0),
      max_highld_dist_bp = integer(0),
      n_pairs_ge_r2 = integer(0)))
  }

  use_chr_ld <- !(is.na(chra) || is.na(chrb))
  if (use_chr_ld) {
    base <- dplyr::transmute(
      ld_df,
      focal_a = .data[[snpa]], focal_b = .data[[snpb]],
      BPA = as.integer(.data[[bpa]]), BPb = as.integer(.data[[bpb]]),
      chr_at_a = as.character(.data[[chra]]),
      chr_bt_b = as.character(.data[[chrb]]),
      RR = suppressWarnings(as.numeric(.data[[rcol]])),
      dist_bp = abs(BPb - BPA))
  } else {
    base <- dplyr::transmute(
      ld_df,
      focal_a = .data[[snpa]], focal_b = .data[[snpb]],
      BPA = as.integer(.data[[bpa]]), BPb = as.integer(.data[[bpb]]),
      RR = suppressWarnings(as.numeric(.data[[rcol]])),
      dist_bp = abs(BPb - BPA))
  }

  hi <- dplyr::filter(
    base,
    is.finite(.data$RR) & .data$RR >= !!r2_min &
      nzchar(trimws(as.character(.data$focal_a))) &
      nzchar(trimws(as.character(.data$focal_b))) &
      .data$focal_a != .data$focal_b
  )

  if (use_chr_ld) {
    pa <- dplyr::transmute(hi, focal_id = .data$focal_a, chr_tag = .data$chr_at_a, BP = .data$BPA, dist_bp = .data$dist_bp)
    pb <- dplyr::transmute(hi, focal_id = .data$focal_b, chr_tag = .data$chr_bt_b, BP = .data$BPb, dist_bp = .data$dist_bp)
  } else {
    pa <- dplyr::transmute(hi, focal_id = .data$focal_a, chr_tag = NA_character_, BP = .data$BPA, dist_bp = .data$dist_bp)
    pb <- dplyr::transmute(hi, focal_id = .data$focal_b, chr_tag = NA_character_, BP = .data$BPb, dist_bp = .data$dist_bp)
  }

  long_part <- dplyr::bind_rows(pa, pb)

  if (nrow(long_part) == 0L) {
    return(dplyr::tibble(
      focal_id          = character(0),
      BP                = integer(0),
      chr_tag_ld        = character(0),
      max_highld_dist_bp = integer(0),
      n_pairs_ge_r2     = integer(0)))
  }

  long_part %>%
    dplyr::filter(is.finite(.data$BP)) %>%
    dplyr::group_by(.data$focal_id, .data$BP) %>%
    dplyr::summarise(
      chr_tag_ld = {
        ux <- stats::na.omit(unique(as.character(.data$chr_tag)))
        if (length(ux)) as.character(ux[[1]]) else NA_character_
      },
      max_highld_dist_bp = suppressWarnings(as.integer(max(.data$dist_bp, na.rm = TRUE))),
      n_pairs_ge_r2 = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      max_highld_dist_bp = dplyr::if_else(
        is.na(.data$max_highld_dist_bp) | .data$max_highld_dist_bp < 0L,
        0L,
        .data$max_highld_dist_bp
      )
    )
}

# One BED row per focal in snp_tbl: symmetric LD extent + sensible floor/fallback widths.
make_ld_windows_bed <- function(
    snp_pos_tbl,
    ld_summary_tbl,
    outfile,
    focal_half_bp,
    floor_half_bp
  ) {
  st <- dplyr::rename_with(as.data.frame(snp_pos_tbl), ~ toupper(trimws(.)))
  ncol_s <- colnames(st)
  sci <- intersect(ncol_s, c("SNP", "CHR", "BP"))
  if (!all(c("SNP", "CHR", "BP") %in% sci)) {
    stop("Position table needs SNP, CHR, BP (any case); got: ", paste(ncol_s, collapse = ", "))
  }

  bd <- dplyr::transmute(
    st,
    SNP = trimws(as.character(.data[["SNP"]])),
    CHR = trimws(as.character(.data[["CHR"]])),
    BP = suppressWarnings(as.integer(.data[["BP"]]))
  ) %>%
    dplyr::filter(nzchar(.data$SNP) & nzchar(.data$CHR) & is.finite(.data$BP))

  ld <- dplyr::rename_with(as.data.frame(ld_summary_tbl), ~ tolower(trimws(.)))
  req_ld <- c("focal_id", "bp", "max_highld_dist_bp")
  if (!all(req_ld %in% names(ld))) {
    warning("LD summary table missing focal_id/BP/max_highld_dist_bp; using empty join for ", outfile)
    ld_join <- dplyr::tibble(
      focal_id = character(0),
      bp_ld = integer(0),
      max_highld_dist_bp = integer(0)
    )
  } else {
    ld_join <- dplyr::transmute(
      ld,
      focal_id = trimws(as.character(.data[["focal_id"]])),
      bp_ld = suppressWarnings(as.integer(.data[["bp"]])),
      max_highld_dist_bp = suppressWarnings(as.integer(.data[["max_highld_dist_bp"]]))
    )
  }

  fh <- suppressWarnings(as.integer(focal_half_bp)[1])
  fl <- suppressWarnings(as.integer(floor_half_bp)[1])
  if (!is.finite(fh) || fh < 1000L) fh <- 200000L
  if (!is.finite(fl) || fl < 10L) fl <- 5000L

  w <- dplyr::left_join(bd, ld_join, by = c("SNP" = "focal_id", "BP" = "bp_ld")) %>%
    dplyr::mutate(
      max_hi = dplyr::coalesce(as.integer(.data$max_highld_dist_bp), NA_integer_),
      radius_bp = dplyr::case_when(
        !is.na(.data$max_hi) & .data$max_hi > 0L ~ as.integer(pmax(fl, pmin(as.numeric(.data$max_hi), as.numeric(fh)))),
        TRUE ~ fh
      ),
      start = suppressWarnings(as.integer(pmax(0L, .data$BP - .data$radius_bp))),
      end = suppressWarnings(as.integer(.data$BP + .data$radius_bp))
    ) %>%
    dplyr::filter(is.finite(.data$start), is.finite(.data$end), .data$end > .data$start)

  w <- dplyr::distinct(dplyr::select(w, "CHR", "start", "end", "SNP"))

  data.table::fwrite(
    w %>% dplyr::select(.data$CHR, .data$start, .data$end, .data$SNP),
    outfile,
    sep = "\t",
    col.names = FALSE
  )
}

if (!is.na(ld_top_file)) {
  ld_top <- fread(ld_top_file)
  if (!all(c("BP_A", "BP_B") %in% names(ld_top))) {
    warning("LD file missing BP_A/BP_B — pairwise distance not computed.\n ", paste(names(ld_top), collapse = ", "))
    ld_top <- dplyr::mutate(ld_top, dist_bp = NA_real_)
  } else {
    ld_top <- dplyr::mutate(
      ld_top,
      dist_bp = abs(as.integer(.data[["BP_A"]]) - as.integer(.data[["BP_B"]]))
    )
  }

  data.table::fwrite(ld_top, file.path(results_dir, "LD_top_pairwise.tsv.gz"))

  r_plot <- intersect(c("R2", "R"), names(ld_top))[1]
    if (!is.na(r_plot)) {
    lt_plot <- if (nrow(ld_top) > 500000L) dplyr::slice_sample(ld_top, n = 500000L) else ld_top
    lt_plot <- dplyr::filter(lt_plot, is.finite(.data$dist_bp))
    if (nrow(lt_plot) == 0L) {
      message("LD plot skipped: no finite pairwise distances (check BP columns).")
    } else {
    lt_plot <- dplyr::transmute(lt_plot, dist_kb = .data$dist_bp / 1000, rval = suppressWarnings(as.numeric(.data[[r_plot]])))

    p_top <- ggplot(lt_plot, aes(x = .data$dist_kb, y = .data$rval)) +
      geom_point(alpha = 0.2, size = 0.35, stroke = 0) +
      labs(
        title = "Pairwise LD among focal SNPs (--ld-snp-list)",
        subtitle = sprintf("Correlation column \"%s\"; r²/threshold summaries use gea_aggregate_ld_pairs_symmetric()", r_plot),
        x = "Distance (kb)",
        y = if (identical(r_plot, "R2")) expression(italic(r)^2) else expression(italic(r))
      ) +
      theme_classic()

    ggsave(
      filename = file.path(plot_dir, "LD_top_pairwise.png"),
      plot = p_top,
      width = 9,
      height = 5,
      dpi = 200,
      bg = "white"
    )
    }
  }

  ld_blocks <- gea_aggregate_ld_pairs_symmetric(ld_top, ld_block_r2)
  data.table::fwrite(ld_blocks, file.path(results_dir, "LD_top_high_ld_blocks_summary.csv"))

  if (file.exists(pRDA_snp_pos_path) && file.exists(lfmm_snp_pos_path)) {
    pRDA_snp_pos <- fread(pRDA_snp_pos_path)
    lfmm_snp_pos <- fread(lfmm_snp_pos_path)

    make_ld_windows_bed(pRDA_snp_pos, ld_blocks, file.path(out_dir, "pRDA_top_snps_LD_windows.bed"),
      focal_half_bp = focal_ld_halfwidth_bp, floor_half_bp = ld_bed_min_halfwidth_bp)
    make_ld_windows_bed(lfmm_snp_pos, ld_blocks, file.path(out_dir, "LFMM_top_snps_LD_windows.bed"),
      focal_half_bp = focal_ld_halfwidth_bp, floor_half_bp = ld_bed_min_halfwidth_bp)
    message(
      "LD block BED paths (one row/SNP): radius = max(high-LD pairwise distance vs floor ",
      ld_bed_min_halfwidth_bp, " bp), capped ", focal_ld_halfwidth_bp,
      " bp; default window if no mates above r² is ", focal_ld_halfwidth_bp, " bp each side."
    )
  } else {
    message(
      "Skipping LD window BED export (needs SNP position tables):\n",
      "  ", pRDA_snp_pos_path, "\n",
      "  ", lfmm_snp_pos_path
    )
  }
} else {
  message(
    "No focal LD table yet. Expected ",
    basename(ld_top_ldgz),
    " (or ",
    basename(paste0(ld_top_prefix, ".ld")),
    "). After ≥2 focal SNPs in top_snps.txt, check PLINK messages above ",
    "(IDs must occur in ",
    basename(vcf_for_ld), ")."
  )
}

# ---- bedtools: merge LD-window BEDs and intersect annotated genes (GFF) ------------
# Mirrors the former shell one-liners; requires `bedtools` on PATH.
# Override GFF with env LGP_BEDTOOLS_GFF or options(lgp.bedtools_gff = "...").
bed_ld_pRDA <- file.path(out_dir, "pRDA_top_snps_LD_windows.bed")
bed_ld_LFMM <- file.path(out_dir, "LFMM_top_snps_LD_windows.bed")
bed_ld_merged <- file.path(results_dir, "merged_top_snps_LD_windows.bed")
tsv_linked_merged <- file.path(results_dir, "linked_genes.tsv")
tsv_linked_pRDA <- file.path(results_dir, "pRDA_linked_genes.tsv")
tsv_linked_lfmm <- file.path(results_dir, "LFMM_linked_genes.tsv")

gene_go_terms_wide <- file.path(data_dir, "InputData", "gene_go_terms_wide.tsv")
gene_func_wide <- file.path(data_dir, "InputData", "gene_functional_descriptions_wide.tsv")

lgp_shell_run <- function(cmd) {
  system2("/bin/bash", args = c("-lc", paste0("set -euo pipefail; ", cmd)))
}

lgp_interval_bed_nonempty <- function(path) {
  if (!isTRUE(file.exists(path))) {
    return(FALSE)
  }
  lines <- suppressWarnings(readLines(path, warn = FALSE, n = 20000L))
  lines <- trimws(lines)
  lines <- lines[!grepl("^\\s*$", lines)]
  any(grepl("^[^#]", lines))
}

#' Parse gene ID from GFF column 9 style attributes (matches Scripts/06-GO_enrichment.R).
lgp_gene_id_from_gff_attrs_col <- function(attrs_vec) {
  vapply(attrs_vec, function(a) {
    a <- trimws(as.character(a)[1])
    if (!nzchar(a) || is.na(a)) {
      return(NA_character_)
    }
    m <- regexec("gene_id=([^;]+)", a)
    rm <- regmatches(a, m)[[1]]
    if (length(rm) >= 2L) {
      return(trimws(rm[[2]]))
    }
    m2 <- regexec("(?:^|;)ID=([^;]+)", a)
    rm2 <- regmatches(a, m2)[[1]]
    if (length(rm2) >= 2L) {
      return(trimws(rm2[[2]]))
    }
    NA_character_
  }, character(1))
}

#' Append gene_id, gene_go_terms, gene_function_description to bedtools intersect -wao output.
lgp_enrich_linked_genes_tsv <- function(tsv_path, go_dt, desc_dt) {
  if (!isTRUE(file.exists(tsv_path))) {
    return(FALSE)
  }
  line1 <- suppressWarnings(readLines(tsv_path, n = 1L, warn = FALSE))
  hdr_line <- trimws(as.character(line1))
  skip_hdr <- nzchar(hdr_line) && grepl("^bed_chr\t", hdr_line)
  raw <- data.table::fread(
    tsv_path,
    sep = "\t",
    quote = "",
    header = FALSE,
    skip = if (skip_hdr) 1L else 0L,
    stringsAsFactors = FALSE
  )
  if (!ncol(raw) || nrow(raw) == 0L) {
    return(FALSE)
  }
  raw_nc <- ncol(raw)
  attrs_col_idx <- raw_nc - 1L
  if (attrs_col_idx < 1L) {
    warning("[bedtools] Cannot annotate linked genes (unexpected columns): ", basename(tsv_path),
      call. = FALSE)
    return(FALSE)
  }
  std_names <- c(
    "bed_chr", "bed_start", "bed_end", "focal_snp",
    "gff_seqname", "gff_source", "gff_feature", "gff_start", "gff_end",
    "gff_score", "gff_strand", "gff_frame", "gff_attributes",
    "overlap_bp"
  )
  if (raw_nc == length(std_names)) {
    data.table::setnames(raw, std_names)
  } else {
    warning(
      "[bedtools] Expected ",
      length(std_names),
      " bedtools -wao columns; got ",
      raw_nc,
      " (still annotating by attributes column index).",
      call. = FALSE
    )
    data.table::setnames(raw, paste0("col", seq_len(raw_nc)))
    data.table::setnames(raw, attrs_col_idx, "gff_attributes")
  }

  gene_ids <- lgp_gene_id_from_gff_attrs_col(raw[["gff_attributes"]])
  raw[, gene_id := gene_ids]
  raw[, gene_go_terms := go_dt$gene_go_terms[data.table::chmatch(raw$gene_id, go_dt$gene_id)]]
  raw[, gene_function_description :=
        desc_dt$gene_function_description[data.table::chmatch(raw$gene_id, desc_dt$gene_id)]]

  tmp_out <- tempfile(
    pattern = paste0(gsub("[^A-Za-z0-9]", "_", basename(tsv_path)), "_"),
    tmpdir = dirname(tsv_path),
    fileext = ".tsv"
  )
  data.table::fwrite(raw, tmp_out, sep = "\t", quote = "auto", na = "")
  if (isTRUE(file.exists(tsv_path))) {
    unlink(tsv_path)
  }
  if (!file.rename(tmp_out, tsv_path)) {
    file.copy(tmp_out, tsv_path, overwrite = TRUE)
    unlink(tmp_out)
  }
  TRUE
}

gff_ann <- trimws(as.character(Sys.getenv("LGP_BEDTOOLS_GFF", "")[1L]))
if (!nzchar(gff_ann)) {
  ox <- trimws(as.character(getOption("lgp.bedtools_gff", ""))[1L])
  if (length(ox) && nzchar(ox)) {
    gff_ann <- ox
  }
}
if (!nzchar(gff_ann)) {
  gcand <- file.path(data_dir, c("InputData", "Inputs"),
    "a_triloba_filtered_longest_isoform_gene.gff")
  gff_ann <- gcand[file.exists(gcand)][1L]
  if (is.na(gff_ann)) gff_ann <- ""
}

bedtools_bin <- Sys.which("bedtools")

if (!nzchar(bedtools_bin)) {
  message(
    "[bedtools] Executable not found on PATH; skipping LD-window \u2229 gene step.\n",
    "Install bedtools 2.x, then rerun this script."
  )
} else if (!nzchar(gff_ann) || !file.exists(gff_ann)) {
  warning(
    "[bedtools] GFF not resolved (try LGP_BEDTOOLS_GFF or options(lgp.bedtools_gff)). ",
    "Skipping intersect; expected under InputData/ for this pipeline.",
    call. = FALSE
  )
} else {
  qb <- function(x) shQuote(x)

  can_enrich <- file.exists(gene_go_terms_wide) && file.exists(gene_func_wide)
  if (!can_enrich) {
    warning(
      "[bedtools] Missing gene_go_terms_wide.tsv or gene_functional_descriptions_wide.tsv under InputData;\n",
      "  pRDA/LFMM linked-gene tables will omit GO/description columns.",
      call. = FALSE
    )
    enrich_deps <- character()
    go_dt <- NULL
    desc_dt <- NULL
  } else {
    enrich_deps <- c(gene_go_terms_wide, gene_func_wide)
    go_dt <- data.table::fread(gene_go_terms_wide, header = FALSE)
    data.table::setnames(go_dt, c("gene_id", "gene_go_terms"))
    go_dt <- unique(go_dt, by = "gene_id")
    desc_dt <- data.table::fread(gene_func_wide, header = FALSE)
    data.table::setnames(desc_dt, c("gene_id", "gene_function_description"))
    desc_dt <- unique(desc_dt, by = "gene_id")
  }

  merged_inputs <- stats::setNames(
    c(bed_ld_pRDA, bed_ld_LFMM),
    c("pRDA", "LFMM")
  )
  merged_inputs <- merged_inputs[vapply(merged_inputs, lgp_interval_bed_nonempty, logical(1))]
  deps_gff <- c(gff_ann)

  if (length(merged_inputs) >= 1L &&
        lgp_should_rerun_external(bed_ld_merged, c(unname(merged_inputs), deps_gff))) {
    merged_script <- sprintf(
      "cat %s | %s sort | %s merge > %s",
      paste(vapply(unname(merged_inputs), qb, character(1)), collapse = " "),
      qb(bedtools_bin),
      qb(bedtools_bin),
      qb(bed_ld_merged)
    )
    st_merge <- lgp_shell_run(merged_script)
    if (!identical(st_merge, 0L)) {
      warning("[bedtools] merge/sort exited with status ", st_merge, call. = FALSE)
    } else {
      message("[bedtools] wrote ", basename(bed_ld_merged))
    }
  } else if (length(merged_inputs) == 0L) {
    message("[bedtools] No nonempty LD-window BED files yet; skipping merge/intersect.")
  }

  intersect_wao <- function(bed_a, tsv_out, enrich_this_file = FALSE) {
    if (!lgp_interval_bed_nonempty(bed_a)) {
      return(invisible(FALSE))
    }
    deps_core <- c(bed_a, deps_gff)
    deps_all <- if (isTRUE(enrich_this_file) && can_enrich) {
      c(deps_core, enrich_deps)
    } else {
      deps_core
    }

    if (!lgp_should_rerun_external(tsv_out, deps_all)) {
      message("[bedtools] skip (fresh): ", basename(tsv_out))
      return(invisible(TRUE))
    }

    need_bt <- lgp_should_rerun_external(tsv_out, deps_core)

    if (need_bt) {
      int_scr <- sprintf(
        "%s intersect -wao -a %s -b %s > %s",
        qb(bedtools_bin), qb(bed_a), qb(gff_ann), qb(tsv_out)
      )
      st <- lgp_shell_run(int_scr)
      if (!identical(st, 0L)) {
        warning("[bedtools] intersect exited with status ", st, " for ", basename(tsv_out), call. = FALSE)
        return(invisible(FALSE))
      }
      message("[bedtools] intersect -> ", basename(tsv_out))
    }

    if (isTRUE(enrich_this_file) && can_enrich && file.exists(tsv_out)) {
      ok_en <- lgp_enrich_linked_genes_tsv(tsv_out, go_dt, desc_dt)
      if (isTRUE(ok_en)) {
        message("[bedtools] added gene_id / GO terms / descriptions -> ", basename(tsv_out))
      }
    }

    invisible(TRUE)
  }

  intersect_wao(bed_ld_merged, tsv_linked_merged, enrich_this_file = FALSE)
  if (length(merged_inputs) >= 1L && !lgp_interval_bed_nonempty(bed_ld_merged)) {
    message(
      "[bedtools] ",
      basename(bed_ld_merged),
      " is missing or empty after merge step (merge may have failed or been skipped)."
    )
  }

  intersect_wao(bed_ld_pRDA, tsv_linked_pRDA, enrich_this_file = TRUE)
  intersect_wao(bed_ld_LFMM, tsv_linked_lfmm, enrich_this_file = TRUE)
}
