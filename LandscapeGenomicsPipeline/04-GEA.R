#######################################################################
# Genotype-Environment Association (GEA) with RDA and LFMM
# https://thewanglab.github.io/algatr/articles/RDA_vignette.html
#######################################################################

rda_packages()
library(algatr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(stringr)

set.seed(1234)

# Algatr colnames are CHROM_<bp>; split for bedtools / LD windows (CHR + BP).
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
# SNP outlier color for every ggplot in this script (matches LFMM / Manhattan helpers)
gea_outlier_col <- "#4E0707"
gea_nonoutlier_col <- grDevices::rgb(0.7, 0.7, 0.7, 0.5)

gea_scale_color_snp_type <- function() {
  ggplot2::scale_color_manual(
    values = c(`Non-outlier` = gea_nonoutlier_col, Outlier = gea_outlier_col),
    na.translate = FALSE
  )
}
gea_scale_fill_snp_type <- function() {
  ggplot2::scale_fill_manual(
    values = c(`Non-outlier` = gea_nonoutlier_col, Outlier = gea_outlier_col),
    na.translate = FALSE
  )
}

# Split SNP layers so grey non-outliers are drawn first; dark-red outliers on top (visible).
gea_geom_snp_biplot <- function(TAB_snps_sub) {
  sn_base <- dplyr::filter(TAB_snps_sub, as.character(.data$type) == "Non-outlier")
  sn_out  <- dplyr::filter(TAB_snps_sub, as.character(.data$type) == "Outlier")
  list(
    ggplot2::geom_point(
      data = sn_base,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$type),
      size     = 1.25,
      alpha    = 0.45,
      shape    = 16
    ),
    ggplot2::geom_point(
      data = sn_out,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$type),
      size     = 2.05,
      alpha    = 1,
      shape    = 16
    )
  )
}

# algatr::rda_biplot layout, with outliers plotted last so they are not hidden.
rda_biplot_gea <- function(TAB_snps, TAB_var, biplot_axes = c(1, 2)) {
  xax <- paste0("RDA", biplot_axes[1])
  yax <- paste0("RDA", biplot_axes[2])
  TAB_snps_sub <- TAB_snps[, c(xax, yax, "type")]
  colnames(TAB_snps_sub) <- c("x", "y", "type")
  TAB_var_sub <- TAB_var[, c(xax, yax), drop = FALSE]
  colnames(TAB_var_sub) <- c("x", "y")
  TAB_var_sub$x <- TAB_var_sub$x * max(TAB_snps_sub$x, na.rm = TRUE) /
    stats::quantile(TAB_var_sub$x, na.rm = TRUE)[4]
  TAB_var_sub$y <- TAB_var_sub$y * max(TAB_snps_sub$y, na.rm = TRUE) /
    stats::quantile(TAB_var_sub$y, na.rm = TRUE)[4]
  layers <- gea_geom_snp_biplot(TAB_snps_sub)
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.8), linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.8), linewidth = 0.6) +
    layers[[1]] +
    layers[[2]] +
    gea_scale_color_snp_type() +
    ggplot2::geom_segment(
      data = TAB_var_sub,
      ggplot2::aes(xend = .data$x, yend = .data$y, x = 0, y = 0),
      colour    = "black",
      linewidth = 0.15,
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))
    ) +
    ggrepel::geom_text_repel(
      data = TAB_var_sub,
      ggplot2::aes(x = .data$x, y = .data$y, label = rownames(TAB_var_sub)),
      size         = 4,
      max.overlaps = Inf,
      segment.size = 0.2
    ) +
    ggplot2::xlab(xax) +
    ggplot2::ylab(yax) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "SNP type")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.background    = ggplot2::element_blank(),
      legend.background   = ggplot2::element_blank(),
      panel.grid          = ggplot2::element_blank(),
      plot.background     = ggplot2::element_blank(),
      legend.text         = ggplot2::element_text(size = ggplot2::rel(0.8)),
      strip.text          = ggplot2::element_text(size = 11)
    )
}

# Same as algatr::rda_plot() but reapplies SNP type colors so outliers use gea_outlier_col
# (algatr defaults use orange/purple in biplots, histograms, and Manhattan).
rda_plot_gea <- function(mod, rda_snps = NULL, pvalues = NULL, axes = "all",
                         biplot_axes = NULL, sig = 0.05, manhattan = NULL,
                         rdaplot = NULL, binwidth = NULL) {
  if (is.null(rdaplot)) rdaplot <- TRUE
  if (is.null(manhattan)) manhattan <- FALSE
  if (axes == "all") axes <- seq_len(ncol(mod$CCA$v))
  if (is.null(rda_snps)) {
    loadings <- vegan::scores(mod, choices = axes, display = "species")
    loadings <- loadings %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "SNP") %>%
      tidyr::pivot_longer(!SNP, names_to = "axis", values_to = "loading")
    print(rda_hist(loadings, binwidth = binwidth))
  }
  if (!is.null(rda_snps)) {
    tidy_list <- rda_ggtidy(mod, rda_snps, axes = axes)
    TAB_snps <- tidy_list[["TAB_snps"]]
    TAB_var <- tidy_list[["TAB_var"]]
    if (rdaplot) {
      if (length(axes) == 1L) {
        print(rda_hist(TAB_snps, binwidth = binwidth) + gea_scale_fill_snp_type())
      } else if (!is.null(biplot_axes)) {
        if (is.vector(biplot_axes)) {
          print(rda_biplot_gea(TAB_snps, TAB_var, biplot_axes = biplot_axes))
        }
        if (is.list(biplot_axes)) {
          lapply(biplot_axes, function(x) {
            print(rda_biplot_gea(TAB_snps, TAB_var, biplot_axes = x))
          })
        }
      } else {
        cb <- combn(length(axes), 2)
        if (!is.null(dim(cb))) {
          apply(cb, 2, function(x) {
            print(rda_biplot_gea(TAB_snps, TAB_var, biplot_axes = x))
          })
        } else {
          print(rda_biplot_gea(TAB_snps, TAB_var, biplot_axes = cb))
        }
      }
    }
    if (manhattan && !is.null(pvalues)) {
      print(rda_manhattan(TAB_snps, rda_snps, pvalues, sig = sig) + gea_scale_color_snp_type())
    }
  }
}

# LFMM Manhattan with x = cumulative genomic position and axis ticks labeled by chromosome
# (requires snp_map.rds from 00-Preprocessing.R: snp, chr, pos)
# merge_across_vars = TRUE (default): one panel like pRDA_z_manhattan_chr — one point per SNP,
#   y = max(-log10(adj. p)) across predictors, Outlier if any predictor passes sig.
# merge_across_vars = FALSE: facet_wrap(~var) (one mini-Manhattan per environmental variable).
lfmm_manhattanplot_chr <- function(
    df,
    snp_map,
    sig = 0.05,
    var = NULL,
    gap_frac = 0.02,
    merge_across_vars = TRUE
) {
  if (!is.null(var)) df <- df[df$var %in% var, ]
  df <- data.frame(df)
  m <- merge(df, snp_map, by = "snp", all.x = TRUE)
  miss <- is.na(m$chr) | is.na(m$pos)
  if (any(miss)) {
    warning("Dropping ", sum(miss), " SNP(s) with no chr/pos in snp_map.")
    m <- m[!miss, , drop = FALSE]
  }
  if (nrow(m) == 0L) stop("No rows left after merging LFMM results with snp_map.")

  facet_by_var <- FALSE
  n_vars <- if ("var" %in% names(m)) dplyr::n_distinct(m$var, na.rm = TRUE) else 1L
  if (merge_across_vars && "var" %in% names(m) && n_vars > 1L) {
    m <- m %>%
      dplyr::group_by(snp, chr, pos) %>%
      dplyr::summarise(
        neg_log10_p = max(-log10(.data$adjusted.pvalue), na.rm = TRUE),
        type = dplyr::if_else(
          any(.data$adjusted.pvalue < sig, na.rm = TRUE),
          "Outlier",
          "Non-outlier"
        ),
        .groups = "drop"
      )
    # max(-log10(p)) can be Inf if p = 0; keep for display
  } else {
    m$neg_log10_p <- -log10(m$adjusted.pvalue)
    m$type <- "Non-outlier"
    ok <- !is.na(m$adjusted.pvalue) & m$adjusted.pvalue < sig
    m$type[ok] <- "Outlier"
    facet_by_var <- !merge_across_vars && "var" %in% names(m) && dplyr::n_distinct(m$var, na.rm = TRUE) > 1L
  }

  chrs_sorted <- stringr::str_sort(unique(as.character(m$chr)), numeric = TRUE)
  chr_max <- setNames(
    vapply(chrs_sorted, function(ch) {
      max(m$pos[as.character(m$chr) == ch], na.rm = TRUE)
    }, numeric(1)),
    chrs_sorted
  )
  gap <- gap_frac * sum(chr_max, na.rm = TRUE) / max(1L, length(chrs_sorted))
  offset <- 0
  m$x_cum <- NA_real_
  for (ch in chrs_sorted) {
    ix <- as.character(m$chr) == ch
    m$x_cum[ix] <- offset + as.numeric(m$pos[ix])
    offset <- offset + chr_max[ch] + gap
  }
  chr_mid <- vapply(chrs_sorted, function(ch) {
    xs <- m$x_cum[as.character(m$chr) == ch]
    mean(range(xs, na.rm = TRUE))
  }, numeric(1))
  vline_x <- cumsum(chr_max[chrs_sorted] + gap) - gap
  vline_x <- head(vline_x, -1L)

  p <- ggplot2::ggplot(m, ggplot2::aes(x = .data$x_cum, y = .data$neg_log10_p, col = .data$type)) +
    ggplot2::geom_hline(
      yintercept = -log10(sig),
      linetype = "dashed",
      color = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(xintercept = vline_x, linetype = "dotted", color = "grey75", linewidth = 0.25) +
    ggplot2::geom_point(alpha = 0.75, pch = 16, ggplot2::aes(col = .data$type)) +
    ggplot2::scale_color_manual(
      values = c(`Non-outlier` = gea_nonoutlier_col, Outlier = gea_outlier_col),
      na.translate = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = chr_mid,
      labels = chrs_sorted,
      expand = ggplot2::expansion(mult = c(0.01, 0.01)),
      name = "Chromosome"
    ) +
    ggplot2::ylab(expression(-log[10](italic(p)))) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "SNP type")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position       = "right",
      legend.background     = ggplot2::element_blank(),
      panel.grid            = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      plot.background       = ggplot2::element_blank(),
      panel.background      = ggplot2::element_blank(),
      legend.text           = ggplot2::element_text(size = ggplot2::rel(0.8)),
      axis.text.x           = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (facet_by_var) {
    p <- p +
      ggplot2::facet_wrap(~var, nrow = length(unique(m$var))) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 11))
  } else {
    p <- p +
      ggplot2::labs(
        title = "LFMM genome scan",
        subtitle = if (merge_across_vars && "var" %in% names(df) && dplyr::n_distinct(df$var, na.rm = TRUE) > 1L) {
          "Per SNP: strongest signal across environmental variables (max -log10 adj. p)"
        } else {
          NULL
        }
      )
  }
  p
}

# RDA "Manhattan" by chromosome when you used z-score outlier detection (no rdadapt p-values).
# For each SNP: species scores are standardized within each constrained axis (mean 0, sd 1),
# then y = max |Z| across axes. Points in sig_z_obj$rda_snps (from rda_getoutliers(..., z=))
# are colored as outliers. Horizontal line = z_threshold (same rule of thumb as algatr).
rda_z_manhattanplot_chr <- function(mod, snp_map, sig_z_obj, z_threshold = 3, gap_frac = 0.02) {
  naxes <- ncol(mod$CCA$v)
  if (is.null(naxes) || naxes < 1L) {
    stop("No constrained CCA axes in mod (mod$CCA$v is empty).")
  }
  load_rda <- vegan::scores(mod, choices = seq_len(naxes), display = "species", scaling = "none")
  snp_names <- rownames(load_rda)
  Z_mat <- apply(as.matrix(load_rda), 2L, function(x) {
    sdx <- stats::sd(x)
    if (!is.finite(sdx) || sdx == 0) rep(0, length(x)) else (x - mean(x)) / sdx
  })
  z_max <- if (ncol(Z_mat) == 1L) abs(as.vector(Z_mat)) else apply(abs(Z_mat), 1L, max)
  df <- data.frame(snp = snp_names, z_max = as.numeric(z_max), stringsAsFactors = FALSE)
  if (is.data.frame(sig_z_obj) && "rda_snps" %in% names(sig_z_obj)) {
    out_snps <- unique(as.character(sig_z_obj$rda_snps))
  } else if (is.list(sig_z_obj) && !is.null(sig_z_obj$rda_snps)) {
    out_snps <- unique(as.character(sig_z_obj$rda_snps))
  } else {
    out_snps <- character(0)
  }
  df$type <- ifelse(df$snp %in% out_snps, "Outlier", "Non-outlier")
  m <- merge(df, snp_map, by = "snp", all.x = TRUE)
  miss <- is.na(m$chr) | is.na(m$pos)
  if (any(miss)) {
    warning("Dropping ", sum(miss), " SNP(s) with no chr/pos in snp_map.")
    m <- m[!miss, , drop = FALSE]
  }
  if (nrow(m) == 0L) stop("No rows left after merging RDA Z-scores with snp_map.")
  chrs_sorted <- stringr::str_sort(unique(as.character(m$chr)), numeric = TRUE)
  chr_max <- setNames(
    vapply(chrs_sorted, function(ch) {
      max(m$pos[as.character(m$chr) == ch], na.rm = TRUE)
    }, numeric(1)),
    chrs_sorted
  )
  gap <- gap_frac * sum(chr_max, na.rm = TRUE) / max(1L, length(chrs_sorted))
  offset <- 0
  m$x_cum <- NA_real_
  for (ch in chrs_sorted) {
    ix <- as.character(m$chr) == ch
    m$x_cum[ix] <- offset + as.numeric(m$pos[ix])
    offset <- offset + chr_max[ch] + gap
  }
  chr_mid <- vapply(chrs_sorted, function(ch) {
    xs <- m$x_cum[as.character(m$chr) == ch]
    mean(range(xs, na.rm = TRUE))
  }, numeric(1))
  vline_x <- cumsum(chr_max[chrs_sorted] + gap) - gap
  vline_x <- head(vline_x, -1L)

  ggplot2::ggplot(m, ggplot2::aes(x = .data$x_cum, y = .data$z_max, col = .data$type)) +
    ggplot2::geom_hline(yintercept = z_threshold, linetype = "dashed", color = "black", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = vline_x, linetype = "dotted", color = "grey75", linewidth = 0.25) +
    ggplot2::geom_point(alpha = 0.75, pch = 16, ggplot2::aes(col = .data$type)) +
    ggplot2::scale_color_manual(
      values = c(`Non-outlier` = gea_nonoutlier_col, Outlier = gea_outlier_col),
      na.translate = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = chr_mid,
      labels = chrs_sorted,
      expand = ggplot2::expansion(mult = c(0.01, 0.01)),
      name = "Chromosome"
    ) +
    ggplot2::ylab("max |Z| (RDA loading, scaled per axis)") +
    ggplot2::labs(
      title = "RDA genome scan (Z-scores)",
      subtitle = paste0("Outliers = SNPs in z-based rda_getoutliers; dashed y = ", z_threshold)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "SNP type")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position       = "right",
      legend.background     = ggplot2::element_blank(),
      panel.grid            = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      plot.background       = ggplot2::element_blank(),
      panel.background      = ggplot2::element_blank(),
      legend.text           = ggplot2::element_text(size = ggplot2::rel(0.8)),
      axis.text.x           = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

# Readable labels for WorldClim bioclim variables (constraint / biplot arrows in RDA)
bioclim_biolabels <- c(
  bio1  = "BIO1: annual mean temperature",
  bio2  = "BIO2: mean diurnal range",
  bio3  = "BIO3: isothermality",
  bio4  = "BIO4: temperature seasonality",
  bio5  = "BIO5: max temp. of warmest month",
  bio6  = "BIO6: min temp. of coldest month",
  bio7  = "BIO7: temperature annual range",
  bio8  = "BIO8: mean temp. of wettest quarter",
  bio9  = "BIO9: mean temp. of driest quarter",
  bio10 = "BIO10: mean temp. of warmest quarter",
  bio11 = "BIO11: mean temp. of coldest quarter",
  bio12 = "BIO12: annual precipitation",
  bio13 = "BIO13: precip. of wettest month",
  bio14 = "BIO14: precip. of driest month",
  bio15 = "BIO15: precip. seasonality (CV)",
  bio16 = "BIO16: precip. of wettest quarter",
  bio17 = "BIO17: precip. of driest quarter",
  bio18 = "BIO18: precip. of warmest quarter",
  bio19 = "BIO19: precip. of coldest quarter"
)

# MMRR coefficient / importance panels (data from Scripts/03-IBD-IBE.R, not hardcoded)
mmrr_var_label <- function(var, labels_map = bioclim_biolabels) {
  v <- as.character(var)
  if (v %in% names(labels_map)) return(unname(labels_map[[v]]))
  if (v == "geodist") return("Geographic distance")
  v
}

load_mmrr_coeff_table <- function(ibd_dir) {
  cache_path <- file.path(ibd_dir, "_step_cache", "mmrr_YX_models.rds")
  if (file.exists(cache_path)) {
    bundle <- readRDS(cache_path)
    if (is.list(bundle) && !is.null(bundle$value)) bundle <- bundle$value
    coeff <- bundle$results_full$coeff_df
    if (is.list(coeff) && !is.null(coeff$var)) {
      return(data.frame(
        var = as.character(coeff$var),
        estimate = as.numeric(coeff$estimate),
        lower_ci = as.numeric(coeff[["95% Lower"]]),
        upper_ci = as.numeric(coeff[["95% Upper"]]),
        p_value = as.numeric(coeff$p),
        stringsAsFactors = FALSE
      ))
    }
  }
  txt_path <- file.path(ibd_dir, "results", "mmrr_results_full.txt")
  if (!file.exists(txt_path)) {
    stop(
      "MMRR results not found. Run Scripts/03-IBD-IBE.R first.\n  Tried: ",
      cache_path, "\n  and: ", txt_path
    )
  }
  raw <- utils::read.delim(txt_path, check.names = FALSE, stringsAsFactors = FALSE)
  names(raw) <- tolower(gsub("[^a-z0-9]+", "_", names(raw)))
  var_col <- intersect(names(raw), c("var", "variable"))[1]
  est_col <- intersect(names(raw), c("estimate", "coefficient"))[1]
  p_col <- intersect(names(raw), c("p", "p_value"))[1]
  lo_col <- grep("lower", names(raw), value = TRUE)[1]
  hi_col <- grep("upper", names(raw), value = TRUE)[1]
  if (any(is.na(c(var_col, est_col, p_col, lo_col, hi_col)))) {
    stop("Could not parse columns in ", txt_path)
  }
  out <- data.frame(
    var = raw[[var_col]],
    estimate = as.numeric(gsub("\u2212", "-", raw[[est_col]], fixed = TRUE)),
    lower_ci = as.numeric(gsub("\u2212", "-", raw[[lo_col]], fixed = TRUE)),
    upper_ci = as.numeric(gsub("\u2212", "-", raw[[hi_col]], fixed = TRUE)),
    p_value = as.numeric(raw[[p_col]]),
    stringsAsFactors = FALSE
  )
  out[!is.na(out$estimate) & !grepl("^R-Squared|^F-", out$var, ignore.case = TRUE), , drop = FALSE]
}

plot_mmrr_coefficients <- function(mmrr_df, out_path) {
  df <- mmrr_df[mmrr_df$var != "Intercept", , drop = FALSE]
  if (!nrow(df)) stop("No MMRR predictors to plot (only Intercept?).")
  df$var_label <- vapply(df$var, mmrr_var_label, character(1))
  df$var_label <- factor(df$var_label, levels = df$var_label[order(abs(df$estimate))])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$var_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$lower_ci, xmax = .data$upper_ci),
      height = 0.2
    ) +
    ggplot2::geom_point(size = 3, color = "steelblue") +
    ggplot2::labs(
      x = "Standardized coefficient (95% CI)",
      y = "Predictor variable",
      title = "MMRR standardized coefficients"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10)
    )
  ggplot2::ggsave(out_path, p, width = 6, height = max(4, 0.35 * nrow(df)), dpi = 300)
  invisible(p)
}

plot_mmrr_relative_importance <- function(mmrr_df, out_path) {
  df <- mmrr_df[mmrr_df$var != "Intercept", , drop = FALSE]
  if (!nrow(df)) stop("No MMRR predictors to plot (only Intercept?).")
  df$abs_estimate <- abs(df$estimate)
  df$rel_importance <- df$abs_estimate / sum(df$abs_estimate) * 100
  df$var_label <- vapply(df$var, mmrr_var_label, character(1))
  df$var_label <- factor(df$var_label, levels = df$var_label[order(df$rel_importance)])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rel_importance, y = .data$var_label, fill = .data$var_label)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      x = "Relative importance (%)",
      y = "",
      title = "MMRR relative predictor importance"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10)
    )
  ggplot2::ggsave(out_path, p, width = 6, height = max(4, 0.35 * nrow(df)), dpi = 300)
  invisible(p)
}

# Biplot matching algatr::rda_biplot scaling but with bioclim_labels on constraint arrows
rda_biplot_bioclim_labeled <- function(mod, rda_snps, biplot_axes = c(1, 2),
                                       labels_map = bioclim_biolabels) {
  axes <- seq_len(max(max(biplot_axes), 2L))
  tidy_list <- rda_ggtidy(mod, rda_snps = rda_snps, axes = axes)
  TAB_snps <- tidy_list[["TAB_snps"]]
  TAB_var  <- tidy_list[["TAB_var"]]
  xax <- paste0("RDA", biplot_axes[1])
  yax <- paste0("RDA", biplot_axes[2])
  TAB_snps_sub <- TAB_snps[, c(xax, yax, "type")]
  colnames(TAB_snps_sub) <- c("x", "y", "type")
  TAB_var_sub <- TAB_var[, c(xax, yax), drop = FALSE]
  colnames(TAB_var_sub) <- c("x", "y")
  TAB_var_sub$x <- TAB_var_sub$x * max(TAB_snps_sub$x) / stats::quantile(TAB_var_sub$x, na.rm = TRUE)[4]
  TAB_var_sub$y <- TAB_var_sub$y * max(TAB_snps_sub$y) / stats::quantile(TAB_var_sub$y, na.rm = TRUE)[4]
  rn <- rownames(TAB_var_sub)
  var_lab <- rn
  if (length(labels_map)) {
    ii <- match(rn, names(labels_map))
    var_lab <- ifelse(!is.na(ii), unname(labels_map)[ii], rn)
  }
  sn_layers <- gea_geom_snp_biplot(TAB_snps_sub)
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.8), linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.8), linewidth = 0.6) +
    sn_layers[[1]] +
    sn_layers[[2]] +
    gea_scale_color_snp_type() +
    geom_segment(
      data = TAB_var_sub,
      aes(xend = x, yend = y, x = 0, y = 0),
      colour = "black",
      linewidth = 0.15,
      arrow = arrow(length = unit(0.02, "npc"))
    ) +
    geom_text_repel(
      data = TAB_var_sub,
      aes(x = x, y = y, label = var_lab),
      size = 3.2,
      max.overlaps = Inf,
      segment.size = 0.2
    ) +
    xlab(xax) + ylab(yax) +
    guides(color = guide_legend(title = "SNP type")) +
    theme_bw(base_size = 11) +
    theme(
      panel.background   = element_blank(),
      legend.background  = element_blank(),
      panel.grid         = element_blank(),
      plot.background    = element_blank(),
      legend.text        = element_text(size = rel(0.8)),
      strip.text         = element_text(size = 11)
    )
}

## Pipeline cache + Paths (Outputs/…) — align with Scripts/00-Preprocessing.R
._lgpr <- Sys.getenv("LGP_PROJECT_ROOT", "~/Desktop/LandscapeGenomicsPipeline")
options(lgp.project_root = sub("/+$", "", path.expand(getOption("lgp.project_root", ._lgpr))))
suppressPackageStartupMessages(base::source(
  base::file.path(getOption("lgp.project_root"), "Scripts", "lgp_pipeline_cache.R"),
  encoding = "UTF-8"
))
rm(._lgpr)

## ------------------------------------------------
## 0. Checks and paths
## ------------------------------------------------

data_dir <- lgp_project_root()
in_dir <- lgp_preprocess_dir()
out_dir <- lgp_outputs_step_dir("04-gea")
plot_dir <- file.path(out_dir, "plots")
results_dir <- file.path(out_dir, "results")
gea_step_cache <- file.path(out_dir, "_step_cache")
dir.create(gea_step_cache, recursive = TRUE, showWarnings = FALSE)
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
load_if_missing("str_dos", "str_dos.rds")

snp_map_path <- file.path(in_dir, "snp_map.rds")
if (!file.exists(snp_map_path)) {
  alt <- file.path(data_dir, "Outputs", "00-preprocessing", "snp_map.rds")
  if (file.exists(alt)) snp_map_path <- alt
}
if (!file.exists(snp_map_path)) {
  stop(
    "Missing snp_map.rds (re-run 00-Preprocessing.R to create it). Tried:\n  ",
    file.path(in_dir, "snp_map.rds"), "\n  ", file.path(data_dir, "Outputs", "00-preprocessing", "snp_map.rds")
  )
}
snp_map <- readRDS(snp_map_path)

#######################################################################
## 01. Redundancy analysis (RDA)
#######################################################################

# Scale the genomic data using Hellinger transformation
gen_h <- decostand(str_dos, "hellinger")

# Run a simple RDA with no variable selection
mod_full <- rda_run(gen_h, env_scaled, model = "full")
writeLines(capture.output(mod_full$call, summary(mod_full), RsquareAdj(mod_full)),
           file.path(results_dir, "rda_full_summary.txt"))

# Full RDA: per-term (sequential) tests for each environmental predictor
rda_full_terms_perm <- 999L
a_full_terms <- anova(mod_full, by = "terms", permutations = rda_full_terms_perm)
writeLines(
  capture.output(a_full_terms),
  file.path(results_dir, "rda_full_anova_by_terms.txt")
)
constr_inertia <- sum(mod_full$CCA$eig, na.rm = TRUE)
term_only <- !grepl("^Residual$", rownames(a_full_terms), ignore.case = TRUE)
prop_constr <- rep(NA_real_, nrow(a_full_terms))
prop_constr[term_only] <- a_full_terms$Variance[term_only] / constr_inertia
rda_full_terms_df <- data.frame(
  Term                 = rownames(a_full_terms),
  Df                   = a_full_terms$Df,
  Variance             = a_full_terms$Variance,
  F                    = a_full_terms$F,
  Pr_F                 = a_full_terms$`Pr(>F)`,
  Prop_of_constrained  = prop_constr,
  row.names            = NULL
)
write.csv(
  rda_full_terms_df,
  file.path(results_dir, "rda_full_anova_terms_table.csv"),
  row.names = FALSE
)

# Run RDA with variable selection
mod_best <- rda_run(gen_h, env_scaled,
                    model = "best",
                    Pin = 0.05,
                    R2permutations = 1000,
                    R2scope = TRUE
)
writeLines(capture.output(mod_best$call, mod_best$anova, RsquareAdj(mod_best)),
           file.path(results_dir, "rda_best_summary.txt"))

# Partial RDA with geography as covariable
mod_pRDA_geo <- rda_run(gen_h, env_scaled, coord,
                        model = "full",
                        correctGEO = TRUE,
                        correctPC = FALSE
)
writeLines(capture.output(anova(mod_pRDA_geo), RsquareAdj(mod_pRDA_geo), summary(mod_pRDA_geo)),
           file.path(results_dir, "rda_partial_geo_summary.txt"))

# Partial RDA with genetic structure and geography as covariables
mod_pRDA_gs <- rda_run(gen_h, env_scaled, coord,
                       model = "full",
                       correctGEO = TRUE,
                       correctPC = TRUE,
                       nPC = 2
)
writeLines(capture.output(summary(mod_pRDA_gs)), file.path(results_dir, "rda_partial_gs_summary.txt"))

png(file.path(plot_dir, "rda_partial_full_gs_plot.png"), width = 1200, height = 1200, res = 200)
plot(mod_pRDA_gs)
dev.off()

# Variance partitioning
varpart <- rda_varpart(gen_h, env_scaled, coord,
                       Pin = 0.05, R2permutations = 1000,
                       R2scope = TRUE, nPC = 2
)
writeLines(capture.output(str(varpart)),
           file.path(results_dir, "rda_varpart_table.txt"))

rda_varpart_table(varpart, call_col = TRUE)

# Partial RDA with variable selection (for outlier detection)
mod_pRDA <- rda_run(gen_h, env_scaled, model = "best", correctPC = TRUE, nPC = 2)
writeLines(capture.output(mod_pRDA$anova), file.path(results_dir, "rda_partial_best_anova.txt"))

# Identifying outliers using the Z-scores method
png(file.path(plot_dir, "prda_best_plot_axes.png"), width = 1600, height = 1200, res = 200)
rda_plot_gea(mod_pRDA, axes = "all", binwidth = 20)
dev.off()

rda_sig_z <- rda_getoutliers(mod_pRDA, naxes = "all", outlier_method = "z", z = 3, plot = FALSE)
writeLines(paste0("Number of RDA outlier SNPs (z): ", length(rda_sig_z$rda_snps)),
           file.path(results_dir, "prda_best_outlier_count.txt"))

write.csv(rda_sig_z, file.path(results_dir, "best_pRDA_sigZ_SNPS.csv"))

##Formatting for later analyses 
# split into CHR and BP
sp <- strsplit(rda_sig_z$rda_snps, "_")

# build data.frame
df <- data.frame(
  SNP = rda_sig_z$rda_snps,
  CHR = sapply(sp, `[`, 1),
  BP  = sapply(sp, `[`, 2),
  stringsAsFactors = FALSE
)

# write to CSV
write.csv(df, "best_pRDA_sigZ_SNPS.csv", row.names = FALSE)

pRDA_snp_positions <- gea_chr_bp_from_snp_markers(rda_sig_z$rda_snps)
utils::write.table(
  pRDA_snp_positions,
  file = file.path(results_dir, "pRDA_snp_positions.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rda_outlier_snps <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x) && "rda_snps" %in% names(x)) return(x$rda_snps)
  if (is.list(x) && !is.null(x$rda_snps)) return(x$rda_snps)
  NULL
}
rda_outlier_pvals <- function(x) {
  if (is.null(x) || !is.list(x)) return(NULL)
  x$pvalues
}

# Per-SNP p-values for rda_manhattan(). algatr::rda_plot(..., manhattan = TRUE) only draws
# when pvalues is non-NULL (`if (manhattan & !is.null(pvalues)) print(rda_manhattan(...))`).
# rda_outlier_pvals() is NULL when: sig object is NULL; or not a list; or $pvalues missing
# (e.g. z-score path). This rebuilds FDR p-values from $rdadapt when present, else rdadapt().
rda_pvals_for_manhattan <- function(mod, sig_p, p_adj = "fdr") {
  if (is.list(sig_p) && !is.null(sig_p$pvalues)) {
    return(sig_p$pvalues)
  }
  if (is.list(sig_p) && !is.null(sig_p$rdadapt)) {
    return(stats::p.adjust(sig_p$rdadapt$p.values, method = p_adj))
  }
  ncca <- ncol(mod$CCA$v)
  if (is.null(ncca) || ncca < 2L) {
    warning(
      "RDA Manhattan needs rdadapt p-values; requires >= 2 constrained CCA axes. Found: ",
      if (is.null(ncca)) 0L else ncca
    )
    return(NULL)
  }
  ra <- rdadapt(mod, ncca)
  stats::p.adjust(ra$p.values, method = p_adj)
}

rda_sig_p_best <- rda_getoutliers(mod_best, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.05, plot = FALSE)
rda_sig_p_gs <- rda_getoutliers(mod_pRDA_gs, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.05, plot = FALSE)
rda_sig_p_prda <- rda_getoutliers(mod_pRDA, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.05, plot = FALSE)

# Visualizing RDA results with rda_plot_gea() (algatr rda_plot + #4E0707 outliers)
png(file.path(plot_dir, "fullRDA_bioplot_best.png"), width = 1600, height = 1200, res = 200)
rda_plot_gea(mod_best, rda_outlier_snps(rda_sig_p_best), biplot_axes = c(1, 2), rdaplot = TRUE, manhattan = FALSE)
dev.off()

png(file.path(plot_dir, "pRDA_bioplot_gs.png"), width = 1600, height = 1200, res = 200)
rda_plot_gea(mod_pRDA_gs, rda_outlier_snps(rda_sig_p_gs), biplot_axes = c(1, 2), rdaplot = TRUE, manhattan = FALSE)
dev.off()

# RDA plot with outlier SNPs highlighted
png(file.path(plot_dir, "rda_plot_outliers.png"), width = 1600, height = 1200, res = 200)
rda_plot_gea(mod_pRDA, rda_sig_z$rda_snps, rdaplot = TRUE, manhattan = FALSE, binwidth = 0.01)
dev.off()

#png(file.path(plot_dir, "prda_best_biplot_bioclim_labels_partial_pc.png"), width = 1600, height = 1200, res = 200)
#print(rda_biplot_bioclim_labeled(mod_pRDA, rda_sig_z$rda_snps, biplot_axes = c(1, 2)))
#dev.off()

# Manhattan plot (per-SNP p-values from rdadapt). Do not rely on rda_plot_gea(..., manhattan=TRUE)
# alone: it draws nothing unless pvalues is non-NULL (see algatr::rda_plot).
png(file.path(plot_dir, "pRDA_best_manhattan.png"), width = 1600, height = 1200, res = 200)
{
  p_mv <- rda_pvals_for_manhattan(mod_pRDA, rda_sig_p_prda, p_adj = "fdr")
  sn_mv <- rda_outlier_snps(rda_sig_p_prda)
  if (is.null(sn_mv)) sn_mv <- character(0)
  if (!is.null(p_mv) && length(p_mv) > 0L) {
    axes <- seq_len(ncol(mod_pRDA$CCA$v))
    tl <- rda_ggtidy(mod_pRDA, sn_mv, axes = axes)
    if (length(p_mv) != nrow(tl$TAB_snps)) {
      warning(
        "p-value length (", length(p_mv), ") does not match TAB_snps (", nrow(tl$TAB_snps),
        "); not drawing RDA Manhattan."
      )
      plot.new()
      title(main = "RDA Manhattan: length mismatch (p-values vs SNP table)")
    } else {
      print(rda_manhattan(tl$TAB_snps, sn_mv, p_mv, sig = 0.05) + gea_scale_color_snp_type())
    }
  } else {
    plot.new()
    title(main = "RDA Manhattan: no p-values (constrained axes < 2?)")
  }
}
dev.off()

# Chromosome-position plot for z-score RDA outliers (no p-values / rdadapt required).
# Uses the same cumulative-x layout as LFMM Manhattan; y = max |Z| across constrained axes.
png(file.path(plot_dir, "pRDA_z_manhattan_chr.png"), width = 1600, height = 800, res = 200)
print(rda_z_manhattanplot_chr(mod_pRDA, snp_map, rda_sig_z, z_threshold = 3))
dev.off()

# Interpreting RDA results
# Extract genotypes for outlier SNPs
rda_snps_z <- rda_sig_z$rda_snps
rda_gen_z <- gen_h[, rda_snps_z]
write.csv(rda_gen_z, file.path(results_dir, "pRDA_genotypes.csv"), row.names = TRUE)

# Run correlation test
cor_df_z <- rda_cor(rda_gen_z, env_scaled)

# Save RDA correlation table (p < 0.05) to file
rda_sig_table <- rda_table(cor_df_z, sig = 0.05, order = TRUE)
write.csv(rda_sig_table, file.path(results_dir, "rda_p0.05.csv"), row.names = TRUE)

#######################################################################
## 02. Latent factor mixed models (LFMM)
#######################################################################

# K selection (run for reference; use chosen K in lfmm_run)
# select_K(gen_h, K_selection = "tracy_widom", criticalpoint = 2.0234)
# select_K(gen_h, K_selection = "quick_elbow", low = 0.08, max.pc = 0.90)
# select_K(gen_h, K_selection = "tess", coords = coord, Kvals = 1:10)
# select_K(gen_h, K_selection = "find_clusters", perc.pca = 90, max.n.clust = 10)

lfmm_ridge_cache <- file.path(gea_step_cache, "lfmm_ridge_fit_wrapped.rds")
str_pre <- file.path(in_dir, "str_dos.rds")
lfmm_force <- isTRUE(getOption("lgp.gea_force_lfmm", FALSE))

ridge_results <- if (exists("ridge_results", inherits = TRUE) && !lfmm_force && !lgp_force_rerun()) {
  get("ridge_results", inherits = TRUE)
} else if (lfmm_force || lgp_force_rerun()) {
  rr <- lfmm_run(gen_h, env_scaled, K = 2, lfmm_method = "ridge")
  lgp_atomic_saveRDS(rr, lfmm_ridge_cache)
  rr
} else {
  lgp_rds_cached(
    lfmm_ridge_cache,
    dep_paths = str_pre,
    meta      = list(K = 2L, lfmm_method = "ridge"),
    label     = "LFMM ridge_results",
    compute   = function() lfmm_run(gen_h, env_scaled, K = 2, lfmm_method = "ridge")
  )
}

# lasso_results <- lfmm_run(gen_h, env_scaled, K = 2, lfmm_method = "lasso")

if (!exists("ridge_results", inherits = TRUE) || is.null(ridge_results[["df"]])) {
  stop(
    paste0(
      "LFMM ridge_results is missing or has no $df component.\n",
      "Run the lfmm_run() block above in this script, or after a full run restore with:\n",
      "  ridge_results <- readRDS(\"", lfmm_ridge_cache, "\")"
    ),
    call. = FALSE
  )
}

# Save LFMM tables to files
# Significant associations: algatr::lfmm_table() renders a {gt} object, not a plain
# data.frame, so we subset ridge_results$df the same way lfmm_table(sig_only = TRUE) does.
lfmm_tab_adj_p <- suppressWarnings(as.numeric(getOption("lgp.lfmm_adj_p_threshold", 0.05))[1])
if (!is.finite(lfmm_tab_adj_p) || lfmm_tab_adj_p <= 0 || lfmm_tab_adj_p > 1) {
  lfmm_tab_adj_p <- 0.05
}
lfmm_ridge_tab <- dplyr::as_tibble(ridge_results$df) %>%
  dplyr::filter(!is.na(.data$adjusted.pvalue), .data$adjusted.pvalue < lfmm_tab_adj_p) %>%
  dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.)))
if (nrow(lfmm_ridge_tab) == 0L) {
  warning(
    "No LFMM ridge rows with adjusted.pvalue < ",
    lfmm_tab_adj_p,
    " — sign_lfmm_ridge_table will be empty (lfmm_table default sig_only = TRUE).",
    call. = FALSE
  )
}

write.table(
  lfmm_ridge_tab,
  file = file.path(results_dir, "sign_lfmm_ridge_table.csv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

snm_lf <- intersect(names(lfmm_ridge_tab), c("snp", "SNP"))
if (length(snm_lf) != 1L) {
  stop(
    "sign_lfmm_ridge_table: need column snp in ridge_results$df (have: ",
    paste(names(ridge_results$df), collapse = ", "),
    "); cannot build lfmm_snp_positions.",
    call. = FALSE
  )
}
if (nrow(lfmm_ridge_tab) > 0L) {
  lfmm_snp_positions <- gea_chr_bp_from_snp_markers(unique(lfmm_ridge_tab[[snm_lf]]))
  utils::write.csv(lfmm_snp_positions, file.path(results_dir, "lfmm_snp_positions.csv"), row.names = FALSE)
}
# Save  SNP tables as CSV for downstream use
write.csv(ridge_results$df, file.path(results_dir, "lfmm_ridge_results.csv"), row.names = TRUE)

# split into CHR and BP
sp <- strsplit(ridge_results$df$snp, "_")

# build data.frame
df <- data.frame(
  SNP = ridge_results$df$snp,
  CHR = sapply(sp, `[`, 1),
  BP  = sapply(sp, `[`, 2),
  stringsAsFactors = FALSE
)

# write to CSV
write.csv(df, "sign_lfmm_ridge_list.csv", row.names = FALSE)


# QQ plot
png(file.path(plot_dir, "lfmm_qqplot_ridge.png"), width = 1200, height = 1200, res = 200)
lfmm_qqplot(ridge_results$df)
dev.off()

# Manhattan plots (x-axis: cumulative position with chromosome labels; needs snp_map.rds)
png(file.path(plot_dir, "lfmm_manhattan_ridge.png"), width = 1600, height = 800, res = 200)
p_lfmm_manhattan <- lfmm_manhattanplot_chr(ridge_results$df, snp_map, sig = 0.05)
print(p_lfmm_manhattan)
dev.off()

#######################################################################
## 03. MMRR coefficient and importance panels (from 03-IBD-IBE.R)
#######################################################################

ibd_dir <- lgp_outputs_step_dir("03-ibd-ibe", root = data_dir)
mmrr_results <- load_mmrr_coeff_table(ibd_dir)
mmrr_results$significant <- ifelse(
  mmrr_results$p_value < 0.01, "**",
  ifelse(mmrr_results$p_value < 0.05, "*", "")
)
utils::write.table(
  mmrr_results,
  file.path(results_dir, "mmrr_coefficient_table.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
plot_mmrr_coefficients(
  mmrr_results,
  file.path(plot_dir, "mmrr_panel_coefficients.png")
)
plot_mmrr_relative_importance(
  mmrr_results,
  file.path(plot_dir, "mmrr_panel_importance.png")
)
message(
  "[04-gea] MMRR panels from ", ibd_dir,
  " -> mmrr_panel_coefficients.png, mmrr_panel_importance.png"
)
