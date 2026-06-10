#######################################################################
# GO enrichment (topGO) from LD-linked genes — pRDA vs LFMM windows
#
# Inputs:
#   • Outputs/…/05-linkage-decay/results/pRDA_linked_genes.tsv
#   • Outputs/…/05-linkage-decay/results/LFMM_linked_genes.tsv
#     (bedtools intersect -wao; gene IDs parsed from GFF attributes column)
#   • InputData/gene_go_terms_wide.tsv       → gene → GO mapping (topGO / annFUN.gene2GO)
#   • InputData/gene_functional_descriptions_wide.tsv (tab, no header: id \\t description)
#
# Requires Bioconductor packages topGO and GO.db:
#   install.packages("BiocManager")
#   BiocManager::install(c("topGO", "GO.db"))
# Heatmap (optional): ggplot2
#   install.packages("ggplot2")
#######################################################################

# Cached paths (Scripts/lgp_pipeline_cache.R).
._lgpr <- Sys.getenv("LGP_PROJECT_ROOT", "~/Desktop/LandscapeGenomicsPipeline")
options(lgp.project_root = sub("/+$", "", path.expand(getOption("lgp.project_root", ._lgpr))))
suppressPackageStartupMessages(base::source(
  base::file.path(getOption("lgp.project_root"), "Scripts", "lgp_pipeline_cache.R"),
  encoding = "UTF-8"
))
rm(._lgpr)

data_dir <- lgp_project_root()
ld_results_dir <- file.path(lgp_outputs_base(data_dir), "05-linkage-decay", "results")

go_step_dir <- lgp_outputs_step_dir("06-go-enrichment")
results_dir <- file.path(go_step_dir, "results")
go_plot_dir <- file.path(go_step_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_plot_dir, recursive = TRUE, showWarnings = FALSE)

linked_prda <- file.path(ld_results_dir, "pRDA_linked_genes.tsv")
linked_lfmm <- file.path(ld_results_dir, "LFMM_linked_genes.tsv")

go_terms_wide <- file.path(data_dir, "InputData", "gene_go_terms_wide.tsv")
go_desc_wide <- file.path(data_dir, "InputData", "gene_functional_descriptions_wide.tsv")

# Minimum genes (after intersect with GO universe) before running Fisher tests.
go_min_genes <- suppressWarnings(as.integer(getOption("lgp.go_min_genes_of_interest", 5L))[1L])
if (!is.finite(go_min_genes) || go_min_genes < 2L) go_min_genes <- 5L

go_ontologies <- getOption("lgp.go_topgo_ontologies", c("MF", "BP", "CC"))
go_ontologies <- toupper(trimws(as.character(go_ontologies)))
go_ontologies <- unique(go_ontologies[go_ontologies %in% c("MF", "BP", "CC")])

## topGO creates GOMFTerm / GOBPTerm / GOCCTerm when the package is *attached*
## (see ?topGO — .onAttach runs groupGOTerms()). Using only topGO::... loads the namespace
## but does not attach the package, so internal get("GOMFTerm") fails.
if (!requireNamespace("GO.db", quietly = TRUE)) {
  stop(
    "Package GO.db is required with topGO.\n",
    "  BiocManager::install(\"GO.db\")\n",
    "Or: BiocManager::install(c(\"topGO\", \"GO.db\"))",
    call. = FALSE
  )
}
if (!requireNamespace("topGO", quietly = TRUE)) {
  stop(
    "Package topGO is required.\n",
    "  BiocManager::install(\"topGO\")",
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  suppressWarnings(library(topGO, quietly = TRUE, warn.conflicts = FALSE))
})
if (!exists("GOMFTerm", inherits = TRUE)) {
  topGO::groupGOTerms()
}
if (!exists("GOMFTerm", inherits = TRUE)) {
  stop(
    "topGO GO term environments missing after attach (GOMFTerm). Reinstall/update:\n",
    "  BiocManager::install(c(\"topGO\", \"GO.db\"))",
    call. = FALSE
  )
}

#' Pull unique gene IDs from bedtools intersect -wao output (BED + GFF + overlap).
#' Uses the GFF attributes column (second-to-last column).
lgp_linked_genes_bedtools_gene_ids <- function(tsv_path) {
  if (!isTRUE(file.exists(tsv_path))) {
    return(character(0))
  }
  hdr_line <- trimws(readLines(tsv_path, n = 1L, warn = FALSE))
  has_hdr <- nzchar(hdr_line) && grepl("^bed_chr\t", hdr_line)
  dt <- utils::read.delim(
    tsv_path,
    header = has_hdr,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  if (!ncol(dt) || nrow(dt) == 0L) {
    return(character(0))
  }
  if ("gene_id" %in% names(dt)) {
    ids <- trimws(as.character(dt[["gene_id"]]))
    ids <- unique(ids[!is.na(ids) & nzchar(ids)])
    return(ids)
  }
  attr_col_idx <- ncol(dt) - 1L
  if (attr_col_idx < 1L) {
    return(character(0))
  }
  attrs <- dt[[attr_col_idx]]
  attrs <- trimws(as.character(attrs))
  attrs <- attrs[!is.na(attrs) & nzchar(attrs)]
  ids <- vapply(attrs, function(a) {
    m <- regexec("gene_id=([^;]+)", a)
    rm <- regmatches(a, m)[[1]]
    if (length(rm) >= 2L) {
      return(trimws(rm[[2]]))
    }
    m2 <- regexec("(?:^|;)ID=([^;]+)", a)
    rm2 <- regmatches(a, m2)[[1]]
    if (length(rm2) >= 2L) trimws(rm2[[2]]) else ""
  }, character(1))
  ids <- ids[nzchar(ids)]
  unique(ids)
}

#' Same mapping format topGO expects (gene ID → GO IDs character vector).
lgp_read_gene2go_wide <- function(path) {
  df <- utils::read.delim(path, header = FALSE, sep = "\t", quote = "",
                          stringsAsFactors = FALSE, comment.char = "")
  if (ncol(df) < 2L) {
    stop("gene_go_terms_wide.tsv: need >= 2 columns (gene \\t GO list).", call. = FALSE)
  }
  gids <- trimws(as.character(df[[1]]))
  lst <- strsplit(trimws(as.character(df[[2]])), ",", fixed = TRUE)
  names(lst) <- gids
  lapply(lst, function(z) unique(trimws(z[nzchar(trimws(z))])))
}

#' Functional descriptions: tab-separated, no header (gene \\t text).
lgp_read_func_desc_wide <- function(path) {
  df <- utils::read.delim(path, header = FALSE, sep = "\t", quote = "",
                          stringsAsFactors = FALSE, comment.char = "")
  if (ncol(df) < 2L) {
    stop("gene_functional_descriptions_wide.tsv: need >= 2 columns.", call. = FALSE)
  }
  data.frame(
    id = trimws(as.character(df[[1]])),
    description = trimws(as.character(df[[2]])),
    stringsAsFactors = FALSE
  )
}

#' Numeric helper for GenTable Fisher column ("0.015", "< 1e-30").
lgp_parse_topgo_fisher_col <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

#' Run topGO Fisher test + tables for one ontology.
lgp_run_topgo_ontology <- function(
    analysis_label,
    ontology,
    genes_of_interest,
    gene_id2go,
    annot_df,
    results_dir,
    dep_paths_for_skip
  ) {
  ontology <- toupper(trimws(as.character(ontology)[1L]))
  stopifnot(ontology %in% c("MF", "BP", "CC"))

  safe_lab <- gsub("[^A-Za-z0-9._-]+", "_", analysis_label)
  out_terms <- file.path(results_dir, paste0(safe_lab, "_", ontology, "_go_terms.csv"))
  out_genes <- file.path(results_dir, paste0(safe_lab, "_", ontology, "_go_term_genes.csv"))

  if (
    !lgp_should_rerun_external(out_terms, dep_paths_for_skip) &&
      !lgp_should_rerun_external(out_genes, dep_paths_for_skip)
  ) {
    message("[lgp-cache] Skipping topGO ", ontology, " for ", analysis_label, " (fresh outputs).")
    return(invisible(NULL))
  }

  gene_universe <- names(gene_id2go)
  goi <- intersect(unique(trimws(as.character(genes_of_interest))), gene_universe)

  if (length(goi) < go_min_genes) {
    warning(
      "[GO] ", analysis_label, " ", ontology, ": only ",
      length(goi),
      " gene(s) with GO annotations (minimum ",
      go_min_genes,
      "); skipping.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  gene_list <- factor(as.integer(gene_universe %in% goi))
  names(gene_list) <- gene_universe

  my_go_data <- methods::new(
    "topGOdata",
    description = paste(analysis_label, ontology),
    ontology = ontology,
    allGenes = gene_list,
    annot = topGO::annFUN.gene2GO,
    gene2GO = gene_id2go
  )

  result_test <- tryCatch(
    topGO::runTest(my_go_data, algorithm = "parentchild", statistic = "fisher"),
    error = function(e) {
      message("[GO] parentchild failed (", conditionMessage(e), "); using classic.")
      topGO::runTest(my_go_data, algorithm = "classic", statistic = "fisher")
    }
  )

  sc <- topGO::score(result_test)
  n_sig <- sum(is.finite(sc) & sc <= 0.05, na.rm = TRUE)
  top_n <- max(20L, min(n_sig + 10L, length(sc)))

  all_res <- topGO::GenTable(
    my_go_data,
    classicFisher = result_test,
    orderBy = "classicFisher",
    ranksOf = "classicFisher",
    topNodes = top_n
  )

  all_res$fisher_numeric <- lgp_parse_topgo_fisher_col(all_res$classicFisher)
  sig_terms_df <- all_res[is.finite(all_res$fisher_numeric) & all_res$fisher_numeric <= 0.05, ,
    drop = FALSE
  ]

  utils::write.csv(sig_terms_df, out_terms, row.names = FALSE)

  my_terms <- sig_terms_df$GO.ID
  if (!length(my_terms)) {
    message("[GO] ", analysis_label, " ", ontology, ": no terms at Fisher <= 0.05 (topNodes=", top_n, ").")
    utils::write.csv(data.frame(note = "no_significant_terms"), out_genes, row.names = FALSE)
    return(invisible(NULL))
  }

  my_genes <- topGO::genesInTerm(my_go_data, my_terms)
  rows <- list()
  for (term in my_terms) {
    gi <- unique(as.character(my_genes[[term]]))
    gi <- gi[gi %in% goi]
    if (!length(gi)) next
    rows[[length(rows) + 1L]] <- data.frame(GO.ID = term, gene_id = gi, stringsAsFactors = FALSE)
  }

  if (!length(rows)) {
    utils::write.csv(data.frame(note = "no_mapped_genes_in_terms"), out_genes, row.names = FALSE)
    return(invisible(NULL))
  }

  m <- do.call(rbind, rows)
  m <- merge(m, annot_df, by.x = "gene_id", by.y = "id", all.x = TRUE)
  utils::write.csv(m, out_genes, row.names = FALSE)
  message("[GO] ", analysis_label, " ", ontology, ": wrote ", nrow(sig_terms_df), " term(s); gene mapping -> ", basename(out_genes))
  invisible(NULL)
}

#' Read one topGO GenTable export; return NULL if empty or placeholder.
lgp_read_go_terms_result_csv <- function(path) {
  if (!isTRUE(file.exists(path))) {
    return(NULL)
  }
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(d)) {
    return(NULL)
  }
  if ("note" %in% names(d)) {
    return(NULL)
  }
  if (!all(c("GO.ID", "Term") %in% names(d))) {
    return(NULL)
  }
  if ("fisher_numeric" %in% names(d)) {
    d$p_value <- suppressWarnings(as.numeric(d$fisher_numeric))
  } else if ("classicFisher" %in% names(d)) {
    d$p_value <- lgp_parse_topgo_fisher_col(d$classicFisher)
  } else {
    return(NULL)
  }
  d <- d[is.finite(d$p_value) & d$p_value > 0 & d$p_value <= 1, , drop = FALSE]
  if (!nrow(d)) {
    return(NULL)
  }
  data.frame(
    GO.ID = trimws(as.character(d$GO.ID)),
    Term = trimws(as.character(d$Term)),
    p_value = d$p_value,
    stringsAsFactors = FALSE
  )
}

#' Heatmap of -log10(Fisher p) for enriched terms (rows) vs assay x ontology (columns).
lgp_plot_go_enrichment_heatmap <- function(
    results_dir,
    go_plot_dir,
    ontologies,
    analyses = c("pRDA", "LFMM"),
    max_terms = NULL,
    nlp_cap = NULL
  ) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning(
      "[GO] Skipping heatmap: install ggplot2: install.packages(\"ggplot2\")",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  mt <- suppressWarnings(as.integer(getOption("lgp.go_heatmap_max_terms", 80L))[1L])
  if (!is.finite(mt) || mt < 5L) {
    mt <- 80L
  }
  max_terms_use <- if (!is.null(max_terms)) max_terms else mt

  cap <- suppressWarnings(as.numeric(getOption("lgp.go_heatmap_neglog10_cap", 12))[1L])
  if (!is.finite(cap) || cap < 3) {
    cap <- 12
  }
  nlp_cap_use <- if (!is.null(nlp_cap)) nlp_cap else cap

  long_lst <- list()
  facet_order <- character(0)

  onts <- toupper(trimws(unique(ontologies[ontologies %in% c("MF", "BP", "CC")])))

  for (a in analyses) {
    safe_lab <- gsub("[^A-Za-z0-9._-]+", "_", a)
    for (o in onts) {
      path <- file.path(results_dir, paste0(safe_lab, "_", o, "_go_terms.csv"))
      blk <- lgp_read_go_terms_result_csv(path)
      if (is.null(blk)) next
      facet <- paste(a, o, sep = "_")
      facet_order <- unique(c(facet_order, facet))
      blk$facet <- facet
      long_lst[[length(long_lst) + 1L]] <- blk
    }
  }

  if (!length(long_lst)) {
    message("[GO] No significant GO term tables found for heatmap (expected *_*_go_terms.csv).")
    return(invisible(NULL))
  }

  long <- do.call(rbind, long_lst)
  rownames(long) <- NULL

  uniq_go <- unique(long$GO.ID)
  ncol_m <- length(facet_order)
  mat <- matrix(NA_real_,
    nrow = length(uniq_go),
    ncol = ncol_m,
    dimnames = list(uniq_go, facet_order))

  long$nlp_uncapped <- -log10(pmax(long$p_value, .Machine$double.xmin))
  long$nlp <- pmin(long$nlp_uncapped, nlp_cap_use)

  for (k in seq_len(nrow(long))) {
    gi <- long$GO.ID[k]
    fj <- long$facet[k]
    v <- long$nlp[k]
    ov <- mat[gi, fj]
    if (!is.finite(ov)) {
      mat[gi, fj] <- v
    } else if (is.finite(v)) {
      mat[gi, fj] <- max(ov, v, na.rm = TRUE)
    }
  }

  dup_term <- aggregate(Term ~ GO.ID, data = long, function(x) as.character(utils::head(x, 1)))
  colnames(dup_term) <- c("GO.ID", "Term")

  sc <- apply(mat, 1L, function(z) suppressWarnings(max(z, na.rm = TRUE)))
  sc[!is.finite(sc)] <- 0

  nk <- names(sort.int(sc, decreasing = TRUE))[seq_len(min(max_terms_use, length(sc)))]
  nk <- nk[sc[nk] > 0 & is.finite(sc[nk])]
  if (!length(nk)) {
    message("[GO] Heatmap skipped: no terms with finite enrichment scores.")
    return(invisible(NULL))
  }

  nk <- nk[seq_len(min(length(nk), max_terms_use))]
  mat_sub <- mat[nk, , drop = FALSE]

  nr <- nrow(mat_sub)
  nc <- ncol(mat_sub)
  plot_df <- data.frame(
    GO.ID = rep(rownames(mat_sub), times = nc),
    facet = rep(colnames(mat_sub), each = nr),
    nlp = as.vector(mat_sub),
    stringsAsFactors = FALSE
  )
  plot_df <- merge(plot_df, dup_term, by = "GO.ID", sort = FALSE, all.x = TRUE)

  go_levels <- rownames(mat_sub)[order(apply(mat_sub, 1L, max, na.rm = TRUE))]
  plot_df$GO.fac <- factor(plot_df$GO.ID, levels = go_levels)

  id2lbl <- dup_term[!duplicated(dup_term$GO.ID), , drop = FALSE]

  tn <- trimws(as.character(id2lbl$Term))
  id_lab <- ifelse(
    nchar(tn) > 58L,
    paste0(trimws(substr(tn, 1L, 55L)), "..."),
    tn
  )
  names(id_lab) <- trimws(as.character(id2lbl$GO.ID))

  plot_df$facet <- factor(plot_df$facet, levels = facet_order)
  lvl <- levels(plot_df$GO.fac)

  ylab_txt <- vapply(lvl, function(g) {
    lab <- suppressWarnings(trimws(as.character(id_lab[g][1])))
    if (length(lab) != 1L || is.na(lab) || identical(lab, "NA") || !nzchar(lab)) {
      paste0(g, " (no Term)")
    } else {
      lab
    }
  }, FUN.VALUE = character(1L))
  ylab_map <- structure(as.character(ylab_txt), names = as.character(lvl))

  heat_png <- file.path(go_plot_dir, "go_enrichment_heatmap_minusLog10_Fisher.png")
  heat_csv <- file.path(go_plot_dir, "go_enrichment_heatmap_matrix.csv")

  deps_heat <- c(
    Sys.glob(file.path(results_dir, "*_*_go_terms.csv"))
  )

  util_write <- TRUE
  if (length(deps_heat) && isTRUE(file.exists(heat_png)) && !lgp_should_rerun_external(heat_png, deps_heat)) {
    message("[GO] Heatmap PNG up to date: ", basename(heat_png))
    util_write <- FALSE
  }

  if (util_write) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(
      x = .data[["facet"]],
      y = .data[["GO.fac"]],
      fill = .data[["nlp"]]
    )) +
      ggplot2::geom_tile(color = "#e8e8e8") +
      ggplot2::scale_fill_gradient(
        limits = c(0, nlp_cap_use),
        na.value = "#f7f7f7",
        low = "#fffff0",
        high = "#084594",
        name = "-log10 Fisher P\n(capped fill)"
      ) +
      ggplot2::scale_y_discrete(
        breaks = lvl,
        labels = ylab_map[lvl],
      ) +
      ggplot2::labs(
        title = "Enriched GO terms (Fisher)",
        subtitle = sprintf(
          paste0(
            "Top %s terms by strongest column; fill = min(-log10(p), %s) per assay_ontology;",
            "\ncolumns: pRDA/LFMM x MF/BP/CC (only assays with enrichment output)."
          ),
          nrow(mat_sub),
          format(signif(nlp_cap_use, 3))
        ),
        x = NULL,
        y = "GO term"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(size = 7)
      )

    h_in <- min(52, max(5, nr * 0.18))
    ggplot2::ggsave(heat_png, p, width = 9, height = h_in, dpi = 220, limitsize = FALSE, bg = "white")
    utils::write.csv(cbind(mat_sub, Term = ylab_map[row.names(mat_sub)]), heat_csv, row.names = TRUE)
    message("[GO] Wrote GO heatmap -> ", basename(heat_png))
  }

  invisible(heat_png)
}

# ---- Load annotation backbone -------------------------------------------------

if (!file.exists(go_terms_wide)) {
  stop("Missing gene GO mapping file:\n  ", go_terms_wide, call. = FALSE)
}
if (!file.exists(go_desc_wide)) {
  stop("Missing gene descriptions file:\n  ", go_desc_wide, call. = FALSE)
}

gene_id2go <- lgp_read_gene2go_wide(go_terms_wide)
annot_tbl <- lgp_read_func_desc_wide(go_desc_wide)

deps_go <- c(go_terms_wide, go_desc_wide)

# ---- pRDA-linked genes ---------------------------------------------------------

genes_prda <- lgp_linked_genes_bedtools_gene_ids(linked_prda)
if (!length(genes_prda)) {
  warning(
    "No genes parsed from ", linked_prda,
    "\nRun Scripts/05-LinkageDecay.R (bedtools intersect) first.",
    call. = FALSE
  )
} else {
  utils::write.csv(
    data.frame(gene_id = genes_prda, stringsAsFactors = FALSE),
    file.path(results_dir, "pRDA_linked_unique_genes.csv"),
    row.names = FALSE
  )
  message("[GO] pRDA linked genes (unique, annotated genome): ", length(genes_prda))
  for (ont in go_ontologies) {
    lgp_run_topgo_ontology(
      "pRDA",
      ont,
      genes_prda,
      gene_id2go,
      annot_tbl,
      results_dir,
      c(deps_go, linked_prda)
    )
  }
}

# ---- LFMM-linked genes ---------------------------------------------------------

genes_lfmm <- lgp_linked_genes_bedtools_gene_ids(linked_lfmm)
if (!length(genes_lfmm)) {
  warning(
    "No genes parsed from ", linked_lfmm,
    "\nRun Scripts/05-LinkageDecay.R (bedtools intersect) first.",
    call. = FALSE
  )
} else {
  utils::write.csv(
    data.frame(gene_id = genes_lfmm, stringsAsFactors = FALSE),
    file.path(results_dir, "LFMM_linked_unique_genes.csv"),
    row.names = FALSE
  )
  message("[GO] LFMM linked genes (unique, annotated genome): ", length(genes_lfmm))
  for (ont in go_ontologies) {
    lgp_run_topgo_ontology(
      "LFMM",
      ont,
      genes_lfmm,
      gene_id2go,
      annot_tbl,
      results_dir,
      c(deps_go, linked_lfmm)
    )
  }
}

lgp_plot_go_enrichment_heatmap(
  results_dir = results_dir,
  go_plot_dir = go_plot_dir,
  ontologies = go_ontologies,
  analyses = c("pRDA", "LFMM")
)

message("[GO] Results directory: ", results_dir)
