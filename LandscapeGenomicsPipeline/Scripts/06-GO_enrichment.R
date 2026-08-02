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
# Heatmaps (optional; one PNG per ontology MF/BP/CC): ggplot2
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

#' Official GO term names from GO.db (repairs topGO GenTable truncation).
lgp_go_term_names <- function(go_ids, fallback = NULL) {
  go_ids <- trimws(as.character(go_ids))
  n <- length(go_ids)
  out <- rep(NA_character_, n)
  if (!is.null(fallback) && length(fallback) == n) {
    out <- trimws(as.character(fallback))
  }
  if (!n) {
    return(out)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("GO.db", quietly = TRUE)) {
    return(out)
  }
  all_terms <- tryCatch(
    AnnotationDbi::Term(GO.db::GOTERM),
    error = function(e) NULL
  )
  if (is.null(all_terms) || !length(all_terms)) {
    return(out)
  }
  looked <- unname(all_terms[go_ids])
  ok <- !is.na(looked) & nzchar(looked)
  out[ok] <- looked[ok]
  out
}

#' Wrap labels for ggplot axes without truncating text.
lgp_wrap_axis_label <- function(x, width = 48L) {
  width <- as.integer(width)[1L]
  if (!is.finite(width) || width < 20L) width <- 48L
  x <- as.character(x)
  vapply(x, function(s) {
    if (length(s) != 1L || is.na(s) || !nzchar(s)) {
      return(as.character(s))
    }
    paste(strwrap(s, width = width), collapse = "\n")
  }, FUN.VALUE = character(1L), USE.NAMES = FALSE)
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
  if (nrow(sig_terms_df) && "Term" %in% names(sig_terms_df) && "GO.ID" %in% names(sig_terms_df)) {
    sig_terms_df$Term <- lgp_go_term_names(sig_terms_df$GO.ID, fallback = sig_terms_df$Term)
  }

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
    Term = lgp_go_term_names(d$GO.ID, fallback = d$Term),
    p_value = d$p_value,
    stringsAsFactors = FALSE
  )
}

#' Flag GO terms plausibly related to climate adaptation in a temperate tree
#' (e.g. Asimina triloba): water/drought, ABA, stomata, light/photosynthesis,
#' temperature-related lipids, oxidative/abiotic stress, root foraging, etc.
#'
#' Uses keyword matches on Term plus a curated GO.ID allowlist (helps when
#' topGO truncates long Term strings). Toggle with
#' options(lgp.go_heatmap_climate_filter = FALSE) to disable.
lgp_go_climate_adaptation_relevant <- function(go_id, term) {
  go_id <- toupper(trimws(as.character(go_id)))
  term <- trimws(as.character(term))
  n <- max(length(go_id), length(term))
  if (!n) {
    return(logical(0))
  }
  go_id <- rep_len(go_id, n)
  term <- rep_len(term, n)

  # Explicit IDs kept even if Term is truncated in GenTable exports.
  allow_ids <- c(
    "GO:1902265", # abscisic acid homeostasis
    "GO:0042631", # cellular response to water deprivation
    "GO:0071462", # cellular response to water stimulus
    "GO:0009414", # response to water deprivation
    "GO:0009415", # response to water
    "GO:0104004", # cellular response to environmental stimulus
    "GO:0071496", # cellular response to external stimulus
    "GO:0009581", # detection of external stimulus
    "GO:0009720", # detection of hormone stimulus
    "GO:0010375", # stomatal complex patterning
    "GO:0010376", # stomatal complex formation
    "GO:0010440", # stomatal lineage progression
    "GO:2000037", # regulation of stomatal complex patterning
    "GO:0010444", # guard mother cell differentiation
    "GO:0009640", # photomorphogenesis
    "GO:0009639", # response to red or far red light
    "GO:2000030", # regulation of response to red or far red light
    "GO:0042548", # regulation of photosynthesis, light reaction
    "GO:0009723", # response to ethylene
    "GO:0010540", # basipetal auxin transport
    "GO:0009735", # response to cytokinin
    "GO:0009884", # cytokinin receptor activity
    "GO:0098754", # detoxification
    "GO:0009636", # response to toxic substance
    "GO:0046686", # response to cadmium ion
    "GO:0009635", # response to herbicide
    "GO:0072756", # cellular response to paraquat (ROS)
    "GO:0071731", # response to nitric oxide
    "GO:0071732", # cellular response to nitric oxide
    "GO:0006749", # glutathione metabolic process
    "GO:0043295", # glutathione binding
    "GO:0009698", # phenylpropanoid metabolic process
    "GO:0010023", # proanthocyanidin biosynthetic process
    "GO:0030968", # ER unfolded protein response
    "GO:0031990", # mRNA export in response to heat stress (often truncated)
    "GO:0010086", # embryonic root morphogenesis
    "GO:0080022", # primary root development
    "GO:0008643", # carbohydrate transport
    "GO:0034219", # carbohydrate transmembrane transport
    "GO:0015144", # carbohydrate transmembrane transporter activity
    "GO:1901569", # fatty acid derivative catabolic process
    "GO:0036115", # fatty-acyl-CoA catabolic process
    "GO:0046459", # short-chain fatty acid metabolic process
    "GO:0005452", # solute:inorganic anion antiporter activity
    "GO:0006855", # xenobiotic transmembrane transport
    "GO:0042908"  # xenobiotic transport
  )

  pos <- paste0(
    "(?i)",
    paste(
      c(
        "water", "drought", "dehydrat", "desiccat",
        "abscisic", "\\bABA\\b",
        "osmotic", "osmolyte", "salinity", "salt stress",
        "cold", "chill", "freez", "\\bheat\\b", "temperature", "therm[ao]",
        "\\bstress\\b",
        "stomatal", "\\bstomata\\b", "guard (mother )?cell",
        "photosynth", "light reaction", "far[- ]?red", "photomorph",
        "ethylene", "jasmon", "salicylic",
        "antioxid", "oxidative", "reactive oxygen", "glutathione", "detox",
        "phenylpropanoid", "flavonoid", "proanthocyanidin", "anthocyanin",
        "cuticle", "\\bwax\\b", "suberin",
        "environmental stim", "external stimulus", "\\babiotic\\b",
        "cadmium", "toxic substance", "xenobiotic", "heavy metal", "paraquat", "herbicide",
        "unfolded protein", "heat shock",
        "nitric oxide",
        "basipetal auxin",
        "root morphogenesis", "primary root", "lateral root",
        "carbohydrate transport", "carbohydrate transmembrane",
        "fatty acid", "fatty-acyl",
        "anion antiporter"
      ),
      collapse = "|"
    )
  )
  neg <- paste0(
    "(?i)",
    paste(
      c(
        "cartilage", "muscle cell", "limbic", "orbitofrontal",
        "bronchodilator", "cytokine production", "heparin",
        "glycosaminoglycan", "connective tissue",
        "\\bX-ray\\b", "bleomycin", "vitamin B1", "monosaccharide"
      ),
      collapse = "|"
    )
  )

  hit_id <- go_id %in% allow_ids
  hit_kw <- grepl(pos, term, perl = TRUE)
  hit_neg <- grepl(neg, term, perl = TRUE)
  (hit_id | hit_kw) & !hit_neg
}

#' Heatmaps of -log10(Fisher p): one plot per ontology (MF / BP / CC),
#' with pRDA and LFMM as columns within each plot.
#' Color scale is shared across ontologies so panels are comparable.
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

  onts <- toupper(trimws(unique(ontologies[ontologies %in% c("MF", "BP", "CC")])))
  ont_labels <- c(
    MF = "Molecular Function (MF)",
    BP = "Biological Process (BP)",
    CC = "Cellular Component (CC)"
  )

  assay_order <- unique(as.character(analyses))
  panels <- list()
  obs_all <- numeric(0)
  climate_filter <- {
    opt <- getOption("lgp.go_heatmap_climate_filter", TRUE)
    if (is.logical(opt)) isTRUE(opt[1L]) else !identical(tolower(trimws(as.character(opt[1L]))), "false")
  }
  kept_terms_rows <- list()

  for (o in onts) {
    long_lst <- list()
    for (a in assay_order) {
      safe_lab <- gsub("[^A-Za-z0-9._-]+", "_", a)
      path <- file.path(results_dir, paste0(safe_lab, "_", o, "_go_terms.csv"))
      blk <- lgp_read_go_terms_result_csv(path)
      if (is.null(blk)) next
      blk$facet <- a
      long_lst[[length(long_lst) + 1L]] <- blk
    }

    if (!length(long_lst)) {
      message("[GO] No significant ", o, " term tables found for heatmap.")
      next
    }

    long <- do.call(rbind, long_lst)
    rownames(long) <- NULL

    climate_filter_local <- climate_filter
    if (isTRUE(climate_filter_local)) {
      keep_cli <- lgp_go_climate_adaptation_relevant(long$GO.ID, long$Term)
      n_before <- length(unique(long$GO.ID))
      long <- long[keep_cli, , drop = FALSE]
      message(
        "[GO] Climate-adaptation filter (", o, "): kept ",
        length(unique(long$GO.ID)), " / ", n_before, " enriched term(s)."
      )
      if (!nrow(long)) {
        message("[GO] No climate-relevant ", o, " terms after filter; skipping heatmap.")
        next
      }
    }

    facet_order <- assay_order[assay_order %in% unique(long$facet)]
    uniq_go <- unique(long$GO.ID)
    mat <- matrix(NA_real_,
      nrow = length(uniq_go),
      ncol = length(facet_order),
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
      message("[GO] Heatmap skipped for ", o, ": no terms with finite enrichment scores.")
      next
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
    id2lbl$Term <- lgp_go_term_names(id2lbl$GO.ID, fallback = id2lbl$Term)

    # Full names for CSVs / inventory; wrapped copy only for axis display.
    tn_full <- trimws(as.character(id2lbl$Term))
    names(tn_full) <- trimws(as.character(id2lbl$GO.ID))
    wrap_w <- suppressWarnings(as.integer(getOption("lgp.go_heatmap_label_wrap", 48L))[1L])
    if (!is.finite(wrap_w) || wrap_w < 20L) wrap_w <- 48L
    tn_wrap <- lgp_wrap_axis_label(tn_full, width = wrap_w)
    names(tn_wrap) <- names(tn_full)

    plot_df$facet <- factor(plot_df$facet, levels = facet_order)
    lvl <- levels(plot_df$GO.fac)

    ylab_txt <- vapply(lvl, function(g) {
      lab <- suppressWarnings(trimws(as.character(tn_full[g][1])))
      if (length(lab) != 1L || is.na(lab) || identical(lab, "NA") || !nzchar(lab)) {
        paste0(g, " (no Term)")
      } else {
        lab
      }
    }, FUN.VALUE = character(1L))
    ylab_map <- structure(as.character(ylab_txt), names = as.character(lvl))

    ylab_display <- vapply(lvl, function(g) {
      lab <- suppressWarnings(as.character(tn_wrap[g][1]))
      if (length(lab) != 1L || is.na(lab) || identical(lab, "NA") || !nzchar(lab)) {
        ylab_map[[g]]
      } else {
        lab
      }
    }, FUN.VALUE = character(1L))
    ylab_display <- structure(as.character(ylab_display), names = as.character(lvl))

    heat_png <- file.path(go_plot_dir, paste0("go_enrichment_heatmap_", o, "_minusLog10_Fisher.png"))
    heat_csv <- file.path(go_plot_dir, paste0("go_enrichment_heatmap_", o, "_matrix.csv"))

    obs <- as.numeric(mat_sub)
    obs <- obs[is.finite(obs)]
    obs_all <- c(obs_all, obs)

    kept_terms_rows[[length(kept_terms_rows) + 1L]] <- data.frame(
      ontology = o,
      GO.ID = rownames(mat_sub),
      Term = unname(ylab_map[rownames(mat_sub)]),
      stringsAsFactors = FALSE
    )

    panels[[o]] <- list(
      ontology = o,
      mat_sub = mat_sub,
      plot_df = plot_df,
      ylab_map = ylab_map,
      ylab_display = ylab_display,
      lvl = lvl,
      nr = nr,
      nc = nc,
      heat_png = heat_png,
      heat_csv = heat_csv
    )
  }

  if (!length(panels)) {
    message("[GO] No significant GO term tables found for heatmap (expected *_*_go_terms.csv).")
    return(invisible(NULL))
  }

  if (length(kept_terms_rows)) {
    kept_df <- do.call(rbind, kept_terms_rows)
    rownames(kept_df) <- NULL
    kept_csv <- file.path(
      go_plot_dir,
      if (isTRUE(climate_filter)) {
        "go_climate_adaptation_terms_in_heatmaps.csv"
      } else {
        "go_terms_in_heatmaps.csv"
      }
    )
    utils::write.csv(kept_df, kept_csv, row.names = FALSE)
    message("[GO] Wrote term inventory -> ", basename(kept_csv), " (", nrow(kept_df), " row(s)).")
  }

  # Shared fill scale across MF/BP/CC so panels are visually comparable.
  scale_mode <- tolower(trimws(as.character(
    getOption("lgp.go_heatmap_scale_mode", "data")[1L]
  )))
  if (identical(scale_mode, "fixed") || !length(obs_all)) {
    fill_lo <- 0
    fill_hi <- nlp_cap_use
  } else {
    fill_lo <- min(obs_all)
    fill_hi <- max(obs_all)
    pad <- max(0.05, (fill_hi - fill_lo) * 0.08)
    fill_lo <- max(0, fill_lo - pad)
    fill_hi <- min(nlp_cap_use, fill_hi + pad)
    if (!is.finite(fill_hi) || fill_hi <= fill_lo) {
      fill_lo <- 0
      fill_hi <- nlp_cap_use
    }
  }
  fill_brks <- pretty(c(fill_lo, fill_hi), n = 5)
  fill_brks <- fill_brks[fill_brks >= fill_lo - 1e-9 & fill_brks <= fill_hi + 1e-9]
  if (!length(fill_brks)) {
    fill_brks <- c(fill_lo, fill_hi)
  }

  n_cols <- 6L
  pal_cols <- c(
    "#f7f4ef",
    "#ffe08a",
    "#ff9f4a",
    "#ef5d5d",
    "#b83280",
    "#3d1a5c"
  )
  if (length(obs_all) >= 2L && fill_hi > fill_lo) {
    q_stops <- as.numeric(stats::quantile(
      obs_all,
      probs = seq(0, 1, length.out = n_cols),
      names = FALSE,
      type = 7
    ))
    fill_vals <- (q_stops - fill_lo) / (fill_hi - fill_lo)
    fill_vals <- pmax(0, pmin(1, fill_vals))
    for (i in seq_along(fill_vals)[-1L]) {
      if (fill_vals[i] <= fill_vals[i - 1L]) {
        fill_vals[i] <- min(1, fill_vals[i - 1L] + 1e-4)
      }
    }
    fill_vals[1L] <- 0
    fill_vals[length(fill_vals)] <- 1
  } else {
    fill_vals <- seq(0, 1, length.out = n_cols)
  }

  message(
    "[GO] Shared heatmap fill scale across ",
    paste(names(panels), collapse = "/"),
    ": ",
    format(signif(fill_lo, 3)),
    "–",
    format(signif(fill_hi, 3))
  )

  # Always refresh per-ontology matrices; PNG is a stacked BP (top) + MF (bottom) figure.
  for (panel in panels) {
    utils::write.csv(
      cbind(panel$mat_sub, Term = panel$ylab_map[row.names(panel$mat_sub)]),
      panel$heat_csv,
      row.names = TRUE
    )
  }

  stack_order <- c("BP", "MF")
  stack_order <- stack_order[stack_order %in% names(panels)]
  # Any remaining ontologies (e.g. CC) get their own single-panel PNG.
  other_onts <- setdiff(names(panels), stack_order)

  combo_png <- file.path(
    go_plot_dir,
    if (length(stack_order) >= 2L) {
      paste0("go_enrichment_heatmap_", paste(stack_order, collapse = "_"), "_minusLog10_Fisher.png")
    } else if (length(stack_order) == 1L) {
      panels[[stack_order]]$heat_png
    } else {
      NA_character_
    }
  )
  other_pngs <- vapply(other_onts, function(o) panels[[o]]$heat_png, character(1L))
  out_pngs <- c(if (is.character(combo_png) && !is.na(combo_png)) combo_png, other_pngs)

  deps_heat <- Sys.glob(file.path(results_dir, "*_*_go_terms.csv"))
  need_write <- !length(deps_heat) || !length(out_pngs) || any(vapply(out_pngs, function(p) {
    !isTRUE(file.exists(p)) || isTRUE(lgp_should_rerun_external(p, deps_heat))
  }, logical(1L)))

  written <- character(0)
  if (!need_write) {
    for (pn in out_pngs) {
      message("[GO] Heatmap PNG up to date: ", basename(pn))
      written <- c(written, pn)
    }
    return(invisible(written))
  }

  lgp_draw_go_heatmap_stacked <- function(plot_df, y_levels, y_labels, title, subtitle) {
    ggplot2::ggplot(plot_df, ggplot2::aes(
      x = .data[["facet"]],
      y = .data[["ykey"]],
      fill = .data[["nlp"]]
    )) +
      ggplot2::geom_tile(color = "#f4f1ec", linewidth = 0.25) +
      ggplot2::scale_fill_gradientn(
        colours = pal_cols,
        values = fill_vals,
        limits = c(fill_lo, fill_hi),
        breaks = fill_brks,
        na.value = "#f0eee9",
        name = "-log10 Fisher P"
      ) +
      ggplot2::scale_y_discrete(breaks = y_levels, labels = y_labels) +
      ggplot2::facet_grid(ontology ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 9),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(size = 7, lineheight = 0.95),
        plot.margin = ggplot2::margin(8, 12, 8, 8),
        strip.placement = "outside",
        strip.text.y.left = ggplot2::element_text(angle = 90, face = "bold", size = 13),
        strip.background.y = ggplot2::element_rect(fill = "#efeae3", colour = NA),
        panel.spacing.y = grid::unit(0.6, "lines")
      )
  }

  if (length(stack_order)) {
    stack_rows <- lapply(stack_order, function(o) {
      pn <- panels[[o]]
      df <- pn$plot_df
      df$ontology <- o
      df$ykey <- paste(o, as.character(df$GO.ID), sep = "||")
      df
    })
    stack_df <- do.call(rbind, stack_rows)
    rownames(stack_df) <- NULL
    stack_df$ontology <- factor(stack_df$ontology, levels = stack_order)

    y_levels <- unlist(lapply(stack_order, function(o) {
      paste(o, panels[[o]]$lvl, sep = "||")
    }), use.names = FALSE)
    y_labels <- unlist(lapply(stack_order, function(o) {
      pn <- panels[[o]]
      lab <- pn$ylab_display
      if (is.null(lab)) lab <- pn$ylab_map
      unname(lab[pn$lvl])
    }), use.names = FALSE)
    names(y_labels) <- y_levels
    stack_df$ykey <- factor(stack_df$ykey, levels = y_levels)

    n_terms <- vapply(stack_order, function(o) panels[[o]]$nr, integer(1L))
    sub_txt <- paste0(
      if (length(stack_order) >= 2L) {
        "BP (top) and MF (bottom); "
      } else {
        paste0(stack_order[[1L]], "; ")
      },
      if (isTRUE(climate_filter)) {
        "climate-adaptation filter; "
      } else {
        ""
      },
      sprintf("shared scale %.2f–%.2f; columns: pRDA vs LFMM.", fill_lo, fill_hi)
    )

    p_stack <- lgp_draw_go_heatmap_stacked(
      stack_df,
      y_levels,
      y_labels,
      title = "Enriched GO terms (Fisher)",
      subtitle = sub_txt
    )

    n_label_lines <- sum(vapply(y_labels, function(s) {
      length(strsplit(as.character(s), "\n", fixed = TRUE)[[1L]])
    }, integer(1L)))
    h_in <- min(72, max(6, n_label_lines * 0.22 + 2.2))
    nc <- max(vapply(stack_order, function(o) panels[[o]]$nc, integer(1L)))
    w_in <- max(8.5, min(12, 6.8 + nc * 1.2))
    ggplot2::ggsave(combo_png, p_stack, width = w_in, height = h_in, dpi = 220, limitsize = FALSE, bg = "white")
    message("[GO] Wrote stacked GO heatmap -> ", basename(combo_png),
            " [", paste(stack_order, collapse = "+"), "]")
    written <- c(written, combo_png)

    # Remove obsolete single-ontology PNGs for stacked panels.
    for (o in stack_order) {
      old_png <- panels[[o]]$heat_png
      if (isTRUE(file.exists(old_png)) && !identical(normalizePath(old_png), normalizePath(combo_png))) {
        unlink(old_png)
      }
    }
  }

  for (o in other_onts) {
    panel <- panels[[o]]
    ont_title <- if (!is.na(ont_labels[o])) ont_labels[o] else o
    df <- panel$plot_df
    df$ykey <- factor(as.character(df$GO.ID), levels = panel$lvl)
    y_labels <- panel$ylab_display
    if (is.null(y_labels)) y_labels <- panel$ylab_map
    y_labels <- unname(y_labels[panel$lvl])
    names(y_labels) <- panel$lvl

    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data[["facet"]],
      y = .data[["ykey"]],
      fill = .data[["nlp"]]
    )) +
      ggplot2::geom_tile(color = "#f4f1ec", linewidth = 0.25) +
      ggplot2::scale_fill_gradientn(
        colours = pal_cols,
        values = fill_vals,
        limits = c(fill_lo, fill_hi),
        breaks = fill_brks,
        na.value = "#f0eee9",
        name = "-log10 Fisher P"
      ) +
      ggplot2::scale_y_discrete(breaks = panel$lvl, labels = y_labels) +
      ggplot2::labs(
        title = paste0("Enriched GO terms — ", ont_title),
        subtitle = sprintf(
          "Shared color scale (%.2f–%.2f). Columns: pRDA vs LFMM.",
          fill_lo, fill_hi
        ),
        x = NULL,
        y = "GO term"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(size = 7, lineheight = 0.95)
      )

    n_label_lines <- sum(vapply(y_labels, function(s) {
      length(strsplit(as.character(s), "\n", fixed = TRUE)[[1L]])
    }, integer(1L)))
    h_in <- min(60, max(5, n_label_lines * 0.22 + 1.8))
    w_in <- max(8, min(12, 6.5 + panel$nc * 1.2))
    ggplot2::ggsave(panel$heat_png, p, width = w_in, height = h_in, dpi = 220, limitsize = FALSE, bg = "white")
    message("[GO] Wrote GO heatmap (", o, ") -> ", basename(panel$heat_png))
    written <- c(written, panel$heat_png)
  }

  invisible(written)
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
