# LandscapeGenomicsPipeline - shared incremental cache helpers
#
# Behaviour:
#   - By default recomputation is skipped when a cache RDS exists and none of its
#     dependency paths are strictly newer than the cache file's mtime.
#   - Set options(lgp.use_cache = FALSE) to disable reading/writing RDS caches only
#     (PNG/CSV skips are separate flags where used).
#   - Set options(lgp.force_rerun = TRUE) OR env LGP_FORCE_RERUN=1 to rebuild caches.
#
# Recommended project root discovery (handled in helpers):
#   - options(lgp.project_root = "...") OR Sys.getenv("LGP_PROJECT_ROOT", unset = "~/Desktop/LandscapeGenomicsPipeline")

.lgp_cache_helpers_version <- "1.1"

## ---- Option defaults -----------------------------------------------------

.init_lgp_options <- function() {
  oo <- options()
  defs <- list(
    lgp.use_cache = TRUE,
    lgp.force_rerun = FALSE,
    ## Set FALSE to skip "output already newer than deps" shortcuts (e.g. PLINK skips)
    lgp.smart_skip_external = TRUE,
    ## Prefer "Outputs" (matches Scripts/00-Preprocessing.R); fall back if only "outputs" exists
    lgp.outputs_dir_name = NULL
  )
  to_set <- defs[vapply(names(defs), function(nm) is.null(oo[[nm]]), NA)]
  if (length(to_set)) do.call(options, to_set)
  invisible(TRUE)
}

## ---- Guards --------------------------------------------------------------

lgp_use_cache <- function() {
  isTRUE(getOption("lgp.use_cache", TRUE)) &&
    !lgp_force_rerun() &&
    !identical(Sys.getenv("LGP_FORCE_RERUN", ""), "1") &&
    !identical(Sys.getenv("LGP_FORCE_RERUN", ""), "true")
}

lgp_force_rerun <- function() {
  isTRUE(getOption("lgp.force_rerun", FALSE)) ||
    identical(Sys.getenv("LGP_FORCE_RERUN", ""), "1") ||
    identical(Sys.getenv("LGP_FORCE_RERUN", ""), "true")
}

lgp_smart_skip_external <- function() {
  isTRUE(getOption("lgp.smart_skip_external", TRUE))
}

## ---- Paths ---------------------------------------------------------------

#' Default pipeline root directory (normalized, no trailing slash).
lgp_project_root <- function() {
  pr <- Sys.getenv("LGP_PROJECT_ROOT", unset = "")
  if (!nzchar(pr)) pr <- getOption("lgp.project_root", "~/Desktop/LandscapeGenomicsPipeline")
  normalizePath(sub("/+$", "", path.expand(pr)), winslash = "/", mustWork = FALSE)
}

#' Top-level Outputs directory (capital O preferred for consistency).
lgp_outputs_base <- function(root = NULL) {
  root <- root %||% lgp_project_root()
  dn <- getOption("lgp.outputs_dir_name", NULL)
  if (!is.null(dn) && nzchar(dn)) {
    p <- file.path(root, dn)
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  cand <- c(file.path(root, "Outputs"), file.path(root, "outputs"))
  for (p in cand) {
    if (dir.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  p <- cand[[1]]
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

#' Directory containing 00-preprocessing RDS artifacts - whichever exists first.
#' @param create_default If neither exists yet, creates `Outputs/00-preprocessing`.
lgp_preprocess_dir <- function(root = NULL, create_default = TRUE) {
  root <- root %||% lgp_project_root()
  cand <- unique(c(
    file.path(lgp_outputs_base(root), "00-preprocessing"),
    file.path(root, "Outputs", "00-preprocessing"),
    file.path(root, "outputs", "00-preprocessing")
  ))
  hit <- cand[vapply(cand, function(p) file.exists(file.path(p, "coord.rds")), logical(1))]
  if (length(hit)) hit[[1]] else if (create_default) {
    p <- cand[[1]]
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    p
  } else cand[[1]]
}

#' Convenience: path to subdirectory under Outputs (step folder).
#' @param step e.g. "04-gea", "05-linkage-decay"
lgp_outputs_step_dir <- function(step, root = NULL) {
  p <- file.path(lgp_outputs_base(root), step)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---- Staleness / atomic IO ---------------------------------------------

#' TRUE when the cache record is stale: missing cache, unreadable cache, listed
#' dependency missing, or any existing dependency strictly newer than the cache mtime.
lgp_deps_newer_than_cache <- function(dep_paths, cache_path) {
  if (!file.exists(cache_path)) return(TRUE)

  ct <- suppressWarnings(file.info(cache_path)$mtime[[1]])
  dep_paths <- dep_paths[vapply(dep_paths, function(x) nzchar(trimws(as.character(x))), logical(1))]
  dep_paths <- unique(dep_paths)
  if (!length(dep_paths)) return(FALSE)

  for (p in dep_paths) {
    if (!file.exists(p)) return(TRUE)
    mt <- suppressWarnings(file.info(p)$mtime[[1]])
    if (!is.na(mt) && mt > ct) return(TRUE)
  }
  FALSE
}

lgp_atomic_saveRDS <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(.basename_safe(path), ".lgp."), tmpdir = dirname(path), fileext = ".rds.partial")
  on.exit(try(unlink(tmp), silent = TRUE), add = TRUE)
  saveRDS(object, file = tmp, compress = compress)
  if (.Platform$OS.type == "windows") {
    if (file.exists(path)) file.remove(path)
  }
  if (!file.rename(tmp, path)) {
    ## Fallback (e.g., cross-device)
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

.basename_safe <- function(x) sub("[.][^.]+$", "", basename(x))

## ---- terra SpatRaster (saveRDS needs wrap/unwrap) -----------------------

.lgp_has_terra <- function() requireNamespace("terra", quietly = TRUE)

.lgp_is_spatraster_like <- function(x) {
  .lgp_has_terra() &&
    (inherits(x, "SpatRaster") || inherits(x, "PackedSpatRaster"))
}

.lgp_spatraster_usable <- function(x) {
  if (!.lgp_is_spatraster_like(x)) return(TRUE)
  tryCatch({
    suppressMessages(terra::nlyr(x))
    TRUE
  }, error = function(e) FALSE)
}

.lgp_terra_serialize <- function(x) {
  if (!.lgp_has_terra()) return(x)
  if (inherits(x, "SpatRaster")) return(terra::wrap(x))
  if (is.list(x)) {
    out <- lapply(x, .lgp_terra_serialize)
    if (!is.null(names(x))) names(out) <- names(x)
    return(out)
  }
  x
}

.lgp_terra_deserialize <- function(x) {
  if (!.lgp_has_terra()) return(x)
  if (inherits(x, "PackedSpatRaster")) return(terra::unwrap(x))
  if (is.list(x)) {
    out <- lapply(x, .lgp_terra_deserialize)
    if (!is.null(names(x))) names(out) <- names(x)
    return(out)
  }
  x
}

.lgp_terra_objects_valid <- function(x) {
  if (!.lgp_has_terra()) return(TRUE)
  if (.lgp_is_spatraster_like(x)) return(.lgp_spatraster_usable(x))
  if (is.list(x)) return(all(vapply(x, .lgp_terra_objects_valid, logical(1))))
  TRUE
}

#' Pack terra objects for plain saveRDS (unwrap on read with \code{lgp_read_rds}).
lgp_pack_for_rds <- function(x) .lgp_terra_serialize(x)

#' Restore values written by \code{lgp_pack_for_rds} / cache helpers.
lgp_read_rds <- function(path) {
  if (!file.exists(path)) {
    stop("RDS not found: ", path, call. = FALSE)
  }
  raw <- readRDS(path)
  val <- if (is.list(raw$.lgp) && !is.null(raw$value)) raw$value else raw
  val <- .lgp_terra_deserialize(val)
  if (!.lgp_terra_objects_valid(val)) {
    stop(
      "Invalid or stale terra raster in ", path,
      " (delete cache or set LGP_FORCE_RERUN=1).",
      call. = FALSE
    )
  }
  val
}

## ---- Main API ------------------------------------------------------------

#' Cached RDS: return `compute()` unless cache exists with fresh deps and matching meta.
#'
#' Wrapped value is saved as list(.lgp = list(version, meta), value = OBJ).
#'
#' @param path Cache file (.rds)
#' @param dep_paths Character paths to inputs (mtime compared to cache).
#' @param compute No-arg expression function or function() OBJ
#' @param meta Optional list; mismatched meta forces recomputation even if deps \"fresh\".
#' @param label Human-readable prefix for messages
#' @param compress Passed to saveRDS
lgp_rds_cached <- function(path,
                           dep_paths = character(),
                           compute,
                           meta = NULL,
                           label = basename(path),
                           compress = "xz") {
  if (!is.character(path) || !nzchar(path)) stop("`path` must be a non-empty string.", call. = FALSE)
  compute_fn <- compute
  if (!is.function(compute_fn)) {
    compute_fn <- function() compute
  }

  dep_paths <- unique(as.character(dep_paths))

  force_re <- lgp_force_rerun()

  if (!force_re && lgp_use_cache() && file.exists(path) && !lgp_deps_newer_than_cache(dep_paths, path)) {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(obj) && is.list(obj$.lgp) && !is.null(obj$value)) {
      lg_ok <- TRUE
      if (!is.null(meta)) {
        if (is.null(obj$.lgp$meta) || !identical(obj$.lgp$meta, meta)) lg_ok <- FALSE
      }
      if (lg_ok) {
        val <- .lgp_terra_deserialize(obj$value)
        if (.lgp_terra_objects_valid(val)) {
          message("[lgp-cache] Loading ", label, " <- ", path)
          return(val)
        }
        message("[lgp-cache] Stale terra raster; recomputing: ", label)
      }
    }
  }

  if (force_re) {
    message("[lgp-cache] Forced recompute: ", label)
  } else {
    message("[lgp-cache] Computing & saving ", label, " -> ", path)
  }

  val <- compute_fn()
  lg <- list(version = .lgp_cache_helpers_version, meta = meta)
  lgp_atomic_saveRDS(
    list(.lgp = lg, value = .lgp_terra_serialize(val)),
    path,
    compress = compress
  )
  val
}

#' Cache a richer object bundle (often a plain list saved without .lgp wrapper when you
#' control version yourself via meta).
#'
#' Convenience wrapper identical to \code{lgp_rds_cached()} but emphasizes list outputs.
lgp_list_cached <- function(path, dep_paths = character(), compute_list, meta = NULL,
                           label = basename(path)) {
  out <- lgp_rds_cached(
    path = path,
    dep_paths = dep_paths,
    compute = compute_list,
    meta = meta,
    label = label
  )
  if (!is.list(out)) warning(
    "`lgp_list_cached` compute did not return a list (got ",
    paste(class(out), collapse = ", "), ")",
    call. = FALSE
  )
  out
}

## ---- Assign list into caller environment ---------------------------------

#' Populate an environment from a named list (typical workflow: unwrap cached bundles).
#'
#' @param lst Named list
#' @param envir Assign here (defaults to `.GlobalEnv` for scripted pipelines).
lgp_assign_list <- function(lst, names = names(lst), envir = globalenv()) {
  if (!length(lst)) invisible(return(NULL))
  if (length(names)) {
    lst <- lst[names]
  }
  for (nm in names(lst)) assign(nm, lst[[nm]], envir = envir)
  invisible(nm)
}

## ---- Raster graphics skip (PNG / JPEG) -----------------------------------

#' Open a raster plot device unless the output exists and deps are newer.
#'
#' @param path Output path; extension selects device (`.png` default, `.jpg`/`.jpeg` -> `jpeg`).
#' Returns `TRUE` if a new graphic was rendered (caller should `force(plot_expr)` and `dev.off()`).
#'
#' Typical pattern:
#'
#' ```
#' plot_fn <- function() { ... plotting code ... }
#' if (lgp_begin_graphics_if_stale(path, dep_paths)) { plot_fn(); dev.off() }
#' ```
lgp_begin_graphics_if_stale <- function(path,
                                        plot_expr,
                                        dep_paths = character(),
                                        width = 1600,
                                        height = 1200,
                                        res = 200,
                                        ...) {
  rewrite <- TRUE
  if (lgp_use_cache() && file.exists(path)) {
    rewrite <- lgp_deps_newer_than_cache(dep_paths, path)
  }
  if (!rewrite) {
    message("[lgp-cache] Skip raster graphic (fresh): ", path)
    return(FALSE)
  }

  ext <- tools::file_ext(path)
  el <- if (tolower(ext) %in% c("jpeg", "jpg")) "jpeg" else "png"

  if (tolower(el) == "jpeg") {
    jargs <- list(filename = path, width = width, height = height, ...)
    if (!is.na(res)) jargs$res <- res
    do.call(grDevices::jpeg, jargs)
  } else {
    pn_args <- list(filename = path, width = width, height = height, ...)
    if (!is.na(res)) pn_args$res <- res
    do.call(grDevices::png, pn_args)
  }
  force(plot_expr)
  grDevices::dev.off()
  TRUE
}

#' PNG shorthand (same as [`lgp_begin_graphics_if_stale`]).
lgp_png_if_stale <- function(path, plot_expr, dep_paths = character(),
                             width = 1600, height = 1200, res = 200, ...) {
  lgp_begin_graphics_if_stale(
    path, plot_expr, dep_paths,
    width = width, height = height, res = res, ...
  )
}

## ---- External command skip ------------------------------------------------

#' Decide whether PLINK/other external output can be reused.
#'
#' Uses mtime comparisons only (not argument fingerprinting unless caller bundles args in dep_paths).
#'
#' @param prod Primary output expected (must exist).
#' @param dep_paths Dependencies that should not be newer than `prod`.
lgp_should_rerun_external <- function(prod, dep_paths) {
  if (!lgp_use_cache()) return(TRUE)
  force <- !lgp_smart_skip_external() || lgp_force_rerun() ||
    identical(Sys.getenv("LGP_FORCE_EXTERNAL", ""), "1")
  if (force) return(TRUE)
  if (!file.exists(prod)) return(TRUE)
  lgp_deps_newer_than_cache(dep_paths, prod)
}

## ---- Initialise ----------------------------------------------------------

.init_lgp_options()
