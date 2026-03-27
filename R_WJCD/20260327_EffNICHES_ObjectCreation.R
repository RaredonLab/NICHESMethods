# Spatial PlotDev
# WJCD

# pack
library(Seurat)
library(SeuratObject)
library(arrow)
library(cowplot)
library(ggplot2)
library(NICHES)
library(patchwork)
library(stringr)
library(scales)
library(dplyr)
library(colorspace)
library(tidyr)
library(ggrepel)

#### Discover Spatial Source ####

#' Discover and select a spatial coordinate source (returns only an index spec)
#'
#' This function scans meta.data, Images(), and reductions for coordinate data,
#' presents numbered options, and returns ONLY a compact "index spec" that a
#' downstream function can use to reliably re-fetch coordinates (and polygon
#' vertices when available) from the Seurat object.
#'
#' The returned object is intentionally minimal:
#' - source$type in c("meta", "image", "reduction")
#' - source$meta_xy OR source$image OR source$reduction
#' - source$has_polygons (TRUE/FALSE)
#' - source$xy_cols for image sources (and optional dims for reduction sources)
#'
#' Downstream should use `spatial_fetch()` (provided below) to materialize
#' centroids/polygons from this spec.
#'
#' @param obj Seurat object
#' @param meta_xy Candidate centroid coordinate column names in meta.data (default c("x","y"))
#' @param reduction_patterns Heuristic patterns for spatial-like reductions
#' @param scan_all_reductions If TRUE, scan all reductions for XY-like embeddings (default FALSE: spatial-like only)
#' @param polygon_rows_per_cell_thresh Threshold to classify image coordinate tables as polygons (default 1.5)
#' @param interactive If TRUE, prompt user to choose an option (default TRUE)
#' @param selection Integer option number to select programmatically (default NULL)
#' @param allow_cancel If TRUE, user can enter 0 to cancel and return NULL (default TRUE)
#'
#' @return A minimal index spec list (or NULL if canceled).
discover_spatial_source <- function(
    obj,
    meta_xy = c("x", "y"),
    reduction_patterns = c("spatial", "space", "image", "coord", "xenium", "visium"),
    scan_all_reductions = FALSE,
    polygon_rows_per_cell_thresh = 1.5,
    interactive = TRUE,
    selection = NULL,
    allow_cancel = TRUE
) {
  stopifnot(inherits(obj, "Seurat"))
  
  # -------------------------
  # helpers
  # -------------------------
  safe_num <- function(x) suppressWarnings(as.numeric(x))
  fmt_num  <- function(x, digits = 3) {
    if (length(x) == 0 || all(!is.finite(x))) return("NA")
    formatC(x, format = "f", digits = digits)
  }
  fmt_rng  <- function(r) paste0(fmt_num(r[1]), " .. ", fmt_num(r[2]))
  n_distinct <- function(x) length(unique(x))
  
  pick_xy_cols <- function(cn) {
    cn_l <- tolower(cn)
    pick_one <- function(cands) {
      hit <- intersect(cands, cn_l)
      if (length(hit) > 0) return(cn[match(hit[1], cn_l)])
      return(NULL)
    }
    xcol <- pick_one(c("x", "row", "imagerow", "pxl_row_in_fullres", "pxl_row", "xcoord", "x_coord"))
    ycol <- pick_one(c("y", "col", "imagecol", "pxl_col_in_fullres", "pxl_col", "ycoord", "y_coord"))
    
    if (is.null(xcol) || is.null(ycol)) {
      if ("row" %in% cn && "col" %in% cn) {
        xcol <- "row"; ycol <- "col"
      } else if (length(cn) >= 2) {
        if (is.null(xcol)) xcol <- cn[1]
        if (is.null(ycol)) ycol <- cn[min(2, length(cn))]
      }
    }
    list(xcol = xcol, ycol = ycol)
  }
  
  pick_cell_ids <- function(coords) {
    cn <- colnames(coords)
    if ("cell" %in% cn) return(as.character(coords$cell))
    rn <- rownames(coords)
    if (!is.null(rn) && length(rn) == nrow(coords)) return(as.character(rn))
    as.character(seq_len(nrow(coords)))
  }
  
  classify_coords_table <- function(df) {
    n_rows  <- nrow(df)
    n_cells <- n_distinct(df$cell)
    rows_per_cell <- n_rows / max(n_cells, 1)
    tab <- table(df$cell)
    frac_multi <- mean(tab > 1)
    
    is_polygon <- (rows_per_cell > polygon_rows_per_cell_thresh) ||
      (frac_multi >= 0.25 && rows_per_cell > 1.1)
    
    list(is_polygon = is_polygon, n_rows = n_rows, n_cells = n_cells, rows_per_cell = rows_per_cell)
  }
  
  scan_reductions_for_xy <- function(obj, red_names, spatial_like, scan_all) {
    # Return minimal reduction candidate descriptors (NOT embeddings)
    target <- if (isTRUE(scan_all)) red_names else spatial_like
    candidates <- list()
    if (length(target) == 0) return(candidates)
    
    for (rn in target) {
      emb <- tryCatch(Seurat::Embeddings(obj, reduction = rn), error = function(e) NULL)
      if (is.null(emb) || nrow(emb) == 0 || ncol(emb) < 2) next
      
      cn <- colnames(emb)
      cn_l <- tolower(cn)
      
      x_idx <- which(cn_l %in% c("x", "xcoord", "x_coord", "imagerow", "row", "pxl_row_in_fullres", "pxl_row"))
      y_idx <- which(cn_l %in% c("y", "ycoord", "y_coord", "imagecol", "col", "pxl_col_in_fullres", "pxl_col"))
      
      used_first_two <- FALSE
      if (length(x_idx) == 0 || length(y_idx) == 0) {
        used_first_two <- TRUE
        x_idx <- 1
        y_idx <- 2
      } else {
        x_idx <- x_idx[1]
        y_idx <- y_idx[1]
      }
      
      x <- safe_num(emb[, x_idx])
      y <- safe_num(emb[, y_idx])
      keep <- is.finite(x) & is.finite(y)
      if (!any(keep)) next
      
      # lightweight summary for menu
      xr <- range(x[keep], na.rm = TRUE)
      yr <- range(y[keep], na.rm = TRUE)
      
      candidates[[paste0("reduction:", rn)]] <- list(
        type = "reduction",
        reduction = rn,
        x_dim = cn[x_idx],
        y_dim = cn[y_idx],
        x_idx = x_idx,
        y_idx = y_idx,
        used_first_two_dims = used_first_two,
        n_cells = sum(keep),
        x_range = xr,
        y_range = yr
      )
    }
    
    candidates
  }
  
  # -------------------------
  # gather candidates (minimal specs + light summaries)
  # -------------------------
  candidates <- list()
  
  # meta.data candidate
  md <- obj@meta.data
  if (all(meta_xy %in% colnames(md))) {
    xv <- safe_num(md[[meta_xy[1]]])
    yv <- safe_num(md[[meta_xy[2]]])
    keep <- is.finite(xv) & is.finite(yv)
    if (sum(keep) > 0) {
      candidates[["meta.data"]] <- list(
        type = "meta",
        meta_xy = meta_xy,
        has_polygons = FALSE,
        n_cells = sum(keep),
        x_range = range(xv[keep], na.rm = TRUE),
        y_range = range(yv[keep], na.rm = TRUE)
      )
    }
  }
  
  # image candidates
  imgs <- tryCatch(Seurat::Images(obj), error = function(e) character(0))
  if (length(imgs) > 0) {
    for (im in imgs) {
      coords <- tryCatch(Seurat::GetTissueCoordinates(obj, image = im), error = function(e) NULL)
      if (is.null(coords) || nrow(coords) == 0) next
      
      cn <- colnames(coords)
      xy <- pick_xy_cols(cn)
      if (is.null(xy$xcol) || is.null(xy$ycol)) next
      
      cell_id <- pick_cell_ids(coords)
      df <- data.frame(
        cell = as.character(cell_id),
        x = safe_num(coords[[xy$xcol]]),
        y = safe_num(coords[[xy$ycol]]),
        stringsAsFactors = FALSE
      )
      df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
      if (nrow(df) == 0) next
      
      cls <- classify_coords_table(df)
      df_xr <- range(df$x, na.rm = TRUE)
      df_yr <- range(df$y, na.rm = TRUE)
      
      candidates[[paste0("image:", im)]] <- list(
        type = "image",
        image = im,
        xy_cols = c(xy$xcol, xy$ycol),
        has_polygons = isTRUE(cls$is_polygon),
        # for menu only:
        n_cells = cls$n_cells,
        approx_vertices_per_cell = cls$rows_per_cell,
        x_range = df_xr,
        y_range = df_yr
      )
    }
  }
  
  # reduction candidates
  red_names <- tryCatch(names(obj@reductions), error = function(e) character(0))
  spatial_like <- red_names[
    vapply(
      red_names,
      function(rn) any(grepl(paste0(reduction_patterns, collapse = "|"), rn, ignore.case = TRUE)),
      logical(1)
    )
  ]
  red_cands <- scan_reductions_for_xy(obj, red_names, spatial_like, scan_all_reductions)
  candidates <- c(candidates, red_cands)
  
  # -------------------------
  # print menu + select
  # -------------------------
  cat("\n==============================\n")
  cat("Spatial Coordinate Source Picker\n")
  cat("==============================\n\n")
  
  if (length(candidates) == 0) {
    cat("No coordinate sources detected.\n")
    cat("Tip: try scan_all_reductions=TRUE if coordinates are stored in non-spatial reductions.\n")
    cat("\n==============================\n")
    return(invisible(NULL))
  }
  
  keys <- names(candidates)
  
  cat("## Options (choose one)\n")
  for (i in seq_along(keys)) {
    k <- keys[i]
    cnd <- candidates[[k]]
    
    poly_flag <- if (isTRUE(cnd$has_polygons)) " [HAS PERIMETER VERTICES]" else ""
    
    if (identical(cnd$type, "meta")) {
      cat(sprintf(
        "%d) %s%s\n   - meta.data cols: %s / %s\n   - cells: %s | X: %s | Y: %s\n",
        i, k, poly_flag,
        cnd$meta_xy[1], cnd$meta_xy[2],
        format(cnd$n_cells, big.mark = ","),
        fmt_rng(cnd$x_range), fmt_rng(cnd$y_range)
      ))
    } else if (identical(cnd$type, "image")) {
      extra <- paste0("image=", cnd$image, "; cols=", cnd$xy_cols[1], "/", cnd$xy_cols[2])
      if (isTRUE(cnd$has_polygons)) {
        extra <- paste0(extra, "; approx_vertices/cell=", fmt_num(cnd$approx_vertices_per_cell, 2))
      }
      cat(sprintf(
        "%d) %s%s\n   - %s\n   - cells: %s | X: %s | Y: %s\n",
        i, k, poly_flag,
        extra,
        format(cnd$n_cells, big.mark = ","),
        fmt_rng(cnd$x_range), fmt_rng(cnd$y_range)
      ))
    } else if (identical(cnd$type, "reduction")) {
      extra <- paste0("reduction=", cnd$reduction, "; dims=", cnd$x_dim, "/", cnd$y_dim,
                      if (isTRUE(cnd$used_first_two_dims)) " (used first two dims)" else "")
      cat(sprintf(
        "%d) %s%s\n   - %s\n   - cells: %s | X: %s | Y: %s\n",
        i, k, poly_flag,
        extra,
        format(cnd$n_cells, big.mark = ","),
        fmt_rng(cnd$x_range), fmt_rng(cnd$y_range)
      ))
    }
  }
  
  if (isTRUE(allow_cancel)) cat("0) Cancel / return NULL\n")
  
  chosen_idx <- NULL
  if (!is.null(selection)) {
    chosen_idx <- as.integer(selection)
  } else if (isTRUE(interactive)) {
    ans <- readline(prompt = "\nEnter option number: ")
    chosen_idx <- suppressWarnings(as.integer(ans))
  } else {
    cat("\nNo selection made (interactive=FALSE and selection=NULL). Returning NULL.\n")
    return(invisible(NULL))
  }
  
  if (!is.finite(chosen_idx)) stop("Invalid selection: must be an integer option number.")
  if (chosen_idx == 0 && isTRUE(allow_cancel)) return(NULL)
  if (chosen_idx < 1 || chosen_idx > length(keys)) stop("Selection out of range.")
  
  chosen_key <- keys[chosen_idx]
  chosen <- candidates[[chosen_key]]
  
  # -------------------------
  # return MINIMAL index spec
  # -------------------------
  if (identical(chosen$type, "meta")) {
    return(list(
      type = "meta",
      meta_xy = chosen$meta_xy,
      has_polygons = FALSE
    ))
  }
  
  if (identical(chosen$type, "image")) {
    return(list(
      type = "image",
      image = chosen$image,
      xy_cols = chosen$xy_cols,
      has_polygons = isTRUE(chosen$has_polygons)
    ))
  }
  
  # reduction
  return(list(
    type = "reduction",
    reduction = chosen$reduction,
    x_idx = chosen$x_idx,
    y_idx = chosen$y_idx,
    x_dim = chosen$x_dim,
    y_dim = chosen$y_dim,
    used_first_two_dims = isTRUE(chosen$used_first_two_dims),
    has_polygons = FALSE
  ))
}


#' Materialize spatial coordinates (and optional polygon vertices) from a spec
#'
#' Downstream functions should call this to fetch data from the Seurat object
#' using the minimal spec returned by discover_spatial_source().
#'
#' @param obj Seurat object
#' @param spec List returned by discover_spatial_source()
#'
#' @return list(centroids = data.frame(cell,x,y), polygons = data.frame(cell,x,y) or NULL)
spatial_fetch <- function(obj, spec) {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(is.list(spec), !is.null(spec$type))
  
  safe_num <- function(x) suppressWarnings(as.numeric(x))
  
  pick_cell_ids <- function(coords) {
    cn <- colnames(coords)
    if ("cell" %in% cn) return(as.character(coords$cell))
    rn <- rownames(coords)
    if (!is.null(rn) && length(rn) == nrow(coords)) return(as.character(rn))
    as.character(seq_len(nrow(coords)))
  }
  
  if (identical(spec$type, "meta")) {
    md <- obj@meta.data
    if (!all(spec$meta_xy %in% colnames(md))) stop("meta_xy columns not found in obj@meta.data.")
    xv <- safe_num(md[[spec$meta_xy[1]]])
    yv <- safe_num(md[[spec$meta_xy[2]]])
    keep <- is.finite(xv) & is.finite(yv)
    centroids <- data.frame(cell = rownames(md)[keep], x = xv[keep], y = yv[keep], stringsAsFactors = FALSE)
    return(list(centroids = centroids, polygons = NULL))
  }
  
  if (identical(spec$type, "image")) {
    coords <- Seurat::GetTissueCoordinates(obj, image = spec$image)
    if (is.null(coords) || nrow(coords) == 0) stop("No coordinates returned by GetTissueCoordinates() for that image.")
    if (!all(spec$xy_cols %in% colnames(coords))) stop("Specified xy_cols not found in image coordinates.")
    cell_id <- pick_cell_ids(coords)
    
    df <- data.frame(
      cell = as.character(cell_id),
      x = safe_num(coords[[spec$xy_cols[1]]]),
      y = safe_num(coords[[spec$xy_cols[2]]]),
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
    if (nrow(df) == 0) stop("Image coordinate table had no finite x/y values.")
    
    if (isTRUE(spec$has_polygons)) {
      # polygons are raw vertex rows; centroids are per-cell vertex-means
      cent <- stats::aggregate(cbind(x, y) ~ cell, data = df, FUN = function(z) mean(z, na.rm = TRUE))
      return(list(centroids = cent, polygons = df))
    } else {
      # ensure unique per cell
      df_cent <- df[!duplicated(df$cell), , drop = FALSE]
      return(list(centroids = df_cent, polygons = NULL))
    }
  }
  
  if (identical(spec$type, "reduction")) {
    emb <- Seurat::Embeddings(obj, reduction = spec$reduction)
    if (is.null(emb) || nrow(emb) == 0 || ncol(emb) < 2) stop("Reduction embeddings not available.")
    if (spec$x_idx > ncol(emb) || spec$y_idx > ncol(emb)) stop("x_idx/y_idx out of bounds for embeddings.")
    
    x <- safe_num(emb[, spec$x_idx])
    y <- safe_num(emb[, spec$y_idx])
    keep <- is.finite(x) & is.finite(y)
    
    centroids <- data.frame(
      cell = rownames(emb)[keep],
      x = x[keep],
      y = y[keep],
      stringsAsFactors = FALSE
    )
    return(list(centroids = centroids, polygons = NULL))
  }
  
  stop("Unknown spec$type. Expected one of: meta, image, reduction.")
}

#### Discover Expression Source ####

#' Discover and select an expression data source (returns only an index spec)
#'
#' This is the expression analogue of discover_spatial_source(): it scans assays and
#' layers, summarizes what is available + heuristic type, presents numbered options,
#' and returns ONLY a compact "index spec" that downstream code can use to reliably
#' re-fetch the matrix.
#'
#' The returned spec is intentionally minimal:
#' - $assay: assay name
#' - $layer: layer name (e.g. "counts", "data", "scale.data")
#' - $kind: heuristic kind ("raw_like", "normalized_like", "scaled_like", "control_like", etc.)
#' - $control_like: TRUE/FALSE (semantic)
#' - $counts_data_identical: TRUE/FALSE (informational; NOT disqualifying)
#'
#' Downstream should use `expression_fetch()` (provided below) to materialize the matrix.
#'
#' @param obj Seurat object
#' @param assays Character vector of assays to scan (default NULL = all assays)
#' @param layers Layers to scan (default c("counts","data","scale.data"))
#' @param sample_n Number of sampled values for heuristics (default 2e5; sparse samples from @x)
#' @param show_quantiles If TRUE, prints sampled 1/50/99% quantiles (default TRUE)
#' @param integer_tol Tolerance for integer-ish check (default 1e-8)
#' @param raw_integer_frac_thresh Fraction integer-ish to call raw_like (default 0.98)
#' @param control_assay_patterns Regex patterns for known control assays
#' @param control_value_max Fallback control-like detection: max value threshold when low-entropy + counts==data (default 2)
#' @param control_unique_max Fallback control-like detection: unique value count threshold when low-entropy + counts==data (default 3)
#' @param interactive If TRUE, prompt user to choose an option (default TRUE)
#' @param selection Integer option number to select programmatically (default NULL)
#' @param allow_cancel If TRUE, user can enter 0 to cancel and return NULL (default TRUE)
#'
#' @return A minimal index spec list (or NULL if canceled).
discover_expression_source <- function(
    obj,
    assays = NULL,
    layers = c("counts", "data", "scale.data"),
    sample_n = 2e5,
    show_quantiles = TRUE,
    integer_tol = 1e-8,
    raw_integer_frac_thresh = 0.98,
    control_assay_patterns = c("^ControlCodeword$", "^ControlProbe$", "^BlankCodeword$", "^BlankProbe$", "^GenomicControl$"),
    control_value_max = 2,
    control_unique_max = 3,
    interactive = TRUE,
    selection = NULL,
    allow_cancel = TRUE
) {
  stopifnot(inherits(obj, "Seurat"))
  
  # -------------------------
  # helpers
  # -------------------------
  fmt_num <- function(x, digits = 4) {
    if (length(x) == 0 || all(!is.finite(x))) return("NA")
    formatC(x, format = "g", digits = digits)
  }
  fmt_rng <- function(r) paste0(fmt_num(r[1]), " .. ", fmt_num(r[2]))
  
  frac_integerish <- function(x, tol = 1e-8) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    mean(abs(x - round(x)) <= tol)
  }
  
  is_control_assay <- function(assay, patterns) {
    any(vapply(patterns, function(p) grepl(p, assay, ignore.case = FALSE), logical(1)))
  }
  
  # Unified matrix fetcher across Seurat v4/v5
  get_mat <- function(obj, assay, layer) {
    a <- obj[[assay]]
    
    # Seurat v5: Assay5 layers
    m <- tryCatch(SeuratObject::LayerData(a, layer = layer), error = function(e) NULL)
    if (!is.null(m)) return(m)
    
    # Seurat v4: classic slots (slot arg is defunct in v5, but fallback for older objects)
    m2 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = layer), error = function(e) NULL)
    if (!is.null(m2)) return(m2)
    
    NULL
  }
  
  sample_values <- function(m, n = 2e5) {
    if (is.null(m)) return(numeric(0))
    if (inherits(m, "dgCMatrix") || inherits(m, "dgTMatrix") || inherits(m, "dgRMatrix")) {
      v <- m@x
      if (length(v) == 0) return(numeric(0))
      if (n <= 0 || length(v) <= n) return(as.numeric(v))
      return(as.numeric(sample(v, size = n)))
    }
    v <- as.numeric(m)
    v <- v[is.finite(v)]
    if (length(v) == 0) return(numeric(0))
    if (n <= 0 || length(v) <= n) return(v)
    sample(v, size = n)
  }
  
  matrices_identical <- function(a, b) {
    if (is.null(a) || is.null(b)) return(FALSE)
    if (!identical(dim(a), dim(b))) return(FALSE)
    
    if (inherits(a, "dgCMatrix") && inherits(b, "dgCMatrix")) {
      return(
        identical(a@i, b@i) &&
          identical(a@p, b@p) &&
          identical(a@x, b@x) &&
          identical(a@Dim, b@Dim)
      )
    }
    
    if (is.matrix(a) && is.matrix(b)) {
      return(isTRUE(all.equal(a, b, tolerance = 0)))
    }
    
    FALSE
  }
  
  classify_matrix <- function(m, sampled_vals, assay, layer, counts_data_identical = FALSE) {
    if (is.null(m)) return(list(kind = "absent", usable = FALSE, control_like = FALSE))
    
    sv <- sampled_vals
    sv <- sv[is.finite(sv)]
    
    d <- dim(m)
    if (length(d) == 2 && (d[1] == 0 || d[2] == 0)) {
      return(list(kind = "present_empty", usable = FALSE, control_like = FALSE))
    }
    if (length(sv) == 0) {
      return(list(kind = "present_empty", usable = FALSE, control_like = FALSE))
    }
    
    fi <- frac_integerish(sv, tol = integer_tol)
    has_neg <- any(sv < 0, na.rm = TRUE)
    
    # Control assays by name
    if (is_control_assay(assay, control_assay_patterns)) {
      return(list(kind = "control_like", usable = TRUE, control_like = TRUE, frac_integerish = fi))
    }
    
    # Fallback control-like: low-entropy discrete + counts==data signature
    u <- unique(sv)
    if (layer %in% c("counts", "data") &&
        isTRUE(counts_data_identical) &&
        length(u) <= control_unique_max &&
        max(u, na.rm = TRUE) <= control_value_max) {
      return(list(kind = "control_like", usable = TRUE, control_like = TRUE, frac_integerish = fi))
    }
    
    if (has_neg) {
      return(list(kind = "scaled_like", usable = TRUE, control_like = FALSE, frac_integerish = fi))
    }
    
    if (is.finite(fi) && fi >= raw_integer_frac_thresh) {
      return(list(kind = "raw_like", usable = TRUE, control_like = FALSE, frac_integerish = fi))
    }
    
    list(kind = "normalized_like", usable = TRUE, control_like = FALSE, frac_integerish = fi)
  }
  
  # -------------------------
  # decide assays to scan
  # -------------------------
  all_assays <- tryCatch(names(obj@assays), error = function(e) character(0))
  if (is.null(assays)) assays <- all_assays
  assays <- intersect(assays, all_assays)
  
  # -------------------------
  # scan and build options
  # -------------------------
  options <- list()
  
  if (length(assays) == 0) {
    cat("\n==============================\n")
    cat("Expression Data Source Picker\n")
    cat("==============================\n\n")
    cat("No assays found to scan.\n")
    cat("\n==============================\n")
    return(invisible(NULL))
  }
  
  for (assay in assays) {
    mats <- list()
    for (layer in layers) mats[[layer]] <- get_mat(obj, assay, layer)
    
    cd_ident <- matrices_identical(mats[["counts"]], mats[["data"]])
    assay_is_ctrl <- is_control_assay(assay, control_assay_patterns)
    
    for (layer in layers) {
      m <- mats[[layer]]
      if (is.null(m)) next
      
      # sample
      sv <- sample_values(m, n = sample_n)
      if (length(sv) == 0) next
      
      cls <- classify_matrix(
        m = m,
        sampled_vals = sv,
        assay = assay,
        layer = layer,
        counts_data_identical = cd_ident
      )
      if (!isTRUE(cls$usable)) next
      
      rng <- range(sv, na.rm = TRUE)
      qs <- NULL
      if (isTRUE(show_quantiles)) {
        qs <- tryCatch(stats::quantile(sv, probs = c(0.01, 0.5, 0.99), na.rm = TRUE), error = function(e) NULL)
      }
      
      key <- paste0("assay:", assay, " | layer:", layer)
      
      options[[key]] <- list(
        # minimal spec fields (what we will return)
        spec = list(
          assay = assay,
          layer = layer,
          kind = cls$kind,
          control_like = isTRUE(cls$control_like) || isTRUE(assay_is_ctrl),
          counts_data_identical = isTRUE(cd_ident)
        ),
        # lightweight menu info only
        menu = list(
          assay = assay,
          layer = layer,
          kind = cls$kind,
          control_like = isTRUE(cls$control_like) || isTRUE(assay_is_ctrl),
          counts_data_identical = isTRUE(cd_ident),
          class = class(m),
          dims = dim(m),
          value_range = rng,
          quantiles = qs,
          frac_integerish = cls$frac_integerish
        )
      )
    }
  }
  
  # -------------------------
  # print menu + select
  # -------------------------
  cat("\n==============================\n")
  cat("Expression Data Source Picker\n")
  cat("==============================\n\n")
  
  cat("- DefaultAssay(obj): ", tryCatch(Seurat::DefaultAssay(obj), error = function(e) NA_character_), "\n", sep = "")
  cat("- Assays scanned: ", paste(assays, collapse = ", "), "\n", sep = "")
  cat("- Layers checked: ", paste(layers, collapse = ", "), "\n\n", sep = "")
  
  if (length(options) == 0) {
    cat("No usable expression matrices detected.\n")
    cat("\n==============================\n")
    return(invisible(NULL))
  }
  
  keys <- names(options)
  
  cat("## Options (choose one)\n")
  for (i in seq_along(keys)) {
    k <- keys[i]
    m <- options[[k]]$menu
    
    flags <- c()
    if (isTRUE(m$control_like)) flags <- c(flags, "CONTROL-LIKE/QC")
    if (isTRUE(m$counts_data_identical)) flags <- c(flags, "counts==data")
    flag_str <- if (length(flags) > 0) paste0(" [", paste(flags, collapse = "; "), "]") else ""
    
    qline <- ""
    if (isTRUE(show_quantiles) && !is.null(m$quantiles)) {
      q <- m$quantiles
      qline <- paste0(
        "   - sampled quantiles (1/50/99%): ",
        fmt_num(unname(q[1])), " / ", fmt_num(unname(q[2])), " / ", fmt_num(unname(q[3])), "\n"
      )
    }
    
    fi_line <- ""
    if (!is.null(m$frac_integerish) && is.finite(m$frac_integerish)) {
      fi_line <- paste0("   - frac integer-ish: ", fmt_num(m$frac_integerish, 3), "\n")
    }
    
    cat(sprintf(
      "%d) %s%s\n   - assay: %s | layer: %s | kind: %s\n   - dims (features x cells): %s x %s\n   - sampled range: %s\n%s%s",
      i, k, flag_str,
      m$assay, m$layer, m$kind,
      m$dims[1], m$dims[2],
      fmt_rng(m$value_range),
      qline,
      fi_line
    ))
  }
  
  if (isTRUE(allow_cancel)) cat("0) Cancel / return NULL\n")
  
  chosen_idx <- NULL
  if (!is.null(selection)) {
    chosen_idx <- as.integer(selection)
  } else if (isTRUE(interactive)) {
    ans <- readline(prompt = "\nEnter option number: ")
    chosen_idx <- suppressWarnings(as.integer(ans))
  } else {
    cat("\nNo selection made (interactive=FALSE and selection=NULL). Returning NULL.\n")
    return(invisible(NULL))
  }
  
  if (!is.finite(chosen_idx)) stop("Invalid selection: must be an integer option number.")
  if (chosen_idx == 0 && isTRUE(allow_cancel)) return(NULL)
  if (chosen_idx < 1 || chosen_idx > length(keys)) stop("Selection out of range.")
  
  chosen_key <- keys[chosen_idx]
  chosen_spec <- options[[chosen_key]]$spec
  
  # return MINIMAL index spec only
  return(chosen_spec)
}

#' Materialize an expression matrix from a spec
#'
#' Downstream functions should call this to fetch the matrix from the Seurat object
#' using the minimal spec returned by discover_expression_source().
#'
#' @param obj Seurat object
#' @param spec List returned by discover_expression_source()
#'
#' @return Matrix-like object (dgCMatrix or matrix) for the chosen assay+layer.
expression_fetch <- function(obj, spec) {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(is.list(spec), !is.null(spec$assay), !is.null(spec$layer))
  
  assay <- spec$assay
  layer <- spec$layer
  
  a <- obj[[assay]]
  
  # Seurat v5: Assay5 layers
  m <- tryCatch(SeuratObject::LayerData(a, layer = layer), error = function(e) NULL)
  if (!is.null(m)) return(m)
  
  # Seurat v4: classic slots
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = layer), error = function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  stop("Could not fetch expression matrix for assay='", assay, "' layer='", layer, "'.")
}

#### Data Extraction and Object Creation Function ####

# EffNICHES: Seurat-independent container + NICHES runner/extractor
#
# Assumes you already have:
#   - discover_spatial_source()
#   - spatial_fetch()
#   - discover_expression_source()
#   - expression_fetch()
#
# Requirements:
#   - Matrix
#   - Seurat / SeuratObject
#   - NICHES (RunNICHES)
#   - fst (for fst-style serialization helpers below)
#   - qs  (for sparse matrices + integer vectors)

.safe_meta_cols <- function(obj, include = NULL) {
  md <- obj@meta.data
  n  <- nrow(md)
  
  ok <- vapply(names(md), function(nm) {
    v <- md[[nm]]
    # disallow list-like and 2D structures
    if (is.list(v) && !is.factor(v)) return(FALSE)
    if (is.matrix(v) || is.data.frame(v)) return(FALSE)
    
    # must be length n (some broken columns are length 0)
    if (length(v) != n) return(FALSE)
    
    # allow atomic, factor, Date/POSIXct, etc.
    is.atomic(v) || is.factor(v) || inherits(v, c("Date", "POSIXct", "POSIXt"))
  }, logical(1))
  
  cols <- names(md)[ok]
  
  # Always keep include (if present and valid length); otherwise warn
  if (!is.null(include)) {
    include <- intersect(include, names(md))
    cols <- union(include, cols)
  }
  
  cols
}
.safe_num <- function(x) suppressWarnings(as.numeric(x))

# Build ragged polygon store:
# polygons_df: data.frame(cell, x, y) with potentially many rows per cell
# ids: canonical node ids, length N
.polygons_to_ragged <- function(polygons_df, ids) {
  if (is.null(polygons_df) || nrow(polygons_df) == 0) return(NULL)
  
  # Ensure columns exist
  stopifnot(all(c("cell", "x", "y") %in% colnames(polygons_df)))
  
  # Keep only cells in ids
  polygons_df$cell <- as.character(polygons_df$cell)
  keep <- polygons_df$cell %in% ids
  polygons_df <- polygons_df[keep, , drop = FALSE]
  if (nrow(polygons_df) == 0) return(NULL)
  
  # Map cell -> index
  id_to_idx <- setNames(seq_along(ids), ids)
  idx <- unname(id_to_idx[polygons_df$cell])
  ok <- is.finite(idx) & !is.na(idx)
  polygons_df <- polygons_df[ok, , drop = FALSE]
  idx <- idx[ok]
  
  # Order by idx to build CSR-style pointers
  o <- order(idx)
  idx <- idx[o]
  xy  <- cbind(.safe_num(polygons_df$x[o]), .safe_num(polygons_df$y[o]))
  storage.mode(xy) <- "double"
  
  N <- length(ids)
  counts <- tabulate(idx, nbins = N)
  start <- integer(N + 1)
  start[1] <- 1L
  if (N > 0) start[2:(N + 1)] <- 1L + cumsum(counts)
  
  list(
    xy = xy,              # [V x 2]
    start = start         # [N+1], 1-based ranges: start[i]:(start[i+1]-1)
  )
}

# Set coordinates for NICHES in meta.data as temporary columns, restore after
.with_temp_xy <- function(obj, centroids_df, x_col = "eff_x", y_col = "eff_y") {
  md <- obj@meta.data
  n <- nrow(md)
  
  # Preserve existing columns if present
  had_x <- x_col %in% colnames(md)
  had_y <- y_col %in% colnames(md)
  old_x <- if (had_x) md[[x_col]] else NULL
  old_y <- if (had_y) md[[y_col]] else NULL
  
  centroids_df$cell <- as.character(centroids_df$cell)
  row_idx <- match(centroids_df$cell, rownames(md))
  ok <- is.finite(row_idx) & !is.na(row_idx)
  
  if (!any(ok)) {
    stop("No centroid cell IDs matched obj@meta.data rownames. ",
         "This usually means your spatial source cell IDs don't match Seurat cell names.")
  }
  
  # Force full-length numeric columns (robust for data.frame and S4Vectors::DataFrame)
  md[[x_col]] <- rep(NA_real_, n)
  md[[y_col]] <- rep(NA_real_, n)
  
  md[[x_col]][row_idx[ok]] <- suppressWarnings(as.numeric(centroids_df$x[ok]))
  md[[y_col]][row_idx[ok]] <- suppressWarnings(as.numeric(centroids_df$y[ok]))
  
  # Sanity check: must have some finite coords
  if (!any(is.finite(md[[x_col]]) & is.finite(md[[y_col]]))) {
    stop("Injected XY columns are all non-finite. Check your spatial source and parsing.")
  }
  
  obj@meta.data <- md
  
  restore <- function(obj2) {
    md2 <- obj2@meta.data
    if (had_x) md2[[x_col]] <- old_x else md2[[x_col]] <- NULL
    if (had_y) md2[[y_col]] <- old_y else md2[[y_col]] <- NULL
    obj2@meta.data <- md2
    obj2
  }
  
  list(obj = obj, x_col = x_col, y_col = y_col, restore = restore)
}  

# Temporarily force NICHES to use a chosen layer for an assay, restore after
# Strategy:
# - Seurat v5 Assay5: try DefaultLayer(assay) <- layer
# - Seurat v4: if layer != "data", copy chosen slot into "data" and restore
.with_temp_assay_layer <- function(obj, assay, layer) {
  stopifnot(assay %in% names(obj@assays))
  a <- obj[[assay]]
  
  # v5 path: DefaultLayer exists
  restore <- function(obj2) obj2
  did <- FALSE
  
  # Try SeuratObject::DefaultLayer setter
  did <- tryCatch({
    cur <- SeuratObject::DefaultLayer(a)
    if (!identical(cur, layer)) {
      SeuratObject::DefaultLayer(a) <- layer
      obj[[assay]] <- a
      restore <- function(obj2) {
        a2 <- obj2[[assay]]
        SeuratObject::DefaultLayer(a2) <- cur
        obj2[[assay]] <- a2
        obj2
      }
    }
    TRUE
  }, error = function(e) FALSE)
  
  if (isTRUE(did)) return(list(obj = obj, restore = restore))
  
  # v4 fallback: copy requested slot into data
  if (!layer %in% c("counts", "data", "scale.data")) {
    stop("Layer '", layer, "' not recognized for v4 fallback.")
  }
  if (layer == "data") return(list(obj = obj, restore = restore))
  
  old_data <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "data"), error = function(e) NULL)
  new_data <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = layer), error = function(e) NULL)
  if (is.null(new_data)) stop("Could not fetch assay='", assay, "' slot='", layer, "' for v4 fallback.")
  
  obj <- Seurat::SetAssayData(obj, assay = assay, slot = "data", new.data = new_data)
  
  restore <- function(obj2) {
    if (!is.null(old_data)) {
      obj2 <- Seurat::SetAssayData(obj2, assay = assay, slot = "data", new.data = old_data)
    }
    obj2
  }
  
  list(obj = obj, restore = restore)
}

# Robust assay/layer fetcher (v4/v5)
.get_mat <- function(obj, assay, layer) {
  a <- obj[[assay]]
  
  m <- tryCatch(SeuratObject::LayerData(a, layer = layer), error = function(e) NULL)
  if (!is.null(m)) return(m)
  
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = layer), error = function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  NULL
}

# Detect sender/receiver id columns in NICHES edge meta.data
.detect_edge_id_cols <- function(md) {
  cn <- colnames(md)
  
  pick <- function(cands) {
    hit <- intersect(cands, cn)
    if (length(hit) > 0) return(hit[1])
    # case-insensitive fallback
    cn_l <- tolower(cn)
    c_l  <- tolower(cands)
    hit2 <- intersect(c_l, cn_l)
    if (length(hit2) > 0) return(cn[match(hit2[1], cn_l)])
    NULL
  }
  
  sender <- pick(c("Sending", "sending", "sender", "Sender", "cell.Sending", "Cell.Sending", "from", "From"))
  recv   <- pick(c("Receiving", "receiving", "receiver", "Receiver", "cell.Receiving", "Cell.Receiving", "to", "To"))
  
  list(sender = sender, receiver = recv)
}

# Parse edge ids from edge "cell names" if no sender/receiver columns exist
.parse_edge_from_rownames <- function(rn) {
  rn <- as.character(rn)
  # common separators
  seps <- c("\\|", "___", "__", "->", "—", "-", ":")
  for (s in seps) {
    parts <- strsplit(rn, s, perl = TRUE)
    ok <- vapply(parts, length, integer(1)) >= 2
    if (any(ok)) {
      sender <- vapply(parts, function(z) z[1], character(1))
      receiver <- vapply(parts, function(z) z[2], character(1))
      return(list(sender = sender, receiver = receiver, sep = s))
    }
  }
  NULL
}

# Convert a NICHES edge object into EffNICHES edge store
.extract_edges_E_by_M <- function(niches_obj, eff_ids, edge_assay = "CellToCellSpatial", layer = "data") {
  stopifnot(edge_assay %in% names(niches_obj@assays))
  
  md <- niches_obj@meta.data
  cols <- .detect_edge_id_cols(md)
  
  if (!is.null(cols$sender) && !is.null(cols$receiver)) {
    sender_id <- as.character(md[[cols$sender]])
    recv_id   <- as.character(md[[cols$receiver]])
  } else {
    parsed <- .parse_edge_from_rownames(rownames(md))
    if (is.null(parsed)) {
      stop("Could not find sender/receiver identifiers in NICHES edge meta.data (no id columns; could not parse rownames).")
    }
    sender_id <- parsed$sender
    recv_id   <- parsed$receiver
  }
  
  # Map to node indices
  id_to_idx <- setNames(seq_along(eff_ids), eff_ids)
  s_idx <- unname(id_to_idx[sender_id])
  r_idx <- unname(id_to_idx[recv_id])
  
  ok <- is.finite(s_idx) & is.finite(r_idx) & !is.na(s_idx) & !is.na(r_idx)
  if (!any(ok)) stop("No edges could be mapped onto EffNICHES node IDs (sender/receiver mismatch).")
  
  # Edge ordering: keep only mapped edges
  s_idx <- as.integer(s_idx[ok])
  r_idx <- as.integer(r_idx[ok])
  
  ij <- cbind(s_idx, r_idx)
  colnames(ij) <- c("sender", "receiver")
  
  # Fetch mechanism-by-edge matrix: features x edges
  m_mech_by_edge <- .get_mat(niches_obj, assay = edge_assay, layer = layer)
  if (is.null(m_mech_by_edge)) stop("Could not fetch NICHES edge assay matrix.")
  
  # Subset columns to mapped edges (edges are "cells" in niches_obj)
  # We subset by ok in the *same order* we derived ids from meta.data rownames (aligned to columns via colnames)
  # NICHES usually aligns meta.data rows to columns; if not, we enforce match.
  edge_cols <- colnames(m_mech_by_edge)
  md_rn <- rownames(md)
  
  # Try to align meta.data rows to edge columns
  if (!is.null(edge_cols) && length(edge_cols) == ncol(m_mech_by_edge) && all(edge_cols %in% md_rn)) {
    # reorder meta-derived vectors to match matrix columns
    ord <- match(edge_cols, md_rn)
    # recompute ok, s_idx, r_idx using ordered meta rows
    sender_id2 <- sender_id[ord]
    recv_id2   <- recv_id[ord]
    s2 <- unname(id_to_idx[sender_id2])
    r2 <- unname(id_to_idx[recv_id2])
    ok2 <- is.finite(s2) & is.finite(r2) & !is.na(s2) & !is.na(r2)
    
    ij <- cbind(as.integer(s2[ok2]), as.integer(r2[ok2]))
    colnames(ij) <- c("sender", "receiver")
    m_mech_by_edge <- m_mech_by_edge[, ok2, drop = FALSE]
  } else {
    # fallback: assume ok corresponds to the matrix column order (best-effort)
    m_mech_by_edge <- m_mech_by_edge[, ok, drop = FALSE]
  }
  
  # Convert to edges x mechanisms sparse matrix
  w_E_by_M <- Matrix::t(m_mech_by_edge)  # dgCMatrix if input is sparse
  if (!inherits(w_E_by_M, "dgCMatrix")) w_E_by_M <- as(w_E_by_M, "dgCMatrix")
  
  list(
    ij = ij,                                # [E x 2]
    w  = w_E_by_M,                          # [E x M]
    mechanisms = colnames(w_E_by_M)
  )
}

# Public API: build_effniches()

#' Build a Seurat-independent EffNICHES container (nodes + expression)
#'
#' @param obj Seurat object
#' @param spatial_spec If provided, a spec from discover_spatial_source(); if NULL, discovery will run
#' @param expr_spec If provided, a spec from discover_expression_source(); if NULL, discovery will run
#' @param meta_cols Character vector of meta.data columns to carry into eff$nodes$meta (default NULL = none)
#' @param interactive If TRUE, allow interactive picker when specs are NULL
#' @param spatial_selection If using discovery non-interactively, select option number
#' @param expr_selection If using discovery non-interactively, select option number
#'
#' @return EffNICHES list with nodes + expr + sources (no NICHES results yet)
build_effniches <- function(
    obj,
    spatial_spec = NULL,
    expr_spec = NULL,
    meta_cols = NULL,
    interactive = TRUE,
    spatial_selection = NULL,
    expr_selection = NULL
) {
  stopifnot(inherits(obj, "Seurat"))
  
  # Resolve specs if not provided
  if (is.null(spatial_spec)) {
    spatial_spec <- discover_spatial_source(
      obj,
      interactive = interactive,
      selection = spatial_selection
    )
  }
  if (is.null(expr_spec)) {
    expr_spec <- discover_expression_source(
      obj,
      interactive = interactive,
      selection = expr_selection
    )
  }
  if (is.null(spatial_spec) || is.null(expr_spec)) {
    stop("Both spatial_spec and expr_spec must be resolved (not NULL).")
  }
  
  # Fetch materialized data from specs
  sp <- spatial_fetch(obj, spatial_spec)         # list(centroids df, polygons df or NULL)
  expr_mat <- expression_fetch(obj, expr_spec)   # features x cells
  
  if (is.null(sp$centroids) || nrow(sp$centroids) == 0) stop("No centroids returned by spatial_fetch().")
  if (is.null(expr_mat) || ncol(expr_mat) == 0 || nrow(expr_mat) == 0) stop("No expression matrix returned by expression_fetch().")
  
  # Canonical IDs: intersection of centroids and expression columns
  cent_ids <- as.character(sp$centroids$cell)
  expr_ids <- colnames(expr_mat)
  if (is.null(expr_ids)) stop("Expression matrix must have colnames (cell/spot IDs).")
  
  ids <- intersect(expr_ids, cent_ids)
  if (length(ids) == 0) stop("No overlapping IDs between spatial centroids and expression columns.")
  
  # Build node xy aligned to ids
  cent_map <- sp$centroids
  cent_map$cell <- as.character(cent_map$cell)
  o <- match(ids, cent_map$cell)
  if (any(!is.finite(o) | is.na(o))) stop("Internal error: failed to align centroids to ids.")
  xy <- cbind(.safe_num(cent_map$x[o]), .safe_num(cent_map$y[o]))
  storage.mode(xy) <- "double"
  colnames(xy) <- c("x", "y")
  
  # Subset & reorder expression to ids
  expr_mat <- expr_mat[, ids, drop = FALSE]
  if (!inherits(expr_mat, "dgCMatrix") && !is.matrix(expr_mat)) {
    # try coerce
    expr_mat <- tryCatch(as(expr_mat, "dgCMatrix"), error = function(e) expr_mat)
  }
  
  # Node metadata
  meta <- NULL
  if (!is.null(meta_cols)) {
    md <- obj@meta.data
    bad <- setdiff(meta_cols, colnames(md))
    if (length(bad) > 0) stop("Requested meta_cols not found in obj@meta.data: ", paste(bad, collapse = ", "))
    meta <- md[ids, meta_cols, drop = FALSE]
  }
  
  # Polygons ragged (if present)
  polygons <- NULL
  if (!is.null(sp$polygons) && nrow(sp$polygons) > 0) {
    polygons <- .polygons_to_ragged(sp$polygons, ids)
  }
  
  eff <- list(
    version = list(format_version = "0.1.0"),
    nodes = list(
      ids = ids,
      xy = xy,
      meta = meta,
      polygons = polygons,
      index = list(id_to_idx = setNames(seq_along(ids), ids))
    ),
    expr = list(
      assay = expr_spec$assay,
      layer = expr_spec$layer,
      kind  = expr_spec$kind,
      control_like = isTRUE(expr_spec$control_like),
      counts_data_identical = isTRUE(expr_spec$counts_data_identical),
      features = rownames(expr_mat),
      mat = expr_mat
    ),
    niches = NULL,
    sources = list(
      spatial_spec = spatial_spec,
      expr_spec = expr_spec
    )
  )
  
  class(eff) <- c("EffNICHES", "list")
  eff
}

# Public API: run_effniches() [BROKEN - USE import_effniches_ctc() TO IMPORT EDGE DATA]

#' Run NICHES using chosen sources and extract results into EffNICHES (Seurat-independent)
#'
#' @param obj Seurat object (input to RunNICHES)
#' @param eff EffNICHES object from build_effniches()
#' @param label_col Metadata column in obj@meta.data used for NICHES cell_types
#' @param LR.database Passed to RunNICHES (default "fantom5")
#' @param species Passed to RunNICHES (e.g. "mouse","rat","human")
#' @param k Passed to RunNICHES (default NULL)
#' @param rad.set Passed to RunNICHES (default NULL)
#' @param CellToCellSpatial default TRUE
#' @param CellToNeighborhood default FALSE
#' @param NeighborhoodToCell default FALSE
#' @param meta.data.to.map Which metadata columns to map into NICHES output (default NULL = names(obj@meta.data))
#' @param niches_assay_layer Layer to read from NICHES output assays (default "data")
#'
#' @return EffNICHES with $niches populated
run_effniches <- function(
    obj,
    eff,
    label_col,
    LR.database = "fantom5",
    species,
    k = NULL,
    rad.set = NULL,
    CellToCellSpatial = TRUE,
    CellToNeighborhood = FALSE,
    NeighborhoodToCell = FALSE,
    meta.data.to.map = NULL,
    niches_assay_layer = "data"
) {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(inherits(eff, "EffNICHES"))
  stopifnot(is.character(label_col), length(label_col) == 1)
  
  if (!label_col %in% colnames(obj@meta.data)) stop("label_col not found in obj@meta.data: ", label_col)
  if (is.null(meta.data.to.map)) meta.data.to.map <- names(obj@meta.data)
  if (is.null(meta.data.to.map)) {
    meta.data.to.map <- .safe_meta_cols(obj, include = label_col)
  }
  
  # Ensure chosen expression assay/layer is what NICHES uses
  expr_spec <- eff$sources$expr_spec
  spatial_spec <- eff$sources$spatial_spec
  
  # Materialize centroids for XY injection (must match Seurat cell IDs)
  sp <- spatial_fetch(obj, spatial_spec)
  
  # Temporarily set coords in meta.data
  xy_ctx <- .with_temp_xy(obj, sp$centroids, x_col = ".eff_x", y_col = ".eff_y")
  obj2 <- xy_ctx$obj
  
  # Temporarily set assay layer for NICHES
  layer_ctx <- .with_temp_assay_layer(obj2, assay = expr_spec$assay, layer = expr_spec$layer)
  obj3 <- layer_ctx$obj
  
  # after xy_ctx is created, before RunNICHES()
  md <- obj3@meta.data
  
  use_cells <- rownames(md)[!is.na(md[[label_col]])]
  n_use <- length(use_cells)
  
  finite_xy <- is.finite(md[use_cells, xy_ctx$x_col]) & is.finite(md[use_cells, xy_ctx$y_col])
  
  if (!any(finite_xy)) {
    stop("For the cells NICHES will analyze (non-NA ", label_col, "), none have finite xy. ",
         "Your spatial centroids are not being mapped onto those cell IDs.")
  }
  
  if (sum(finite_xy) < n_use) {
    warning("Only ", sum(finite_xy), " / ", n_use, " analyzable cells have finite xy; NICHES may drop others.")
  }
  
  # Run NICHES
  niches_obj <- NULL
  niches_obj <- tryCatch(
    RunNICHES(
      object = obj3,
      assay = expr_spec$assay,
      cell_types = label_col,
      LR.database = LR.database,
      species = species,
      position.x = ".eff_x",
      position.y = ".eff_y",
      k = k,
      rad.set = rad.set,
      meta.data.to.map = meta.data.to.map,
      CellToCellSpatial = isTRUE(CellToCellSpatial),
      CellToCell = FALSE,
      SystemToCell = FALSE,
      CellToSystem = FALSE,
      CellToNeighborhood = isTRUE(CellToNeighborhood),
      NeighborhoodToCell = isTRUE(NeighborhoodToCell)
    ),
    error = function(e) {
      # restore before stopping
      obj3r <- tryCatch(layer_ctx$restore(obj3), error = function(e2) obj3)
      obj2r <- tryCatch(xy_ctx$restore(obj3r), error = function(e2) obj3r)
      stop("RunNICHES failed: ", conditionMessage(e))
    }
  )
  
  # Restore Seurat object (best-effort)
  obj3 <- tryCatch(layer_ctx$restore(obj3), error = function(e) obj3)
  obj3 <- tryCatch(xy_ctx$restore(obj3), error = function(e) obj3)
  
  # Extract results into eff$niches
  out <- list(
    params = list(
      LR.database = LR.database,
      species = species,
      k = k,
      rad.set = rad.set,
      CellToCellSpatial = isTRUE(CellToCellSpatial),
      CellToNeighborhood = isTRUE(CellToNeighborhood),
      NeighborhoodToCell = isTRUE(NeighborhoodToCell),
      expression_source = expr_spec,
      spatial_source = spatial_spec
    ),
    mechanisms = NULL,
    CellToCellSpatial = NULL,
    CellToNeighborhood = NULL,
    NeighborhoodToCell = NULL
  )
  
  # CellToCellSpatial edges: store separately as required
  if (isTRUE(CellToCellSpatial)) {
    if (!"CellToCellSpatial" %in% names(niches_obj@assays)) {
      stop("Requested CellToCellSpatial=TRUE but NICHES output has no 'CellToCellSpatial' assay.")
    }
    edges <- .extract_edges_E_by_M(niches_obj, eff$nodes$ids, edge_assay = "CellToCellSpatial", layer = niches_assay_layer)
    out$CellToCellSpatial <- list(
      ij = edges$ij,                         # [E x 2]
      w  = edges$w,                          # [E x M] dgCMatrix (zeros implicitly filtered)
      mechanisms = edges$mechanisms,
      directed = TRUE
    )
    out$mechanisms <- edges$mechanisms
  }
  
  # Node-level outputs: stored as N x M dgCMatrix (nodes rows, mechanisms cols)
  if (isTRUE(CellToNeighborhood)) {
    if (!"CellToNeighborhood" %in% names(niches_obj@assays)) {
      stop("Requested CellToNeighborhood=TRUE but NICHES output has no 'CellToNeighborhood' assay.")
    }
    m <- .get_mat(niches_obj, assay = "CellToNeighborhood", layer = niches_assay_layer)
    if (is.null(m)) stop("Could not fetch CellToNeighborhood matrix from NICHES output.")
    nm <- Matrix::t(m)  # edges: actually nodes x mechanisms after transpose (features x cells)
    nm <- as(nm, "dgCMatrix")
    # align node ordering
    node_ids <- rownames(nm)
    if (!is.null(node_ids)) {
      o <- match(eff$nodes$ids, node_ids)
      keep <- is.finite(o) & !is.na(o)
      nm2 <- nm[o[keep], , drop = FALSE]
      # if not all nodes present, you may want to pad; for now, require full coverage
      if (nrow(nm2) != length(eff$nodes$ids)) {
        stop("CellToNeighborhood nodes do not fully cover eff$nodes$ids (coverage mismatch).")
      }
      rownames(nm2) <- eff$nodes$ids
      nm <- nm2
    }
    out$CellToNeighborhood <- nm
    if (is.null(out$mechanisms)) out$mechanisms <- colnames(nm)
  }
  
  if (isTRUE(NeighborhoodToCell)) {
    if (!"NeighborhoodToCell" %in% names(niches_obj@assays)) {
      stop("Requested NeighborhoodToCell=TRUE but NICHES output has no 'NeighborhoodToCell' assay.")
    }
    m <- .get_mat(niches_obj, assay = "NeighborhoodToCell", layer = niches_assay_layer)
    if (is.null(m)) stop("Could not fetch NeighborhoodToCell matrix from NICHES output.")
    nm <- Matrix::t(m)
    nm <- as(nm, "dgCMatrix")
    node_ids <- rownames(nm)
    if (!is.null(node_ids)) {
      o <- match(eff$nodes$ids, node_ids)
      keep <- is.finite(o) & !is.na(o)
      nm2 <- nm[o[keep], , drop = FALSE]
      if (nrow(nm2) != length(eff$nodes$ids)) {
        stop("NeighborhoodToCell nodes do not fully cover eff$nodes$ids (coverage mismatch).")
      }
      rownames(nm2) <- eff$nodes$ids
      nm <- nm2
    }
    out$NeighborhoodToCell <- nm
    if (is.null(out$mechanisms)) out$mechanisms <- colnames(nm)
  }
  
  eff$niches <- out
  eff
}

# Optional: minimal validator
validate_effniches <- function(eff) {
  stopifnot(inherits(eff, "EffNICHES"))
  N <- length(eff$nodes$ids)
  stopifnot(nrow(eff$nodes$xy) == N, ncol(eff$nodes$xy) == 2)
  stopifnot(ncol(eff$expr$mat) == N)
  
  if (!is.null(eff$niches) && !is.null(eff$niches$CellToCellSpatial)) {
    e <- eff$niches$CellToCellSpatial
    stopifnot(nrow(e$ij) == nrow(e$w))
    stopifnot(ncol(e$w) == length(e$mechanisms))
  }
  TRUE
}

# fst-style serialization helpers
# Uses fst for node tables, qs for sparse matrices + ragged arrays.
save_effniches <- function(eff, dir, prefix = "effniches") {
  stopifnot(inherits(eff, "EffNICHES"))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  manifest <- list(
    version = eff$version,
    prefix = prefix,
    files = list()
  )
  
  # nodes: ids, xy, meta (fst)
  nodes_df <- data.frame(
    id = eff$nodes$ids,
    x = eff$nodes$xy[, 1],
    y = eff$nodes$xy[, 2],
    stringsAsFactors = FALSE
  )
  if (!is.null(eff$nodes$meta)) {
    nodes_df <- cbind(nodes_df, eff$nodes$meta)
  }
  f_nodes <- file.path(dir, paste0(prefix, "_nodes.fst"))
  fst::write_fst(nodes_df, f_nodes)
  manifest$files$nodes <- basename(f_nodes)
  
  # polygons (qs)
  if (!is.null(eff$nodes$polygons)) {
    f_poly_xy <- file.path(dir, paste0(prefix, "_poly_xy.qs"))
    f_poly_st <- file.path(dir, paste0(prefix, "_poly_start.qs"))
    qs::qsave(eff$nodes$polygons$xy, f_poly_xy)
    qs::qsave(eff$nodes$polygons$start, f_poly_st)
    manifest$files$polygons <- list(xy = basename(f_poly_xy), start = basename(f_poly_st))
  }
  
  # expression matrix (qs)
  f_expr <- file.path(dir, paste0(prefix, "_expr_mat.qs"))
  qs::qsave(eff$expr$mat, f_expr)
  manifest$files$expr_mat <- basename(f_expr)
  
  # expr metadata + sources (rds)
  f_info <- file.path(dir, paste0(prefix, "_info.rds"))
  info <- list(expr = eff$expr[setdiff(names(eff$expr), "mat")], sources = eff$sources)
  saveRDS(info, f_info)
  manifest$files$info <- basename(f_info)
  
  # niches (optional) (qs for matrices, rds for params)
  if (!is.null(eff$niches)) {
    # params
    f_nparams <- file.path(dir, paste0(prefix, "_niches_params.rds"))
    saveRDS(eff$niches$params, f_nparams)
    manifest$files$niches_params <- basename(f_nparams)
    
    # edges
    if (!is.null(eff$niches$CellToCellSpatial)) {
      f_ij <- file.path(dir, paste0(prefix, "_edges_ij.qs"))
      f_w  <- file.path(dir, paste0(prefix, "_edges_w.qs"))
      f_me <- file.path(dir, paste0(prefix, "_edges_mechanisms.qs"))
      qs::qsave(eff$niches$CellToCellSpatial$ij, f_ij)
      qs::qsave(eff$niches$CellToCellSpatial$w,  f_w)
      qs::qsave(eff$niches$CellToCellSpatial$mechanisms, f_me)
      manifest$files$CellToCellSpatial <- list(
        ij = basename(f_ij), w = basename(f_w), mechanisms = basename(f_me)
      )
    }
    
    # node-level
    if (!is.null(eff$niches$CellToNeighborhood)) {
      f <- file.path(dir, paste0(prefix, "_CellToNeighborhood.qs"))
      qs::qsave(eff$niches$CellToNeighborhood, f)
      manifest$files$CellToNeighborhood <- basename(f)
    }
    if (!is.null(eff$niches$NeighborhoodToCell)) {
      f <- file.path(dir, paste0(prefix, "_NeighborhoodToCell.qs"))
      qs::qsave(eff$niches$NeighborhoodToCell, f)
      manifest$files$NeighborhoodToCell <- basename(f)
    }
  }
  
  # manifest
  f_manifest <- file.path(dir, paste0(prefix, "_manifest.rds"))
  saveRDS(manifest, f_manifest)
  
  invisible(f_manifest)
}

load_effniches <- function(dir, prefix = "effniches") {
  f_manifest <- file.path(dir, paste0(prefix, "_manifest.rds"))
  man <- readRDS(f_manifest)
  
  nodes_df <- fst::read_fst(file.path(dir, man$files$nodes), as.data.table = FALSE)
  ids <- nodes_df$id
  xy <- as.matrix(nodes_df[, c("x", "y"), drop = FALSE])
  storage.mode(xy) <- "double"
  
  meta <- NULL
  meta_cols <- setdiff(colnames(nodes_df), c("id", "x", "y"))
  if (length(meta_cols) > 0) {
    meta <- nodes_df[, meta_cols, drop = FALSE]
    rownames(meta) <- ids
  }
  
  polygons <- NULL
  if (!is.null(man$files$polygons)) {
    polygons <- list(
      xy = qs::qread(file.path(dir, man$files$polygons$xy)),
      start = qs::qread(file.path(dir, man$files$polygons$start))
    )
  }
  
  expr_mat <- qs::qread(file.path(dir, man$files$expr_mat))
  info <- readRDS(file.path(dir, man$files$info))
  
  eff <- list(
    version = man$version,
    nodes = list(
      ids = ids,
      xy = xy,
      meta = meta,
      polygons = polygons,
      index = list(id_to_idx = setNames(seq_along(ids), ids))
    ),
    expr = c(info$expr, list(mat = expr_mat)),
    niches = NULL,
    sources = info$sources
  )
  
  # niches (optional)
  if (!is.null(man$files$niches_params)) {
    params <- readRDS(file.path(dir, man$files$niches_params))
    niches <- list(params = params, mechanisms = NULL)
    
    if (!is.null(man$files$CellToCellSpatial)) {
      ij <- qs::qread(file.path(dir, man$files$CellToCellSpatial$ij))
      w  <- qs::qread(file.path(dir, man$files$CellToCellSpatial$w))
      me <- qs::qread(file.path(dir, man$files$CellToCellSpatial$mechanisms))
      niches$CellToCellSpatial <- list(ij = ij, w = w, mechanisms = me, directed = TRUE)
      niches$mechanisms <- me
    }
    if (!is.null(man$files$CellToNeighborhood)) {
      niches$CellToNeighborhood <- qs::qread(file.path(dir, man$files$CellToNeighborhood))
      if (is.null(niches$mechanisms)) niches$mechanisms <- colnames(niches$CellToNeighborhood)
    }
    if (!is.null(man$files$NeighborhoodToCell)) {
      niches$NeighborhoodToCell <- qs::qread(file.path(dir, man$files$NeighborhoodToCell))
      if (is.null(niches$mechanisms)) niches$mechanisms <- colnames(niches$NeighborhoodToCell)
    }
    
    eff$niches <- niches
  }
  
  class(eff) <- c("EffNICHES", "list")
  eff
}

#### EffNICHES Add Polygon Data ####

# Allows for addition of segmentation polygon data after build_effniches
effniches_addpoly <- function(
    eff,
    obj,
    images = NULL,
    interactive = TRUE,
    selection = NULL,
    polygon_rows_per_cell_thresh = 1.5,
    allow_none = TRUE,
    overwrite = TRUE,
    verbose = TRUE
) {
  stopifnot(inherits(eff, "EffNICHES"))
  stopifnot(inherits(obj, "Seurat"))
  
  # -------------------------
  # local helpers
  # -------------------------
  .safe_num <- function(x) suppressWarnings(as.numeric(x))
  
  pick_xy_cols <- function(cn) {
    cn_l <- tolower(cn)
    
    pick_one <- function(cands) {
      hit <- intersect(cands, cn_l)
      if (length(hit) > 0) return(cn[match(hit[1], cn_l)])
      NULL
    }
    
    xcol <- pick_one(c("x", "row", "imagerow", "pxl_row_in_fullres", "pxl_row", "xcoord", "x_coord"))
    ycol <- pick_one(c("y", "col", "imagecol", "pxl_col_in_fullres", "pxl_col", "ycoord", "y_coord"))
    
    if (is.null(xcol) || is.null(ycol)) {
      if ("row" %in% cn && "col" %in% cn) {
        xcol <- "row"
        ycol <- "col"
      } else if (length(cn) >= 2) {
        if (is.null(xcol)) xcol <- cn[1]
        if (is.null(ycol)) ycol <- cn[min(2, length(cn))]
      }
    }
    
    list(xcol = xcol, ycol = ycol)
  }
  
  pick_cell_ids <- function(coords) {
    cn <- colnames(coords)
    if ("cell" %in% cn) return(as.character(coords$cell))
    rn <- rownames(coords)
    if (!is.null(rn) && length(rn) == nrow(coords)) return(as.character(rn))
    as.character(seq_len(nrow(coords)))
  }
  
  classify_coords_table <- function(df) {
    n_rows <- nrow(df)
    n_cells <- length(unique(df$cell))
    rows_per_cell <- n_rows / max(n_cells, 1)
    tab <- table(df$cell)
    frac_multi <- mean(tab > 1)
    
    is_polygon <- (rows_per_cell > polygon_rows_per_cell_thresh) ||
      (frac_multi >= 0.25 && rows_per_cell > 1.1)
    
    list(
      is_polygon = is_polygon,
      n_rows = n_rows,
      n_cells = n_cells,
      rows_per_cell = rows_per_cell,
      frac_multi = frac_multi
    )
  }
  
  polygons_to_ragged <- function(polygons_df, ids) {
    if (is.null(polygons_df) || nrow(polygons_df) == 0) return(NULL)
    
    stopifnot(all(c("cell", "x", "y") %in% colnames(polygons_df)))
    
    polygons_df$cell <- as.character(polygons_df$cell)
    keep <- polygons_df$cell %in% ids
    polygons_df <- polygons_df[keep, , drop = FALSE]
    if (nrow(polygons_df) == 0) return(NULL)
    
    id_to_idx <- setNames(seq_along(ids), ids)
    idx <- unname(id_to_idx[polygons_df$cell])
    ok <- is.finite(idx) & !is.na(idx)
    polygons_df <- polygons_df[ok, , drop = FALSE]
    idx <- idx[ok]
    
    o <- order(idx)
    idx <- idx[o]
    xy <- cbind(.safe_num(polygons_df$x[o]), .safe_num(polygons_df$y[o]))
    storage.mode(xy) <- "double"
    
    N <- length(ids)
    counts <- tabulate(idx, nbins = N)
    start <- integer(N + 1)
    start[1] <- 1L
    if (N > 0) start[2:(N + 1)] <- 1L + cumsum(counts)
    
    list(
      xy = xy,
      start = start
    )
  }
  
  fmt_num <- function(x, digits = 2) {
    if (length(x) == 0 || all(!is.finite(x))) return("NA")
    formatC(x, format = "f", digits = digits)
  }
  
  # -------------------------
  # discover image candidates
  # -------------------------
  imgs <- tryCatch(Seurat::Images(obj), error = function(e) character(0))
  if (length(imgs) == 0) {
    stop("No images found in Seurat object.")
  }
  
  candidates <- list()
  
  for (im in imgs) {
    coords <- tryCatch(Seurat::GetTissueCoordinates(obj, image = im), error = function(e) NULL)
    if (is.null(coords) || nrow(coords) == 0) next
    
    cn <- colnames(coords)
    xy <- pick_xy_cols(cn)
    if (is.null(xy$xcol) || is.null(xy$ycol)) next
    
    cell_id <- pick_cell_ids(coords)
    df <- data.frame(
      cell = as.character(cell_id),
      x = .safe_num(coords[[xy$xcol]]),
      y = .safe_num(coords[[xy$ycol]]),
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
    if (nrow(df) == 0) next
    
    cls <- classify_coords_table(df)
    
    candidates[[im]] <- list(
      image = im,
      xy_cols = c(xy$xcol, xy$ycol),
      is_polygon = isTRUE(cls$is_polygon),
      n_rows = cls$n_rows,
      n_cells = cls$n_cells,
      rows_per_cell = cls$rows_per_cell,
      frac_multi = cls$frac_multi
    )
  }
  
  if (length(candidates) == 0) {
    stop("No usable image coordinate tables detected.")
  }
  
  polygon_candidates <- candidates[vapply(candidates, function(x) isTRUE(x$is_polygon), logical(1))]
  if (length(polygon_candidates) == 0) {
    stop("No image entries with perimeter-vertex style coordinates were detected.")
  }
  
  # -------------------------
  # resolve selected images
  # -------------------------
  selected_images <- NULL
  
  if (!is.null(images)) {
    bad <- setdiff(images, names(candidates))
    if (length(bad) > 0) {
      stop("Requested images not found in object: ", paste(bad, collapse = ", "))
    }
    selected_images <- intersect(images, names(polygon_candidates))
    if (length(selected_images) == 0) {
      stop("None of the requested images appear to contain perimeter vertices.")
    }
    
  } else if (!is.null(selection)) {
    # programmatic selection support:
    # - numeric indices into polygon candidate list
    # - "all" or "all_polygons"
    if (is.character(selection) && length(selection) == 1 && selection %in% c("all", "all_polygons")) {
      selected_images <- names(polygon_candidates)
    } else {
      idx <- as.integer(selection)
      keys <- names(polygon_candidates)
      if (any(!is.finite(idx)) || any(idx < 1) || any(idx > length(keys))) {
        stop("Invalid numeric selection for polygon images.")
      }
      selected_images <- keys[idx]
    }
    
  } else if (isTRUE(interactive)) {
    cat("\n==============================\n")
    cat("EffNICHES Polygon Image Picker\n")
    cat("==============================\n\n")
    cat("Polygon-capable image entries:\n")
    
    pkeys <- names(polygon_candidates)
    cat("A) All images with perimeter vertices\n")
    for (i in seq_along(pkeys)) {
      x <- polygon_candidates[[pkeys[i]]]
      cat(sprintf(
        "%d) %s\n   - rows: %s | cells: %s | approx vertices/cell: %s | frac multi-row cells: %s\n",
        i,
        pkeys[i],
        format(x$n_rows, big.mark = ","),
        format(x$n_cells, big.mark = ","),
        fmt_num(x$rows_per_cell, 2),
        fmt_num(x$frac_multi, 3)
      ))
    }
    if (isTRUE(allow_none)) cat("0) Cancel / return unchanged EffNICHES\n")
    
    ans <- readline(prompt = "\nEnter option number (or A for all): ")
    ans <- trimws(ans)
    
    if (identical(toupper(ans), "A")) {
      selected_images <- pkeys
    } else {
      idx <- suppressWarnings(as.integer(ans))
      if (!is.finite(idx)) stop("Invalid selection.")
      if (idx == 0 && isTRUE(allow_none)) return(eff)
      if (idx < 1 || idx > length(pkeys)) stop("Selection out of range.")
      selected_images <- pkeys[idx]
    }
    
  } else {
    stop("No image selection provided. Use images=, selection=, or interactive=TRUE.")
  }
  
  # -------------------------
  # gather polygons from selected images
  # -------------------------
  ids <- eff$nodes$ids
  poly_list <- vector("list", length(selected_images))
  names(poly_list) <- selected_images
  
  for (im in selected_images) {
    coords <- Seurat::GetTissueCoordinates(obj, image = im)
    cnd <- polygon_candidates[[im]]
    
    cell_id <- pick_cell_ids(coords)
    df <- data.frame(
      image = im,
      cell = as.character(cell_id),
      x = .safe_num(coords[[cnd$xy_cols[1]]]),
      y = .safe_num(coords[[cnd$xy_cols[2]]]),
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
    
    # keep only cells already in eff
    df <- df[df$cell %in% ids, , drop = FALSE]
    poly_list[[im]] <- df
  }
  
  poly_df <- do.call(rbind, poly_list)
  
  if (is.null(poly_df) || nrow(poly_df) == 0) {
    warning("No matching polygon vertices found for eff$nodes$ids in the selected image(s).")
    return(eff)
  }
  
  # -------------------------
  # check for cells appearing in multiple selected images
  # -------------------------
  cell_image_tab <- table(poly_df$cell, poly_df$image)
  cells_in_multiple_images <- rownames(cell_image_tab)[rowSums(cell_image_tab > 0) > 1]
  
  if (length(cells_in_multiple_images) > 0) {
    stop(
      "Some eff node IDs appear in more than one selected image, which would create ambiguous polygon mappings. ",
      "Example cell(s): ", paste(head(cells_in_multiple_images, 5), collapse = ", ")
    )
  }
  
  # remove image column before packing
  poly_df_out <- poly_df[, c("cell", "x", "y"), drop = FALSE]
  poly_ragged <- polygons_to_ragged(poly_df_out, ids = ids)
  
  if (is.null(poly_ragged)) {
    warning("Polygon extraction completed but no polygons could be mapped onto eff$nodes$ids.")
    return(eff)
  }
  
  # -------------------------
  # assign polygons only; do not alter centroids
  # -------------------------
  if (!isTRUE(overwrite) && !is.null(eff$nodes$polygons)) {
    stop("eff already has polygons and overwrite=FALSE.")
  }
  
  eff$nodes$polygons <- poly_ragged
  
  if (isTRUE(verbose)) {
    matched_cells <- sum((poly_ragged$start[-1] - poly_ragged$start[-length(poly_ragged$start)]) > 0)
    total_vertices <- nrow(poly_ragged$xy)
    
    message("Added polygons from ", length(selected_images), " image(s).")
    message("Eff nodes with polygons: ", format(matched_cells, big.mark = ","),
            " / ", format(length(ids), big.mark = ","))
    message("Total stored vertices: ", format(total_vertices, big.mark = ","))
    message("Centroid coordinates in eff$nodes$xy were left unchanged.")
  }
  
  eff
}

#### Import NICHES Data ####

# Helper: robust assay/layer fetcher
.get_mat <- function(obj, assay, layer) {
  a <- obj[[assay]]
  
  m <- tryCatch(SeuratObject::LayerData(a, layer = layer), error = function(e) NULL)
  if (!is.null(m)) return(m)
  
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = layer), error = function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  NULL
}

# Helper: recover sender / receiver barcodes from ctc edge names
.get_ctc_edge_endpoints <- function(ctc, edge_delim = "—") {
  stopifnot(inherits(ctc, "Seurat"))
  
  edge_names <- Seurat::Cells(ctc)
  if (is.null(edge_names) || length(edge_names) == 0) {
    stop("ctc has no Seurat cell names / edge names.")
  }
  
  parts <- strsplit(edge_names, split = edge_delim, fixed = TRUE)
  bad <- vapply(parts, length, integer(1)) != 2L
  if (any(bad)) {
    stop(
      "Some ctc edge names did not split into exactly 2 parts using delimiter '",
      edge_delim, "'. Example bad edge: ", edge_names[which(bad)[1]]
    )
  }
  
  sender_ids   <- vapply(parts, `[[`, character(1), 1)
  receiver_ids <- vapply(parts, `[[`, character(1), 2)
  
  list(
    edge_names = edge_names,
    sender_ids = sender_ids,
    receiver_ids = receiver_ids
  )
}

# Helper: convert ctc into the EffNICHES edge store format
.extract_edges_from_ctc_E_by_M <- function(
    ctc,
    eff_ids,
    ctc_assay = "CellToCellSpatial",
    ctc_layer = c("counts", "data"),
    edge_delim = "—",
    verbose = TRUE
) {
  ctc_layer <- match.arg(ctc_layer)
  
  if (!inherits(ctc, "Seurat")) {
    stop("'ctc' must be a Seurat object.")
  }
  if (!ctc_assay %in% names(ctc@assays)) {
    stop("Assay '", ctc_assay, "' not found in ctc.")
  }
  
  ep <- .get_ctc_edge_endpoints(ctc, edge_delim = edge_delim)
  
  edge_names   <- ep$edge_names
  sender_ids   <- ep$sender_ids
  receiver_ids <- ep$receiver_ids
  
  id_to_idx <- setNames(seq_along(eff_ids), eff_ids)
  s_idx <- unname(id_to_idx[sender_ids])
  r_idx <- unname(id_to_idx[receiver_ids])
  
  keep <- is.finite(s_idx) & !is.na(s_idx) & is.finite(r_idx) & !is.na(r_idx)
  
  if (isTRUE(verbose)) {
    message("Total ctc edges: ", format(length(edge_names), big.mark = ","))
    message("Retained edges:  ", format(sum(keep), big.mark = ","))
    message("Dropped edges:   ", format(sum(!keep), big.mark = ","))
  }
  
  if (!any(keep)) {
    stop("No ctc edges remained after restricting to eff$nodes$ids.")
  }
  
  m <- .get_mat(ctc, assay = ctc_assay, layer = ctc_layer)
  if (is.null(m)) {
    stop("Could not fetch ctc matrix for assay='", ctc_assay, "', layer='", ctc_layer, "'.")
  }
  
  if (is.null(rownames(m))) {
    stop("ctc matrix is missing mechanism names (rownames).")
  }
  if (is.null(colnames(m))) {
    stop("ctc matrix is missing edge names (colnames).")
  }
  
  if (!identical(colnames(m), edge_names)) {
    if (all(edge_names %in% colnames(m))) {
      m <- m[, edge_names, drop = FALSE]
    } else {
      stop("ctc matrix column names do not align with Cells(ctc).")
    }
  }
  
  # mechanism-by-edge -> subset edges -> transpose to edge-by-mechanism
  m_sub <- m[, keep, drop = FALSE]
  w_E_by_M <- Matrix::t(m_sub)
  if (!inherits(w_E_by_M, "dgCMatrix")) {
    w_E_by_M <- methods::as(w_E_by_M, "dgCMatrix")
  }
  
  mechanisms <- rownames(m)
  colnames(w_E_by_M) <- mechanisms
  
  ij <- cbind(
    sender   = as.integer(s_idx[keep]),
    receiver = as.integer(r_idx[keep])
  )
  
  list(
    ij = ij,
    w = w_E_by_M,
    mechanisms = mechanisms,
    edge_names = edge_names[keep],
    sender_ids = sender_ids[keep],
    receiver_ids = receiver_ids[keep]
  )
}

# Public API: import precomputed CellToCellSpatial from ctc Seurat into an EffNICHES object
import_effniches_ctc <- function(
    eff,
    ctc,
    ctc_assay = "CellToCellSpatial",
    ctc_layer = c("counts", "data"),
    edge_delim = "—",
    replace_existing = TRUE,
    keep_edge_names = FALSE,
    verbose = TRUE
) {
  ctc_layer <- match.arg(ctc_layer)
  
  stopifnot(inherits(eff, "EffNICHES"))
  stopifnot(inherits(ctc, "Seurat"))
  
  if (!isTRUE(replace_existing) && !is.null(eff$niches)) {
    stop("eff$niches is already populated and replace_existing=FALSE.")
  }
  
  edges <- .extract_edges_from_ctc_E_by_M(
    ctc = ctc,
    eff_ids = eff$nodes$ids,
    ctc_assay = ctc_assay,
    ctc_layer = ctc_layer,
    edge_delim = edge_delim,
    verbose = verbose
  )
  
  cell_to_cell <- list(
    ij = edges$ij,
    w = edges$w,
    mechanisms = edges$mechanisms,
    directed = TRUE
  )
  
  if (isTRUE(keep_edge_names)) {
    cell_to_cell$edge_names <- edges$edge_names
  }
  
  eff$niches <- list(
    params = list(
      imported_from_existing_ctc = TRUE,
      ctc_assay = ctc_assay,
      ctc_layer = ctc_layer,
      edge_delim = edge_delim,
      expression_source = eff$sources$expr_spec,
      spatial_source = eff$sources$spatial_spec
    ),
    mechanisms = edges$mechanisms,
    CellToCellSpatial = cell_to_cell,
    CellToNeighborhood = NULL,
    NeighborhoodToCell = NULL
  )
  
  eff
}

#### Subset NICHES ####

subsetEffNICHES <- function(eff,
                            cells = NULL,
                            subset,
                            invert = FALSE,
                            edge_delim = "\u2014",
                            drop_empty_mechanisms = FALSE,
                            verbose = TRUE) {
  if (!inherits(eff, "EffNICHES") && !is.list(eff)) {
    stop("'eff' must be an EffNICHES object (or compatible list).")
  }
  if (is.null(eff$nodes) || is.null(eff$nodes$ids)) {
    stop("eff$nodes$ids is missing.")
  }
  if (is.null(eff$nodes$xy)) {
    stop("eff$nodes$xy is missing.")
  }
  if (is.null(eff$expr) || is.null(eff$expr$mat)) {
    stop("eff$expr$mat is missing.")
  }
  if (is.null(eff$niches) || is.null(eff$niches$CellToCellSpatial)) {
    stop("eff$niches$CellToCellSpatial is missing.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  
  .safe_num <- function(x) suppressWarnings(as.numeric(x))
  
  .get_node_meta <- function(eff, node_ids) {
    if (!is.null(eff$nodes$meta)) {
      meta <- eff$nodes$meta
      if (!is.data.frame(meta)) {
        stop("eff$nodes$meta exists but is not a data.frame.")
      }
      if (is.null(rownames(meta))) {
        stop("eff$nodes$meta exists but has no rownames.")
      }
      if (!all(node_ids %in% rownames(meta))) {
        stop("Not all eff$nodes$ids are present in rownames(eff$nodes$meta).")
      }
      return(meta[node_ids, , drop = FALSE])
    }
    data.frame(row.names = node_ids)
  }
  
  .parse_edges_from_names <- function(edge_ids, delim) {
    edge_ids <- as.character(edge_ids)
    sender <- sub(paste0(delim, ".*$"), "", edge_ids, perl = TRUE)
    receiver <- sub(paste0("^.*", delim), "", edge_ids, perl = TRUE)
    data.frame(
      edge_id = edge_ids,
      sender = sender,
      receiver = receiver,
      stringsAsFactors = FALSE,
      row.names = edge_ids
    )
  }
  
  .get_edge_df <- function(ctc, edge_ids, delim) {
    candidate_names <- c("edges", "edge_df", "meta", "pairs")
    
    for (nm in candidate_names) {
      obj <- ctc[[nm]]
      if (!is.data.frame(obj) || nrow(obj) == 0L) next
      
      nms <- colnames(obj)
      sender_col <- intersect(nms, c("sender", "from", "source", "cell_sender", "sending_cell"))
      receiver_col <- intersect(nms, c("receiver", "to", "target", "cell_receiver", "receiving_cell"))
      
      if (length(sender_col) < 1L || length(receiver_col) < 1L) next
      
      sender_col <- sender_col[1]
      receiver_col <- receiver_col[1]
      
      out <- obj
      if (is.null(rownames(out))) {
        if (nrow(out) != length(edge_ids)) {
          stop("Found explicit edge table but could not align it to edge IDs.")
        }
        rownames(out) <- edge_ids
      }
      if (!all(edge_ids %in% rownames(out))) {
        stop("Explicit edge table does not contain all edge IDs.")
      }
      
      out <- out[edge_ids, , drop = FALSE]
      out$sender <- as.character(out[[sender_col]])
      out$receiver <- as.character(out[[receiver_col]])
      out$edge_id <- rownames(out)
      return(out)
    }
    
    .parse_edges_from_names(edge_ids, delim = delim)
  }
  
  .rebuild_polygon_store <- function(polys, keep_idx) {
    if (is.null(polys)) return(NULL)
    if (!is.list(polys)) return(polys)
    if (is.null(polys$start) || is.null(polys$xy)) return(polys)
    
    old_start <- as.integer(polys$start)
    old_xy <- polys$xy
    
    N_old <- length(old_start) - 1L
    if (N_old < 0L) {
      stop("Invalid polygon store: length(start) must be at least 1.")
    }
    if (length(keep_idx) == 0L) {
      return(list(
        xy = matrix(numeric(0), ncol = 2),
        start = 1L
      ))
    }
    if (any(keep_idx < 1L | keep_idx > N_old)) {
      stop("Polygon keep_idx out of bounds.")
    }
    
    counts_new <- integer(length(keep_idx))
    vertex_ranges <- vector("list", length(keep_idx))
    
    for (j in seq_along(keep_idx)) {
      i <- keep_idx[j]
      a <- old_start[i]
      b <- old_start[i + 1L] - 1L
      
      if (is.na(a) || is.na(b) || b < a) {
        vertex_ranges[[j]] <- integer(0)
        counts_new[j] <- 0L
      } else {
        idx <- seq.int(a, b)
        vertex_ranges[[j]] <- idx
        counts_new[j] <- length(idx)
      }
    }
    
    vertex_idx <- unlist(vertex_ranges, use.names = FALSE)
    
    xy_new <- old_xy[vertex_idx, , drop = FALSE]
    if (!is.matrix(xy_new)) {
      xy_new <- as.matrix(xy_new)
    }
    storage.mode(xy_new) <- "double"
    
    start_new <- integer(length(keep_idx) + 1L)
    start_new[1] <- 1L
    start_new[-1] <- 1L + cumsum(counts_new)
    
    out <- list(
      xy = xy_new,
      start = start_new
    )
    
    # Preserve any additional polygon-level components only when their alignment
    # is unambiguous.
    extra_names <- setdiff(names(polys), c("xy", "start"))
    for (nm in extra_names) {
      x <- polys[[nm]]
      
      if (is.atomic(x) && !is.list(x) && length(x) == N_old) {
        out[[nm]] <- x[keep_idx]
      } else if (is.atomic(x) && !is.list(x) && length(x) == nrow(old_xy)) {
        out[[nm]] <- x[vertex_idx]
      } else if (is.data.frame(x) && nrow(x) == N_old) {
        out[[nm]] <- x[keep_idx, , drop = FALSE]
      } else if ((is.matrix(x) || inherits(x, "Matrix")) && nrow(x) == N_old) {
        out[[nm]] <- x[keep_idx, , drop = FALSE]
      } else if (is.data.frame(x) && nrow(x) == nrow(old_xy)) {
        out[[nm]] <- x[vertex_idx, , drop = FALSE]
      } else if ((is.matrix(x) || inherits(x, "Matrix")) && nrow(x) == nrow(old_xy)) {
        out[[nm]] <- x[vertex_idx, , drop = FALSE]
      } else {
        out[[nm]] <- x
      }
    }
    
    out
  }
  
  .drop_empty_mechs <- function(ctc) {
    if (is.null(ctc$w)) return(ctc)
    w <- ctc$w
    if (!(is.matrix(w) || inherits(w, "Matrix"))) return(ctc)
    
    keep_mech <- Matrix::colSums(abs(w) != 0) > 0
    keep_mech[is.na(keep_mech)] <- FALSE
    
    ctc$w <- w[, keep_mech, drop = FALSE]
    
    if (!is.null(ctc$mechanisms) && length(ctc$mechanisms) == ncol(w)) {
      ctc$mechanisms <- ctc$mechanisms[keep_mech]
    }
    
    ctc
  }
  
  all_nodes <- as.character(eff$nodes$ids)
  N_old <- length(all_nodes)
  keep_nodes <- all_nodes
  
  meta <- .get_node_meta(eff, all_nodes)
  
  if (!missing(subset)) {
    subset_expr <- substitute(subset)
    r <- eval(subset_expr, envir = meta, enclos = parent.frame())
    if (!is.logical(r) || length(r) != nrow(meta)) {
      stop("'subset' must evaluate to a logical vector of length equal to the number of nodes.")
    }
    r[is.na(r)] <- FALSE
    keep_nodes <- rownames(meta)[r]
  }
  
  if (!is.null(cells)) {
    cells <- intersect(as.character(cells), all_nodes)
    keep_nodes <- intersect(keep_nodes, cells)
  }
  
  if (isTRUE(invert)) {
    keep_nodes <- setdiff(all_nodes, keep_nodes)
  }
  
  keep_nodes <- all_nodes[all_nodes %in% keep_nodes]
  
  if (length(keep_nodes) == 0L) {
    stop("Subsetting removed all nodes.")
  }
  
  keep_idx <- match(keep_nodes, all_nodes)
  if (any(is.na(keep_idx))) {
    stop("Internal error: failed to map keep_nodes onto all_nodes.")
  }
  
  new_id_to_idx <- setNames(seq_along(keep_nodes), keep_nodes)
  
  out <- eff
  
  ## -------------------------
  ## nodes
  ## -------------------------
  out$nodes$ids <- keep_nodes
  
  if (!is.null(eff$nodes$xy)) {
    xy <- eff$nodes$xy
    if (!(is.matrix(xy) || is.data.frame(xy))) {
      stop("eff$nodes$xy must be a matrix or data.frame.")
    }
    
    if (!is.null(rownames(xy)) && all(all_nodes %in% rownames(xy))) {
      out$nodes$xy <- xy[keep_nodes, , drop = FALSE]
    } else {
      if (nrow(xy) != N_old) {
        stop("eff$nodes$xy does not align to eff$nodes$ids.")
      }
      out$nodes$xy <- xy[keep_idx, , drop = FALSE]
      rownames(out$nodes$xy) <- keep_nodes
    }
  }
  
  if (!is.null(eff$nodes$meta)) {
    out$nodes$meta <- meta[keep_nodes, , drop = FALSE]
  }
  
  if (!is.null(eff$nodes$polygons)) {
    out$nodes$polygons <- .rebuild_polygon_store(eff$nodes$polygons, keep_idx = keep_idx)
  }
  
  out$nodes$index <- list(id_to_idx = new_id_to_idx)
  
  if (!is.null(out$N) && length(out$N) == 1L) {
    out$N <- length(keep_nodes)
  }
  if (!is.null(out$nodes$N) && length(out$nodes$N) == 1L) {
    out$nodes$N <- length(keep_nodes)
  }
  
  ## -------------------------
  ## expression
  ## -------------------------
  if (!is.null(eff$expr)) {
    out$expr <- eff$expr
    
    if (!is.null(eff$expr$mat)) {
      x <- eff$expr$mat
      if (!(is.matrix(x) || inherits(x, "Matrix"))) {
        stop("eff$expr$mat must be a matrix-like object.")
      }
      
      if (!is.null(colnames(x)) && all(keep_nodes %in% colnames(x))) {
        out$expr$mat <- x[, keep_nodes, drop = FALSE]
      } else if (ncol(x) == N_old) {
        out$expr$mat <- x[, keep_idx, drop = FALSE]
        colnames(out$expr$mat) <- keep_nodes
      } else {
        stop("eff$expr$mat does not align to eff$nodes$ids.")
      }
    }
    
    if (!is.null(eff$expr$cell_ids) &&
        is.atomic(eff$expr$cell_ids) &&
        length(eff$expr$cell_ids) == N_old) {
      out$expr$cell_ids <- keep_nodes
    }
  }
  
  ## -------------------------
  ## niches
  ## -------------------------
  out$niches <- eff$niches
  ctc <- eff$niches$CellToCellSpatial
  
  edge_ids <- NULL
  if (!is.null(ctc$w) && (is.matrix(ctc$w) || inherits(ctc$w, "Matrix"))) {
    edge_ids <- rownames(ctc$w)
  }
  if (is.null(edge_ids) || length(edge_ids) == 0L) {
    edge_ids <- NULL
    for (nm in c("edges", "edge_df", "meta")) {
      obj <- ctc[[nm]]
      if (is.data.frame(obj) && !is.null(rownames(obj)) && nrow(obj) > 0L) {
        edge_ids <- rownames(obj)
        break
      }
    }
  }
  if (is.null(edge_ids) || length(edge_ids) == 0L) {
    stop("Could not determine edge IDs in eff$niches$CellToCellSpatial.")
  }
  
  edge_ids <- as.character(edge_ids)
  edge_df <- .get_edge_df(ctc, edge_ids = edge_ids, delim = edge_delim)
  
  keep_edge_mask <- edge_df$sender %in% keep_nodes & edge_df$receiver %in% keep_nodes
  keep_edge_mask[is.na(keep_edge_mask)] <- FALSE
  
  keep_edges <- edge_ids[keep_edge_mask]
  keep_edge_idx <- which(keep_edge_mask)
  
  ctc_out <- ctc
  
  ## edge x mechanism matrix
  if (!is.null(ctc$w)) {
    w <- ctc$w
    if (!(is.matrix(w) || inherits(w, "Matrix"))) {
      stop("eff$niches$CellToCellSpatial$w must be matrix-like.")
    }
    ctc_out$w <- w[keep_edge_idx, , drop = FALSE]
    if (!is.null(rownames(ctc_out$w))) {
      rownames(ctc_out$w) <- keep_edges
    }
  }
  
  ## remap ij to new node indexing
  if (!is.null(ctc$ij)) {
    ij <- ctc$ij
    if (!is.matrix(ij) && !is.data.frame(ij)) {
      stop("eff$niches$CellToCellSpatial$ij must be matrix/data.frame.")
    }
    if (nrow(ij) != length(edge_ids)) {
      stop("nrow(ctc$ij) must equal number of edges.")
    }
    
    ij_sub <- ij[keep_edge_idx, , drop = FALSE]
    if (ncol(ij_sub) < 2L) {
      stop("ctc$ij must have at least two columns (sender, receiver).")
    }
    
    sender_old_idx <- as.integer(ij_sub[, 1])
    receiver_old_idx <- as.integer(ij_sub[, 2])
    
    sender_ids <- all_nodes[sender_old_idx]
    receiver_ids <- all_nodes[receiver_old_idx]
    
    sender_new_idx <- unname(new_id_to_idx[sender_ids])
    receiver_new_idx <- unname(new_id_to_idx[receiver_ids])
    
    if (any(is.na(sender_new_idx) | is.na(receiver_new_idx))) {
      stop("Internal error remapping ctc$ij to subset node indices.")
    }
    
    ij_new <- cbind(
      sender = as.integer(sender_new_idx),
      receiver = as.integer(receiver_new_idx)
    )
    ctc_out$ij <- ij_new
  }
  
  ## edge metadata-like components
  for (nm in c("edges", "edge_df", "meta", "pairs")) {
    x <- ctc[[nm]]
    if (is.data.frame(x) && nrow(x) == length(edge_ids)) {
      ctc_out[[nm]] <- x[keep_edge_idx, , drop = FALSE]
      if (!is.null(rownames(ctc_out[[nm]]))) {
        rownames(ctc_out[[nm]]) <- keep_edges
      }
    }
  }
  
  if (!is.null(ctc$ids) && is.atomic(ctc$ids) && length(ctc$ids) == length(edge_ids)) {
    ctc_out$ids <- keep_edges
  }
  if (!is.null(ctc$edge_names) && is.atomic(ctc$edge_names) && length(ctc$edge_names) == length(edge_ids)) {
    ctc_out$edge_names <- ctc$edge_names[keep_edge_idx]
  }
  if (!is.null(ctc$N) && length(ctc$N) == 1L) {
    ctc_out$N <- length(keep_edges)
  }
  
  if (isTRUE(drop_empty_mechanisms)) {
    ctc_out <- .drop_empty_mechs(ctc_out)
  }
  
  out$niches$CellToCellSpatial <- ctc_out
  
  ## node-level niche matrices, if present
  for (nm in c("CellToNeighborhood", "NeighborhoodToCell")) {
    x <- eff$niches[[nm]]
    if (is.null(x)) next
    
    if (is.matrix(x) || inherits(x, "Matrix")) {
      if (!is.null(rownames(x)) && all(keep_nodes %in% rownames(x))) {
        out$niches[[nm]] <- x[keep_nodes, , drop = FALSE]
      } else if (nrow(x) == N_old) {
        out$niches[[nm]] <- x[keep_idx, , drop = FALSE]
        rownames(out$niches[[nm]]) <- keep_nodes
      } else {
        out$niches[[nm]] <- x
      }
    } else {
      out$niches[[nm]] <- x
    }
  }
  
  if (verbose) {
    message(
      "subsetEffNICHES complete: kept ",
      length(keep_nodes), " / ", length(all_nodes), " nodes and ",
      length(keep_edges), " / ", length(edge_ids), " edges."
    )
  }
  
  out
}
