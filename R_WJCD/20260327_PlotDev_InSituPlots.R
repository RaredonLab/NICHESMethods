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

#### Helper Functions ####

## Reusable LR helper functions

# Normalize common dash variants so user input can match stored
# mechanism names more reliably.
norm_lr_dash <- function(s) {
  s <- gsub("\u2014", "-", s)
  s <- gsub("\u2013", "-", s)
  s <- gsub("\u2212", "-", s)
  s
}

# Resolve one or more LR mechanism specifications against a
# mechanism name vector. Supports numeric indices and character
# names / partial matches.
resolveLRindices <- function(lr, mech_names, deduplicate = TRUE) {
  if (length(lr) < 1) {
    stop("'lr' must contain at least one mechanism name or index.")
  }
  if (is.null(mech_names) || !is.character(mech_names) || length(mech_names) < 1) {
    stop("'mech_names' must be a non-empty character vector.")
  }
  
  mn_n <- norm_lr_dash(mech_names)
  out <- integer(0)
  
  for (one_lr in lr) {
    if (is.numeric(one_lr) && length(one_lr) == 1) {
      ii <- as.integer(one_lr)
      if (!is.finite(ii) || ii < 1L || ii > length(mech_names)) {
        stop("Numeric 'lr' values must be between 1 and ", length(mech_names), ".")
      }
      out <- c(out, ii)
      next
    }
    
    if (!is.character(one_lr) || length(one_lr) != 1) {
      stop("'lr' must be a character vector of mechanism names and/or numeric indices.")
    }
    
    lr_n <- norm_lr_dash(one_lr)
    
    hit <- which(mn_n == lr_n)
    if (length(hit) == 1) {
      out <- c(out, hit)
      next
    }
    if (length(hit) > 1) {
      stop("Could not uniquely resolve lr = '", one_lr, "' by exact match.")
    }
    
    hit <- which(grepl(lr_n, mn_n, fixed = TRUE))
    if (length(hit) == 1) {
      out <- c(out, hit)
      next
    }
    if (length(hit) > 1) {
      stop("Could not uniquely resolve lr = '", one_lr, "'; multiple partial matches found.")
    }
    
    stop("Could not resolve lr = '", one_lr, "'.")
  }
  
  if (deduplicate) out <- unique(out)
  if (length(out) < 1) stop("No LR mechanisms resolved.")
  out
}

# Normalize an edges x mechanisms matrix column-wise so different
# mechanisms can be compared or combined on more equal footing.
normalizeLRcolumns <- function(mat,
                               method = c("none", "max", "percentile", "zscore")) {
  method <- match.arg(method)
  
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  
  if (ncol(mat) == 0L) return(mat)
  if (method == "none") return(mat)
  
  out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  colnames(out) <- colnames(mat)
  
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]
    x[!is.finite(x)] <- NA_real_
    
    if (all(is.na(x))) {
      out[, j] <- 0
      next
    }
    
    if (method == "max") {
      mx <- suppressWarnings(max(x, na.rm = TRUE))
      if (!is.finite(mx) || mx <= 0) {
        out[, j] <- 0
      } else {
        out[, j] <- x / mx
      }
    }
    
    if (method == "percentile") {
      ok <- is.finite(x)
      y <- numeric(length(x))
      if (any(ok)) {
        r <- rank(x[ok], ties.method = "average")
        n_ok <- sum(ok)
        if (n_ok <= 1L) {
          y[ok] <- 1
        } else {
          y[ok] <- (r - 1) / (n_ok - 1)
        }
      }
      out[, j] <- y
    }
    
    if (method == "zscore") {
      mu <- suppressWarnings(mean(x, na.rm = TRUE))
      s  <- suppressWarnings(stats::sd(x, na.rm = TRUE))
      if (!is.finite(s) || s == 0) {
        out[, j] <- 0
      } else {
        out[, j] <- (x - mu) / s
      }
    }
  }
  
  out[!is.finite(out)] <- 0
  out
}

# Aggregate an edges x mechanisms matrix into one score per edge.
aggregateLRcolumns <- function(mat,
                               method = c("sum", "mean", "max")) {
  method <- match.arg(method)
  
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  
  if (ncol(mat) == 0L) stop("'mat' must have at least one column.")
  if (ncol(mat) == 1L) return(as.numeric(mat[, 1]))
  
  if (method == "sum")  return(rowSums(mat, na.rm = TRUE))
  if (method == "mean") return(rowMeans(mat, na.rm = TRUE))
  if (method == "max")  return(apply(mat, 1, max, na.rm = TRUE))
}

## LR edge normalization and aggregation helper

# High-level engine:
# Resolve LR mechanisms, extract edge weights, normalize selected
# columns, aggregate them, and optionally transform the final
# edge-level signal.
computeEdgeLRsignal <- function(eff,
                                lr,
                                lr_normalization = c("none", "max", "percentile", "zscore"),
                                lr_aggregation   = c("sum", "mean", "max"),
                                deduplicate_lr   = TRUE,
                                transform = identity,
                                warn_many_lr = TRUE,
                                max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  
  if (is.null(eff$niches$CellToCellSpatial)) {
    stop("eff$niches$CellToCellSpatial is missing.")
  }
  
  c2c <- eff$niches$CellToCellSpatial
  
  if (is.null(c2c$ij) || !is.matrix(c2c$ij) || ncol(c2c$ij) != 2) {
    stop("CellToCellSpatial$ij must be a matrix with 2 columns (node indices).")
  }
  if (is.null(c2c$w) || !inherits(c2c$w, c("dgCMatrix", "Matrix", "matrix"))) {
    stop("CellToCellSpatial$w must be a matrix-like object (edges x mechanisms).")
  }
  if (is.null(c2c$mechanisms) || !is.character(c2c$mechanisms)) {
    stop("CellToCellSpatial$mechanisms must be a character vector.")
  }
  if (nrow(c2c$w) != nrow(c2c$ij)) {
    stop("Row mismatch: nrow(CellToCellSpatial$w) must equal nrow(CellToCellSpatial$ij).")
  }
  if (ncol(c2c$w) != length(c2c$mechanisms)) {
    stop("Column mismatch: ncol(CellToCellSpatial$w) must equal length(c2c$mechanisms).")
  }
  if (nrow(c2c$ij) > max_edges) {
    stop("Edge count (", nrow(c2c$ij), ") exceeds max_edges (", max_edges, ").")
  }
  
  idx <- resolveLRindices(
    lr = lr,
    mech_names = c2c$mechanisms,
    deduplicate = deduplicate_lr
  )
  
  selected_mechanisms <- c2c$mechanisms[idx]
  
  if (isTRUE(warn_many_lr) && length(idx) > 6L) {
    warning(
      "Aggregating ", length(idx),
      " LR mechanisms may obscure spatial structure and mechanistic interpretation."
    )
  }
  
  if (identical(lr_normalization, "zscore") && lr_aggregation %in% c("sum", "mean")) {
    warning(
      "zscore normalization can produce negative values; aggregated edge scores may reflect cancellation rather than total signaling."
    )
  }
  
  w_sub <- as.matrix(c2c$w[, idx, drop = FALSE])
  storage.mode(w_sub) <- "double"
  colnames(w_sub) <- selected_mechanisms
  
  w_norm <- normalizeLRcolumns(w_sub, method = lr_normalization)
  w_edge_raw <- aggregateLRcolumns(w_norm, method = lr_aggregation)
  
  if (length(w_edge_raw) != nrow(c2c$ij)) {
    stop("Aggregated edge signal length mismatch.")
  }
  
  resolve_transform <- function(transform) {
    if (is.function(transform)) {
      return(transform)
    }
    
    if (is.character(transform) && length(transform) == 1L && !is.na(transform)) {
      key <- tolower(trimws(transform))
      
      alias_map <- c(
        "identity" = "identity",
        "none"     = "identity",
        "linear"   = "identity",
        "log1p"    = "log1p",
        "logp1"    = "log1p",
        "log"      = "log",
        "sqrt"     = "sqrt",
        "square"   = "square"
      )
      
      if (!key %in% names(alias_map)) {
        stop(
          "Unknown transform string: '", transform, "'. ",
          "Use a function or one of: ",
          paste(sprintf("'%s'", unique(names(alias_map))), collapse = ", "),
          "."
        )
      }
      
      key <- alias_map[[key]]
      
      fn <- switch(
        key,
        "identity" = identity,
        "log1p"    = base::log1p,
        "log"      = base::log,
        "sqrt"     = base::sqrt,
        "square"   = function(x) x^2
      )
      
      return(fn)
    }
    
    stop(
      "'transform' must be either a function or a single character string ",
      "(for example identity, 'log1p', 'sqrt')."
    )
  }
  
  transform_fun <- resolve_transform(transform)
  
  w_edge <- tryCatch(
    transform_fun(w_edge_raw),
    error = function(e) {
      stop("transform failed on aggregated edge signal: ", conditionMessage(e))
    }
  )
  
  if (is.null(w_edge)) {
    stop("transform returned NULL; it must return a numeric vector.")
  }
  if (is.list(w_edge) || is.data.frame(w_edge)) {
    stop("transform returned a list/data.frame; it must return an atomic numeric vector.")
  }
  if (!is.atomic(w_edge)) {
    stop("transform returned a non-atomic object; it must return an atomic numeric vector.")
  }
  
  w_edge <- as.numeric(w_edge)
  
  if (length(w_edge) != length(w_edge_raw)) {
    stop(
      "transform returned length ", length(w_edge),
      " but expected length ", length(w_edge_raw), "."
    )
  }
  
  w_edge[!is.finite(w_edge)] <- NA_real_
  
  list(
    edge_signal = w_edge,
    edge_signal_raw = w_edge_raw,
    mechanisms = selected_mechanisms,
    mechanism_idx = idx,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    transform = transform
  )
}

#' Build a plot title from selected LR mechanisms
makeLRPlotTitle <- function(selected_mechanisms,
                            plot_title = NULL,
                            suffix = NULL) {
  if (!is.null(plot_title)) {
    moi <- as.character(plot_title)[1]
  } else {
    if (length(selected_mechanisms) == 1L) {
      moi <- selected_mechanisms
    } else if (length(selected_mechanisms) <= 4L) {
      moi <- paste(selected_mechanisms, collapse = ", ")
    } else {
      moi <- paste0(length(selected_mechanisms), " LR mechanisms")
    }
  }
  
  if (!is.null(suffix) && nzchar(suffix)) {
    moi <- paste0(moi, suffix)
  }
  
  moi
}

#' Resolve plotting coordinates aligned to eff$nodes$ids
resolveEffLayoutCoords <- function(eff,
                                   layout = c("spatial", "umap"),
                                   coords = NULL) {
  stopifnot(inherits(eff, "EffNICHES"))
  layout <- match.arg(layout)
  
  xy_spatial <- eff$nodes$xy
  ids <- eff$nodes$ids
  
  if (is.null(xy_spatial) || !is.matrix(xy_spatial) || ncol(xy_spatial) != 2) {
    stop("eff$nodes$xy must be an N x 2 matrix.")
  }
  if (is.null(ids) || length(ids) != nrow(xy_spatial)) {
    stop("eff$nodes$ids must match eff$nodes$xy.")
  }
  
  if (is.null(colnames(xy_spatial)) || !all(c("x", "y") %in% colnames(xy_spatial))) {
    colnames(xy_spatial) <- c("x", "y")
  }
  
  if (layout == "spatial") {
    xy <- xy_spatial
    storage.mode(xy) <- "numeric"
    if (anyNA(xy)) stop("Spatial coordinates contain NA values.")
    return(list(
      xy = xy,
      layout_name = "spatial",
      node_ids = ids
    ))
  }
  
  # -------------------------
  # UMAP resolution
  # -------------------------
  out <- coords
  
  if (is.null(out)) {
    if (!is.null(eff$nodes$umap)) {
      out <- eff$nodes$umap
    } else if (!is.null(eff$nodes$reductions$umap)) {
      out <- eff$nodes$reductions$umap
    } else if (!is.null(eff$umap)) {
      out <- eff$umap
    } else if (!is.null(eff$reductions$umap)) {
      out <- eff$reductions$umap
    }
  }
  
  if (is.null(out)) {
    stop(
      "UMAP layout requested, but no UMAP coordinates were found. ",
      "Provide coords=, or store them in eff$nodes$umap, eff$nodes$reductions$umap, ",
      "eff$umap, or eff$reductions$umap."
    )
  }
  
  if (inherits(out, c("DimReduc", "dimreduc"))) {
    if (!is.null(out@cell.embeddings)) {
      out <- out@cell.embeddings
    }
  }
  
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  if (ncol(out) < 2) {
    stop("UMAP coordinates must have at least 2 columns.")
  }
  out <- out[, 1:2, drop = FALSE]
  
  rn <- rownames(out)
  if (!is.null(rn) && all(ids %in% rn)) {
    out <- out[ids, , drop = FALSE]
  } else {
    if (nrow(out) != length(ids)) {
      stop(
        "UMAP coordinates could not be aligned to eff$nodes$ids. ",
        "Supply rownames matching eff$nodes$ids, or provide exactly one row per node ",
        "in the same order as eff$nodes$ids."
      )
    }
  }
  
  xy <- as.matrix(out)
  storage.mode(xy) <- "numeric"
  colnames(xy) <- c("x", "y")
  
  if (anyNA(xy)) stop("UMAP coordinates contain NA values.")
  
  list(
    xy = xy,
    layout_name = "umap",
    node_ids = ids
  )
}

#' Prepare one LR plot panel worth of data without rendering
prepareLRPanelData <- function(eff,
                               lr,
                               layout = c("spatial", "umap"),
                               coords = NULL,
                               transform = identity,
                               # LR combination controls
                               lr_normalization = c("none", "max", "percentile", "zscore"),
                               lr_aggregation   = c("sum", "mean", "max"),
                               deduplicate_lr   = TRUE,
                               warn_many_lr     = TRUE,
                               # title controls
                               plot_title       = NULL,
                               title_suffix     = NULL,
                               # visualization filtering
                               edge_fraction = 1,
                               edge_top = NULL,
                               # safety
                               max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  layout <- match.arg(layout)
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  
  if (is.null(eff$niches$CellToCellSpatial)) {
    stop("eff$niches$CellToCellSpatial is missing.")
  }
  
  c2c <- eff$niches$CellToCellSpatial
  ij  <- c2c$ij
  
  if (is.null(ij) || !is.matrix(ij) || ncol(ij) != 2) {
    stop("CellToCellSpatial$ij must be a matrix with 2 columns (node indices).")
  }
  if (nrow(ij) > max_edges) {
    stop("Edge count (", nrow(ij), ") exceeds max_edges (", max_edges, ").")
  }
  
  coords_out <- resolveEffLayoutCoords(
    eff = eff,
    layout = layout,
    coords = coords
  )
  
  lr_out <- computeEdgeLRsignal(
    eff = eff,
    lr = lr,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    transform = transform,
    warn_many_lr = warn_many_lr,
    max_edges = max_edges
  )
  
  edge_df <- buildLREdgePlotData(
    ij = ij,
    xy = coords_out$xy,
    edge_signal = lr_out$edge_signal,
    edge_fraction = edge_fraction,
    edge_top = edge_top
  )
  
  node_df <- buildNodePlotData(coords_out$xy)
  
  title_use <- makeLRPlotTitle(
    selected_mechanisms = lr_out$selected_mechanisms,
    plot_title = plot_title,
    suffix = title_suffix
  )
  
  list(
    edge_df = edge_df,
    node_df = node_df,
    edge_signal = lr_out$edge_signal,
    edge_signal_raw = lr_out$edge_signal_raw,
    selected_idx = lr_out$selected_idx,
    selected_mechanisms = lr_out$selected_mechanisms,
    normalization = lr_out$normalization,
    aggregation = lr_out$aggregation,
    layout_name = coords_out$layout_name,
    node_ids = coords_out$node_ids,
    plot_title = title_use
  )
}

#' Convert node coordinates to a node plotting data.frame
buildNodePlotData <- function(xy) {
  if (is.null(xy) || !is.matrix(xy) || ncol(xy) != 2) {
    stop("xy must be an N x 2 matrix.")
  }
  
  if (is.null(colnames(xy)) || !all(c("x", "y") %in% colnames(xy))) {
    colnames(xy) <- c("x", "y")
  }
  
  df_nodes <- data.frame(
    x = suppressWarnings(as.numeric(xy[, "x"])),
    y = suppressWarnings(as.numeric(xy[, "y"])),
    stringsAsFactors = FALSE
  )
  
  if (anyNA(df_nodes$x) || anyNA(df_nodes$y)) {
    stop("NA detected in node coordinates.")
  }
  
  df_nodes
}


#' Build edge plotting data from indices, coordinates, and edge signal
buildLREdgePlotData <- function(ij,
                                xy,
                                edge_signal,
                                edge_fraction = 1,
                                edge_top = NULL) {
  if (is.null(ij) || !is.matrix(ij) || ncol(ij) != 2) {
    stop("ij must be a matrix with 2 columns.")
  }
  if (is.null(xy) || !is.matrix(xy) || ncol(xy) != 2) {
    stop("xy must be an N x 2 matrix.")
  }
  if (length(edge_signal) != nrow(ij)) {
    stop("edge_signal must have length equal to nrow(ij).")
  }
  
  n_nodes <- nrow(xy)
  ij_num <- suppressWarnings(matrix(as.integer(ij), ncol = 2))
  
  if (anyNA(ij_num)) stop("NA detected in ij.")
  if (any(ij_num < 1L) || any(ij_num > n_nodes)) {
    stop("ij contains node indices out of bounds.")
  }
  
  if (!is.numeric(edge_fraction) || length(edge_fraction) != 1L ||
      !is.finite(edge_fraction) || edge_fraction <= 0 || edge_fraction > 1) {
    stop("edge_fraction must be a single number in (0, 1].")
  }
  
  if (!is.null(edge_top)) {
    if (!is.numeric(edge_top) || length(edge_top) != 1L ||
        !is.finite(edge_top) || edge_top < 1) {
      stop("edge_top must be NULL or a single positive number.")
    }
    edge_top <- as.integer(edge_top)
  }
  
  from_idx <- ij_num[, 1]
  to_idx   <- ij_num[, 2]
  
  x1 <- suppressWarnings(as.numeric(xy[from_idx, "x"]))
  y1 <- suppressWarnings(as.numeric(xy[from_idx, "y"]))
  x2 <- suppressWarnings(as.numeric(xy[to_idx,   "x"]))
  y2 <- suppressWarnings(as.numeric(xy[to_idx,   "y"]))
  
  if (anyNA(x1) || anyNA(y1) || anyNA(x2) || anyNA(y2)) {
    stop("NA detected in constructed edge coordinates.")
  }
  
  keep <- !is.na(edge_signal) & is.finite(edge_signal)
  if (!any(keep)) stop("No finite values to plot after transform().")
  
  finite_idx <- which(keep)
  if (length(finite_idx) > 0L && (edge_fraction < 1 || !is.null(edge_top))) {
    ord_desc <- order(edge_signal[finite_idx], decreasing = TRUE)
    
    n_keep_fraction <- ceiling(length(finite_idx) * edge_fraction)
    n_keep_final <- n_keep_fraction
    
    if (!is.null(edge_top)) {
      n_keep_final <- min(n_keep_final, edge_top)
    }
    n_keep_final <- max(1L, min(length(finite_idx), n_keep_final))
    
    selected_idx <- finite_idx[ord_desc[seq_len(n_keep_final)]]
    keep <- rep(FALSE, length(edge_signal))
    keep[selected_idx] <- TRUE
  }
  
  if (!any(keep)) stop("No edges remained after filtering.")
  
  data.frame(
    x1 = x1[keep],
    y1 = y1[keep],
    x2 = x2[keep],
    y2 = y2[keep],
    for.plotting = edge_signal[keep],
    stringsAsFactors = FALSE
  )
}

#' Shared color limits across prepared LR panels
computeSharedLRColorLimits <- function(panel_data_list,
                                       color_limits = NULL) {
  if (is.null(color_limits)) {
    all_vals <- unlist(
      lapply(panel_data_list, function(x) x$edge_signal),
      use.names = FALSE
    )
    all_vals <- all_vals[is.finite(all_vals)]
    
    if (length(all_vals) == 0L) {
      stop("No finite edge_signal values were found across panels.")
    }
    
    scale_limits <- range(all_vals, na.rm = TRUE)
  } else {
    if (!is.numeric(color_limits) || length(color_limits) != 2 || any(!is.finite(color_limits))) {
      stop("color_limits must be NULL or a numeric vector of length 2.")
    }
    if (color_limits[1] >= color_limits[2]) {
      stop("color_limits must satisfy color_limits[1] < color_limits[2].")
    }
    scale_limits <- color_limits
  }
  
  if (diff(scale_limits) == 0) {
    eps <- max(abs(scale_limits[1]), 1) * 1e-9
    scale_limits <- scale_limits + c(-eps, eps)
  }
  
  scale_limits
}

#' Validate / derive color limits
validateColorLimits <- function(values, color_limits = NULL) {
  observed_range <- range(values, na.rm = TRUE)
  
  if (is.null(color_limits)) {
    scale_limits <- observed_range
  } else {
    if (!is.numeric(color_limits) || length(color_limits) != 2 || any(!is.finite(color_limits))) {
      stop("color_limits must be NULL or a numeric vector of length 2.")
    }
    if (color_limits[1] >= color_limits[2]) {
      stop("color_limits must satisfy color_limits[1] < color_limits[2].")
    }
    scale_limits <- color_limits
  }
  
  if (diff(scale_limits) == 0) {
    eps <- max(abs(scale_limits[1]), 1) * 1e-9
    scale_limits <- scale_limits + c(-eps, eps)
  }
  
  list(
    observed_range = observed_range,
    scale_limits = scale_limits
  )
}

#' Remove axes from a rendered LR panel
stripLRPanelAxes <- function(p) {
  p + ggplot2::theme(
    axis.title   = ggplot2::element_blank(),
    axis.text    = ggplot2::element_blank(),
    axis.ticks   = ggplot2::element_blank(),
    axis.line    = ggplot2::element_blank(),
    axis.line.x  = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )
}

#' Render a prepared LR panel
renderPreparedLRPanel <- function(panel_data,
                                  plot_title = NULL,
                                  show_axes = TRUE,
                                  reverse_y = FALSE,
                                  edge_width_bg = 0.10,
                                  edge_alpha_bg = 0.075,
                                  edge_color_bg = "white",
                                  edge_width_fg = 0.20,
                                  pt_size = 0,
                                  pt_color = "white",
                                  alpha_range = c(0.05, 1),
                                  alpha_limits = NULL,
                                  color_option = "C",
                                  color_limits = NULL,
                                  oob_squish = TRUE,
                                  x_limits = NULL,
                                  y_limits = NULL) {
  if (!is.list(panel_data) ||
      is.null(panel_data$edge_df) ||
      is.null(panel_data$node_df)) {
    stop("panel_data must be a list returned by prepareLRPanelData().")
  }
  
  title_use <- if (is.null(plot_title)) panel_data$plot_title else plot_title
  
  p <- renderLRNetworkPlot(
    edge_df = panel_data$edge_df,
    node_df = panel_data$node_df,
    plot_title = title_use,
    reverse_y = reverse_y,
    edge_width_bg = edge_width_bg,
    edge_alpha_bg = edge_alpha_bg,
    edge_color_bg = edge_color_bg,
    edge_width_fg = edge_width_fg,
    pt_size = pt_size,
    pt_color = pt_color,
    alpha_range = alpha_range,
    alpha_limits = alpha_limits,
    color_option = color_option,
    color_limits = color_limits,
    oob_squish = oob_squish,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  if (!isTRUE(show_axes)) {
    p <- stripLRPanelAxes(p)
  }
  
  p
}

#' Internal LR panel-grid engine
.plotLRPanelGrid <- function(
    eff_list,
    lr_list,
    layout = c("spatial", "umap"),
    coords_list = NULL,
    plot_titles = NULL,
    ncol = NULL,
    nrow = NULL,
    guides = "collect",
    share_color_limits = TRUE,
    share_alpha_limits = TRUE,
    show_axes = TRUE,
    transform = identity,
    # visual controls
    edge_width_bg = 0.10,
    edge_alpha_bg = 0.075,
    edge_color_bg = "white",
    edge_width_fg = 0.20,
    pt_size = 0.0,
    pt_color = "white",
    reverse_y = FALSE,
    # scaling controls
    alpha_range = c(0.05, 1),
    alpha_limits = NULL,
    color_option = "C",
    color_limits = NULL,
    oob_squish = TRUE,
    # coordinate controls
    x_limits = NULL,
    y_limits = NULL,
    # LR combination controls
    lr_normalization = c("none", "max", "percentile", "zscore"),
    lr_aggregation   = c("sum", "mean", "max"),
    deduplicate_lr   = TRUE,
    warn_many_lr     = TRUE,
    title_suffix     = NULL,
    # visualization filtering
    edge_fraction = 1,
    edge_top = NULL,
    # safety
    max_edges = 2e6) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  
  layout <- match.arg(layout)
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  
  if (!is.list(eff_list) || length(eff_list) == 0L ||
      !all(vapply(eff_list, inherits, logical(1), "EffNICHES"))) {
    stop("eff_list must be a non-empty list of EffNICHES objects.")
  }
  
  if (!is.list(lr_list) || length(lr_list) == 0L) {
    stop("lr_list must be a non-empty list of LR specifications.")
  }
  
  n_panels <- length(eff_list)
  
  if (length(lr_list) != n_panels) {
    stop("lr_list must have the same length as eff_list.")
  }
  
  if (is.null(coords_list)) {
    coords_list <- rep(list(NULL), n_panels)
  } else {
    if (!is.list(coords_list)) {
      coords_list <- rep(list(coords_list), n_panels)
    }
    if (length(coords_list) != n_panels) {
      stop("coords_list must be NULL, a single coords object, or a list with length equal to eff_list.")
    }
  }
  
  if (is.null(plot_titles)) {
    plot_titles <- rep(NA_character_, n_panels)
  } else {
    if (length(plot_titles) != n_panels) {
      stop("plot_titles must have the same length as eff_list.")
    }
    plot_titles <- as.character(plot_titles)
  }
  
  panel_data_list <- vector("list", n_panels)
  
  for (i in seq_len(n_panels)) {
    panel_data_list[[i]] <- prepareLRPanelData(
      eff = eff_list[[i]],
      lr = lr_list[[i]],
      layout = layout,
      coords = coords_list[[i]],
      transform = transform,
      lr_normalization = lr_normalization,
      lr_aggregation = lr_aggregation,
      deduplicate_lr = deduplicate_lr,
      warn_many_lr = warn_many_lr,
      plot_title = NULL,
      title_suffix = title_suffix,
      edge_fraction = edge_fraction,
      edge_top = edge_top,
      max_edges = max_edges
    )
  }
  
  if (isTRUE(share_color_limits)) {
    color_limits_use <- computeSharedLRColorLimits(
      panel_data_list = panel_data_list,
      color_limits = color_limits
    )
  } else {
    color_limits_use <- color_limits
  }
  
  if (isTRUE(share_alpha_limits)) {
    alpha_limits_use <- computeSharedLRAlphaLimits(
      panel_data_list = panel_data_list,
      alpha_limits = alpha_limits
    )
  } else {
    alpha_limits_use <- alpha_limits
  }
  
  plot_list <- vector("list", n_panels)
  
  for (i in seq_len(n_panels)) {
    title_i <- if (is.na(plot_titles[i])) NULL else plot_titles[i]
    
    plot_list[[i]] <- renderPreparedLRPanel(
      panel_data = panel_data_list[[i]],
      plot_title = title_i,
      show_axes = show_axes,
      reverse_y = reverse_y,
      edge_width_bg = edge_width_bg,
      edge_alpha_bg = edge_alpha_bg,
      edge_color_bg = edge_color_bg,
      edge_width_fg = edge_width_fg,
      pt_size = pt_size,
      pt_color = pt_color,
      alpha_range = alpha_range,
      alpha_limits = alpha_limits_use,
      color_option = color_option,
      color_limits = color_limits_use,
      oob_squish = oob_squish,
      x_limits = x_limits,
      y_limits = y_limits
    )
  }
  
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n_panels))
  }
  
  patchwork::wrap_plots(
    plot_list,
    ncol = ncol,
    nrow = nrow,
    guides = guides
  ) &
    ggplot2::theme(
      legend.position = "right",
      plot.background  = ggplot2::element_rect(fill = "black", colour = NA),
      panel.background = ggplot2::element_rect(fill = "black", colour = NA)
    )
}

#' Compute alpha values for edge plotting
computeEdgeAlpha <- function(values,
                             alpha_range = c(0.05, 1),
                             from = NULL) {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }
  
  a_rng <- alpha_range
  if (length(a_rng) != 2 || any(!is.finite(a_rng)) ||
      a_rng[1] < 0 || a_rng[2] <= 0 || a_rng[1] >= a_rng[2]) {
    stop("alpha_range must be something like c(0.05, 1).")
  }
  
  observed_range <- range(values, na.rm = TRUE)
  
  if (is.null(from)) {
    scale_range <- observed_range
  } else {
    if (!is.numeric(from) || length(from) != 2 || any(!is.finite(from))) {
      stop("from must be NULL or a numeric vector of length 2.")
    }
    if (from[1] >= from[2]) {
      stop("from must satisfy from[1] < from[2].")
    }
    scale_range <- from
  }
  
  if (diff(scale_range) == 0) {
    alpha_val <- rep(a_rng[2], length(values))
  } else {
    alpha_val <- scales::rescale(
      values,
      to = a_rng,
      from = scale_range
    )
    alpha_val[!is.finite(alpha_val)] <- a_rng[1]
  }
  
  alpha_val
}

#' Shared alpha limits across prepared LR panels
computeSharedLRAlphaLimits <- function(panel_data_list,
                                       alpha_limits = NULL) {
  if (is.null(alpha_limits)) {
    all_vals <- unlist(
      lapply(panel_data_list, function(x) x$edge_signal),
      use.names = FALSE
    )
    all_vals <- all_vals[is.finite(all_vals)]
    
    if (length(all_vals) == 0L) {
      stop("No finite edge_signal values were found across panels.")
    }
    
    scale_limits <- range(all_vals, na.rm = TRUE)
  } else {
    if (!is.numeric(alpha_limits) || length(alpha_limits) != 2L || any(!is.finite(alpha_limits))) {
      stop("alpha_limits must be NULL or a numeric vector of length 2.")
    }
    if (alpha_limits[1] > alpha_limits[2]) {
      stop("alpha_limits must satisfy alpha_limits[1] <= alpha_limits[2].")
    }
    scale_limits <- alpha_limits
  }
  
  scale_limits
}

#' Dark theme used by LR plot renderers
makeDarkTheme <- function() {
  if (exists("DarkTheme", mode = "function")) {
    DarkTheme()
  } else {
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = "black", color = NA),
      panel.background  = ggplot2::element_rect(fill = "black", color = NA),
      panel.grid        = ggplot2::element_blank(),
      axis.text         = ggplot2::element_text(color = "white"),
      axis.title        = ggplot2::element_text(color = "white"),
      axis.ticks        = ggplot2::element_line(color = "white"),
      axis.line         = ggplot2::element_line(color = "white"),
      plot.title        = ggplot2::element_text(color = "white"),
      legend.background = ggplot2::element_rect(fill = "black", color = NA),
      legend.key        = ggplot2::element_rect(fill = "black", color = NA),
      legend.text       = ggplot2::element_text(color = "white"),
      legend.title      = ggplot2::element_text(color = "white")
    )
  }
}


#' Render LR network plot from prepared data
renderLRNetworkPlot <- function(edge_df,
                                node_df = NULL,
                                plot_title = NULL,
                                reverse_y = FALSE,
                                edge_width_bg = 0.10,
                                edge_alpha_bg = 0.075,
                                edge_color_bg = "white",
                                edge_width_fg = 0.20,
                                pt_size = 0,
                                pt_color = "white",
                                alpha_range = c(0.05, 1),
                                alpha_limits = NULL,
                                color_option = "C",
                                color_limits = NULL,
                                oob_squish = TRUE,
                                x_limits = NULL,
                                y_limits = NULL) {
  if (!is.data.frame(edge_df) || !all(c("x1", "y1", "x2", "y2", "for.plotting") %in% names(edge_df))) {
    stop("edge_df must contain x1, y1, x2, y2, and for.plotting.")
  }
  
  if (!is.null(x_limits)) {
    if (!is.numeric(x_limits) || length(x_limits) != 2L || any(!is.finite(x_limits)) || x_limits[1] >= x_limits[2]) {
      stop("x_limits must be NULL or a numeric vector of length 2 with x_limits[1] < x_limits[2].")
    }
  }
  
  if (!is.null(y_limits)) {
    if (!is.numeric(y_limits) || length(y_limits) != 2L || any(!is.finite(y_limits)) || y_limits[1] >= y_limits[2]) {
      stop("y_limits must be NULL or a numeric vector of length 2 with y_limits[1] < y_limits[2].")
    }
  }
  
  lims <- validateColorLimits(edge_df$for.plotting, color_limits = color_limits)
  edge_df$alpha_val <- computeEdgeAlpha(
    edge_df$for.plotting,
    alpha_range = alpha_range,
    from = alpha_limits
  )
  edge_df <- edge_df[order(edge_df$for.plotting), , drop = FALSE]
  
  color_oob_fun <- if (isTRUE(oob_squish)) scales::squish else scales::censor
  
  p <- ggplot2::ggplot()
  
  if (!is.null(node_df) && isTRUE(pt_size > 0)) {
    p <- p + ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      size = pt_size,
      color = pt_color,
      inherit.aes = FALSE
    )
  }
  
  p <- p +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = .data$x1, y = .data$y1,
                   xend = .data$x2, yend = .data$y2),
      linewidth = edge_width_bg,
      alpha = edge_alpha_bg,
      color = edge_color_bg
    ) +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = .data$x1, y = .data$y1,
                   xend = .data$x2, yend = .data$y2,
                   color = .data$for.plotting,
                   alpha = .data$alpha_val),
      linewidth = edge_width_fg
    ) +
    ggplot2::scale_color_viridis_c(
      option = color_option,
      limits = lims$scale_limits,
      oob = color_oob_fun
    ) +
    ggplot2::scale_alpha_continuous(range = alpha_range, guide = "none") +
    ggplot2::coord_fixed(
      xlim = x_limits,
      ylim = y_limits,
      expand = FALSE
    ) +
    ggplot2::theme_classic() +
    makeDarkTheme() +
    ggplot2::ggtitle(plot_title)
  
  if (isTRUE(reverse_y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  
  p
}

#' Single-panel LR plotting engine
.plotLRLayout <- function(eff,
                          lr,
                          layout = c("spatial", "umap"),
                          coords = NULL,
                          transform = identity,
                          # visual controls
                          edge_width_bg = 0.10,
                          edge_alpha_bg = 0.075,
                          edge_color_bg = "white",
                          edge_width_fg = 0.20,
                          pt_size = 0,
                          pt_color = "white",
                          reverse_y = FALSE,
                          # scaling controls
                          alpha_range = c(0.05, 1),
                          color_option = "C",
                          color_limits = NULL,
                          oob_squish = TRUE,
                          # LR combination controls
                          lr_normalization = c("none", "max", "percentile", "zscore"),
                          lr_aggregation   = c("sum", "mean", "max"),
                          deduplicate_lr   = TRUE,
                          warn_many_lr     = TRUE,
                          plot_title       = NULL,
                          title_suffix     = NULL,
                          # visualization filtering
                          edge_fraction = 1,
                          edge_top = NULL,
                          # safety
                          max_edges = 2e6) {
  panel_data <- prepareLRPanelData(
    eff = eff,
    lr = lr,
    layout = layout,
    coords = coords,
    transform = transform,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    warn_many_lr = warn_many_lr,
    plot_title = plot_title,
    title_suffix = title_suffix,
    edge_fraction = edge_fraction,
    edge_top = edge_top,
    max_edges = max_edges
  )
  
  renderPreparedLRPanel(
    panel_data = panel_data,
    plot_title = NULL,
    show_axes = TRUE,
    reverse_y = reverse_y,
    edge_width_bg = edge_width_bg,
    edge_alpha_bg = edge_alpha_bg,
    edge_color_bg = edge_color_bg,
    edge_width_fg = edge_width_fg,
    pt_size = pt_size,
    pt_color = pt_color,
    alpha_range = alpha_range,
    color_option = color_option,
    color_limits = color_limits,
    oob_squish = oob_squish
  )
}

#' Pull one node-level grouping variable aligned to eff$nodes$ids
getEffNodeGrouping <- function(eff, group_by) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  if (!is.character(group_by) || length(group_by) != 1L || !nzchar(group_by)) {
    stop("group_by must be a single non-empty character string.")
  }
  
  ids <- eff$nodes$ids
  if (is.null(ids) || length(ids) < 1L) {
    stop("eff$nodes$ids is missing or empty.")
  }
  
  meta <- eff$nodes$meta
  if (is.null(meta) || !is.data.frame(meta)) {
    stop("eff$nodes$meta must be a data.frame containing group_by.")
  }
  if (!(group_by %in% colnames(meta))) {
    stop("group_by not found in eff$nodes$meta: ", group_by)
  }
  
  if (!is.null(rownames(meta))) {
    o <- match(ids, rownames(meta))
    if (anyNA(o)) {
      if (nrow(meta) != length(ids)) {
        stop("eff$nodes$meta rownames do not cover eff$nodes$ids and nrow(meta) != N.")
      }
      vals <- meta[[group_by]]
    } else {
      vals <- meta[o, group_by, drop = TRUE]
    }
  } else {
    if (nrow(meta) != length(ids)) {
      stop("eff$nodes$meta has no rownames and nrow(meta) != N; cannot align.")
    }
    vals <- meta[[group_by]]
  }
  
  vals <- as.character(vals)
  
  if (anyNA(vals) || any(!nzchar(vals))) {
    stop("group_by contains NA or empty values after alignment: ", group_by)
  }
  
  vals
}

#' Build typed edge drawing data.
#' If edge_gradient = FALSE, each edge is split into two hard color stepoffs at midpoint.
#' If edge_gradient = TRUE, each edge is approximated by multiple short colored segments.
buildTypedLREdgeDrawData <- function(edge_df,
                                     palette_map,
                                     edge_gradient = FALSE,
                                     edge_n_segments = 9L) {
  if (!is.data.frame(edge_df) ||
      !all(c("x1", "y1", "x2", "y2",
             "for.plotting", "alpha_val",
             "sender_type", "receiver_type",
             "edge_order") %in% names(edge_df))) {
    stop("edge_df is missing required columns.")
  }
  
  used_types <- sort(unique(c(edge_df$sender_type, edge_df$receiver_type)))
  missing_cols <- setdiff(used_types, names(palette_map))
  if (length(missing_cols) > 0L) {
    stop(
      "palette is missing colors for: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  sender_col <- unname(palette_map[edge_df$sender_type])
  receiver_col <- unname(palette_map[edge_df$receiver_type])
  
  # --------------------------------
  # Fast mode: hard color step at midpoint
  # --------------------------------
  if (!isTRUE(edge_gradient)) {
    xm <- (edge_df$x1 + edge_df$x2) / 2
    ym <- (edge_df$y1 + edge_df$y2) / 2
    
    n <- nrow(edge_df)
    
    out1 <- data.frame(
      x = edge_df$x1,
      y = edge_df$y1,
      xend = xm,
      yend = ym,
      alpha_val = edge_df$alpha_val,
      segment_color = sender_col,
      edge_order = edge_df$edge_order * 2L - 1L,
      stringsAsFactors = FALSE
    )
    
    out2 <- data.frame(
      x = xm,
      y = ym,
      xend = edge_df$x2,
      yend = edge_df$y2,
      alpha_val = edge_df$alpha_val,
      segment_color = receiver_col,
      edge_order = edge_df$edge_order * 2L,
      stringsAsFactors = FALSE
    )
    
    out <- rbind(out1, out2)
    return(out[order(out$edge_order), , drop = FALSE])
  }
  
  # --------------------------------
  # Gradient mode: segmented approximation
  # --------------------------------
  edge_n_segments <- as.integer(edge_n_segments)
  if (!is.finite(edge_n_segments) || edge_n_segments < 3L) {
    stop("edge_n_segments must be an integer >= 3.")
  }
  
  pair_df <- unique(edge_df[, c("sender_type", "receiver_type"), drop = FALSE])
  pair_key <- paste(pair_df$sender_type, pair_df$receiver_type, sep = " -> ")
  
  color_profiles <- vector("list", length(pair_key))
  names(color_profiles) <- pair_key
  
  for (i in seq_len(nrow(pair_df))) {
    s_type <- pair_df$sender_type[i]
    r_type <- pair_df$receiver_type[i]
    
    s_col <- unname(palette_map[[s_type]])
    r_col <- unname(palette_map[[r_type]])
    
    ramp_fun <- grDevices::colorRampPalette(c(s_col, r_col), space = "Lab")
    color_profiles[[i]] <- ramp_fun(edge_n_segments)
  }
  
  n_edges <- nrow(edge_df)
  out_n <- n_edges * edge_n_segments
  
  out <- data.frame(
    x = numeric(out_n),
    y = numeric(out_n),
    xend = numeric(out_n),
    yend = numeric(out_n),
    alpha_val = numeric(out_n),
    segment_color = character(out_n),
    edge_order = integer(out_n),
    stringsAsFactors = FALSE
  )
  
  idx_start <- 1L
  
  for (i in seq_len(n_edges)) {
    idx_end <- idx_start + edge_n_segments - 1L
    
    x_seq <- seq(edge_df$x1[i], edge_df$x2[i], length.out = edge_n_segments + 1L)
    y_seq <- seq(edge_df$y1[i], edge_df$y2[i], length.out = edge_n_segments + 1L)
    
    pair_name <- paste(edge_df$sender_type[i], edge_df$receiver_type[i], sep = " -> ")
    seg_cols <- color_profiles[[pair_name]]
    
    out$x[idx_start:idx_end] <- x_seq[-(edge_n_segments + 1L)]
    out$y[idx_start:idx_end] <- y_seq[-(edge_n_segments + 1L)]
    out$xend[idx_start:idx_end] <- x_seq[-1L]
    out$yend[idx_start:idx_end] <- y_seq[-1L]
    out$alpha_val[idx_start:idx_end] <- edge_df$alpha_val[i]
    out$segment_color[idx_start:idx_end] <- seg_cols
    out$edge_order[idx_start:idx_end] <- edge_df$edge_order * edge_n_segments + seq_len(edge_n_segments) - edge_n_segments
    
    idx_start <- idx_end + 1L
  }
  
  out[order(out$edge_order), , drop = FALSE]
}

#' Internal helper: collect drawn cell types for one typed LR panel
.collectTypedLRUsedTypes <- function(eff,
                                     lr,
                                     group_by,
                                     transform = identity,
                                     cell_types = NULL,
                                     cell_types_mode = c("either", "both"),
                                     sender_types = NULL,
                                     receiver_types = NULL,
                                     lr_normalization = c("none", "max", "percentile", "zscore"),
                                     lr_aggregation   = c("sum", "mean", "max"),
                                     deduplicate_lr   = TRUE,
                                     warn_many_lr     = TRUE,
                                     edge_fraction = 1,
                                     edge_top = NULL,
                                     max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  cell_types_mode  <- match.arg(cell_types_mode)
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  
  if (is.null(eff$niches$CellToCellSpatial)) {
    stop("eff$niches$CellToCellSpatial is missing.")
  }
  
  c2c <- eff$niches$CellToCellSpatial
  ij  <- c2c$ij
  
  if (is.null(ij) || !is.matrix(ij) || ncol(ij) != 2) {
    stop("CellToCellSpatial$ij must be a matrix with 2 columns (node indices).")
  }
  if (nrow(ij) > max_edges) {
    stop("Edge count (", nrow(ij), ") exceeds max_edges (", max_edges, ").")
  }
  
  coords_out <- resolveEffLayoutCoords(
    eff = eff,
    layout = "spatial",
    coords = NULL
  )
  
  lr_out <- computeEdgeLRsignal(
    eff = eff,
    lr = lr,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    transform = transform,
    warn_many_lr = warn_many_lr,
    max_edges = max_edges
  )
  
  node_groups <- getEffNodeGrouping(eff, group_by = group_by)
  xy <- coords_out$xy
  
  from_idx <- as.integer(ij[, 1])
  to_idx   <- as.integer(ij[, 2])
  
  edge_df <- data.frame(
    for.plotting = as.numeric(lr_out$edge_signal),
    sender_type = node_groups[from_idx],
    receiver_type = node_groups[to_idx],
    stringsAsFactors = FALSE
  )
  
  keep <- is.finite(edge_df$for.plotting) &
    !is.na(edge_df$for.plotting) &
    (edge_df$for.plotting > 0)
  
  if (!any(keep)) return(character(0))
  edge_df <- edge_df[keep, , drop = FALSE]
  
  if (!is.numeric(edge_fraction) || length(edge_fraction) != 1L ||
      !is.finite(edge_fraction) || edge_fraction <= 0 || edge_fraction > 1) {
    stop("edge_fraction must be a single number in (0, 1].")
  }
  
  if (!is.null(edge_top)) {
    if (!is.numeric(edge_top) || length(edge_top) != 1L ||
        !is.finite(edge_top) || edge_top < 1) {
      stop("edge_top must be NULL or a single positive number.")
    }
    edge_top <- as.integer(edge_top)
  }
  
  if (nrow(edge_df) > 0L && (edge_fraction < 1 || !is.null(edge_top))) {
    ord_desc <- order(edge_df$for.plotting, decreasing = TRUE)
    n_keep_fraction <- ceiling(nrow(edge_df) * edge_fraction)
    n_keep_final <- n_keep_fraction
    
    if (!is.null(edge_top)) {
      n_keep_final <- min(n_keep_final, edge_top)
    }
    n_keep_final <- max(1L, min(nrow(edge_df), n_keep_final))
    
    edge_df <- edge_df[ord_desc[seq_len(n_keep_final)], , drop = FALSE]
  }
  
  normalize_type_filter <- function(x, arg_name) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- unique(trimws(x))
    x <- x[nzchar(x)]
    if (length(x) < 1L) {
      stop(arg_name, " must contain at least one non-empty cell type or be NULL.")
    }
    x
  }
  
  cell_types     <- normalize_type_filter(cell_types, "cell_types")
  sender_types   <- normalize_type_filter(sender_types, "sender_types")
  receiver_types <- normalize_type_filter(receiver_types, "receiver_types")
  
  role_keep <- rep(TRUE, nrow(edge_df))
  
  if (!is.null(cell_types)) {
    if (cell_types_mode == "either") {
      role_keep <- role_keep & (
        edge_df$sender_type %in% cell_types |
          edge_df$receiver_type %in% cell_types
      )
    } else {
      role_keep <- role_keep & (
        edge_df$sender_type %in% cell_types &
          edge_df$receiver_type %in% cell_types
      )
    }
  }
  
  if (!is.null(sender_types)) {
    role_keep <- role_keep & (edge_df$sender_type %in% sender_types)
  }
  
  if (!is.null(receiver_types)) {
    role_keep <- role_keep & (edge_df$receiver_type %in% receiver_types)
  }
  
  edge_df <- edge_df[role_keep, , drop = FALSE]
  if (nrow(edge_df) < 1L) return(character(0))
  
  sort(unique(c(edge_df$sender_type, edge_df$receiver_type)))
}

#' Internal helper: remove axes from typed plot panels
stripTypedLRPanelAxes <- function(p) {
  p + ggplot2::theme(
    axis.title   = ggplot2::element_blank(),
    axis.text    = ggplot2::element_blank(),
    axis.ticks   = ggplot2::element_blank(),
    axis.line    = ggplot2::element_blank(),
    axis.line.x  = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )
}

#' Compute per-object spatial extents
getSpatialExtent <- function(eff) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  coords_out <- resolveEffLayoutCoords(
    eff = eff,
    layout = "spatial",
    coords = NULL
  )
  
  xy <- coords_out$xy
  
  xr <- range(xy[, "x"], na.rm = TRUE)
  yr <- range(xy[, "y"], na.rm = TRUE)
  
  list(
    x_range = xr,
    y_range = yr,
    x_mid = mean(xr),
    y_mid = mean(yr),
    width = diff(xr),
    height = diff(yr)
  )
}

#' Compute shared panel dimensions for across-object spatial plotting
computeSharedSpatialPanelSize <- function(eff_list) {
  if (!is.list(eff_list) || length(eff_list) == 0L ||
      !all(vapply(eff_list, inherits, logical(1), "EffNICHES"))) {
    stop("eff_list must be a non-empty list of EffNICHES objects.")
  }
  
  extents <- lapply(eff_list, getSpatialExtent)
  
  shared_width  <- max(vapply(extents, `[[`, numeric(1), "width"),  na.rm = TRUE)
  shared_height <- max(vapply(extents, `[[`, numeric(1), "height"), na.rm = TRUE)
  
  if (!is.finite(shared_width) || shared_width < 0) {
    stop("Could not compute shared spatial panel width.")
  }
  if (!is.finite(shared_height) || shared_height < 0) {
    stop("Could not compute shared spatial panel height.")
  }
  
  list(
    width = shared_width,
    height = shared_height
  )
}

#' Center one object's coordinates inside a shared plotting window
computeCenteredSpatialLimits <- function(eff,
                                         shared_width,
                                         shared_height) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  if (!is.numeric(shared_width) || length(shared_width) != 1L ||
      !is.finite(shared_width) || shared_width < 0) {
    stop("shared_width must be a single non-negative finite number.")
  }
  
  if (!is.numeric(shared_height) || length(shared_height) != 1L ||
      !is.finite(shared_height) || shared_height < 0) {
    stop("shared_height must be a single non-negative finite number.")
  }
  
  ext <- getSpatialExtent(eff)
  
  x_limits <- ext$x_mid + c(-0.5, 0.5) * shared_width
  y_limits <- ext$y_mid + c(-0.5, 0.5) * shared_height
  
  if (diff(x_limits) == 0) {
    eps <- max(abs(ext$x_mid), 1) * 1e-9
    x_limits <- x_limits + c(-eps, eps)
  }
  
  if (diff(y_limits) == 0) {
    eps <- max(abs(ext$y_mid), 1) * 1e-9
    y_limits <- y_limits + c(-eps, eps)
  }
  
  list(
    x_limits = x_limits,
    y_limits = y_limits
  )
}

#### Plot Node ####

plotNode <- function(eff,
                     group_by = NULL,
                     pt_size = 0.3,
                     palette = NULL,
                     reverse_y = FALSE,
                     legend_square_size = 4,
                     use_polygons = FALSE,
                     polygon_outline = NA,
                     polygon_alpha = 1,
                     log_transform = FALSE,
                     color_limits = NULL,
                     oob_squish = TRUE) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  # -------------------------
  # Pull xy / ids
  # -------------------------
  if (is.null(eff$nodes$xy) || ncol(eff$nodes$xy) != 2) {
    stop("eff$nodes$xy must be an N x 2 matrix with columns x,y.")
  }
  xy <- eff$nodes$xy
  if (is.null(colnames(xy)) || !all(c("x", "y") %in% colnames(xy))) {
    colnames(xy) <- c("x", "y")
  }
  
  ids <- eff$nodes$ids
  if (is.null(ids) || length(ids) != nrow(xy)) {
    stop("eff$nodes$ids must exist and match nrow(eff$nodes$xy).")
  }
  
  # -------------------------
  # Base centroid df (includes meta)
  # -------------------------
  df <- data.frame(
    id = ids,
    x = suppressWarnings(as.numeric(xy[, "x"])),
    y = suppressWarnings(as.numeric(xy[, "y"])),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(eff$nodes$meta)) {
    meta <- eff$nodes$meta
    
    # Align meta to ids
    if (!is.null(rownames(meta))) {
      o <- match(ids, rownames(meta))
      if (anyNA(o)) {
        if (nrow(meta) != length(ids)) {
          stop("eff$nodes$meta rownames do not cover eff$nodes$ids and nrow(meta) != N.")
        }
        df <- cbind(df, meta)
      } else {
        df <- cbind(df, meta[o, , drop = FALSE])
      }
    } else {
      if (nrow(meta) != length(ids)) {
        stop("eff$nodes$meta has no rownames and nrow(meta) != N; cannot align.")
      }
      df <- cbind(df, meta)
    }
  }
  
  # -------------------------
  # Helper: find an expression matrix inside eff
  #   Accept matrix/dgCMatrix with rownames=features, colnames=ids/cells
  # -------------------------
  .is_matlike <- function(x) {
    is.matrix(x) || inherits(x, "dgCMatrix") || inherits(x, "Matrix")
  }
  
  .candidate_mats <- function(eff) {
    out <- list()
    
    paths <- list(
      expr            = eff$expr,
      expr_mat        = eff$expr$mat,
      expression      = eff$expression,
      expression_mat  = eff$expression$mat,
      counts          = eff$counts,
      data            = eff$data,
      X               = eff$X,
      mat             = eff$mat
    )
    
    for (nm in names(paths)) {
      x <- paths[[nm]]
      if (.is_matlike(x)) out[[nm]] <- x
    }
    out
  }
  
  .pick_expr_mat <- function(eff, ids) {
    mats <- .candidate_mats(eff)
    if (length(mats) == 0) return(NULL)
    
    best <- NULL
    best_score <- -Inf
    for (nm in names(mats)) {
      m <- mats[[nm]]
      cn <- colnames(m)
      rn <- rownames(m)
      if (is.null(cn) || is.null(rn)) next
      
      score <- sum(ids %in% cn)
      if (score > best_score) {
        best <- m
        best_score <- score
      }
    }
    best
  }
  
  .fetch_values <- function(eff, df, ids, group_by) {
    if (is.null(group_by)) return(NULL)
    
    # Case 1: numeric vector provided directly
    if (is.numeric(group_by)) {
      if (length(group_by) != length(ids)) {
        stop("If group_by is a numeric vector, it must have length N = length(eff$nodes$ids).")
      }
      v <- as.numeric(group_by)
      names(v) <- NULL
      return(v)
    }
    
    # Case 2: string: meta column OR gene/feature from expression matrix
    if (!is.character(group_by) || length(group_by) != 1L) {
      stop("group_by must be NULL, a single column/feature name, or a numeric vector of length N.")
    }
    
    # 2a: meta column
    if (group_by %in% colnames(df)) {
      return(df[[group_by]])
    }
    
    # 2b: gene/feature from expression matrix
    m <- .pick_expr_mat(eff, ids)
    if (is.null(m)) {
      stop(
        "group_by not found in eff$nodes$meta (as a column), and no expression matrix was detected in eff.\n",
        "Tried: eff$expr, eff$expr$mat, eff$expression, eff$expression$mat, eff$counts, eff$data, eff$X, eff$mat."
      )
    }
    if (is.null(rownames(m))) stop("Detected expression matrix has no rownames (feature names).")
    if (is.null(colnames(m))) stop("Detected expression matrix has no colnames (cell/node ids).")
    
    if (!(group_by %in% rownames(m))) {
      stop("group_by not found as meta column or feature in detected expression matrix: ", group_by)
    }
    
    # align matrix columns to ids
    ci <- match(ids, colnames(m))
    if (anyNA(ci)) {
      stop("Expression matrix colnames do not cover all eff$nodes$ids; cannot align expression to nodes.")
    }
    
    v <- m[group_by, ci]
    as.numeric(v)
  }
  
  # -------------------------
  # Resolve values + legend title
  # -------------------------
  values <- .fetch_values(eff, df, ids, group_by)
  
  legend_title <- NULL
  if (!is.null(group_by)) {
    if (is.character(group_by) && length(group_by) == 1L) {
      legend_title <- group_by
    } else {
      legend_title <- "value"
    }
  }
  
  # -------------------------
  # Apply transform (numeric only)
  # -------------------------
  if (!is.null(values) && is.numeric(values)) {
    if (isTRUE(log_transform)) {
      values <- log1p(values)
      legend_title <- paste0("log1p(", legend_title, ")")
    } else if (is.function(log_transform)) {
      values <- log_transform(values)
      legend_title <- paste0("transformed(", legend_title, ")")
    } else if (!identical(log_transform, FALSE)) {
      stop("log_transform must be FALSE, TRUE (log1p), or a function.")
    }
  }
  
  if (!is.null(values)) {
    df[[".color_value"]] <- values
  }
  
  # -------------------------
  # Determine quantitative vs discrete
  # -------------------------
  is_quantitative <- !is.null(values) && is.numeric(values)
  
  # -------------------------
  # Validate color limits (quantitative only)
  # -------------------------
  scale_limits <- NULL
  if (is_quantitative) {
    observed_range <- range(df[[".color_value"]], na.rm = TRUE)
    
    if (is.null(color_limits)) {
      scale_limits <- observed_range
    } else {
      if (!is.numeric(color_limits) || length(color_limits) != 2 || any(!is.finite(color_limits))) {
        stop("color_limits must be NULL or a numeric vector of length 2.")
      }
      if (color_limits[1] >= color_limits[2]) {
        stop("color_limits must satisfy color_limits[1] < color_limits[2].")
      }
      scale_limits <- color_limits
    }
    
    if (!all(is.finite(scale_limits))) {
      stop("Could not determine finite color limits for quantitative plotting.")
    }
    if (diff(scale_limits) == 0) {
      eps <- max(abs(scale_limits[1]), 1) * 1e-9
      scale_limits <- scale_limits + c(-eps, eps)
    }
  }
  
  # -------------------------
  # Build polygon vertex df if requested; replicate meta by node id
  # -------------------------
  poly_df <- NULL
  if (isTRUE(use_polygons)) {
    poly <- eff$nodes$polygons
    if (is.null(poly) || is.null(poly$xy) || is.null(poly$start)) {
      stop("use_polygons=TRUE but eff$nodes$polygons is missing (needs $xy and $start).")
    }
    
    xyv <- poly$xy
    st  <- poly$start
    
    if (!is.matrix(xyv) || ncol(xyv) != 2) stop("eff$nodes$polygons$xy must be a V x 2 matrix.")
    if (!is.integer(st)) st <- as.integer(st)
    if (length(st) != length(ids) + 1L) {
      stop("eff$nodes$polygons$start must have length N+1, where N=length(eff$nodes$ids).")
    }
    
    counts <- diff(st)
    keep_cells <- which(counts >= 3L)
    if (length(keep_cells) == 0) {
      stop("use_polygons=TRUE but no cells have >= 3 vertices in eff$nodes$polygons.")
    }
    
    idx_list <- lapply(keep_cells, function(i) {
      a <- st[i]; b <- st[i + 1L] - 1L
      if (b < a) integer(0) else a:b
    })
    vidx <- unlist(idx_list, use.names = FALSE)
    rep_ids <- rep(ids[keep_cells], times = counts[keep_cells])
    
    poly_df <- data.frame(
      id = rep_ids,
      x  = suppressWarnings(as.numeric(xyv[vidx, 1])),
      y  = suppressWarnings(as.numeric(xyv[vidx, 2])),
      stringsAsFactors = FALSE
    )
    
    # replicate per-node columns (including .color_value if present) onto vertices
    meta_cols <- setdiff(colnames(df), c("id", "x", "y"))
    if (length(meta_cols) > 0) {
      mi <- match(poly_df$id, df$id)
      if (anyNA(mi)) stop("Internal error: polygon ids not found in centroid ids.")
      poly_df[meta_cols] <- df[mi, meta_cols, drop = FALSE]
    }
  }
  
  # -------------------------
  # Plot base
  # -------------------------
  if (isTRUE(use_polygons)) {
    p <- ggplot2::ggplot(poly_df, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id))
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y))
  }
  
  color_oob_fun <- if (isTRUE(oob_squish)) scales::squish else scales::censor
  
  # -------------------------
  # Layers + scales
  # -------------------------
  if (!is.null(values)) {
    v <- df[[".color_value"]]
    
    if (is_quantitative) {
      # continuous
      if (isTRUE(use_polygons)) {
        p <- p +
          ggplot2::geom_polygon(
            ggplot2::aes(fill = .data[[".color_value"]]),
            color = polygon_outline,
            alpha = polygon_alpha,
            linewidth = 0.05
          ) +
          ggplot2::guides(fill = ggplot2::guide_colorbar(title = legend_title))
        
        if (!is.null(palette)) {
          p <- p + ggplot2::scale_fill_gradientn(
            colors = palette,
            name = legend_title,
            limits = scale_limits,
            oob = color_oob_fun
          )
        } else {
          p <- p + ggplot2::scale_fill_viridis_c(
            name = legend_title,
            limits = scale_limits,
            oob = color_oob_fun
          )
        }
      } else {
        p <- p +
          ggplot2::geom_point(
            ggplot2::aes(color = .data[[".color_value"]]),
            size = pt_size
          ) +
          ggplot2::guides(color = ggplot2::guide_colorbar(title = legend_title))
        
        if (!is.null(palette)) {
          p <- p + ggplot2::scale_color_gradientn(
            colors = palette,
            name = legend_title,
            limits = scale_limits,
            oob = color_oob_fun
          )
        } else {
          p <- p + ggplot2::scale_color_viridis_c(
            name = legend_title,
            limits = scale_limits,
            oob = color_oob_fun
          )
        }
      }
      
    } else {
      # discrete
      if (isTRUE(use_polygons)) {
        poly_df[[".color_value"]] <- as.factor(poly_df[[".color_value"]])
        
        p <- p +
          ggplot2::geom_polygon(
            ggplot2::aes(fill = .data[[".color_value"]]),
            color = polygon_outline,
            alpha = polygon_alpha,
            linewidth = 0.05
          ) +
          ggplot2::guides(fill = ggplot2::guide_legend(
            title = legend_title,
            override.aes = list(shape = 15, size = legend_square_size)
          ))
        
        if (!is.null(palette)) {
          p <- p + ggplot2::scale_fill_manual(values = palette, name = legend_title)
        } else {
          p <- p + ggplot2::scale_fill_discrete(name = legend_title)
        }
        
      } else {
        df[[".color_value"]] <- as.factor(df[[".color_value"]])
        
        p <- p +
          ggplot2::geom_point(
            ggplot2::aes(color = .data[[".color_value"]]),
            size = pt_size
          ) +
          ggplot2::guides(color = ggplot2::guide_legend(
            title = legend_title,
            override.aes = list(shape = 15, size = legend_square_size)
          ))
        
        if (!is.null(palette)) {
          p <- p + ggplot2::scale_color_manual(values = palette, name = legend_title)
        } else {
          p <- p + ggplot2::scale_color_discrete(name = legend_title)
        }
      }
    }
    
  } else {
    # no group_by
    if (isTRUE(use_polygons)) {
      p <- p + ggplot2::geom_polygon(
        fill = "white",
        color = polygon_outline,
        alpha = polygon_alpha,
        linewidth = 0.05
      )
    } else {
      p <- p + ggplot2::geom_point(size = pt_size, color = "white")
    }
  }
  
  # -------------------------
  # Coordinate handling
  # -------------------------
  p <- p + ggplot2::coord_fixed()
  if (isTRUE(reverse_y)) p <- p + ggplot2::scale_y_reverse()
  
  # -------------------------
  # Black background theme
  # -------------------------
  p <- p +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill = "black", color = NA),
      plot.background   = ggplot2::element_rect(fill = "black", color = NA),
      panel.grid        = ggplot2::element_blank(),
      axis.text         = ggplot2::element_text(color = "white"),
      axis.title        = ggplot2::element_text(color = "white"),
      axis.line         = ggplot2::element_line(color = "white"),
      axis.ticks        = ggplot2::element_line(color = "white"),
      legend.background = ggplot2::element_rect(fill = "black", color = NA),
      legend.key        = ggplot2::element_rect(fill = "black", color = NA),
      legend.text       = ggplot2::element_text(color = "white"),
      legend.title      = ggplot2::element_text(color = "white")
    )
  
  p
}

plotNodeMulti <- function(effs,
                          group_by = NULL,
                          combine = FALSE,
                          ncol = NULL,
                          plot_titles = names(effs),
                          ...) {
  
  # ----------------------------
  # Validate eff input
  # ----------------------------
  if (inherits(effs, "EffNICHES")) {
    effs <- list(effs)
  }
  
  if (!is.list(effs) || length(effs) == 0) {
    stop("'effs' must be a non-empty list of EffNICHES objects.")
  }
  
  if (is.null(names(effs))) {
    names(effs) <- paste0("eff_", seq_along(effs))
  } else {
    blank_names <- is.na(names(effs)) | names(effs) == ""
    names(effs)[blank_names] <- paste0("eff_", which(blank_names))
  }
  
  if (length(plot_titles) != length(effs)) {
    stop("'plot_titles' must have the same length as 'effs'.")
  }
  
  # ----------------------------
  # Helpers copied from plotNode
  # ----------------------------
  .is_matlike <- function(x) {
    is.matrix(x) || inherits(x, "dgCMatrix") || inherits(x, "Matrix")
  }
  
  .candidate_mats <- function(eff) {
    out <- list()
    
    paths <- list(
      expr            = eff$expr,
      expr_mat        = eff$expr$mat,
      expression      = eff$expression,
      expression_mat  = eff$expression$mat,
      counts          = eff$counts,
      data            = eff$data,
      X               = eff$X,
      mat             = eff$mat
    )
    
    for (nm in names(paths)) {
      x <- paths[[nm]]
      if (.is_matlike(x)) out[[nm]] <- x
    }
    out
  }
  
  .pick_expr_mat <- function(eff, ids) {
    mats <- .candidate_mats(eff)
    if (length(mats) == 0) return(NULL)
    
    best <- NULL
    best_score <- -Inf
    for (nm in names(mats)) {
      m <- mats[[nm]]
      cn <- colnames(m)
      rn <- rownames(m)
      if (is.null(cn) || is.null(rn)) next
      
      score <- sum(ids %in% cn)
      if (score > best_score) {
        best <- m
        best_score <- score
      }
    }
    best
  }
  
  .fetch_values <- function(eff, group_by) {
    if (is.null(group_by)) return(NULL)
    
    if (is.null(eff$nodes$xy) || ncol(eff$nodes$xy) != 2) {
      stop("Each eff must have eff$nodes$xy with 2 columns.")
    }
    ids <- eff$nodes$ids
    if (is.null(ids) || length(ids) != nrow(eff$nodes$xy)) {
      stop("Each eff must have eff$nodes$ids matching nrow(eff$nodes$xy).")
    }
    
    df <- data.frame(id = ids, stringsAsFactors = FALSE)
    
    if (!is.null(eff$nodes$meta)) {
      meta <- eff$nodes$meta
      if (!is.null(rownames(meta))) {
        o <- match(ids, rownames(meta))
        if (anyNA(o)) {
          if (nrow(meta) != length(ids)) {
            stop("eff$nodes$meta rownames do not cover eff$nodes$ids and nrow(meta) != N.")
          }
          df <- cbind(df, meta)
        } else {
          df <- cbind(df, meta[o, , drop = FALSE])
        }
      } else {
        if (nrow(meta) != length(ids)) {
          stop("eff$nodes$meta has no rownames and nrow(meta) != N; cannot align.")
        }
        df <- cbind(df, meta)
      }
    }
    
    if (is.numeric(group_by)) {
      if (length(group_by) != length(ids)) {
        stop("If group_by is a numeric vector, it must have length N = length(eff$nodes$ids).")
      }
      return(as.numeric(group_by))
    }
    
    if (!is.character(group_by) || length(group_by) != 1L) {
      stop("group_by must be NULL, a single column/feature name, or a numeric vector of length N.")
    }
    
    if (group_by %in% colnames(df)) {
      return(df[[group_by]])
    }
    
    m <- .pick_expr_mat(eff, ids)
    if (is.null(m)) {
      stop(
        "group_by not found in eff$nodes$meta (as a column), and no expression matrix was detected in eff.\n",
        "Tried: eff$expr, eff$expr$mat, eff$expression, eff$expression$mat, eff$counts, eff$data, eff$X, eff$mat."
      )
    }
    if (is.null(rownames(m))) stop("Detected expression matrix has no rownames (feature names).")
    if (is.null(colnames(m))) stop("Detected expression matrix has no colnames (cell/node ids).")
    
    if (!(group_by %in% rownames(m))) {
      stop("group_by not found as meta column or feature in detected expression matrix: ", group_by)
    }
    
    ci <- match(ids, colnames(m))
    if (anyNA(ci)) {
      stop("Expression matrix colnames do not cover all eff$nodes$ids; cannot align expression to nodes.")
    }
    
    as.numeric(m[group_by, ci])
  }
  
  .apply_transform <- function(values, log_transform) {
    if (is.null(values)) return(NULL)
    if (!is.numeric(values)) return(values)
    
    if (isTRUE(log_transform)) {
      return(log1p(values))
    } else if (is.function(log_transform)) {
      return(log_transform(values))
    } else if (identical(log_transform, FALSE)) {
      return(values)
    } else {
      stop("log_transform must be FALSE, TRUE (log1p), or a function.")
    }
  }
  
  # extract log_transform if supplied in ...
  dots <- list(...)
  log_transform <- if ("log_transform" %in% names(dots)) dots$log_transform else FALSE
  
  # ----------------------------
  # Determine whether shared quantitative scale is needed
  # ----------------------------
  global_limits <- NULL
  is_quantitative_global <- FALSE
  
  if (!is.null(group_by)) {
    vals1 <- .fetch_values(effs[[1]], group_by)
    vals1 <- .apply_transform(vals1, log_transform)
    is_quantitative_global <- is.numeric(vals1)
    
    if (is_quantitative_global) {
      global_min <- Inf
      global_max <- -Inf
      
      for (k in seq_along(effs)) {
        eff <- effs[[k]]
        if (!inherits(eff, "EffNICHES")) {
          stop("All elements of 'effs' must inherit from 'EffNICHES'. Problem at element ", k, ".")
        }
        
        vals <- .fetch_values(eff, group_by)
        vals <- .apply_transform(vals, log_transform)
        vals <- vals[is.finite(vals)]
        
        if (length(vals) == 0) {
          stop("No finite quantitative values found in effs[[", k, "]].")
        }
        
        global_min <- min(global_min, min(vals, na.rm = TRUE))
        global_max <- max(global_max, max(vals, na.rm = TRUE))
      }
      
      global_limits <- c(global_min, global_max)
      if (global_limits[1] == global_limits[2]) {
        eps <- max(abs(global_limits[1]), 1) * 1e-9
        global_limits <- global_limits + c(-eps, eps)
      }
    }
  }
  
  # ----------------------------
  # Build plots
  # ----------------------------
  plots <- vector("list", length(effs))
  names(plots) <- names(effs)
  
  for (k in seq_along(effs)) {
    eff <- effs[[k]]
    
    if (!inherits(eff, "EffNICHES")) {
      stop("All elements of 'effs' must inherit from 'EffNICHES'. Problem at element ", k, ".")
    }
    
    call_args <- c(
      list(
        eff = eff,
        group_by = group_by
      ),
      dots
    )
    
    if (is_quantitative_global) {
      call_args$color_limits <- global_limits
    }
    
    p <- do.call(plotNode, call_args)
    
    if (!is.null(plot_titles)) {
      p <- p + ggplot2::ggtitle(plot_titles[k])
    }
    
    plots[[k]] <- p
  }
  
  attr(plots, "color_limits") <- global_limits
  
  if (isTRUE(combine)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' is required when combine = TRUE.")
    }
    return(
      patchwork::wrap_plots(plots, ncol = ncol) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "right")
    )
  }
  
  plots
}

#### Plot Network Density ####

plotNetworkDensity <- function(eff,
                               set.rad = NULL,
                               k = NULL,
                               pt_size = 0.3,
                               edge_width = 0.15,
                               edge_alpha = 0.33,
                               reverse_y = FALSE,
                               max_edges = 2e6,
                               report = TRUE) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  # --- enforce exactly one of set.rad or k ---
  if (is.null(set.rad) && is.null(k)) {
    stop("Specify exactly one of set.rad (radius) or k (kNN). Both are NULL.")
  }
  if (!is.null(set.rad) && !is.null(k)) {
    stop("Specify exactly one of set.rad (radius) or k (kNN). Both were provided.")
  }
  
  # --- pull xy ---
  xy <- eff$nodes$xy
  if (is.null(xy) || !is.matrix(xy) || ncol(xy) != 2) {
    stop("eff$nodes$xy must be an N x 2 matrix.")
  }
  
  ids <- eff$nodes$ids
  if (is.null(ids) || length(ids) != nrow(xy)) {
    stop("eff$nodes$ids must exist and have length equal to nrow(eff$nodes$xy).")
  }
  
  x <- suppressWarnings(as.numeric(xy[, 1]))
  y <- suppressWarnings(as.numeric(xy[, 2]))
  if (anyNA(x) || anyNA(y)) stop("NA detected in x/y after numeric coercion.")
  
  n <- length(ids)
  coords <- cbind(x, y)
  
  df <- data.frame(
    id = ids,
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )
  
  # -------------------------
  # Edge builders return integer index pairs only
  # -------------------------
  
  make_edges_radius_idx <- function(coords, rad, max_edges) {
    rad <- as.numeric(rad)
    if (!is.finite(rad) || rad <= 0) {
      stop("set.rad must be a positive finite number.")
    }
    
    # Fast and memory-safe radius search
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      stop(
        "Radius mode now requires the 'dbscan' package for memory-safe neighbor search.\n",
        "Install it with: install.packages('dbscan')"
      )
    }
    
    fr <- dbscan::frNN(coords, eps = rad, sort = FALSE)
    
    # Count unique undirected edges before allocating
    edge_n <- 0L
    for (i in seq_len(n)) {
      nei <- fr$id[[i]]
      if (length(nei)) edge_n <- edge_n + sum(nei > i)
      if (edge_n > max_edges) {
        stop("Edge count exceeds max_edges (", max_edges, "). ",
             "Use a smaller radius or increase max_edges.")
      }
    }
    
    if (edge_n == 0L) {
      return(list(a = integer(0), b = integer(0)))
    }
    
    a <- integer(edge_n)
    b <- integer(edge_n)
    pos <- 1L
    
    for (i in seq_len(n)) {
      nei <- fr$id[[i]]
      if (!length(nei)) next
      nei <- nei[nei > i]   # unique undirected edges only
      m <- length(nei)
      if (!m) next
      
      rng <- pos:(pos + m - 1L)
      a[rng] <- i
      b[rng] <- nei
      pos <- pos + m
    }
    
    list(a = a, b = b)
  }
  
  make_edges_knn_idx <- function(coords, k, max_edges) {
    k <- as.integer(k)
    if (!is.finite(k) || k < 1L) stop("k must be a positive integer.")
    if (k >= n) stop("k must be < number of nodes (", n, ").")
    
    if (!requireNamespace("FNN", quietly = TRUE)) {
      stop(
        "kNN mode requires the 'FNN' package.\n",
        "Install it with: install.packages('FNN')"
      )
    }
    
    idx <- FNN::get.knn(coords, k = k)$nn.index
    
    i_rep <- rep.int(seq_len(n), times = k)
    j_rep <- as.vector(idx)
    
    keep <- !is.na(j_rep) & i_rep != j_rep
    i_rep <- i_rep[keep]
    j_rep <- j_rep[keep]
    
    a <- pmin.int(i_rep, j_rep)
    b <- pmax.int(i_rep, j_rep)
    
    # deduplicate unordered pairs without paste()
    ord <- order(a, b)
    a <- a[ord]
    b <- b[ord]
    
    keep_unique <- c(TRUE, (a[-1L] != a[-length(a)]) | (b[-1L] != b[-length(b)]))
    a <- a[keep_unique]
    b <- b[keep_unique]
    
    if (length(a) > max_edges) {
      stop("Edge count (", length(a), ") exceeds max_edges (", max_edges, "). ",
           "Reduce k or increase max_edges.")
    }
    
    list(a = a, b = b)
  }
  
  edge_idx <- if (!is.null(set.rad)) {
    make_edges_radius_idx(coords, set.rad, max_edges)
  } else {
    make_edges_knn_idx(coords, k, max_edges)
  }
  
  a <- edge_idx$a
  b <- edge_idx$b
  
  edges <- data.frame(
    from = ids[a],
    to   = ids[b],
    x1   = x[a],
    y1   = y[a],
    x2   = x[b],
    y2   = y[b],
    stringsAsFactors = FALSE
  )
  
  # -------------------------
  # Report
  # -------------------------
  report_text <- NULL
  if (isTRUE(report)) {
    method_line <- if (!is.null(set.rad)) {
      paste0("Radius = ", format(set.rad, digits = 6))
    } else {
      paste0("k = ", as.integer(k))
    }
    
    deg <- integer(n)
    if (length(a)) {
      deg <- deg + tabulate(a, nbins = n)
      deg <- deg + tabulate(b, nbins = n)
    }
    
    qs <- as.numeric(stats::quantile(
      deg,
      probs = c(0, 0.25, 0.5, 0.75, 1),
      names = FALSE,
      type = 7
    ))
    
    report_text <- paste(
      method_line,
      paste0("Edges drawn = ", format(length(a), big.mark = ",")),
      paste0("Min / Max = ", qs[1], " / ", qs[5]),
      paste0("Q1 / Q2 / Q3 = ", qs[2], " / ", qs[3], " / ", qs[4]),
      paste0("Mean = ", round(mean(deg), 3)),
      sep = "\n"
    )
  }
  
  # -------------------------
  # Plot
  # -------------------------
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x1, y = .data$y1, xend = .data$x2, yend = .data$y2),
      linewidth = edge_width,
      alpha = edge_alpha,
      color = "black"
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data$x, y = .data$y),
      size = pt_size,
      color = "black"
    ) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 5.5, r = 90, b = 5.5, l = 5.5, unit = "pt")
    )
  
  if (isTRUE(reverse_y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  
  if (isTRUE(report) && !is.null(report_text)) {
    xr <- range(x, na.rm = TRUE)
    xpad <- if (diff(xr) == 0) 1 else 0.05 * diff(xr)
    
    y_anchor <- if (isTRUE(reverse_y)) min(y, na.rm = TRUE) else max(y, na.rm = TRUE)
    x_anchor <- max(x, na.rm = TRUE) + xpad
    
    p <- p + ggplot2::annotate(
      "text",
      x = x_anchor,
      y = y_anchor,
      label = report_text,
      hjust = 0, vjust = 1,
      size = 3.3,
      family = "mono"
    )
  }
  
  p
}

#### LR Edge Plot ####

#' Plot LR network in spatial tissue coordinates
plotLRInSitu <- function(eff,
                         lr,
                         transform = identity,
                         # visual controls
                         edge_width_bg = 0.10,
                         edge_alpha_bg = 0.075,
                         edge_color_bg = "white",
                         edge_width_fg = 0.20,
                         pt_size = 0.0,
                         pt_color = "white",
                         reverse_y = FALSE,
                         # scaling controls
                         alpha_range = c(0.05, 1),
                         color_option = "C",
                         color_limits = NULL,
                         oob_squish = TRUE,
                         # LR combination controls
                         lr_normalization = c("none", "max", "percentile", "zscore"),
                         lr_aggregation   = c("sum", "mean", "max"),
                         deduplicate_lr   = TRUE,
                         warn_many_lr     = TRUE,
                         plot_title       = NULL,
                         # visualization filtering
                         edge_fraction = 1,
                         edge_top = NULL,
                         # safety
                         max_edges = 2e6) {
  
  .plotLRLayout(
    eff = eff,
    lr = lr,
    layout = "spatial",
    coords = NULL,
    transform = transform,
    edge_width_bg = edge_width_bg,
    edge_alpha_bg = edge_alpha_bg,
    edge_color_bg = edge_color_bg,
    edge_width_fg = edge_width_fg,
    pt_size = pt_size,
    pt_color = pt_color,
    reverse_y = reverse_y,
    alpha_range = alpha_range,
    color_option = color_option,
    color_limits = color_limits,
    oob_squish = oob_squish,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    warn_many_lr = warn_many_lr,
    plot_title = plot_title,
    title_suffix = NULL,
    edge_fraction = edge_fraction,
    edge_top = edge_top,
    max_edges = max_edges
  )
}

#' Plot LR network in UMAP coordinates
plotLRNetwork <- function(eff,
                          lr,
                          umap_coords = NULL,
                          transform = identity,
                          # visual controls
                          edge_width_bg = 0.10,
                          edge_alpha_bg = 0.075,
                          edge_color_bg = "white",
                          edge_width_fg = 0.20,
                          pt_size = 0.20,
                          pt_color = "grey55",
                          reverse_y = FALSE,
                          # scaling controls
                          alpha_range = c(0.05, 1),
                          color_option = "C",
                          color_limits = NULL,
                          oob_squish = TRUE,
                          # LR combination controls
                          lr_normalization = c("none", "max", "percentile", "zscore"),
                          lr_aggregation   = c("sum", "mean", "max"),
                          deduplicate_lr   = TRUE,
                          warn_many_lr     = TRUE,
                          plot_title       = NULL,
                          # visualization filtering
                          edge_fraction = 1,
                          edge_top = NULL,
                          # safety
                          max_edges = 2e6) {
  
  .plotLRLayout(
    eff = eff,
    lr = lr,
    layout = "umap",
    coords = umap_coords,
    transform = transform,
    edge_width_bg = edge_width_bg,
    edge_alpha_bg = edge_alpha_bg,
    edge_color_bg = edge_color_bg,
    edge_width_fg = edge_width_fg,
    pt_size = pt_size,
    pt_color = pt_color,
    reverse_y = reverse_y,
    alpha_range = alpha_range,
    color_option = color_option,
    color_limits = color_limits,
    oob_squish = oob_squish,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    warn_many_lr = warn_many_lr,
    plot_title = plot_title,
    title_suffix = NULL,
    edge_fraction = edge_fraction,
    edge_top = edge_top,
    max_edges = max_edges
  )
}

#' Plot multiple LR networks in spatial tissue coordinates
plotLRInSituMulti <- function(eff,
                              lr,
                              plot_titles = NULL,
                              ncol = NULL,
                              nrow = NULL,
                              guides = "collect",
                              share_color_limits = TRUE,
                              transform = identity,
                              # visual controls
                              edge_width_bg = 0.10,
                              edge_alpha_bg = 0.075,
                              edge_color_bg = "white",
                              edge_width_fg = 0.20,
                              pt_size = 0.0,
                              pt_color = "white",
                              reverse_y = FALSE,
                              # scaling controls
                              alpha_range = c(0.05, 1),
                              color_option = "C",
                              color_limits = NULL,
                              oob_squish = TRUE,
                              # LR combination controls
                              lr_normalization = c("none", "max", "percentile", "zscore"),
                              lr_aggregation   = c("sum", "mean", "max"),
                              deduplicate_lr   = TRUE,
                              warn_many_lr     = TRUE,
                              # visualization filtering
                              edge_fraction = 1,
                              edge_top = NULL,
                              # safety
                              max_edges = 2e6) {
  if (missing(lr) || is.null(lr) || length(lr) == 0) {
    stop("lr must be a non-empty character vector or list.")
  }
  
  if (is.character(lr)) {
    lr_list <- as.list(lr)
    lr_labels <- lr
  } else if (is.list(lr)) {
    lr_list <- lr
    lr_labels <- vapply(
      lr_list,
      function(x) {
        if (length(x) == 1L) {
          as.character(x)
        } else {
          paste(as.character(x), collapse = ", ")
        }
      },
      character(1)
    )
  } else {
    stop("lr must be either a character vector or a list.")
  }
  
  n_panels <- length(lr_list)
  eff_list <- rep(list(eff), n_panels)
  
  if (is.null(plot_titles)) {
    plot_titles <- lr_labels
  } else {
    if (length(plot_titles) != n_panels) {
      stop("plot_titles must have the same length as lr.")
    }
    plot_titles <- as.character(plot_titles)
  }
  
  .plotLRPanelGrid(
    eff_list = eff_list,
    lr_list = lr_list,
    layout = "spatial",
    coords_list = rep(list(NULL), n_panels),
    plot_titles = plot_titles,
    ncol = ncol,
    nrow = nrow,
    guides = guides,
    share_color_limits = share_color_limits,
    show_axes = TRUE,
    transform = transform,
    edge_width_bg = edge_width_bg,
    edge_alpha_bg = edge_alpha_bg,
    edge_color_bg = edge_color_bg,
    edge_width_fg = edge_width_fg,
    pt_size = pt_size,
    pt_color = pt_color,
    reverse_y = reverse_y,
    alpha_range = alpha_range,
    color_option = color_option,
    color_limits = color_limits,
    oob_squish = oob_squish,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    warn_many_lr = warn_many_lr,
    title_suffix = NULL,
    edge_fraction = edge_fraction,
    edge_top = edge_top,
    max_edges = max_edges
  )
}

#' Plot multiple LR networks in UMAP coordinates
plotLRNetworkMulti <- function(eff,
                               lr,
                               umap_coords = NULL,
                               plot_titles = NULL,
                               ncol = NULL,
                               nrow = NULL,
                               guides = "collect",
                               share_color_limits = TRUE,
                               transform = identity,
                               # visual controls
                               edge_width_bg = 0.10,
                               edge_alpha_bg = 0.075,
                               edge_color_bg = "white",
                               edge_width_fg = 0.20,
                               pt_size = 0.20,
                               pt_color = "grey55",
                               reverse_y = FALSE,
                               # scaling controls
                               alpha_range = c(0.05, 1),
                               color_option = "C",
                               color_limits = NULL,
                               oob_squish = TRUE,
                               # LR combination controls
                               lr_normalization = c("none", "max", "percentile", "zscore"),
                               lr_aggregation   = c("sum", "mean", "max"),
                               deduplicate_lr   = TRUE,
                               warn_many_lr     = TRUE,
                               # visualization filtering
                               edge_fraction = 1,
                               edge_top = NULL,
                               # safety
                               max_edges = 2e6) {
  if (missing(lr) || is.null(lr) || length(lr) == 0) {
    stop("lr must be a non-empty character vector or list.")
  }
  
  if (is.character(lr)) {
    lr_list <- as.list(lr)
    lr_labels <- lr
  } else if (is.list(lr)) {
    lr_list <- lr
    lr_labels <- vapply(
      lr_list,
      function(x) {
        if (length(x) == 1L) {
          as.character(x)
        } else {
          paste(as.character(x), collapse = ", ")
        }
      },
      character(1)
    )
  } else {
    stop("lr must be either a character vector or a list.")
  }
  
  n_panels <- length(lr_list)
  eff_list <- rep(list(eff), n_panels)
  coords_list <- rep(list(umap_coords), n_panels)
  
  if (is.null(plot_titles)) {
    plot_titles <- lr_labels
  } else {
    if (length(plot_titles) != n_panels) {
      stop("plot_titles must have the same length as lr.")
    }
    plot_titles <- as.character(plot_titles)
  }
  
  .plotLRPanelGrid(
    eff_list = eff_list,
    lr_list = lr_list,
    layout = "umap",
    coords_list = coords_list,
    plot_titles = plot_titles,
    ncol = ncol,
    nrow = nrow,
    guides = guides,
    share_color_limits = share_color_limits,
    show_axes = TRUE,
    transform = transform,
    edge_width_bg = edge_width_bg,
    edge_alpha_bg = edge_alpha_bg,
    edge_color_bg = edge_color_bg,
    edge_width_fg = edge_width_fg,
    pt_size = pt_size,
    pt_color = pt_color,
    reverse_y = reverse_y,
    alpha_range = alpha_range,
    color_option = color_option,
    color_limits = color_limits,
    oob_squish = oob_squish,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    warn_many_lr = warn_many_lr,
    title_suffix = NULL,
    edge_fraction = edge_fraction,
    edge_top = edge_top,
    max_edges = max_edges
  )
}

plotLRInSituAcrossObjects <- function(
    eff_list,
    lr,
    plot_titles = names(eff_list),
    ncol = 4,
    nrow = 2,
    transform = identity,
    # visual controls
    edge_width = NULL,
    edge_width_bg = 0.10,
    edge_alpha_bg = 0.075,
    edge_color_bg = "white",
    edge_width_fg = 0.20,
    pt_size = 0.0,
    pt_color = "white",
    reverse_y = FALSE,
    # scaling controls
    alpha_range = c(0.05, 1),
    color_option = "C",
    color_limits = NULL,
    oob_squish = TRUE,
    # LR combination controls
    lr_normalization = c("none", "max", "percentile", "zscore"),
    lr_aggregation   = c("sum", "mean", "max"),
    deduplicate_lr   = TRUE,
    warn_many_lr     = TRUE,
    # visualization filtering
    edge_fraction = 1,
    edge_top = NULL,
    # safety
    max_edges = 2e6,
    guides = "collect",
    share_color_limits = TRUE,
    share_alpha_limits = TRUE
) {
  if (!is.list(eff_list) || length(eff_list) == 0L ||
      !all(vapply(eff_list, inherits, logical(1), "EffNICHES"))) {
    stop("All elements of eff_list must inherit from 'EffNICHES'.")
  }
  
  if (length(plot_titles) != length(eff_list)) {
    stop("plot_titles must match length of eff_list.")
  }
  
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  
  # Backward-compatible alias:
  # if user supplies edge_width, use it for the foreground edges
  if (!is.null(edge_width)) {
    edge_width_fg <- edge_width
  }
  
  shared_panel <- computeSharedSpatialPanelSize(eff_list)
  
  panel_data_list <- vector("list", length(eff_list))
  plot_list <- vector("list", length(eff_list))
  
  for (i in seq_along(eff_list)) {
    panel_data_list[[i]] <- prepareLRPanelData(
      eff = eff_list[[i]],
      lr = lr,
      layout = "spatial",
      coords = NULL,
      transform = transform,
      lr_normalization = lr_normalization,
      lr_aggregation = lr_aggregation,
      deduplicate_lr = deduplicate_lr,
      warn_many_lr = warn_many_lr,
      plot_title = NULL,
      title_suffix = NULL,
      edge_fraction = edge_fraction,
      edge_top = edge_top,
      max_edges = max_edges
    )
  }
  
  if (isTRUE(share_color_limits)) {
    color_limits_use <- computeSharedLRColorLimits(
      panel_data_list = panel_data_list,
      color_limits = color_limits
    )
  } else {
    color_limits_use <- color_limits
  }
  
  if (isTRUE(share_alpha_limits)) {
    alpha_limits_use <- computeSharedLRAlphaLimits(
      panel_data_list = panel_data_list,
      alpha_limits = NULL
    )
  } else {
    alpha_limits_use <- NULL
  }
  
  for (i in seq_along(eff_list)) {
    lims_i <- computeCenteredSpatialLimits(
      eff = eff_list[[i]],
      shared_width = shared_panel$width,
      shared_height = shared_panel$height
    )
    
    plot_list[[i]] <- renderPreparedLRPanel(
      panel_data = panel_data_list[[i]],
      plot_title = as.character(plot_titles[i]),
      show_axes = FALSE,
      reverse_y = reverse_y,
      edge_width_bg = edge_width_bg,
      edge_alpha_bg = edge_alpha_bg,
      edge_color_bg = edge_color_bg,
      edge_width_fg = edge_width_fg,
      pt_size = pt_size,
      pt_color = pt_color,
      alpha_range = alpha_range,
      alpha_limits = alpha_limits_use,
      color_option = color_option,
      color_limits = color_limits_use,
      oob_squish = oob_squish,
      x_limits = lims_i$x_limits,
      y_limits = lims_i$y_limits
    )
  }
  
  patchwork::wrap_plots(
    plot_list,
    ncol = ncol,
    nrow = nrow,
    guides = guides
  ) &
    ggplot2::theme(
      legend.position = "right",
      plot.background  = ggplot2::element_rect(fill = "black", colour = NA),
      panel.background = ggplot2::element_rect(fill = "black", colour = NA)
    )
}

#' Plot LR network in spatial tissue coordinates with sender/receiver identity coloring
plotLRInSituTyped <- function(eff,
                              lr,
                              group_by,
                              transform = identity,
                              # identity colors
                              palette = NULL,
                              legend_types = NULL,
                              # edge identity filtering
                              cell_types = NULL,
                              cell_types_mode = c("either", "both"),
                              sender_types = NULL,
                              receiver_types = NULL,
                              # contextual scaffold
                              show_edge_scaffold = TRUE,
                              scaffold_color = "grey20",
                              scaffold_alpha = 0.35,
                              scaffold_width = 0.05,
                              # rendering
                              edge_gradient = FALSE,
                              edge_width = 0.20,
                              edge_n_segments = 9L,
                              pt_size = 0.0,
                              pt_color = "white",
                              reverse_y = FALSE,
                              # alpha controls
                              alpha_range = c(0.0, 1),
                              alpha_limits = NULL,
                              # shared spatial limits
                              x_limits = NULL,
                              y_limits = NULL,
                              # LR combination controls
                              lr_normalization = c("none", "max", "percentile", "zscore"),
                              lr_aggregation   = c("sum", "mean", "max"),
                              deduplicate_lr   = TRUE,
                              warn_many_lr     = TRUE,
                              plot_title       = NULL,
                              # visualization filtering
                              edge_fraction = 1,
                              edge_top = NULL,
                              # safety
                              max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  cell_types_mode  <- match.arg(cell_types_mode)
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  
  if (!is.logical(show_edge_scaffold) || length(show_edge_scaffold) != 1L || is.na(show_edge_scaffold)) {
    stop("show_edge_scaffold must be TRUE or FALSE.")
  }
  
  if (!is.logical(edge_gradient) || length(edge_gradient) != 1L || is.na(edge_gradient)) {
    stop("edge_gradient must be TRUE or FALSE.")
  }
  
  if (!is.numeric(edge_width) || length(edge_width) != 1L ||
      !is.finite(edge_width) || edge_width <= 0) {
    stop("edge_width must be a single positive number.")
  }
  
  if (is.null(scaffold_width)) {
    scaffold_width <- edge_width
  }
  if (!is.numeric(scaffold_width) || length(scaffold_width) != 1L ||
      !is.finite(scaffold_width) || scaffold_width <= 0) {
    stop("scaffold_width must be NULL or a single positive number.")
  }
  
  if (!is.numeric(scaffold_alpha) || length(scaffold_alpha) != 1L ||
      !is.finite(scaffold_alpha) || scaffold_alpha < 0 || scaffold_alpha > 1) {
    stop("scaffold_alpha must be a single number in [0, 1].")
  }
  
  edge_n_segments <- as.integer(edge_n_segments)
  if (!is.finite(edge_n_segments) || edge_n_segments < 3L) {
    stop("edge_n_segments must be an integer >= 3.")
  }
  
  if (!is.numeric(edge_fraction) || length(edge_fraction) != 1L ||
      !is.finite(edge_fraction) || edge_fraction <= 0 || edge_fraction > 1) {
    stop("edge_fraction must be a single number in (0, 1].")
  }
  
  if (!is.null(edge_top)) {
    if (!is.numeric(edge_top) || length(edge_top) != 1L ||
        !is.finite(edge_top) || edge_top < 1) {
      stop("edge_top must be NULL or a single positive number.")
    }
    edge_top <- as.integer(edge_top)
  }
  
  if (!is.null(x_limits)) {
    if (!is.numeric(x_limits) || length(x_limits) != 2L || any(!is.finite(x_limits)) || x_limits[1] >= x_limits[2]) {
      stop("x_limits must be NULL or a numeric vector of length 2 with x_limits[1] < x_limits[2].")
    }
  }
  
  if (!is.null(y_limits)) {
    if (!is.numeric(y_limits) || length(y_limits) != 2L || any(!is.finite(y_limits)) || y_limits[1] >= y_limits[2]) {
      stop("y_limits must be NULL or a numeric vector of length 2 with y_limits[1] < y_limits[2].")
    }
  }
  
  if (!is.null(alpha_limits)) {
    if (!is.numeric(alpha_limits) || length(alpha_limits) != 2L || any(!is.finite(alpha_limits)) || alpha_limits[1] > alpha_limits[2]) {
      stop("alpha_limits must be NULL or a numeric vector of length 2 with alpha_limits[1] <= alpha_limits[2].")
    }
  }
  
  if (is.null(eff$niches$CellToCellSpatial)) {
    stop("eff$niches$CellToCellSpatial is missing.")
  }
  
  c2c <- eff$niches$CellToCellSpatial
  ij  <- c2c$ij
  
  if (is.null(ij) || !is.matrix(ij) || ncol(ij) != 2) {
    stop("CellToCellSpatial$ij must be a matrix with 2 columns (node indices).")
  }
  if (nrow(ij) > max_edges) {
    stop("Edge count (", nrow(ij), ") exceeds max_edges (", max_edges, ").")
  }
  
  coords_out <- resolveEffLayoutCoords(
    eff = eff,
    layout = "spatial",
    coords = NULL
  )
  
  lr_out <- computeEdgeLRsignal(
    eff = eff,
    lr = lr,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    transform = transform,
    warn_many_lr = warn_many_lr,
    max_edges = max_edges
  )
  
  node_groups <- getEffNodeGrouping(eff, group_by = group_by)
  xy <- coords_out$xy
  
  x_limits_use <- if (is.null(x_limits)) range(xy[, "x"], na.rm = TRUE) else x_limits
  y_limits_use <- if (is.null(y_limits)) range(xy[, "y"], na.rm = TRUE) else y_limits
  
  from_idx <- as.integer(ij[, 1])
  to_idx   <- as.integer(ij[, 2])
  
  if (anyNA(from_idx) || anyNA(to_idx)) {
    stop("NA detected in CellToCellSpatial$ij.")
  }
  if (any(from_idx < 1L) || any(from_idx > nrow(xy)) ||
      any(to_idx   < 1L) || any(to_idx   > nrow(xy))) {
    stop("CellToCellSpatial$ij contains indices out of bounds.")
  }
  
  all_edge_df <- data.frame(
    x1 = as.numeric(xy[from_idx, "x"]),
    y1 = as.numeric(xy[from_idx, "y"]),
    x2 = as.numeric(xy[to_idx,   "x"]),
    y2 = as.numeric(xy[to_idx,   "y"]),
    for.plotting = as.numeric(lr_out$edge_signal),
    sender_type = node_groups[from_idx],
    receiver_type = node_groups[to_idx],
    stringsAsFactors = FALSE
  )
  
  normalize_type_filter <- function(x, arg_name) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- unique(trimws(x))
    x <- x[nzchar(x)]
    if (length(x) < 1L) {
      stop(arg_name, " must contain at least one non-empty cell type or be NULL.")
    }
    x
  }
  
  cell_types     <- normalize_type_filter(cell_types, "cell_types")
  sender_types   <- normalize_type_filter(sender_types, "sender_types")
  receiver_types <- normalize_type_filter(receiver_types, "receiver_types")
  
  role_keep <- rep(TRUE, nrow(all_edge_df))
  
  if (!is.null(cell_types)) {
    if (cell_types_mode == "either") {
      role_keep <- role_keep & (
        all_edge_df$sender_type %in% cell_types |
          all_edge_df$receiver_type %in% cell_types
      )
    } else {
      role_keep <- role_keep & (
        all_edge_df$sender_type %in% cell_types &
          all_edge_df$receiver_type %in% cell_types
      )
    }
  }
  
  if (!is.null(sender_types)) {
    role_keep <- role_keep & (all_edge_df$sender_type %in% sender_types)
  }
  
  if (!is.null(receiver_types)) {
    role_keep <- role_keep & (all_edge_df$receiver_type %in% receiver_types)
  }
  
  scaffold_df <- all_edge_df[
    ,
    c("x1", "y1", "x2", "y2", "sender_type", "receiver_type"),
    drop = FALSE
  ]
  
  edge_df <- all_edge_df[role_keep, , drop = FALSE]
  
  if (nrow(edge_df) > 0L) {
    keep_realized <- is.finite(edge_df$for.plotting) &
      !is.na(edge_df$for.plotting) &
      (edge_df$for.plotting > 0)
    
    edge_df <- edge_df[keep_realized, , drop = FALSE]
  } else {
    edge_df <- edge_df[0, , drop = FALSE]
  }
  
  if (nrow(edge_df) > 0L && (edge_fraction < 1 || !is.null(edge_top))) {
    ord_desc <- order(edge_df$for.plotting, decreasing = TRUE)
    n_keep_fraction <- ceiling(nrow(edge_df) * edge_fraction)
    n_keep_final <- n_keep_fraction
    
    if (!is.null(edge_top)) {
      n_keep_final <- min(n_keep_final, edge_top)
    }
    n_keep_final <- max(1L, min(nrow(edge_df), n_keep_final))
    
    edge_df <- edge_df[ord_desc[seq_len(n_keep_final)], , drop = FALSE]
  }
  
  has_realized_edges <- nrow(edge_df) > 0L
  
  if (has_realized_edges) {
    edge_df$alpha_val <- computeEdgeAlpha(
      edge_df$for.plotting,
      alpha_range = alpha_range,
      from = alpha_limits
    )
    
    edge_df <- edge_df[order(edge_df$for.plotting), , drop = FALSE]
    edge_df$edge_order <- seq_len(nrow(edge_df))
  }
  
  used_types_realized <- if (has_realized_edges) {
    sort(unique(c(edge_df$sender_type, edge_df$receiver_type)))
  } else {
    character(0)
  }
  
  used_types_scaffold <- if (nrow(scaffold_df) > 0L) {
    sort(unique(c(scaffold_df$sender_type, scaffold_df$receiver_type)))
  } else {
    character(0)
  }
  
  used_types <- sort(unique(c(used_types_realized, used_types_scaffold)))
  
  if (length(used_types) < 1L) {
    used_types <- sort(unique(as.character(node_groups[!is.na(node_groups)])))
    used_types <- used_types[nzchar(used_types)]
  }
  
  if (is.null(legend_types)) {
    legend_types_use <- used_types
  } else {
    legend_types_use <- as.character(legend_types)
    legend_types_use <- legend_types_use[!is.na(legend_types_use)]
    legend_types_use <- unique(trimws(legend_types_use))
    legend_types_use <- legend_types_use[nzchar(legend_types_use)]
    if (length(legend_types_use) < 1L) {
      stop("legend_types must contain at least one non-empty cell type or be NULL.")
    }
    if (!all(used_types %in% legend_types_use)) {
      stop("legend_types must include all drawn cell types in this panel.")
    }
    legend_types_use <- sort(legend_types_use)
  }
  
  if (is.null(palette)) {
    pal_vals <- scales::hue_pal()(length(legend_types_use))
    palette_map <- stats::setNames(pal_vals, legend_types_use)
  } else {
    if (!is.character(palette) || is.null(names(palette))) {
      stop("palette must be NULL or a named character vector keyed by cell type.")
    }
    if (anyNA(names(palette)) || any(!nzchar(names(palette)))) {
      stop("palette must have non-empty names for all supplied colors.")
    }
    
    missing_cols <- setdiff(legend_types_use, names(palette))
    if (length(missing_cols) > 0L) {
      stop(
        "palette is missing colors for legend cell types: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    palette_map <- palette[legend_types_use]
  }
  
  if (has_realized_edges) {
    draw_df <- buildTypedLREdgeDrawData(
      edge_df = edge_df,
      palette_map = palette_map,
      edge_gradient = edge_gradient,
      edge_n_segments = edge_n_segments
    )
  } else {
    draw_df <- data.frame(
      x = numeric(0),
      y = numeric(0),
      xend = numeric(0),
      yend = numeric(0),
      alpha_val = numeric(0),
      segment_color = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  title_use <- makeLRPlotTitle(
    selected_mechanisms = lr_out$selected_mechanisms,
    plot_title = plot_title,
    suffix = NULL
  )
  
  legend_df <- data.frame(
    cell_type = legend_types_use,
    x = rep(x_limits_use[1], length(legend_types_use)),
    y = rep(y_limits_use[1], length(legend_types_use)),
    stringsAsFactors = FALSE
  )
  
  p <- ggplot2::ggplot()
  
  if (isTRUE(pt_size > 0)) {
    node_df <- buildNodePlotData(xy)
    p <- p +
      ggplot2::geom_point(
        data = node_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        size = pt_size,
        color = pt_color,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  
  if (isTRUE(show_edge_scaffold) && nrow(scaffold_df) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = scaffold_df,
        ggplot2::aes(
          x = .data$x1, y = .data$y1,
          xend = .data$x2, yend = .data$y2
        ),
        inherit.aes = FALSE,
        color = scaffold_color,
        alpha = scaffold_alpha,
        linewidth = scaffold_width,
        lineend = "butt",
        show.legend = FALSE
      )
  }
  
  if (nrow(draw_df) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = draw_df,
        ggplot2::aes(
          x = .data$x, y = .data$y,
          xend = .data$xend, yend = .data$yend,
          alpha = .data$alpha_val
        ),
        color = draw_df$segment_color,
        linewidth = edge_width,
        lineend = "butt",
        show.legend = FALSE
      ) +
      ggplot2::scale_alpha_identity(guide = "none")
  } else {
    p <- p + ggplot2::scale_alpha_identity(guide = "none")
  }
  
  p <- p +
    ggplot2::geom_point(
      data = legend_df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$cell_type),
      shape = 22,
      size = 3,
      alpha = 0,
      color = NA,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = palette_map[legend_types_use],
      name = group_by,
      breaks = legend_types_use,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        override.aes = list(alpha = 1, size = 5, color = NA, shape = 22)
      )
    ) +
    ggplot2::coord_fixed(
      xlim = x_limits_use,
      ylim = y_limits_use,
      expand = FALSE
    ) +
    ggplot2::theme_classic() +
    makeDarkTheme() +
    ggplot2::ggtitle(title_use)
  
  if (isTRUE(reverse_y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  
  p
}

#' Plot multiple typed LR networks in spatial tissue coordinates for one object
plotLRInSituTypedMulti <- function(eff,
                                   lr,
                                   group_by,
                                   plot_titles = NULL,
                                   ncol = NULL,
                                   nrow = NULL,
                                   guides = "collect",
                                   show_axes = TRUE,
                                   transform = identity,
                                   # identity colors
                                   palette = NULL,
                                   # edge identity filtering
                                   cell_types = NULL,
                                   cell_types_mode = c("either", "both"),
                                   sender_types = NULL,
                                   receiver_types = NULL,
                                   # rendering
                                   show_edge_scaffold = TRUE,
                                   scaffold_color = "grey20",
                                   scaffold_alpha = 0.35,
                                   scaffold_width = 0.05,
                                   edge_gradient = FALSE,
                                   edge_width = 0.20,
                                   edge_n_segments = 9L,
                                   pt_size = 0.0,
                                   pt_color = "white",
                                   reverse_y = FALSE,
                                   # alpha controls
                                   alpha_range = c(0.0, 1),
                                   # LR combination controls
                                   lr_normalization = c("none", "max", "percentile", "zscore"),
                                   lr_aggregation   = c("sum", "mean", "max"),
                                   deduplicate_lr   = TRUE,
                                   warn_many_lr     = TRUE,
                                   # visualization filtering
                                   edge_fraction = 1,
                                   edge_top = NULL,
                                   # safety
                                   max_edges = 2e6) {
  if (missing(lr) || is.null(lr) || length(lr) == 0) {
    stop("lr must be a non-empty character vector or list.")
  }
  
  cell_types_mode  <- match.arg(cell_types_mode)
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  
  if (is.character(lr)) {
    lr_list <- as.list(lr)
    lr_labels <- lr
  } else if (is.list(lr)) {
    lr_list <- lr
    lr_labels <- vapply(
      lr_list,
      function(x) {
        if (length(x) == 1L) {
          as.character(x)
        } else {
          paste(as.character(x), collapse = ", ")
        }
      },
      character(1)
    )
  } else {
    stop("lr must be either a character vector or a list.")
  }
  
  n_panels <- length(lr_list)
  
  if (is.null(plot_titles)) {
    plot_titles <- lr_labels
  } else {
    if (length(plot_titles) != n_panels) {
      stop("plot_titles must have the same length as lr.")
    }
    plot_titles <- as.character(plot_titles)
  }
  
  # Shared palette across all panels if user did not supply one
  palette_use <- palette
  if (is.null(palette_use)) {
    all_used_types <- unique(unlist(
      lapply(lr_list, function(one_lr) {
        .collectTypedLRUsedTypes(
          eff = eff,
          lr = one_lr,
          group_by = group_by,
          transform = transform,
          cell_types = cell_types,
          cell_types_mode = cell_types_mode,
          sender_types = sender_types,
          receiver_types = receiver_types,
          lr_normalization = lr_normalization,
          lr_aggregation = lr_aggregation,
          deduplicate_lr = deduplicate_lr,
          warn_many_lr = warn_many_lr,
          edge_fraction = edge_fraction,
          edge_top = edge_top,
          max_edges = max_edges
        )
      }),
      use.names = FALSE
    ))
    
    all_used_types <- sort(unique(all_used_types))
    if (length(all_used_types) > 0L) {
      palette_use <- stats::setNames(scales::hue_pal()(length(all_used_types)), all_used_types)
    }
  }
  
  plot_list <- vector("list", n_panels)
  
  for (i in seq_len(n_panels)) {
    p <- plotLRInSituTyped(
      eff = eff,
      lr = lr_list[[i]],
      group_by = group_by,
      transform = transform,
      palette = palette_use,
      cell_types = cell_types,
      cell_types_mode = cell_types_mode,
      sender_types = sender_types,
      receiver_types = receiver_types,
      edge_gradient = edge_gradient,
      edge_width = edge_width,
      edge_n_segments = edge_n_segments,
      pt_size = pt_size,
      pt_color = pt_color,
      reverse_y = reverse_y,
      alpha_range = alpha_range,
      lr_normalization = lr_normalization,
      lr_aggregation = lr_aggregation,
      deduplicate_lr = deduplicate_lr,
      warn_many_lr = warn_many_lr,
      plot_title = plot_titles[i],
      edge_fraction = edge_fraction,
      edge_top = edge_top,
      max_edges = max_edges
    )
    
    if (!isTRUE(show_axes)) {
      p <- stripTypedLRPanelAxes(p)
    }
    
    plot_list[[i]] <- p
  }
  
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n_panels))
  }
  
  patchwork::wrap_plots(
    plot_list,
    ncol = ncol,
    nrow = nrow,
    guides = guides
  ) &
    ggplot2::theme(
      legend.position = "right",
      plot.background  = ggplot2::element_rect(fill = "black", colour = NA),
      panel.background = ggplot2::element_rect(fill = "black", colour = NA)
    )
}

#' Plot one typed LR network across multiple EffNICHES objects
plotLRInSituTypedAcrossObjects <- function(eff_list,
                                           lr,
                                           group_by,
                                           plot_titles = names(eff_list),
                                           ncol = 4,
                                           nrow = 2,
                                           guides = "collect",
                                           show_axes = FALSE,
                                           transform = identity,
                                           # identity colors
                                           palette = NULL,
                                           # edge identity filtering
                                           cell_types = NULL,
                                           cell_types_mode = c("either", "both"),
                                           sender_types = NULL,
                                           receiver_types = NULL,
                                           # contextual scaffold
                                           show_edge_scaffold = TRUE,
                                           scaffold_color = "grey20",
                                           scaffold_alpha = 0.35,
                                           scaffold_width = 0.05,
                                           # rendering
                                           edge_gradient = FALSE,
                                           edge_width = 0.20,
                                           edge_n_segments = 9L,
                                           pt_size = 0.0,
                                           pt_color = "white",
                                           reverse_y = FALSE,
                                           # alpha controls
                                           alpha_range = c(0.0, 1),
                                           # LR combination controls
                                           lr_normalization = c("none", "max", "percentile", "zscore"),
                                           lr_aggregation   = c("sum", "mean", "max"),
                                           deduplicate_lr   = TRUE,
                                           warn_many_lr     = TRUE,
                                           # visualization filtering
                                           edge_fraction = 1,
                                           edge_top = NULL,
                                           # safety
                                           max_edges = 2e6) {
  if (!is.list(eff_list) || length(eff_list) == 0L ||
      !all(vapply(eff_list, inherits, logical(1), "EffNICHES"))) {
    stop("All elements of eff_list must inherit from 'EffNICHES'.")
  }
  
  cell_types_mode  <- match.arg(cell_types_mode)
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  
  if (length(plot_titles) != length(eff_list)) {
    stop("plot_titles must match length of eff_list.")
  }
  
  shared_panel <- computeSharedSpatialPanelSize(eff_list)
  
  collectAllGroupTypes <- function(one_eff, group_by) {
    meta <- one_eff$nodes$meta
    if (is.null(meta) || !group_by %in% colnames(meta)) {
      return(character(0))
    }
    vals <- as.character(meta[[group_by]])
    vals <- vals[!is.na(vals)]
    vals <- trimws(vals)
    vals <- vals[nzchar(vals)]
    unique(vals)
  }
  
  all_legend_types <- sort(unique(unlist(
    lapply(eff_list, collectAllGroupTypes, group_by = group_by),
    use.names = FALSE
  )))
  
  if (length(all_legend_types) < 1L) {
    stop("No cell types were found across the supplied objects for group_by = '", group_by, "'.")
  }
  
  if (is.null(palette)) {
    palette_use <- stats::setNames(
      scales::hue_pal()(length(all_legend_types)),
      all_legend_types
    )
  } else {
    if (!is.character(palette) || is.null(names(palette))) {
      stop("palette must be NULL or a named character vector keyed by cell type.")
    }
    if (anyNA(names(palette)) || any(!nzchar(names(palette)))) {
      stop("palette must have non-empty names for all supplied colors.")
    }
    
    missing_cols <- setdiff(all_legend_types, names(palette))
    if (length(missing_cols) > 0L) {
      stop(
        "palette is missing colors for across-object legend cell types: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    palette_use <- palette[all_legend_types]
  }
  
  # Shared alpha limits from true edge signals across all objects
  panel_data_list <- lapply(eff_list, function(one_eff) {
    prepareLRPanelData(
      eff = one_eff,
      lr = lr,
      layout = "spatial",
      coords = NULL,
      transform = transform,
      lr_normalization = lr_normalization,
      lr_aggregation = lr_aggregation,
      deduplicate_lr = deduplicate_lr,
      warn_many_lr = warn_many_lr,
      plot_title = NULL,
      title_suffix = NULL,
      edge_fraction = 1,
      edge_top = NULL,
      max_edges = max_edges
    )
  })
  
  alpha_limits_use <- computeSharedLRAlphaLimits(panel_data_list)
  
  make_empty_panel <- function(title_text = NULL,
                               subtitle_text = "No positive edges",
                               x_limits,
                               y_limits) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_blank(ggplot2::aes(x = 0, y = 0)) +
      ggplot2::annotate(
        "text",
        x = mean(x_limits),
        y = mean(y_limits),
        label = subtitle_text,
        colour = "white",
        size = 4
      ) +
      ggplot2::coord_fixed(
        xlim = x_limits,
        ylim = y_limits,
        expand = FALSE
      ) +
      ggplot2::theme_classic() +
      makeDarkTheme() +
      ggplot2::ggtitle(title_text) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank()
      )
    
    p
  }
  
  is_empty_edge_error <- function(msg) {
    pats <- c(
      "No edges remained after cell-type filtering",
      "No positive finite edge values to plot after transform\\(\\)",
      "No positive edges remained after visualization filtering"
    )
    any(vapply(pats, function(p) grepl(p, msg), logical(1)))
  }
  
  plot_list <- vector("list", length(eff_list))
  
  for (i in seq_along(eff_list)) {
    lims_i <- computeCenteredSpatialLimits(
      eff = eff_list[[i]],
      shared_width = shared_panel$width,
      shared_height = shared_panel$height
    )
    
    plot_list[[i]] <- tryCatch(
      {
        p <- plotLRInSituTyped(
          eff = eff_list[[i]],
          lr = lr,
          group_by = group_by,
          transform = transform,
          palette = palette_use,
          legend_types = all_legend_types,
          cell_types = cell_types,
          cell_types_mode = cell_types_mode,
          sender_types = sender_types,
          receiver_types = receiver_types,
          show_edge_scaffold = show_edge_scaffold,
          scaffold_color = scaffold_color,
          scaffold_alpha = scaffold_alpha,
          scaffold_width = scaffold_width,
          edge_gradient = edge_gradient,
          edge_width = edge_width,
          edge_n_segments = edge_n_segments,
          pt_size = pt_size,
          pt_color = pt_color,
          reverse_y = reverse_y,
          alpha_range = alpha_range,
          alpha_limits = alpha_limits_use,
          x_limits = lims_i$x_limits,
          y_limits = lims_i$y_limits,
          lr_normalization = lr_normalization,
          lr_aggregation = lr_aggregation,
          deduplicate_lr = deduplicate_lr,
          warn_many_lr = warn_many_lr,
          plot_title = as.character(plot_titles[i]),
          edge_fraction = edge_fraction,
          edge_top = edge_top,
          max_edges = max_edges
        )
        
        if (!isTRUE(show_axes)) {
          p <- stripTypedLRPanelAxes(p)
        }
        
        p
      },
      error = function(e) {
        msg <- conditionMessage(e)
        
        if (is_empty_edge_error(msg)) {
          p <- make_empty_panel(
            title_text = as.character(plot_titles[i]),
            subtitle_text = "No positive edges",
            x_limits = lims_i$x_limits,
            y_limits = lims_i$y_limits
          )
          
          if (!isTRUE(show_axes)) {
            p <- stripTypedLRPanelAxes(p)
          }
          
          return(p)
        }
        
        stop(e)
      }
    )
  }
  
  patchwork::wrap_plots(
    plot_list,
    ncol = ncol,
    nrow = nrow,
    guides = guides
  ) &
    ggplot2::theme(
      legend.position = "right",
      plot.background  = ggplot2::element_rect(fill = "black", colour = NA),
      panel.background = ggplot2::element_rect(fill = "black", colour = NA)
    )
}

#### Grouped LR Mechanisms ####

# Canonical LR family groupings curated from the catalog.
# Non-canonical / weakly supported direct-binding pairs were omitted.
# Strings are standardized as "Lig-Rec" with ASCII hyphen.

### Hedgehog
lr_hh_full <- c(
  "Shh-Ptch1",
  "Dhh-Ptch1",
  "Shh-Boc",
  "Shh-Cdon",
  "Dhh-Boc",
  "Dhh-Cdon",
  "Shh-Lrp2",
  "Shh-Hhip",
  "Dhh-Hhip"
)

lr_hh_activators <- c(
  "Shh-Ptch1",
  "Dhh-Ptch1"
)

lr_hh_potentiators <- c(
  "Shh-Boc",
  "Shh-Cdon",
  "Dhh-Boc",
  "Dhh-Cdon",
  "Shh-Lrp2"
)

lr_hh_inhibitors <- c(
  "Shh-Hhip",
  "Dhh-Hhip"
)


### FGF
lr_fgf_full <- c(
  "Fgf10-Fgfr2",
  "Fgf7-Fgfr2",
  "Fgf1-Fgfr1",
  "Fgf1-Fgfr2",
  "Fgf1-Fgfr3",
  "Fgf2-Fgfr1",
  "Fgf2-Fgfr2",
  "Fgf2-Fgfr3",
  "Fgf9-Fgfr1",
  "Fgf9-Fgfr2",
  "Fgf9-Fgfr3",
  "Fgf18-Fgfr2",
  "Fgf18-Fgfr3",
  "Fgf1-Cd44",
  "Fgf2-Cd44",
  "Fgf1-Nrp1",
  "Fgf2-Nrp1",
  "Fgf7-Nrp1",
  "Fgf1-Fgfrl1",
  "Fgf2-Fgfrl1"
)

lr_fgf_activators <- c(
  "Fgf10-Fgfr2",
  "Fgf7-Fgfr2",
  "Fgf1-Fgfr1",
  "Fgf1-Fgfr2",
  "Fgf1-Fgfr3",
  "Fgf2-Fgfr1",
  "Fgf2-Fgfr2",
  "Fgf2-Fgfr3",
  "Fgf9-Fgfr1",
  "Fgf9-Fgfr2",
  "Fgf9-Fgfr3",
  "Fgf18-Fgfr2",
  "Fgf18-Fgfr3"
)

lr_fgf_potentiators <- c(
  "Fgf1-Cd44",
  "Fgf2-Cd44",
  "Fgf1-Nrp1",
  "Fgf2-Nrp1",
  "Fgf7-Nrp1"
)

lr_fgf_inhibitors <- c(
  "Fgf1-Fgfrl1",
  "Fgf2-Fgfrl1"
)


### WNT / RSPO
lr_wnt_full <- c(
  "Wnt2-Fzd1",
  "Wnt2-Fzd5",
  "Wnt2-Fzd7",
  "Wnt3a-Fzd1",
  "Wnt3a-Fzd2",
  "Wnt3a-Fzd5",
  "Wnt3a-Fzd6",
  "Wnt3a-Fzd7",
  "Wnt4-Fzd2",
  "Wnt4-Fzd6",
  "Wnt7a-Fzd5",
  "Wnt7b-Fzd1",
  "Wnt11-Fzd7",
  "Wnt3a-Ryk",
  "Wnt5a-Ryk",
  "Wnt5a-Ror1",
  "Wnt2-Lrp6",
  "Wnt3a-Lrp6",
  "Wnt5a-Lrp5",
  "Wnt7a-Lrp6",
  "Wnt7b-Lrp5",
  "Rspo1-Lgr4",
  "Rspo1-Lgr5",
  "Rspo1-Lgr6",
  "Rspo2-Lgr4",
  "Rspo2-Lgr5",
  "Rspo2-Lgr6",
  "Rspo3-Lgr4",
  "Rspo3-Lgr5",
  "Rspo3-Lgr6",
  "Cthrc1-Fzd3",
  "Cthrc1-Fzd5",
  "Cthrc1-Fzd6",
  "Dkk2-Lrp6",
  "Igfbp4-Lrp6",
  "Sfrp1-Fzd2",
  "Sfrp1-Fzd6"
)

lr_wnt_activators <- c(
  "Wnt2-Fzd1",
  "Wnt2-Fzd5",
  "Wnt2-Fzd7",
  "Wnt3a-Fzd1",
  "Wnt3a-Fzd2",
  "Wnt3a-Fzd5",
  "Wnt3a-Fzd6",
  "Wnt3a-Fzd7",
  "Wnt4-Fzd2",
  "Wnt4-Fzd6",
  "Wnt7a-Fzd5",
  "Wnt7b-Fzd1",
  "Wnt11-Fzd7",
  "Wnt3a-Ryk",
  "Wnt5a-Ryk",
  "Wnt5a-Ror1"
)

lr_wnt_potentiators <- c(
  "Wnt2-Lrp6",
  "Wnt3a-Lrp6",
  "Wnt5a-Lrp5",
  "Wnt7a-Lrp6",
  "Wnt7b-Lrp5",
  "Rspo1-Lgr4",
  "Rspo1-Lgr5",
  "Rspo1-Lgr6",
  "Rspo2-Lgr4",
  "Rspo2-Lgr5",
  "Rspo2-Lgr6",
  "Rspo3-Lgr4",
  "Rspo3-Lgr5",
  "Rspo3-Lgr6",
  "Cthrc1-Fzd3",
  "Cthrc1-Fzd5",
  "Cthrc1-Fzd6"
)

lr_wnt_inhibitors <- c(
  "Dkk2-Lrp6",
  "Igfbp4-Lrp6",
  "Sfrp1-Fzd2",
  "Sfrp1-Fzd6"
)


### BMP
lr_bmp_full <- c(
  "Bmp2-Acvr1",
  "Bmp2-Bmpr1a",
  "Bmp2-Bmpr1b",
  "Bmp2-Bmpr2",
  "Bmp2-Acvr2a",
  "Bmp2-Acvr2b",
  "Bmp4-Acvr1",
  "Bmp4-Bmpr1a",
  "Bmp4-Bmpr1b",
  "Bmp4-Bmpr2",
  "Bmp4-Acvr2a",
  "Bmp4-Acvr2b",
  "Bmp5-Acvr1",
  "Bmp5-Bmpr1a",
  "Bmp5-Bmpr1b",
  "Bmp5-Bmpr2",
  "Bmp5-Acvr2a",
  "Bmp5-Acvr2b",
  "Bmp6-Acvr1",
  "Bmp6-Bmpr1a",
  "Bmp6-Bmpr1b",
  "Bmp6-Bmpr2",
  "Bmp6-Acvr2a",
  "Bmp6-Acvr2b",
  "Bmp7-Acvr1",
  "Bmp7-Bmpr1a",
  "Bmp7-Bmpr1b",
  "Bmp7-Bmpr2",
  "Bmp7-Acvr2a",
  "Bmp7-Acvr2b",
  "Bmp3-Acvr2b",
  "Bmp3-Acvr2a"
)

lr_bmp_activators <- c(
  "Bmp2-Acvr1",
  "Bmp2-Bmpr1a",
  "Bmp2-Bmpr1b",
  "Bmp2-Bmpr2",
  "Bmp2-Acvr2a",
  "Bmp2-Acvr2b",
  "Bmp4-Acvr1",
  "Bmp4-Bmpr1a",
  "Bmp4-Bmpr1b",
  "Bmp4-Bmpr2",
  "Bmp4-Acvr2a",
  "Bmp4-Acvr2b",
  "Bmp5-Acvr1",
  "Bmp5-Bmpr1a",
  "Bmp5-Bmpr1b",
  "Bmp5-Bmpr2",
  "Bmp5-Acvr2a",
  "Bmp5-Acvr2b",
  "Bmp6-Acvr1",
  "Bmp6-Bmpr1a",
  "Bmp6-Bmpr1b",
  "Bmp6-Bmpr2",
  "Bmp6-Acvr2a",
  "Bmp6-Acvr2b",
  "Bmp7-Acvr1",
  "Bmp7-Bmpr1a",
  "Bmp7-Bmpr1b",
  "Bmp7-Bmpr2",
  "Bmp7-Acvr2a",
  "Bmp7-Acvr2b"
)

lr_bmp_potentiators <- character(0)

lr_bmp_inhibitors <- c(
  "Bmp3-Acvr2b",
  "Bmp3-Acvr2a"
)


### Activin / Inhibin
lr_activin_full <- c(
  "Inhba-Acvr2a",
  "Inhba-Acvr2b",
  "Inhba-Acvr1b",
  "Inhbb-Acvr2a",
  "Inhbb-Acvr2b",
  "Inhbb-Acvr1b",
  "Inhba-Acvr1",
  "Inha-Tgfbr3",
  "Inha-Acvr2a",
  "Inha-Acvr2b"
)

lr_activin_activators <- c(
  "Inhba-Acvr2a",
  "Inhba-Acvr2b",
  "Inhba-Acvr1b",
  "Inhbb-Acvr2a",
  "Inhbb-Acvr2b",
  "Inhbb-Acvr1b"
)

lr_activin_potentiators <- c(
  "Inha-Tgfbr3"
)

lr_activin_inhibitors <- c(
  "Inhba-Acvr1",
  "Inha-Acvr2a",
  "Inha-Acvr2b"
)


### TGF-beta
lr_tgfb_full <- c(
  "Tgfb1-Tgfbr2",
  "Tgfb1-Tgfbr1",
  "Tgfb2-Tgfbr2",
  "Tgfb2-Tgfbr1",
  "Tgfb3-Tgfbr2",
  "Tgfb3-Tgfbr1",
  "Tgfb1-Tgfbr3",
  "Tgfb2-Tgfbr3",
  "Tgfb3-Tgfbr3",
  "Tgfb1-Eng",
  "Tgfb2-Eng",
  "Tgfb3-Eng",
  "Tgfb1-Itgb6",
  "Tgfb1-Cd109"
)

lr_tgfb_activators <- c(
  "Tgfb1-Tgfbr2",
  "Tgfb1-Tgfbr1",
  "Tgfb2-Tgfbr2",
  "Tgfb2-Tgfbr1",
  "Tgfb3-Tgfbr2",
  "Tgfb3-Tgfbr1"
)

lr_tgfb_potentiators <- c(
  "Tgfb1-Tgfbr3",
  "Tgfb2-Tgfbr3",
  "Tgfb3-Tgfbr3",
  "Tgfb1-Eng",
  "Tgfb2-Eng",
  "Tgfb3-Eng",
  "Tgfb1-Itgb6"
)

lr_tgfb_inhibitors <- c(
  "Tgfb1-Cd109"
)


### RGM-BMP modulators
lr_rgm_full <- c(
  "Rgma-Neo1",
  "Rgmb-Neo1",
  "Rgma-Bmpr1b",
  "Rgma-Bmpr2",
  "Rgmb-Bmpr1b",
  "Rgmb-Bmpr2"
)

lr_rgm_activators <- c(
  "Rgma-Neo1",
  "Rgmb-Neo1"
)

lr_rgm_potentiators <- c(
  "Rgma-Bmpr1b",
  "Rgma-Bmpr2",
  "Rgmb-Bmpr1b",
  "Rgmb-Bmpr2"
)

lr_rgm_inhibitors <- character(0)


### Notch
lr_notch_full <- c(
  "Dll4-Notch1",
  "Dll4-Notch2",
  "Dll4-Notch3",
  "Dll4-Notch4",
  "Jag1-Notch1",
  "Jag1-Notch2",
  "Jag1-Notch3",
  "Jag1-Notch4",
  "Jag2-Notch1",
  "Jag2-Notch2",
  "Jag2-Notch3",
  "Jag2-Notch4"
)

lr_notch_activators <- lr_notch_full
lr_notch_potentiators <- character(0)
lr_notch_inhibitors <- character(0)


### EGFR / ERBB
lr_egfr_full <- c(
  "Areg-Egfr",
  "Egf-Egfr",
  "Tgfa-Egfr",
  "Ereg-Egfr",
  "Hbegf-Egfr",
  "Nrg2-Erbb3",
  "Hbegf-Cd9",
  "Hbegf-Cd44"
)

lr_egfr_activators <- c(
  "Areg-Egfr",
  "Egf-Egfr",
  "Tgfa-Egfr",
  "Ereg-Egfr",
  "Hbegf-Egfr",
  "Nrg2-Erbb3"
)

lr_egfr_potentiators <- c(
  "Hbegf-Cd9",
  "Hbegf-Cd44"
)

lr_egfr_inhibitors <- character(0)


### Semaphorin
lr_sema_full <- c(
  "Sema3a-Nrp1",
  "Sema3a-Nrp2",
  "Sema3b-Nrp1",
  "Sema3b-Nrp2",
  "Sema3c-Nrp1",
  "Sema3c-Nrp2",
  "Sema3f-Nrp1",
  "Sema3f-Nrp2",
  "Sema3a-Plxna1",
  "Sema3a-Plxna3",
  "Sema3f-Plxna1",
  "Sema3f-Plxna3",
  "Sema3e-Nrp1",
  "Sema4d-Met",
  "Sema6d-Plxna1",
  "Sema6d-Kdr"
)

lr_sema_activators <- c(
  "Sema3a-Nrp1",
  "Sema3a-Nrp2",
  "Sema3b-Nrp1",
  "Sema3b-Nrp2",
  "Sema3c-Nrp1",
  "Sema3c-Nrp2",
  "Sema3f-Nrp1",
  "Sema3f-Nrp2",
  "Sema3a-Plxna1",
  "Sema3a-Plxna3",
  "Sema3f-Plxna1",
  "Sema3f-Plxna3",
  "Sema6d-Plxna1"
)

lr_sema_potentiators <- c(
  "Sema3e-Nrp1",
  "Sema4d-Met",
  "Sema6d-Kdr"
)

lr_sema_inhibitors <- character(0)


### Netrin
lr_netrin_full <- c(
  "Ntn1-Neo1",
  "Ntn3-Neo1",
  "Ntn3-Cdon"
)

lr_netrin_activators <- c(
  "Ntn1-Neo1",
  "Ntn3-Neo1"
)

lr_netrin_potentiators <- c(
  "Ntn3-Cdon"
)

lr_netrin_inhibitors <- character(0)


### Slit-Robo
lr_slit_full <- c(
  "Slit2-Robo2",
  "Slit3-Robo2",
  "Slit2-Gpc1"
)

lr_slit_activators <- c(
  "Slit2-Robo2",
  "Slit3-Robo2"
)

lr_slit_potentiators <- c(
  "Slit2-Gpc1"
)

lr_slit_inhibitors <- character(0)


### Ephrin-Eph
lr_eph_full <- c(
  "Efna1-Epha1",
  "Efna1-Epha2",
  "Efna1-Epha4",
  "Efna1-Epha7",
  "Efna2-Epha1",
  "Efna2-Epha2",
  "Efna2-Epha4",
  "Efna2-Epha7",
  "Efna4-Epha1",
  "Efna4-Epha2",
  "Efna4-Epha4",
  "Efna4-Epha7",
  "Efna5-Epha1",
  "Efna5-Epha2",
  "Efna5-Epha4",
  "Efna5-Epha7",
  "Efna5-Ephb2",
  "Efnb1-Epha4",
  "Efnb1-Ephb1",
  "Efnb1-Ephb2",
  "Efnb1-Ephb3",
  "Efnb1-Ephb4",
  "Efnb1-Ephb6",
  "Efnb2-Epha4",
  "Efnb2-Ephb1",
  "Efnb2-Ephb2",
  "Efnb2-Ephb3",
  "Efnb2-Ephb4",
  "Efnb2-Ephb6",
  "Efnb3-Epha4",
  "Efnb3-Ephb1",
  "Efnb3-Ephb2",
  "Efnb3-Ephb3",
  "Efnb3-Ephb4",
  "Efnb3-Ephb6"
)

lr_eph_activators <- lr_eph_full
lr_eph_potentiators <- character(0)
lr_eph_inhibitors <- character(0)


### Neurotrophin / CNTF
lr_neurotrophin_full <- c(
  "Bdnf-Ntrk2",
  "Ntf3-Ntrk3",
  "Ntf3-Ntrk2",
  "Ntf4-Ntrk2",
  "Ntf4-Ntrk3",
  "Cntf-Cntfr",
  "Cntf-Lifr",
  "Cntf-Il6st"
)

lr_neurotrophin_activators <- lr_neurotrophin_full
lr_neurotrophin_potentiators <- character(0)
lr_neurotrophin_inhibitors <- character(0)


### Synaptic adhesion
lr_synaptic_full <- c(
  "Agrn-Lrp4",
  "Nlgn2-Nrxn1",
  "Nlgn2-Nrxn2",
  "Nlgn3-Nrxn1",
  "Nlgn3-Nrxn2",
  "Nxph1-Nrxn1",
  "Nxph1-Nrxn2"
)

lr_synaptic_activators <- lr_synaptic_full
lr_synaptic_potentiators <- character(0)
lr_synaptic_inhibitors <- character(0)


### Angiopoietin / Apelin
lr_angio_tie_full <- c(
  "Angpt1-Tek",
  "Angpt2-Tek",
  "Angpt4-Tek",
  "Angpt1-Itga5",
  "Apln-Aplnr"
)

lr_angio_tie_activators <- c(
  "Angpt1-Tek",
  "Angpt2-Tek",
  "Angpt4-Tek",
  "Apln-Aplnr"
)

lr_angio_tie_potentiators <- c(
  "Angpt1-Itga5"
)

lr_angio_tie_inhibitors <- character(0)


### PDGF
lr_pdgf_full <- c(
  "Pdgfa-Pdgfra",
  "Pdgfb-Pdgfra",
  "Pdgfb-Pdgfrb",
  "Pdgfc-Pdgfra",
  "Pdgfc-Pdgfrb",
  "Pdgfd-Pdgfrb",
  "Pdgfc-Kdr"
)

lr_pdgf_activators <- c(
  "Pdgfa-Pdgfra",
  "Pdgfb-Pdgfra",
  "Pdgfb-Pdgfrb",
  "Pdgfc-Pdgfra",
  "Pdgfc-Pdgfrb",
  "Pdgfd-Pdgfrb"
)

lr_pdgf_potentiators <- c(
  "Pdgfc-Kdr"
)

lr_pdgf_inhibitors <- character(0)


### VEGF
lr_vegf_full <- c(
  "Vegfa-Flt1",
  "Vegfa-Kdr",
  "Vegfa-Nrp1",
  "Vegfa-Nrp2",
  "Vegfa-Gpc1",
  "Vegfc-Flt4",
  "Vegfc-Kdr",
  "Vegfc-Nrp2",
  "Pgf-Flt1",
  "Pgf-Nrp1",
  "Pgf-Nrp2",
  "Col18a1-Kdr"
)

lr_vegf_activators <- c(
  "Vegfa-Flt1",
  "Vegfa-Kdr",
  "Vegfc-Flt4",
  "Vegfc-Kdr",
  "Pgf-Flt1"
)

lr_vegf_potentiators <- c(
  "Vegfa-Nrp1",
  "Vegfa-Nrp2",
  "Vegfa-Gpc1",
  "Vegfc-Nrp2",
  "Pgf-Nrp1",
  "Pgf-Nrp2"
)

lr_vegf_inhibitors <- c(
  "Col18a1-Kdr"
)


### HGF
lr_hgf_full <- c(
  "Hgf-Met",
  "Hgf-Cd44",
  "Hgf-Sdc2"
)

lr_hgf_activators <- c(
  "Hgf-Met"
)

lr_hgf_potentiators <- c(
  "Hgf-Cd44",
  "Hgf-Sdc2"
)

lr_hgf_inhibitors <- character(0)


### IGF
lr_igf_full <- c(
  "Igf1-Igf1r",
  "Igf1-Insr",
  "Igf2-Igf1r",
  "Igf2-Insr"
)

lr_igf_activators <- lr_igf_full
lr_igf_potentiators <- character(0)
lr_igf_inhibitors <- character(0)


### KIT
lr_kit_full <- c(
  "Kitlg-Kit"
)

lr_kit_activators <- c(
  "Kitlg-Kit"
)

lr_kit_potentiators <- character(0)
lr_kit_inhibitors <- character(0)


### Secretoglobin / Megalin
lr_scgb_full <- c(
  "Scgb1a1-Lrp2"
)

lr_scgb_activators <- character(0)
lr_scgb_potentiators <- character(0)
lr_scgb_inhibitors <- c(
  "Scgb1a1-Lrp2"
)


### Secretin
lr_secretin_full <- c(
  "Sct-Sctr"
)

lr_secretin_activators <- c(
  "Sct-Sctr"
)

lr_secretin_potentiators <- character(0)
lr_secretin_inhibitors <- character(0)


### SP-D / innate regulation
lr_spd_full <- c(
  "Sftpd-Tlr4"
)

lr_spd_activators <- character(0)
lr_spd_potentiators <- character(0)
lr_spd_inhibitors <- c(
  "Sftpd-Tlr4"
)


### Chemokines
lr_chemokine_full <- c(
  "Cxcl11-Cxcr3",
  "Cxcl12-Cxcr4",
  "Cxcl16-Cxcr6",
  "Cxcl1-Cxcr2",
  "Cxcl2-Cxcr2",
  "Cxcl3-Cxcr2",
  "Cxcl6-Cxcr2",
  "Ccl3-Ccr1",
  "Ccl3-Ccr5",
  "Ccl11-Ccr3",
  "Ccl11-Ccr5",
  "Ccl7-Ccr1",
  "Ccl7-Ccr3",
  "Ccl7-Ccr5",
  "Cxcl11-Ackr3",
  "Cxcl12-Ackr3"
)

lr_chemokine_activators <- c(
  "Cxcl11-Cxcr3",
  "Cxcl12-Cxcr4",
  "Cxcl16-Cxcr6",
  "Cxcl1-Cxcr2",
  "Cxcl2-Cxcr2",
  "Cxcl3-Cxcr2",
  "Cxcl6-Cxcr2",
  "Ccl3-Ccr1",
  "Ccl3-Ccr5",
  "Ccl11-Ccr3",
  "Ccl11-Ccr5",
  "Ccl7-Ccr1",
  "Ccl7-Ccr3",
  "Ccl7-Ccr5"
)

lr_chemokine_potentiators <- character(0)

lr_chemokine_inhibitors <- c(
  "Cxcl11-Ackr3",
  "Cxcl12-Ackr3"
)


### Cytokines / IL-6 family / CSF
lr_cytokine_full <- c(
  "Csf1-Csf1r",
  "Il34-Csf1r",
  "Il6-Il6r",
  "Il6-Il6st",
  "Il11-Il11ra1",
  "Il11-Il6st",
  "Il2-Il2ra",
  "Il2-Il2rb",
  "Il2-Il2rg",
  "Il7-Il7r",
  "Il7-Il2rg",
  "Il17b-Il17rb",
  "Il25-Il17rb",
  "Il1rn-Il1r1",
  "Il1rn-Il1r2"
)

lr_cytokine_activators <- c(
  "Csf1-Csf1r",
  "Il34-Csf1r",
  "Il6-Il6r",
  "Il6-Il6st",
  "Il11-Il11ra1",
  "Il11-Il6st",
  "Il2-Il2ra",
  "Il2-Il2rb",
  "Il2-Il2rg",
  "Il7-Il7r",
  "Il7-Il2rg",
  "Il17b-Il17rb",
  "Il25-Il17rb"
)

lr_cytokine_potentiators <- character(0)

lr_cytokine_inhibitors <- c(
  "Il1rn-Il1r1",
  "Il1rn-Il1r2"
)


### TNF superfamily
lr_tnf_full <- c(
  "Faslg-Fas",
  "Tnf-Tnfrsf1a",
  "Tnfsf10-Tnfrsf11b",
  "Tnfsf15-Tnfrsf25"
)

lr_tnf_activators <- c(
  "Faslg-Fas",
  "Tnf-Tnfrsf1a",
  "Tnfsf15-Tnfrsf25"
)

lr_tnf_potentiators <- character(0)

lr_tnf_inhibitors <- c(
  "Tnfsf10-Tnfrsf11b"
)


### ECM / fibrosis / mechanoregulation
lr_ecm_full <- c(
  "Bgn-Tlr2",
  "Bgn-Tlr4",
  "Col1a1-Ddr1",
  "Tnc-Egfr",
  "Thbs1-Cd47",
  "Col4a3-Cd47",
  "Mmp9-Cd44",
  "Dcn-Met",
  "Dcn-Egfr",
  "Thbs2-Notch3",
  "Vim-Cd44"
)

lr_ecm_activators <- c(
  "Bgn-Tlr2",
  "Bgn-Tlr4",
  "Col1a1-Ddr1",
  "Tnc-Egfr",
  "Thbs1-Cd47",
  "Col4a3-Cd47"
)

lr_ecm_potentiators <- c(
  "Mmp9-Cd44",
  "Thbs2-Notch3",
  "Vim-Cd44"
)

lr_ecm_inhibitors <- c(
  "Dcn-Met",
  "Dcn-Egfr"
)


### Master list
lr_family_groups <- list(
  hh = list(full = lr_hh_full, activators = lr_hh_activators, potentiators = lr_hh_potentiators, inhibitors = lr_hh_inhibitors),
  fgf = list(full = lr_fgf_full, activators = lr_fgf_activators, potentiators = lr_fgf_potentiators, inhibitors = lr_fgf_inhibitors),
  wnt = list(full = lr_wnt_full, activators = lr_wnt_activators, potentiators = lr_wnt_potentiators, inhibitors = lr_wnt_inhibitors),
  bmp = list(full = lr_bmp_full, activators = lr_bmp_activators, potentiators = lr_bmp_potentiators, inhibitors = lr_bmp_inhibitors),
  activin = list(full = lr_activin_full, activators = lr_activin_activators, potentiators = lr_activin_potentiators, inhibitors = lr_activin_inhibitors),
  tgfb = list(full = lr_tgfb_full, activators = lr_tgfb_activators, potentiators = lr_tgfb_potentiators, inhibitors = lr_tgfb_inhibitors),
  rgm = list(full = lr_rgm_full, activators = lr_rgm_activators, potentiators = lr_rgm_potentiators, inhibitors = lr_rgm_inhibitors),
  notch = list(full = lr_notch_full, activators = lr_notch_activators, potentiators = lr_notch_potentiators, inhibitors = lr_notch_inhibitors),
  egfr = list(full = lr_egfr_full, activators = lr_egfr_activators, potentiators = lr_egfr_potentiators, inhibitors = lr_egfr_inhibitors),
  sema = list(full = lr_sema_full, activators = lr_sema_activators, potentiators = lr_sema_potentiators, inhibitors = lr_sema_inhibitors),
  netrin = list(full = lr_netrin_full, activators = lr_netrin_activators, potentiators = lr_netrin_potentiators, inhibitors = lr_netrin_inhibitors),
  slit = list(full = lr_slit_full, activators = lr_slit_activators, potentiators = lr_slit_potentiators, inhibitors = lr_slit_inhibitors),
  eph = list(full = lr_eph_full, activators = lr_eph_activators, potentiators = lr_eph_potentiators, inhibitors = lr_eph_inhibitors),
  neurotrophin = list(full = lr_neurotrophin_full, activators = lr_neurotrophin_activators, potentiators = lr_neurotrophin_potentiators, inhibitors = lr_neurotrophin_inhibitors),
  synaptic = list(full = lr_synaptic_full, activators = lr_synaptic_activators, potentiators = lr_synaptic_potentiators, inhibitors = lr_synaptic_inhibitors),
  angio_tie = list(full = lr_angio_tie_full, activators = lr_angio_tie_activators, potentiators = lr_angio_tie_potentiators, inhibitors = lr_angio_tie_inhibitors),
  pdgf = list(full = lr_pdgf_full, activators = lr_pdgf_activators, potentiators = lr_pdgf_potentiators, inhibitors = lr_pdgf_inhibitors),
  vegf = list(full = lr_vegf_full, activators = lr_vegf_activators, potentiators = lr_vegf_potentiators, inhibitors = lr_vegf_inhibitors),
  hgf = list(full = lr_hgf_full, activators = lr_hgf_activators, potentiators = lr_hgf_potentiators, inhibitors = lr_hgf_inhibitors),
  igf = list(full = lr_igf_full, activators = lr_igf_activators, potentiators = lr_igf_potentiators, inhibitors = lr_igf_inhibitors),
  kit = list(full = lr_kit_full, activators = lr_kit_activators, potentiators = lr_kit_potentiators, inhibitors = lr_kit_inhibitors),
  scgb = list(full = lr_scgb_full, activators = lr_scgb_activators, potentiators = lr_scgb_potentiators, inhibitors = lr_scgb_inhibitors),
  secretin = list(full = lr_secretin_full, activators = lr_secretin_activators, potentiators = lr_secretin_potentiators, inhibitors = lr_secretin_inhibitors),
  sftpd = list(full = lr_spd_full, activators = lr_spd_activators, potentiators = lr_spd_potentiators, inhibitors = lr_spd_inhibitors),
  chemokine = list(full = lr_chemokine_full, activators = lr_chemokine_activators, potentiators = lr_chemokine_potentiators, inhibitors = lr_chemokine_inhibitors),
  cytokine = list(full = lr_cytokine_full, activators = lr_cytokine_activators, potentiators = lr_cytokine_potentiators, inhibitors = lr_cytokine_inhibitors),
  tnf = list(full = lr_tnf_full, activators = lr_tnf_activators, potentiators = lr_tnf_potentiators, inhibitors = lr_tnf_inhibitors),
  ecm = list(full = lr_ecm_full, activators = lr_ecm_activators, potentiators = lr_ecm_potentiators, inhibitors = lr_ecm_inhibitors)
)

