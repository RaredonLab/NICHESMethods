# PlotDev
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
  s <- gsub("\u2014", "-", s, fixed = TRUE)
  s <- gsub("\u2013", "-", s, fixed = TRUE)
  s <- gsub("\u2212", "-", s, fixed = TRUE)
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
  
  normalize_one <- switch(
    method,
    max = function(x) {
      x[!is.finite(x)] <- NA_real_
      if (all(is.na(x))) return(rep(0, length(x)))
      mx <- suppressWarnings(max(x, na.rm = TRUE))
      if (!is.finite(mx) || mx <= 0) {
        rep(0, length(x))
      } else {
        x / mx
      }
    },
    percentile = function(x) {
      x[!is.finite(x)] <- NA_real_
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
      y
    },
    zscore = function(x) {
      x[!is.finite(x)] <- NA_real_
      if (all(is.na(x))) return(rep(0, length(x)))
      mu <- suppressWarnings(mean(x, na.rm = TRUE))
      s  <- suppressWarnings(stats::sd(x, na.rm = TRUE))
      if (!is.finite(s) || s == 0) {
        rep(0, length(x))
      } else {
        (x - mu) / s
      }
    }
  )
  
  for (j in seq_len(ncol(mat))) {
    out[, j] <- normalize_one(mat[, j])
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
  
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::rowMaxs(mat, na.rm = TRUE))
  }
  
  apply(mat, 1, max, na.rm = TRUE)
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
    stop("Column mismatch: ncol(CellToCellSpatial$w) must equal length(CellToCellSpatial$mechanisms).")
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
  
  w_edge <- transform(w_edge_raw)
  w_edge[!is.finite(w_edge)] <- NA_real_
  
  list(
    edge_signal = w_edge,
    edge_signal_raw = w_edge_raw,
    selected_idx = idx,
    selected_mechanisms = selected_mechanisms,
    normalization = lr_normalization,
    aggregation = lr_aggregation
  )
}

#' Build a plot title from selected LR mechanisms
makeLRPlotTitle <- function(selected_mechanisms,
                            plot_title = NULL,
                            suffix = NULL) {
  if (!is.null(plot_title)) {
    moi <- plot_title
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

#' Convert node coordinates to a node plotting data.frame
buildNodePlotData <- function(xy) {
  if (is.null(xy) || !is.matrix(xy) || ncol(xy) != 2) {
    stop("xy must be an N x 2 matrix.")
  }
  
  if (is.null(colnames(xy)) || !all(c("x", "y") %in% colnames(xy))) {
    colnames(xy) <- c("x", "y")
  }
  
  df_nodes <- data.frame(
    x = as.numeric(xy[, "x"]),
    y = as.numeric(xy[, "y"]),
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
  ij_num <- matrix(as.integer(ij), ncol = 2)
  
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
  
  x1 <- as.numeric(xy[from_idx, "x"])
  y1 <- as.numeric(xy[from_idx, "y"])
  x2 <- as.numeric(xy[to_idx,   "x"])
  y2 <- as.numeric(xy[to_idx,   "y"])
  
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

#' Compute alpha values for edge plotting
computeEdgeAlpha <- function(values,
                             alpha_range = c(0.05, 1)) {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }
  
  a_rng <- alpha_range
  if (length(a_rng) != 2 || any(!is.finite(a_rng)) ||
      a_rng[1] < 0 || a_rng[2] <= 0 || a_rng[1] >= a_rng[2]) {
    stop("alpha_range must be something like c(0.05, 1).")
  }
  
  observed_range <- range(values, na.rm = TRUE)
  
  if (diff(observed_range) == 0) {
    alpha_val <- rep(a_rng[2], length(values))
  } else {
    alpha_val <- scales::rescale(
      values,
      to = a_rng,
      from = observed_range
    )
    alpha_val[!is.finite(alpha_val)] <- a_rng[1]
  }
  
  alpha_val
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
                                color_option = "C",
                                color_limits = NULL,
                                oob_squish = TRUE) {
  if (!is.data.frame(edge_df) || !all(c("x1", "y1", "x2", "y2", "for.plotting") %in% names(edge_df))) {
    stop("edge_df must contain x1, y1, x2, y2, and for.plotting.")
  }
  
  lims <- validateColorLimits(edge_df$for.plotting, color_limits = color_limits)
  edge_df$alpha_val <- computeEdgeAlpha(edge_df$for.plotting, alpha_range = alpha_range)
  edge_df <- edge_df[order(edge_df$for.plotting), , drop = FALSE]
  
  color_oob_fun <- if (isTRUE(oob_squish)) scales::squish else scales::censor
  
  p <- ggplot2::ggplot()
  
  # Nodes first so they stay behind edges
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
    ggplot2::coord_fixed() +
    ggplot2::theme_classic() +
    makeDarkTheme() +
    ggplot2::ggtitle(plot_title)
  
  if (isTRUE(reverse_y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  
  p
}

.plotLRMulti <- function(plot_fun,
                         eff,
                         lr,
                         plot_titles = NULL,
                         ncol = NULL,
                         nrow = NULL,
                         guides = "collect",
                         ...) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for multi-plot wrappers.")
  }
  
  if (missing(lr) || is.null(lr) || length(lr) == 0) {
    stop("lr must be a non-empty character vector or list.")
  }
  
  # Standardize LR input
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
  
  n_plots <- length(lr_list)
  
  # Standardize plot titles
  if (is.null(plot_titles)) {
    plot_titles <- lr_labels
  } else {
    if (length(plot_titles) != n_plots) {
      stop("plot_titles must have the same length as lr.")
    }
    plot_titles <- as.character(plot_titles)
  }
  
  # Build plots
  plot_list <- vector("list", n_plots)
  for (i in seq_len(n_plots)) {
    plot_list[[i]] <- plot_fun(
      eff = eff,
      lr = lr_list[[i]],
      plot_title = plot_titles[[i]],
      ...
    )
  }
  
  # Sensible default layout if neither ncol nor nrow is provided
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n_plots))
  }
  
  patchwork::wrap_plots(plot_list, ncol = ncol, nrow = nrow, guides = guides) &
    ggplot2::theme(legend.position = "right")
}

#' Generic helper to attach stacked-lane metadata to curve geometry
#'
#' This standardizes the contract used by family and compare stacked-lane
#' builders. Mode-specific wrappers should validate their own required edge
#' columns, then delegate here.
attachStackedCurveMetadata <- function(curve_df,
                                       edge_row,
                                       id_col,
                                       id_value,
                                       extra_cols = character(0),
                                       default_self_arrow_side = -1L,
                                       default_nonself_arrow_side = 1L) {
  req_curve <- c("x", "y", "point_id", "from_group", "to_group")
  
  if (!is.data.frame(curve_df) || !all(req_curve %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req_curve, collapse = ", "))
  }
  if (!is.data.frame(edge_row) || nrow(edge_row) != 1L) {
    stop("edge_row must be a one-row data.frame.")
  }
  if (!is.character(id_col) || length(id_col) != 1L || !nzchar(id_col)) {
    stop("id_col must be a single non-empty character string.")
  }
  if (!is.numeric(id_value) || length(id_value) != 1L || !is.finite(id_value)) {
    stop("id_value must be a single finite numeric value.")
  }
  if (!is.character(extra_cols)) {
    stop("extra_cols must be a character vector.")
  }
  
  miss_extra <- setdiff(extra_cols, names(edge_row))
  if (length(miss_extra) > 0L) {
    stop(
      "edge_row is missing required metadata column(s): ",
      paste(miss_extra, collapse = ", "), "."
    )
  }
  
  out <- curve_df
  
  is_self_i <- identical(
    as.character(edge_row$from_group),
    as.character(edge_row$to_group)
  )
  
  if (!"arrow_side" %in% names(out)) {
    out$arrow_side <- if (is_self_i) default_self_arrow_side else default_nonself_arrow_side
  }
  if (!"is_self" %in% names(out)) {
    out$is_self <- is_self_i
  }
  
  for (nm in extra_cols) {
    out[[nm]] <- edge_row[[nm]]
  }
  
  out[[id_col]] <- as.integer(id_value)
  rownames(out) <- NULL
  out
}

#' Create one base directional lane per retained from_group -> to_group edge
#'
#' Generic alias for the shared stacked-circuit base-lane geometry builder.
buildStackedCircuitBaseCurves <- function(node_tbl,
                                          edge_tbl,
                                          edge_gap = 0.015,
                                          bidirectional_offset = NULL,
                                          self_loop_spread = 0.5,
                                          self_loop_height = 0.5,
                                          curve_strength = 0.16,
                                          n_pts = 401,
                                          center = c(0, 0)) {
  req_node <- c("group", "x", "y", "theta", "radius_plot")
  req_edge <- c("from_group", "to_group")
  
  if (!is.data.frame(node_tbl) || !all(req_node %in% names(node_tbl))) {
    stop("node_tbl must contain: ", paste(req_node, collapse = ", "))
  }
  if (!is.data.frame(edge_tbl) || !all(req_edge %in% names(edge_tbl))) {
    stop("edge_tbl must contain: ", paste(req_edge, collapse = ", "))
  }
  
  if (nrow(edge_tbl) == 0L) {
    return(data.frame())
  }
  
  base_edges <- unique(edge_tbl[, c("from_group", "to_group"), drop = FALSE])
  base_edges$weight_raw <- 1
  base_edges$weight_plot <- 1
  
  curve_df <- buildCircuitCurves(
    node_tbl = node_tbl,
    edge_tbl = base_edges,
    edge_gap = edge_gap,
    bidirectional_offset = bidirectional_offset,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    curve_strength = curve_strength,
    n_pts = n_pts,
    center = center
  )
  
  if (nrow(curve_df) == 0L) {
    return(curve_df)
  }
  
  curve_df$base_edge_id <- interaction(curve_df$from_group, curve_df$to_group, drop = TRUE, lex.order = TRUE)
  curve_df
}

#' Offset a base lane by a constant distance in data units
#'
#' Generic alias for the shared stacked-circuit offset helper.
offsetStackedCircuitCurve <- function(curve_df,
                                      offset_distance,
                                      center = c(0, 0),
                                      probe_offset = NULL) {
  req <- c("x", "y", "point_id")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  if (!is.numeric(offset_distance) || length(offset_distance) != 1L ||
      !is.finite(offset_distance) || offset_distance < 0) {
    stop("offset_distance must be a single non-negative number.")
  }
  
  d <- curve_df[order(curve_df$point_id), , drop = FALSE]
  if (nrow(d) < 2L || offset_distance == 0) {
    return(d)
  }
  
  xy <- as.matrix(d[, c("x", "y"), drop = FALSE])
  storage.mode(xy) <- "double"
  
  tang <- .circuit_estimate_tangent(xy)
  n_left <- cbind(x = -tang[, "y"], y = tang[, "x"])
  
  if ("arrow_side" %in% names(d)) {
    side <- unique(d$arrow_side)
    side <- side[is.finite(side)]
    if (length(side) == 1L && side != 0) {
      side_use <- sign(side)
    } else {
      side_use <- 1
    }
    
    if ("is_self" %in% names(d)) {
      is_self_i <- unique(d$is_self)
      is_self_i <- is_self_i[!is.na(is_self_i)]
      if (length(is_self_i) == 1L && isTRUE(is_self_i)) {
        side_use <- -side_use
      }
    }
  } else {
    if (is.null(probe_offset)) {
      probe_offset <- max(offset_distance * 0.25, 1e-4)
    }
    
    side_use <- inferStackedCircuitOutwardSide(
      curve_df = d,
      center = center,
      probe_offset = probe_offset
    )
  }
  
  d$x <- xy[, 1] + side_use * offset_distance * n_left[, "x"]
  d$y <- xy[, 2] + side_use * offset_distance * n_left[, "y"]
  d
}

#' Build one self-loop track for stacked-circuit rendering
#'
#' Generic alias for the shared stacked-circuit self-loop geometry helper.
buildStackedCircuitSelfLoopCurve <- function(node_row,
                                             offset_index,
                                             n_tracks_on_edge = 1L,
                                             edge_gap = 0.015,
                                             self_loop_spread = 0.5,
                                             self_loop_height = 0.5,
                                             track_spacing = 0.0035,
                                             n_pts = 401,
                                             base_half_angle = 0.30 * pi,
                                             angle_step_scale = 1.15,
                                             max_half_angle = 0.80 * pi,
                                             loop_compact_threshold = 7L,
                                             loop_compact_step = 0.022 * pi,
                                             min_base_half_angle = 0.12 * pi) {
  req_node <- c("group", "x", "y", "theta", "radius_plot")
  
  if (!is.data.frame(node_row) || nrow(node_row) != 1L || !all(req_node %in% names(node_row))) {
    stop("node_row must be a one-row data.frame containing: ", paste(req_node, collapse = ", "))
  }
  if (!is.numeric(offset_index) || length(offset_index) != 1L ||
      !is.finite(offset_index) || offset_index < 0) {
    stop("offset_index must be a single non-negative number.")
  }
  if (!is.numeric(n_tracks_on_edge) || length(n_tracks_on_edge) != 1L ||
      !is.finite(n_tracks_on_edge) || n_tracks_on_edge < 1) {
    stop("n_tracks_on_edge must be a single positive number.")
  }
  if (!is.numeric(edge_gap) || length(edge_gap) != 1L || !is.finite(edge_gap) || edge_gap < 0) {
    stop("edge_gap must be a single non-negative number.")
  }
  if (!is.numeric(self_loop_spread) || length(self_loop_spread) != 1L ||
      !is.finite(self_loop_spread) || self_loop_spread <= 0) {
    stop("self_loop_spread must be a single positive number.")
  }
  if (!is.numeric(self_loop_height) || length(self_loop_height) != 1L ||
      !is.finite(self_loop_height) || self_loop_height <= 0) {
    stop("self_loop_height must be a single positive number.")
  }
  if (!is.numeric(track_spacing) || length(track_spacing) != 1L ||
      !is.finite(track_spacing) || track_spacing < 0) {
    stop("track_spacing must be a single non-negative number.")
  }
  if (!is.numeric(n_pts) || length(n_pts) != 1L || !is.finite(n_pts) || n_pts < 50) {
    stop("n_pts must be a single integer >= 50.")
  }
  if (!is.numeric(base_half_angle) || length(base_half_angle) != 1L ||
      !is.finite(base_half_angle) || base_half_angle <= 0 || base_half_angle >= pi) {
    stop("base_half_angle must be a single number in (0, pi).")
  }
  if (!is.numeric(angle_step_scale) || length(angle_step_scale) != 1L ||
      !is.finite(angle_step_scale) || angle_step_scale <= 0) {
    stop("angle_step_scale must be a single positive number.")
  }
  if (!is.numeric(max_half_angle) || length(max_half_angle) != 1L ||
      !is.finite(max_half_angle) || max_half_angle <= 0 || max_half_angle >= pi) {
    stop("max_half_angle must be a single number in (0, pi).")
  }
  if (!is.numeric(loop_compact_threshold) || length(loop_compact_threshold) != 1L ||
      !is.finite(loop_compact_threshold) || loop_compact_threshold < 1) {
    stop("loop_compact_threshold must be a single positive number.")
  }
  if (!is.numeric(loop_compact_step) || length(loop_compact_step) != 1L ||
      !is.finite(loop_compact_step) || loop_compact_step < 0) {
    stop("loop_compact_step must be a single non-negative number.")
  }
  if (!is.numeric(min_base_half_angle) || length(min_base_half_angle) != 1L ||
      !is.finite(min_base_half_angle) || min_base_half_angle <= 0 || min_base_half_angle >= pi) {
    stop("min_base_half_angle must be a single number in (0, pi).")
  }
  
  n_pts <- as.integer(n_pts)
  n_tracks_on_edge <- as.integer(round(n_tracks_on_edge))
  loop_compact_threshold <- as.integer(round(loop_compact_threshold))
  
  .cubic_bezier <- function(p0, p1, p2, p3, t) {
    omt <- 1 - t
    cbind(
      x = omt^3 * p0[1] +
        3 * omt^2 * t * p1[1] +
        3 * omt * t^2 * p2[1] +
        t^3 * p3[1],
      y = omt^3 * p0[2] +
        3 * omt^2 * t * p1[2] +
        3 * omt * t^2 * p2[2] +
        t^3 * p3[2]
    )
  }
  
  s <- node_row
  
  u_out <- c(cos(s$theta), sin(s$theta))
  t_hat <- c(-u_out[2], u_out[1])
  
  r_contact <- s$radius_plot + edge_gap
  
  n_extra <- max(0L, n_tracks_on_edge - loop_compact_threshold)
  base_half_angle_use <- max(
    min_base_half_angle,
    base_half_angle - n_extra * loop_compact_step
  )
  
  angle_step <- angle_step_scale * track_spacing / max(r_contact, 1e-8)
  half_angle <- min(base_half_angle_use + offset_index * angle_step, max_half_angle)
  
  n_start <- cos(half_angle) * u_out + sin(half_angle) * t_hat
  n_end   <- cos(half_angle) * u_out - sin(half_angle) * t_hat
  
  p0 <- c(s$x, s$y) + r_contact * n_start
  p6 <- c(s$x, s$y) + r_contact * n_end
  
  apex_dist <- r_contact +
    self_loop_height * s$radius_plot +
    offset_index * track_spacing * 1.35
  
  apex <- c(s$x, s$y) + apex_dist * u_out
  
  h_contact <- 0.95 * (
    self_loop_height * s$radius_plot +
      offset_index * track_spacing * 0.85
  )
  
  p1 <- p0 + h_contact * n_start
  p5 <- p6 + h_contact * n_end
  
  shoulder <- 1.35 * self_loop_spread * s$radius_plot +
    offset_index * track_spacing
  
  p2 <- apex + shoulder * t_hat
  p4 <- apex - shoulder * t_hat
  
  tt1 <- seq(0, 1, length.out = ceiling(n_pts / 2))
  tt2 <- seq(0, 1, length.out = n_pts - length(tt1) + 1L)
  
  seg1 <- .cubic_bezier(p0, p1, p2, apex, tt1)
  seg2 <- .cubic_bezier(apex, p4, p5, p6, tt2)
  
  loop_xy <- rbind(
    seg1,
    seg2[-1, , drop = FALSE]
  )
  
  out <- data.frame(
    x = loop_xy[, "x"],
    y = loop_xy[, "y"],
    point_id = seq_len(nrow(loop_xy)),
    from_group = s$group,
    to_group = s$group,
    arrow_side = -1,
    is_self = TRUE,
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- NULL
  out
}

#' Expand base directional lanes into stacked tracks
#'
#' Generic engine used by family and compare circuit plotting.
#'
#' Required edge-table columns:
#' - from_group, to_group
#' - one rank/order column
#' - one offset column
#' - one count column giving total tracks on that directional edge
#'
#' Mode-specific metadata are attached through attach_metadata_fun().
buildStackedCircuitCurves <- function(base_curve_df,
                                      edge_tbl,
                                      node_tbl,
                                      track_spacing,
                                      attach_metadata_fun,
                                      track_order_cols,
                                      offset_col,
                                      count_col,
                                      id_start = 0L,
                                      edge_gap = 0.015,
                                      self_loop_spread = 0.62,
                                      self_loop_height = 0.5,
                                      n_pts = 401,
                                      self_anchor_frac = 0.16) {
  req_base <- c("from_group", "to_group", "point_id", "x", "y")
  req_edge_core <- c("from_group", "to_group")
  req_node <- c("group", "x", "y", "theta", "radius_plot")
  
  if (!is.data.frame(base_curve_df) || !all(req_base %in% names(base_curve_df))) {
    stop("base_curve_df must contain: ", paste(req_base, collapse = ", "))
  }
  if (!is.data.frame(edge_tbl) || !all(req_edge_core %in% names(edge_tbl))) {
    stop("edge_tbl must contain: ", paste(req_edge_core, collapse = ", "))
  }
  if (!is.data.frame(node_tbl) || !all(req_node %in% names(node_tbl))) {
    stop("node_tbl must contain: ", paste(req_node, collapse = ", "))
  }
  if (!is.numeric(track_spacing) || length(track_spacing) != 1L ||
      !is.finite(track_spacing) || track_spacing < 0) {
    stop("track_spacing must be a single non-negative number.")
  }
  if (!is.function(attach_metadata_fun)) {
    stop("attach_metadata_fun must be a function.")
  }
  if (!is.character(track_order_cols) || length(track_order_cols) < 1L) {
    stop("track_order_cols must be a non-empty character vector.")
  }
  if (!all(track_order_cols %in% names(edge_tbl))) {
    stop(
      "edge_tbl is missing track_order_cols: ",
      paste(setdiff(track_order_cols, names(edge_tbl)), collapse = ", ")
    )
  }
  if (!is.character(offset_col) || length(offset_col) != 1L || !(offset_col %in% names(edge_tbl))) {
    stop("offset_col must name a column in edge_tbl.")
  }
  if (!is.character(count_col) || length(count_col) != 1L || !(count_col %in% names(edge_tbl))) {
    stop("count_col must name a column in edge_tbl.")
  }
  if (!is.numeric(id_start) || length(id_start) != 1L || !is.finite(id_start)) {
    stop("id_start must be a single finite numeric value.")
  }
  if (!is.numeric(edge_gap) || length(edge_gap) != 1L ||
      !is.finite(edge_gap) || edge_gap < 0) {
    stop("edge_gap must be a single non-negative number.")
  }
  if (!is.numeric(self_loop_spread) || length(self_loop_spread) != 1L ||
      !is.finite(self_loop_spread) || self_loop_spread <= 0) {
    stop("self_loop_spread must be a single positive number.")
  }
  if (!is.numeric(self_loop_height) || length(self_loop_height) != 1L ||
      !is.finite(self_loop_height) || self_loop_height <= 0) {
    stop("self_loop_height must be a single positive number.")
  }
  if (!is.numeric(n_pts) || length(n_pts) != 1L ||
      !is.finite(n_pts) || n_pts < 50) {
    stop("n_pts must be a single integer >= 50.")
  }
  if (!is.numeric(self_anchor_frac) || length(self_anchor_frac) != 1L ||
      !is.finite(self_anchor_frac) || self_anchor_frac <= 0 || self_anchor_frac >= 0.5) {
    stop("self_anchor_frac must be a single number in (0, 0.5).")
  }
  n_pts <- as.integer(n_pts)
  
  if (nrow(edge_tbl) == 0L || nrow(base_curve_df) == 0L) {
    return(data.frame())
  }
  
  curve_key <- interaction(base_curve_df$from_group, base_curve_df$to_group,
                           drop = TRUE, lex.order = TRUE)
  curve_split <- split(base_curve_df, curve_key)
  
  node_split <- split(node_tbl, node_tbl$group)
  
  edge_tbl$base_edge_id <- interaction(edge_tbl$from_group, edge_tbl$to_group,
                                       drop = TRUE, lex.order = TRUE)
  edge_split <- split(edge_tbl, edge_tbl$base_edge_id)
  
  out <- vector("list", length = 0L)
  track_id <- as.integer(id_start)
  
  for (key in names(edge_split)) {
    ed <- edge_split[[key]]
    
    ord_args <- unname(ed[track_order_cols])
    ord_args <- c(ord_args, list(method = "radix"))
    ord_idx <- do.call(order, ord_args)
    ed <- ed[ord_idx, , drop = FALSE]
    
    fg <- as.character(ed$from_group[1])
    tg <- as.character(ed$to_group[1])
    is_self <- identical(fg, tg)
    
    if (!is_self) {
      base_curve <- curve_split[[key]]
      if (is.null(base_curve) || nrow(base_curve) < 2L) next
      
      for (i in seq_len(nrow(ed))) {
        erow <- ed[i, , drop = FALSE]
        
        d_off <- track_spacing * erow[[offset_col]]
        geom_df <- offsetStackedCircuitCurve(
          curve_df = base_curve,
          offset_distance = d_off
        )
        
        keep_cols <- intersect(
          c("x", "y", "point_id", "from_group", "to_group", "arrow_side", "is_self"),
          names(geom_df)
        )
        geom_df <- geom_df[, keep_cols, drop = FALSE]
        rownames(geom_df) <- NULL
        
        track_id <- track_id + 1L
        out[[length(out) + 1L]] <- attach_metadata_fun(
          curve_df = geom_df,
          edge_row = erow,
          track_id = track_id
        )
      }
      
    } else {
      node_row <- node_split[[fg]]
      if (is.null(node_row) || nrow(node_row) != 1L) next
      
      prev_geom <- NULL
      
      for (i in seq_len(nrow(ed))) {
        erow <- ed[i, , drop = FALSE]
        oi <- erow[[offset_col]]
        n_tracks <- erow[[count_col]]
        
        target_geom <- buildStackedCircuitSelfLoopCurve(
          node_row = node_row,
          offset_index = oi,
          n_tracks_on_edge = n_tracks,
          edge_gap = edge_gap,
          self_loop_spread = self_loop_spread,
          self_loop_height = self_loop_height,
          track_spacing = track_spacing,
          n_pts = n_pts
        )
        
        if (is.null(prev_geom)) {
          geom_df <- target_geom
        } else {
          geom_df <- offsetStackedCircuitCurve(
            curve_df = prev_geom,
            offset_distance = track_spacing
          )
          
          keep_cols <- intersect(
            c("x", "y", "point_id", "from_group", "to_group", "arrow_side", "is_self"),
            names(geom_df)
          )
          geom_df <- geom_df[, keep_cols, drop = FALSE]
          rownames(geom_df) <- NULL
          
          geom_df <- anchorCurveToTargetEndpoints(
            curve_df = geom_df,
            target_df = target_geom,
            anchor_frac = self_anchor_frac
          )
          
          if (!"arrow_side" %in% names(geom_df)) geom_df$arrow_side <- -1
          if (!"is_self" %in% names(geom_df)) geom_df$is_self <- TRUE
        }
        
        track_id <- track_id + 1L
        out[[length(out) + 1L]] <- attach_metadata_fun(
          curve_df = geom_df,
          edge_row = erow,
          track_id = track_id
        )
        
        prev_geom <- geom_df
      }
    }
  }
  
  if (length(out) == 0L) {
    return(data.frame())
  }
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Build a manual identity legend using colored / linetyped line exemplars
#'
#' Generic builder for legends like:
#' - mechanism legend
#' - object legend
buildCircuitIdentityLegend <- function(df,
                                       label_col,
                                       color_col,
                                       linetype_col,
                                       legend_title = "Identity",
                                       edge_linewidth = 0.7,
                                       text_size = 3.4,
                                       title_size = 11,
                                       x_label = 1.35,
                                       xlim = c(0, 3.2),
                                       sort_labels = TRUE) {
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  
  req <- c(label_col, color_col, linetype_col)
  if (!all(req %in% names(df))) {
    stop("df must contain: ", paste(req, collapse = ", "))
  }
  
  leg_df <- unique(df[, req, drop = FALSE])
  
  if (isTRUE(sort_labels)) {
    ord <- order(as.character(leg_df[[label_col]]))
    leg_df <- leg_df[ord, , drop = FALSE]
  }
  
  rownames(leg_df) <- NULL
  leg_df$y <- rev(seq_len(nrow(leg_df)))
  
  ggplot2::ggplot(leg_df) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = 1.1, y = .data$y, yend = .data$y),
      color = leg_df[[color_col]],
      linetype = leg_df[[linetype_col]],
      linewidth = edge_linewidth,
      lineend = "round"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_label, y = .data$y, label = .data[[label_col]]),
      hjust = 0,
      color = "black",
      size = text_size
    ) +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(0.5, max(leg_df$y) + 0.8) +
    ggplot2::labs(title = legend_title) +
    makeCircuitLegendTheme(title_size = title_size, text_size = text_size)
}

#' Normalize a set of identity labels
#'
#' Removes NA / empty entries while preserving first appearance order.
normalizeIdentityLabels <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  unique(x)
}

#' Resolve a named or unnamed aesthetic vector against identity labels
#'
#' Supports:
#' - NULL -> caller handles defaults separately
#' - named vector with one entry per identity label
#' - unnamed vector with length equal to n_labels
resolveIdentityAestheticVector <- function(labels,
                                           values,
                                           arg_name = "values") {
  labels <- normalizeIdentityLabels(labels)
  
  if (length(labels) < 1L) {
    stop("No identity labels available for ", arg_name, " resolution.")
  }
  
  if (is.null(values)) {
    return(NULL)
  }
  
  if (!is.character(values)) {
    stop(arg_name, " must be NULL or a character vector.")
  }
  
  if (!is.null(names(values))) {
    miss <- setdiff(labels, names(values))
    if (length(miss) > 0L) {
      stop(
        "Named ", arg_name, " is missing entries for: ",
        paste(miss, collapse = ", "), "."
      )
    }
    out <- values[labels]
    return(out)
  }
  
  if (length(values) != length(labels)) {
    stop(
      "Unnamed ", arg_name, " must have length equal to the number of identities (",
      length(labels), ")."
    )
  }
  
  out <- values
  names(out) <- labels
  out
}

#' Resolve an aesthetic vector that may be provided at either group level
#' or identity level
#'
#' This is useful when identities are nested within groups, e.g.:
#' - mechanisms within ligands
#' - mechanisms within receptor families
resolveGroupedIdentityAestheticVector <- function(labels,
                                                  groups,
                                                  values,
                                                  arg_name = "values") {
  labels <- as.character(labels)
  groups <- as.character(groups)
  
  if (length(labels) != length(groups)) {
    stop("labels and groups must have the same length.")
  }
  
  keep <- !is.na(labels) & nzchar(labels)
  labels_keep <- labels[keep]
  groups_keep <- groups[keep]
  
  if (length(labels_keep) < 1L) {
    stop("No identity labels available for ", arg_name, " resolution.")
  }
  
  # Preserve first appearance order of labels
  first_idx <- !duplicated(labels_keep)
  lab_full <- labels_keep[first_idx]
  grp_full <- groups_keep[first_idx]
  
  # Sanity check: each identity label should map to exactly one group
  grp_by_label <- split(groups_keep, labels_keep)
  bad_labels <- names(grp_by_label)[vapply(
    grp_by_label,
    function(x) length(unique(x[!is.na(x)])) > 1L,
    logical(1)
  )]
  
  if (length(bad_labels) > 0L) {
    stop(
      "Each identity label must map to exactly one group. Ambiguous group mapping for: ",
      paste(bad_labels, collapse = ", "), "."
    )
  }
  
  group_levels <- unique(grp_full)
  
  if (is.null(values)) {
    return(NULL)
  }
  
  if (!is.character(values)) {
    stop(arg_name, " must be NULL or a character vector.")
  }
  
  if (!is.null(names(values))) {
    # First try full identity-level mapping
    miss_labels <- setdiff(lab_full, names(values))
    if (length(miss_labels) == 0L) {
      out <- values[lab_full]
      names(out) <- lab_full
      return(out)
    }
    
    # Then try group-level mapping
    miss_groups <- setdiff(group_levels, names(values))
    if (length(miss_groups) == 0L) {
      group_map <- values[group_levels]
      names(group_map) <- group_levels
      
      out <- group_map[grp_full]
      names(out) <- lab_full
      return(out)
    }
    
    stop(
      "Named ", arg_name, " must cover either all identities or all groups."
    )
  }
  
  if (length(values) == length(group_levels)) {
    group_map <- values
    names(group_map) <- group_levels
    
    out <- group_map[grp_full]
    names(out) <- lab_full
    return(out)
  }
  
  if (length(values) == length(lab_full)) {
    out <- values
    names(out) <- lab_full
    return(out)
  }
  
  stop(
    "Unnamed ", arg_name, " must have length equal to either ",
    "the number of identities or the number of groups."
  )
}

#' Default discrete color map for identity labels
defaultIdentityColors <- function(labels,
                                  palette = "Set 2") {
  labels <- normalizeIdentityLabels(labels)
  
  if (length(labels) < 1L) {
    stop("No identity labels available for default color generation.")
  }
  
  cols <- grDevices::hcl.colors(length(labels), palette = palette)
  names(cols) <- labels
  cols
}

#' Default object-comparison color map
defaultCompareObjectColors <- function(object_names) {
  object_names <- normalizeIdentityLabels(object_names)
  
  default_cols <- c("#1F3A93", "#8B1E1E")
  if (length(object_names) <= 2L) {
    cols <- default_cols[seq_along(object_names)]
  } else {
    extra <- grDevices::hcl.colors(length(object_names) - 2L, palette = "Dark 3")
    cols <- c(default_cols, extra)
  }
  names(cols) <- object_names
  cols
}

#' Default object-comparison linetype map
defaultCompareObjectLinetypes <- function(object_names) {
  object_names <- normalizeIdentityLabels(object_names)
  ltys <- rep("solid", length(object_names))
  names(ltys) <- object_names
  ltys
}

prepareStackedCircuitRenderGeometry <- function(curve_df,
                                                order_cols,
                                                track_id_col = "track_id",
                                                arrow_head_length = 0.060,
                                                arrow_head_width = 0.022,
                                                taper_length = 0.060,
                                                taper_base_halfwidth = 0.0035) {
  if (!is.data.frame(curve_df)) {
    stop("curve_df must be a data.frame.")
  }
  if (!is.character(track_id_col) || length(track_id_col) != 1L || !nzchar(track_id_col)) {
    stop("track_id_col must be a single non-empty character string.")
  }
  if (!is.character(order_cols) || length(order_cols) < 1L) {
    stop("order_cols must be a non-empty character vector.")
  }
  if (!is.numeric(arrow_head_length) || length(arrow_head_length) != 1L ||
      !is.finite(arrow_head_length) || arrow_head_length <= 0) {
    stop("arrow_head_length must be a single positive number.")
  }
  if (!is.numeric(arrow_head_width) || length(arrow_head_width) != 1L ||
      !is.finite(arrow_head_width) || arrow_head_width <= 0) {
    stop("arrow_head_width must be a single positive number.")
  }
  if (!is.numeric(taper_length) || length(taper_length) != 1L ||
      !is.finite(taper_length) || taper_length < 0) {
    stop("taper_length must be a single non-negative number.")
  }
  if (!is.numeric(taper_base_halfwidth) || length(taper_base_halfwidth) != 1L ||
      !is.finite(taper_base_halfwidth) || taper_base_halfwidth <= 0) {
    stop("taper_base_halfwidth must be a single positive number.")
  }
  
  req <- c(track_id_col, "point_id", "x", "y", "weight_plot", "weight_raw")
  miss <- setdiff(req, names(curve_df))
  if (length(miss) > 0L) {
    stop("curve_df is missing required column(s): ", paste(miss, collapse = ", "))
  }
  
  miss_ord <- setdiff(order_cols, names(curve_df))
  if (length(miss_ord) > 0L) {
    stop("curve_df is missing order_cols: ", paste(miss_ord, collapse = ", "))
  }
  
  if (nrow(curve_df) == 0L) {
    empty <- curve_df[0, , drop = FALSE]
    empty$edge_id <- integer(0)
    empty$head_length_use <- numeric(0)
    empty$head_width_use <- numeric(0)
    empty$taper_length_use <- numeric(0)
    empty$taper_base_halfwidth_use <- numeric(0)
    
    return(list(
      curve_geom = empty,
      edge_order_df = data.frame(),
      shaft_df = empty,
      taper_poly = data.frame(),
      arrow_poly = data.frame()
    ))
  }
  
  ord_df <- curve_df[, order_cols, drop = FALSE]
  ord_args <- c(unname(ord_df), list(method = "radix"))
  ord_idx <- do.call(order, ord_args)
  curve_df <- curve_df[ord_idx, , drop = FALSE]
  rownames(curve_df) <- NULL
  
  curve_geom <- curve_df
  curve_geom$edge_id <- curve_geom[[track_id_col]]
  curve_geom$head_length_use <- arrow_head_length
  curve_geom$head_width_use <- arrow_head_width
  curve_geom$taper_length_use <- taper_length
  curve_geom$taper_base_halfwidth_use <- taper_base_halfwidth
  
  edge_order_df <- unique(
    curve_geom[, c("edge_id", "weight_plot", "weight_raw"), drop = FALSE]
  )
  edge_order_df <- edge_order_df[order(edge_order_df$weight_plot,
                                       edge_order_df$weight_raw,
                                       edge_order_df$edge_id), , drop = FALSE]
  rownames(edge_order_df) <- NULL
  edge_order_df$draw_rank <- seq_len(nrow(edge_order_df))
  
  curve_split <- split(curve_geom, curve_geom$edge_id)
  
  shaft_list <- lapply(curve_split, function(d) {
    d <- d[order(d$point_id), , drop = FALSE]
    
    trim_i <- unique(d$taper_length_use)
    trim_i <- trim_i[is.finite(trim_i)]
    if (length(trim_i) != 1L) trim_i <- taper_length
    
    out <- trimCircuitCurvesForArrowhead(
      curve_df = d,
      trim_length = trim_i
    )
    
    if (is.null(out) || nrow(out) < 2L) {
      return(NULL)
    }
    
    out
  })
  shaft_list <- shaft_list[!vapply(shaft_list, is.null, logical(1))]
  
  if (length(shaft_list) > 0L) {
    shaft_df <- do.call(rbind, shaft_list)
    rownames(shaft_df) <- NULL
    shaft_df$draw_rank <- edge_order_df$draw_rank[match(shaft_df$edge_id, edge_order_df$edge_id)]
    shaft_df <- shaft_df[order(shaft_df$draw_rank,
                               shaft_df$edge_id,
                               shaft_df$point_id), , drop = FALSE]
    rownames(shaft_df) <- NULL
  } else {
    shaft_df <- curve_geom[0, , drop = FALSE]
    shaft_df$draw_rank <- numeric(0)
  }
  
  taper_poly <- buildCircuitTerminalTaperPolygons(curve_df = curve_geom)
  arrow_poly <- buildCircuitArrowPolygons(curve_df = curve_geom)
  
  if (nrow(taper_poly) > 0L) {
    taper_poly$draw_rank <- edge_order_df$draw_rank[match(taper_poly$edge_id, edge_order_df$edge_id)]
    taper_poly$poly_group <- interaction(taper_poly$edge_id, taper_poly$part_id, drop = TRUE)
    taper_poly <- taper_poly[order(taper_poly$draw_rank,
                                   taper_poly$edge_id,
                                   taper_poly$part_id,
                                   taper_poly$point_id), , drop = FALSE]
    rownames(taper_poly) <- NULL
  }
  
  if (nrow(arrow_poly) > 0L) {
    arrow_poly$draw_rank <- edge_order_df$draw_rank[match(arrow_poly$edge_id, edge_order_df$edge_id)]
    arrow_poly$poly_group <- interaction(arrow_poly$edge_id, arrow_poly$part_id, drop = TRUE)
    arrow_poly <- arrow_poly[order(arrow_poly$draw_rank,
                                   arrow_poly$edge_id,
                                   arrow_poly$part_id,
                                   arrow_poly$point_id), , drop = FALSE]
    rownames(arrow_poly) <- NULL
  }
  
  list(
    curve_geom = curve_geom,
    edge_order_df = edge_order_df,
    shaft_df = shaft_df,
    taper_poly = taper_poly,
    arrow_poly = arrow_poly
  )
}

renderStackedCircuitTracksFixed <- function(p,
                                            edge_order_df,
                                            shaft_df,
                                            taper_poly,
                                            arrow_poly,
                                            color_col = "identity_color",
                                            linetype_col = "identity_linetype",
                                            edge_linewidth = 0.7) {
  if (!inherits(p, "ggplot")) {
    stop("p must be a ggplot object.")
  }
  if (!is.data.frame(edge_order_df) || !"edge_id" %in% names(edge_order_df)) {
    stop("edge_order_df must contain edge_id.")
  }
  if (!is.character(color_col) || length(color_col) != 1L || !nzchar(color_col)) {
    stop("color_col must be a single non-empty character string.")
  }
  if (!is.character(linetype_col) || length(linetype_col) != 1L || !nzchar(linetype_col)) {
    stop("linetype_col must be a single non-empty character string.")
  }
  if (!is.numeric(edge_linewidth) || length(edge_linewidth) != 1L ||
      !is.finite(edge_linewidth) || edge_linewidth <= 0) {
    stop("edge_linewidth must be a single positive number.")
  }
  
  req_shaft <- c("edge_id", "x", "y", "point_id", color_col, linetype_col)
  if (!all(req_shaft %in% names(shaft_df))) {
    stop("shaft_df must contain: ", paste(req_shaft, collapse = ", "))
  }
  
  if (nrow(taper_poly) > 0L && !all(c("edge_id", "x", "y", "poly_group") %in% names(taper_poly))) {
    stop("taper_poly must contain edge_id, x, y, and poly_group when non-empty.")
  }
  if (nrow(arrow_poly) > 0L && !all(c("edge_id", "x", "y", "poly_group") %in% names(arrow_poly))) {
    stop("arrow_poly must contain edge_id, x, y, and poly_group when non-empty.")
  }
  
  shaft_split <- if (nrow(shaft_df) > 0L) split(shaft_df, shaft_df$edge_id) else list()
  taper_split <- if (nrow(taper_poly) > 0L) split(taper_poly, taper_poly$edge_id) else list()
  arrow_split <- if (nrow(arrow_poly) > 0L) split(arrow_poly, arrow_poly$edge_id) else list()
  
  for (eid in edge_order_df$edge_id) {
    d_shaft <- shaft_split[[as.character(eid)]]
    if (!is.null(d_shaft) && nrow(d_shaft) > 1L) {
      lt <- unique(as.character(d_shaft[[linetype_col]]))
      lt <- lt[!is.na(lt) & nzchar(lt)]
      if (length(lt) != 1L) lt <- "solid"
      
      col_i <- unique(as.character(d_shaft[[color_col]]))
      col_i <- col_i[!is.na(col_i) & nzchar(col_i)]
      if (length(col_i) != 1L) col_i <- "black"
      
      p <- p +
        ggplot2::geom_path(
          data = d_shaft,
          ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id),
          color = col_i,
          linetype = lt,
          linewidth = edge_linewidth,
          lineend = "butt",
          linejoin = "round",
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
    
    d_taper <- taper_split[[as.character(eid)]]
    if (!is.null(d_taper) && nrow(d_taper) > 0L) {
      fill_i <- if (!is.null(d_shaft) && nrow(d_shaft) > 0L) {
        col_i
      } else {
        "black"
      }
      
      p <- p +
        ggplot2::geom_polygon(
          data = d_taper,
          ggplot2::aes(x = .data$x, y = .data$y, group = .data$poly_group),
          fill = fill_i,
          color = NA,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
    
    d_arrow <- arrow_split[[as.character(eid)]]
    if (!is.null(d_arrow) && nrow(d_arrow) > 0L) {
      fill_i <- if (!is.null(d_shaft) && nrow(d_shaft) > 0L) {
        col_i
      } else {
        "black"
      }
      
      p <- p +
        ggplot2::geom_polygon(
          data = d_arrow,
          ggplot2::aes(x = .data$x, y = .data$y, group = .data$poly_group),
          fill = fill_i,
          color = NA,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
  }
  
  p
}

renderStackedCircuitTracksMapped <- function(p,
                                             edge_order_df,
                                             shaft_df,
                                             taper_poly,
                                             arrow_poly,
                                             color_col = "identity_color",
                                             linetype_col = "identity_linetype",
                                             edge_aes = c("alpha", "width"),
                                             edge_alpha_range = c(0.01, 1),
                                             alpha_gamma = 2.2,
                                             edge_width_range = c(0.25, 2.6),
                                             edge_linewidth = 0.7) {
  edge_aes <- match.arg(edge_aes)
  
  if (!inherits(p, "ggplot")) {
    stop("p must be a ggplot object.")
  }
  if (!is.data.frame(edge_order_df) || !all(c("edge_id", "weight_plot") %in% names(edge_order_df))) {
    stop("edge_order_df must contain edge_id and weight_plot.")
  }
  if (!is.character(color_col) || length(color_col) != 1L || !nzchar(color_col)) {
    stop("color_col must be a single non-empty character string.")
  }
  if (!is.character(linetype_col) || length(linetype_col) != 1L || !nzchar(linetype_col)) {
    stop("linetype_col must be a single non-empty character string.")
  }
  if (!is.numeric(edge_alpha_range) || length(edge_alpha_range) != 2L ||
      any(!is.finite(edge_alpha_range)) || edge_alpha_range[1] < 0 ||
      edge_alpha_range[2] > 1 || edge_alpha_range[1] >= edge_alpha_range[2]) {
    stop("edge_alpha_range must be a numeric vector like c(0.01, 1).")
  }
  if (!is.numeric(alpha_gamma) || length(alpha_gamma) != 1L ||
      !is.finite(alpha_gamma) || alpha_gamma <= 0) {
    stop("alpha_gamma must be a single positive number.")
  }
  if (!is.numeric(edge_width_range) || length(edge_width_range) != 2L ||
      any(!is.finite(edge_width_range)) || edge_width_range[1] <= 0 ||
      edge_width_range[1] >= edge_width_range[2]) {
    stop("edge_width_range must be a numeric vector like c(0.25, 2.6).")
  }
  if (!is.numeric(edge_linewidth) || length(edge_linewidth) != 1L ||
      !is.finite(edge_linewidth) || edge_linewidth <= 0) {
    stop("edge_linewidth must be a single positive number.")
  }
  
  req_shaft <- c("edge_id", "x", "y", "point_id", "weight_plot", color_col, linetype_col)
  if (!all(req_shaft %in% names(shaft_df))) {
    stop("shaft_df must contain: ", paste(req_shaft, collapse = ", "))
  }
  
  if (nrow(taper_poly) > 0L && !all(c("edge_id", "x", "y", "weight_plot", "poly_group") %in% names(taper_poly))) {
    stop("taper_poly must contain edge_id, x, y, weight_plot, and poly_group when non-empty.")
  }
  if (nrow(arrow_poly) > 0L && !all(c("edge_id", "x", "y", "weight_plot", "poly_group") %in% names(arrow_poly))) {
    stop("arrow_poly must contain edge_id, x, y, weight_plot, and poly_group when non-empty.")
  }
  
  vals <- edge_order_df$weight_plot[is.finite(edge_order_df$weight_plot)]
  if (length(vals) < 1L) vals <- c(0, 1)
  scale_limits <- range(vals, na.rm = TRUE)
  if (diff(scale_limits) == 0) {
    eps <- max(abs(scale_limits[1]), 1) * 1e-9
    scale_limits <- scale_limits + c(-eps, eps)
  }
  
  shaft_split <- if (nrow(shaft_df) > 0L) split(shaft_df, shaft_df$edge_id) else list()
  taper_split <- if (nrow(taper_poly) > 0L) split(taper_poly, taper_poly$edge_id) else list()
  arrow_split <- if (nrow(arrow_poly) > 0L) split(arrow_poly, arrow_poly$edge_id) else list()
  
  if (edge_aes == "alpha") {
    remap_alpha_value <- function(x, limits, gamma = 2.2) {
      z <- scales::rescale(x, to = c(0, 1), from = limits)
      z[!is.finite(z)] <- 0
      z <- pmax(0, pmin(1, z))
      z ^ gamma
    }
    
    for (eid in edge_order_df$edge_id) {
      d_shaft <- shaft_split[[as.character(eid)]]
      d_taper <- taper_split[[as.character(eid)]]
      d_arrow <- arrow_split[[as.character(eid)]]
      
      a_i <- edge_order_df$weight_plot[match(eid, edge_order_df$edge_id)]
      a_i <- remap_alpha_value(a_i, limits = scale_limits, gamma = alpha_gamma)
      a_i <- edge_alpha_range[1] + a_i * diff(edge_alpha_range)
      
      col_i <- "black"
      lt <- "solid"
      
      if (!is.null(d_shaft) && nrow(d_shaft) > 1L) {
        lt <- unique(as.character(d_shaft[[linetype_col]]))
        lt <- lt[!is.na(lt) & nzchar(lt)]
        if (length(lt) != 1L) lt <- "solid"
        
        col_i <- unique(as.character(d_shaft[[color_col]]))
        col_i <- col_i[!is.na(col_i) & nzchar(col_i)]
        if (length(col_i) != 1L) col_i <- "black"
        
        p <- p +
          ggplot2::geom_path(
            data = d_shaft,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id),
            color = col_i,
            alpha = a_i,
            linetype = lt,
            linewidth = edge_linewidth,
            lineend = "butt",
            linejoin = "round",
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
      
      if (!is.null(d_taper) && nrow(d_taper) > 0L) {
        p <- p +
          ggplot2::geom_polygon(
            data = d_taper,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$poly_group),
            fill = col_i,
            alpha = a_i,
            color = NA,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
      
      if (!is.null(d_arrow) && nrow(d_arrow) > 0L) {
        p <- p +
          ggplot2::geom_polygon(
            data = d_arrow,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$poly_group),
            fill = col_i,
            alpha = a_i,
            color = NA,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
    }
  } else {
    for (eid in edge_order_df$edge_id) {
      d_shaft <- shaft_split[[as.character(eid)]]
      d_taper <- taper_split[[as.character(eid)]]
      d_arrow <- arrow_split[[as.character(eid)]]
      
      lw_i <- edge_order_df$weight_plot[match(eid, edge_order_df$edge_id)]
      lw_i <- scales::rescale(lw_i, to = edge_width_range, from = scale_limits)
      
      col_i <- "black"
      lt <- "solid"
      
      if (!is.null(d_shaft) && nrow(d_shaft) > 1L) {
        lt <- unique(as.character(d_shaft[[linetype_col]]))
        lt <- lt[!is.na(lt) & nzchar(lt)]
        if (length(lt) != 1L) lt <- "solid"
        
        col_i <- unique(as.character(d_shaft[[color_col]]))
        col_i <- col_i[!is.na(col_i) & nzchar(col_i)]
        if (length(col_i) != 1L) col_i <- "black"
        
        p <- p +
          ggplot2::geom_path(
            data = d_shaft,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id),
            color = col_i,
            linetype = lt,
            linewidth = lw_i,
            lineend = "butt",
            linejoin = "round",
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
      
      if (!is.null(d_taper) && nrow(d_taper) > 0L) {
        p <- p +
          ggplot2::geom_polygon(
            data = d_taper,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$poly_group),
            fill = col_i,
            color = NA,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
      
      if (!is.null(d_arrow) && nrow(d_arrow) > 0L) {
        p <- p +
          ggplot2::geom_polygon(
            data = d_arrow,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$poly_group),
            fill = col_i,
            color = NA,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
    }
  }
  
  p
}

validateStackedRenderColumns <- function(curve_df) {
  req <- c(
    "track_id",
    "identity",
    "identity_color",
    "identity_linetype",
    "weight_plot",
    "weight_raw",
    "point_id",
    "x",
    "y",
    "from_group",
    "to_group",
    "arrow_side",
    "is_self"
  )
  
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  
  invisible(TRUE)
}

renderCircuitNodesBasic <- function(p,
                                    node_tbl,
                                    node_palette = NULL,
                                    node_outline = "black",
                                    label_size = 4) {
  req <- c("group", "x", "y", "radius_plot", "label_x", "label_y", "hjust")
  if (!inherits(p, "ggplot")) {
    stop("p must be a ggplot object.")
  }
  if (!is.data.frame(node_tbl) || !all(req %in% names(node_tbl))) {
    stop("node_tbl must contain: ", paste(req, collapse = ", "))
  }
  
  pal <- resolveCircuitNodePalette(node_tbl$group, node_palette = node_palette)
  node_df <- node_tbl
  node_df$fill_col <- unname(pal[node_df$group])
  
  p +
    ggforce::geom_circle(
      data = node_df,
      ggplot2::aes(
        x0 = .data$x,
        y0 = .data$y,
        r = .data$radius_plot,
        fill = .data$fill_col
      ),
      color = node_outline,
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(
      data = node_df,
      ggplot2::aes(
        x = .data$label_x,
        y = .data$label_y,
        label = .data$group,
        hjust = .data$hjust
      ),
      color = "black",
      size = label_size,
      show.legend = FALSE
    )
}

renderCircuitNodesCompareGlyph <- function(p,
                                           node_tbl,
                                           object_aesthetics = NULL,
                                           node_palette = NULL,
                                           node_outline = "black",
                                           label_size = 4) {
  req <- c(
    "group", "x", "y", "label_x", "label_y", "hjust",
    "radius_outer", "radius_inner", "outer_object", "inner_object"
  )
  if (!inherits(p, "ggplot")) {
    stop("p must be a ggplot object.")
  }
  if (!is.data.frame(node_tbl) || !all(req %in% names(node_tbl))) {
    stop("node_tbl must contain: ", paste(req, collapse = ", "))
  }
  
  pal <- resolveCircuitNodePalette(node_tbl$group, node_palette = node_palette)
  node_df <- node_tbl
  node_df$fill_col <- unname(pal[node_df$group])
  
  obj_names_needed <- unique(c(node_df$outer_object, node_df$inner_object))
  obj_names_needed <- obj_names_needed[
    !is.na(obj_names_needed) &
      nzchar(obj_names_needed) &
      obj_names_needed != "neutral"
  ]
  
  if (is.null(object_aesthetics)) {
    object_aesthetics <- resolveCompareObjectAesthetics(obj_names_needed)
  }
  
  node_df$outer_col <- ifelse(
    node_df$outer_object == "neutral",
    node_outline,
    unname(object_aesthetics$colors[as.character(node_df$outer_object)])
  )
  
  node_df$inner_col <- unname(object_aesthetics$colors[as.character(node_df$inner_object)])
  
  outer_nodes <- node_df[
    is.finite(node_df$radius_outer) &
      node_df$radius_outer > 0 &
      !is.na(node_df$outer_col) &
      nzchar(node_df$outer_col),
    , drop = FALSE
  ]
  
  inner_nodes <- node_df[
    is.finite(node_df$radius_inner) &
      node_df$radius_inner > 0 &
      !is.na(node_df$inner_col) &
      nzchar(node_df$inner_col) &
      (
        !("is_neutral" %in% names(node_df)) |
          is.na(node_df$is_neutral) |
          !node_df$is_neutral
      ),
    , drop = FALSE
  ]
  
  if (nrow(outer_nodes) > 0L) {
    p <- p +
      ggforce::geom_circle(
        data = outer_nodes,
        ggplot2::aes(
          x0 = .data$x,
          y0 = .data$y,
          r = .data$radius_outer,
          fill = .data$fill_col,
          color = .data$outer_col
        ),
        linewidth = 0.45,
        show.legend = FALSE
      )
  }
  
  if (nrow(inner_nodes) > 0L) {
    p <- p +
      ggforce::geom_circle(
        data = inner_nodes,
        ggplot2::aes(
          x0 = .data$x,
          y0 = .data$y,
          r = .data$radius_inner,
          color = .data$inner_col
        ),
        fill = NA,
        linewidth = 0.45,
        show.legend = FALSE
      )
  }
  
  p +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::geom_text(
      data = node_df,
      ggplot2::aes(
        x = .data$label_x,
        y = .data$label_y,
        label = .data$group,
        hjust = .data$hjust
      ),
      color = "black",
      size = label_size,
      show.legend = FALSE
    )
}

useCompareCircuitNodeGlyph <- function(node_tbl) {
  req <- c(
    "radius_obj1", "radius_obj2",
    "radius_outer", "radius_inner",
    "outer_object", "inner_object"
  )
  
  if (!is.data.frame(node_tbl) || !all(req %in% names(node_tbl))) {
    return(FALSE)
  }
  
  outer_ok <- is.finite(node_tbl$radius_outer)
  inner_ok <- is.finite(node_tbl$radius_inner)
  obj_ok <- !is.na(node_tbl$outer_object) & nzchar(node_tbl$outer_object) &
    !is.na(node_tbl$inner_object) & nzchar(node_tbl$inner_object)
  
  any(outer_ok & inner_ok & obj_ok)
}

assembleStackedCircuitLegends <- function(node_tbl,
                                          node_palette = NULL,
                                          node_outline = "black",
                                          signal_legend = NULL,
                                          identity_legend = NULL,
                                          show_node_size_legend = NULL,
                                          node_size_legend_title = "Group abundance",
                                          node_size_legend_range = c(4, 10),
                                          group_legend_title = "Group",
                                          heights_no_size = c(1, 1.25, 1),
                                          heights_with_size = c(1, 1.25, 1, 1)) {
  if (!is.data.frame(node_tbl) || !"group" %in% names(node_tbl)) {
    stop("node_tbl must contain a 'group' column.")
  }
  
  if (is.null(show_node_size_legend)) {
    if (!"radius_plot" %in% names(node_tbl)) {
      show_node_size_legend <- FALSE
    } else {
      show_node_size_legend <- length(unique(round(node_tbl$radius_plot, 10))) > 1L
    }
  }
  
  if (is.null(signal_legend)) {
    signal_legend <- ggplot2::ggplot() + makeCircuitLegendTheme()
  }
  if (is.null(identity_legend)) {
    identity_legend <- ggplot2::ggplot() + makeCircuitLegendTheme()
  }
  
  group_legend <- buildCircuitGroupLegend(
    node_tbl = node_tbl,
    node_palette = node_palette,
    node_outline = node_outline,
    legend_title = group_legend_title
  )
  
  size_legend <- NULL
  if (isTRUE(show_node_size_legend)) {
    size_legend <- buildCircuitNodeSizeLegend(
      node_tbl = node_tbl,
      node_size_legend_title = node_size_legend_title,
      node_size_legend_range = node_size_legend_range
    )
  }
  
  assembleCircuitLegendColumn(
    signal_legend = signal_legend,
    identity_legend = identity_legend,
    group_legend = group_legend,
    size_legend = size_legend,
    show_node_size_legend = show_node_size_legend,
    heights_no_size = heights_no_size,
    heights_with_size = heights_with_size
  )
}

makeEmptyCircuitLegend <- function() {
  ggplot2::ggplot() + makeCircuitLegendTheme()
}

#### New Legend Helpers (LEGENDS STILL NEED WORK) ####

#' Shared theme for manual legend panels
makeCircuitLegendTheme <- function(title_size = 11,
                                   text_size = 3.3) {
  ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(
        hjust = 0,
        color = "black",
        size = title_size
      ),
      text = ggplot2::element_text(color = "black", size = text_size),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

#' Build manual group legend matching circuit node appearance
buildCircuitGroupLegend <- function(node_tbl,
                                    node_palette = NULL,
                                    node_outline = "black",
                                    legend_title = "Group",
                                    point_size = 5.5,
                                    text_size = 3.3,
                                    title_size = 11,
                                    x_label = 1.15,
                                    xlim = c(0, 2.8)) {
  if (!is.data.frame(node_tbl) || !"group" %in% names(node_tbl)) {
    stop("node_tbl must contain a 'group' column.")
  }
  
  pal <- resolveCircuitNodePalette(node_tbl$group, node_palette = node_palette)
  
  grp_df <- unique(node_tbl[, "group", drop = FALSE])
  grp_df <- grp_df[order(grp_df$group), , drop = FALSE]
  rownames(grp_df) <- NULL
  
  grp_df$fill <- unname(pal[grp_df$group])
  grp_df$y <- rev(seq_len(nrow(grp_df)))
  
  ggplot2::ggplot(grp_df) +
    ggplot2::geom_point(
      ggplot2::aes(x = 0.5, y = .data$y),
      shape = 21,
      size = point_size,
      fill = grp_df$fill,
      color = node_outline,
      stroke = 0.45
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_label, y = .data$y, label = .data$group),
      hjust = 0,
      color = "black",
      size = text_size
    ) +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(0.5, max(grp_df$y) + 0.8) +
    ggplot2::labs(title = legend_title) +
    makeCircuitLegendTheme(title_size = title_size, text_size = text_size)
}

#' Build manual node-size legend
#'
#' Updated behavior:
#' - always shows exactly 3 circles
#' - smallest node in the plot
#' - largest node in the plot
#' - arithmetic mean of those two node sizes
#' - labels are the corresponding population percentages
#'   reported to 3 significant figures
#' - legend circles are drawn with ggforce::geom_circle() using the
#'   same radius_plot scale as the main plot, so size ratios match
#'   the actual nodes
#'
#' Note:
#' - node_size_legend_range is retained only for backward compatibility
#'   and is not used, because legend circles are now drawn directly in
#'   data units rather than point-size units.
buildCircuitNodeSizeLegend <- function(node_tbl,
                                       node_size_legend_title = "Group abundance",
                                       node_size_legend_range = c(4, 10),
                                       fill = "grey75",
                                       color = "black",
                                       stroke = 0.45,
                                       text_size = 3.2,
                                       title_size = 11,
                                       x_label = NULL,
                                       x_padding_left = 0.15,
                                       x_padding_right = 1.35,
                                       y_gap_factor = 0.55) {
  req <- c("radius_plot", "frac_cells")
  if (!is.data.frame(node_tbl) || !all(req %in% names(node_tbl))) {
    stop("node_tbl must contain: ", paste(req, collapse = ", "))
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Package 'ggforce' is required.")
  }
  
  # retained for backward compatibility; not used
  invisible(node_size_legend_range)
  
  rp <- node_tbl$radius_plot
  fp <- node_tbl$frac_cells
  
  ok <- is.finite(rp) & is.finite(fp)
  if (!any(ok)) {
    stop("node_tbl must contain at least one finite radius_plot / frac_cells pair.")
  }
  
  rp <- rp[ok]
  fp <- fp[ok]
  
  i_min <- which.min(rp)
  i_max <- which.max(rp)
  
  radius_min <- rp[i_min]
  radius_max <- rp[i_max]
  frac_min <- fp[i_min]
  frac_max <- fp[i_max]
  
  radius_mid <- mean(c(radius_min, radius_max))
  frac_mid <- mean(c(frac_min, frac_max))
  
  legend_radius <- c(radius_min, radius_mid, radius_max)
  legend_frac <- c(frac_min, frac_mid, frac_max)
  
  # Display largest at top, midpoint in middle, smallest at bottom
  ord <- c(3, 2, 1)
  legend_radius <- legend_radius[ord]
  legend_frac <- legend_frac[ord]
  
  r_max <- max(legend_radius, na.rm = TRUE)
  if (!is.finite(r_max) || r_max <= 0) {
    stop("Could not determine a valid maximum legend radius.")
  }
  
  # Vertical placement of circle centers with adequate spacing
  y_bot <- r_max
  y_mid <- y_bot + r_max + y_gap_factor * r_max
  y_top <- y_mid + r_max + y_gap_factor * r_max
  
  legend_df <- data.frame(
    x0 = rep(r_max + x_padding_left, 3),
    y0 = c(y_top, y_mid, y_bot),
    r = legend_radius,
    label = paste0(signif(100 * legend_frac, 3), "%"),
    stringsAsFactors = FALSE
  )
  
  if (is.null(x_label)) {
    x_label <- 2 * r_max + x_padding_right
  }
  
  x_min <- min(legend_df$x0 - legend_df$r) - 0.05 * r_max
  x_max <- x_label + 0.90 * r_max
  y_min <- min(legend_df$y0 - legend_df$r) - 0.25 * r_max
  y_max <- max(legend_df$y0 + legend_df$r) + 0.25 * r_max
  
  ggplot2::ggplot(legend_df) +
    ggforce::geom_circle(
      ggplot2::aes(x0 = .data$x0, y0 = .data$y0, r = .data$r),
      fill = fill,
      color = color,
      linewidth = stroke
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_label, y = .data$y0, label = .data$label),
      hjust = 0,
      color = "black",
      size = text_size
    ) +
    ggplot2::coord_equal(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      clip = "off"
    ) +
    ggplot2::labs(title = node_size_legend_title) +
    makeCircuitLegendTheme(title_size = title_size, text_size = text_size)
}

#' Build manual signal legend for alpha-mode LR circuit plots
#'
#' Shows grayscale swatches matching the whitening behavior used on-plot.
buildCircuitSignalLegendAlpha <- function(values,
                                          legend_title = "Signal strength",
                                          alpha_range = c(0.01, 1),
                                          alpha_gamma = 2.2,
                                          n = 5,
                                          text_size = 3.2,
                                          title_size = 11,
                                          x_label = 1.15,
                                          xlim = c(0, 2.3)) {
  if (!is.numeric(values) || length(values) < 1L) {
    stop("values must be a non-empty numeric vector.")
  }
  vals <- values[is.finite(values)]
  if (length(vals) < 1L) vals <- c(0, 1)
  
  vr <- range(vals, na.rm = TRUE)
  if (!all(is.finite(vr))) vr <- c(0, 1)
  if (diff(vr) == 0) vr <- vr + c(-1e-9, 1e-9)
  
  if (!is.numeric(alpha_range) || length(alpha_range) != 2L ||
      any(!is.finite(alpha_range)) || alpha_range[1] < 0 ||
      alpha_range[2] > 1 || alpha_range[1] >= alpha_range[2]) {
    stop("alpha_range must be a numeric vector like c(0.01, 1).")
  }
  
  remap_alpha_value <- function(x, limits, gamma = 2.2) {
    z <- scales::rescale(x, to = c(0, 1), from = limits)
    z[!is.finite(z)] <- 0
    z <- pmax(0, pmin(1, z))
    z ^ gamma
  }
  
  sig_vals <- seq(vr[1], vr[2], length.out = n)
  z <- remap_alpha_value(sig_vals, limits = vr, gamma = alpha_gamma)
  alpha_vals <- alpha_range[1] + z * diff(alpha_range)
  alpha_vals <- pmax(alpha_range[1], pmin(alpha_range[2], alpha_vals))
  
  sig_cols <- vapply(
    alpha_vals,
    function(a) grDevices::gray(1 - a),
    character(1)
  )
  
  sig_df <- data.frame(
    y = rev(seq_along(sig_vals)),
    label = signif(rev(sig_vals), 3),
    col = rev(sig_cols),
    stringsAsFactors = FALSE
  )
  
  ggplot2::ggplot(sig_df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 0.5, y = .data$y),
      width = 0.7,
      height = 0.8,
      fill = sig_df$col,
      color = "black",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_label, y = .data$y, label = .data$label),
      hjust = 0,
      color = "black",
      size = text_size
    ) +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(0.5, max(sig_df$y) + 0.8) +
    ggplot2::labs(title = legend_title) +
    makeCircuitLegendTheme(title_size = title_size, text_size = text_size)
}

#' Build manual signal legend for width-mode LR circuit plots
buildCircuitSignalLegendWidth <- function(values,
                                          legend_title = "Signal strength",
                                          edge_width_range = c(0.25, 2.6),
                                          edge_color = "black",
                                          edge_alpha = 0.9,
                                          n = 5,
                                          text_size = 3.2,
                                          title_size = 11,
                                          x_label = 1.35,
                                          xlim = c(0, 2.6)) {
  if (!is.numeric(values) || length(values) < 1L) {
    stop("values must be a non-empty numeric vector.")
  }
  vals <- values[is.finite(values)]
  if (length(vals) < 1L) vals <- c(0, 1)
  
  vr <- range(vals, na.rm = TRUE)
  if (!all(is.finite(vr))) vr <- c(0, 1)
  if (diff(vr) == 0) vr <- vr + c(-1e-9, 1e-9)
  
  sig_vals <- seq(vr[1], vr[2], length.out = n)
  widths <- scales::rescale(sig_vals, to = edge_width_range, from = vr)
  
  sig_df <- data.frame(
    y = rev(seq_along(sig_vals)),
    label = signif(rev(sig_vals), 3),
    linewidth = rev(widths),
    stringsAsFactors = FALSE
  )
  
  ggplot2::ggplot(sig_df) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0.15, xend = 1.0,
        y = .data$y, yend = .data$y,
        linewidth = .data$linewidth
      ),
      color = edge_color,
      alpha = edge_alpha,
      lineend = "butt"
    ) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::geom_text(
      ggplot2::aes(x = x_label, y = .data$y, label = .data$label),
      hjust = 0,
      color = "black",
      size = text_size
    ) +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(0.5, max(sig_df$y) + 0.8) +
    ggplot2::labs(title = legend_title) +
    makeCircuitLegendTheme(title_size = title_size, text_size = text_size)
}

#' Build manual edge-style legend for LR circuit plots
#'
#' For alpha mode, shows a representative curve in grayscale.
#' For width mode, shows a representative curve at median width.
buildCircuitEdgeStyleLegend <- function(edge_aes = c("alpha", "width"),
                                        legend_title = "Edge encoding",
                                        edge_color = "black",
                                        edge_alpha = 0.9,
                                        edge_linewidth = 0.7,
                                        edge_alpha_range = c(0.01, 1),
                                        alpha_gamma = 2.2,
                                        edge_width_range = c(0.25, 2.6),
                                        text_size = 3.2,
                                        title_size = 11) {
  edge_aes <- match.arg(edge_aes)
  
  base_df <- data.frame(
    x = c(0.15, 0.45, 0.75, 1.05),
    y = c(1.00, 1.22, 0.86, 1.00),
    stringsAsFactors = FALSE
  )
  
  if (edge_aes == "alpha") {
    # representative mid-high signal whitening
    z <- 0.75 ^ alpha_gamma
    a <- edge_alpha_range[1] + z * diff(edge_alpha_range)
    col_i <- grDevices::gray(1 - a)
    
    p <- ggplot2::ggplot(base_df) +
      ggplot2::geom_path(
        ggplot2::aes(x = .data$x, y = .data$y),
        color = col_i,
        linewidth = edge_linewidth,
        lineend = "butt",
        linejoin = "round"
      )
  } else {
    width_i <- mean(edge_width_range)
    
    p <- ggplot2::ggplot(base_df) +
      ggplot2::geom_path(
        ggplot2::aes(x = .data$x, y = .data$y),
        color = edge_color,
        alpha = edge_alpha,
        linewidth = width_i,
        lineend = "butt",
        linejoin = "round"
      )
  }
  
  p +
    ggplot2::geom_text(
      data = data.frame(x = 1.35, y = 1.0, lab = if (edge_aes == "alpha") "Grayscale" else "Line width"),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
      hjust = 0,
      color = "black",
      size = text_size
    ) +
    ggplot2::xlim(0, 2.4) +
    ggplot2::ylim(0.55, 1.45) +
    ggplot2::labs(title = legend_title) +
    makeCircuitLegendTheme(title_size = title_size, text_size = text_size)
}

#' Build manual mechanism legend for family circuit plots
buildFamilyMechanismLegend <- function(curve_df,
                                       mechanism_legend_title = "Mechanism",
                                       edge_linewidth = 0.7,
                                       text_size = 3.4,
                                       title_size = 11,
                                       x_label = 1.35,
                                       xlim = c(0, 3.2)) {
  req <- c("identity", "base_color", "identity_linetype")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  
  leg_df <- curve_df
  leg_df$identity_linetype <- as.character(leg_df$identity_linetype)
  
  buildCircuitIdentityLegend(
    df = leg_df,
    label_col = "identity",
    color_col = "base_color",
    linetype_col = "identity_linetype",
    legend_title = mechanism_legend_title,
    edge_linewidth = edge_linewidth,
    text_size = text_size,
    title_size = title_size,
    x_label = x_label,
    xlim = xlim,
    sort_labels = TRUE
  )
}

#' Build manual signal legend for family circuit plots
buildFamilySignalLegend <- function(curve_df,
                                    signal_legend_title = "Signal strength",
                                    signal_legend_n = 5,
                                    text_size = 3.2,
                                    title_size = 11,
                                    x_label = 1.15,
                                    xlim = c(0, 2.3)) {
  req <- c("weight_plot", "intensity_plot")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  
  sig_rng_raw <- range(curve_df$weight_plot, na.rm = TRUE)
  if (!all(is.finite(sig_rng_raw))) sig_rng_raw <- c(0, 1)
  if (diff(sig_rng_raw) == 0) sig_rng_raw <- sig_rng_raw + c(-1e-9, 1e-9)
  
  int_rng <- range(curve_df$intensity_plot, na.rm = TRUE)
  if (!all(is.finite(int_rng))) int_rng <- c(0.12, 1)
  if (diff(int_rng) == 0) int_rng <- int_rng + c(-1e-9, 1e-9)
  
  sig_vals <- seq(sig_rng_raw[1], sig_rng_raw[2], length.out = signal_legend_n)
  sig_int <- scales::rescale(sig_vals, to = int_rng, from = sig_rng_raw)
  sig_cols <- whitenFamilyColors(rep("black", length(sig_int)), sig_int)
  
  sig_df <- data.frame(
    y = rev(seq_along(sig_vals)),
    label = signif(rev(sig_vals), 3),
    col = rev(sig_cols),
    stringsAsFactors = FALSE
  )
  
  ggplot2::ggplot(sig_df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 0.5, y = .data$y),
      width = 0.7,
      height = 0.8,
      fill = sig_df$col,
      color = "black",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_label, y = .data$y, label = .data$label),
      hjust = 0,
      color = "black",
      size = text_size
    ) +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(0.5, max(sig_df$y) + 0.8) +
    ggplot2::labs(title = signal_legend_title) +
    makeCircuitLegendTheme(title_size = title_size, text_size = text_size)
}

#' Assemble a standard circuit legend column
#'
#' Expected order:
#' - signal legend
#' - identity legend
#' - group legend
#' - optional node-size legend
assembleCircuitLegendColumn <- function(signal_legend,
                                        identity_legend,
                                        group_legend,
                                        size_legend = NULL,
                                        show_node_size_legend = FALSE,
                                        heights_no_size = c(1, 1.25, 1),
                                        heights_with_size = c(1, 1.25, 1, 1)) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  
  if (isTRUE(show_node_size_legend) && !is.null(size_legend)) {
    signal_legend / identity_legend / group_legend / size_legend +
      patchwork::plot_layout(heights = heights_with_size)
  } else {
    signal_legend / identity_legend / group_legend +
      patchwork::plot_layout(heights = heights_no_size)
  }
}

#### plotCircuit Helpers ####

#' Resolve circuit grouping aligned to eff node order
resolveCircuitGroups <- function(eff,
                                 group_by,
                                 drop_na = TRUE,
                                 drop_empty = TRUE) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  ids  <- eff$nodes$ids
  meta <- eff$nodes$meta
  
  if (is.null(ids) || length(ids) < 1L) {
    stop("eff$nodes$ids is missing or empty.")
  }
  if (is.null(meta) || !is.data.frame(meta) || nrow(meta) != length(ids)) {
    stop("eff$nodes$meta must be a data.frame with one row per node.")
  }
  
  if (is.character(group_by) && length(group_by) == 1L) {
    if (!(group_by %in% colnames(meta))) {
      stop("group_by = '", group_by, "' was not found in eff$nodes$meta.")
    }
    grp <- meta[[group_by]]
  } else {
    if (length(group_by) != length(ids)) {
      stop(
        "If group_by is supplied as a vector, it must have length equal to ",
        "length(eff$nodes$ids)."
      )
    }
    grp <- group_by
  }
  
  grp_factor <- if (is.factor(grp)) grp else factor(as.character(grp))
  grp_chr    <- as.character(grp_factor)
  
  keep <- rep(TRUE, length(grp_chr))
  if (isTRUE(drop_na)) {
    keep <- keep & !is.na(grp_chr)
  }
  if (isTRUE(drop_empty)) {
    keep <- keep & nzchar(grp_chr)
  }
  
  grp_chr[!keep] <- NA_character_
  
  levels_use <- levels(grp_factor)
  levels_use <- levels_use[levels_use %in% unique(grp_chr[!is.na(grp_chr)])]
  
  list(
    group = grp_chr,
    levels = levels_use,
    keep = keep
  )
}

#' Compute node-level circuit table and circular layout
computeCircuitNodeTable <- function(groups,
                                    levels = NULL,
                                    start_angle = pi / 2,
                                    clockwise = TRUE,
                                    circle_radius = 1,
                                    label_radius = 1.20,
                                    size_nodes_by_abundance = FALSE,
                                    node_radius = 0.06,
                                    node_radius_range = c(0.025, 0.11)) {
  if (length(groups) < 1L) {
    stop("groups must be non-empty.")
  }
  
  grp <- as.character(groups)
  grp <- grp[!is.na(grp)]
  
  if (length(grp) < 1L) {
    stop("No non-missing groups remained.")
  }
  
  if (is.null(levels)) {
    ord <- sort(unique(grp))
  } else {
    ord <- levels[levels %in% unique(grp)]
  }
  
  tab <- table(factor(grp, levels = ord))
  n_cells <- as.integer(tab)
  frac_cells <- n_cells / sum(n_cells)
  
  n <- length(ord)
  step <- 2 * pi / n
  theta <- start_angle + if (isTRUE(clockwise)) -(0:(n - 1)) * step else (0:(n - 1)) * step
  
  x <- circle_radius * cos(theta)
  y <- circle_radius * sin(theta)
  
  if (isTRUE(size_nodes_by_abundance)) {
    # Scale by area, then convert to radius
    frac01 <- if (diff(range(frac_cells)) == 0) {
      rep(1, length(frac_cells))
    } else {
      (frac_cells - min(frac_cells)) / diff(range(frac_cells))
    }
    radius_plot <- node_radius_range[1] + frac01 * diff(node_radius_range)
  } else {
    radius_plot <- rep(node_radius, length(ord))
  }
  
  label_x <- (label_radius + radius_plot) * cos(theta)
  label_y <- (label_radius + radius_plot) * sin(theta)
  
  hjust <- ifelse(label_x >= 0, 0, 1)
  
  data.frame(
    group = ord,
    n_cells = n_cells,
    frac_cells = frac_cells,
    theta = theta,
    x = x,
    y = y,
    radius_plot = radius_plot,
    label_x = label_x,
    label_y = label_y,
    hjust = hjust,
    stringsAsFactors = FALSE
  )
}

#' Aggregate edge-level LR signal to group-level circuit weights
#' Aggregate edge-level LR signal to group-level circuit weights
computeCircuitEdgeTable <- function(eff,
                                    groups,
                                    lr,
                                    transform = identity,
                                    lr_normalization = c("none", "max", "percentile", "zscore"),
                                    lr_aggregation   = c("sum", "mean", "max"),
                                    deduplicate_lr   = TRUE,
                                    warn_many_lr     = TRUE,
                                    edge_group_aggregate = c("mean", "mean_per_sender", "sum", "mean_realized", "mean_per_receiver"),
                                    include_self_loops = TRUE,
                                    include_same_cell_edges = TRUE,
                                    max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  edge_group_aggregate <- match.arg(edge_group_aggregate)
  
  if (is.null(eff$niches$CellToCellSpatial)) {
    stop("eff$niches$CellToCellSpatial is missing.")
  }
  
  c2c <- eff$niches$CellToCellSpatial
  ij  <- c2c$ij
  
  if (is.null(ij) || !is.matrix(ij) || ncol(ij) != 2) {
    stop("CellToCellSpatial$ij must be a matrix with 2 columns.")
  }
  if (nrow(ij) > max_edges) {
    stop("Edge count (", nrow(ij), ") exceeds max_edges (", max_edges, ").")
  }
  if (length(groups) != length(eff$nodes$ids)) {
    stop("groups must have length equal to length(eff$nodes$ids).")
  }
  
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
  
  sender_idx_all <- as.integer(ij[, 1])
  receiver_idx_all <- as.integer(ij[, 2])
  
  if (anyNA(sender_idx_all) || anyNA(receiver_idx_all)) {
    stop("NA detected in CellToCellSpatial$ij.")
  }
  if (any(sender_idx_all < 1L | sender_idx_all > length(groups)) ||
      any(receiver_idx_all < 1L | receiver_idx_all > length(groups))) {
    stop("CellToCellSpatial$ij contains node indices out of bounds.")
  }
  
  from_group_all <- groups[sender_idx_all]
  to_group_all   <- groups[receiver_idx_all]
  edge_signal_all <- lr_out$edge_signal
  
  all_groups <- unique(groups[!is.na(groups)])
  if (is.factor(groups)) {
    all_groups <- levels(groups)[levels(groups) %in% all_groups]
  } else {
    all_groups <- sort(all_groups)
  }
  
  n_groups <- length(all_groups)
  if (n_groups < 1L) {
    stop("No non-missing groups available for circuit aggregation.")
  }
  
  node_group_sizes <- as.numeric(table(factor(groups[!is.na(groups)], levels = all_groups)))
  storage.mode(node_group_sizes) <- "double"
  names(node_group_sizes) <- all_groups
  
  fg_id_all <- match(from_group_all, all_groups)
  tg_id_all <- match(to_group_all, all_groups)
  
  keep_structural <- !is.na(fg_id_all) & !is.na(tg_id_all)
  if (!isTRUE(include_same_cell_edges)) {
    keep_structural <- keep_structural & (sender_idx_all != receiver_idx_all)
  }
  if (!isTRUE(include_self_loops)) {
    keep_structural <- keep_structural & (fg_id_all != tg_id_all)
  }
  
  total_constructed_edges <- sum(keep_structural)
  keep_signal <- keep_structural & is.finite(edge_signal_all)
  
  all_pairs <- expand.grid(
    from_group = all_groups,
    to_group   = all_groups,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  all_pairs$from_id <- match(all_pairs$from_group, all_groups)
  all_pairs$to_id   <- match(all_pairs$to_group, all_groups)
  
  if (!isTRUE(include_self_loops)) {
    all_pairs <- all_pairs[all_pairs$from_id != all_pairs$to_id, , drop = FALSE]
  }
  
  all_pairs$pair_id <- (all_pairs$from_id - 1L) * n_groups + all_pairs$to_id
  all_pairs$n_sender_total <- node_group_sizes[all_pairs$from_id]
  all_pairs$n_receiver_total <- node_group_sizes[all_pairs$to_id]
  
  same_group <- all_pairs$from_id == all_pairs$to_id
  all_pairs$n_possible_edges <- ifelse(
    same_group & !isTRUE(include_same_cell_edges),
    all_pairs$n_sender_total * pmax(all_pairs$n_receiver_total - 1, 0),
    all_pairs$n_sender_total * all_pairs$n_receiver_total
  )
  
  n_pair_bins <- n_groups * n_groups
  
  if (any(keep_structural)) {
    pair_id_struct <- (fg_id_all[keep_structural] - 1L) * n_groups + tg_id_all[keep_structural]
    n_constructed_all <- tabulate(pair_id_struct, nbins = n_pair_bins)
  } else {
    n_constructed_all <- numeric(n_pair_bins)
  }
  all_pairs$n_constructed_edges <- as.numeric(n_constructed_all[all_pairs$pair_id])
  
  if (!any(keep_signal)) {
    edge_tbl <- all_pairs
    edge_tbl$weight_raw <- 0
    edge_tbl$n_sender_participating <- 0
    edge_tbl$n_receiver_participating <- 0
    edge_tbl$n_realized_edges <- 0
    edge_tbl <- edge_tbl[order(edge_tbl$weight_raw), , drop = FALSE]
    edge_tbl$from_id <- NULL
    edge_tbl$to_id <- NULL
    edge_tbl$pair_id <- NULL
    
    return(list(
      edge_table = edge_tbl,
      lr_out = lr_out
    ))
  }
  
  sender_idx <- sender_idx_all[keep_signal]
  receiver_idx <- receiver_idx_all[keep_signal]
  fg_id <- fg_id_all[keep_signal]
  tg_id <- tg_id_all[keep_signal]
  edge_signal <- edge_signal_all[keep_signal]
  
  pair_id_sig <- (fg_id - 1L) * n_groups + tg_id
  
  sum_by_pair_all <- rowsum(edge_signal, group = pair_id_sig, reorder = FALSE)
  sum_by_pair <- numeric(n_pair_bins)
  sum_by_pair[as.integer(rownames(sum_by_pair_all))] <- as.numeric(sum_by_pair_all[, 1])
  
  n_realized_all <- tabulate(pair_id_sig, nbins = n_pair_bins)
  
  sender_pair_key <- paste(pair_id_sig, sender_idx, sep = "\r")
  n_sender_participating_all <- tabulate(
    pair_id_sig[!duplicated(sender_pair_key)],
    nbins = n_pair_bins
  )
  
  receiver_pair_key <- paste(pair_id_sig, receiver_idx, sep = "\r")
  n_receiver_participating_all <- tabulate(
    pair_id_sig[!duplicated(receiver_pair_key)],
    nbins = n_pair_bins
  )
  
  weight_all <- switch(
    edge_group_aggregate,
    sum = sum_by_pair,
    mean = {
      if (!is.finite(total_constructed_edges) || total_constructed_edges <= 0) {
        numeric(n_pair_bins)
      } else {
        sum_by_pair / total_constructed_edges
      }
    },
    mean_realized = {
      out <- numeric(n_pair_bins)
      ok <- n_realized_all > 0
      out[ok] <- sum_by_pair[ok] / n_realized_all[ok]
      out
    },
    mean_per_sender = {
      out <- numeric(n_pair_bins)
      sender_totals_by_pair <- node_group_sizes[((seq_len(n_pair_bins) - 1L) %/% n_groups) + 1L]
      ok <- is.finite(sender_totals_by_pair) & sender_totals_by_pair > 0
      out[ok] <- sum_by_pair[ok] / sender_totals_by_pair[ok]
      out
    },
    mean_per_receiver = {
      out <- numeric(n_pair_bins)
      receiver_totals_by_pair <- node_group_sizes[((seq_len(n_pair_bins) - 1L) %% n_groups) + 1L]
      ok <- is.finite(receiver_totals_by_pair) & receiver_totals_by_pair > 0
      out[ok] <- sum_by_pair[ok] / receiver_totals_by_pair[ok]
      out
    }
  )
  
  idx_pair <- all_pairs$pair_id
  edge_tbl <- all_pairs
  edge_tbl$weight_raw <- as.numeric(weight_all[idx_pair])
  edge_tbl$n_sender_participating <- as.numeric(n_sender_participating_all[idx_pair])
  edge_tbl$n_receiver_participating <- as.numeric(n_receiver_participating_all[idx_pair])
  edge_tbl$n_realized_edges <- as.numeric(n_realized_all[idx_pair])
  
  edge_tbl$weight_raw[!is.finite(edge_tbl$weight_raw)] <- 0
  edge_tbl$n_sender_participating[!is.finite(edge_tbl$n_sender_participating)] <- 0
  edge_tbl$n_receiver_participating[!is.finite(edge_tbl$n_receiver_participating)] <- 0
  edge_tbl$n_realized_edges[!is.finite(edge_tbl$n_realized_edges)] <- 0
  
  edge_tbl <- edge_tbl[order(edge_tbl$weight_raw), , drop = FALSE]
  edge_tbl$from_id <- NULL
  edge_tbl$to_id <- NULL
  edge_tbl$pair_id <- NULL
  
  list(
    edge_table = edge_tbl,
    lr_out = lr_out
  )
}

#' Normalize group-level circuit weights for plotting
normalizeCircuitWeights <- function(edge_tbl,
                                    edge_group_normalize = c("global_max", "none", "row_sum", "col_sum")) {
  edge_group_normalize <- match.arg(edge_group_normalize)
  
  if (!is.data.frame(edge_tbl) || !"weight_raw" %in% colnames(edge_tbl)) {
    stop("edge_tbl must contain a 'weight_raw' column.")
  }
  
  out <- edge_tbl
  w <- out$weight_raw
  
  if (edge_group_normalize == "none") {
    out$weight_plot <- w
    return(out)
  }
  
  if (edge_group_normalize == "global_max") {
    mx <- suppressWarnings(max(w, na.rm = TRUE))
    out$weight_plot <- if (is.finite(mx) && mx > 0) w / mx else 0
    return(out)
  }
  
  if (edge_group_normalize == "row_sum") {
    rs <- tapply(w, out$from_group, sum, na.rm = TRUE)
    out$weight_plot <- ifelse(
      rs[out$from_group] > 0,
      w / rs[out$from_group],
      0
    )
    out$weight_plot[!is.finite(out$weight_plot)] <- 0
    return(out)
  }
  
  if (edge_group_normalize == "col_sum") {
    cs <- tapply(w, out$to_group, sum, na.rm = TRUE)
    out$weight_plot <- ifelse(
      cs[out$to_group] > 0,
      w / cs[out$to_group],
      0
    )
    out$weight_plot[!is.finite(out$weight_plot)] <- 0
    return(out)
  }
  
  out
}

# Shared circuit geometry helpers
.circuit_path_cumlen <- function(xy) {
  dx <- diff(xy[, 1])
  dy <- diff(xy[, 2])
  c(0, cumsum(sqrt(dx^2 + dy^2)))
}

.circuit_interp_point_on_path <- function(xy, s, cumlen = NULL) {
  if (is.null(cumlen)) cumlen <- .circuit_path_cumlen(xy)
  
  total_len <- tail(cumlen, 1)
  if (!is.finite(total_len) || total_len <= 0) {
    return(as.numeric(xy[1, ]))
  }
  
  s <- max(0, min(total_len, s))
  
  j <- which(cumlen >= s)[1]
  if (is.na(j)) return(as.numeric(xy[nrow(xy), ]))
  if (j == 1L) return(as.numeric(xy[1, ]))
  
  i <- j - 1L
  seg_len <- cumlen[j] - cumlen[i]
  
  if (!is.finite(seg_len) || seg_len <= 0) {
    return(as.numeric(xy[j, ]))
  }
  
  u <- (s - cumlen[i]) / seg_len
  as.numeric(xy[i, ] + u * (xy[j, ] - xy[i, ]))
}

.circuit_resample_terminal_path <- function(xy, terminal_length, n = 25L) {
  cumlen <- .circuit_path_cumlen(xy)
  total_len <- tail(cumlen, 1)
  
  if (!is.finite(total_len) || total_len <= 0) {
    return(matrix(
      as.numeric(xy[nrow(xy), ]),
      nrow = 1,
      dimnames = list(NULL, c("x", "y"))
    ))
  }
  
  if (total_len <= terminal_length) {
    s_seq <- seq(0, total_len, length.out = n)
  } else {
    s_seq <- seq(total_len - terminal_length, total_len, length.out = n)
  }
  
  out <- t(vapply(
    s_seq,
    function(s) .circuit_interp_point_on_path(xy, s, cumlen = cumlen),
    FUN.VALUE = c(x = 0, y = 0)
  ))
  
  out <- as.matrix(out)
  colnames(out) <- c("x", "y")
  out
}

.circuit_estimate_tangent <- function(xy) {
  n <- nrow(xy)
  if (n < 2L) {
    return(matrix(c(1, 0), nrow = 1, dimnames = list(NULL, c("x", "y"))))
  }
  
  dx <- numeric(n)
  dy <- numeric(n)
  
  if (n == 2L) {
    dx[] <- xy[2, 1] - xy[1, 1]
    dy[] <- xy[2, 2] - xy[1, 2]
  } else {
    dx[1] <- xy[2, 1] - xy[1, 1]
    dy[1] <- xy[2, 2] - xy[1, 2]
    
    dx[n] <- xy[n, 1] - xy[n - 1, 1]
    dy[n] <- xy[n, 2] - xy[n - 1, 2]
    
    for (i in 2:(n - 1)) {
      dx[i] <- (xy[i + 1, 1] - xy[i - 1, 1]) / 2
      dy[i] <- (xy[i + 1, 2] - xy[i - 1, 2]) / 2
    }
  }
  
  dlen <- sqrt(dx^2 + dy^2)
  dlen[!is.finite(dlen) | dlen == 0] <- 1
  
  cbind(x = dx / dlen, y = dy / dlen)
}

.circuit_trim_path_by_lengths <- function(xy, trim_start, trim_end) {
  cumlen <- .circuit_path_cumlen(xy)
  total_len <- tail(cumlen, 1)
  
  if (!is.finite(total_len) || total_len <= 0) {
    return(xy[1, , drop = FALSE])
  }
  
  s0 <- max(0, trim_start)
  s1 <- min(total_len, total_len - trim_end)
  
  if (!is.finite(s0)) s0 <- 0
  if (!is.finite(s1)) s1 <- total_len
  
  if (s1 <= s0) {
    mid_s <- total_len / 2
    pt <- .circuit_interp_point_on_path(xy, mid_s, cumlen)
    return(matrix(pt, nrow = 1, dimnames = list(NULL, c("x", "y"))))
  }
  
  keep_mid <- cumlen > s0 & cumlen < s1
  
  start_pt <- .circuit_interp_point_on_path(xy, s0, cumlen)
  end_pt   <- .circuit_interp_point_on_path(xy, s1, cumlen)
  
  out <- rbind(
    start_pt,
    xy[keep_mid, , drop = FALSE],
    end_pt
  )
  
  out <- as.matrix(out)
  colnames(out) <- c("x", "y")
  out
}

#' Build circuit edge paths with exact endpoint truncation on the rendered lane
#'
#' Non-self edges:
#' - build one shared outward-bowing quadratic guide arc per unordered pair
#' - offset to sender/receiver-specific lanes
#' - orient each lane from sender -> receiver
#' - truncate the ACTUAL rendered lane by:
#'     sender radius + edge_gap   from the sending end
#'     receiver radius + edge_gap from the receiving end
#' - return the truncated sampled path directly
#'
#' Self edges:
#' - retain a compact loop path outside the node boundary
buildCircuitCurves <- function(node_tbl,
                               edge_tbl,
                               edge_width_plot = 0.03,
                               edge_gap = 0.015,
                               bidirectional_offset = NULL,
                               self_loop_spread = 0.62,
                               self_loop_height = 0.72,
                               curve_strength = 0.13,
                               n_pts = 401,
                               center = c(0, 0)) {
  if (!is.data.frame(node_tbl) ||
      !all(c("group", "x", "y", "theta", "radius_plot") %in% names(node_tbl))) {
    stop("node_tbl must contain group, x, y, theta, and radius_plot.")
  }
  if (!is.data.frame(edge_tbl) ||
      !all(c("from_group", "to_group", "weight_raw", "weight_plot") %in% names(edge_tbl))) {
    stop("edge_tbl must contain from_group, to_group, weight_raw, and weight_plot.")
  }
  
  if (is.null(bidirectional_offset)) {
    med_r <- stats::median(node_tbl$radius_plot, na.rm = TRUE)
    if (!is.finite(med_r) || med_r <= 0) med_r <- 0.06
    bidirectional_offset <- med_r * 0.20
  }
  
  if (!is.numeric(n_pts) || length(n_pts) != 1L || !is.finite(n_pts) || n_pts < 50) {
    stop("n_pts must be a single integer >= 50.")
  }
  n_pts <- as.integer(n_pts)
  
  node_map <- split(node_tbl, node_tbl$group)
  n_nodes_total <- nrow(node_tbl)
  if (!is.finite(n_nodes_total) || n_nodes_total < 2L) {
    stop("node_tbl must contain at least two nodes.")
  }
  
  group_index <- seq_len(n_nodes_total)
  names(group_index) <- node_tbl$group
  
  .quad_bezier <- function(p0, p1, p2, t) {
    omt <- 1 - t
    cbind(
      x = omt^2 * p0[1] + 2 * omt * t * p1[1] + t^2 * p2[1],
      y = omt^2 * p0[2] + 2 * omt * t * p1[2] + t^2 * p2[2]
    )
  }
  
  .quad_bezier_deriv <- function(p0, p1, p2, t) {
    cbind(
      x = 2 * (1 - t) * (p1[1] - p0[1]) + 2 * t * (p2[1] - p1[1]),
      y = 2 * (1 - t) * (p1[2] - p0[2]) + 2 * t * (p2[2] - p1[2])
    )
  }
  
  .cubic_bezier <- function(p0, p1, p2, p3, t) {
    omt <- 1 - t
    cbind(
      x = omt^3 * p0[1] +
        3 * omt^2 * t * p1[1] +
        3 * omt * t^2 * p2[1] +
        t^3 * p3[1],
      y = omt^3 * p0[2] +
        3 * omt^2 * t * p1[2] +
        3 * omt * t^2 * p2[2] +
        t^3 * p3[2]
    )
  }
  
  make_nonself <- function(i, row) {
    fg <- as.character(row$from_group)
    tg <- as.character(row$to_group)
    
    s <- node_map[[fg]]
    t <- node_map[[tg]]
    
    pair_groups <- sort(c(fg, tg))
    a <- node_map[[pair_groups[1]]]
    b <- node_map[[pair_groups[2]]]
    
    p0 <- c(a$x, a$y)
    p2 <- c(b$x, b$y)
    
    gv <- p2 - p0
    glen <- sqrt(sum(gv^2))
    if (!is.finite(glen) || glen == 0) glen <- 1
    gu <- gv / glen
    
    n1 <- c(-gu[2],  gu[1])
    n2 <- -n1
    
    mid <- (p0 + p2) / 2
    outward <- mid - center
    
    gn <- if (sum(n1 * outward) >= sum(n2 * outward)) n2 else n1
    
    ia <- unname(group_index[pair_groups[1]])
    ib <- unname(group_index[pair_groups[2]])
    
    step_dist <- abs(ia - ib)
    circ_dist <- min(step_dist, n_nodes_total - step_dist)
    
    curve_scale <- if (circ_dist <= 1L) {
      -0.40
    } else if (circ_dist == 2L) {
      0.90
    } else if (circ_dist == 3L) {
      0.85
    } else {
      0.80
    }
    
    p1 <- mid + (curve_strength * curve_scale) * gn
    
    tt <- seq(0, 1, length.out = n_pts)
    base_xy <- .quad_bezier(p0, p1, p2, tt)
    base_d  <- .quad_bezier_deriv(p0, p1, p2, tt)
    
    dlen <- sqrt(base_d[, "x"]^2 + base_d[, "y"]^2)
    dlen[!is.finite(dlen) | dlen == 0] <- 1
    
    tx <- base_d[, "x"] / dlen
    ty <- base_d[, "y"] / dlen
    
    nx <- -ty
    ny <-  tx
    
    forward_dir <- (fg == pair_groups[1] && tg == pair_groups[2])
    lane_side <- if (forward_dir) 1 else -1
    
    lane_xy <- cbind(
      x = base_xy[, "x"] + lane_side * bidirectional_offset * nx,
      y = base_xy[, "y"] + lane_side * bidirectional_offset * ny
    )
    
    if (!forward_dir) {
      lane_xy <- lane_xy[nrow(lane_xy):1, , drop = FALSE]
    }
    
    trim_start <- as.numeric(s$radius_plot + edge_gap)
    trim_end   <- as.numeric(t$radius_plot + edge_gap)
    
    lane_trim <- .circuit_trim_path_by_lengths(
      xy = lane_xy,
      trim_start = trim_start,
      trim_end = trim_end
    )
    
    arrow_side <- if (forward_dir) lane_side else -lane_side
    
    data.frame(
      edge_id = i,
      x = lane_trim[, "x"],
      y = lane_trim[, "y"],
      point_id = seq_len(nrow(lane_trim)),
      from_group = fg,
      to_group = tg,
      weight_raw = row$weight_raw,
      weight_plot = row$weight_plot,
      lane_side = lane_side,
      arrow_side = arrow_side,
      is_self = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  make_self <- function(i, row) {
    fg <- as.character(row$from_group)
    s <- node_map[[fg]]
    
    # Local outward radial direction for this node
    u_out <- c(cos(s$theta), sin(s$theta))
    t_hat <- c(-u_out[2], u_out[1])
    
    # Contact radius just outside the node
    r_contact <- s$radius_plot + edge_gap
    
    # Place the two contact points on either side of the outward pole.
    # Their normals are the true node-circle normals at the contact points.
    contact_half_angle <- 0.28 * pi
    
    n_start <- cos(contact_half_angle) * u_out + sin(contact_half_angle) * t_hat
    n_end   <- cos(contact_half_angle) * u_out - sin(contact_half_angle) * t_hat
    
    p0 <- c(s$x, s$y) + r_contact * n_start
    p6 <- c(s$x, s$y) + r_contact * n_end
    
    # Apex of the loop, farther outward than the contact circle
    apex_dist <- r_contact + self_loop_height * s$radius_plot
    apex <- c(s$x, s$y) + apex_dist * u_out
    
    # Control handle lengths at the node contacts:
    # these enforce exact normal-aligned entry/exit
    h_contact <- 0.95 * self_loop_height * s$radius_plot
    
    p1 <- p0 + h_contact * n_start
    p5 <- p6 + h_contact * n_end
    
    # Lateral shoulder controls around the apex.
    # Using a shoulder wider than the contact spread forces the loop
    # to wrap past a semicircle (>180°) rather than form a simple bump.
    shoulder <- 1.35 * self_loop_spread * s$radius_plot
    
    p2 <- apex + shoulder * t_hat
    p4 <- apex - shoulder * t_hat
    
    # Build two cubic segments:
    # p0 -> apex and apex -> p6
    tt1 <- seq(0, 1, length.out = ceiling(n_pts / 2))
    tt2 <- seq(0, 1, length.out = n_pts - length(tt1) + 1L)
    
    seg1 <- .cubic_bezier(p0, p1, p2, apex, tt1)
    seg2 <- .cubic_bezier(apex, p4, p5, p6, tt2)
    
    # Join without duplicating the apex row
    loop_xy <- rbind(
      seg1,
      seg2[-1, , drop = FALSE]
    )
    
    data.frame(
      edge_id = i,
      x = loop_xy[, "x"],
      y = loop_xy[, "y"],
      point_id = seq_len(nrow(loop_xy)),
      from_group = fg,
      to_group = fg,
      weight_raw = row$weight_raw,
      weight_plot = row$weight_plot,
      lane_side = 1,
      arrow_side = -1,
      is_self = TRUE,
      stringsAsFactors = FALSE
    )
  }
  
  if (nrow(edge_tbl) == 0L) {
    return(data.frame())
  }
  
  out <- lapply(seq_len(nrow(edge_tbl)), function(i) {
    row <- edge_tbl[i, , drop = FALSE]
    if (row$from_group == row$to_group) {
      make_self(i, row)
    } else {
      make_nonself(i, row)
    }
  })
  
  do.call(rbind, out)
}

#' Render LR circuit plot from prepared node and edge tables
renderLRCircuitPlot <- function(node_tbl,
                                edge_tbl,
                                plot_title = NULL,
                                node_palette = NULL,
                                node_outline = "black",
                                edge_aes = c("alpha", "width"),
                                edge_color = "black",
                                edge_alpha = 0.9,
                                edge_linewidth = 0.7,
                                edge_alpha_range = c(0.01, 1),
                                alpha_gamma = 2.2,
                                edge_width_range = c(0.25, 2.6),
                                edge_gap = 0.015,
                                bidirectional_offset = NULL,
                                self_loop_spread = 0.62,
                                self_loop_height = 0.72,
                                curve_strength = 0.16,
                                label_size = 4,
                                xlim = c(-1.40, 1.40),
                                ylim = c(-1.40, 1.40),
                                show_node_size_legend = NULL,
                                node_size_legend_title = "Group abundance",
                                node_size_legend_range = c(4, 10),
                                signal_legend_title = "Signal strength",
                                edge_encoding_legend_title = "Edge encoding",
                                group_legend_title = "Group",
                                legend_width_ratio = 1.45) {
  edge_aes <- match.arg(edge_aes)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Package 'ggforce' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }
  
  if (!is.numeric(alpha_gamma) || length(alpha_gamma) != 1L ||
      !is.finite(alpha_gamma) || alpha_gamma <= 0) {
    stop("alpha_gamma must be a single positive number.")
  }
  
  pal <- resolveCircuitNodePalette(node_tbl$group, node_palette = node_palette)
  node_tbl$fill_col <- unname(pal[node_tbl$group])
  
  if (is.null(show_node_size_legend)) {
    show_node_size_legend <- length(unique(round(node_tbl$radius_plot, 10))) > 1L
  }
  
  curve_df <- buildCircuitCurves(
    node_tbl = node_tbl,
    edge_tbl = edge_tbl,
    edge_gap = edge_gap,
    bidirectional_offset = bidirectional_offset,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    curve_strength = curve_strength
  )
  
  p_main <- ggplot2::ggplot()
  
  if (nrow(curve_df) > 0L) {
    wvals <- curve_df$weight_plot[is.finite(curve_df$weight_plot)]
    if (length(wvals) < 1L) {
      stop("No finite edge weights available for plotting.")
    }
    
    scale_limits <- range(wvals, na.rm = TRUE)
    if (!all(is.finite(scale_limits))) {
      stop("Could not determine finite scale limits for edge weights.")
    }
    if (diff(scale_limits) == 0) {
      eps <- max(abs(scale_limits[1]), 1) * 1e-9
      scale_limits <- scale_limits + c(-eps, eps)
    }
    
    clip_weights <- function(x) {
      pmax(scale_limits[1], pmin(scale_limits[2], x))
    }
    
    curve_df$weight_plot_clipped <- clip_weights(curve_df$weight_plot)
    
    edge_param_df <- unique(
      curve_df[, c("edge_id", "weight_plot_clipped"), drop = FALSE]
    )
    edge_param_df <- edge_param_df[order(edge_param_df$weight_plot_clipped,
                                         edge_param_df$edge_id), , drop = FALSE]
    
    if (edge_aes == "width") {
      edge_param_df$head_width_use <- scales::rescale(
        edge_param_df$weight_plot_clipped,
        to = c(0.012, 0.040),
        from = scale_limits
      )
      edge_param_df$head_length_use <- scales::rescale(
        edge_param_df$weight_plot_clipped,
        to = c(0.040, 0.090),
        from = scale_limits
      )
      
      edge_param_df$taper_length_use <- edge_param_df$head_length_use
      edge_param_df$taper_base_halfwidth_use <- scales::rescale(
        edge_param_df$weight_plot_clipped,
        to = c(0.0025, 0.0075),
        from = scale_limits
      )
    } else {
      edge_param_df$head_width_use <- 0.022
      edge_param_df$head_length_use <- 0.060
      edge_param_df$taper_length_use <- 0.060
      edge_param_df$taper_base_halfwidth_use <- 0.0035
    }
    
    m_edge <- match(curve_df$edge_id, edge_param_df$edge_id)
    curve_df$head_width_use <- edge_param_df$head_width_use[m_edge]
    curve_df$head_length_use <- edge_param_df$head_length_use[m_edge]
    curve_df$taper_length_use <- edge_param_df$taper_length_use[m_edge]
    curve_df$taper_base_halfwidth_use <- edge_param_df$taper_base_halfwidth_use[m_edge]
    
    curve_df <- curve_df[order(curve_df$weight_plot_clipped,
                               curve_df$edge_id,
                               curve_df$point_id), , drop = FALSE]
    
    curve_split <- split(curve_df, curve_df$edge_id)
    
    shaft_list <- lapply(curve_split, function(d) {
      trim_i <- unique(d$taper_length_use)
      trim_i <- trim_i[is.finite(trim_i)]
      if (length(trim_i) != 1L) trim_i <- 0.060
      
      out <- trimCircuitCurvesForArrowhead(
        curve_df = d,
        trim_length = trim_i
      )
      
      if (is.null(out) || nrow(out) < 2L) {
        return(NULL)
      }
      
      out
    })
    
    shaft_list <- shaft_list[!vapply(shaft_list, is.null, logical(1))]
    
    if (length(shaft_list) > 0L) {
      curve_df_shaft <- do.call(rbind, shaft_list)
      rownames(curve_df_shaft) <- NULL
      curve_df_shaft$weight_plot_clipped <- clip_weights(curve_df_shaft$weight_plot)
    } else {
      curve_df_shaft <- curve_df[0, , drop = FALSE]
      curve_df_shaft$weight_plot_clipped <- numeric(0)
    }
    
    taper_poly <- buildCircuitTerminalTaperPolygons(curve_df = curve_df)
    arrow_poly <- buildCircuitArrowPolygons(curve_df = curve_df)
    
    if (nrow(taper_poly) > 0L) {
      taper_poly$weight_plot_clipped <- clip_weights(taper_poly$weight_plot)
      taper_poly$poly_group <- interaction(taper_poly$edge_id, taper_poly$part_id, drop = TRUE)
      taper_poly <- taper_poly[order(taper_poly$weight_plot_clipped,
                                     taper_poly$edge_id,
                                     taper_poly$part_id,
                                     taper_poly$point_id), , drop = FALSE]
    }
    
    if (nrow(arrow_poly) > 0L) {
      arrow_poly$weight_plot_clipped <- clip_weights(arrow_poly$weight_plot)
      arrow_poly$poly_group <- interaction(arrow_poly$edge_id, arrow_poly$part_id, drop = TRUE)
      arrow_poly <- arrow_poly[order(arrow_poly$weight_plot_clipped,
                                     arrow_poly$edge_id,
                                     arrow_poly$part_id,
                                     arrow_poly$point_id), , drop = FALSE]
    }
    
    if (edge_aes == "alpha") {
      remap_alpha_value <- function(x, limits, gamma = 2.2) {
        z <- scales::rescale(x, to = c(0, 1), from = limits)
        z[!is.finite(z)] <- 0
        z <- pmax(0, pmin(1, z))
        z ^ gamma
      }
      
      if (nrow(curve_df_shaft) > 0L) {
        curve_df_shaft$weight_plot_alpha <- remap_alpha_value(
          curve_df_shaft$weight_plot_clipped,
          limits = scale_limits,
          gamma = alpha_gamma
        )
      } else {
        curve_df_shaft$weight_plot_alpha <- numeric(0)
      }
      
      if (nrow(taper_poly) > 0L) {
        taper_poly$weight_plot_alpha <- remap_alpha_value(
          taper_poly$weight_plot_clipped,
          limits = scale_limits,
          gamma = alpha_gamma
        )
      }
      
      if (nrow(arrow_poly) > 0L) {
        arrow_poly$weight_plot_alpha <- remap_alpha_value(
          arrow_poly$weight_plot_clipped,
          limits = scale_limits,
          gamma = alpha_gamma
        )
      }
      
      if (nrow(curve_df_shaft) > 0L) {
        edge_order_df <- unique(
          curve_df_shaft[, c("edge_id", "weight_plot_clipped", "weight_plot_alpha"), drop = FALSE]
        )
      } else {
        edge_order_df <- unique(
          curve_df[, c("edge_id", "weight_plot_clipped"), drop = FALSE]
        )
        edge_order_df$weight_plot_alpha <- remap_alpha_value(
          edge_order_df$weight_plot_clipped,
          limits = scale_limits,
          gamma = alpha_gamma
        )
      }
      
      edge_order_df <- edge_order_df[order(edge_order_df$weight_plot_clipped,
                                           edge_order_df$edge_id), , drop = FALSE]
      edge_ids_in_order <- edge_order_df$edge_id
      
      shaft_split <- if (nrow(curve_df_shaft) > 0L) {
        split(curve_df_shaft, curve_df_shaft$edge_id)
      } else {
        list()
      }
      taper_split <- if (nrow(taper_poly) > 0L) split(taper_poly, taper_poly$edge_id) else list()
      arrow_split <- if (nrow(arrow_poly) > 0L) split(arrow_poly, arrow_poly$edge_id) else list()
      
      p_main <- p_main
      
      for (eid in edge_ids_in_order) {
        d_shaft <- shaft_split[[as.character(eid)]]
        if (!is.null(d_shaft) && nrow(d_shaft) > 1L) {
          p_main <- p_main +
            ggplot2::geom_path(
              data = d_shaft,
              ggplot2::aes(
                x = .data$x,
                y = .data$y,
                group = .data$edge_id,
                color = .data$weight_plot_alpha
              ),
              linewidth = edge_linewidth,
              alpha = 1,
              lineend = "butt",
              linejoin = "round",
              inherit.aes = FALSE,
              show.legend = FALSE
            )
        }
        
        d_taper <- taper_split[[as.character(eid)]]
        if (!is.null(d_taper) && nrow(d_taper) > 0L) {
          d_taper$poly_group <- interaction(d_taper$edge_id, d_taper$part_id, drop = TRUE)
          
          p_main <- p_main +
            ggplot2::geom_polygon(
              data = d_taper,
              ggplot2::aes(
                x = .data$x,
                y = .data$y,
                group = .data$poly_group,
                color = .data$weight_plot_alpha,
                fill = after_scale(colour)
              ),
              linewidth = 0,
              inherit.aes = FALSE,
              show.legend = FALSE
            )
        }
        
        d_arrow <- arrow_split[[as.character(eid)]]
        if (!is.null(d_arrow) && nrow(d_arrow) > 0L) {
          d_arrow$poly_group <- interaction(d_arrow$edge_id, d_arrow$part_id, drop = TRUE)
          
          p_main <- p_main +
            ggplot2::geom_polygon(
              data = d_arrow,
              ggplot2::aes(
                x = .data$x,
                y = .data$y,
                group = .data$poly_group,
                color = .data$weight_plot_alpha,
                fill = after_scale(colour)
              ),
              linewidth = 0,
              inherit.aes = FALSE,
              show.legend = FALSE
            )
        }
      }
      
      grey_cols <- alphaRangeToGrey(edge_alpha_range)
      p_main <- p_main +
        ggplot2::scale_color_gradient(
          low = grey_cols$low,
          high = grey_cols$high,
          limits = c(0, 1),
          guide = "none"
        )
    } else {
      if (nrow(curve_df_shaft) > 0L) {
        p_main <- p_main +
          ggplot2::geom_path(
            data = curve_df_shaft,
            ggplot2::aes(
              x = .data$x,
              y = .data$y,
              group = .data$edge_id,
              linewidth = .data$weight_plot_clipped
            ),
            color = edge_color,
            alpha = edge_alpha,
            lineend = "butt",
            linejoin = "round",
            show.legend = FALSE
          )
      }
      
      if (nrow(taper_poly) > 0L) {
        p_main <- p_main +
          ggplot2::geom_polygon(
            data = taper_poly,
            ggplot2::aes(
              x = .data$x,
              y = .data$y,
              group = .data$poly_group
            ),
            fill = edge_color,
            alpha = edge_alpha,
            color = NA,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
      
      if (nrow(arrow_poly) > 0L) {
        p_main <- p_main +
          ggplot2::geom_polygon(
            data = arrow_poly,
            ggplot2::aes(
              x = .data$x,
              y = .data$y,
              group = .data$poly_group
            ),
            fill = edge_color,
            alpha = edge_alpha,
            color = NA,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
      }
      
      p_main <- p_main +
        ggplot2::scale_linewidth_continuous(
          range = edge_width_range,
          limits = scale_limits,
          guide = "none"
        )
    }
  }
  
  p_main <- p_main +
    ggforce::geom_circle(
      data = node_tbl,
      ggplot2::aes(
        x0 = .data$x,
        y0 = .data$y,
        r = .data$radius_plot,
        fill = .data$fill_col
      ),
      color = node_outline,
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(
      data = node_tbl,
      ggplot2::aes(
        x = .data$label_x,
        y = .data$label_y,
        label = .data$group,
        hjust = .data$hjust
      ),
      color = "black",
      size = label_size,
      show.legend = FALSE
    ) +
    ggplot2::coord_equal(xlim = xlim, ylim = ylim, clip = "off") +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(hjust = 0.5, color = "black")
    )
  
  signal_values <- if (nrow(edge_tbl) > 0L) edge_tbl$weight_plot else c(0, 1)
  
  if (edge_aes == "alpha") {
    p_signal_legend <- buildCircuitSignalLegendAlpha(
      values = signal_values,
      legend_title = signal_legend_title,
      alpha_range = edge_alpha_range,
      alpha_gamma = alpha_gamma,
      n = 5
    )
  } else {
    p_signal_legend <- buildCircuitSignalLegendWidth(
      values = signal_values,
      legend_title = signal_legend_title,
      edge_width_range = edge_width_range,
      edge_color = edge_color,
      edge_alpha = edge_alpha,
      n = 5
    )
  }
  
  p_edge_encoding_legend <- buildCircuitEdgeStyleLegend(
    edge_aes = edge_aes,
    legend_title = edge_encoding_legend_title,
    edge_color = edge_color,
    edge_alpha = edge_alpha,
    edge_linewidth = edge_linewidth,
    edge_alpha_range = edge_alpha_range,
    alpha_gamma = alpha_gamma,
    edge_width_range = edge_width_range
  )
  
  p_group_legend <- buildCircuitGroupLegend(
    node_tbl = node_tbl,
    node_palette = node_palette,
    node_outline = node_outline,
    legend_title = group_legend_title
  )
  
  if (isTRUE(show_node_size_legend)) {
    p_size_legend <- buildCircuitNodeSizeLegend(
      node_tbl = node_tbl,
      node_size_legend_title = node_size_legend_title,
      node_size_legend_range = node_size_legend_range
    )
    
    legend_col <- p_signal_legend / p_edge_encoding_legend / p_group_legend / p_size_legend +
      patchwork::plot_layout(heights = c(1, 0.85, 1, 1))
  } else {
    legend_col <- p_signal_legend / p_edge_encoding_legend / p_group_legend +
      patchwork::plot_layout(heights = c(1, 0.85, 1))
  }
  
  out <- p_main | legend_col
  out + patchwork::plot_layout(widths = c(3.8, legend_width_ratio))
}

#' Resolve node palette for circuit groups
resolveCircuitNodePalette <- function(groups, node_palette = NULL) {
  groups <- normalizeIdentityLabels(groups)
  
  if (length(groups) < 1L) {
    stop("No groups available for node palette resolution.")
  }
  
  if (is.null(node_palette)) {
    return(defaultIdentityColors(groups, palette = "Set 2"))
  }
  
  resolveIdentityAestheticVector(
    labels = groups,
    values = node_palette,
    arg_name = "node_palette"
  )
}

#' Convert an alpha range on a white background into equivalent grayscale colors
alphaRangeToGrey <- function(alpha_range = c(0.08, 1)) {
  if (length(alpha_range) != 2 || any(!is.finite(alpha_range)) ||
      alpha_range[1] < 0 || alpha_range[2] <= 0 ||
      alpha_range[1] >= alpha_range[2] || alpha_range[2] > 1) {
    stop("alpha_range must be something like c(0.08, 1) with values in [0, 1].")
  }
  
  list(
    low  = grDevices::gray(1 - alpha_range[1]),
    high = grDevices::gray(1 - alpha_range[2])
  )
}

buildCircuitArrowPolygons <- function(curve_df,
                                      head_length = 0.060,
                                      head_width = 0.022,
                                      n_resample = 25) {
  req <- c("edge_id", "x", "y", "point_id", "weight_plot", "arrow_side")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  
  if (!is.numeric(head_length) || length(head_length) != 1L ||
      !is.finite(head_length) || head_length <= 0) {
    stop("head_length must be a single positive number.")
  }
  if (!is.numeric(head_width) || length(head_width) != 1L ||
      !is.finite(head_width) || head_width <= 0) {
    stop("head_width must be a single positive number.")
  }
  if (!is.numeric(n_resample) || length(n_resample) != 1L ||
      !is.finite(n_resample) || n_resample < 5) {
    stop("n_resample must be a single integer >= 5.")
  }
  n_resample <- as.integer(n_resample)
  
  edge_split <- split(curve_df, curve_df$edge_id)
  
  out <- lapply(edge_split, function(d) {
    d <- d[order(d$point_id), , drop = FALSE]
    if (nrow(d) < 2L) return(NULL)
    
    head_length_i <- if ("head_length_use" %in% names(d)) unique(d$head_length_use) else head_length
    head_length_i <- head_length_i[is.finite(head_length_i)]
    if (length(head_length_i) != 1L) head_length_i <- head_length
    
    head_width_i <- if ("head_width_use" %in% names(d)) unique(d$head_width_use) else head_width
    head_width_i <- head_width_i[is.finite(head_width_i)]
    if (length(head_width_i) != 1L) head_width_i <- head_width
    
    xy_full <- as.matrix(d[, c("x", "y"), drop = FALSE])
    storage.mode(xy_full) <- "double"
    
    cumlen_full <- .circuit_path_cumlen(xy_full)
    total_len <- tail(cumlen_full, 1)
    if (!is.finite(total_len) || total_len <= 0) return(NULL)
    
    xy_term <- .circuit_resample_terminal_path(
      xy = xy_full,
      terminal_length = head_length_i,
      n = n_resample
    )
    
    if (nrow(xy_term) < 2L) return(NULL)
    
    tang <- .circuit_estimate_tangent(xy_term)
    norm_left <- cbind(x = -tang[, "y"], y = tang[, "x"])
    
    side <- unique(d$arrow_side)
    side <- side[is.finite(side)]
    if (length(side) != 1L) side <- 1
    side <- sign(side)
    if (side == 0) side <- 1
    
    t01 <- seq(0, 1, length.out = nrow(xy_term))
    width_profile <- head_width_i * (1 - t01)^1.15
    
    outer_xy <- cbind(
      x = xy_term[, "x"] + side * width_profile * norm_left[, "x"],
      y = xy_term[, "y"] + side * width_profile * norm_left[, "y"]
    )
    
    poly_xy <- rbind(
      xy_term[nrow(xy_term):1, , drop = FALSE],
      outer_xy
    )
    
    data.frame(
      edge_id = d$edge_id[1],
      part_id = 1L,
      x = poly_xy[, "x"],
      y = poly_xy[, "y"],
      point_id = seq_len(nrow(poly_xy)),
      weight_plot = d$weight_plot[nrow(d)],
      stringsAsFactors = FALSE
    )
  })
  
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) return(data.frame())
  
  poly_df <- do.call(rbind, out)
  rownames(poly_df) <- NULL
  poly_df
}

trimCircuitCurvesForArrowhead <- function(curve_df,
                                          trim_length = 0.060) {
  req <- c("edge_id", "x", "y", "point_id")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  
  if (!is.numeric(trim_length) || length(trim_length) != 1L ||
      !is.finite(trim_length) || trim_length < 0) {
    stop("trim_length must be a single non-negative number.")
  }
  
  if (trim_length == 0) {
    out <- curve_df[order(curve_df$edge_id, curve_df$point_id), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  
  edge_split <- split(curve_df, curve_df$edge_id)
  
  out <- lapply(edge_split, function(d) {
    d <- d[order(d$point_id), , drop = FALSE]
    if (nrow(d) < 2L) return(NULL)
    
    xy <- as.matrix(d[, c("x", "y"), drop = FALSE])
    storage.mode(xy) <- "double"
    
    cumlen <- .circuit_path_cumlen(xy)
    total_len <- tail(cumlen, 1)
    
    if (!is.finite(total_len) || total_len <= 0) return(NULL)
    
    keep_to <- total_len - trim_length
    
    # If the path is too short, keep only a minimal stub at the start.
    if (!is.finite(keep_to) || keep_to <= 0) {
      xy_trim <- matrix(
        as.numeric(xy[1, ]),
        nrow = 1,
        dimnames = list(NULL, c("x", "y"))
      )
    } else {
      keep_mid <- cumlen < keep_to
      end_pt <- .circuit_interp_point_on_path(xy, keep_to, cumlen = cumlen)
      
      xy_trim <- rbind(
        xy[keep_mid, , drop = FALSE],
        end_pt
      )
      
      xy_trim <- as.matrix(xy_trim)
      colnames(xy_trim) <- c("x", "y")
      
      # Remove exact duplicate consecutive rows if interpolation lands on an existing point
      if (nrow(xy_trim) > 1L) {
        dup <- c(FALSE, rowSums(abs(diff(xy_trim)) < .Machine$double.eps^0.5) == 2)
        xy_trim <- xy_trim[!dup, , drop = FALSE]
      }
    }
    
    out_d <- d[rep(1, nrow(xy_trim)), , drop = FALSE]
    out_d$x <- xy_trim[, "x"]
    out_d$y <- xy_trim[, "y"]
    out_d$point_id <- seq_len(nrow(xy_trim))
    out_d
  })
  
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) return(curve_df[0, , drop = FALSE])
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

buildCircuitTerminalTaperPolygons <- function(curve_df,
                                              taper_length = 0.060,
                                              taper_base_halfwidth = 0.010,
                                              n_resample = 41) {
  req <- c("edge_id", "x", "y", "point_id", "weight_plot")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  
  if (!is.numeric(taper_length) || length(taper_length) != 1L ||
      !is.finite(taper_length) || taper_length <= 0) {
    stop("taper_length must be a single positive number.")
  }
  if (!is.numeric(taper_base_halfwidth) || length(taper_base_halfwidth) != 1L ||
      !is.finite(taper_base_halfwidth) || taper_base_halfwidth <= 0) {
    stop("taper_base_halfwidth must be a single positive number.")
  }
  if (!is.numeric(n_resample) || length(n_resample) != 1L ||
      !is.finite(n_resample) || n_resample < 7) {
    stop("n_resample must be a single integer >= 7.")
  }
  n_resample <- as.integer(n_resample)
  
  edge_split <- split(curve_df, curve_df$edge_id)
  
  out <- lapply(edge_split, function(d) {
    d <- d[order(d$point_id), , drop = FALSE]
    if (nrow(d) < 2L) return(NULL)
    
    taper_length_i <- if ("taper_length_use" %in% names(d)) unique(d$taper_length_use) else taper_length
    taper_length_i <- taper_length_i[is.finite(taper_length_i)]
    if (length(taper_length_i) != 1L) taper_length_i <- taper_length
    
    taper_halfwidth_i <- if ("taper_base_halfwidth_use" %in% names(d)) unique(d$taper_base_halfwidth_use) else taper_base_halfwidth
    taper_halfwidth_i <- taper_halfwidth_i[is.finite(taper_halfwidth_i)]
    if (length(taper_halfwidth_i) != 1L) taper_halfwidth_i <- taper_base_halfwidth
    
    xy_full <- as.matrix(d[, c("x", "y"), drop = FALSE])
    storage.mode(xy_full) <- "double"
    
    cumlen_full <- .circuit_path_cumlen(xy_full)
    total_len <- tail(cumlen_full, 1)
    if (!is.finite(total_len) || total_len <= 0) return(NULL)
    
    xy_term <- .circuit_resample_terminal_path(
      xy = xy_full,
      terminal_length = taper_length_i,
      n = n_resample
    )
    if (nrow(xy_term) < 2L) return(NULL)
    
    tang <- .circuit_estimate_tangent(xy_term)
    norm_left <- cbind(x = -tang[, "y"], y = tang[, "x"])
    
    t01 <- seq(0, 1, length.out = nrow(xy_term))
    
    # Full ribbon taper: wide at the start, narrows smoothly to 0 at the tip.
    halfwidth_profile <- taper_halfwidth_i * (1 - t01)^1.10
    
    left_xy <- cbind(
      x = xy_term[, "x"] + halfwidth_profile * norm_left[, "x"],
      y = xy_term[, "y"] + halfwidth_profile * norm_left[, "y"]
    )
    
    right_xy <- cbind(
      x = xy_term[, "x"] - halfwidth_profile * norm_left[, "x"],
      y = xy_term[, "y"] - halfwidth_profile * norm_left[, "y"]
    )
    
    poly_xy <- rbind(
      left_xy,
      right_xy[nrow(right_xy):1, , drop = FALSE]
    )
    
    data.frame(
      edge_id = d$edge_id[1],
      part_id = 1L,
      x = poly_xy[, "x"],
      y = poly_xy[, "y"],
      point_id = seq_len(nrow(poly_xy)),
      weight_plot = d$weight_plot[nrow(d)],
      stringsAsFactors = FALSE
    )
  })
  
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) return(data.frame())
  
  poly_df <- do.call(rbind, out)
  rownames(poly_df) <- NULL
  poly_df
}

#### plotCircuit Functions ####

# Compute group-level LR circuit data
computeLRCircuitData <- function(eff,
                                 group_by,
                                 lr,
                                 transform = identity,
                                 lr_normalization = c("none", "max", "percentile", "zscore"),
                                 lr_aggregation   = c("sum", "mean", "max"),
                                 deduplicate_lr   = TRUE,
                                 warn_many_lr     = TRUE,
                                 edge_group_aggregate = c("sum", "mean", "mean_realized", "mean_per_sender", "mean_per_receiver"),
                                 edge_group_normalize = c("global_max", "none", "row_sum", "col_sum"),
                                 include_self_loops = TRUE,
                                 include_same_cell_edges = TRUE,
                                 drop_na_groups = TRUE,
                                 drop_empty_groups = TRUE,
                                 start_angle = pi / 2,
                                 clockwise = TRUE,
                                 size_nodes_by_abundance = FALSE,
                                 max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  edge_group_aggregate <- match.arg(edge_group_aggregate)
  edge_group_normalize <- match.arg(edge_group_normalize)
  
  grp_out <- resolveCircuitGroups(
    eff = eff,
    group_by = group_by,
    drop_na = drop_na_groups,
    drop_empty = drop_empty_groups
  )
  
  node_tbl <- computeCircuitNodeTable(
    groups = grp_out$group,
    levels = grp_out$levels,
    start_angle = start_angle,
    clockwise = clockwise,
    size_nodes_by_abundance = size_nodes_by_abundance
  )
  
  edge_out <- computeCircuitEdgeTable(
    eff = eff,
    groups = grp_out$group,
    lr = lr,
    transform = transform,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    warn_many_lr = warn_many_lr,
    edge_group_aggregate = edge_group_aggregate,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    max_edges = max_edges
  )
  
  edge_tbl <- normalizeCircuitWeights(
    edge_tbl = edge_out$edge_table,
    edge_group_normalize = edge_group_normalize
  )
  
  list(
    nodes = node_tbl,
    edges = edge_tbl,
    lr_out = edge_out$lr_out,
    params = list(
      group_by = group_by,
      lr = lr,
      lr_normalization = lr_normalization,
      lr_aggregation = lr_aggregation,
      deduplicate_lr = deduplicate_lr,
      warn_many_lr = warn_many_lr,
      edge_group_aggregate = edge_group_aggregate,
      edge_group_normalize = edge_group_normalize,
      include_self_loops = include_self_loops,
      include_same_cell_edges = include_same_cell_edges,
      drop_na_groups = drop_na_groups,
      drop_empty_groups = drop_empty_groups,
      start_angle = start_angle,
      clockwise = clockwise,
      max_edges = max_edges
    )
  )
}

#' Plot group-level LR circuit
plotLRCircuit <- function(eff,
                          group_by,
                          lr,
                          transform = identity,
                          # LR combination controls
                          lr_normalization = c("none", "max", "percentile", "zscore"),
                          lr_aggregation   = c("sum", "mean", "max"),
                          deduplicate_lr   = TRUE,
                          warn_many_lr     = TRUE,
                          # group-edge controls
                          edge_group_aggregate = c("sum", "mean", "mean_realized", "mean_per_sender", "mean_per_receiver"),
                          edge_group_normalize = c("global_max", "none", "row_sum", "col_sum"),
                          include_self_loops = TRUE,
                          include_same_cell_edges = TRUE,
                          # grouping controls
                          drop_na_groups = TRUE,
                          drop_empty_groups = TRUE,
                          # node / edge visuals
                          size_nodes_by_abundance = FALSE,
                          node_palette = NULL,
                          node_outline = "black",
                          edge_aes = c("alpha", "width"),
                          edge_color = "black",
                          edge_gap = 0.015,
                          bidirectional_offset = NULL,
                          self_loop_spread = 0.62,
                          self_loop_height = 0.72,
                          edge_alpha = 0.9,
                          edge_linewidth = 0.7,
                          edge_alpha_range = c(0.01, 1),
                          alpha_gamma = 2.2,
                          edge_width_range = c(0.25, 2.6),
                          min_edge = 0,
                          # layout / geometry
                          start_angle = pi / 2,
                          clockwise = TRUE,
                          curve_strength = 0.18,
                          # labels / title
                          label_size = 4,
                          plot_title = NULL,
                          title_suffix = NULL,
                          # legend controls
                          signal_legend_title = "Signal strength",
                          edge_encoding_legend_title = "Edge encoding",
                          group_legend_title = "Group",
                          show_node_size_legend = NULL,
                          node_size_legend_title = "Group abundance",
                          node_size_legend_range = c(4, 10),
                          legend_width_ratio = 1.45,
                          # safety
                          max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  lr_normalization <- match.arg(lr_normalization)
  lr_aggregation   <- match.arg(lr_aggregation)
  edge_group_aggregate <- match.arg(edge_group_aggregate)
  edge_group_normalize <- match.arg(edge_group_normalize)
  edge_aes <- match.arg(edge_aes)
  
  circ <- computeLRCircuitData(
    eff = eff,
    group_by = group_by,
    lr = lr,
    transform = transform,
    lr_normalization = lr_normalization,
    lr_aggregation = lr_aggregation,
    deduplicate_lr = deduplicate_lr,
    warn_many_lr = warn_many_lr,
    edge_group_aggregate = edge_group_aggregate,
    edge_group_normalize = edge_group_normalize,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    drop_na_groups = drop_na_groups,
    drop_empty_groups = drop_empty_groups,
    size_nodes_by_abundance = size_nodes_by_abundance,
    start_angle = start_angle,
    clockwise = clockwise,
    max_edges = max_edges
  )
  
  edge_tbl <- circ$edges
  
  if (!is.numeric(min_edge) || length(min_edge) != 1L || !is.finite(min_edge)) {
    stop("min_edge must be a single finite numeric value.")
  }
  edge_tbl <- edge_tbl[edge_tbl$weight_plot > min_edge, , drop = FALSE]
  
  title_use <- makeLRPlotTitle(
    selected_mechanisms = circ$lr_out$selected_mechanisms,
    plot_title = plot_title,
    suffix = title_suffix
  )
  
  renderLRCircuitPlot(
    node_tbl = circ$nodes,
    edge_tbl = edge_tbl,
    plot_title = title_use,
    node_palette = node_palette,
    node_outline = node_outline,
    edge_aes = edge_aes,
    edge_color = edge_color,
    edge_alpha = edge_alpha,
    edge_linewidth = edge_linewidth,
    edge_alpha_range = edge_alpha_range,
    alpha_gamma = alpha_gamma,
    edge_width_range = edge_width_range,
    edge_gap = edge_gap,
    bidirectional_offset = bidirectional_offset,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    curve_strength = curve_strength,
    label_size = label_size,
    show_node_size_legend = show_node_size_legend,
    node_size_legend_title = node_size_legend_title,
    node_size_legend_range = node_size_legend_range,
    signal_legend_title = signal_legend_title,
    edge_encoding_legend_title = edge_encoding_legend_title,
    group_legend_title = group_legend_title,
    legend_width_ratio = legend_width_ratio
  )
}

#' Multi-plot wrapper for LR circuit plots
plotLRCircuitMulti <- function(eff,
                               group_by,
                               lr,
                               plot_titles = NULL,
                               ncol = NULL,
                               nrow = NULL,
                               guides = "collect",
                               ...) {
  .plotLRMulti(
    plot_fun = function(eff, lr, plot_title, ...) {
      plotLRCircuit(
        eff = eff,
        group_by = group_by,
        lr = lr,
        plot_title = plot_title,
        ...
      )
    },
    eff = eff,
    lr = lr,
    plot_titles = plot_titles,
    ncol = ncol,
    nrow = nrow,
    guides = guides,
    ...
  )
}
#### plotFamilyCircuit Helpers ####

#' Resolve and validate LR family mechanisms
resolveFamilyMechanisms <- function(eff,
                                    lr,
                                    deduplicate_lr = TRUE) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  if (missing(lr) || is.null(lr) || length(lr) < 1L) {
    stop("'lr' must contain at least one LR mechanism.")
  }
  
  if (is.null(eff$niches$CellToCellSpatial)) {
    stop("eff$niches$CellToCellSpatial is missing.")
  }
  
  mech_names <- eff$niches$CellToCellSpatial$mechanisms
  if (is.null(mech_names) || !is.character(mech_names) || length(mech_names) < 1L) {
    stop("CellToCellSpatial$mechanisms must be a non-empty character vector.")
  }
  
  idx <- resolveLRindices(
    lr = lr,
    mech_names = mech_names,
    deduplicate = deduplicate_lr
  )
  
  selected_mechanisms <- mech_names[idx]
  
  if (length(selected_mechanisms) < 1L) {
    stop("No LR mechanisms were resolved.")
  }
  
  list(
    selected_idx = idx,
    selected_mechanisms = selected_mechanisms
  )
}

#' Compute edge x mechanism matrix for a family, with per-mechanism max normalization
computeFamilyEdgeSignalMatrix <- function(eff,
                                          lr,
                                          deduplicate_lr = TRUE,
                                          max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  if (is.null(eff$niches$CellToCellSpatial)) {
    stop("eff$niches$CellToCellSpatial is missing.")
  }
  
  c2c <- eff$niches$CellToCellSpatial
  
  if (is.null(c2c$ij) || !is.matrix(c2c$ij) || ncol(c2c$ij) != 2) {
    stop("CellToCellSpatial$ij must be a matrix with 2 columns.")
  }
  if (is.null(c2c$w) || !inherits(c2c$w, c("dgCMatrix", "Matrix", "matrix"))) {
    stop("CellToCellSpatial$w must be a matrix-like object.")
  }
  if (nrow(c2c$ij) != nrow(c2c$w)) {
    stop("Row mismatch: nrow(CellToCellSpatial$ij) must equal nrow(CellToCellSpatial$w).")
  }
  if (nrow(c2c$ij) > max_edges) {
    stop("Edge count (", nrow(c2c$ij), ") exceeds max_edges (", max_edges, ").")
  }
  
  fam <- resolveFamilyMechanisms(
    eff = eff,
    lr = lr,
    deduplicate_lr = deduplicate_lr
  )
  
  w_sub <- as.matrix(c2c$w[, fam$selected_idx, drop = FALSE])
  storage.mode(w_sub) <- "double"
  colnames(w_sub) <- fam$selected_mechanisms
  
  # Fixed behavior for plotFamilyCircuit:
  # each mechanism is unity-normalized independently
  w_norm <- normalizeLRcolumns(w_sub, method = "max")
  w_norm[!is.finite(w_norm)] <- 0
  
  list(
    ij = c2c$ij,
    w_mech = w_norm,
    selected_idx = fam$selected_idx,
    selected_mechanisms = fam$selected_mechanisms
  )
}

#' Aggregate family mechanism signal to group-level directional circuit edges
computeFamilyCircuitEdgeTable <- function(eff,
                                          groups,
                                          lr,
                                          include_self_loops = TRUE,
                                          include_same_cell_edges = TRUE,
                                          deduplicate_lr = TRUE,
                                          max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  if (length(groups) != length(eff$nodes$ids)) {
    stop("groups must have length equal to length(eff$nodes$ids).")
  }
  
  sig_out <- computeFamilyEdgeSignalMatrix(
    eff = eff,
    lr = lr,
    deduplicate_lr = deduplicate_lr,
    max_edges = max_edges
  )
  
  ij <- sig_out$ij
  w_mech <- sig_out$w_mech
  mechs <- sig_out$selected_mechanisms
  
  sender_idx_all <- as.integer(ij[, 1])
  receiver_idx_all <- as.integer(ij[, 2])
  
  if (anyNA(sender_idx_all) || anyNA(receiver_idx_all)) {
    stop("NA detected in CellToCellSpatial$ij.")
  }
  if (any(sender_idx_all < 1L | sender_idx_all > length(groups)) ||
      any(receiver_idx_all < 1L | receiver_idx_all > length(groups))) {
    stop("CellToCellSpatial$ij contains node indices out of bounds.")
  }
  
  from_group_all <- groups[sender_idx_all]
  to_group_all   <- groups[receiver_idx_all]
  
  all_groups <- unique(groups[!is.na(groups)])
  if (is.factor(groups)) {
    all_groups <- levels(groups)[levels(groups) %in% all_groups]
  } else {
    all_groups <- sort(all_groups)
  }
  
  n_groups <- length(all_groups)
  if (n_groups < 1L) {
    stop("No non-missing groups available for circuit aggregation.")
  }
  
  fg_id_all <- match(from_group_all, all_groups)
  tg_id_all <- match(to_group_all, all_groups)
  
  keep_structural <- !is.na(fg_id_all) & !is.na(tg_id_all)
  if (!isTRUE(include_same_cell_edges)) {
    keep_structural <- keep_structural & (sender_idx_all != receiver_idx_all)
  }
  if (!isTRUE(include_self_loops)) {
    keep_structural <- keep_structural & (fg_id_all != tg_id_all)
  }
  
  all_pairs <- expand.grid(
    from_group = all_groups,
    to_group   = all_groups,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  all_pairs$from_id <- match(all_pairs$from_group, all_groups)
  all_pairs$to_id   <- match(all_pairs$to_group, all_groups)
  
  if (!isTRUE(include_self_loops)) {
    all_pairs <- all_pairs[all_pairs$from_id != all_pairs$to_id, , drop = FALSE]
  }
  
  all_pairs$pair_id <- (all_pairs$from_id - 1L) * n_groups + all_pairs$to_id
  n_pair_bins <- n_groups * n_groups
  
  pair_id_struct <- if (any(keep_structural)) {
    (fg_id_all[keep_structural] - 1L) * n_groups + tg_id_all[keep_structural]
  } else {
    integer(0)
  }
  
  n_constructed_all <- if (length(pair_id_struct) > 0L) {
    tabulate(pair_id_struct, nbins = n_pair_bins)
  } else {
    numeric(n_pair_bins)
  }
  
  all_pairs$n_constructed_edges <- as.numeric(n_constructed_all[all_pairs$pair_id])
  
  if (!any(keep_structural)) {
    out <- all_pairs[rep(seq_len(nrow(all_pairs)), each = length(mechs)), , drop = FALSE]
    out$mechanism <- rep(mechs, times = nrow(all_pairs))
    out$weight_raw <- 0
    out <- out[order(out$weight_raw, out$pair_id, out$mechanism), , drop = FALSE]
    rownames(out) <- NULL
    
    return(list(
      edge_table = out[, c("from_group", "to_group", "mechanism", "weight_raw", "n_constructed_edges"), drop = FALSE],
      selected_mechanisms = mechs
    ))
  }
  
  fg_id <- fg_id_all[keep_structural]
  tg_id <- tg_id_all[keep_structural]
  pair_id <- (fg_id - 1L) * n_groups + tg_id
  
  w_keep <- w_mech[keep_structural, , drop = FALSE]
  
  sum_by_pair_and_mech <- lapply(seq_len(ncol(w_keep)), function(j) {
    tmp <- rowsum(w_keep[, j], group = pair_id, reorder = FALSE)
    out_j <- numeric(n_pair_bins)
    out_j[as.integer(rownames(tmp))] <- as.numeric(tmp[, 1])
    out_j
  })
  
  sum_mat <- do.call(cbind, sum_by_pair_and_mech)
  colnames(sum_mat) <- mechs
  
  out <- all_pairs[rep(seq_len(nrow(all_pairs)), each = length(mechs)), , drop = FALSE]
  out$mechanism <- rep(mechs, times = nrow(all_pairs))
  out$weight_raw <- as.numeric(sum_mat[cbind(out$pair_id, match(out$mechanism, mechs))])
  out$weight_raw[!is.finite(out$weight_raw)] <- 0
  
  out <- out[order(out$weight_raw, out$pair_id, out$mechanism), , drop = FALSE]
  rownames(out) <- NULL
  
  list(
    edge_table = out[, c("from_group", "to_group", "mechanism", "weight_raw", "n_constructed_edges"), drop = FALSE],
    selected_mechanisms = mechs
  )
}

#' Filter and normalize mechanism-specific group edges for family plotting
#'
#' Updated behavior for plotFamilyCircuit only:
#' - threshold is applied independently within each mechanism
#' - retained edges are re-normalized within that mechanism to unity
#' - signal visibility is encoded linearly in alpha
prepareFamilyCircuitEdgeTable <- function(edge_tbl,
                                          min_mechanism_edge = 0,
                                          mechanism_relative_threshold = 0.25,
                                          visible_floor = 0.12,
                                          visible_ceiling = 1) {
  req <- c("from_group", "to_group", "mechanism", "weight_raw")
  if (!is.data.frame(edge_tbl) || !all(req %in% names(edge_tbl))) {
    stop("edge_tbl must contain: ", paste(req, collapse = ", "))
  }
  
  if (!is.numeric(min_mechanism_edge) || length(min_mechanism_edge) != 1L || !is.finite(min_mechanism_edge)) {
    stop("min_mechanism_edge must be a single finite numeric value.")
  }
  if (!is.numeric(mechanism_relative_threshold) || length(mechanism_relative_threshold) != 1L ||
      !is.finite(mechanism_relative_threshold) ||
      mechanism_relative_threshold < 0 || mechanism_relative_threshold > 1) {
    stop("mechanism_relative_threshold must be a single number in [0, 1].")
  }
  if (!is.numeric(visible_floor) || length(visible_floor) != 1L || !is.finite(visible_floor) ||
      visible_floor < 0 || visible_floor > 1) {
    stop("visible_floor must be a single number in [0, 1].")
  }
  if (!is.numeric(visible_ceiling) || length(visible_ceiling) != 1L || !is.finite(visible_ceiling) ||
      visible_ceiling <= 0 || visible_ceiling > 1 || visible_ceiling < visible_floor) {
    stop("visible_ceiling must be in (0, 1] and >= visible_floor.")
  }

  if (nrow(edge_tbl) == 0L) {
    edge_tbl$weight_plot <- numeric(0)
    edge_tbl$intensity_plot <- numeric(0)
    edge_tbl$alpha_plot <- numeric(0)
    return(edge_tbl)
  }
  
  mech_split <- split(edge_tbl, as.character(edge_tbl$mechanism))
  
  out <- lapply(mech_split, function(d) {
    w <- d$weight_raw
    w[!is.finite(w)] <- NA_real_
    
    mech_max <- suppressWarnings(max(w, na.rm = TRUE))
    if (!is.finite(mech_max) || mech_max <= 0) {
      return(NULL)
    }
    
    abs_threshold <- max(min_mechanism_edge, mechanism_relative_threshold * mech_max)
    
    keep <- is.finite(w) & (w >= abs_threshold)
    d <- d[keep, , drop = FALSE]
    if (nrow(d) == 0L) {
      return(NULL)
    }
    
    mech_max_retained <- suppressWarnings(max(d$weight_raw, na.rm = TRUE))
    if (!is.finite(mech_max_retained) || mech_max_retained <= 0) {
      d$weight_plot <- 0
    } else {
      d$weight_plot <- d$weight_raw / mech_max_retained
    }
    d$weight_plot[!is.finite(d$weight_plot)] <- 0
    
    z <- d$weight_plot
    z[!is.finite(z)] <- 0
    z <- pmax(0, pmin(1, z))
    
    # linear alpha for retained edges
    d$alpha_plot <- visible_floor + z * (visible_ceiling - visible_floor)
    d$alpha_plot[!is.finite(d$alpha_plot)] <- visible_floor
    
    # retain this column name for downstream compatibility
    d$intensity_plot <- d$alpha_plot
    
    d
  })
  
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) {
    out <- edge_tbl[0, , drop = FALSE]
    out$weight_plot <- numeric(0)
    out$intensity_plot <- numeric(0)
    out$alpha_plot <- numeric(0)
    return(out)
  }
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Rank retained mechanisms within each directional group edge
rankFamilyMechanismsWithinEdges <- function(edge_tbl) {
  req <- c("from_group", "to_group", "mechanism", "weight_raw", "weight_plot", "intensity_plot")
  if (!is.data.frame(edge_tbl) || !all(req %in% names(edge_tbl))) {
    stop("edge_tbl must contain: ", paste(req, collapse = ", "))
  }
  
  if (!"alpha_plot" %in% names(edge_tbl)) {
    edge_tbl$alpha_plot <- edge_tbl$intensity_plot
  }
  
  if (nrow(edge_tbl) == 0L) {
    edge_tbl$rank_within_edge <- integer(0)
    edge_tbl$offset_index <- integer(0)
    edge_tbl$n_mechanisms_on_edge <- integer(0)
    return(edge_tbl)
  }
  
  split_key <- interaction(edge_tbl$from_group, edge_tbl$to_group, drop = TRUE, lex.order = TRUE)
  edge_split <- split(edge_tbl, split_key)
  
  out <- lapply(edge_split, function(d) {
    d <- d[order(-d$weight_plot, -d$weight_raw, d$mechanism), , drop = FALSE]
    
    n <- nrow(d)
    d$rank_within_edge <- seq_len(n)
    d$offset_index <- seq_len(n) - 1L
    d$n_mechanisms_on_edge <- n
    
    d$alpha_plot <- cummin(d$alpha_plot)
    d$intensity_plot <- d$alpha_plot
    
    d$alpha_plot[!is.finite(d$alpha_plot)] <- 0
    d$intensity_plot[!is.finite(d$intensity_plot)] <- 0
    
    d
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  
  out <- out[order(out$weight_plot, out$weight_raw, out$from_group, out$to_group, out$mechanism), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Build node table for family circuit plotting
computeFamilyCircuitNodeTable <- function(eff,
                                          group_by,
                                          drop_na_groups = TRUE,
                                          drop_empty_groups = TRUE,
                                          start_angle = pi / 2,
                                          clockwise = TRUE,
                                          size_nodes_by_abundance = FALSE,
                                          circle_radius = 1,
                                          label_radius = 1.20,
                                          node_radius = 0.06,
                                          node_radius_range = c(0.025, 0.11)) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  grp_out <- resolveCircuitGroups(
    eff = eff,
    group_by = group_by,
    drop_na = drop_na_groups,
    drop_empty = drop_empty_groups
  )
  
  node_tbl <- computeCircuitNodeTable(
    groups = grp_out$group,
    levels = grp_out$levels,
    start_angle = start_angle,
    clockwise = clockwise,
    circle_radius = circle_radius,
    label_radius = label_radius,
    size_nodes_by_abundance = size_nodes_by_abundance,
    node_radius = node_radius,
    node_radius_range = node_radius_range
  )
  
  list(
    node_table = node_tbl,
    groups = grp_out$group,
    levels = grp_out$levels
  )
}

#' Determine which side of a base lane is the convex/outward side
#'
#' Returns +1 if left-normal is outward, -1 if right-normal is outward.
inferStackedCircuitOutwardSide <- function(curve_df,
                                   center = c(0, 0),
                                   probe_offset = 0.01) {
  req <- c("x", "y", "point_id")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  
  d <- curve_df[order(curve_df$point_id), , drop = FALSE]
  if (nrow(d) < 2L) return(1)
  
  xy <- as.matrix(d[, c("x", "y"), drop = FALSE])
  storage.mode(xy) <- "double"
  
  tang <- .circuit_estimate_tangent(xy)
  n_left <- cbind(x = -tang[, "y"], y = tang[, "x"])
  
  p_left <- cbind(
    x = xy[, 1] + probe_offset * n_left[, "x"],
    y = xy[, 2] + probe_offset * n_left[, "y"]
  )
  p_right <- cbind(
    x = xy[, 1] - probe_offset * n_left[, "x"],
    y = xy[, 2] - probe_offset * n_left[, "y"]
  )
  
  d0 <- sqrt((xy[, 1] - center[1])^2 + (xy[, 2] - center[2])^2)
  d_left <- sqrt((p_left[, 1] - center[1])^2 + (p_left[, 2] - center[2])^2)
  d_right <- sqrt((p_right[, 1] - center[1])^2 + (p_right[, 2] - center[2])^2)
  
  gain_left <- mean(d_left - d0, na.rm = TRUE)
  gain_right <- mean(d_right - d0, na.rm = TRUE)
  
  if (!is.finite(gain_left)) gain_left <- -Inf
  if (!is.finite(gain_right)) gain_right <- -Inf
  
  if (gain_left >= gain_right) 1 else -1
}

#' Expand base directional lanes into stacked mechanism-specific lanes
#'
#' Thin wrapper over buildStackedCircuitCurves().
buildFamilyStackedCurves <- function(base_curve_df,
                                     edge_tbl,
                                     node_tbl,
                                     track_spacing,
                                     edge_gap = 0.015,
                                     self_loop_spread = 0.62,
                                     self_loop_height = 0.5,
                                     n_pts = 401,
                                     self_anchor_frac = 0.16) {
  req_base <- c("from_group", "to_group", "point_id", "x", "y")
  req_edge <- c(
    "from_group", "to_group", "mechanism",
    "weight_raw", "weight_plot", "intensity_plot",
    "rank_within_edge", "offset_index", "n_mechanisms_on_edge"
  )
  req_node <- c("group", "x", "y", "theta", "radius_plot")
  
  if (!is.data.frame(base_curve_df) || !all(req_base %in% names(base_curve_df))) {
    stop("base_curve_df must contain: ", paste(req_base, collapse = ", "))
  }
  if (!is.data.frame(edge_tbl) || !all(req_edge %in% names(edge_tbl))) {
    stop("edge_tbl must contain: ", paste(req_edge, collapse = ", "))
  }
  if (!is.data.frame(node_tbl) || !all(req_node %in% names(node_tbl))) {
    stop("node_tbl must contain: ", paste(req_node, collapse = ", "))
  }
  if (!is.numeric(track_spacing) || length(track_spacing) != 1L ||
      !is.finite(track_spacing) || track_spacing < 0) {
    stop("track_spacing must be a single non-negative number.")
  }
  if (!is.numeric(edge_gap) || length(edge_gap) != 1L ||
      !is.finite(edge_gap) || edge_gap < 0) {
    stop("edge_gap must be a single non-negative number.")
  }
  if (!is.numeric(self_loop_spread) || length(self_loop_spread) != 1L ||
      !is.finite(self_loop_spread) || self_loop_spread <= 0) {
    stop("self_loop_spread must be a single positive number.")
  }
  if (!is.numeric(self_loop_height) || length(self_loop_height) != 1L ||
      !is.finite(self_loop_height) || self_loop_height <= 0) {
    stop("self_loop_height must be a single positive number.")
  }
  if (!is.numeric(n_pts) || length(n_pts) != 1L ||
      !is.finite(n_pts) || n_pts < 50) {
    stop("n_pts must be a single integer >= 50.")
  }
  if (!is.numeric(self_anchor_frac) || length(self_anchor_frac) != 1L ||
      !is.finite(self_anchor_frac) || self_anchor_frac <= 0 || self_anchor_frac >= 0.5) {
    stop("self_anchor_frac must be a single number in (0, 0.5).")
  }
  
  edge_tbl$.ord_weight_plot_desc <- -edge_tbl$weight_plot
  edge_tbl$.ord_weight_raw_desc  <- -edge_tbl$weight_raw
  
  out <- buildStackedCircuitCurves(
    base_curve_df = base_curve_df,
    edge_tbl = edge_tbl,
    node_tbl = node_tbl,
    track_spacing = track_spacing,
    attach_metadata_fun = function(curve_df, edge_row, track_id) {
      attachFamilyCurveMetadata(
        curve_df = curve_df,
        edge_row = edge_row,
        family_edge_id = track_id
      )
    },
    track_order_cols = c("offset_index", ".ord_weight_plot_desc", ".ord_weight_raw_desc", "mechanism"),
    offset_col = "offset_index",
    count_col = "n_mechanisms_on_edge",
    id_start = 0L,
    edge_gap = edge_gap,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    n_pts = n_pts,
    self_anchor_frac = self_anchor_frac
  )
  
  if (nrow(out) == 0L) {
    return(out)
  }
  
  out <- out[order(out$weight_raw, out$track_id, out$point_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

anchorCurveToTargetEndpoints <- function(curve_df,
                                         target_df,
                                         anchor_frac = 0.16) {
  req <- c("x", "y", "point_id")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  if (!is.data.frame(target_df) || !all(req %in% names(target_df))) {
    stop("target_df must contain: ", paste(req, collapse = ", "))
  }
  
  d <- curve_df[order(curve_df$point_id), , drop = FALSE]
  t <- target_df[order(target_df$point_id), , drop = FALSE]
  
  if (nrow(d) != nrow(t)) {
    stop("curve_df and target_df must have the same number of rows.")
  }
  if (!is.numeric(anchor_frac) || length(anchor_frac) != 1L ||
      !is.finite(anchor_frac) || anchor_frac <= 0 || anchor_frac >= 0.5) {
    stop("anchor_frac must be a single number in (0, 0.5).")
  }
  
  n <- nrow(d)
  if (n < 3L) return(d)
  
  k <- max(2L, floor(anchor_frac * n))
  k <- min(k, floor(n / 2))
  
  delta_start <- c(
    t$x[1] - d$x[1],
    t$y[1] - d$y[1]
  )
  delta_end <- c(
    t$x[n] - d$x[n],
    t$y[n] - d$y[n]
  )
  
  w_start <- numeric(n)
  w_end <- numeric(n)
  
  # cosine taper: 1 -> 0 across first k points
  u1 <- seq(0, 1, length.out = k)
  w_start[seq_len(k)] <- 0.5 * (1 + cos(pi * u1))
  
  # cosine taper: 0 -> 1 across last k points
  u2 <- seq(0, 1, length.out = k)
  idx_end <- (n - k + 1L):n
  w_end[idx_end] <- 0.5 * (1 - cos(pi * u2))
  
  d$x <- d$x + w_start * delta_start[1] + w_end * delta_end[1]
  d$y <- d$y + w_start * delta_start[2] + w_end * delta_end[2]
  
  # enforce exact endpoint anchoring
  d$x[1] <- t$x[1]
  d$y[1] <- t$y[1]
  d$x[n] <- t$x[n]
  d$y[n] <- t$y[n]
  
  d
}


attachFamilyCurveMetadata <- function(curve_df, edge_row, family_edge_id) {
  req_edge <- c(
    "from_group", "to_group", "mechanism",
    "weight_raw", "weight_plot", "intensity_plot",
    "rank_within_edge", "offset_index", "n_mechanisms_on_edge"
  )
  
  if (!is.data.frame(edge_row) || nrow(edge_row) != 1L || !all(req_edge %in% names(edge_row))) {
    stop("edge_row must be a one-row data.frame containing: ", paste(req_edge, collapse = ", "))
  }
  
  out <- attachStackedCurveMetadata(
    curve_df = curve_df,
    edge_row = edge_row,
    id_col = "family_edge_id",
    id_value = family_edge_id,
    extra_cols = req_edge
  )
  
  out$track_id <- out$family_edge_id
  
  out <- out[, c(
    "x", "y", "point_id",
    "from_group", "to_group",
    "arrow_side", "is_self",
    "mechanism",
    "weight_raw", "weight_plot", "intensity_plot",
    "rank_within_edge", "offset_index", "n_mechanisms_on_edge",
    "family_edge_id", "track_id"
  ), drop = FALSE]
  
  rownames(out) <- NULL
  out
}

#' Resolve family mechanism identity aesthetics
#'
#' Updated behavior:
#' - high-saturation colors by default
#' - mechanisms sharing the same ligand share a color when possible
#' - black is allowed and used first when helpful
#' - linetypes prioritize visually dense styles with limited whitespace
#'
#' Identity is still the combination of color + linetype.
resolveFamilyMechanismAesthetics <- function(
    mechanisms,
    mechanism_colors = NULL,
    mechanism_linetypes = NULL,
    linetype_cycle = c(
      "solid",
      "longdash",
      "dashed",
      "dotdash",
      "twodash",
      "dotted"
    ),
    solid_color_palette = c(
      "black",
      "#8B1E1E", # dark red
      "#0B3C8C", # dark blue
      "#0B6E4F", # dark green
      "#6A1B73", # dark purple
      "#8C5A00", # dark gold / ochre
      "#5A3E36", # dark brown
      "#3F3F3F"  # dark gray
    )
) {
  mechanisms <- normalizeIdentityLabels(mechanisms)
  
  if (length(mechanisms) < 1L) {
    stop("No mechanisms available for aesthetic resolution.")
  }
  
  n_mech <- length(mechanisms)
  
  parse_ligand <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("\u2014", "-", x, fixed = TRUE)
    x <- gsub("\u2013", "-", x, fixed = TRUE)
    x <- gsub("\u2212", "-", x, fixed = TRUE)
    
    out <- sub("-.*$", "", x)
    out[!nzchar(out)] <- x[!nzchar(out)]
    out
  }
  
  ligand <- parse_ligand(mechanisms)
  ligand_levels <- unique(ligand)
  n_lig <- length(ligand_levels)
  
  # ------------------------------------------------------------
  # USER-SUPPLIED AESTHETICS
  # ------------------------------------------------------------
  if (!is.null(mechanism_colors)) {
    cols <- resolveGroupedIdentityAestheticVector(
      labels = mechanisms,
      groups = ligand,
      values = mechanism_colors,
      arg_name = "mechanism_colors"
    )
  } else {
    cols <- NULL
  }
  
  if (!is.null(mechanism_linetypes)) {
    ltys <- resolveGroupedIdentityAestheticVector(
      labels = mechanisms,
      groups = ligand,
      values = mechanism_linetypes,
      arg_name = "mechanism_linetypes"
    )
  } else {
    ltys <- NULL
  }
  
  if (!is.null(cols) && !is.null(ltys)) {
    return(list(
      colors = cols,
      linetypes = ltys
    ))
  }
  
  # normalize palette
  solid_color_palette <- as.character(solid_color_palette)
  solid_color_palette <- solid_color_palette[!is.na(solid_color_palette) & nzchar(solid_color_palette)]
  solid_color_palette <- unique(solid_color_palette)
  
  if (length(solid_color_palette) < 1L) {
    stop("solid_color_palette must contain at least one valid color.")
  }
  
  non_solid_cycle <- setdiff(linetype_cycle, "solid")
  if (length(non_solid_cycle) < 1L) {
    non_solid_cycle <- c("longdash", "dashed", "dotdash", "twodash", "dotted")
  }
  
  # ------------------------------------------------------------
  # CASE 1: SINGLE LIGAND
  # PRIORITY:
  #   1. use as many distinct SOLID colors as possible
  #   2. only then begin reusing colors with non-solid linetypes
  # ------------------------------------------------------------
  if (n_lig == 1L) {
    if (is.null(cols)) {
      cols <- rep(NA_character_, n_mech)
      names(cols) <- mechanisms
      
      n_solid_capacity <- length(solid_color_palette)
      
      if (n_mech <= n_solid_capacity) {
        cols[] <- solid_color_palette[seq_len(n_mech)]
      } else {
        combo_colors <- rep(
          solid_color_palette,
          times = ceiling(n_mech / n_solid_capacity)
        )[seq_len(n_mech)]
        cols[] <- combo_colors
      }
    }
    
    if (is.null(ltys)) {
      ltys <- rep(NA_character_, n_mech)
      names(ltys) <- mechanisms
      
      n_solid_capacity <- length(solid_color_palette)
      
      if (n_mech <= n_solid_capacity) {
        ltys[] <- "solid"
      } else {
        combo_ltys <- c(
          rep("solid", n_solid_capacity),
          rep(
            rep(non_solid_cycle, each = n_solid_capacity),
            length.out = max(0, n_mech - n_solid_capacity)
          )
        )
        ltys[] <- combo_ltys[seq_len(n_mech)]
      }
    }
    
    return(list(
      colors = cols,
      linetypes = ltys
    ))
  }
  
  # ------------------------------------------------------------
  # CASE 2: MULTIPLE LIGANDS
  # PRIORITY:
  #   1. preserve shared color within ligand
  #   2. use solid first within each ligand
  #   3. only use non-solid if a ligand has many mechanisms
  # ------------------------------------------------------------
  if (is.null(cols)) {
    ligand_palette <- solid_color_palette
    
    if (length(ligand_levels) <= length(ligand_palette)) {
      ligand_cols <- ligand_palette[seq_len(length(ligand_levels))]
    } else {
      extra_needed <- length(ligand_levels) - length(ligand_palette)
      ligand_cols <- c(
        ligand_palette,
        grDevices::hcl.colors(extra_needed, palette = "Dark 3", rev = FALSE)
      )
    }
    
    names(ligand_cols) <- ligand_levels
    cols <- ligand_cols[ligand]
    names(cols) <- mechanisms
  }
  
  if (is.null(ltys)) {
    ltys <- character(n_mech)
    names(ltys) <- mechanisms
    
    mech_by_lig <- split(mechanisms, ligand)
    
    for (lig in names(mech_by_lig)) {
      mech_i <- mech_by_lig[[lig]]
      n_i <- length(mech_i)
      
      if (n_i == 1L) {
        ltys[mech_i] <- "solid"
      } else {
        ltys_i <- c(
          "solid",
          rep(non_solid_cycle, length.out = max(0, n_i - 1L))
        )
        ltys[mech_i] <- ltys_i[seq_len(n_i)]
      }
    }
  }
  
  list(
    colors = cols,
    linetypes = ltys
  )
}

#' Mix one or more colors toward white by intensity in [0, 1]
#'
#' intensity = 1 retains the original color
#' intensity = 0 becomes pure white
whitenFamilyColors <- function(colors,
                               intensity) {
  if (!is.character(colors)) {
    stop("colors must be a character vector.")
  }
  if (!is.numeric(intensity)) {
    stop("intensity must be numeric.")
  }
  
  if (length(intensity) == 1L) {
    intensity <- rep(intensity, length(colors))
  }
  if (length(colors) != length(intensity)) {
    stop("colors and intensity must have the same length, or intensity must have length 1.")
  }
  
  intensity[!is.finite(intensity)] <- 0
  intensity <- pmax(0, pmin(1, intensity))
  
  rgb_in <- grDevices::col2rgb(colors) / 255
  out <- vapply(seq_along(colors), function(i) {
    rgb_i <- rgb_in[, i]
    z <- intensity[i]
    rgb_out <- (1 - z) * c(1, 1, 1) + z * rgb_i
    grDevices::rgb(rgb_out[1], rgb_out[2], rgb_out[3])
  }, character(1))
  
  out
}

#' Add mechanism identity aesthetics and plotting colors to stacked curve data
prepareFamilyCurveAesthetics <- function(curve_df,
                                         mechanism_aesthetics) {
  req <- c("mechanism", "intensity_plot", "track_id")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  if (!is.list(mechanism_aesthetics) ||
      is.null(mechanism_aesthetics$colors) ||
      is.null(mechanism_aesthetics$linetypes)) {
    stop("mechanism_aesthetics must be a list with $colors and $linetypes.")
  }
  
  cols <- mechanism_aesthetics$colors
  ltys <- mechanism_aesthetics$linetypes
  
  miss_col <- setdiff(unique(curve_df$mechanism), names(cols))
  miss_lty <- setdiff(unique(curve_df$mechanism), names(ltys))
  
  if (length(miss_col) > 0L) {
    stop("Missing colors for mechanism(s): ", paste(miss_col, collapse = ", "), ".")
  }
  if (length(miss_lty) > 0L) {
    stop("Missing linetypes for mechanism(s): ", paste(miss_lty, collapse = ", "), ".")
  }
  
  curve_df$base_color <- unname(cols[curve_df$mechanism])
  curve_df$linetype_mech <- unname(ltys[curve_df$mechanism])
  
  if (!"alpha_plot" %in% names(curve_df)) {
    curve_df$alpha_plot <- curve_df$intensity_plot
  }
  curve_df$alpha_plot[!is.finite(curve_df$alpha_plot)] <- 0
  curve_df$alpha_plot <- pmax(0, pmin(1, curve_df$alpha_plot))
  
  curve_df$identity <- curve_df$mechanism
  curve_df$identity_color <- whitenFamilyColors(
    colors = curve_df$base_color,
    intensity = curve_df$alpha_plot
  )
  curve_df$identity_linetype <- curve_df$linetype_mech
  
  curve_df
}

#' Approximate data-unit track spacing from user line width
#'
#' This is the conversion layer needed because ggplot line widths live in device space
#' while stacked curve offsets live in data coordinates.
#'
#' The returned spacing is intended to make adjacent mechanism tracks visually tangent.
estimateStackedCircuitTrackSpacing <- function(edge_linewidth = 0.7,
                                       xlim = c(-1.40, 1.40),
                                       ylim = c(-1.40, 1.40),
                                       spacing_factor = 0.00525) {
  if (!is.numeric(edge_linewidth) || length(edge_linewidth) != 1L ||
      !is.finite(edge_linewidth) || edge_linewidth <= 0) {
    stop("edge_linewidth must be a single positive number.")
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || any(!is.finite(xlim)) || xlim[1] >= xlim[2]) {
    stop("xlim must be a numeric vector of length 2 with xlim[1] < xlim[2].")
  }
  if (!is.numeric(ylim) || length(ylim) != 2L || any(!is.finite(ylim)) || ylim[1] >= ylim[2]) {
    stop("ylim must be a numeric vector of length 2 with ylim[1] < ylim[2].")
  }
  if (!is.numeric(spacing_factor) || length(spacing_factor) != 1L ||
      !is.finite(spacing_factor) || spacing_factor <= 0) {
    stop("spacing_factor must be a single positive number.")
  }
  
  span <- max(diff(xlim), diff(ylim))
  edge_linewidth * span * spacing_factor
}

#' High-level family circuit data engine (helpers only; no plotting yet)
computeFamilyCircuitData <- function(eff,
                                     group_by,
                                     lr,
                                     deduplicate_lr = TRUE,
                                     include_self_loops = TRUE,
                                     include_same_cell_edges = TRUE,
                                     drop_na_groups = TRUE,
                                     drop_empty_groups = TRUE,
                                     size_nodes_by_abundance = FALSE,
                                     start_angle = pi / 2,
                                     clockwise = TRUE,
                                     min_mechanism_edge = 0,
                                     mechanism_relative_threshold = 0.25,
                                     visible_floor = 0.12,
                                     visible_ceiling = 1,
                                     max_edges = 2e6) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  node_out <- computeFamilyCircuitNodeTable(
    eff = eff,
    group_by = group_by,
    drop_na_groups = drop_na_groups,
    drop_empty_groups = drop_empty_groups,
    start_angle = start_angle,
    clockwise = clockwise,
    size_nodes_by_abundance = size_nodes_by_abundance
  )
  
  edge_out <- computeFamilyCircuitEdgeTable(
    eff = eff,
    groups = node_out$groups,
    lr = lr,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    deduplicate_lr = deduplicate_lr,
    max_edges = max_edges
  )
  
  edge_tbl <- prepareFamilyCircuitEdgeTable(
    edge_tbl = edge_out$edge_table,
    min_mechanism_edge = min_mechanism_edge,
    mechanism_relative_threshold = mechanism_relative_threshold,
    visible_floor = visible_floor,
    visible_ceiling = visible_ceiling
  )
  
  edge_tbl <- rankFamilyMechanismsWithinEdges(edge_tbl)
  
  list(
    nodes = node_out$node_table,
    edges = edge_tbl,
    groups = node_out$groups,
    levels = node_out$levels,
    selected_mechanisms = edge_out$selected_mechanisms,
    params = list(
      group_by = group_by,
      lr = lr,
      deduplicate_lr = deduplicate_lr,
      include_self_loops = include_self_loops,
      include_same_cell_edges = include_same_cell_edges,
      drop_na_groups = drop_na_groups,
      drop_empty_groups = drop_empty_groups,
      size_nodes_by_abundance = size_nodes_by_abundance,
      start_angle = start_angle,
      clockwise = clockwise,
      min_mechanism_edge = min_mechanism_edge,
      mechanism_relative_threshold = mechanism_relative_threshold,
      visible_floor = visible_floor,
      visible_ceiling = visible_ceiling,
      max_edges = max_edges
    )
  )
}


#' Render LR family circuit plot from prepared node and stacked mechanism-curve tables
#'
#' Updated legend design:
#' - manual legend column now includes:
#'     * signal strength
#'     * mechanism
#'     * group
#'     * optional group abundance
#' - all legend sections share one visual grammar
renderLRFamilyCircuitPlot <- function(node_tbl,
                                      curve_df,
                                      plot_title = NULL,
                                      node_palette = NULL,
                                      node_outline = "black",
                                      edge_linewidth = 0.7,
                                      label_size = 4,
                                      xlim = c(-1.40, 1.40),
                                      ylim = c(-1.40, 1.40),
                                      show_node_size_legend = NULL,
                                      node_size_legend_title = "Group abundance",
                                      node_size_legend_range = c(4, 10),
                                      signal_legend_title = "Signal strength",
                                      signal_legend_n = 5,
                                      mechanism_legend_title = "Mechanism",
                                      group_legend_title = "Group",
                                      legend_width_ratio = 1.45,
                                      arrow_head_length = 0.060,
                                      arrow_head_width = 0.022,
                                      taper_length = 0.060,
                                      taper_base_halfwidth = 0.0035) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Package 'ggforce' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }
  
  req_node <- c("group", "x", "y", "radius_plot", "label_x", "label_y", "hjust")
  if (!is.data.frame(node_tbl) || !all(req_node %in% names(node_tbl))) {
    stop("node_tbl must contain: ", paste(req_node, collapse = ", "))
  }
  
  req_curve <- c(
    "track_id", "point_id", "x", "y",
    "from_group", "to_group",
    "weight_raw", "weight_plot", "intensity_plot",
    "identity", "identity_color", "identity_linetype",
    "base_color",
    "arrow_side", "is_self"
  )
  
  if (!is.data.frame(curve_df) || !all(req_curve %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req_curve, collapse = ", "))
  }
  
  if (is.null(show_node_size_legend)) {
    show_node_size_legend <- length(unique(round(node_tbl$radius_plot, 10))) > 1L
  }
  
  p_main <- ggplot2::ggplot()
  
  if (nrow(curve_df) > 0L) {
    validateStackedRenderColumns(curve_df)
    geom_out <- prepareStackedCircuitRenderGeometry(
      curve_df = curve_df,
      order_cols = c("weight_plot", "weight_raw", "track_id", "point_id"),
      arrow_head_length = arrow_head_length,
      arrow_head_width = arrow_head_width,
      taper_length = taper_length,
      taper_base_halfwidth = taper_base_halfwidth
    )
    
    p_main <- renderStackedCircuitTracksFixed(
      p = p_main,
      edge_order_df = geom_out$edge_order_df,
      shaft_df = geom_out$shaft_df,
      taper_poly = geom_out$taper_poly,
      arrow_poly = geom_out$arrow_poly,
      edge_linewidth = edge_linewidth
    )
  }
  
  p_main <- renderCircuitNodesBasic(
    p = p_main,
    node_tbl = node_tbl,
    node_palette = node_palette,
    node_outline = node_outline,
    label_size = label_size
  )
  
  p_main <- p_main +
    ggplot2::coord_equal(xlim = xlim, ylim = ylim, clip = "off") +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(hjust = 0.5, color = "black")
    )
  
  if (nrow(curve_df) > 0L) {
    p_signal_legend <- buildFamilySignalLegend(
      curve_df = curve_df,
      signal_legend_title = signal_legend_title,
      signal_legend_n = signal_legend_n
    )
    
    p_identity_legend <- buildFamilyMechanismLegend(
      curve_df = curve_df,
      mechanism_legend_title = mechanism_legend_title,
      edge_linewidth = edge_linewidth
    )
  } else {
    p_signal_legend <- makeEmptyCircuitLegend()
    p_identity_legend <- makeEmptyCircuitLegend()
  }
  
  legend_col <- assembleStackedCircuitLegends(
    node_tbl = node_tbl,
    node_palette = node_palette,
    node_outline = node_outline,
    signal_legend = p_signal_legend,
    identity_legend = p_identity_legend,
    show_node_size_legend = show_node_size_legend,
    node_size_legend_title = node_size_legend_title,
    node_size_legend_range = node_size_legend_range,
    group_legend_title = group_legend_title,
    heights_no_size = c(1, 1.5, 1),
    heights_with_size = c(1, 1.5, 1, 1)
  )
  
  out <- p_main | legend_col
  out + patchwork::plot_layout(widths = c(3.8, legend_width_ratio))
}

#### plotFamilyCircuit Public Functions ####

prepareFamilyCircuitPlotData <- function(
    eff,
    group_by,
    lr,
    deduplicate_lr = TRUE,
    include_self_loops = TRUE,
    include_same_cell_edges = TRUE,
    drop_na_groups = TRUE,
    drop_empty_groups = TRUE,
    size_nodes_by_abundance = FALSE,
    start_angle = pi / 2,
    clockwise = TRUE,
    min_mechanism_edge = 0,
    mechanism_relative_threshold = 0.25,
    visible_floor = 0.12,
    visible_ceiling = 1,
    node_palette = NULL,
    mechanism_colors = NULL,
    mechanism_linetypes = NULL,
    linetype_cycle = c("solid", "longdash", "dashed", "twodash", "dotdash", "dotted"),
    edge_gap = 0.015,
    bidirectional_offset = NULL,
    self_loop_spread = 0.5,
    self_loop_height = 0.5,
    curve_strength = 0.16,
    edge_linewidth = 0.7,
    xlim = c(-1.40, 1.40),
    ylim = c(-1.40, 1.40),
    track_spacing = NULL,
    track_spacing_factor = 0.00525,
    n_curve_pts = 401,
    center = c(0, 0),
    max_edges = 2e6
) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  fam <- computeFamilyCircuitData(
    eff = eff,
    group_by = group_by,
    lr = lr,
    deduplicate_lr = deduplicate_lr,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    drop_na_groups = drop_na_groups,
    drop_empty_groups = drop_empty_groups,
    size_nodes_by_abundance = size_nodes_by_abundance,
    start_angle = start_angle,
    clockwise = clockwise,
    min_mechanism_edge = min_mechanism_edge,
    mechanism_relative_threshold = mechanism_relative_threshold,
    visible_floor = visible_floor,
    visible_ceiling = visible_ceiling,
    max_edges = max_edges
  )
  
  base_curve_df <- buildStackedCircuitBaseCurves(
    node_tbl = fam$nodes,
    edge_tbl = fam$edges,
    edge_gap = edge_gap,
    bidirectional_offset = bidirectional_offset,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    curve_strength = curve_strength,
    n_pts = n_curve_pts,
    center = center
  )
  
  if (is.null(track_spacing)) {
    track_spacing <- estimateStackedCircuitTrackSpacing(
      edge_linewidth = edge_linewidth,
      xlim = xlim,
      ylim = ylim,
      spacing_factor = track_spacing_factor
    )
  } else {
    if (!is.numeric(track_spacing) || length(track_spacing) != 1L ||
        !is.finite(track_spacing) || track_spacing < 0) {
      stop("track_spacing must be NULL or a single non-negative number.")
    }
  }
  
  curve_df <- buildFamilyStackedCurves(
    base_curve_df = base_curve_df,
    edge_tbl = fam$edges,
    node_tbl = fam$nodes,
    track_spacing = track_spacing,
    edge_gap = edge_gap,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    n_pts = n_curve_pts
  )  
  
  mech_aes <- resolveFamilyMechanismAesthetics(
    mechanisms = fam$selected_mechanisms,
    mechanism_colors = mechanism_colors,
    mechanism_linetypes = mechanism_linetypes,
    linetype_cycle = linetype_cycle
  )
  
  if (nrow(curve_df) > 0L) {
    curve_df <- prepareFamilyCurveAesthetics(
      curve_df = curve_df,
      mechanism_aesthetics = mech_aes
    )
  }
  
  list(
    nodes = fam$nodes,
    edges = fam$edges,
    base_curves = base_curve_df,
    curves = curve_df,
    selected_mechanisms = fam$selected_mechanisms,
    mechanism_aesthetics = mech_aes,
    params = c(
      fam$params,
      list(
        node_palette = node_palette,
        mechanism_colors = mechanism_colors,
        mechanism_linetypes = mechanism_linetypes,
        linetype_cycle = linetype_cycle,
        edge_gap = edge_gap,
        bidirectional_offset = bidirectional_offset,
        self_loop_spread = self_loop_spread,
        self_loop_height = self_loop_height,
        curve_strength = curve_strength,
        edge_linewidth = edge_linewidth,
        xlim = xlim,
        ylim = ylim,
        track_spacing = track_spacing,
        track_spacing_factor = track_spacing_factor,
        n_curve_pts = n_curve_pts,
        center = center
      )
    )
  )
}

#' Plot group-level LR family circuit
plotFamilyCircuit <- function(
    eff,
    group_by,
    lr,
    deduplicate_lr = TRUE,
    include_self_loops = TRUE,
    include_same_cell_edges = TRUE,
    drop_na_groups = TRUE,
    drop_empty_groups = TRUE,
    size_nodes_by_abundance = FALSE,
    node_palette = NULL,
    node_outline = "black",
    mechanism_colors = NULL,
    mechanism_linetypes = NULL,
    linetype_cycle = c("solid", "longdash", "dashed", "twodash", "dotdash", "dotted"),
    min_mechanism_edge = 0,
    mechanism_relative_threshold = 0.25,
    visible_floor = 0.12,
    visible_ceiling = 1,
    edge_gap = 0.015,
    bidirectional_offset = NULL,
    self_loop_spread = 0.5,
    self_loop_height = 0.5,
    curve_strength = 0.16,
    edge_linewidth = 0.7,
    start_angle = pi / 2,
    clockwise = TRUE,
    label_size = 4,
    plot_title = NULL,
    title_suffix = NULL,
    xlim = c(-1.40, 1.40),
    ylim = c(-1.40, 1.40),
    track_spacing = NULL,
    track_spacing_factor = 0.00525,
    n_curve_pts = 401,
    center = c(0, 0),
    show_node_size_legend = NULL,
    node_size_legend_title = "Group abundance",
    node_size_legend_range = c(4, 10),
    signal_legend_title = "Signal strength",
    signal_legend_n = 5,
    mechanism_legend_title = "Mechanism",
    group_legend_title = "Group",
    legend_width_ratio = 1.45,
    arrow_head_length = 0.060,
    arrow_head_width = 0.022,
    taper_length = 0.060,
    taper_base_halfwidth = 0.0035,
    max_edges = 2e6
) {
  stopifnot(inherits(eff, "EffNICHES"))
  
  fam_plot <- prepareFamilyCircuitPlotData(
    eff = eff,
    group_by = group_by,
    lr = lr,
    deduplicate_lr = deduplicate_lr,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    drop_na_groups = drop_na_groups,
    drop_empty_groups = drop_empty_groups,
    size_nodes_by_abundance = size_nodes_by_abundance,
    start_angle = start_angle,
    clockwise = clockwise,
    min_mechanism_edge = min_mechanism_edge,
    mechanism_relative_threshold = mechanism_relative_threshold,
    visible_floor = visible_floor,
    visible_ceiling = visible_ceiling,
    node_palette = node_palette,
    mechanism_colors = mechanism_colors,
    mechanism_linetypes = mechanism_linetypes,
    linetype_cycle = linetype_cycle,
    edge_gap = edge_gap,
    bidirectional_offset = bidirectional_offset,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    curve_strength = curve_strength,
    edge_linewidth = edge_linewidth,
    xlim = xlim,
    ylim = ylim,
    track_spacing = track_spacing,
    track_spacing_factor = track_spacing_factor,
    n_curve_pts = n_curve_pts,
    center = center,
    max_edges = max_edges
  )
  
  title_use <- makeLRPlotTitle(
    selected_mechanisms = fam_plot$selected_mechanisms,
    plot_title = plot_title,
    suffix = title_suffix
  )
  
  renderLRFamilyCircuitPlot(
    node_tbl = fam_plot$nodes,
    curve_df = fam_plot$curves,
    plot_title = title_use,
    node_palette = node_palette,
    node_outline = node_outline,
    edge_linewidth = edge_linewidth,
    label_size = label_size,
    xlim = xlim,
    ylim = ylim,
    show_node_size_legend = show_node_size_legend,
    node_size_legend_title = node_size_legend_title,
    node_size_legend_range = node_size_legend_range,
    signal_legend_title = signal_legend_title,
    signal_legend_n = signal_legend_n,
    mechanism_legend_title = mechanism_legend_title,
    group_legend_title = group_legend_title,
    legend_width_ratio = legend_width_ratio,
    arrow_head_length = arrow_head_length,
    arrow_head_width = arrow_head_width,
    taper_length = taper_length,
    taper_base_halfwidth = taper_base_halfwidth
  )
}

#' Multi-plot wrapper for LR family circuit plots
#'
#' `lr` can be:
#' - a character vector, where each element is treated as its own plot, or
#' - a list, where each element is a vector of LR mechanisms defining one family plot
plotFamilyCircuitMulti <- function(
    eff,
    group_by,
    lr,
    plot_titles = NULL,
    ncol = NULL,
    nrow = NULL,
    guides = "collect",
    ...
) {
  .plotLRMulti(
    plot_fun = function(eff, lr, plot_title, ...) {
      plotFamilyCircuit(
        eff = eff,
        group_by = group_by,
        lr = lr,
        plot_title = plot_title,
        ...
      )
    },
    eff = eff,
    lr = lr,
    plot_titles = plot_titles,
    ncol = ncol,
    nrow = nrow,
    guides = guides,
    ...
  )
}

#### plotCompareCircuit Helpers ####

#' Resolve a common union of groups across multiple EffNICHES objects
resolveCompareCircuitGroups <- function(eff_list,
                                        group_by,
                                        drop_na = TRUE,
                                        drop_empty = TRUE,
                                        warn_if_nonidentical = TRUE) {
  if (!is.list(eff_list) || length(eff_list) < 2L) {
    stop("eff_list must be a list containing at least 2 EffNICHES objects.")
  }
  if (is.null(names(eff_list))) {
    names(eff_list) <- paste0("Object ", seq_along(eff_list))
  } else {
    bad <- !nzchar(names(eff_list))
    names(eff_list)[bad] <- paste0("Object ", which(bad))
  }
  
  grp_each <- vector("list", length(eff_list))
  lev_each <- vector("list", length(eff_list))
  
  for (i in seq_along(eff_list)) {
    eff <- eff_list[[i]]
    stopifnot(inherits(eff, "EffNICHES"))
    
    out_i <- resolveCircuitGroups(
      eff = eff,
      group_by = group_by,
      drop_na = drop_na,
      drop_empty = drop_empty
    )
    grp_each[[i]] <- out_i$group
    lev_each[[i]] <- unique(out_i$group[!is.na(out_i$group)])
  }
  
  union_levels <- unique(unlist(lev_each, use.names = FALSE))
  union_levels <- sort(union_levels)
  
  grp_sets <- lapply(lev_each, function(x) sort(unique(as.character(x))))
  identical_sets <- all(vapply(
    grp_sets,
    function(x) identical(x, grp_sets[[1]]),
    logical(1)
  ))
  
  if (isTRUE(warn_if_nonidentical) && !identical_sets) {
    warning(
      "Objects do not have identical group sets under group_by = '",
      group_by,
      "'. Using the union of groups across objects."
    )
  }
  
  list(
    groups_each = grp_each,
    levels_union = union_levels,
    object_names = names(eff_list),
    identical_group_sets = identical_sets
  )
}

#' Compute a common node table from the union of groups across objects
computeCompareCircuitNodeTable <- function(eff_list,
                                           group_by,
                                           start_angle = pi / 2,
                                           clockwise = TRUE,
                                           circle_radius = 1,
                                           label_radius = 1.20,
                                           size_nodes_by_abundance = FALSE,
                                           node_radius = 0.06,
                                           node_radius_range = c(0.025, 0.11),
                                           node_abundance_source = c("mean", "sum", "first"),
                                           drop_na_groups = TRUE,
                                           drop_empty_groups = TRUE,
                                           warn_if_nonidentical = TRUE) {
  node_abundance_source <- match.arg(node_abundance_source)
  
  grp_out <- resolveCompareCircuitGroups(
    eff_list = eff_list,
    group_by = group_by,
    drop_na = drop_na_groups,
    drop_empty = drop_empty_groups,
    warn_if_nonidentical = warn_if_nonidentical
  )
  
  levels_union <- grp_out$levels_union
  if (length(levels_union) < 1L) {
    stop("No non-missing groups were found across eff_list.")
  }
  
  counts_mat <- sapply(grp_out$groups_each, function(g) {
    as.numeric(table(factor(g[!is.na(g)], levels = levels_union)))
  })
  if (is.null(dim(counts_mat))) {
    counts_mat <- matrix(counts_mat, ncol = 1L)
  }
  rownames(counts_mat) <- levels_union
  colnames(counts_mat) <- grp_out$object_names
  
  # Object-specific within-object proportions
  totals_by_object <- colSums(counts_mat)
  frac_mat <- sweep(counts_mat, 2, totals_by_object, FUN = "/")
  frac_mat[!is.finite(frac_mat)] <- 0
  
  # Always initialize these so downstream code is safe
  delta_frac <- rep(NA_real_, length(levels_union))
  is_neutral <- rep(NA, length(levels_union))
  
  # Node radii
  if (isTRUE(size_nodes_by_abundance)) {
    if (length(grp_out$object_names) != 2L) {
      stop(
        "The compare-mode node abundance glyph currently requires exactly 2 objects ",
        "when size_nodes_by_abundance = TRUE."
      )
    }
    
    all_frac_vals <- as.numeric(frac_mat)
    all_frac_vals <- all_frac_vals[is.finite(all_frac_vals)]
    
    if (length(all_frac_vals) < 1L) {
      radius_mat <- matrix(
        node_radius_range[1],
        nrow = nrow(frac_mat),
        ncol = ncol(frac_mat),
        dimnames = dimnames(frac_mat)
      )
    } else {
      frac_rng <- range(all_frac_vals, na.rm = TRUE)
      
      if (diff(frac_rng) == 0) {
        radius_mat <- matrix(
          node_radius_range[2],
          nrow = nrow(frac_mat),
          ncol = ncol(frac_mat),
          dimnames = dimnames(frac_mat)
        )
      } else {
        radius_mat <- apply(frac_mat, 2, function(v) {
          node_radius_range[1] +
            ((v - frac_rng[1]) / diff(frac_rng)) * diff(node_radius_range)
        })
        radius_mat <- as.matrix(radius_mat)
        rownames(radius_mat) <- rownames(frac_mat)
        colnames(radius_mat) <- colnames(frac_mat)
      }
    }
    
    obj1 <- grp_out$object_names[1]
    obj2 <- grp_out$object_names[2]
    
    r1 <- as.numeric(radius_mat[, obj1])
    r2 <- as.numeric(radius_mat[, obj2])
    
    outer_is_obj1 <- r1 >= r2
    outer_object <- ifelse(outer_is_obj1, obj1, obj2)
    inner_object <- ifelse(outer_is_obj1, obj2, obj1)
    
    radius_outer <- pmax(r1, r2)
    radius_inner <- pmin(r1, r2)
    
    # ------------------------------------------------------------
    # Determine dominance + neutral threshold behavior
    # ------------------------------------------------------------
    neutral_threshold <- 0.04  # default
    
    frac1 <- frac_mat[, obj1]
    frac2 <- frac_mat[, obj2]
    
    delta_frac <- abs(frac1 - frac2)
    
    is_neutral <- delta_frac < neutral_threshold
    
    outer_is_obj1 <- frac1 >= frac2
    outer_object <- ifelse(outer_is_obj1, obj1, obj2)
    inner_object <- ifelse(outer_is_obj1, obj2, obj1)
    
    # Override with neutral flag
    outer_object[is_neutral] <- "neutral"
    
    # For backward compatibility and node-size legend:
    n_cells <- radius_outer
    frac_cells <- radius_outer
    
    radius_plot <- radius_outer
  } else {
    # Backward-compatible pooled node abundance behavior
    n_cells <- switch(
      node_abundance_source,
      mean = rowMeans(counts_mat),
      sum = rowSums(counts_mat),
      first = counts_mat[, 1]
    )
    
    n_cells <- as.numeric(n_cells)
    names(n_cells) <- levels_union
    
    frac_cells <- if (sum(n_cells) > 0) n_cells / sum(n_cells) else rep(0, length(n_cells))
    
    if (diff(range(frac_cells)) == 0) {
      radius_plot <- rep(node_radius_range[2], length(frac_cells))
    } else if (isTRUE(size_nodes_by_abundance)) {
      frac01 <- (frac_cells - min(frac_cells)) / diff(range(frac_cells))
      radius_plot <- node_radius_range[1] + frac01 * diff(node_radius_range)
    } else {
      radius_plot <- rep(node_radius, length(levels_union))
    }
    
    radius_outer <- radius_plot
    radius_inner <- NA_real_
    outer_object <- NA_character_
    inner_object <- NA_character_
    
    if (ncol(frac_mat) >= 1L) {
      rtmp <- rep(NA_real_, nrow(frac_mat))
      names(rtmp) <- rownames(frac_mat)
      r1 <- rtmp
      r2 <- rtmp
      if (ncol(frac_mat) >= 1L) r1 <- as.numeric(frac_mat[, 1])
      if (ncol(frac_mat) >= 2L) r2 <- as.numeric(frac_mat[, 2])
    } else {
      r1 <- rep(NA_real_, length(levels_union))
      r2 <- rep(NA_real_, length(levels_union))
    }
  }
  
  n <- length(levels_union)
  step <- 2 * pi / n
  theta <- start_angle + if (isTRUE(clockwise)) -(0:(n - 1)) * step else (0:(n - 1)) * step
  
  x <- circle_radius * cos(theta)
  y <- circle_radius * sin(theta)
  
  label_x <- (label_radius + radius_plot) * cos(theta)
  label_y <- (label_radius + radius_plot) * sin(theta)
  hjust <- ifelse(label_x >= 0, 0, 1)
  
  node_tbl <- data.frame(
    group = levels_union,
    n_cells = as.numeric(n_cells),
    frac_cells = as.numeric(frac_cells),
    theta = theta,
    x = x,
    y = y,
    radius_plot = as.numeric(radius_plot),
    label_x = label_x,
    label_y = label_y,
    hjust = hjust,
    stringsAsFactors = FALSE
  )
  
  # Add compare-mode abundance fields when exactly 2 objects are present
  if (length(grp_out$object_names) >= 2L) {
    node_tbl$count_obj1 <- as.numeric(counts_mat[, 1])
    node_tbl$count_obj2 <- as.numeric(counts_mat[, 2])
    node_tbl$frac_obj1  <- as.numeric(frac_mat[, 1])
    node_tbl$frac_obj2  <- as.numeric(frac_mat[, 2])
    node_tbl$radius_obj1 <- if (exists("r1")) as.numeric(if (all(is.na(r1))) NA_real_ else if (isTRUE(size_nodes_by_abundance)) radius_mat[, 1] else NA_real_) else NA_real_
    node_tbl$radius_obj2 <- if (exists("r2")) as.numeric(if (all(is.na(r2))) NA_real_ else if (isTRUE(size_nodes_by_abundance)) radius_mat[, 2] else NA_real_) else NA_real_
    node_tbl$radius_outer <- as.numeric(radius_outer)
    node_tbl$radius_inner <- as.numeric(radius_inner)
    node_tbl$outer_object <- as.character(outer_object)
    node_tbl$inner_object <- as.character(inner_object)
    node_tbl$delta_frac <- as.numeric(delta_frac)
    node_tbl$is_neutral <- as.logical(is_neutral)
  }
  
  list(
    node_table = node_tbl,
    groups_each = grp_out$groups_each,
    object_names = grp_out$object_names,
    identical_group_sets = grp_out$identical_group_sets
  )
}

#' Relevel an edge table to a supplied common group ordering
relevelCircuitEdgeTableToUnion <- function(edge_tbl,
                                           levels_union,
                                           include_self_loops = TRUE,
                                           include_same_cell_edges = TRUE,
                                           node_group_sizes = NULL) {
  if (!is.data.frame(edge_tbl) ||
      !all(c("from_group", "to_group", "weight_raw") %in% names(edge_tbl))) {
    stop("edge_tbl must contain from_group, to_group, and weight_raw.")
  }
  
  all_pairs <- expand.grid(
    from_group = levels_union,
    to_group = levels_union,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!isTRUE(include_self_loops)) {
    all_pairs <- all_pairs[all_pairs$from_group != all_pairs$to_group, , drop = FALSE]
  }
  
  out <- merge(
    all_pairs,
    edge_tbl,
    by = c("from_group", "to_group"),
    all.x = TRUE,
    sort = FALSE
  )
  
  if (!"n_possible_edges" %in% names(out)) out$n_possible_edges <- NA_real_
  if (!"n_constructed_edges" %in% names(out)) out$n_constructed_edges <- NA_real_
  if (!"n_realized_edges" %in% names(out)) out$n_realized_edges <- 0
  if (!"n_sender_participating" %in% names(out)) out$n_sender_participating <- 0
  if (!"n_receiver_participating" %in% names(out)) out$n_receiver_participating <- 0
  
  out$weight_raw[is.na(out$weight_raw)] <- 0
  out$n_realized_edges[is.na(out$n_realized_edges)] <- 0
  out$n_sender_participating[is.na(out$n_sender_participating)] <- 0
  out$n_receiver_participating[is.na(out$n_receiver_participating)] <- 0
  
  if (!is.null(node_group_sizes)) {
    out$n_sender_total <- unname(node_group_sizes[out$from_group])
    out$n_receiver_total <- unname(node_group_sizes[out$to_group])
    
    same_group <- out$from_group == out$to_group
    out$n_possible_edges <- ifelse(
      same_group & !isTRUE(include_same_cell_edges),
      out$n_sender_total * pmax(out$n_receiver_total - 1, 0),
      out$n_sender_total * out$n_receiver_total
    )
  }
  
  out
}

#' Compute one-object / one-mechanism circuit edge table on a common group union
computeSingleObjectCompareEdgeTable <- function(eff,
                                                groups,
                                                levels_union,
                                                lr,
                                                transform = identity,
                                                lr_normalization = c("none", "max", "percentile", "zscore"),
                                                edge_group_aggregate = c(
                                                  "sum",
                                                  "mean",
                                                  "mean_realized",
                                                  "mean_per_sender",
                                                  "mean_per_receiver"
                                                ),
                                                include_self_loops = TRUE,
                                                include_same_cell_edges = TRUE,
                                                max_edges = 2e6) {
  lr_normalization <- match.arg(lr_normalization)
  edge_group_aggregate <- match.arg(edge_group_aggregate)
  
  edge_out <- computeCircuitEdgeTable(
    eff = eff,
    groups = groups,
    lr = lr,
    transform = transform,
    lr_normalization = lr_normalization,
    lr_aggregation = "sum",
    deduplicate_lr = TRUE,
    warn_many_lr = FALSE,
    edge_group_aggregate = edge_group_aggregate,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    max_edges = max_edges
  )
  
  node_group_sizes <- as.numeric(table(factor(groups[!is.na(groups)], levels = levels_union)))
  names(node_group_sizes) <- levels_union
  
  relevelCircuitEdgeTableToUnion(
    edge_tbl = edge_out$edge_table,
    levels_union = levels_union,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    node_group_sizes = node_group_sizes
  )
}

#' Normalize comparison edge weights across objects
normalizeCompareCircuitWeights <- function(edge_tbl,
                                           object_normalization = c("shared_raw", "shared_max", "per_object_max")) {
  object_normalization <- match.arg(object_normalization)
  
  req <- c("object_name", "weight_raw")
  if (!is.data.frame(edge_tbl) || !all(req %in% names(edge_tbl))) {
    stop("edge_tbl must contain: ", paste(req, collapse = ", "))
  }
  
  out <- edge_tbl
  
  if (object_normalization == "shared_raw") {
    out$weight_plot <- out$weight_raw
    return(out)
  }
  
  if (object_normalization == "shared_max") {
    mx <- suppressWarnings(max(out$weight_raw, na.rm = TRUE))
    out$weight_plot <- if (is.finite(mx) && mx > 0) out$weight_raw / mx else 0
    out$weight_plot[!is.finite(out$weight_plot)] <- 0
    return(out)
  }
  
  if (object_normalization == "per_object_max") {
    obj_max <- tapply(out$weight_raw, out$object_name, max, na.rm = TRUE)
    out$weight_plot <- ifelse(
      obj_max[out$object_name] > 0,
      out$weight_raw / obj_max[out$object_name],
      0
    )
    out$weight_plot[!is.finite(out$weight_plot)] <- 0
    return(out)
  }
  
  out
}

#' Build comparison circuit data for one mechanism across multiple objects
computeCompareCircuitData <- function(eff_list,
                                      group_by,
                                      lr,
                                      transform = identity,
                                      lr_normalization = c("none", "max", "percentile", "zscore"),
                                      edge_group_aggregate = c(
                                        "sum",
                                        "mean",
                                        "mean_realized",
                                        "mean_per_sender",
                                        "mean_per_receiver"
                                      ),
                                      object_normalization = c("shared_raw", "shared_max", "per_object_max"),
                                      include_self_loops = TRUE,
                                      include_same_cell_edges = TRUE,
                                      drop_na_groups = TRUE,
                                      drop_empty_groups = TRUE,
                                      start_angle = pi / 2,
                                      clockwise = TRUE,
                                      size_nodes_by_abundance = FALSE,
                                      node_abundance_source = c("mean", "sum", "first"),
                                      warn_if_nonidentical = TRUE,
                                      max_edges = 2e6) {
  if (!is.list(eff_list) || length(eff_list) < 2L) {
    stop("eff_list must contain at least 2 EffNICHES objects.")
  }
  if (length(lr) != 1L) {
    stop("This comparison plot is for one LR mechanism at a time; lr must have length 1.")
  }
  
  lr_normalization <- match.arg(lr_normalization)
  edge_group_aggregate <- match.arg(edge_group_aggregate)
  object_normalization <- match.arg(object_normalization)
  node_abundance_source <- match.arg(node_abundance_source)
  
  node_out <- computeCompareCircuitNodeTable(
    eff_list = eff_list,
    group_by = group_by,
    start_angle = start_angle,
    clockwise = clockwise,
    size_nodes_by_abundance = size_nodes_by_abundance,
    node_abundance_source = node_abundance_source,
    drop_na_groups = drop_na_groups,
    drop_empty_groups = drop_empty_groups,
    warn_if_nonidentical = warn_if_nonidentical
  )
  
  node_tbl <- node_out$node_table
  levels_union <- node_tbl$group
  object_names <- node_out$object_names
  
  edge_list <- vector("list", length(eff_list))
  
  for (i in seq_along(eff_list)) {
    edge_i <- computeSingleObjectCompareEdgeTable(
      eff = eff_list[[i]],
      groups = node_out$groups_each[[i]],
      levels_union = levels_union,
      lr = lr,
      transform = transform,
      lr_normalization = lr_normalization,
      edge_group_aggregate = edge_group_aggregate,
      include_self_loops = include_self_loops,
      include_same_cell_edges = include_same_cell_edges,
      max_edges = max_edges
    )
    
    edge_i$object_name <- object_names[i]
    edge_i$object_index <- i
    edge_list[[i]] <- edge_i
  }
  
  edge_tbl <- do.call(rbind, edge_list)
  rownames(edge_tbl) <- NULL
  
  edge_tbl <- normalizeCompareCircuitWeights(
    edge_tbl = edge_tbl,
    object_normalization = object_normalization
  )
  
  edge_tbl <- edge_tbl[order(
    edge_tbl$from_group,
    edge_tbl$to_group,
    edge_tbl$object_index
  ), , drop = FALSE]
  rownames(edge_tbl) <- NULL
  
  list(
    nodes = node_tbl,
    edges = edge_tbl,
    object_names = object_names,
    identical_group_sets = node_out$identical_group_sets,
    params = list(
      group_by = group_by,
      lr = lr,
      lr_normalization = lr_normalization,
      edge_group_aggregate = edge_group_aggregate,
      object_normalization = object_normalization,
      include_self_loops = include_self_loops,
      include_same_cell_edges = include_same_cell_edges,
      drop_na_groups = drop_na_groups,
      drop_empty_groups = drop_empty_groups,
      start_angle = start_angle,
      clockwise = clockwise,
      size_nodes_by_abundance = size_nodes_by_abundance,
      node_abundance_source = node_abundance_source,
      max_edges = max_edges
    )
  )
}

#' Add fixed lane order metadata for object-comparison edges
rankCompareObjectsWithinEdges <- function(edge_tbl,
                                          object_names) {
  req <- c("from_group", "to_group", "object_name", "object_index", "weight_raw", "weight_plot")
  if (!is.data.frame(edge_tbl) || !all(req %in% names(edge_tbl))) {
    stop("edge_tbl must contain: ", paste(req, collapse = ", "))
  }
  
  out <- edge_tbl
  out$object_name <- factor(out$object_name, levels = object_names)
  out <- out[order(out$from_group, out$to_group, out$object_index), , drop = FALSE]
  
  split_key <- interaction(out$from_group, out$to_group, drop = TRUE, lex.order = TRUE)
  edge_split <- split(out, split_key)
  
  out2 <- lapply(edge_split, function(d) {
    d <- d[order(d$object_index), , drop = FALSE]
    d$rank_within_edge <- seq_len(nrow(d))
    d$offset_index <- seq_len(nrow(d)) - 1L
    d$n_objects_on_edge <- nrow(d)
    d
  })
  
  out2 <- do.call(rbind, out2)
  rownames(out2) <- NULL
  out2
}

#' Resolve aesthetics for object-comparison tracks
resolveCompareObjectAesthetics <- function(object_names,
                                           object_colors = NULL,
                                           object_linetypes = NULL) {
  object_names <- normalizeIdentityLabels(object_names)
  
  if (length(object_names) < 1L) {
    stop("No object names available for aesthetic resolution.")
  }
  
  cols <- if (is.null(object_colors)) {
    defaultCompareObjectColors(object_names)
  } else {
    resolveIdentityAestheticVector(
      labels = object_names,
      values = object_colors,
      arg_name = "object_colors"
    )
  }
  
  ltys <- if (is.null(object_linetypes)) {
    defaultCompareObjectLinetypes(object_names)
  } else {
    resolveIdentityAestheticVector(
      labels = object_names,
      values = object_linetypes,
      arg_name = "object_linetypes"
    )
  }
  
  list(colors = cols, linetypes = ltys)
}

#' Prepare compare-curve aesthetics
prepareCompareCurveAesthetics <- function(curve_df,
                                          object_aesthetics) {
  req <- c("object_name", "weight_plot", "track_id")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain object_name, weight_plot, and track_id.")
  }
  
  cols <- object_aesthetics$colors
  ltys <- object_aesthetics$linetypes
  
  miss_col <- setdiff(unique(as.character(curve_df$object_name)), names(cols))
  miss_lty <- setdiff(unique(as.character(curve_df$object_name)), names(ltys))
  
  if (length(miss_col) > 0L) {
    stop("Missing object colors for: ", paste(miss_col, collapse = ", "))
  }
  if (length(miss_lty) > 0L) {
    stop("Missing object linetypes for: ", paste(miss_lty, collapse = ", "))
  }
  
  curve_df$base_color <- unname(cols[as.character(curve_df$object_name)])
  curve_df$linetype_obj <- unname(ltys[as.character(curve_df$object_name)])
  
  curve_df$identity <- as.character(curve_df$object_name)
  curve_df$identity_color <- curve_df$base_color
  curve_df$identity_linetype <- curve_df$linetype_obj
  
  curve_df
}

#' Attach object-comparison metadata to curve geometry
attachCompareCurveMetadata <- function(curve_df, edge_row, compare_edge_id) {
  req_edge <- c(
    "from_group", "to_group", "object_name", "object_index",
    "weight_raw", "weight_plot",
    "rank_within_edge", "offset_index", "n_objects_on_edge"
  )
  
  if (!is.data.frame(edge_row) || nrow(edge_row) != 1L || !all(req_edge %in% names(edge_row))) {
    stop("edge_row must be a one-row data.frame containing: ", paste(req_edge, collapse = ", "))
  }
  
  out <- attachStackedCurveMetadata(
    curve_df = curve_df,
    edge_row = edge_row,
    id_col = "compare_edge_id",
    id_value = compare_edge_id,
    extra_cols = req_edge
  )
  
  out$track_id <- out$compare_edge_id
  
  rownames(out) <- NULL
  out
}

#' Expand base directional lanes into stacked object-comparison lanes
#'
#' Thin wrapper over buildStackedCircuitCurves().
buildCompareStackedCurves <- function(base_curve_df,
                                      edge_tbl,
                                      node_tbl,
                                      track_spacing,
                                      edge_gap = 0.015,
                                      self_loop_spread = 0.62,
                                      self_loop_height = 0.5,
                                      n_pts = 401,
                                      self_anchor_frac = 0.16) {
  req_base <- c("from_group", "to_group", "point_id", "x", "y")
  req_edge <- c(
    "from_group", "to_group", "object_name", "object_index",
    "weight_raw", "weight_plot",
    "rank_within_edge", "offset_index", "n_objects_on_edge"
  )
  req_node <- c("group", "x", "y", "theta", "radius_plot")
  
  if (!is.data.frame(base_curve_df) || !all(req_base %in% names(base_curve_df))) {
    stop("base_curve_df must contain: ", paste(req_base, collapse = ", "))
  }
  if (!is.data.frame(edge_tbl) || !all(req_edge %in% names(edge_tbl))) {
    stop("edge_tbl must contain: ", paste(req_edge, collapse = ", "))
  }
  if (!is.data.frame(node_tbl) || !all(req_node %in% names(node_tbl))) {
    stop("node_tbl must contain: ", paste(req_node, collapse = ", "))
  }
  
  out <- buildStackedCircuitCurves(
    base_curve_df = base_curve_df,
    edge_tbl = edge_tbl,
    node_tbl = node_tbl,
    track_spacing = track_spacing,
    attach_metadata_fun = function(curve_df, edge_row, track_id) {
      attachCompareCurveMetadata(
        curve_df = curve_df,
        edge_row = edge_row,
        compare_edge_id = track_id
      )
    },
    track_order_cols = c("object_index"),
    offset_col = "offset_index",
    count_col = "n_objects_on_edge",
    id_start = 0L,
    edge_gap = edge_gap,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    n_pts = n_pts,
    self_anchor_frac = self_anchor_frac
  )
  
  if (nrow(out) == 0L) {
    return(out)
  }
  
  out <- out[order(out$weight_raw, out$track_id, out$point_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Build all prepared data needed for comparison circuit rendering
prepareCompareCircuitPlotData <- function(eff_list,
                                          group_by,
                                          lr,
                                          transform = identity,
                                          lr_normalization = c("none", "max", "percentile", "zscore"),
                                          edge_group_aggregate = c(
                                            "sum",
                                            "mean",
                                            "mean_realized",
                                            "mean_per_sender",
                                            "mean_per_receiver"
                                          ),
                                          object_normalization = c("shared_raw", "shared_max", "per_object_max"),
                                          include_self_loops = TRUE,
                                          include_same_cell_edges = TRUE,
                                          drop_na_groups = TRUE,
                                          drop_empty_groups = TRUE,
                                          size_nodes_by_abundance = FALSE,
                                          node_abundance_source = c("mean", "sum", "first"),
                                          start_angle = pi / 2,
                                          clockwise = TRUE,
                                          node_palette = NULL,
                                          object_colors = NULL,
                                          object_linetypes = NULL,
                                          edge_gap = 0.015,
                                          bidirectional_offset = NULL,
                                          self_loop_spread = 0.5,
                                          self_loop_height = 0.5,
                                          curve_strength = 0.16,
                                          edge_linewidth = 0.7,
                                          xlim = c(-1.40, 1.40),
                                          ylim = c(-1.40, 1.40),
                                          track_spacing = NULL,
                                          track_spacing_factor = 0.00525,
                                          n_curve_pts = 401,
                                          center = c(0, 0),
                                          warn_if_nonidentical = TRUE,
                                          max_edges = 2e6) {
  lr_normalization <- match.arg(lr_normalization)
  edge_group_aggregate <- match.arg(edge_group_aggregate)
  object_normalization <- match.arg(object_normalization)
  node_abundance_source <- match.arg(node_abundance_source)
  
  cmp <- computeCompareCircuitData(
    eff_list = eff_list,
    group_by = group_by,
    lr = lr,
    transform = transform,
    lr_normalization = lr_normalization,
    edge_group_aggregate = edge_group_aggregate,
    object_normalization = object_normalization,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    drop_na_groups = drop_na_groups,
    drop_empty_groups = drop_empty_groups,
    start_angle = start_angle,
    clockwise = clockwise,
    size_nodes_by_abundance = size_nodes_by_abundance,
    node_abundance_source = node_abundance_source,
    warn_if_nonidentical = warn_if_nonidentical,
    max_edges = max_edges
  )
  
  cmp$edges <- rankCompareObjectsWithinEdges(
    edge_tbl = cmp$edges,
    object_names = cmp$object_names
  )
  
  base_curve_df <- buildStackedCircuitBaseCurves(
    node_tbl = cmp$nodes,
    edge_tbl = cmp$edges,
    edge_gap = edge_gap,
    bidirectional_offset = bidirectional_offset,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    curve_strength = curve_strength,
    n_pts = n_curve_pts,
    center = center
  )
  
  if (is.null(track_spacing)) {
    track_spacing <- estimateStackedCircuitTrackSpacing(
      edge_linewidth = edge_linewidth,
      xlim = xlim,
      ylim = ylim,
      spacing_factor = track_spacing_factor
    )
  }
  
  curve_df <- buildCompareStackedCurves(
    base_curve_df = base_curve_df,
    edge_tbl = cmp$edges,
    node_tbl = cmp$nodes,
    track_spacing = track_spacing,
    edge_gap = edge_gap,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    n_pts = n_curve_pts
  )
  
  obj_aes <- resolveCompareObjectAesthetics(
    object_names = cmp$object_names,
    object_colors = object_colors,
    object_linetypes = object_linetypes
  )
  
  if (nrow(curve_df) > 0L) {
    curve_df <- prepareCompareCurveAesthetics(
      curve_df = curve_df,
      object_aesthetics = obj_aes
    )
  }
  
  list(
    nodes = cmp$nodes,
    edges = cmp$edges,
    base_curves = base_curve_df,
    curves = curve_df,
    object_names = cmp$object_names,
    object_aesthetics = obj_aes,
    identical_group_sets = cmp$identical_group_sets,
    params = c(
      cmp$params,
      list(
        node_palette = node_palette,
        object_colors = object_colors,
        object_linetypes = object_linetypes,
        edge_gap = edge_gap,
        bidirectional_offset = bidirectional_offset,
        self_loop_spread = self_loop_spread,
        self_loop_height = self_loop_height,
        curve_strength = curve_strength,
        edge_linewidth = edge_linewidth,
        xlim = xlim,
        ylim = ylim,
        track_spacing = track_spacing,
        track_spacing_factor = track_spacing_factor,
        n_curve_pts = n_curve_pts,
        center = center
      )
    )
  )
}

buildCompareObjectLegend <- function(curve_df,
                                     object_legend_title = "Object",
                                     edge_linewidth = 0.7,
                                     text_size = 3.4,
                                     title_size = 11,
                                     x_label = 1.35,
                                     xlim = c(0, 3.2)) {
  req <- c("identity", "base_color", "identity_linetype")
  if (!is.data.frame(curve_df) || !all(req %in% names(curve_df))) {
    stop("curve_df must contain: ", paste(req, collapse = ", "))
  }
  
  leg_df <- curve_df
  leg_df$identity_linetype <- as.character(leg_df$identity_linetype)
  
  buildCircuitIdentityLegend(
    df = leg_df,
    label_col = "identity",
    color_col = "base_color",
    linetype_col = "identity_linetype",
    legend_title = object_legend_title,
    edge_linewidth = edge_linewidth,
    text_size = text_size,
    title_size = title_size,
    x_label = x_label,
    xlim = xlim,
    sort_labels = TRUE
  )
}

renderCompareCircuitPlot <- function(node_tbl,
                                     curve_df,
                                     plot_title = NULL,
                                     node_palette = NULL,
                                     node_outline = "black",
                                     object_aesthetics = NULL,
                                     edge_aes = c("alpha", "width"),
                                     edge_alpha_range = c(0.01, 1),
                                     alpha_gamma = 2.2,
                                     edge_width_range = c(0.25, 2.6),
                                     edge_linewidth = 0.7,
                                     label_size = 4,
                                     xlim = c(-1.40, 1.40),
                                     ylim = c(-1.40, 1.40),
                                     show_node_size_legend = NULL,
                                     node_size_legend_title = "Group abundance",
                                     node_size_legend_range = c(4, 10),
                                     signal_legend_title = "Signal strength",
                                     object_legend_title = "Object",
                                     group_legend_title = "Group",
                                     legend_width_ratio = 1.45,
                                     arrow_head_length = 0.060,
                                     arrow_head_width = 0.022,
                                     taper_length = 0.060,
                                     taper_base_halfwidth = 0.0035) {
  edge_aes <- match.arg(edge_aes)
  
  if (is.null(show_node_size_legend)) {
    show_node_size_legend <- length(unique(round(node_tbl$radius_plot, 10))) > 1L
  }
  
  use_compare_node_glyph <- useCompareCircuitNodeGlyph(node_tbl)
  
  p_main <- ggplot2::ggplot()
  
  if (nrow(curve_df) > 0L) {
    validateStackedRenderColumns(curve_df)
    req_curve_extra <- c("base_color")
    if (!all(req_curve_extra %in% names(curve_df))) {
      stop("curve_df must contain: ", paste(req_curve_extra, collapse = ", "))
    }
    geom_out <- prepareStackedCircuitRenderGeometry(
      curve_df = curve_df,
      order_cols = c("weight_plot", "weight_raw", "track_id", "point_id"),
      arrow_head_length = arrow_head_length,
      arrow_head_width = arrow_head_width,
      taper_length = taper_length,
      taper_base_halfwidth = taper_base_halfwidth
    )
    
    p_main <- renderStackedCircuitTracksMapped(
      p = p_main,
      edge_order_df = geom_out$edge_order_df,
      shaft_df = geom_out$shaft_df,
      taper_poly = geom_out$taper_poly,
      arrow_poly = geom_out$arrow_poly,
      edge_aes = edge_aes,
      edge_alpha_range = edge_alpha_range,
      alpha_gamma = alpha_gamma,
      edge_width_range = edge_width_range,
      edge_linewidth = edge_linewidth
    )
  }
  
  if (isTRUE(use_compare_node_glyph)) {
    p_main <- renderCircuitNodesCompareGlyph(
      p = p_main,
      node_tbl = node_tbl,
      object_aesthetics = object_aesthetics,
      node_palette = node_palette,
      node_outline = node_outline,
      label_size = label_size
    )
  } else {
    p_main <- renderCircuitNodesBasic(
      p = p_main,
      node_tbl = node_tbl,
      node_palette = node_palette,
      node_outline = node_outline,
      label_size = label_size
    )
  }
  
  p_main <- p_main +
    ggplot2::coord_equal(xlim = xlim, ylim = ylim, clip = "off") +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(hjust = 0.5, color = "black")
    )
  
  signal_values <- if (nrow(curve_df) > 0L) curve_df$weight_plot else c(0, 1)
  
  if (edge_aes == "alpha") {
    p_signal_legend <- buildCircuitSignalLegendAlpha(
      values = signal_values,
      legend_title = signal_legend_title,
      alpha_range = edge_alpha_range,
      alpha_gamma = alpha_gamma,
      n = 5
    )
  } else {
    p_signal_legend <- buildCircuitSignalLegendWidth(
      values = signal_values,
      legend_title = signal_legend_title,
      edge_width_range = edge_width_range,
      edge_color = "black",
      edge_alpha = 0.9,
      n = 5
    )
  }
  
  p_identity_legend <- if (nrow(curve_df) > 0L) {
    buildCompareObjectLegend(
      curve_df = curve_df,
      object_legend_title = object_legend_title,
      edge_linewidth = edge_linewidth
    )
  } else {
    makeEmptyCircuitLegend()
  }
  
  legend_col <- assembleStackedCircuitLegends(
    node_tbl = node_tbl,
    node_palette = node_palette,
    node_outline = node_outline,
    signal_legend = p_signal_legend,
    identity_legend = p_identity_legend,
    show_node_size_legend = show_node_size_legend,
    node_size_legend_title = node_size_legend_title,
    node_size_legend_range = node_size_legend_range,
    group_legend_title = group_legend_title,
    heights_no_size = c(1, 1.1, 1),
    heights_with_size = c(1, 1.1, 1, 1)
  )
  
  (p_main | legend_col) + patchwork::plot_layout(widths = c(3.8, legend_width_ratio))
}

#### plotCompareCircuit Public ####

plotCompareCircuit <- function(eff_list,
                               group_by,
                               lr,
                               transform = identity,
                               lr_normalization = c("none", "max", "percentile", "zscore"),
                               edge_group_aggregate = c(
                                 "sum",
                                 "mean",
                                 "mean_realized",
                                 "mean_per_sender",
                                 "mean_per_receiver"
                               ),
                               object_normalization = c("shared_raw", "shared_max", "per_object_max"),
                               include_self_loops = TRUE,
                               include_same_cell_edges = TRUE,
                               drop_na_groups = TRUE,
                               drop_empty_groups = TRUE,
                               size_nodes_by_abundance = FALSE,
                               node_abundance_source = c("mean", "sum", "first"),
                               node_palette = NULL,
                               node_outline = "black",
                               object_colors = NULL,
                               object_linetypes = NULL,
                               edge_aes = c("alpha", "width"),
                               edge_gap = 0.015,
                               bidirectional_offset = NULL,
                               self_loop_spread = 0.5,
                               self_loop_height = 0.5,
                               edge_alpha_range = c(0.01, 1),
                               alpha_gamma = 2.2,
                               edge_width_range = c(0.25, 2.6),
                               edge_linewidth = 0.7,
                               min_edge = 0,
                               start_angle = pi / 2,
                               clockwise = TRUE,
                               curve_strength = 0.16,
                               label_size = 4,
                               plot_title = NULL,
                               title_suffix = NULL,
                               xlim = c(-1.40, 1.40),
                               ylim = c(-1.40, 1.40),
                               track_spacing = NULL,
                               track_spacing_factor = 0.00525,
                               n_curve_pts = 401,
                               center = c(0, 0),
                               warn_if_nonidentical = TRUE,
                               show_node_size_legend = NULL,
                               node_size_legend_title = "Group abundance",
                               node_size_legend_range = c(4, 10),
                               signal_legend_title = "Signal strength",
                               object_legend_title = "Object",
                               group_legend_title = "Group",
                               legend_width_ratio = 1.45,
                               arrow_head_length = 0.060,
                               arrow_head_width = 0.022,
                               taper_length = 0.060,
                               taper_base_halfwidth = 0.0035,
                               max_edges = 2e6) {
  lr_normalization <- match.arg(lr_normalization)
  edge_group_aggregate <- match.arg(edge_group_aggregate)
  object_normalization <- match.arg(object_normalization)
  node_abundance_source <- match.arg(node_abundance_source)
  edge_aes <- match.arg(edge_aes)
  
  cmp_plot <- prepareCompareCircuitPlotData(
    eff_list = eff_list,
    group_by = group_by,
    lr = lr,
    transform = transform,
    lr_normalization = lr_normalization,
    edge_group_aggregate = edge_group_aggregate,
    object_normalization = object_normalization,
    include_self_loops = include_self_loops,
    include_same_cell_edges = include_same_cell_edges,
    drop_na_groups = drop_na_groups,
    drop_empty_groups = drop_empty_groups,
    size_nodes_by_abundance = size_nodes_by_abundance,
    node_abundance_source = node_abundance_source,
    start_angle = start_angle,
    clockwise = clockwise,
    node_palette = node_palette,
    object_colors = object_colors,
    object_linetypes = object_linetypes,
    edge_gap = edge_gap,
    bidirectional_offset = bidirectional_offset,
    self_loop_spread = self_loop_spread,
    self_loop_height = self_loop_height,
    curve_strength = curve_strength,
    edge_linewidth = edge_linewidth,
    xlim = xlim,
    ylim = ylim,
    track_spacing = track_spacing,
    track_spacing_factor = track_spacing_factor,
    n_curve_pts = n_curve_pts,
    center = center,
    warn_if_nonidentical = warn_if_nonidentical,
    max_edges = max_edges
  )
  
  if (!is.numeric(min_edge) || length(min_edge) != 1L || !is.finite(min_edge)) {
    stop("min_edge must be a single finite numeric value.")
  }
  
  curve_df <- cmp_plot$curves
  if (nrow(curve_df) > 0L) {
    curve_df <- curve_df[curve_df$weight_plot > min_edge, , drop = FALSE]
  }
  
  title_use <- makeLRPlotTitle(
    selected_mechanisms = lr,
    plot_title = plot_title,
    suffix = title_suffix
  )
  
  renderCompareCircuitPlot(
    node_tbl = cmp_plot$nodes,
    curve_df = curve_df,
    plot_title = title_use,
    node_palette = node_palette,
    node_outline = node_outline,
    object_aesthetics = cmp_plot$object_aesthetics,
    edge_aes = edge_aes,
    edge_alpha_range = edge_alpha_range,
    alpha_gamma = alpha_gamma,
    edge_width_range = edge_width_range,
    edge_linewidth = edge_linewidth,
    label_size = label_size,
    xlim = xlim,
    ylim = ylim,
    show_node_size_legend = show_node_size_legend,
    node_size_legend_title = node_size_legend_title,
    node_size_legend_range = node_size_legend_range,
    signal_legend_title = signal_legend_title,
    object_legend_title = object_legend_title,
    group_legend_title = group_legend_title,
    legend_width_ratio = legend_width_ratio,
    arrow_head_length = arrow_head_length,
    arrow_head_width = arrow_head_width,
    taper_length = taper_length,
    taper_base_halfwidth = taper_base_halfwidth
  )
}

