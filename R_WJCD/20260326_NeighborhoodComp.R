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

#### Neighborhood Comp Plot Dev ####

# Creation of Neighborhood Comp Matrix - no plot yet
computeNeighborhoodCompositionMatrix <- function(
    eff,
    subject_col,
    neighbor_col,
    edge_name = "CellToCellSpatial",
    zero_neighbor = c("drop", "zero_row"),
    include_self = TRUE,
    replicate_col = NULL,
    verbose = TRUE
) {
  zero_neighbor <- match.arg(zero_neighbor)
  
  # -------------------------
  # Normalize input: single EffNICHES or list of EffNICHES
  # -------------------------
  if (inherits(eff, "EffNICHES")) {
    eff_list <- list(eff)
    if (is.null(names(eff_list))) names(eff_list) <- "eff1"
  } else if (is.list(eff) && length(eff) > 0 &&
             all(vapply(eff, inherits, logical(1), "EffNICHES"))) {
    eff_list <- eff
    if (is.null(names(eff_list))) {
      names(eff_list) <- paste0("eff", seq_along(eff_list))
    } else {
      blank_names <- is.na(names(eff_list)) | names(eff_list) == ""
      names(eff_list)[blank_names] <- paste0("eff", which(blank_names))
    }
  } else {
    stop("eff must be either a single EffNICHES object or a non-empty list of EffNICHES objects.")
  }
  
  # -------------------------
  # Helpers
  # -------------------------
  .resolve_node_ids <- function(obj, obj_name) {
    node_ids <- NULL
    if (!is.null(obj$nodes$ids)) {
      node_ids <- as.character(obj$nodes$ids)
    } else if (!is.null(obj$nodes$xy) && !is.null(rownames(obj$nodes$xy))) {
      node_ids <- rownames(obj$nodes$xy)
    }
    if (is.null(node_ids)) {
      stop(sprintf("[%s] Could not resolve node IDs. Need eff$nodes$ids or rownames(eff$nodes$xy).", obj_name))
    }
    node_ids
  }
  
  .resolve_metadata <- function(obj, obj_name, node_ids, subject_col, neighbor_col, replicate_col = NULL) {
    md <- obj$nodes$meta
    if (is.null(md) || !is.data.frame(md)) {
      stop(sprintf("[%s] eff$nodes$meta must be a data.frame containing subject_col and neighbor_col.", obj_name))
    }
    if (is.null(rownames(md))) {
      stop(sprintf("[%s] eff$nodes$meta must have rownames matching node IDs for safe alignment.", obj_name))
    }
    
    missing_md <- setdiff(node_ids, rownames(md))
    if (length(missing_md) > 0) {
      stop(sprintf(
        "[%s] Metadata rownames do not cover all node IDs. Missing %d IDs (e.g. %s).",
        obj_name, length(missing_md), paste(head(missing_md, 5), collapse = ", ")
      ))
    }
    md <- md[node_ids, , drop = FALSE]
    
    if (!(subject_col %in% colnames(md))) {
      stop(sprintf("[%s] subject_col '%s' not found in eff$nodes$meta.", obj_name, subject_col))
    }
    if (!(neighbor_col %in% colnames(md))) {
      stop(sprintf("[%s] neighbor_col '%s' not found in eff$nodes$meta.", obj_name, neighbor_col))
    }
    if (!is.null(replicate_col) && !(replicate_col %in% colnames(md))) {
      stop(sprintf("[%s] replicate_col '%s' not found in eff$nodes$meta.", obj_name, replicate_col))
    }
    
    md
  }
  
  .resolve_ij <- function(obj, obj_name, edge_name, n_nodes) {
    if (is.null(obj$niches) || is.null(obj$niches[[edge_name]])) {
      stop(sprintf("[%s] Could not find eff$niches[['%s']].", obj_name, edge_name))
    }
    edge_obj <- obj$niches[[edge_name]]
    if (is.null(edge_obj$ij) || !is.matrix(edge_obj$ij) || ncol(edge_obj$ij) != 2) {
      stop(sprintf("[%s] eff$niches[['%s']]$ij must be an E x 2 matrix.", obj_name, edge_name))
    }
    
    ij <- edge_obj$ij
    if (!is.numeric(ij)) ij <- apply(ij, 2, as.integer)
    storage.mode(ij) <- "integer"
    
    if (any(is.na(ij))) {
      stop(sprintf("[%s] Edge index matrix ij contains NA values.", obj_name))
    }
    if (any(ij < 1L) || any(ij > n_nodes)) {
      stop(sprintf("[%s] Edge index matrix ij contains indices outside [1, n_nodes].", obj_name))
    }
    
    ij
  }
  
  .subset_and_remap_ij <- function(ij, keep_idx, n_nodes_total) {
    keep_flag <- rep(FALSE, n_nodes_total)
    keep_flag[keep_idx] <- TRUE
    
    if (nrow(ij) == 0) {
      return(matrix(integer(0), ncol = 2))
    }
    
    edge_keep <- keep_flag[ij[, 1]] & keep_flag[ij[, 2]]
    ij_sub <- ij[edge_keep, , drop = FALSE]
    
    if (nrow(ij_sub) == 0) {
      return(matrix(integer(0), ncol = 2))
    }
    
    map <- integer(n_nodes_total)
    map[keep_idx] <- seq_along(keep_idx)
    
    cbind(
      map[ij_sub[, 1]],
      map[ij_sub[, 2]]
    )
  }
  
  .compute_one_matrix <- function(
    ij,
    subject_chr,
    neighbor_chr,
    all_subject_levels,
    all_neighbor_levels,
    zero_neighbor,
    include_self
  ) {
    n_nodes <- length(subject_chr)
    n_edges <- if (length(ij) == 0) 0L else nrow(ij)
    
    subj_to_row <- setNames(seq_along(all_subject_levels), all_subject_levels)
    neigh_to_col <- setNames(seq_along(all_neighbor_levels), all_neighbor_levels)
    
    mat_sum <- matrix(
      0,
      nrow = length(all_subject_levels),
      ncol = length(all_neighbor_levels),
      dimnames = list(all_subject_levels, all_neighbor_levels)
    )
    
    n_cells_used <- setNames(integer(length(all_subject_levels)), all_subject_levels)
    n_cells_total <- setNames(integer(length(all_subject_levels)), all_subject_levels)
    n_cells_zero_neighbors <- setNames(integer(length(all_subject_levels)), all_subject_levels)
    
    # ---- Build undirected adjacency list (unweighted)
    if (n_edges > 0) {
      from <- ij[, 1]
      to   <- ij[, 2]
      
      from2 <- c(from, to)
      to2   <- c(to, from)
      neigh_list <- split(to2, from2)
    } else {
      neigh_list <- list()
    }
    
    empty <- vector("list", n_nodes)
    names(empty) <- as.character(seq_len(n_nodes))
    neigh_list_full <- empty
    if (length(neigh_list) > 0) {
      neigh_list_full[names(neigh_list)] <- neigh_list
    }
    
    # ---- Include/exclude self
    for (i in seq_len(n_nodes)) {
      nbr_idx <- neigh_list_full[[i]]
      if (is.null(nbr_idx)) nbr_idx <- integer(0)
      nbr_idx <- unique(as.integer(nbr_idx))
      
      if (isTRUE(include_self)) {
        nbr_idx <- unique(c(nbr_idx, i))
      } else {
        nbr_idx <- setdiff(nbr_idx, i)
      }
      
      neigh_list_full[[i]] <- nbr_idx
    }
    
    # ---- Aggregate per-cell compositions
    for (i in seq_len(n_nodes)) {
      s_lab <- subject_chr[i]
      if (is.na(s_lab)) next
      if (!(s_lab %in% all_subject_levels)) next
      
      n_cells_total[s_lab] <- n_cells_total[s_lab] + 1L
      
      nbr_idx <- neigh_list_full[[i]]
      if (length(nbr_idx) == 0) {
        n_cells_zero_neighbors[s_lab] <- n_cells_zero_neighbors[s_lab] + 1L
        if (zero_neighbor == "drop") {
          next
        } else {
          n_cells_used[s_lab] <- n_cells_used[s_lab] + 1L
          next
        }
      }
      
      nbr_lab <- neighbor_chr[nbr_idx]
      nbr_lab <- nbr_lab[!is.na(nbr_lab)]
      
      if (length(nbr_lab) == 0) {
        n_cells_zero_neighbors[s_lab] <- n_cells_zero_neighbors[s_lab] + 1L
        if (zero_neighbor == "drop") {
          next
        } else {
          n_cells_used[s_lab] <- n_cells_used[s_lab] + 1L
          next
        }
      }
      
      nbr_lab <- nbr_lab[nbr_lab %in% all_neighbor_levels]
      if (length(nbr_lab) == 0) {
        n_cells_zero_neighbors[s_lab] <- n_cells_zero_neighbors[s_lab] + 1L
        if (zero_neighbor == "drop") {
          next
        } else {
          n_cells_used[s_lab] <- n_cells_used[s_lab] + 1L
          next
        }
      }
      
      tab <- table(nbr_lab)
      frac <- as.numeric(tab) / sum(tab)
      cols <- neigh_to_col[names(tab)]
      r <- subj_to_row[[s_lab]]
      
      mat_sum[r, cols] <- mat_sum[r, cols] + frac
      n_cells_used[s_lab] <- n_cells_used[s_lab] + 1L
    }
    
    mat_mean <- mat_sum
    for (s in all_subject_levels) {
      denom <- n_cells_used[s]
      if (denom > 0) {
        mat_mean[s, ] <- mat_sum[s, ] / denom
      } else {
        mat_mean[s, ] <- NA_real_
      }
    }
    
    list(
      mat = mat_mean,
      stats = list(
        n_nodes = n_nodes,
        n_edges = n_edges,
        n_cells_total_by_subject = n_cells_total,
        n_cells_used_by_subject = n_cells_used,
        n_cells_zero_neighbors_by_subject = n_cells_zero_neighbors,
        n_cells_na_subject = sum(is.na(subject_chr)),
        n_cells_na_neighbor = sum(is.na(neighbor_chr)),
        subject_levels_present = sort(unique(subject_chr[!is.na(subject_chr)])),
        neighbor_levels_present = sort(unique(neighbor_chr[!is.na(neighbor_chr)]))
      )
    )
  }
  
  # -------------------------
  # First pass: collect global subject/neighbor levels
  # -------------------------
  all_subject_levels <- character(0)
  all_neighbor_levels <- character(0)
  
  for (obj_name in names(eff_list)) {
    obj <- eff_list[[obj_name]]
    node_ids <- .resolve_node_ids(obj, obj_name)
    md <- .resolve_metadata(
      obj, obj_name, node_ids,
      subject_col = subject_col,
      neighbor_col = neighbor_col,
      replicate_col = replicate_col
    )
    
    subject_chr  <- as.character(md[[subject_col]])
    neighbor_chr <- as.character(md[[neighbor_col]])
    
    all_subject_levels  <- union(all_subject_levels, sort(unique(subject_chr[!is.na(subject_chr)])))
    all_neighbor_levels <- union(all_neighbor_levels, sort(unique(neighbor_chr[!is.na(neighbor_chr)])))
  }
  
  all_subject_levels  <- sort(all_subject_levels)
  all_neighbor_levels <- sort(all_neighbor_levels)
  
  if (length(all_subject_levels) == 0) stop("No non-NA values found in subject_col across input object(s).")
  if (length(all_neighbor_levels) == 0) stop("No non-NA values found in neighbor_col across input object(s).")
  
  # -------------------------
  # Initialize pooled accumulators
  # -------------------------
  mat_sum <- matrix(
    0,
    nrow = length(all_subject_levels),
    ncol = length(all_neighbor_levels),
    dimnames = list(all_subject_levels, all_neighbor_levels)
  )
  
  n_cells_used <- setNames(integer(length(all_subject_levels)), all_subject_levels)
  n_cells_total <- setNames(integer(length(all_subject_levels)), all_subject_levels)
  n_cells_zero_neighbors <- setNames(integer(length(all_subject_levels)), all_subject_levels)
  
  n_cells_na_subject <- 0L
  n_cells_na_neighbor <- 0L
  
  total_nodes <- 0L
  total_edges <- 0L
  
  object_stats <- vector("list", length(eff_list))
  names(object_stats) <- names(eff_list)
  
  by_replicate <- list()
  replicate_stats <- list()
  
  # -------------------------
  # Main pass over objects
  # -------------------------
  for (obj_name in names(eff_list)) {
    obj <- eff_list[[obj_name]]
    
    node_ids <- .resolve_node_ids(obj, obj_name)
    n_nodes <- length(node_ids)
    md <- .resolve_metadata(
      obj, obj_name, node_ids,
      subject_col = subject_col,
      neighbor_col = neighbor_col,
      replicate_col = replicate_col
    )
    ij <- .resolve_ij(obj, obj_name, edge_name, n_nodes)
    
    n_edges <- nrow(ij)
    total_nodes <- total_nodes + n_nodes
    total_edges <- total_edges + n_edges
    
    subject_chr  <- as.character(md[[subject_col]])
    neighbor_chr <- as.character(md[[neighbor_col]])
    
    # ---- pooled object contribution
    pooled_res <- .compute_one_matrix(
      ij = ij,
      subject_chr = subject_chr,
      neighbor_chr = neighbor_chr,
      all_subject_levels = all_subject_levels,
      all_neighbor_levels = all_neighbor_levels,
      zero_neighbor = zero_neighbor,
      include_self = include_self
    )
    
    for (s in all_subject_levels) {
      denom <- pooled_res$stats$n_cells_used_by_subject[s]
      if (denom > 0) {
        mat_sum[s, ] <- mat_sum[s, ] + pooled_res$mat[s, ] * denom
      }
    }
    
    n_cells_used <- n_cells_used + pooled_res$stats$n_cells_used_by_subject
    n_cells_total <- n_cells_total + pooled_res$stats$n_cells_total_by_subject
    n_cells_zero_neighbors <- n_cells_zero_neighbors + pooled_res$stats$n_cells_zero_neighbors_by_subject
    n_cells_na_subject <- n_cells_na_subject + pooled_res$stats$n_cells_na_subject
    n_cells_na_neighbor <- n_cells_na_neighbor + pooled_res$stats$n_cells_na_neighbor
    
    object_stats[[obj_name]] <- list(
      n_nodes = n_nodes,
      n_edges = n_edges,
      subject_levels_present = pooled_res$stats$subject_levels_present,
      neighbor_levels_present = pooled_res$stats$neighbor_levels_present,
      n_cells_total_by_subject = pooled_res$stats$n_cells_total_by_subject,
      n_cells_used_by_subject = pooled_res$stats$n_cells_used_by_subject,
      n_cells_zero_neighbors_by_subject = pooled_res$stats$n_cells_zero_neighbors_by_subject,
      n_cells_na_subject = pooled_res$stats$n_cells_na_subject,
      n_cells_na_neighbor = pooled_res$stats$n_cells_na_neighbor
    )
    
    # ---- replicate-resolved computation
    if (!is.null(replicate_col)) {
      rep_chr <- as.character(md[[replicate_col]])
      rep_vals <- sort(unique(rep_chr[!is.na(rep_chr)]))
      
      for (rep_val in rep_vals) {
        keep_idx <- which(rep_chr == rep_val)
        if (length(keep_idx) == 0) next
        
        ij_rep <- .subset_and_remap_ij(
          ij = ij,
          keep_idx = keep_idx,
          n_nodes_total = n_nodes
        )
        
        rep_res <- .compute_one_matrix(
          ij = ij_rep,
          subject_chr = subject_chr[keep_idx],
          neighbor_chr = neighbor_chr[keep_idx],
          all_subject_levels = all_subject_levels,
          all_neighbor_levels = all_neighbor_levels,
          zero_neighbor = zero_neighbor,
          include_self = include_self
        )
        
        rep_name <- paste(obj_name, rep_val, sep = "::")
        by_replicate[[rep_name]] <- rep_res$mat
        
        replicate_stats[[rep_name]] <- list(
          object_name = obj_name,
          replicate_id = rep_val,
          replicate_name = rep_name,
          n_nodes = rep_res$stats$n_nodes,
          n_edges = rep_res$stats$n_edges,
          subject_levels_present = rep_res$stats$subject_levels_present,
          neighbor_levels_present = rep_res$stats$neighbor_levels_present,
          n_cells_total_by_subject = rep_res$stats$n_cells_total_by_subject,
          n_cells_used_by_subject = rep_res$stats$n_cells_used_by_subject,
          n_cells_zero_neighbors_by_subject = rep_res$stats$n_cells_zero_neighbors_by_subject,
          n_cells_na_subject = rep_res$stats$n_cells_na_subject,
          n_cells_na_neighbor = rep_res$stats$n_cells_na_neighbor
        )
      }
    }
  }
  
  # -------------------------
  # Convert pooled sums to mean-of-fractions across all cells
  # -------------------------
  mat_mean <- mat_sum
  for (s in all_subject_levels) {
    denom <- n_cells_used[s]
    if (denom > 0) {
      mat_mean[s, ] <- mat_sum[s, ] / denom
    } else {
      mat_mean[s, ] <- NA_real_
    }
  }
  
  diag <- list(
    edge_name = edge_name,
    n_objects = length(eff_list),
    object_names = names(eff_list),
    n_nodes = total_nodes,
    n_edges = total_edges,
    subject_col = subject_col,
    neighbor_col = neighbor_col,
    replicate_col = replicate_col,
    zero_neighbor = zero_neighbor,
    include_self = include_self,
    subject_levels = all_subject_levels,
    neighbor_levels = all_neighbor_levels,
    n_cells_total_by_subject = n_cells_total,
    n_cells_used_by_subject = n_cells_used,
    n_cells_zero_neighbors_by_subject = n_cells_zero_neighbors,
    n_cells_na_subject = n_cells_na_subject,
    n_cells_na_neighbor = n_cells_na_neighbor,
    by_object = object_stats,
    by_replicate = replicate_stats
  )
  
  if (verbose) {
    msg <- c(
      sprintf(
        "Neighborhood composition computed (edge-derived, undirected, unweighted) from '%s'.",
        edge_name
      ),
      sprintf(
        "Objects: %d | Nodes (total): %d | Edges (raw total): %d | Subject groups: %d | Neighbor groups: %d",
        length(eff_list), total_nodes, total_edges, length(all_subject_levels), length(all_neighbor_levels)
      ),
      sprintf("include_self: %s", include_self),
      sprintf("zero_neighbor handling: %s", zero_neighbor),
      sprintf(
        "replicate mode: %s",
        if (is.null(replicate_col)) "off" else paste0("on (metadata column: ", replicate_col, ")")
      )
    )
    if (!is.null(replicate_col)) {
      msg <- c(msg, sprintf("Replicates detected: %d", length(by_replicate)))
    }
    message(paste(msg, collapse = "\n"))
  }
  
  out <- list(
    mat = mat_mean,
    diag = diag
  )
  
  if (!is.null(replicate_col)) {
    out$by_replicate <- by_replicate
  }
  
  return(out)
}

# Plotting Neighborhood Composition Matrix
plotNeighborhoodComp <- function(
    eff,
    subject_col,
    neighbor_col,
    edge_name = "CellToCellSpatial",
    zero_neighbor = "drop",
    include_self = TRUE,
    replicate_col = NULL,
    summary_mode = c("pooled_cells", "mean_replicates"),
    title = NULL,
    palette = "inferno",
    na_fill = "grey20",
    show_values = FALSE,
    value_digits = 2,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    row_order = NULL,
    col_order = NULL,
    cap = NULL,
    min_display = NULL,
    text_size = 10,
    axis_text_x_angle = 45,
    grid_color = "white",
    grid_linewidth = 0.15,
    verbose = TRUE
) {
  summary_mode <- match.arg(summary_mode)
  
  comp <- computeNeighborhoodCompositionMatrix(
    eff = eff,
    subject_col = subject_col,
    neighbor_col = neighbor_col,
    edge_name = edge_name,
    zero_neighbor = zero_neighbor,
    include_self = include_self,
    replicate_col = replicate_col,
    verbose = verbose
  )
  
  # -------------------------
  # Choose matrix to display
  # -------------------------
  if (summary_mode == "pooled_cells") {
    mat <- as.matrix(comp$mat)
  } else {
    if (is.null(replicate_col)) {
      stop("summary_mode = 'mean_replicates' requires replicate_col to be provided.")
    }
    if (is.null(comp$by_replicate) || length(comp$by_replicate) == 0) {
      stop("No replicate-level matrices were available to average.")
    }
    
    rep_arr <- simplify2array(comp$by_replicate)
    if (length(dim(rep_arr)) == 2) {
      mat <- rep_arr
    } else {
      mat <- apply(rep_arr, c(1, 2), function(x) {
        if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      })
    }
    
    rownames(mat) <- rownames(comp$mat)
    colnames(mat) <- colnames(comp$mat)
  }
  
  # -------------------------
  # Optional display modifications
  # -------------------------
  mat_disp <- mat
  
  if (!is.null(cap)) {
    mat_disp <- pmin(mat_disp, cap)
  }
  if (!is.null(min_display)) {
    mat_disp[!is.na(mat_disp) & mat_disp < min_display] <- min_display
  }
  
  # -------------------------
  # Determine dynamic max for scaling
  # -------------------------
  max_val <- suppressWarnings(max(mat_disp, na.rm = TRUE))
  if (!is.finite(max_val) || max_val == 0) {
    max_val <- 1
  }
  
  # -------------------------
  # Ordering / clustering
  # -------------------------
  if (!is.null(row_order)) {
    mat_disp <- mat_disp[row_order, , drop = FALSE]
  } else if (isTRUE(cluster_rows) && nrow(mat_disp) > 1) {
    rr <- hclust(dist(mat_disp))
    mat_disp <- mat_disp[rr$order, , drop = FALSE]
  }
  
  if (!is.null(col_order)) {
    mat_disp <- mat_disp[, col_order, drop = FALSE]
  } else if (isTRUE(cluster_cols) && ncol(mat_disp) > 1) {
    cc <- hclust(dist(t(mat_disp)))
    mat_disp <- mat_disp[, cc$order, drop = FALSE]
  }
  
  # Reverse row order so first alphabetical appears at top
  mat_disp <- mat_disp[rev(rownames(mat_disp)), , drop = FALSE]
  
  # -------------------------
  # Long format
  # -------------------------
  df <- as.data.frame(as.table(mat_disp), stringsAsFactors = FALSE)
  colnames(df) <- c("subject", "neighbor", "fraction")
  
  df$subject  <- factor(df$subject, levels = rownames(mat_disp))
  df$neighbor <- factor(df$neighbor, levels = colnames(mat_disp))
  
  if (is.null(title)) {
    title <- sprintf("Neighborhood composition: %s vs %s", subject_col, neighbor_col)
    if (!is.null(replicate_col) && summary_mode == "mean_replicates") {
      title <- paste0(title, " (mean across replicates)")
    }
  }
  
  # -------------------------
  # Plot
  # -------------------------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = neighbor, y = subject, fill = fraction)) +
    ggplot2::geom_tile(
      color = grid_color,
      linewidth = grid_linewidth
    )
  
  n_rows <- nrow(mat_disp)
  
  p <- p +
    ggplot2::geom_hline(
      yintercept = seq(0.5, n_rows + 0.5, by = 1),
      color = grid_color,
      linewidth = grid_linewidth * 3
    )
  
  p <- p +
    viridis::scale_fill_viridis(
      option = palette,
      na.value = na_fill,
      limits = c(0, max_val),
      name = "Fraction"
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(color = "black", face = "bold"),
      axis.text.x.top = ggplot2::element_text(
        color = "black",
        angle = axis_text_x_angle,
        hjust = 0,
        vjust = 0
      ),
      axis.text.x.bottom = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      panel.grid = ggplot2::element_blank()
    )
  
  # -------------------------
  # Optional numeric overlay
  # -------------------------
  if (isTRUE(show_values)) {
    df_raw <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df_raw) <- c("subject", "neighbor", "fraction_raw")
    
    df_raw$subject  <- factor(df_raw$subject, levels = levels(df$subject))
    df_raw$neighbor <- factor(df_raw$neighbor, levels = levels(df$neighbor))
    
    df_plot <- merge(df, df_raw, by = c("subject", "neighbor"), all.x = TRUE)
    
    p <- p + ggplot2::geom_text(
      data = df_plot,
      ggplot2::aes(label = ifelse(
        is.na(fraction_raw),
        "",
        sprintf(paste0("%.", value_digits, "f"), fraction_raw)
      )),
      color = "black",
      size = (text_size / 3.2)
    )
  }
  
  return(p)
}

#### Comparitive Neighborhood Comp Plot Dev ####

# Compare neighborhood composition between two biological groups at the replicate level.
#
# This function performs statistically appropriate group comparisons for
# neighborhood composition matrices derived from Xenium / spatial transcriptomic
# EffNICHES objects. To avoid pseudoreplication, the unit of replication is the
# biological replicate specified by `replicate_col` (for example, a TMA core),
# not the individual cell.
#
# Supported input modes:
# 1. A single EffNICHES object or list of EffNICHES objects in `eff`, where
#    group identity is defined in node metadata via `group_col`.
# 2. Two separate group inputs, `eff_group1` and `eff_group2`, each of which
#    can be a single EffNICHES object or a list of EffNICHES objects.
#
# Workflow:
# - Computes one neighborhood composition matrix per replicate using
#   `computeNeighborhoodCompositionMatrix()`.
# - Assigns each replicate to one of two groups.
# - For each subject-neighbor matrix cell, compares replicate-level values
#   between groups using either Welch's t-test or Wilcoxon rank-sum test.
# - Adjusts p-values across all tested cells using the specified multiple-testing
#   correction method (default: Benjamini-Hochberg FDR).
#
# Output includes:
# - Mean matrix for group 1
# - Mean matrix for group 2
# - Difference matrix (group 2 - group 1)
# - Raw p-value matrix
# - FDR-adjusted p-value matrix
# - Test statistic matrix
# - Replicate count matrices per group
# - Long-format results table
# - Long-format replicate-level values table
# - Replicate-level matrices and diagnostic metadata
#
# This function is intended as the statistical comparison backend for downstream
# neighborhood composition comparison plots and summaries.
compareNeighborhoodComp <- function(
    eff = NULL,
    eff_group1 = NULL,
    eff_group2 = NULL,
    subject_col,
    neighbor_col,
    group_col = NULL,
    replicate_col,
    edge_name = "CellToCellSpatial",
    zero_neighbor = c("drop", "zero_row"),
    include_self = TRUE,
    test = c("welch", "wilcox"),
    p_adjust_method = "BH",
    group_order = NULL,
    min_reps_per_group = 2,
    verbose = TRUE
) {
  zero_neighbor <- match.arg(zero_neighbor)
  test <- match.arg(test)
  
  # -------------------------
  # Helpers
  # -------------------------
  .normalize_eff_list <- function(x, arg_name) {
    if (inherits(x, "EffNICHES")) {
      out <- list(x)
      names(out) <- arg_name
      return(out)
    }
    
    if (is.list(x) && length(x) > 0 &&
        all(vapply(x, inherits, logical(1), "EffNICHES"))) {
      out <- x
      if (is.null(names(out))) {
        names(out) <- paste0(arg_name, "_", seq_along(out))
      } else {
        blank_names <- is.na(names(out)) | names(out) == ""
        names(out)[blank_names] <- paste0(arg_name, "_", which(blank_names))
      }
      return(out)
    }
    
    stop(sprintf(
      "%s must be either a single EffNICHES object or a non-empty list of EffNICHES objects.",
      arg_name
    ))
  }
  
  .prefix_list_names <- function(x, prefix) {
    names(x) <- paste0(prefix, names(x))
    x
  }
  
  .resolve_node_ids <- function(obj, obj_name) {
    node_ids <- NULL
    if (!is.null(obj$nodes$ids)) {
      node_ids <- as.character(obj$nodes$ids)
    } else if (!is.null(obj$nodes$xy) && !is.null(rownames(obj$nodes$xy))) {
      node_ids <- rownames(obj$nodes$xy)
    }
    if (is.null(node_ids)) {
      stop(sprintf("[%s] Could not resolve node IDs. Need eff$nodes$ids or rownames(eff$nodes$xy).", obj_name))
    }
    node_ids
  }
  
  .resolve_metadata <- function(obj, obj_name, node_ids, cols_needed) {
    md <- obj$nodes$meta
    if (is.null(md) || !is.data.frame(md)) {
      stop(sprintf("[%s] eff$nodes$meta must be a data.frame.", obj_name))
    }
    if (is.null(rownames(md))) {
      stop(sprintf("[%s] eff$nodes$meta must have rownames matching node IDs for safe alignment.", obj_name))
    }
    
    missing_md <- setdiff(node_ids, rownames(md))
    if (length(missing_md) > 0) {
      stop(sprintf(
        "[%s] Metadata rownames do not cover all node IDs. Missing %d IDs (e.g. %s).",
        obj_name, length(missing_md), paste(head(missing_md, 5), collapse = ", ")
      ))
    }
    
    md <- md[node_ids, , drop = FALSE]
    
    missing_cols <- setdiff(cols_needed, colnames(md))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "[%s] Missing required metadata column(s): %s",
        obj_name, paste(missing_cols, collapse = ", ")
      ))
    }
    
    md
  }
  
  .safe_test <- function(x, y, test) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    
    out <- list(
      p_value = NA_real_,
      statistic = NA_real_,
      estimate_group1 = if (length(x) > 0) mean(x) else NA_real_,
      estimate_group2 = if (length(y) > 0) mean(y) else NA_real_
    )
    
    if (length(x) == 0 || length(y) == 0) return(out)
    
    if (test == "welch") {
      tt <- try(stats::t.test(x, y, var.equal = FALSE), silent = TRUE)
      if (!inherits(tt, "try-error")) {
        out$p_value <- unname(tt$p.value)
        out$statistic <- unname(tt$statistic)
      }
    } else if (test == "wilcox") {
      wt <- try(stats::wilcox.test(x, y, exact = FALSE), silent = TRUE)
      if (!inherits(wt, "try-error")) {
        out$p_value <- unname(wt$p.value)
        out$statistic <- unname(wt$statistic)
      }
    }
    
    out
  }
  
  .mean_matrix_list <- function(mat_list, template_mat) {
    if (length(mat_list) == 0) {
      out <- template_mat
      out[,] <- NA_real_
      return(out)
    }
    
    arr <- simplify2array(mat_list)
    
    if (length(dim(arr)) == 2) {
      out <- arr
    } else {
      out <- apply(arr, c(1, 2), function(x) {
        if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      })
    }
    
    rownames(out) <- rownames(template_mat)
    colnames(out) <- colnames(template_mat)
    out
  }
  
  # -------------------------
  # Determine comparison mode
  # -------------------------
  using_internal_groups <- !is.null(eff)
  using_separate_groups <- !is.null(eff_group1) || !is.null(eff_group2)
  
  if (using_internal_groups && using_separate_groups) {
    stop("Use either eff with group_col, or eff_group1/eff_group2. Do not mix both modes.")
  }
  if (!using_internal_groups && !using_separate_groups) {
    stop("Provide either eff with group_col, or both eff_group1 and eff_group2.")
  }
  
  comparison_mode <- NULL
  combined_eff <- NULL
  rep_group_df <- NULL
  
  # -------------------------
  # Mode 1: groups defined in metadata
  # -------------------------
  if (using_internal_groups) {
    comparison_mode <- "metadata_groups"
    
    if (is.null(group_col)) {
      stop("group_col must be provided when using eff.")
    }
    
    eff_list <- .normalize_eff_list(eff, "eff")
    combined_eff <- eff_list
    
    # ---- Compute replicate-level matrices
    comp <- computeNeighborhoodCompositionMatrix(
      eff = combined_eff,
      subject_col = subject_col,
      neighbor_col = neighbor_col,
      edge_name = edge_name,
      zero_neighbor = zero_neighbor,
      include_self = include_self,
      replicate_col = replicate_col,
      verbose = verbose
    )
    
    if (is.null(comp$by_replicate) || length(comp$by_replicate) == 0) {
      stop("No replicate-level matrices were returned. Check replicate_col and input data.")
    }
    
    # ---- Build replicate -> group map from metadata
    rep_group_list <- list()
    
    for (obj_name in names(combined_eff)) {
      obj <- combined_eff[[obj_name]]
      node_ids <- .resolve_node_ids(obj, obj_name)
      md <- .resolve_metadata(
        obj = obj,
        obj_name = obj_name,
        node_ids = node_ids,
        cols_needed = c(subject_col, neighbor_col, group_col, replicate_col)
      )
      
      rep_chr <- as.character(md[[replicate_col]])
      grp_chr <- as.character(md[[group_col]])
      rep_vals <- sort(unique(rep_chr[!is.na(rep_chr)]))
      
      for (rep_val in rep_vals) {
        idx <- which(rep_chr == rep_val)
        grp_vals <- sort(unique(grp_chr[idx][!is.na(grp_chr[idx])]))
        rep_name <- paste(obj_name, rep_val, sep = "::")
        
        if (length(grp_vals) == 0) {
          stop(sprintf(
            "[%s] Replicate '%s' has no non-NA group values in '%s'.",
            obj_name, rep_val, group_col
          ))
        }
        if (length(grp_vals) > 1) {
          stop(sprintf(
            "[%s] Replicate '%s' maps to multiple groups in '%s': %s",
            obj_name, rep_val, group_col, paste(grp_vals, collapse = ", ")
          ))
        }
        
        rep_group_list[[rep_name]] <- data.frame(
          object_name = obj_name,
          replicate_id = rep_val,
          replicate_name = rep_name,
          group = grp_vals,
          stringsAsFactors = FALSE
        )
      }
    }
    
    rep_group_df <- do.call(rbind, rep_group_list)
    rownames(rep_group_df) <- rep_group_df$replicate_name
    
    missing_map <- setdiff(names(comp$by_replicate), rownames(rep_group_df))
    if (length(missing_map) > 0) {
      stop(sprintf(
        "Could not determine group assignments for %d replicate(s): %s",
        length(missing_map), paste(head(missing_map, 5), collapse = ", ")
      ))
    }
    
    rep_group_df <- rep_group_df[names(comp$by_replicate), , drop = FALSE]
  }
  
  # -------------------------
  # Mode 2: two separate group objects
  # -------------------------
  if (using_separate_groups) {
    comparison_mode <- "separate_objects"
    
    if (is.null(eff_group1) || is.null(eff_group2)) {
      stop("When using separate group mode, both eff_group1 and eff_group2 must be provided.")
    }
    
    eff_list_g1 <- .normalize_eff_list(eff_group1, "group1")
    eff_list_g2 <- .normalize_eff_list(eff_group2, "group2")
    
    eff_list_g1 <- .prefix_list_names(eff_list_g1, "g1__")
    eff_list_g2 <- .prefix_list_names(eff_list_g2, "g2__")
    
    combined_eff <- c(eff_list_g1, eff_list_g2)
    
    # ---- Group names
    if (is.null(group_order)) {
      group_order <- c("group1", "group2")
    } else {
      if (length(group_order) != 2) {
        stop("group_order must contain exactly 2 group names.")
      }
    }
    
    g1_name_tmp <- group_order[1]
    g2_name_tmp <- group_order[2]
    
    # ---- Compute replicate-level matrices on combined object list
    comp <- computeNeighborhoodCompositionMatrix(
      eff = combined_eff,
      subject_col = subject_col,
      neighbor_col = neighbor_col,
      edge_name = edge_name,
      zero_neighbor = zero_neighbor,
      include_self = include_self,
      replicate_col = replicate_col,
      verbose = verbose
    )
    
    if (is.null(comp$by_replicate) || length(comp$by_replicate) == 0) {
      stop("No replicate-level matrices were returned. Check replicate_col and input data.")
    }
    
    # ---- Build replicate -> group map from source object side
    rep_group_list <- list()
    
    for (obj_name in names(combined_eff)) {
      obj <- combined_eff[[obj_name]]
      node_ids <- .resolve_node_ids(obj, obj_name)
      md <- .resolve_metadata(
        obj = obj,
        obj_name = obj_name,
        node_ids = node_ids,
        cols_needed = c(subject_col, neighbor_col, replicate_col)
      )
      
      rep_chr <- as.character(md[[replicate_col]])
      rep_vals <- sort(unique(rep_chr[!is.na(rep_chr)]))
      
      this_group <- if (startsWith(obj_name, "g1__")) g1_name_tmp else g2_name_tmp
      
      for (rep_val in rep_vals) {
        rep_name <- paste(obj_name, rep_val, sep = "::")
        
        rep_group_list[[rep_name]] <- data.frame(
          object_name = obj_name,
          replicate_id = rep_val,
          replicate_name = rep_name,
          group = this_group,
          stringsAsFactors = FALSE
        )
      }
    }
    
    rep_group_df <- do.call(rbind, rep_group_list)
    rownames(rep_group_df) <- rep_group_df$replicate_name
    
    missing_map <- setdiff(names(comp$by_replicate), rownames(rep_group_df))
    if (length(missing_map) > 0) {
      stop(sprintf(
        "Could not determine group assignments for %d replicate(s): %s",
        length(missing_map), paste(head(missing_map, 5), collapse = ", ")
      ))
    }
    
    rep_group_df <- rep_group_df[names(comp$by_replicate), , drop = FALSE]
  }
  
  # -------------------------
  # Resolve group order
  # -------------------------
  groups_present <- unique(rep_group_df$group)
  
  if (is.null(group_order)) {
    if (length(groups_present) != 2) {
      stop(sprintf(
        "Exactly 2 groups are required unless group_order is provided. Found %d group(s): %s",
        length(groups_present), paste(sort(groups_present), collapse = ", ")
      ))
    }
    group_order <- sort(groups_present)
  } else {
    if (length(group_order) != 2) {
      stop("group_order must contain exactly 2 group names.")
    }
    if (!all(group_order %in% groups_present)) {
      stop(sprintf(
        "Not all group_order values were found in the data. Missing: %s",
        paste(setdiff(group_order, groups_present), collapse = ", ")
      ))
    }
  }
  
  g1 <- group_order[1]
  g2 <- group_order[2]
  
  # -------------------------
  # Build long replicate-level table
  # -------------------------
  long_list <- vector("list", length(comp$by_replicate))
  rep_names <- names(comp$by_replicate)
  
  for (i in seq_along(rep_names)) {
    rep_name <- rep_names[i]
    mat <- comp$by_replicate[[rep_name]]
    
    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df) <- c("subject", "neighbor", "fraction")
    
    df$replicate_name <- rep_name
    df$object_name <- rep_group_df[rep_name, "object_name"]
    df$replicate_id <- rep_group_df[rep_name, "replicate_id"]
    df$group <- rep_group_df[rep_name, "group"]
    
    long_list[[i]] <- df[, c(
      "object_name", "replicate_id", "replicate_name",
      "group", "subject", "neighbor", "fraction"
    )]
  }
  
  long_df <- do.call(rbind, long_list)
  rownames(long_df) <- NULL
  
  # -------------------------
  # Group replicate counts
  # -------------------------
  rep_counts <- table(rep_group_df$group)
  
  if (verbose) {
    message(
      paste0(
        "Comparing neighborhood composition across groups.\n",
        "Mode: ", comparison_mode, "\n",
        "Group 1: ", g1, " (n = ", unname(rep_counts[g1]), " replicate(s))\n",
        "Group 2: ", g2, " (n = ", unname(rep_counts[g2]), " replicate(s))\n",
        "Test: ", test
      )
    )
  }
  
  # -------------------------
  # Per-cell testing
  # -------------------------
  subject_levels <- rownames(comp$mat)
  neighbor_levels <- colnames(comp$mat)
  
  mean_g1 <- matrix(
    NA_real_,
    nrow = length(subject_levels),
    ncol = length(neighbor_levels),
    dimnames = list(subject_levels, neighbor_levels)
  )
  mean_g2 <- mean_g1
  diff_mat <- mean_g1
  p_mat <- mean_g1
  fdr_mat <- mean_g1
  n_g1_mat <- mean_g1
  n_g2_mat <- mean_g1
  stat_mat <- mean_g1
  
  results_long <- vector("list", length(subject_levels) * length(neighbor_levels))
  k <- 0L
  
  for (s in subject_levels) {
    for (n in neighbor_levels) {
      k <- k + 1L
      
      sub_df <- long_df[long_df$subject == s & long_df$neighbor == n, , drop = FALSE]
      
      x <- sub_df$fraction[sub_df$group == g1]
      y <- sub_df$fraction[sub_df$group == g2]
      
      n1 <- sum(!is.na(x))
      n2 <- sum(!is.na(y))
      
      m1 <- if (n1 > 0) mean(x, na.rm = TRUE) else NA_real_
      m2 <- if (n2 > 0) mean(y, na.rm = TRUE) else NA_real_
      d  <- m2 - m1
      
      mean_g1[s, n] <- m1
      mean_g2[s, n] <- m2
      diff_mat[s, n] <- d
      n_g1_mat[s, n] <- n1
      n_g2_mat[s, n] <- n2
      
      pval <- NA_real_
      stat <- NA_real_
      
      if (n1 >= min_reps_per_group && n2 >= min_reps_per_group) {
        tst <- .safe_test(x, y, test = test)
        pval <- tst$p_value
        stat <- tst$statistic
      }
      
      p_mat[s, n] <- pval
      stat_mat[s, n] <- stat
      
      results_long[[k]] <- data.frame(
        subject = s,
        neighbor = n,
        group1 = g1,
        group2 = g2,
        mean_group1 = m1,
        mean_group2 = m2,
        difference = d,
        n_group1 = n1,
        n_group2 = n2,
        statistic = stat,
        p_value = pval,
        stringsAsFactors = FALSE
      )
    }
  }
  
  results_long <- do.call(rbind, results_long)
  
  # -------------------------
  # Multiple-testing correction
  # -------------------------
  valid_p <- !is.na(results_long$p_value)
  results_long$fdr <- NA_real_
  if (any(valid_p)) {
    results_long$fdr[valid_p] <- stats::p.adjust(
      results_long$p_value[valid_p],
      method = p_adjust_method
    )
  }
  
  for (i in seq_len(nrow(results_long))) {
    s <- results_long$subject[i]
    n <- results_long$neighbor[i]
    fdr_mat[s, n] <- results_long$fdr[i]
  }
  
  # -------------------------
  # Replicate-level mean matrices by group
  # -------------------------
  reps_g1 <- rep_group_df$replicate_name[rep_group_df$group == g1]
  reps_g2 <- rep_group_df$replicate_name[rep_group_df$group == g2]
  
  mats_g1 <- comp$by_replicate[reps_g1]
  mats_g2 <- comp$by_replicate[reps_g2]
  
  mean_mat_g1 <- .mean_matrix_list(mats_g1, comp$mat)
  mean_mat_g2 <- .mean_matrix_list(mats_g2, comp$mat)
  
  # -------------------------
  # Output
  # -------------------------
  diag <- list(
    comparison_mode = comparison_mode,
    subject_col = subject_col,
    neighbor_col = neighbor_col,
    group_col = group_col,
    replicate_col = replicate_col,
    edge_name = edge_name,
    zero_neighbor = zero_neighbor,
    include_self = include_self,
    test = test,
    p_adjust_method = p_adjust_method,
    min_reps_per_group = min_reps_per_group,
    group_order = group_order,
    n_groups = 2L,
    replicate_counts_by_group = as.list(rep_counts[group_order]),
    replicate_map = rep_group_df,
    comparison = sprintf("%s - %s", g2, g1)
  )
  
  out <- list(
    mean_group1 = mean_mat_g1,
    mean_group2 = mean_mat_g2,
    difference = diff_mat,
    p_value = p_mat,
    fdr = fdr_mat,
    statistic = stat_mat,
    n_group1 = n_g1_mat,
    n_group2 = n_g2_mat,
    results = results_long,
    replicate_values = long_df,
    by_replicate = comp$by_replicate,
    diag = diag
  )
  
  return(out)
}

# Plot comparative neighborhood composition heatmap.
#
# This function is a wrapper around `compareNeighborhoodComp()`. It computes
# replicate-level group differences in neighborhood composition and displays the
# resulting difference matrix as a heatmap in the same overall format as
# `plotNeighborhoodComp()`.
#
# Visual encoding:
# - White = no change
# - Red = increased neighborhood composition in group 2 relative to group 1
# - Blue = decreased neighborhood composition in group 2 relative to group 1
#
# Statistical significance:
# - Cells with significance at or below `sig_threshold` are marked with a black
#   asterisk, using either FDR or raw p-values.
plotNeighborhoodCompComparison <- function(
    eff = NULL,
    eff_group1 = NULL,
    eff_group2 = NULL,
    subject_col,
    neighbor_col,
    group_col = NULL,
    replicate_col,
    edge_name = "CellToCellSpatial",
    zero_neighbor = "drop",
    include_self = TRUE,
    test = "welch",
    p_adjust_method = "BH",
    group_order = NULL,
    min_reps_per_group = 2,
    title = NULL,
    low_color = "blue",
    mid_color = "white",
    high_color = "red",
    midpoint = 0,
    na_fill = "grey20",
    show_values = FALSE,
    value_digits = 2,
    show_sig = TRUE,
    sig_metric = c("fdr", "p_value"),
    sig_threshold = 0.05,
    sig_symbol = "*",
    sig_text_size = 5,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    row_order = NULL,
    col_order = NULL,
    cap = NULL,
    min_display = NULL,
    symmetric_limits = TRUE,
    fill_limits = NULL,
    legend_title = "Difference",
    oob = scales::squish,
    text_size = 10,
    axis_text_x_angle = 45,
    grid_color = "white",
    grid_linewidth = 0.15,
    verbose = TRUE
) {
  sig_metric <- match.arg(sig_metric)
  
  # -------------------------
  # Compute comparison
  # -------------------------
  cmp <- compareNeighborhoodComp(
    eff = eff,
    eff_group1 = eff_group1,
    eff_group2 = eff_group2,
    subject_col = subject_col,
    neighbor_col = neighbor_col,
    group_col = group_col,
    replicate_col = replicate_col,
    edge_name = edge_name,
    zero_neighbor = zero_neighbor,
    include_self = include_self,
    test = test,
    p_adjust_method = p_adjust_method,
    group_order = group_order,
    min_reps_per_group = min_reps_per_group,
    verbose = verbose
  )
  
  diff_mat <- as.matrix(cmp$difference)
  sig_mat  <- as.matrix(cmp[[sig_metric]])
  
  if (!identical(dim(diff_mat), dim(sig_mat))) {
    stop("difference matrix and significance matrix must have identical dimensions.")
  }
  if (!identical(rownames(diff_mat), rownames(sig_mat)) ||
      !identical(colnames(diff_mat), colnames(sig_mat))) {
    stop("difference matrix and significance matrix must have matching row/column names.")
  }
  
  # -------------------------
  # Optional display modifications
  # -------------------------
  mat_disp <- diff_mat
  
  if (!is.null(cap)) {
    mat_disp[!is.na(mat_disp)] <- pmax(pmin(mat_disp[!is.na(mat_disp)], cap), -cap)
  }
  
  if (!is.null(min_display)) {
    idx_small_pos <- !is.na(mat_disp) & mat_disp > 0 & mat_disp < min_display
    idx_small_neg <- !is.na(mat_disp) & mat_disp < 0 & mat_disp > -min_display
    
    mat_disp[idx_small_pos] <- min_display
    mat_disp[idx_small_neg] <- -min_display
  }
  
  # -------------------------
  # Ordering / clustering
  # -------------------------
  if (!is.null(row_order)) {
    mat_disp <- mat_disp[row_order, , drop = FALSE]
    sig_mat  <- sig_mat[row_order, , drop = FALSE]
  } else if (isTRUE(cluster_rows) && nrow(mat_disp) > 1) {
    rr <- stats::hclust(stats::dist(mat_disp))
    mat_disp <- mat_disp[rr$order, , drop = FALSE]
    sig_mat  <- sig_mat[rr$order, , drop = FALSE]
  }
  
  if (!is.null(col_order)) {
    mat_disp <- mat_disp[, col_order, drop = FALSE]
    sig_mat  <- sig_mat[, col_order, drop = FALSE]
  } else if (isTRUE(cluster_cols) && ncol(mat_disp) > 1) {
    cc <- stats::hclust(stats::dist(t(mat_disp)))
    mat_disp <- mat_disp[, cc$order, drop = FALSE]
    sig_mat  <- sig_mat[, cc$order, drop = FALSE]
  }
  
  # Reverse row order so first alphabetical appears at top
  mat_disp <- mat_disp[rev(rownames(mat_disp)), , drop = FALSE]
  sig_mat  <- sig_mat[rownames(mat_disp), colnames(mat_disp), drop = FALSE]
  
  # -------------------------
  # Determine color limits
  # -------------------------
  if (is.null(fill_limits)) {
    max_abs <- suppressWarnings(max(abs(mat_disp), na.rm = TRUE))
    if (!is.finite(max_abs) || max_abs == 0) {
      max_abs <- 1
    }
    
    if (isTRUE(symmetric_limits)) {
      fill_limits_used <- c(-max_abs, max_abs)
    } else {
      fill_limits_used <- range(mat_disp, na.rm = TRUE)
      if (!all(is.finite(fill_limits_used)) || diff(fill_limits_used) == 0) {
        fill_limits_used <- c(-1, 1)
      }
    }
  } else {
    if (!is.numeric(fill_limits) || length(fill_limits) != 2 || any(!is.finite(fill_limits))) {
      stop("'fill_limits' must be a numeric vector of length 2 with finite values.")
    }
    if (fill_limits[1] >= fill_limits[2]) {
      stop("'fill_limits' must satisfy fill_limits[1] < fill_limits[2].")
    }
    fill_limits_used <- fill_limits
  }
  
  # -------------------------
  # Long format
  # -------------------------
  df <- as.data.frame(as.table(mat_disp), stringsAsFactors = FALSE)
  colnames(df) <- c("subject", "neighbor", "difference")
  
  df_sig <- as.data.frame(as.table(sig_mat), stringsAsFactors = FALSE)
  colnames(df_sig) <- c("subject", "neighbor", "sig_value")
  
  df <- merge(df, df_sig, by = c("subject", "neighbor"), all.x = TRUE)
  
  df$subject  <- factor(df$subject, levels = rownames(mat_disp))
  df$neighbor <- factor(df$neighbor, levels = colnames(mat_disp))
  
  df$is_significant <- !is.na(df$sig_value) & df$sig_value <= sig_threshold
  df$sig_label <- ifelse(df$is_significant & show_sig, sig_symbol, "")
  
  if (is.null(title)) {
    comp_label <- NULL
    if (!is.null(cmp$diag$comparison)) {
      comp_label <- cmp$diag$comparison
    } else if (!is.null(cmp$diag$group_order) && length(cmp$diag$group_order) == 2) {
      comp_label <- paste0(cmp$diag$group_order[2], " - ", cmp$diag$group_order[1])
    }
    
    if (is.null(comp_label)) {
      title <- "Neighborhood composition difference"
    } else {
      title <- paste0("Neighborhood composition difference: ", comp_label)
    }
  }
  
  # -------------------------
  # Plot
  # -------------------------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = neighbor, y = subject, fill = difference)) +
    ggplot2::geom_tile(
      color = grid_color,
      linewidth = grid_linewidth
    )
  
  n_rows <- nrow(mat_disp)
  
  p <- p +
    ggplot2::geom_hline(
      yintercept = seq(0.5, n_rows + 0.5, by = 1),
      color = grid_color,
      linewidth = grid_linewidth * 3
    ) +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = midpoint,
      limits = fill_limits_used,
      na.value = na_fill,
      oob = oob,
      name = legend_title
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(color = "black", face = "bold"),
      axis.text.x.top = ggplot2::element_text(
        color = "black",
        angle = axis_text_x_angle,
        hjust = 0,
        vjust = 0
      ),
      axis.text.x.bottom = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      panel.grid = ggplot2::element_blank()
    )
  
  # -------------------------
  # Significance overlay
  # -------------------------
  if (isTRUE(show_sig)) {
    p <- p + ggplot2::geom_text(
      data = df[df$is_significant, , drop = FALSE],
      ggplot2::aes(label = sig_label),
      color = "black",
      size = sig_text_size
    )
  }
  
  # -------------------------
  # Optional numeric overlay
  # -------------------------
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(
      data = df,
      ggplot2::aes(label = ifelse(
        is.na(difference),
        "",
        sprintf(paste0("%.", value_digits, "f"), difference)
      )),
      color = "black",
      size = (text_size / 3.2),
      vjust = if (show_sig) 1.6 else 0.5
    )
  }
  
  attr(p, "fill_limits_used") <- fill_limits_used
  attr(p, "comparison_result") <- cmp
  
  return(p)
}
