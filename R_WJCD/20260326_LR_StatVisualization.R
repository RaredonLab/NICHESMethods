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

#### NICHES Data Distribution Exploration ####

plot_eff_mechanism_distribution <- function(
    eff,
    mechanism,
    include_zeros = TRUE,
    bins = 50
) {
  # Pull CellToCellSpatial weight matrix
  w <- eff$niches$CellToCellSpatial$w
  mech_names <- colnames(w)
  
  # Normalize dash variants for matching
  mech_names_norm <- norm_lr_dash(mech_names)
  mechanism_norm  <- norm_lr_dash(mechanism)
  
  # Match requested mechanism
  hit_idx <- which(mech_names_norm == mechanism_norm)
  
  if (length(hit_idx) == 0) {
    stop("Mechanism not found in eff$niches$CellToCellSpatial$w")
  }
  if (length(hit_idx) > 1) {
    stop("Multiple mechanism matches found after dash normalization")
  }
  
  mechanism_matched <- mech_names[hit_idx]
  
  # Extract values
  x <- as.numeric(w[, hit_idx])
  
  # Optionally restrict to realized edges only
  if (!include_zeros) {
    x <- x[x != 0]
  }
  
  # Basic summaries
  cat("Requested mechanism:", mechanism, "\n")
  cat("Matched mechanism:", mechanism_matched, "\n")
  cat("Length:", length(x), "\n")
  cat("Zeros:", sum(x == 0, na.rm = TRUE), "\n")
  cat("Nonzero:", sum(x != 0, na.rm = TRUE), "\n")
  cat("NA:", sum(is.na(x)), "\n\n")
  print(summary(x))
  
  # Plot
  hist(
    x,
    breaks = bins,
    main = paste0(
      mechanism_matched,
      if (include_zeros) " (all edges)" else " (nonzero edges only)"
    ),
    xlab = "NICHES value",
    col = "grey70",
    border = "white"
  )
  
  invisible(x)
}

plot_eff_mechanism_distribution_two <- function(
    eff1,
    eff2,
    mechanism,
    include_zeros = TRUE,
    bins = 50,
    label1 = "Object 1",
    label2 = "Object 2"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  
  # Internal helper to extract one mechanism vector
  .get_mech_values <- function(eff, mechanism) {
    w <- eff$niches$CellToCellSpatial$w
    mech_names <- colnames(w)
    
    mech_names_norm <- norm_lr_dash(mech_names)
    mechanism_norm  <- norm_lr_dash(mechanism)
    
    hit_idx <- which(mech_names_norm == mechanism_norm)
    
    if (length(hit_idx) == 0) {
      stop("Mechanism not found in eff$niches$CellToCellSpatial$w")
    }
    if (length(hit_idx) > 1) {
      stop("Multiple mechanism matches found after dash normalization")
    }
    
    list(
      x = as.numeric(w[, hit_idx]),
      matched_name = mech_names[hit_idx],
      total_edges = nrow(w)
    )
  }
  
  # Pull data
  d1 <- .get_mech_values(eff1, mechanism)
  d2 <- .get_mech_values(eff2, mechanism)
  
  total_edges1 <- d1$total_edges
  total_edges2 <- d2$total_edges
  
  x1 <- d1$x
  x2 <- d2$x
  
  if (!include_zeros) {
    x1 <- x1[x1 != 0]
    x2 <- x2[x2 != 0]
  }
  
  combined_x <- c(x1, x2)
  combined_x <- combined_x[!is.na(combined_x)]
  
  if (length(combined_x) == 0) {
    stop("No non-NA values available for plotting.")
  }
  
  x_range <- range(combined_x)
  
  if (diff(x_range) == 0) {
    x_range <- x_range + c(-0.5, 0.5)
  }
  
  breaks <- seq(x_range[1], x_range[2], length.out = bins + 1)
  
  h1 <- hist(x1, breaks = breaks, plot = FALSE)
  h2 <- hist(x2, breaks = breaks, plot = FALSE)
  
  y1 <- h1$counts / total_edges1
  y2 <- h2$counts / total_edges2
  
  ymax <- max(c(y1, y2), na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  
  # Build rectangle-based data (true bin edges)
  df1 <- data.frame(
    xmin = breaks[-length(breaks)],
    xmax = breaks[-1],
    proportion = y1
  )
  
  df2 <- data.frame(
    xmin = breaks[-length(breaks)],
    xmax = breaks[-1],
    proportion = y2
  )
  
  # Apply gap
  bin_width <- diff(breaks)[1]
  gap_frac <- 0.10
  inset <- bin_width * gap_frac / 2
  
  df1$xmin_gap <- df1$xmin + inset
  df1$xmax_gap <- df1$xmax - inset
  
  df2$xmin_gap <- df2$xmin + inset
  df2$xmax_gap <- df2$xmax - inset
  
  # Plot 1
  p1 <- ggplot2::ggplot(df1) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin_gap,
        xmax = xmax_gap,
        ymin = 0,
        ymax = proportion
      ),
      fill = "black",
      color = NA
    ) +
    ggplot2::coord_cartesian(
      xlim = range(breaks),
      ylim = c(0, ymax)
    ) +
    ggplot2::labs(
      title = paste0(label1, ": ", d1$matched_name),
      x = "NICHES value",
      y = "Realized edges / total edges"
    ) +
    ggplot2::theme_classic()
  
  # Plot 2
  p2 <- ggplot2::ggplot(df2) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin_gap,
        xmax = xmax_gap,
        ymin = 0,
        ymax = proportion
      ),
      fill = "red",
      color = NA
    ) +
    ggplot2::coord_cartesian(
      xlim = range(breaks),
      ylim = c(0, ymax)
    ) +
    ggplot2::labs(
      title = paste0(label2, ": ", d2$matched_name),
      x = "NICHES value",
      y = "Realized edges / total edges"
    ) +
    ggplot2::theme_classic()
  
  # Console summaries
  cat("Requested mechanism:", mechanism, "\n\n")
  
  cat("Plot 1:", label1, "\n")
  cat("Matched mechanism:", d1$matched_name, "\n")
  cat("Total edges:", total_edges1, "\n")
  cat("Zeros:", sum(d1$x == 0, na.rm = TRUE), "\n")
  cat("Nonzero:", sum(d1$x != 0, na.rm = TRUE), "\n")
  cat("NA:", sum(is.na(d1$x)), "\n\n")
  print(summary(if (include_zeros) d1$x else x1))
  cat("\n")
  
  cat("Plot 2:", label2, "\n")
  cat("Matched mechanism:", d2$matched_name, "\n")
  cat("Total edges:", total_edges2, "\n")
  cat("Zeros:", sum(d2$x == 0, na.rm = TRUE), "\n")
  cat("Nonzero:", sum(d2$x != 0, na.rm = TRUE), "\n")
  cat("NA:", sum(is.na(d2$x)), "\n\n")
  print(summary(if (include_zeros) d2$x else x2))
  
  p <- p1 + p2 + patchwork::plot_layout(ncol = 2)
  
  attr(p, "plot1_data") <- list(
    values = x1,
    matched_mechanism = d1$matched_name,
    total_edges = total_edges1,
    hist = h1,
    proportions = y1
  )
  
  attr(p, "plot2_data") <- list(
    values = x2,
    matched_mechanism = d2$matched_name,
    total_edges = total_edges2,
    hist = h2,
    proportions = y2
  )
  
  attr(p, "breaks") <- breaks
  attr(p, "xlim") <- range(breaks)
  attr(p, "ylim") <- c(0, ymax)
  
  return(p)
}

#### Plot quartile pies of one mechanism from an EffNICHES object ####

plot_mechanism_quartile_pies <- function(
    eff,
    mechanism,
    group_col,
    include_zeros = FALSE,
    palette = NULL,
    inside_label_min_prop = 0.05,
    other_pool_thresh = 0.02,
    drop_na_groups = FALSE,
    other_label = "Other",
    pie_radius = 1,
    donut_hole = 0.08,
    repel_force = 1.0,
    repel_box_padding = 0.25,
    repel_point_padding = 0.05,
    repel_segment_size = 0.3,
    inside_label_size = 2.8,
    outside_label_size = 2.6,
    slice_convention = c("S3", "S4")
) {
  slice_convention <- match.arg(slice_convention)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Package 'ggforce' is required.")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.")
  }
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggforce)
    library(ggrepel)
  })
  
  w    <- eff$niches$CellToCellSpatial$w
  ij   <- eff$niches$CellToCellSpatial$ij
  ids  <- eff$nodes$ids
  meta <- eff$nodes$meta
  
  if (is.null(colnames(w))) {
    stop("eff$niches$CellToCellSpatial$w has no column names.")
  }
  if (!group_col %in% colnames(meta)) {
    stop("group_col not found in eff$nodes$meta")
  }
  if (nrow(ij) != nrow(w)) {
    stop("nrow(ij) must equal nrow(w)")
  }
  if (length(ids) < max(ij)) {
    stop("eff$nodes$ids is shorter than the indices in ij")
  }
  if (other_pool_thresh <= 0 || other_pool_thresh >= 1) {
    stop("other_pool_thresh must be between 0 and 1")
  }
  if (inside_label_min_prop <= 0 || inside_label_min_prop >= 1) {
    stop("inside_label_min_prop must be between 0 and 1")
  }
  if (inside_label_min_prop < other_pool_thresh) {
    stop("inside_label_min_prop should be >= other_pool_thresh")
  }
  
  mech_names <- colnames(w)
  mech_names_norm <- trimws(tolower(norm_lr_dash(mech_names)))
  mechanism_norm  <- trimws(tolower(norm_lr_dash(mechanism)))
  hit_idx <- which(mech_names_norm == mechanism_norm)
  
  if (length(hit_idx) == 0) {
    stop("Mechanism not found in eff$niches$CellToCellSpatial$w")
  }
  if (length(hit_idx) > 1) {
    stop("Multiple mechanism matches found after dash normalization")
  }
  
  mechanism_matched <- mech_names[hit_idx]
  x <- as.numeric(w[, hit_idx])
  
  edge_df <- data.frame(
    edge_index   = seq_along(x),
    sender_idx   = ij[, "sender"],
    receiver_idx = ij[, "receiver"],
    weight       = x,
    stringsAsFactors = FALSE
  )
  
  if (!include_zeros) {
    edge_df <- edge_df[!is.na(edge_df$weight) & edge_df$weight > 0, , drop = FALSE]
  } else {
    edge_df <- edge_df[!is.na(edge_df$weight), , drop = FALSE]
  }
  
  if (nrow(edge_df) == 0) {
    stop("No edges remain after filtering.")
  }
  
  edge_df$sender_id   <- ids[edge_df$sender_idx]
  edge_df$receiver_id <- ids[edge_df$receiver_idx]
  
  if (!is.null(rownames(meta)) &&
      all(edge_df$sender_id %in% rownames(meta)) &&
      all(edge_df$receiver_id %in% rownames(meta))) {
    edge_df$sender_group   <- as.character(meta[edge_df$sender_id, group_col])
    edge_df$receiver_group <- as.character(meta[edge_df$receiver_id, group_col])
  } else {
    edge_df$sender_group   <- as.character(meta[[group_col]][match(edge_df$sender_id, ids)])
    edge_df$receiver_group <- as.character(meta[[group_col]][match(edge_df$receiver_id, ids)])
  }
  
  if (drop_na_groups) {
    edge_df <- edge_df[
      !is.na(edge_df$sender_group) & !is.na(edge_df$receiver_group),
      ,
      drop = FALSE
    ]
  } else {
    edge_df$sender_group[is.na(edge_df$sender_group)]     <- "NA"
    edge_df$receiver_group[is.na(edge_df$receiver_group)] <- "NA"
  }
  
  if (nrow(edge_df) == 0) {
    stop("No edges remain after group assignment.")
  }
  
  qs <- as.numeric(stats::quantile(
    edge_df$weight,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
  
  edge_df$quartile <- NA_integer_
  for (i in 1:4) {
    lower <- qs[i]
    upper <- qs[i + 1]
    
    if (i < 4) {
      keep <- edge_df$weight >= lower & edge_df$weight < upper
    } else {
      keep <- edge_df$weight >= lower & edge_df$weight <= upper
    }
    
    edge_df$quartile[keep] <- i
  }
  
  if (any(is.na(edge_df$quartile))) {
    ord <- order(edge_df$weight)
    ranks <- integer(nrow(edge_df))
    ranks[ord] <- seq_len(nrow(edge_df))
    edge_df$quartile <- ceiling(4 * ranks / nrow(edge_df))
    edge_df$quartile[edge_df$quartile < 1] <- 1
    edge_df$quartile[edge_df$quartile > 4] <- 4
  }
  
  long_df <- rbind(
    data.frame(
      quartile = edge_df$quartile,
      role     = "Sender",
      group    = edge_df$sender_group,
      weight   = edge_df$weight,
      stringsAsFactors = FALSE
    ),
    data.frame(
      quartile = edge_df$quartile,
      role     = "Receiver",
      group    = edge_df$receiver_group,
      weight   = edge_df$weight,
      stringsAsFactors = FALSE
    )
  )
  
  raw_plot_df <- long_df %>%
    dplyr::count(quartile, role, group, name = "n") %>%
    dplyr::group_by(quartile, role) %>%
    dplyr::mutate(raw_prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  plot_df <- raw_plot_df %>%
    dplyr::mutate(group_plot = ifelse(raw_prop < other_pool_thresh, other_label, group)) %>%
    dplyr::group_by(quartile, role, group_plot) %>%
    dplyr::summarise(n = sum(n), .groups = "drop") %>%
    dplyr::group_by(quartile, role) %>%
    dplyr::mutate(
      prop = n / sum(n),
      pct  = 100 * prop
    ) %>%
    dplyr::ungroup()
  
  q_info <- edge_df %>%
    dplyr::group_by(quartile) %>%
    dplyr::summarise(
      q_min   = min(weight, na.rm = TRUE),
      q_max   = max(weight, na.rm = TRUE),
      n_edges = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      quartile_lab = paste0(
        "Q", quartile,
        "\n",
        "n=", n_edges,
        "\n[", signif(q_min, 4), ", ", signif(q_max, 4), "]"
      )
    )
  
  plot_df <- plot_df %>%
    dplyr::left_join(q_info[, c("quartile", "quartile_lab")], by = "quartile") %>%
    dplyr::mutate(
      role = factor(role, levels = c("Sender", "Receiver")),
      quartile_lab = factor(
        quartile_lab,
        levels = q_info$quartile_lab[order(q_info$quartile)]
      )
    ) %>%
    dplyr::group_by(quartile_lab, role) %>%
    dplyr::arrange(
      dplyr::if_else(group_plot == other_label, 1L, 0L),
      dplyr::desc(prop),
      group_plot,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      frac       = prop / sum(prop),
      frac_end   = cumsum(frac),
      frac_start = dplyr::lag(frac_end, default = 0),
      frac_mid   = (frac_start + frac_end) / 2,
      
      ## Label geometry: L1
      label_angle = pi/2 - 2*pi*frac_mid,
      
      label_text = paste0(group_plot, "\n", sprintf("%.1f%%", pct))
    ) %>%
    dplyr::ungroup()
  
  ## Slice geometry chosen independently from label geometry
  if (slice_convention == "S3") {
    plot_df$plot_start <- 0 + 2*pi*plot_df$frac_start
    plot_df$plot_end   <- 0 + 2*pi*plot_df$frac_end
  } else if (slice_convention == "S4") {
    plot_df$plot_start <- 0 + 2*pi*plot_df$frac_end
    plot_df$plot_end   <- 0 + 2*pi*plot_df$frac_start
  }
  
  all_groups <- unique(plot_df$group_plot)
  
  if (is.null(palette)) {
    base_groups <- setdiff(sort(all_groups), other_label)
    pal_vals <- grDevices::hcl.colors(length(base_groups), palette = "Set 3")
    group_cols <- stats::setNames(pal_vals, base_groups)
    
    if (other_label %in% all_groups) {
      group_cols <- c(group_cols, stats::setNames("grey70", other_label))
    }
    
    group_cols <- group_cols[unique(c(base_groups, intersect(other_label, all_groups)))]
  } else {
    if (length(palette) < length(all_groups)) {
      stop("palette does not contain enough colors for all groups.")
    }
    ordered_groups <- c(setdiff(sort(all_groups), other_label), intersect(other_label, all_groups))
    group_cols <- stats::setNames(palette[seq_along(ordered_groups)], ordered_groups)
  }
  
  inside_df <- plot_df %>%
    dplyr::filter(prop >= inside_label_min_prop) %>%
    dplyr::mutate(
      x = 0.58 * pie_radius * cos(label_angle),
      y = 0.58 * pie_radius * sin(label_angle)
    )
  
  outside_df <- plot_df %>%
    dplyr::filter(prop < inside_label_min_prop & prop >= other_pool_thresh) %>%
    dplyr::mutate(
      x = 1.03 * pie_radius * cos(label_angle),
      y = 1.03 * pie_radius * sin(label_angle),
      nudge_x = 0.18 * sign(cos(label_angle)),
      nudge_y = 0.05 * sin(label_angle)
    )
  
  x_lim <- c(-1.65 * pie_radius, 1.65 * pie_radius)
  y_lim <- c(-1.35 * pie_radius, 1.35 * pie_radius)
  
  p <- ggplot2::ggplot() +
    ggforce::geom_arc_bar(
      data = plot_df,
      ggplot2::aes(
        x0 = 0,
        y0 = 0,
        r0 = donut_hole * pie_radius,
        r  = pie_radius,
        start = plot_start,
        end   = plot_end,
        fill  = group_plot
      ),
      color = "white",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = inside_df,
      ggplot2::aes(
        x = x,
        y = y,
        label = label_text
      ),
      size = inside_label_size,
      lineheight = 0.9
    ) +
    ggrepel::geom_text_repel(
      data = outside_df,
      ggplot2::aes(
        x = x,
        y = y,
        label = label_text
      ),
      size = outside_label_size,
      lineheight = 0.9,
      force = repel_force,
      box.padding = repel_box_padding,
      point.padding = repel_point_padding,
      min.segment.length = 0,
      segment.size = repel_segment_size,
      max.overlaps = Inf,
      seed = 1,
      direction = "both",
      nudge_x = outside_df$nudge_x,
      nudge_y = outside_df$nudge_y
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(quartile_lab),
      cols = ggplot2::vars(role)
    ) +
    ggplot2::coord_fixed(
      xlim = x_lim,
      ylim = y_lim,
      clip = "off"
    ) +
    ggplot2::scale_fill_manual(
      values = group_cols,
      breaks = names(group_cols),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = mechanism_matched,
      subtitle = paste0(
        "Grouping: ", group_col,
        if (include_zeros) " | all edges" else " | nonzero edges only",
        " | <", 100 * other_pool_thresh, "% pooled as ", other_label,
        " | <", 100 * inside_label_min_prop, "% labeled outside",
        " | slice=", slice_convention, ", label=L1"
      ),
      fill = group_col
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(12, 50, 12, 50)
    )
  
  invisible(list(
    plot = p,
    mechanism = mechanism_matched,
    group_col = group_col,
    include_zeros = include_zeros,
    quartile_breaks = qs,
    edge_table = edge_df,
    plot_table = plot_df,
    colors = group_cols,
    slice_convention = slice_convention,
    label_convention = "L1"
  ))
}

#### Plot Mechanism Histogram Matrix ####

plotMechanismGroupPairHistMatrix <- function(
    eff1,
    eff2,
    group_col,
    mechanism,
    include_zeros = FALSE,
    bins = 50,
    label1 = "Object 1",
    label2 = "Object 2",
    fill1 = "black",
    fill2 = "red",
    alpha1 = 1,
    alpha2 = 1,
    gap_frac = 0.10,
    panel_border_color = "grey40",
    panel_border_linewidth = 0.45,
    axis_line_color = "grey55",
    axis_linewidth = 0.20,
    strip_text_size = 10,
    base_text_size = 10,
    verbose = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  
  # ---------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------
  .norm_long_dash <- function(x) {
    x <- gsub("\u2013", "\u2014", x, fixed = TRUE) # en dash -> em dash
    x <- gsub("\u2212", "\u2014", x, fixed = TRUE) # minus  -> em dash
    x
  }
  
  .resolve_mechanism_col <- function(w, mechanism) {
    mech_names <- colnames(w)
    mech_names_norm <- norm_lr_dash(mech_names)
    mechanism_norm  <- norm_lr_dash(mechanism)
    
    hit_idx <- which(mech_names_norm == mechanism_norm)
    if (length(hit_idx) == 0) {
      stop("Mechanism not found in eff$niches$CellToCellSpatial$w")
    }
    if (length(hit_idx) > 1) {
      stop("Multiple mechanism matches found after dash normalization")
    }
    
    list(
      idx = hit_idx,
      matched_name = mech_names[hit_idx]
    )
  }
  
  .resolve_node_ids <- function(eff, obj_label) {
    node_ids <- NULL
    if (!is.null(eff$nodes$ids)) {
      node_ids <- as.character(eff$nodes$ids)
    } else if (!is.null(eff$nodes$xy) && !is.null(rownames(eff$nodes$xy))) {
      node_ids <- rownames(eff$nodes$xy)
    } else if (!is.null(eff$nodes$meta) && !is.null(rownames(eff$nodes$meta))) {
      node_ids <- rownames(eff$nodes$meta)
    }
    
    if (is.null(node_ids)) {
      stop(sprintf("[%s] Could not resolve node IDs from eff.", obj_label))
    }
    node_ids
  }
  
  .resolve_meta <- function(eff, obj_label, node_ids, group_col) {
    md <- eff$nodes$meta
    if (is.null(md) || !is.data.frame(md)) {
      stop(sprintf("[%s] eff$nodes$meta must be a data.frame.", obj_label))
    }
    if (is.null(rownames(md))) {
      stop(sprintf("[%s] eff$nodes$meta must have rownames matching node IDs.", obj_label))
    }
    if (!(group_col %in% colnames(md))) {
      stop(sprintf("[%s] group_col '%s' not found in eff$nodes$meta.", obj_label, group_col))
    }
    
    missing_ids <- setdiff(node_ids, rownames(md))
    if (length(missing_ids) > 0) {
      stop(sprintf(
        "[%s] Metadata rownames do not cover all node IDs. Missing %d IDs (e.g. %s).",
        obj_label, length(missing_ids), paste(head(missing_ids, 5), collapse = ", ")
      ))
    }
    
    md[node_ids, , drop = FALSE]
  }
  
  .parse_edge_pairs <- function(edge_ids, obj_label) {
    if (is.null(edge_ids)) {
      stop(sprintf("[%s] Row names of w are required to recover sender/receiver IDs.", obj_label))
    }
    
    edge_ids_norm <- .norm_long_dash(edge_ids)
    pieces <- strsplit(edge_ids_norm, "\u2014", fixed = TRUE)
    lens <- lengths(pieces)
    
    bad <- which(lens != 2L)
    if (length(bad) > 0) {
      stop(sprintf(
        "[%s] Could not parse %d edge IDs into sender/receiver pairs using an em/en/minus dash separator. Example bad ID: %s",
        obj_label, length(bad), edge_ids[bad[1]]
      ))
    }
    
    sender <- vapply(pieces, `[`, character(1), 1)
    receiver <- vapply(pieces, `[`, character(1), 2)
    
    data.frame(
      sender_id = sender,
      receiver_id = receiver,
      stringsAsFactors = FALSE
    )
  }
  
  .extract_object_edge_table <- function(eff, mechanism, group_col, obj_label, plot_label) {
    w <- eff$niches$CellToCellSpatial$w
    if (is.null(w)) {
      stop(sprintf("[%s] eff$niches$CellToCellSpatial$w not found.", obj_label))
    }
    
    mech_hit <- .resolve_mechanism_col(w, mechanism)
    node_ids <- .resolve_node_ids(eff, obj_label)
    md <- .resolve_meta(eff, obj_label, node_ids, group_col)
    edge_pairs <- .parse_edge_pairs(rownames(w), obj_label)
    
    sender_group <- as.character(md[edge_pairs$sender_id, group_col])
    receiver_group <- as.character(md[edge_pairs$receiver_id, group_col])
    
    if (length(sender_group) != nrow(w) || length(receiver_group) != nrow(w)) {
      stop(sprintf("[%s] Internal alignment failure when mapping sender/receiver groups.", obj_label))
    }
    
    data.frame(
      object = plot_label,
      sender = sender_group,
      receiver = receiver_group,
      value = as.numeric(w[, mech_hit$idx]),
      stringsAsFactors = FALSE
    )
  }
  
  .build_hist_rects_rotated <- function(df_obj, breaks, total_edges, side, object_label) {
    # side = "right" for object 1, "left" for object 2
    if (nrow(df_obj) == 0) {
      return(data.frame(
        sender = character(0),
        receiver = character(0),
        xmin = numeric(0),
        xmax = numeric(0),
        ymin = numeric(0),
        ymax = numeric(0),
        object = character(0),
        stringsAsFactors = FALSE
      ))
    }
    
    split_key <- paste(df_obj$sender, df_obj$receiver, sep = "___")
    split_list <- split(df_obj$value, split_key, drop = TRUE)
    
    rect_list <- lapply(names(split_list), function(k) {
      vals <- split_list[[k]]
      h <- hist(vals, breaks = breaks, plot = FALSE)
      prop <- h$counts / total_edges
      
      parts <- strsplit(k, "___", fixed = TRUE)[[1]]
      sender_i <- parts[1]
      receiver_i <- parts[2]
      
      out <- data.frame(
        sender = sender_i,
        receiver = receiver_i,
        # rotate -90°: histogram value axis becomes vertical position
        ymin = breaks[-length(breaks)],
        ymax = breaks[-1],
        proportion = prop,
        object = object_label,
        stringsAsFactors = FALSE
      )
      
      if (identical(side, "right")) {
        out$xmin <- 0
        out$xmax <- out$proportion
      } else {
        out$xmin <- -out$proportion
        out$xmax <- 0
      }
      
      out
    })
    
    do.call(rbind, rect_list)
  }
  
  # ---------------------------------------------------------
  # Pull per-edge tables
  # ---------------------------------------------------------
  df1_all <- .extract_object_edge_table(
    eff = eff1,
    mechanism = mechanism,
    group_col = group_col,
    obj_label = "eff1",
    plot_label = label1
  )
  
  df2_all <- .extract_object_edge_table(
    eff = eff2,
    mechanism = mechanism,
    group_col = group_col,
    obj_label = "eff2",
    plot_label = label2
  )
  
  matched_name1 <- .resolve_mechanism_col(eff1$niches$CellToCellSpatial$w, mechanism)$matched_name
  matched_name2 <- .resolve_mechanism_col(eff2$niches$CellToCellSpatial$w, mechanism)$matched_name
  
  total_edges1 <- nrow(eff1$niches$CellToCellSpatial$w)
  total_edges2 <- nrow(eff2$niches$CellToCellSpatial$w)
  
  # ---------------------------------------------------------
  # Group level handling
  # ---------------------------------------------------------
  groups1 <- sort(unique(c(df1_all$sender, df1_all$receiver)))
  groups2 <- sort(unique(c(df2_all$sender, df2_all$receiver)))
  
  groups1 <- groups1[!is.na(groups1)]
  groups2 <- groups2[!is.na(groups2)]
  
  all_groups <- sort(union(groups1, groups2))
  if (length(all_groups) == 0) {
    stop("No non-NA group labels were found.")
  }
  
  if (!setequal(groups1, groups2)) {
    warning("Group labels are not identical across objects. Using the union of labels for rows and columns.")
  }
  
  # ---------------------------------------------------------
  # NA handling / value filtering
  # ---------------------------------------------------------
  df1_plot <- df1_all[!is.na(df1_all$sender) & !is.na(df1_all$receiver) & !is.na(df1_all$value), , drop = FALSE]
  df2_plot <- df2_all[!is.na(df2_all$sender) & !is.na(df2_all$receiver) & !is.na(df2_all$value), , drop = FALSE]
  
  if (!include_zeros) {
    df1_plot <- df1_plot[df1_plot$value != 0, , drop = FALSE]
    df2_plot <- df2_plot[df2_plot$value != 0, , drop = FALSE]
  }
  
  combined_x <- c(df1_plot$value, df2_plot$value)
  combined_x <- combined_x[!is.na(combined_x)]
  
  if (length(combined_x) == 0) {
    stop("No non-NA values available for plotting after filtering.")
  }
  
  value_range <- range(combined_x)
  if (diff(value_range) == 0) {
    value_range <- value_range + c(-0.5, 0.5)
  }
  
  breaks <- seq(value_range[1], value_range[2], length.out = bins + 1)
  
  # ---------------------------------------------------------
  # Histogram rectangles per sender->receiver panel
  # ---------------------------------------------------------
  rect1 <- .build_hist_rects_rotated(
    df_obj = df1_plot,
    breaks = breaks,
    total_edges = total_edges1,
    side = "left",
    object_label = label1
  )
  
  rect2 <- .build_hist_rects_rotated(
    df_obj = df2_plot,
    breaks = breaks,
    total_edges = total_edges2,
    side = "right",
    object_label = label2
  )
  
  rect_df <- rbind(rect1, rect2)
  
  # Preserve gaps between bars, now along the vertical axis
  bin_width <- diff(breaks)[1]
  inset <- bin_width * gap_frac / 2
  
  if (nrow(rect_df) > 0) {
    rect_df$ymin_gap <- rect_df$ymin + inset
    rect_df$ymax_gap <- rect_df$ymax - inset
    
    bad_gap <- rect_df$ymin_gap > rect_df$ymax_gap
    rect_df$ymin_gap[bad_gap] <- rect_df$ymin[bad_gap]
    rect_df$ymax_gap[bad_gap] <- rect_df$ymax[bad_gap]
  }
  
  # ---------------------------------------------------------
  # Shared x-scale across all panels (mirrored proportions)
  # ---------------------------------------------------------
  xmax <- 0
  if (nrow(rect_df) > 0) {
    xmax <- max(abs(c(rect_df$xmin, rect_df$xmax)), na.rm = TRUE)
  }
  if (!is.finite(xmax) || xmax <= 0) xmax <- 1
  
  # ---------------------------------------------------------
  # Force full sender x receiver grid
  # ---------------------------------------------------------
  row_levels <- all_groups
  col_levels <- all_groups
  
  panel_df <- expand.grid(
    sender = row_levels,
    receiver = col_levels,
    stringsAsFactors = FALSE
  )
  
  panel_df$sender <- factor(panel_df$sender, levels = row_levels)
  panel_df$receiver <- factor(panel_df$receiver, levels = col_levels)
  panel_df$x0 <- 0
  panel_df$y0 <- mean(range(breaks))
  
  if (nrow(rect_df) > 0) {
    rect_df$sender <- factor(rect_df$sender, levels = row_levels)
    rect_df$receiver <- factor(rect_df$receiver, levels = col_levels)
  }
  
  # ---------------------------------------------------------
  # Diagnostics
  # ---------------------------------------------------------
  count1_full <- as.data.frame(table(df1_all$sender, df1_all$receiver), stringsAsFactors = FALSE)
  count2_full <- as.data.frame(table(df2_all$sender, df2_all$receiver), stringsAsFactors = FALSE)
  colnames(count1_full) <- c("sender", "receiver", "n_edges")
  colnames(count2_full) <- c("sender", "receiver", "n_edges")
  
  if (verbose) {
    cat("Requested mechanism:", mechanism, "\n\n")
    
    cat("Object 1:", label1, "\n")
    cat("Matched mechanism:", matched_name1, "\n")
    cat("Total edges:", total_edges1, "\n")
    cat("Zeros:", sum(df1_all$value == 0, na.rm = TRUE), "\n")
    cat("Nonzero:", sum(df1_all$value != 0, na.rm = TRUE), "\n")
    cat("NA values:", sum(is.na(df1_all$value)), "\n")
    cat("Edges with NA sender/receiver group:", sum(is.na(df1_all$sender) | is.na(df1_all$receiver)), "\n\n")
    print(summary(if (include_zeros) df1_all$value else df1_plot$value))
    cat("\n")
    
    cat("Object 2:", label2, "\n")
    cat("Matched mechanism:", matched_name2, "\n")
    cat("Total edges:", total_edges2, "\n")
    cat("Zeros:", sum(df2_all$value == 0, na.rm = TRUE), "\n")
    cat("Nonzero:", sum(df2_all$value != 0, na.rm = TRUE), "\n")
    cat("NA values:", sum(is.na(df2_all$value)), "\n")
    cat("Edges with NA sender/receiver group:", sum(is.na(df2_all$sender) | is.na(df2_all$receiver)), "\n\n")
    print(summary(if (include_zeros) df2_all$value else df2_plot$value))
    cat("\n")
  }
  
  # ---------------------------------------------------------
  # Plot
  # ---------------------------------------------------------
  p <- ggplot2::ggplot() +
    ggplot2::geom_blank(
      data = panel_df,
      ggplot2::aes(x = x0, y = y0)
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = axis_line_color,
      linewidth = axis_linewidth
    )
  
  if (nrow(rect_df) > 0) {
    p <- p +
      ggplot2::geom_rect(
        data = rect_df[rect_df$object == label1, , drop = FALSE],
        ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin_gap,
          ymax = ymax_gap
        ),
        fill = fill1,
        color = NA,
        alpha = alpha1
      ) +
      ggplot2::geom_rect(
        data = rect_df[rect_df$object == label2, , drop = FALSE],
        ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin_gap,
          ymax = ymax_gap
        ),
        fill = fill2,
        color = NA,
        alpha = alpha2
      )
  }
  
  title_txt <- paste0(
    matched_name1,
    " | mirrored distributions by ",
    group_col,
    " sender \u2192 ",
    group_col,
    " receiver"
  )
  
  subtitle_txt <- paste0(
    label1, " (left, ", fill1, ") vs ",
    label2, " (right, ", fill2, "); ",
    if (include_zeros) "zeros included" else "zeros excluded",
    "; x = bin count / total edges in object"
  )
  
  p <- p +
    ggplot2::facet_grid(
      rows = ggplot2::vars(sender),
      cols = ggplot2::vars(receiver),
      drop = FALSE,
      switch = "y"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-xmax, xmax),
      ylim = range(breaks),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Mirrored proportion of total object edges",
      y = "NICHES value"
    ) +
    ggplot2::theme_classic(base_size = base_text_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "white", color = panel_border_color, linewidth = panel_border_linewidth),
      strip.text = ggplot2::element_text(size = strip_text_size, face = "bold"),
      strip.placement = "outside",
      panel.border = ggplot2::element_rect(
        color = panel_border_color,
        fill = NA,
        linewidth = panel_border_linewidth
      ),
      panel.spacing = grid::unit(0, "lines"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  
  # ---------------------------------------------------------
  # Attach useful internals
  # ---------------------------------------------------------
  attr(p, "plot1_data") <- list(
    edge_table = df1_all,
    plotted_edge_table = df1_plot,
    matched_mechanism = matched_name1,
    total_edges = total_edges1,
    counts_by_panel = count1_full
  )
  
  attr(p, "plot2_data") <- list(
    edge_table = df2_all,
    plotted_edge_table = df2_plot,
    matched_mechanism = matched_name2,
    total_edges = total_edges2,
    counts_by_panel = count2_full
  )
  
  attr(p, "rect_data") <- rect_df
  attr(p, "panel_data") <- panel_df
  attr(p, "breaks") <- breaks
  attr(p, "xlim") <- c(-xmax, xmax)
  attr(p, "ylim") <- range(breaks)
  attr(p, "group_levels") <- all_groups
  
  return(p)
}
