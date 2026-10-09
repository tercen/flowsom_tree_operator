# The FlowSOM tree from a passed model: rebuild an R FlowSOM object from the JSON model that
# flowsom_operator (>= 2.2.0) writes, and draw FlowSOM's own plots from it.
#
# Why the model and not the labels: FlowSOM's tree is the minimum spanning tree of the map's
# codes, and the stars are the nodes' medians. The per-cell labels give neither. The model
# carries both, so the tree drawn here is the tree FlowSOM itself would draw for that map.
suppressMessages({
  library(FlowSOM)
  library(ggplot2)
  library(jsonlite)
})

# Palettes by name, the same lists Tercen's Plot operator uses, so a population has the same
# colour here as in the workflow's other plots (levels sorted, colours assigned in order).
PALETTES <- list(
  "tab20" = c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
              "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
              "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"),
  "tab10" = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
              "#7f7f7f", "#bcbd22", "#17becf"),
  "Palette-1" = c("#1F78B4", "#E31A1C", "#33A02C", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                  "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
  "Set1" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
             "#F781BF", "#999999"),
  "Dark2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
              "#666666")
)

palette_for <- function(levels, name) {
  pal <- PALETTES[[name]]
  if (is.null(pal)) stop("unknown palette '", name, "'; known: ", paste(names(PALETTES), collapse = ", "))
  setNames(rep_len(pal, length(levels)), levels)
}

# Parse and check the model JSON.
read_model <- function(json) {
  m <- fromJSON(json, simplifyVector = TRUE)
  if (!identical(m$format, "flowsom_model")) stop("not a flowsom_model (format = ", m$format, ")")
  if (m$version != 1) stop("flowsom_model version ", m$version, " is not supported (expected 1)")
  ncodes <- m$xdim * m$ydim
  p <- length(m$markers)
  codes <- matrix(unlist(m$codes), nrow = ncodes, ncol = p, byrow = TRUE)
  # medians: one vector per node, NULL for an empty node.
  med <- matrix(NA_real_, nrow = ncodes, ncol = p)
  meds <- if (is.list(m$medians)) m$medians else lapply(seq_len(nrow(m$medians)), function(k) m$medians[k, ])
  for (k in seq_len(ncodes)) if (length(meds[[k]]) == p) med[k, ] <- unlist(meds[[k]])
  colnames(codes) <- colnames(med) <- m$markers
  m$codes <- codes
  m$medians <- med
  m
}

# Node number from a label like "c007" (or 7).
node_number <- function(x) as.integer(gsub("[^0-9]", "", as.character(x)))

# A FlowSOM object for the model, with `node` the node of each projected cell (1-based). The
# node sizes are the model's own counts (every cell the map assigned); the cells here are used
# only for the pies and the node labels.
rebuild_fsom <- function(m, node) {
  ncodes <- nrow(m$codes)
  # FlowSOM computes metacluster MFIs from `data`, one row per cell. Nothing drawn here uses
  # them (stars and markers read the node medians), so each cell stands in as its node's median.
  stand_in <- m$medians[node, , drop = FALSE]
  stand_in[is.na(stand_in)] <- 0
  fs <- list(data = stand_in, compensate = FALSE, spillover = NULL, transform = FALSE,
             scale = FALSE, prettyColnames = setNames(m$markers, m$markers))
  fs$map <- list(codes = m$codes, grid = expand.grid(seq_len(m$xdim), seq_len(m$ydim)),
                 xdim = m$xdim, ydim = m$ydim, colsUsed = seq_along(m$markers),
                 mapping = cbind(node, 0), distf = 2, nNodes = ncodes)
  counts <- unlist(m$counts)
  fs$map$pctgs <- setNames(counts / sum(counts), seq_len(ncodes))
  fs$map$medianValues <- m$medians
  fs$map$cvValues <- fs$map$sdValues <- fs$map$madValues <- m$medians
  fs$metaclustering <- factor(unlist(m$metaclustering))
  class(fs) <- "FlowSOM"
  BuildMST(fs)
}

# The label of every node: the most common population among its cells; an empty node takes the
# most common population of its metacluster. Without populations, the metacluster.
node_labels <- function(m, node, population = NULL) {
  ncodes <- nrow(m$codes)
  meta <- unlist(m$metaclustering)
  width <- nchar(as.character(max(meta)))
  if (is.null(population)) {
    return(factor(sprintf(paste0("metacluster %0", width, "d"), meta)))
  }
  population <- as.character(population)
  majority <- function(x) if (length(x)) names(which.max(table(x))) else NA_character_
  lab <- vapply(seq_len(ncodes), function(k) majority(population[node == k]), "")
  by_meta <- vapply(seq_len(max(meta)), function(j) majority(population[meta[node] == j]), "")
  empty <- is.na(lab)
  lab[empty] <- by_meta[meta[empty]]
  lab[is.na(lab)] <- "(no cells)"
  factor(lab, levels = sort(unique(lab)))
}

# Each plot as (filename, ggplot-like object). Filenames sort into reading order, because the
# report export orders a step's images by filename and titles each page with it.
tree_plots <- function(fs, bg, bg_colors, group = NULL, group_colors = NULL,
                       markers = character(0), max_node_size = 1,
                       which = c("stars", "numbers", "grid", "pies", "markers")) {
  out <- list()
  # Stars on the tree, background = population (FlowSOM vignette: PlotStars).
  if ("stars" %in% which) out[["01 stars on the tree.png"]] <- PlotStars(fs, backgroundValues = bg,
                                                 backgroundColors = bg_colors,
                                                 maxNodeSize = max_node_size)
  # Nodes numbered and coloured by population (PlotNumbers).
  if ("numbers" %in% which) out[["02 node numbers on the tree.png"]] <- PlotNumbers(fs, level = "clusters",
                                                          maxNodeSize = max_node_size,
                                                          backgroundValues = bg,
                                                          backgroundColors = bg_colors)
  # The same stars on the SOM grid, equal sizes (PlotStars, view = "grid").
  if ("grid" %in% which) out[["03 stars on the SOM grid.png"]] <- PlotStars(fs, backgroundValues = bg,
                                                     backgroundColors = bg_colors,
                                                     view = "grid", equalNodeSize = TRUE)
  # Each node's cells by group (PlotPies).
  if (!is.null(group) && "pies" %in% which) {
    # FlowSOM adds an "Empty" slice for a node none of these cells fall on.
    out[["04 pies by group on the tree.png"]] <- PlotPies(fs, cellTypes = group,
                                                          colorPalette = c(group_colors, Empty = "#FFFFFF"),
                                                          backgroundValues = bg,
                                                          backgroundColors = bg_colors,
                                                          maxNodeSize = max_node_size)
  }
  # One tree per marker, nodes coloured by its median (PlotMarker).
  if (!"markers" %in% which) markers <- character(0)
  for (i in seq_along(markers)) {
    f <- sprintf("05 marker %02d %s.png", i, gsub("[/\\\\:*?\"<>|]", "_", markers[i]))
    out[[f]] <- PlotMarker(fs, markers[i], maxNodeSize = max_node_size)
  }
  out
}

save_png <- function(p, file, width, height, dpi) {
  grDevices::png(file, width = width, height = height, units = "in", res = dpi, type = "cairo",
                 bg = "white")
  on.exit(grDevices::dev.off())
  print(p)
  invisible(file)
}
