# flowsom_tree_operator: FlowSOM's tree, star, pie and marker plots, drawn from the model that
# flowsom_operator (>= 2.2.0) writes.
#
# Projection:
#   rows    = the channel and `<ns>.flowsom_model` (the model is the same on every row)
#   columns = the cells, with the node factor (`<ns>.cluster_id`), and optionally a population
#             factor (colours the nodes) and a group factor (the pies)
#   y       = anything (not read; Tercen needs a y)
# Output: one PNG per plot in an image table, which Tercen shows on the step and the report
# export picks up.
suppressMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(base64enc)
})
source("R/tree.R")

ctx <- tercenCtx()

prop <- function(name, default, type = as.character) ctx$op.value(name, type = type, default = default)
s <- list(
  node_factor = trimws(prop("node_factor", "")),
  population_factor = trimws(prop("population_factor", "")),
  group_factor = trimws(prop("group_factor", "")),
  markers = trimws(prop("markers", "")),
  plots = prop("plots", "stars, numbers, grid, pies, markers"),
  population_palette = prop("population_palette", "tab20"),
  group_palette = prop("group_palette", "Palette-1"),
  max_node_size = prop("max_node_size", 1, as.double),
  width = prop("plot_width", 9, as.double),
  height = prop("plot_height", 7, as.double),
  dpi = prop("dpi", 150, as.double)
)

# The model: the row factor whose name ends in flowsom_model.
rnames <- unlist(ctx$rnames)
model_col <- rnames[grepl("flowsom_model$", rnames)]
if (length(model_col) != 1) {
  stop("Put the FlowSOM model on rows: project the `flowsom_model` column of a FlowSOM (>= 2.2.0) ",
       "step next to the channel. Row factors found: ", paste(rnames, collapse = ", "))
}
models <- unique(ctx$rselect(model_col)[[1]])
models <- models[!is.na(models) & nchar(models) > 0]
if (length(models) != 1) stop("expected one FlowSOM model on rows, found ", length(models))
m <- read_model(models)

# The cells: node, and the optional population and group.
cnames <- unlist(ctx$cnames)
pick <- function(name, what, pattern = NULL) {
  if (nchar(name) > 0) {
    if (!name %in% cnames) stop(what, " '", name, "' is not a column factor; found: ", paste(cnames, collapse = ", "))
    return(name)
  }
  if (is.null(pattern)) return(NULL)
  hit <- cnames[grepl(pattern, cnames)]
  if (length(hit) != 1) stop("set ", what, ": ", length(hit), " column factors match '", pattern, "'")
  hit
}
node_col <- pick(s$node_factor, "node_factor", "cluster_id$")
pop_col <- pick(s$population_factor, "population_factor")
group_col <- pick(s$group_factor, "group_factor")
cells <- ctx$cselect(as.list(c(node_col, pop_col, group_col)))
node <- node_number(cells[[node_col]])
if (anyNA(node) || any(node < 1) || any(node > m$xdim * m$ydim)) {
  stop("'", node_col, "' does not hold node numbers of this ", m$xdim, "x", m$ydim, " map")
}

fs <- rebuild_fsom(m, node)
bg <- node_labels(m, node, if (is.null(pop_col)) NULL else cells[[pop_col]])
group <- if (is.null(group_col)) NULL else factor(cells[[group_col]], levels = sort(unique(cells[[group_col]])))
which <- trimws(strsplit(s$plots, ",")[[1]])
unknown <- setdiff(which, c("stars", "numbers", "grid", "pies", "markers"))
if (length(unknown)) stop("unknown plots: ", paste(unknown, collapse = ", "),
                          "; choose from stars, numbers, grid, pies, markers")
markers <- if (nchar(s$markers) > 0) {
  want <- trimws(strsplit(s$markers, ",")[[1]])
  missing <- setdiff(want, m$markers)
  if (length(missing)) stop("markers not in the model: ", paste(missing, collapse = ", "))
  want
} else m$markers

plots <- tree_plots(fs, bg, palette_for(levels(bg), s$population_palette),
                    group, if (is.null(group)) NULL else palette_for(levels(group), s$group_palette),
                    markers = markers, max_node_size = s$max_node_size, which = which)

dir <- tempfile("flowsom_tree"); dir.create(dir)
out <- do.call(rbind, lapply(names(plots), function(f) {
  path <- save_png(plots[[f]], file.path(dir, f), s$width, s$height, s$dpi)
  data.frame(filename = f, mimetype = "image/png", checksum = as.vector(tools::md5sum(path)),
             .content = base64encode(readBin(path, "raw", file.info(path)$size)),
             stringsAsFactors = FALSE, check.names = FALSE)
}))
unlink(dir, recursive = TRUE)

# An image table joins on nothing: it is a standalone relation beside the crosstab, the shape
# Tercen's Plot operator writes and the report export looks for (a `.content` column).
out %>%
  ctx$addNamespace() %>%
  as_relation() %>%
  as_join_operator(list(), list()) %>%
  save_relation(ctx)
