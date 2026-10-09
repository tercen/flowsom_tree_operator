# Every plot renders from a real model: the one flowsom_operator 2.2.0 writes on its own
# platform test (synthetic data, 3000 cells x 4 markers, 10x10 map, 14 empty nodes). Run in the
# image: docker run --rm -v "$PWD/tests:/t" --entrypoint Rscript <image> /t/smoke.R
source("/operator/R/tree.R")
m <- read_model(paste(readLines("/t/model_synthetic.json", warn = FALSE), collapse = ""))
node <- node_number(read.csv("/t/cells_synthetic.csv", stringsAsFactors = FALSE)$node)
meta <- unlist(m$metaclustering)
pop <- c("pop A", "pop B", "pop C", "pop D", "pop E")[meta[node]]
group <- factor(rep(c("group 1", "group 2"), length.out = length(node)))
fs <- rebuild_fsom(m, node)
stopifnot(nrow(fs$map$codes) == 100, sum(is.na(fs$map$medianValues[, 1])) == 14)
# The tree is built from the codes: one edge fewer than nodes.
stopifnot(igraph::ecount(fs$MST$graph) == 99)
bg <- node_labels(m, node, pop)
stopifnot(length(bg) == 100, !anyNA(bg))
plots <- tree_plots(fs, bg, palette_for(levels(bg), "tab20"), group,
                    palette_for(levels(group), "Palette-1"), markers = m$markers)
stopifnot(length(plots) == 4 + length(m$markers))
d <- tempfile(); dir.create(d)
for (f in names(plots)) {
  save_png(plots[[f]], file.path(d, f), 9, 7, 72)
  stopifnot(file.info(file.path(d, f))$size > 10000)
  cat("ok", f, "\n")
}
