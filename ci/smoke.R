# Every plot renders from a real model: the one flowsom_operator 2.2.0 writes on its own
# platform test (synthetic data, 3000 cells x 4 markers, 10x10 map, 14 empty nodes). Run in the
# image: docker run --rm -v "$PWD/ci:/c" --entrypoint Rscript <image> /c/smoke.R
source("/operator/R/tree.R")
m <- read_model(paste(readLines("/c/model_synthetic.txt", warn = FALSE), collapse = ""))
# The codes are the JSON's per-node vectors, node by node (1.0.x scrambled them).
raw <- jsonlite::fromJSON(paste(readLines("/c/model_synthetic.txt", warn = FALSE), collapse = ""), simplifyVector = FALSE)
for (k in c(1, 7, 100)) stopifnot(isTRUE(all.equal(unname(m$codes[k, ]), unlist(raw$codes[[k]]), tolerance = 0)))
# And the tree is the one FlowSOM builds from those codes.
ref_codes <- do.call(rbind, lapply(raw$codes, unlist))
ref_mst <- igraph::mst(igraph::graph_from_adjacency_matrix(as.matrix(dist(ref_codes)), mode = "undirected", weighted = TRUE))
node <- node_number(read.csv("/c/cells_synthetic.csv", stringsAsFactors = FALSE)$node)
meta <- unlist(m$metaclustering)
pop <- c("pop A", "pop B", "pop C", "pop D", "pop E")[meta[node]]
group <- factor(rep(c("group 1", "group 2"), length.out = length(node)))
m <- relabel_markers(m, c("CD3", "CD19", "CD14", "CD56"))
fs <- rebuild_fsom(m, node)
stopifnot(identical(colnames(fs$map$codes), m$markers))
stopifnot(nrow(fs$map$codes) == 100, sum(is.na(fs$map$medianValues[, 1])) == 14)
# The tree is built from the codes: one edge fewer than nodes.
stopifnot(igraph::ecount(fs$MST$graph) == 99)
edges <- function(g) sort(apply(igraph::as_edgelist(g), 1, function(e) paste(sort(e), collapse = "-")))
stopifnot(identical(edges(fs$MST$graph), edges(ref_mst)))
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
