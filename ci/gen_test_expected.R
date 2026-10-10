# Expected output of tests/test.json, made in the operator image with main.R's own calls on the
# test input (the grid view only: it does not depend on igraph's float layout, so the PNG is
# byte-stable). docker run --rm -v "$PWD:/w" -w /w --entrypoint Rscript <image> ci/gen_test_expected.R
source("/operator/R/tree.R")
x <- read.csv("tests/flowsom_tree_input.csv", stringsAsFactors = FALSE)
m <- read_model(x$flowsom_model[1])
cells <- x[!duplicated(x$cell_id), ]
cells <- cells[order(cells$cell_id), ]
node <- node_number(cells$cluster_id)
fs <- rebuild_fsom(m, node)
bg <- node_labels(m, node, cells$population)
p <- tree_plots(fs, bg, palette_for(levels(bg), "tab20"), which = "grid")
f <- file.path(tempdir(), names(p)[1]); save_png(p[[1]], f, 6, 5, 50)
out <- data.frame(ds0.filename = names(p)[1], ds0.mimetype = "image/png",
                  ds0.checksum = as.vector(tools::md5sum(f)),
                  .content = base64enc::base64encode(readBin(f, "raw", file.info(f)$size)),
                  check.names = FALSE)
write.csv(out, "tests/table1.csv", row.names = FALSE)
cat("bytes", file.info(f)$size, "\n")
