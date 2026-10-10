# FlowSOM tree operator

Draws FlowSOM's own plots in a Tercen workflow — the tree (minimum spanning tree) with stars,
node numbers, pies and one tree per marker, and the stars on the SOM grid — from the **model**
that a FlowSOM step writes (`flowsom_operator` ≥ 2.2.0, column `flowsom_model`).

The tree is FlowSOM's tree for that map: the operator rebuilds the R `FlowSOM` object from the
model's codes, node medians, node counts and metaclustering, and calls FlowSOM 2.4's
`BuildMST`, `PlotStars`, `PlotNumbers`, `PlotPies` and `PlotMarker`. It is not re-estimated from
the cell labels: the MST is built from the codes, and node medians give a different tree (22 of
35 edges shared on a 6×6 test map), and an empty node has no median at all.

It replaces [`flowsom_mst_operator`](https://github.com/tercen/flowsom_mst_operator), which read
the serialised R object that FlowSOM 1.x operators wrote and the Rust FlowSOM (2.x) does not.

## Projection

| | |
|---|---|
| rows | the channel (and a factor to name markers by, e.g. `c.antigen`). The model is **not** projected: Tercen's engine cannot use a string over 4096 bytes as a factor, and a model is tens of KB. The operator reads it from the nearest FlowSOM (≥ 2.2.0) step upstream |
| columns | the cells, with `<ns>.cluster_id` (the SOM node), and optionally a population factor and a group factor |
| y | any value; not read |

Filter the cells to the ones you want counted in the pies and the node colours (for example the
study samples). Node sizes are the model's own counts: every cell the FlowSOM step assigned.

## Output

One PNG per plot, in an image table (`filename`, `mimetype`, `checksum`, `.content`) joined on
nothing, the shape Tercen's Plot operator writes; the step shows the images and the report
export picks them up. Files are named so they sort into reading order:

| file | FlowSOM call |
|---|---|
| `01 stars on the tree.png` | `PlotStars`, background = population |
| `02 node numbers on the tree.png` | `PlotNumbers`, background = population |
| `03 stars on the SOM grid.png` | `PlotStars(view = "grid", equalNodeSize = TRUE)` |
| `04 pies by group on the tree.png` | `PlotPies(cellTypes = group)` (needs `group_factor`) |
| `05 marker NN <marker>.png` | `PlotMarker`, one per marker |

## Properties

| name | default | |
|---|---|---|
| `node_factor` | (auto) | column factor with each cell's node; empty = the one ending in `cluster_id` |
| `population_factor` | | colours each node by its cells' most common population (an empty node: its metacluster's); empty = by metacluster |
| `group_factor` | | the pies; empty = no pie plot |
| `plots` | `stars, numbers, grid, pies, markers` | which plots |
| `markers` | | per-marker trees (by their label); empty = every marker in the model |
| `marker_label_factor` | | a row factor to name markers by (e.g. `c.antigen`) when the FlowSOM step's first row factor is the FCS channel; project it on rows too |
| `population_palette` | `tab20` | colours for the sorted populations, as Tercen's Plot operator assigns them (`tab20`, `tab10`, `Palette-1`, `Set1`, `Dark2`) |
| `group_palette` | `Palette-1` | colours for the sorted groups |
| `model_column` | | which model, when several FlowSOM steps are upstream (e.g. `ds5.flowsom_model`); empty = the nearest |
| `tree_layout` | `unweighted` | `unweighted`: every edge one length, as Python FlowSOM (and run 9) draws it; `weighted`: R `BuildMST`'s layout, edge length = code distance, where a population split over many similar nodes packs into a knot. The tree's edges are FlowSOM's either way |
| `max_node_size` | 1 | FlowSOM's `maxNodeSize` |
| `plot_width`, `plot_height`, `dpi` | 9, 7, 150 | image size in inches and resolution |

## Tests

`ci/smoke.R` renders every plot from the model `flowsom_operator` 2.2.0 writes on its own
platform test (synthetic data), checks the codes are read node by node, and that the tree is
the minimum spanning tree of those codes (1.0.x read the codes scrambled and drew a wrong tree); CI runs it in
the built image.

## Licence

GPL-2.0-or-later (it calls FlowSOM, GPL ≥ 2). See `LICENSE`.

The smoke test lives in `ci/`, not `tests/`: Tercen's installer reads every `tests/*.json` as a
platform test (1.0.0 failed on the model fixture) and fails a `tests/` folder with none (1.0.1).
`tests/test.json` (since 1.0.3) is the platform test: the grid view on a 3x3 model from
flowsom_operator 2.2 on synthetic data; only the grid, because it does not depend on igraph's
floating-point tree layout, so the PNG is byte-stable. `ci/gen_test_expected.R` regenerates it.
