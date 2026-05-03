################################################################################
# 07_scsd.R
# Cross-species analysis of meiotic and post-meiotic germ cells
# (Spermatocytes and Spermatids) across four species
# Species: human, macaque, tree shrew, mouse
#
# Functions:
#   1. Subset and integrate SC + SD cells across four species (CCA)
#   2. Identify fine-grained subtypes:
#      SC: Lept. SC, Zyg. SC, Pach. SC, Dipl. SC
#      SD: Early SD (round), Late SD (elongating)
#   3. Conserved marker gene identification per subtype
#   4. Pseudotime trajectory inference (Monocle3)
#   5. Moran's I analysis for trajectory-dependent genes
#   6. Cross-species comparison of SC/SD gene expression dynamics
#   7. Marker gene heatmap across species and subtypes
#
# Cell type nomenclature (aligned with manuscript):
#   Lept. SC   - Leptotene Spermatocytes
#   Zyg. SC    - Zygotene Spermatocytes
#   Pach. SC   - Pachytene Spermatocytes
#   Dipl. SC   - Diplotene Spermatocytes
#   Early SD   - Early stage Spermatids (round spermatids)
#   Late SD    - Late stage Spermatids (elongating spermatids)
#
# NOTE: In the integrated object, fine subtypes use short names
#       (Lept. SC, Zyg. SC, etc.) while the manuscript cell type column
#       uses the broader categories (SC, Early SD, Late SD).
################################################################################

# =========================== Environment Setup ===========================
setwd("./cross_testis")
# ---------- Load R packages ----------
library(Seurat)            # Single-cell analysis framework
library(tidyverse)         # Data manipulation suite
library(monocle3)          # Pseudotime trajectory analysis
library(pheatmap)          # Heatmap plotting
library(ComplexHeatmap)    # Complex heatmap plotting
library(scales)            # Scale functions for ggplot2

# =========================== Load Pseudotime Data & Subset SC + SD ===========================

load("./data/species.pseudo.RData")

# Add species metadata
human_pseu$Species     <- "human"
macaque_pseu$Species   <- "macaque"
mouse_pseu$Species     <- "mouse"
treeshrew_pseu$Species <- "treeshrew"

# Convert all gene names to uppercase for cross-species consistency
rownames(human_pseu@assays$RNA@features)     <- toupper(rownames(human_pseu@assays$RNA@features))
rownames(macaque_pseu@assays$RNA@features)   <- toupper(rownames(macaque_pseu@assays$RNA@features))
rownames(mouse_pseu@assays$RNA@features)     <- toupper(rownames(mouse_pseu@assays$RNA@features))
rownames(treeshrew_pseu@assays$RNA@features) <- toupper(rownames(treeshrew_pseu@assays$RNA@features))

# Subset to SC + SD cells (meiotic and post-meiotic germ cells)
# Manuscript cell types: SC, Early SD (round SD), Late SD (elongating SD)
human_scsd     <- subset(human_pseu, NEW_CELL_TYPE %in% c("SC", "Early SD", "Late SD"))
macaque_scsd   <- subset(macaque_pseu, NEW_CELL_TYPE %in% c("SC", "Early SD", "Late SD"))
mouse_scsd     <- subset(mouse_pseu, NEW_CELL_TYPE %in% c("SC", "Early SD", "Late SD"))
treeshrew_scsd <- subset(treeshrew_pseu, NEW_CELL_TYPE %in% c("SC", "Early SD", "Late SD"))

# =========================== Merge & Integrate SC+SD Cells (CCA) ===========================

obj_scsd <- list(human_scsd, macaque_scsd, mouse_scsd, treeshrew_scsd)
obj_scsd <- merge(
  x = obj_scsd[[1]],
  y = do.call(c, obj_scsd[-1]),
  add.cell.ids = c("human", "macaque", "mouse", "treeshrew"),
  project = "species"
)

# --- CCA Integration ---
obj_scsd <- obj_scsd %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

obj_scsd <- IntegrateLayers(
  object = obj_scsd,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE,
  dims = 1:20
)

obj_scsd <- RunUMAP(obj_scsd, reduction = "integrated.cca", dims = 1:15,
                    reduction.name = "umap.cca")

# --- Clustering & Fine-Grained Annotation ---
obj_scsd <- FindClusters(obj_scsd, resolution = 0.21, algorithm = 1)

# Annotate fine subtypes (order matches cluster IDs from Seurat)
# These represent the developmental continuum: SC stages -> SD stages
new_cell_type <- c(
  "Early SD",       # Cluster 0: Early round spermatids
  "Late SD",         # Cluster 1: Elongating spermatids
  "Zyg. SC",         # Cluster 2: Zygotene spermatocytes
  "Early SD",        # Cluster 3: Late round spermatids
  "Pach. SC",        # Cluster 4: Pachytene spermatocytes
  "Dipl. SC",        # Cluster 5: Diplotene spermatocytes
  "Lept. SC"         # Cluster 6: Leptotene spermatocytes
)
names(new_cell_type) <- levels(obj_scsd$seurat_clusters)
obj_scsd <- RenameIdents(obj_scsd, new_cell_type)

# Set factor levels in developmental order
obj_scsd$scsd_cell_type <- factor(
  Idents(obj_scsd),
  levels = c("Lept. SC", "Zyg. SC", "Pach. SC", "Dipl. SC",
             "Early SD", "Late SD")
)
Idents(obj_scsd) <- obj_scsd$scsd_cell_type

save(obj_scsd, file = "./data/scsd.RData")

# =========================== Conserved Marker Genes per SC/SD Subtype ===========================

scsd.makers <- data.frame()
for (stage in levels(Idents(obj_scsd))) {
  df <- FindConservedMarkers(
    obj_scsd,
    ident.1 = stage,
    ident.2 = NULL,
    grouping.var = "Species",
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.5,
    only.pos = TRUE
  )
  df$scsd_cell_type <- stage
  df$gene <- rownames(df)
  scsd.makers <- rbind(scsd.makers, df)
}

write.csv(scsd.makers, file = "./figure/SCSD_markers.csv")

# =========================== SC+SD Pseudotime Trajectory (Monocle3) ===========================

load("./data/scsd.RData")

# Convert to Monocle3 CellDataSet
data <- GetAssayData(obj_scsd, assay = 'RNA', layer = 'counts')
cell_metadata <- obj_scsd@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 30)

# Embed CCA UMAP into Monocle3
int.embed <- Embeddings(obj_scsd, reduction = "umap.cca")
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")

# Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE,
  learn_graph_control = list(
    ncenter = 250, minimal_branch_len = 10,
    euclidean_distance_ratio = 0.2, prune_graph = TRUE, maxiter = 50
  ))

plot_cells(cds, show_trajectory_graph = TRUE,
           color_cells_by = "scsd_cell_type",
           group_label_size = 5, label_principal_points = TRUE)

# Order cells along pseudotime
cds <- order_cells(cds)

plot_cells(cds, show_trajectory_graph = TRUE,
           color_cells_by = "pseudotime",
           trajectory_graph_color = "green",
           label_roots = TRUE, label_leaves = FALSE,
           label_branch_points = FALSE, label_cell_groups = FALSE,
           cell_size = 2, cell_stroke = 0, graph_label_size = 5)

ggsave("./figure/figurescsd_trajectory.pdf", width = 10, height = 8)
save(cds, file = "./data/scsd_cds.RData")

# =========================== Moran's I Analysis for Trajectory-Dependent Genes ===========================

load("./data/scsd_cds.RData")

# Global Moran's I on the integrated CDS
track_genes <- graph_test(cds, neighbor_graph = "principal_graph",
                          method = "Moran_I", cores = 20)

# Top trajectory-dependent genes (global)
track_sig_genes <- track_genes %>%
  filter(q_value < 1e-5) %>%
  top_n(20, morans_I) %>%
  pull(gene_short_name)

# --- Per-species Moran's I analysis ---
# NOTE: track_genes_human, track_genes_macaque, etc. are expected to be
# pre-computed per species. If not available, run graph_test on each species subset.

tryCatch({
  track_genes_human %>%
    filter(q_value < 0.05) %>%
    mutate(Z_morans_I = morans_I) -> track_genes_human_sig

  track_genes_macaque %>%
    filter(q_value < 0.05) %>%
    mutate(Z_morans_I = morans_I) -> track_genes_macaque_sig

  track_genes_mouse %>%
    filter(q_value < 0.05) %>%
    mutate(Z_morans_I = morans_I) -> track_genes_mouse_sig

  track_genes_treeshrew %>%
    filter(q_value < 0.05) %>%
    mutate(Z_morans_I = morans_I) -> track_genes_treeshrew_sig

  # Find genes significant in all four species
  track_sig_genes <- Reduce(intersect, list(
    track_genes_human_sig$gene_short_name,
    track_genes_macaque_sig$gene_short_name,
    track_genes_mouse_sig$gene_short_name,
    track_genes_treeshrew_sig$gene_short_name
  ))

  # Build comparison data frame
  df <- data.frame(
    gene               = track_sig_genes,
    human_morans_I     = track_genes_human_sig[track_sig_genes, "Z_morans_I"],
    macaque_morans_I   = track_genes_macaque_sig[track_sig_genes, "Z_morans_I"],
    mouse_morans_I     = track_genes_mouse_sig[track_sig_genes, "Z_morans_I"],
    treeshrew_morans_I = track_genes_treeshrew_sig[track_sig_genes, "Z_morans_I"]
  )

  # Calculate conservation score:
  # Higher score = more divergent in mouse relative to primates + tree shrew
  df$score1 <- log2(df$mouse_morans_I)
  df$score2 <- log2(df$human_morans_I * df$macaque_morans_I *
                    df$treeshrew_morans_I)
  df$score <- df$score1 / df$score2

  write.csv(df, file = "./figure/figuremorans_I_track_genes.v2.csv")

}, error = function(e) {
  cat("NOTE: Moran's I per-species analysis skipped.",
      "Please ensure track_genes_[species] objects are available.\n")
  cat("Error message:", e$message, "\n")
})

# =========================== Heatmap of Top Trajectory-Dependent Genes ===========================

df <- read.csv("./figure/figuremorans_I_track_genes.v2.csv", row.names = 1)

# Select top genes by conservation score
genes <- df %>%
  arrange(desc(score)) %>%
  top_n(23) %>%
  pull(gene)

load("./data/scsd_cds.RData")
genes <- c(genes, "SYCP3")  # Add SYCP3 as key meiotic marker

# Helper function to calculate heatmap size
calc_ht_size <- function(ht, unit = "inch") {
  pdf(NULL)
  ht <- draw(ht)
  w <- ComplexHeatmap:::width(ht)
  w <- grid::convertX(w, unit, valueOnly = TRUE)
  h <- ComplexHeatmap:::height(ht)
  h <- grid::convertY(h, unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}

# Generate smoothed expression heatmaps per species along SC+SD pseudotime
for (species in c("human", "macaque", "mouse", "treeshrew")) {
  cat(">>> Generating SC+SD heatmap for:", species, "\n")

  # Subset CDS to this species
  cells_species <- WhichCells(obj_scsd, expression = Species == species)

  # Extract expression matrix, order by pseudotime
  pt.matrix <- exprs(cds[, cells_species])[
    match(genes, rownames(rowData(cds))),
    order(pseudotime(cds[, cells_species]))
  ]

  # Smooth expression along pseudotime
  pt.matrix <- t(apply(pt.matrix, 1, function(x) {
    smooth.spline(x, df = 10)$y
  }))

  # Z-score normalization per gene
  pt.matrix <- t(apply(pt.matrix, 1, function(x) {
    (x - mean(x)) / sd(x)
  }))

  # Annotation: SC/SD subtype ordered by pseudotime
  anno_col <- colData(cds[, cells_species])$scsd_cell_type[
    order(pseudotime(cds[, cells_species]))
  ]
  anno_col <- as.data.frame(anno_col)
  rownames(anno_col) <- cells_species[
    order(pseudotime(cds[, cells_species]))
  ]
  names(anno_col) <- "CellType"

  # Color mapping
  anno_color <- hue_pal()(length(levels(anno_col$CellType)))
  names(anno_color) <- levels(anno_col$CellType)
  anno_color <- list(CellType = anno_color)

  # Heatmap
  ht <- pheatmap(
    mat = pt.matrix,
    name = "Z score",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = species,
    annotation_col = anno_col,
    annotation_colors = anno_color,
    cellheight = 5
  )

  pdf(paste0("./figure/figurescsd_heatmap_", species, ".pdf"),
      width = calc_ht_size(ht)[1], height = calc_ht_size(ht)[2])
  print(ht)
  dev.off()
}

# =========================== Marker Gene Heatmap (Cross-Species) ===========================
# Heatmap showing canonical SC/SD marker expression across species and subtypes
# (manuscript Figure 7C)

markers <- read.csv("./marker.csv", header = FALSE)
markers <- markers[markers$V2 %in% levels(Idents(obj_scsd)), ]

# Extract scaled expression data
mat <- FetchData(obj_scsd, vars = markers$V1, layer = "data")
mat <- mat[, markers$V1] %>%
  scale() %>%
  t() %>%
  as.matrix()

# Cell metadata for annotation
info <- FetchData(obj_scsd, vars = c("Species", "scsd_cell_type"), layer = "data")

# Color palettes
cols1 <- hue_pal(c = 60, l = 50)(length(unique(info$Species)))
names(cols1) <- unique(info$Species)
cols2 <- hue_pal()(length(unique(info$scsd_cell_type)))
names(cols2) <- levels(info$scsd_cell_type)

# Top annotation
topanno <- HeatmapAnnotation(
  df = info,
  col = list(Species = cols1, scsd_cell_type = cols2),
  which = "column",
  show_annotation_name = FALSE
)

# Create ordered cluster factor for column split
# Order: human -> macaque -> treeshrew -> mouse, within each species by developmental stage
info$cluster <- paste(info$Species, info$scsd_cell_type, sep = "-")
info$cluster <- factor(info$cluster, levels = as.character(outer(
  c("human", "macaque", "treeshrew", "mouse"),
  c("Lept. SC", "Zyg. SC", "Pach. SC", "Dipl. SC", "Early SD", "Late SD"),
  paste, sep = "-"
)))

# Row split by cell type
markers$V2 <- factor(markers$V2,
  levels = c("Lept. SC", "Zyg. SC", "Pach. SC", "Dipl. SC", "Early SD", "Late SD"))

# ComplexHeatmap
ht <- Heatmap(
  matrix = mat,
  use_raster = FALSE,
  name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_split = markers$V2,
  column_split = info$cluster,
  row_gap = unit(1, "mm"),
  column_gap = unit(1, "mm"),
  column_title = "Heatmap of SC/SD Cell Markers",
  row_title_rot = 90,
  top_annotation = topanno,
  height = unit(10, "inch"),
  width = unit(8, "inch")
)

size <- calc_ht_size(ht)
pdf(paste0("./figure/figurescsd.heatmap.pdf"), width = size[1], height = size[2])
print(ht)
dev.off()

cat(">>> 07_scsd.R completed successfully!\n")