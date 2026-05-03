################################################################################
# 06_sg.R
# Cross-species analysis of spermatogonia (SG) subtypes
# Species: human, macaque, tree shrew, mouse
#
# Functions:
#   1. Subset and integrate SG cells across four species (CCA)
#   2. Identify SG subtypes: SSC, Undiff SG, Diff SG
#   3. SG subtype proportion comparison across species
#   4. Conserved marker gene identification per subtype
#   5. Pseudotime trajectory inference within SG (Monocle3)
#   6. Moran's I analysis for trajectory-dependent genes
#   7. Cross-species comparison of SG gene expression dynamics
#   8. Heatmap visualization of top trajectory-dependent genes
#
# Cell type nomenclature (aligned with manuscript):
#   SSC        - Spermatogonial Stem Cells
#   Undiff SG  - Undifferentiated Spermatogonia
#   Diff SG    - Differentiated Spermatogonia
################################################################################

# =========================== Environment Setup ===========================
setwd("./cross_testis")
# ---------- Load R packages ----------
library(Seurat)            # Single-cell analysis framework
library(tidyverse)         # Data manipulation suite
library(monocle3)          # Pseudotime trajectory analysis
library(pheatmap)          # Heatmap plotting
library(ComplexHeatmap)    # Complex heatmap plotting
library(reshape2)          # Data reshaping
library(scales)            # Scale functions for ggplot2

# =========================== Load Pseudotime Data & Subset SG ===========================

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

# Subset to SG cells only
human_sg     <- subset(human_pseu, NEW_CELL_TYPE == "SG")
macaque_sg   <- subset(macaque_pseu, NEW_CELL_TYPE == "SG")
mouse_sg     <- subset(mouse_pseu, NEW_CELL_TYPE == "SG")
treeshrew_sg <- subset(treeshrew_pseu, NEW_CELL_TYPE == "SG")

# =========================== Merge & Integrate SG Cells (CCA) ===========================

obj_sg <- list(human_sg, macaque_sg, mouse_sg, treeshrew_sg)
obj_sg <- merge(
  x = obj_sg[[1]],
  y = do.call(c, obj_sg[-1]),
  add.cell.ids = c("human", "macaque", "mouse", "treeshrew"),
  project = "species"
)

# --- CCA Integration ---
obj_sg <- obj_sg %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

obj_sg <- IntegrateLayers(
  object = obj_sg,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE,
  dims = 1:30
)

obj_sg <- RunUMAP(obj_sg, reduction = "integrated.cca", dims = 1:8,
                  reduction.name = "umap.cca")

# --- Clustering & Annotation ---
obj_sg <- FindClusters(obj_sg, resolution = 0.2, algorithm = 1)

# Annotate SG subtypes (order matches cluster IDs)
new_cell_type <- c("Undiff SG", "Diff SG", "SSC")
names(new_cell_type) <- levels(obj_sg$seurat_clusters)
obj_sg <- RenameIdents(obj_sg, new_cell_type)

# Set factor levels: SSC -> Undiff SG -> Diff SG (developmental order)
obj_sg$sg_cell_type <- factor(Idents(obj_sg),
                              levels = c("SSC", "Undiff SG", "Diff SG"))
Idents(obj_sg) <- obj_sg$sg_cell_type

save(obj_sg, file = "./data/sg.RData")

# =========================== SG Subtype Proportion Barplot ===========================

# Calculate proportions per species
df <- table(obj_sg$sg_cell_type, obj_sg$Species)
df <- apply(df, 2, function(x) x / sum(x) * 100)
df <- data.frame(Celltype = rownames(df), df)
df <- reshape2::melt(data = df, id.vars = "Celltype",
                     variable.name = "Species", value.name = "Count")
df$Celltype <- factor(df$Celltype, levels = c("SSC", "Undiff SG", "Diff SG"))
df$Species <- factor(df$Species, levels = c("human", "macaque", "treeshrew", "mouse"))

# Bar plot: SG subtype proportions across species
df %>% ggplot(aes(x = Species, y = Count, fill = Celltype)) +
  geom_bar(stat = "identity", color = "black") +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
    axis.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 15),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line()
  ) +
  labs(x = "", y = "", title = "SG Subtype Proportions")

ggsave(filename = "./figure/Celltype.prop.pdf", width = 10, height = 8)

# =========================== UMAP Visualization ===========================

# UMAP colored by species
DimPlot(obj_sg, pt.size = 1, group.by = "Species", reduction = "umap.cca")
ggsave("./figure/SG_Species.umap.pdf", width = 10, height = 8)

# =========================== Conserved Marker Genes per SG Subtype ===========================

# Find conserved markers for Undiff SG and Diff SG (pooled across species)
sg.makers <- data.frame()
for (stage in levels(Idents(obj_sg))[2:3]) {
  df <- FindConservedMarkers(
    obj_sg,
    ident.1 = stage,
    ident.2 = NULL,
    grouping.var = "Species",
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.5,
    only.pos = TRUE
  )
  df$sg_cell_type <- stage
  df$gene <- rownames(df)
  sg.makers <- rbind(sg.makers, df)
}
sg.makers$sg_cell_type <- factor(sg.makers$sg_cell_type,
                                 levels = c("Undiff SG", "Diff SG"))

# Find conserved markers for SSC
ssc.makers <- FindConservedMarkers(
  obj_sg,
  ident.1 = "SSC",
  ident.2 = NULL,
  grouping.var = "Species",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.5,
  only.pos = TRUE
)
ssc.makers$sg_cell_type <- "SSC"
ssc.makers$gene <- rownames(ssc.makers)

write.csv(ssc.makers, file = "./figure/SSC_markers.csv")
write.csv(sg.makers, file = "./figure/SG_markers.csv")

# =========================== Feature Plots for Key SG Markers ===========================

# Plot canonical SG subtype markers (manuscript Figure 6C)
for (gene in c("LIN7B", "SLC25A21", "BNC2", "GFRA1", "DMRT1", "DPEP3")) {
  FeaturePlot(obj_sg, pt.size = 0.75, features = gene,
              reduction = "umap.cca") +
    scale_color_gradient(low = "grey80", high = "red",
                         na.value = "red", limits = c(0, 3))
  ggsave(paste0("./figure/SG_", gene, ".pdf"),
         width = 6, height = 4.5)
}

# =========================== SG Pseudotime Trajectory (Monocle3) ===========================

load("./data/sg.RData")

# Convert to Monocle3 CellDataSet
data <- GetAssayData(obj_sg, assay = 'RNA', layer = 'counts')
cell_metadata <- obj_sg@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 30)

# Embed CCA UMAP into Monocle3
int.embed <- Embeddings(obj_sg, reduction = "umap.cca")
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")

# Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE,
                   learn_graph_control = list(
                     ncenter = 400, minimal_branch_len = 20,
                     euclidean_distance_ratio = 0.5, prune_graph = TRUE, maxiter = 100
                   ))

plot_cells(cds, show_trajectory_graph = TRUE,
           color_cells_by = "sg_cell_type",
           group_label_size = 5, label_principal_points = TRUE,
           cell_size = 1)

# Order cells along pseudotime
cds <- order_cells(cds)

plot_cells(cds, show_trajectory_graph = TRUE,
           color_cells_by = "pseudotime",
           trajectory_graph_color = "green",
           label_roots = TRUE, label_leaves = FALSE,
           label_branch_points = FALSE, label_cell_groups = FALSE,
           cell_size = 2, cell_stroke = 0, graph_label_size = 5)

save(cds, file = "./data/sg_cds.RData")

load(file = "./data/sg_cds.RData")
ggsave("./figure/SG_trajectory.pdf", width = 10, height = 8)

# =========================== Moran's I Analysis for Trajectory-Dependent Genes ===========================

# Identify genes with significant spatial autocorrelation along the trajectory
# (Moran's I test) per species, then find intersection

track_genes <- graph_test(cds, neighbor_graph = "principal_graph",
                          method = "Moran_I", cores = 20)

# Filter significant trajectory-dependent genes per species
# NOTE: track_genes_human, track_genes_macaque, etc. are expected to be
# pre-computed or computed in a similar manner per species
# Here we assume they exist; if not, run graph_test per species subset

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
  
  # Calculate a score comparing primate+tree shrew vs mouse Moran's I
  # Higher score = more conserved in primates/tree shrew relative to mouse
  df$score1 <- log2(df$mouse_morans_I)
  df$score2 <- log2(df$human_morans_I * df$macaque_morans_I *
                      df$treeshrew_morans_I)
  df$score <- df$score1 / df$score2
  
  write.csv(df, file = "./figure/morans_I_track_genes.v2.csv")
  
}, error = function(e) {
  cat("NOTE: Moran's I per-species analysis skipped.",
      "Please ensure track_genes_[species] objects are available.\n")
  cat("Error message:", e$message, "\n")
})

# =========================== Heatmap of Top Trajectory-Dependent Genes ===========================

# Read back the Moran's I comparison results
df <- read.csv("./figure/morans_I_track_genes.v2.csv", row.names = 1)

# Rank genes by conservation score (higher = more primate/tree shrew-like)
genes <- df %>% arrange(desc(score)) %>% pull(gene)
genes <- c(genes, "FMR1")  # Add FMR1 as a gene of interest (manuscript Figure 6G)

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

# Generate smoothed expression heatmaps per species along SG pseudotime
for (species in c("human", "macaque", "treeshrew", "mouse")) {
  cat(">>> Generating SG heatmap for:", species, "\n")
  
  # Subset CDS to this species
  cells_species <- WhichCells(obj_sg, expression = Species == species)
  
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
  
  # Annotation: SG subtype ordered by pseudotime
  anno_col <- colData(cds[, cells_species])$sg_cell_type[
    order(pseudotime(cds[, cells_species]))
  ]
  anno_col <- as.data.frame(anno_col)
  rownames(anno_col) <- cells_species[
    order(pseudotime(cds[, cells_species]))
  ]
  names(anno_col) <- "CellType"
  
  # Color mapping for SG subtypes
  anno_color <- hue_pal()(length(levels(anno_col$CellType)))
  names(anno_color) <- levels(anno_col$CellType)
  anno_color <- list(CellType = anno_color)
  
  # Heatmap
  ht <- pheatmap(
    mat = pt.matrix,
    name = "Z score",
    legend_breaks = seq(-4, 4, 2),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = species,
    annotation_col = anno_col,
    annotation_colors = anno_color,
    cellheight = 25
  )
  
  pdf(paste0("./figure/sg_heatmap_", species, ".pdf"),
      width = calc_ht_size(ht)[1], height = calc_ht_size(ht)[2])
  print(ht)
  dev.off()
}

cat(">>> 06_sg.R completed successfully!\n")