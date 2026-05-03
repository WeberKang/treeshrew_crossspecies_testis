################################################################################
# 03_cross_species_integration.R
# Cross-species transcriptomic integration and comparative analysis
# Species: human, macaque, tree shrew, mouse
#
# Functions:
#   1. Gene intersection across species & PCA on pseudo-bulk transcriptomes
#   2. Pearson correlation analysis of cell-type transcriptomes
#   3. Cell type proportion comparison across species
#   4. Species-specific highly expressed gene identification
#   5. Venn diagram of shared/testis-specific genes
#   6. Pseudotime trajectory inference (Monocle3)
#   7. Cross-species trajectory alignment (cellAlign)
#   8. CCA integration & LISI evaluation (Seurat v5)
#
# Cell type nomenclature (aligned with manuscript):
#   SG        - Spermatogonia
#   SC        - Spermatocytes
#   Early SD  - Early stage Spermatids (round spermatids, rSD in code)
#   Late SD   - Late stage Spermatids (elongating spermatids, eSD in code)
#   Sertoli   - Sertoli cells
#   OSC       - Other Somatic Cells
#
# NOTE: In this script, cell types in variable names and factor levels
#       use manuscript nomenclature (Early SD, Late SD, OSC).
#       The objects loaded from previous steps may still have rSD/eSD
#       in their NEW_CELL_TYPE column and will be renamed here.
################################################################################

# =========================== Environment Setup ===========================

setwd("./cross_testis")

# ---------- Load R packages ----------
library(data.table)        # Fast data reading
library(Seurat)            # Single-cell analysis framework (v5)
library(harmony)           # Harmony integration
library(biomaRt)           # Ensembl gene annotation
library(tidyverse)         # Data manipulation suite
library(Matrix)            # Sparse matrix operations
library(ComplexHeatmap)    # Complex heatmap plotting
library(scales)            # Scale functions for ggplot2
library(ggVennDiagram)     # Venn diagram plotting
library(rtracklayer)       # GTF/GFF file I/O
library(RColorBrewer)      # Color palettes
library(ggrepel)           # Repel labels for ggplot2
library(patchwork)         # Plot composition
library(monocle3)          # Pseudotime trajectory analysis
library(cellAlign)         # Cross-species trajectory alignment
library(lisi)              # LISI score for integration evaluation
library(reshape2)          # Data reshaping (melt)

# =========================== Load Processed Species Data ===========================

load("./data/species.RData")

# =========================== Gene Intersection & HVG Filtering ===========================

# Convert mouse gene names to uppercase for cross-species consistency
rownames(mouse@assays$RNA@features) <- toupper(rownames(mouse@assays$RNA@features))

# Find genes shared across all four species
gene_inter <- Reduce(intersect, list(
  Features(human), Features(macaque),
  Features(mouse), Features(treeshrew)
))

# Subset each species to shared genes only
human     <- human[Features(human) %in% gene_inter, ]
macaque   <- macaque[Features(macaque) %in% gene_inter, ]
mouse     <- mouse[Features(mouse) %in% gene_inter, ]
treeshrew <- treeshrew[Features(treeshrew) %in% gene_inter, ]

# Check gene counts after filtering
lapply(list(human, macaque, mouse, treeshrew), function(x) dim(x)[1])

# Add species metadata
human@meta.data["Species"]     <- "human"
macaque@meta.data["Species"]   <- "macaque"
mouse@meta.data["Species"]     <- "mouse"
treeshrew@meta.data["Species"] <- "treeshrew"

save(human, macaque, mouse, treeshrew, file = "./data/cross.RData")

# =========================== PCA on Pseudo-Bulk Transcriptomes ===========================

load("./data/cross.RData")

# --- Identify highly variable genes (HVGs) per species ---
human_hvg     <- human %>% NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>% VariableFeatures() %>% as.character()
macaque_hvg   <- macaque %>% NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>% VariableFeatures() %>% as.character()
mouse_hvg     <- mouse %>% NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>% VariableFeatures() %>% as.character()
treeshrew_hvg <- treeshrew %>% NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>% VariableFeatures() %>% as.character()

# Union of all HVGs
cross_union <- Reduce(union, list(human_hvg, macaque_hvg, mouse_hvg, treeshrew_hvg))

# Subset to HVG union
human     <- human[Features(human) %in% cross_union, ]
macaque   <- macaque[Features(macaque) %in% cross_union, ]
mouse     <- mouse[Features(mouse) %in% cross_union, ]
treeshrew <- treeshrew[Features(treeshrew) %in% cross_union, ]

lapply(list(human, macaque, mouse, treeshrew), function(x) dim(x)[1])

# Normalize
human     <- NormalizeData(human)
macaque   <- NormalizeData(macaque)
mouse     <- NormalizeData(mouse)
treeshrew <- NormalizeData(treeshrew)

# --- Calculate pseudo-bulk average expression per species-celltype combination ---
Idents(human) <- paste(human@meta.data$orig.ident,
                       human@meta.data$Species,
                       human@meta.data$NEW_CELL_TYPE, sep = "-")
human_avg <- as.matrix(AverageExpression(human, assays = "RNA", layer = "data")$RNA)

Idents(macaque) <- paste(macaque@meta.data$orig.ident,
                         macaque@meta.data$Species,
                         macaque@meta.data$NEW_CELL_TYPE, sep = "-")
macaque_avg <- as.matrix(AverageExpression(macaque, assays = "RNA", layer = "data")$RNA)

Idents(mouse) <- paste(mouse@meta.data$orig.ident,
                       mouse@meta.data$Species,
                       mouse@meta.data$NEW_CELL_TYPE, sep = "-")
mouse_avg <- as.matrix(AverageExpression(mouse, assays = "RNA", layer = "data")$RNA)

Idents(treeshrew) <- paste(treeshrew@meta.data$orig.ident,
                           treeshrew@meta.data$Species,
                           treeshrew@meta.data$NEW_CELL_TYPE, sep = "-")
treeshrew_avg <- as.matrix(AverageExpression(treeshrew, assays = "RNA", layer = "data")$RNA)

# Merge across species
cross_avg <- merge(human_avg, macaque_avg, by = 0, all = TRUE) %>%
  merge(mouse_avg, by.x = "Row.names", by.y = 0, all = TRUE) %>%
  merge(treeshrew_avg, by.x = "Row.names", by.y = 0, all = TRUE)

# Extract group annotation (species and cell type)
group <- names(cross_avg) %>%
  .[!grepl("Row.names", .)] %>%
  strsplit(split = "-") %>%
  sapply(function(x) x[2:3]) %>%
  t()
colnames(group) <- c("species", "celltype")

# Convert to matrix and remove batch effect
cross_avg <- cross_avg %>%
  data.frame(row.names = .[, "Row.names"]) %>%
  .[, -which(names(.) == "Row.names")] %>%
  as.matrix() %>%
  limma::removeBatchEffect(batch = c(rep("batch1", 54), rep("batch2", 12)))

write.csv(cross_avg, file = "./figure/way2.remove.batch.avg.csv")

# --- PCA ---
pca <- prcomp(t(cross_avg), scale = TRUE, center = TRUE)

pca.res <- pca$x %>% as.data.frame()
pca.res <- cbind(pca.res, group)
pca.res$celltype <- factor(pca.res$celltype,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC"))

# Calculate variance explained
pca.var <- pca$sdev^2 %>% as.data.frame()
pca.var$var <- round(pca.var$. / sum(pca.var) * 100, 2)
pca.var$pc <- colnames(pca.res)[1:(ncol(pca.res) - 2)]
pca.var$pc <- as.numeric(gsub(x = pca.var$pc, pattern = "PC", replacement = ""))
pca.var$pc <- factor(pca.var$pc, levels = pca.var$pc)

# PCA scatter plot
ggplot(pca.res, aes(PC1, PC2, color = celltype)) +
  geom_point(size = 5, aes(shape = species)) +
  labs(
    x = paste('PC1(', pca.var$var[1], '%)', sep = ''),
    y = paste('PC2(', pca.var$var[2], '%)', sep = '')
  ) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  ggtitle("Species Celltype PCA") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

write.csv(pca.res, file = "./figure/way2.remove.batch.pca.csv")
ggsave(filename = "./figure/way2.remove.batch.pca.pdf", width = 12, height = 9)

# =========================== Pearson Correlation Heatmap ===========================

# Re-calculate HVGs for pseudo-bulk
human_hvg     <- NormalizeData(human) %>% FindVariableFeatures(selection.method = "vst") %>%
  VariableFeatures(human)
macaque_hvg   <- NormalizeData(macaque) %>% FindVariableFeatures(selection.method = "vst") %>%
  VariableFeatures(macaque)
mouse_hvg     <- NormalizeData(mouse) %>% FindVariableFeatures(selection.method = "vst") %>%
  VariableFeatures(mouse)
treeshrew_hvg <- NormalizeData(treeshrew) %>% FindVariableFeatures(selection.method = "vst") %>%
  VariableFeatures(treeshrew)

cross_union <- Reduce(union, list(human_hvg, macaque_hvg, mouse_hvg, treeshrew_hvg))

# Aggregate expression using HVG union
Idents(human) <- paste(human@meta.data$orig.ident,
                       human@meta.data$Species,
                       human@meta.data$NEW_CELL_TYPE, sep = "-")
human_avg <- as.matrix(AggregateExpression(human, assays = "RNA",
                                           features = cross_union)$RNA)

Idents(macaque) <- paste(macaque@meta.data$orig.ident,
                         macaque@meta.data$Species,
                         macaque@meta.data$NEW_CELL_TYPE, sep = "-")
macaque_avg <- as.matrix(AggregateExpression(macaque, assays = "RNA",
                                             features = cross_union)$RNA)

Idents(mouse) <- paste(mouse@meta.data$orig.ident,
                       mouse@meta.data$Species,
                       mouse@meta.data$NEW_CELL_TYPE, sep = "-")
mouse_avg <- as.matrix(AggregateExpression(mouse, assays = "RNA",
                                           features = cross_union)$RNA)

Idents(treeshrew) <- paste(treeshrew@meta.data$orig.ident,
                           treeshrew@meta.data$Species,
                           treeshrew@meta.data$NEW_CELL_TYPE, sep = "-")
treeshrew_avg <- as.matrix(AggregateExpression(treeshrew, assays = "RNA",
                                               features = cross_union)$RNA)

# Merge and remove batch effect
cross_avg <- merge(human_avg, macaque_avg, by = 0, all = TRUE) %>%
  merge(mouse_avg, by.x = "Row.names", by.y = 0, all = TRUE) %>%
  merge(treeshrew_avg, by.x = "Row.names", by.y = 0, all = TRUE)

group <- names(cross_avg) %>%
  .[!grepl("Row.names", .)] %>%
  strsplit(split = "-") %>%
  sapply(function(x) x[1:3]) %>%
  t() %>% as.data.frame()
colnames(group) <- c("ID", "species", "celltype")

cross_avg <- cross_avg %>%
  data.frame(row.names = .[, "Row.names"]) %>%
  .[, -which(names(.) == "Row.names")] %>%
  as.matrix() %>%
  limma::removeBatchEffect(batch = c(rep("batch1", 36), rep("batch2", 8)))

# Pearson correlation matrix
cor_gather <- cor(t(scale(t(cross_avg))), method = "pearson")
rownames(group) <- rownames(cor_gather)

# Color annotations for heatmap
ID <- hue_pal(c = 60, l = 50)(length(unique(group$ID)))
names(ID) <- unique(group$ID)
species <- hue_pal(h.start = 50, c = 200)(length(unique(group$species)))
names(species) <- c("human", "macaque", "mouse", "treeshrew")
celltype <- hue_pal()(length(unique(group$celltype)) + 2)[1:4]
names(celltype) <- c("SG", "SC", "Early SD", "Late SD")

# Column annotation
ha <- HeatmapAnnotation(
  df = group[, c(3, 2, 1)],
  which = "column",
  show_annotation_name = TRUE,
  col = list(ID = ID, species = species, celltype = celltype),
  border = TRUE,
  annotation_name_side = "left",
  show_legend = TRUE
)

# Row annotation
ha2 <- HeatmapAnnotation(
  df = group,
  which = "row",
  show_annotation_name = FALSE,
  col = list(ID = ID, species = species, celltype = celltype),
  border = TRUE,
  show_legend = FALSE
)

# Helper function to calculate heatmap size for PDF output
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

# Plot correlation heatmap
ht <- Heatmap(
  matrix = cor_gather,
  name = "Correlation",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  top_annotation = ha,
  right_annotation = ha2,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_column_dend = FALSE,
  row_dend_width = unit(20, "mm"),
  column_names_side = "top",
  column_names_rot = 45,
  height = unit(5, "mm") * nrow(cor_gather),
  width = unit(5, "mm") * ncol(cor_gather),
  col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
)

size <- calc_ht_size(ht)
pdf(paste0("./figure/Pearson.Correlation.pdf"), width = size[1], height = size[2])
print(ht)
dev.off()

# =========================== Cell Type Proportion Barplot ===========================

load("./data/species.RData")

# Set consistent factor levels
human@meta.data$NEW_CELL_TYPE <- factor(
  human@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC"))
macaque@meta.data$NEW_CELL_TYPE <- factor(
  macaque@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC"))
mouse@meta.data$NEW_CELL_TYPE <- factor(
  mouse@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC"))
treeshrew@meta.data$NEW_CELL_TYPE <- factor(
  treeshrew@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC"))

# Calculate cell type proportions (germ cells only)
df <- cbind(
  human     = table(human@meta.data$NEW_CELL_TYPE),
  macaque   = table(macaque@meta.data$NEW_CELL_TYPE),
  mouse     = table(mouse@meta.data$NEW_CELL_TYPE),
  treeshrew = table(treeshrew@meta.data$NEW_CELL_TYPE)
)
df <- df[c("SG", "SC", "Early SD", "Late SD"), ]
df <- apply(df, 2, function(x) x / sum(x) * 100)
df <- data.frame(Celltype = rownames(df), df)
df <- reshape2::melt(data = df, id.vars = "Celltype",
                     variable.name = "Species", value.name = "Count")
df$Celltype <- factor(df$Celltype, levels = c("SG", "SC", "Early SD", "Late SD"))
df$Species <- factor(df$Species, levels = c("human", "macaque", "treeshrew", "mouse"))

# Bar plot
df %>% ggplot(aes(x = Species, y = Count, fill = Celltype)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0(round(Count, 2), "%")),
            position = position_stack(vjust = 0.5), size = 5) +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
    axis.text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 15),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line()
  ) +
  labs(x = "", y = "", title = "Celltype Proportion (Germ Cells)")

ggsave(filename = "./figure/Celltype.germ.pdf", width = 10, height = 8)

# =========================== Species-Specific Highly Expressed Genes ===========================

load("./data/species.RData")

for (species in c("human", "macaque", "mouse", "treeshrew")) {
  obj <- get(species)
  Idents(obj) <- obj$NEW_CELL_TYPE

  avg <- AverageExpression(obj, assays = "RNA", layer = "counts")$RNA
  avg <- LogNormalize(avg) %>% as.data.frame()

  stages <- c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
  alist <- list()

  for (stage in stages) {
    current_stage_value <- avg[, stage]
    other_stages_sum <- rowSums(avg[, stages[!stages %in% stage]]) / 5
    condition <- current_stage_value > 0.25 &
                 current_stage_value / other_stages_sum > 4
    alist[[stage]] <- avg[condition, ]
  }

  assign(paste0(species, "_high"), do.call(rbind, alist))
}

# Plot marker gene heatmaps per species
for (species in c("human", "macaque", "mouse", "treeshrew")) {
  mat <- get(paste0(species, "_high")) %>% as.matrix()

  Celltype <- c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
  cols <- hue_pal()(length(Celltype))
  col_info <- factor(Celltype, levels = Celltype)
  names(cols) <- col_info
  names(col_info) <- Celltype

  column_ha <- HeatmapAnnotation(
    CellType = col_info,
    col = list(CellType = cols),
    border = TRUE,
    which = "column",
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  ht <- Heatmap(
    matrix = mat,
    name = species,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = column_ha,
    show_row_names = FALSE,
    show_column_names = TRUE,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    column_names_side = "top",
    column_names_rot = 45,
    col = circlize::colorRamp2(
      c(0, max(apply(mat, 2, summary)[2, ]), max(apply(mat, 2, summary)[5, ])),
      c("blue", "white", "red")
    ),
    width = ncol(mat) * unit(30, "mm")
  )

  size <- calc_ht_size(ht)
  pdf(paste0("./figure/", species, ".pdf"), width = size[1], height = size[2])
  print(ht)
  dev.off()

  write.csv(get(paste0(species, "_high")), file = paste0("./figure/", species, ".csv"))
}

# =========================== Venn Diagram of Shared Genes ===========================

load("./data/species.RData")

venn <- Venn(list(
  human     = toupper(Features(human)),
  macaque   = toupper(Features(macaque)),
  mouse     = toupper(Features(mouse)),
  treeshrew = toupper(Features(treeshrew))
))

data <- process_data(venn)

ggplot() +
  # Region fill layer
  geom_polygon(aes(X, Y, fill = id, group = id),
               data = venn_setedge(data), show.legend = FALSE, alpha = 0.5) +
  # Set edge layer
  geom_path(aes(X, Y, color = id, group = id),
            data = venn_setedge(data), linewidth = 2, show.legend = FALSE) +
  # Set label layer
  geom_text(aes(X, Y, label = name),
            data = venn_setlabel(data), size = 8) +
  # Region count label layer
  geom_text(aes(X, Y, label = count),
            data = venn_regionlabel(data), size = 8) +
  coord_equal() +
  scale_fill_manual(values = hue_pal()(4)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Cross-species Gene Overlap") +
  scale_x_continuous(expand = c(0.2, 0, 0.2, 0))

# Save Venn data
df <- data$regionData %>% as.data.frame()
df$item <- sapply(df$item, function(x) paste(unlist(x), collapse = "\n"))
write.csv(df, file = paste0("./figure/venn.csv"))

df <- data$regionData$item
names(df) <- data$regionData$name
df <- as.data.frame(sapply(df, "[", i = 1:max(sapply(df, length))))
write.csv(df, file = paste0("./figure/venn.v2.csv"), na = "")

ggsave(filename = paste0("./figure/Cross_gene.venn.pdf"), width = 8, height = 8)

# =========================== Pseudotime Trajectory Inference (Monocle3) ===========================

load("./data/species.RData")

# Subset to germ cells only
Idents(human)     <- human$NEW_CELL_TYPE
Idents(macaque)   <- macaque$NEW_CELL_TYPE
Idents(mouse)     <- mouse$NEW_CELL_TYPE
Idents(treeshrew) <- treeshrew$NEW_CELL_TYPE

human     <- subset(human, idents = c("SG", "SC", "Early SD", "Late SD"))
macaque   <- subset(macaque, idents = c("SG", "SC", "Early SD", "Late SD"))
mouse     <- subset(mouse, idents = c("SG", "SC", "Early SD", "Late SD"))
treeshrew <- subset(treeshrew, idents = c("SG", "SC", "Early SD", "Late SD"))

# Run pseudotime inference for each species
for (species in c("human", "macaque", "mouse", "treeshrew")) {
  cat(">>> Running pseudotime for:", species, "\n")

  obj <- get(species)

  # Identify HVGs
  hvg <- NormalizeData(obj) %>%
    FindVariableFeatures(selection.method = "vst") %>%
    VariableFeatures(obj)
  assign(paste0(species, "_hvg"), hvg)

  # Convert to Monocle3 CellDataSet
  data <- GetAssayData(obj, assay = 'RNA', layer = 'counts')
  cell_metadata <- obj@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)

  cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)

  cds <- preprocess_cds(cds, num_dim = 30)

  # Embed Seurat UMAP into Monocle3
  int.embed <- Embeddings(obj, reduction = "umap")
  cds@int_colData$reducedDims$UMAP <- int.embed

  cds <- cluster_cells(cds)
  plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")

  # Learn trajectory graph (species-specific parameters)
  if (species == "human") {
    cds <- learn_graph(cds, use_partition = FALSE,
      learn_graph_control = list(ncenter = 150, minimal_branch_len = 8,
                                 euclidean_distance_ratio = 0.5, prune_graph = TRUE,
                                 maxiter = 100))
  } else if (species == "macaque") {
    cds <- learn_graph(cds, use_partition = FALSE,
      learn_graph_control = list(ncenter = 65, minimal_branch_len = 6,
                                 euclidean_distance_ratio = 0.5, prune_graph = TRUE,
                                 maxiter = 100))
  } else if (species == "mouse") {
    cds <- learn_graph(cds, use_partition = FALSE,
      learn_graph_control = list(ncenter = 100, minimal_branch_len = 8,
                                 euclidean_distance_ratio = 0.5, prune_graph = TRUE,
                                 maxiter = 100))
  } else if (species == "treeshrew") {
    cds <- learn_graph(cds, use_partition = FALSE,
      learn_graph_control = list(ncenter = 140, minimal_branch_len = 5,
                                 euclidean_distance_ratio = 0.5, prune_graph = TRUE,
                                 maxiter = 100))
  }

  # Order cells along pseudotime
  cds <- order_cells(cds)

  # Plot trajectory
  plot_cells(cds, show_trajectory_graph = TRUE, color_cells_by = "NEW_CELL_TYPE",
             group_label_size = 5, label_principal_points = TRUE)
  plot_cells(cds, show_trajectory_graph = TRUE, color_cells_by = "pseudotime",
             label_cell_groups = FALSE)

  # Store pseudotime (rescaled to 0-1)
  assign(paste0(species, "_pseudotime"), rescale(pseudotime(cds)))

  # Add pseudotime to Seurat metadata
  pseudotime <- get(paste0(species, "_pseudotime"))
  obj@meta.data$pseudotime <- pseudotime
  assign(paste0(species, "_pseu"), obj)
}

# Save pseudotime results
save(list = c(ls(pattern = "_pseudotime"), ls(pattern = "_hvg"), ls(pattern = "_pseu$")),
     file = paste0("./data/species.pseudo.RData"))

# =========================== Cross-Species Trajectory Alignment (cellAlign) ===========================

load("./data/species.pseudo.RData")

# Interpolate expression along pseudotime for each species
for (species in c("human", "macaque", "mouse", "treeshrew")) {
  cat(">>> Interpolating for:", species, "\n")

  obj <- get(paste0(species, "_pseu"))
  data <- GetAssayData(obj, assay = 'RNA', layer = 'counts')
  pseudotime <- get(paste0(species, "_pseudotime"))

  interGlobal <- interWeights(expDataBatch = data, trajCond = pseudotime,
                              winSz = 0.1, numPts = 200)
  assign(x = paste0("interGlobal_", species), value = interGlobal)

  # Scale interpolated data
  interGlobal <- scaleInterpolate(interGlobal)
  assign(x = paste0("interScaledGlobal_", species), value = interGlobal)
}

save(list = c(ls(pattern = "interGlobal_"), ls(pattern = "interScaledGlobal_")),
     file = "./data/species.align.RData")

# =========================== Pairwise Trajectory Alignment ===========================

load("./data/species.align.RData")

# Generate all pairwise species comparisons (lower triangle)
compare_group <- outer(c("human", "macaque", "mouse", "treeshrew"),
                       c("human", "macaque", "mouse", "treeshrew"),
                       paste, sep = "-")
compare_group <- compare_group[lower.tri(compare_group)]

for (pair in compare_group) {
  specie1 <- strsplit(pair, "-")[[1]][1]
  specie2 <- strsplit(pair, "-")[[1]][2]

  # Get shared markers
  sharedMarkers <- Reduce(intersect, lapply(
    paste0(c(specie1, specie2), "_hvg"),
    function(x) toupper(get(x))
  ))

  # Load and subset interpolated data for specie1
  interGlobal1 <- get(paste0("interGlobal_", specie1))
  rownames(interGlobal1$interpolatedVals) <- toupper(rownames(interGlobal1$interpolatedVals))
  rownames(interGlobal1$error) <- toupper(rownames(interGlobal1$error))
  interGlobal1$interpolatedVals <- interGlobal1$interpolatedVals[sharedMarkers, ]
  interGlobal1$error <- interGlobal1$error[sharedMarkers, ]
  interGlobal1 <- scaleInterpolate(interGlobal1)
  interScaledGlobal1 <- interGlobal1

  # Load and subset interpolated data for specie2
  interGlobal2 <- get(paste0("interGlobal_", specie2))
  rownames(interGlobal2$interpolatedVals) <- toupper(rownames(interGlobal2$interpolatedVals))
  rownames(interGlobal2$error) <- toupper(rownames(interGlobal2$error))
  interGlobal2$interpolatedVals <- interGlobal2$interpolatedVals[sharedMarkers, ]
  interGlobal2$error <- interGlobal2$error[sharedMarkers, ]
  interGlobal2 <- scaleInterpolate(interGlobal2)
  interScaledGlobal2 <- interGlobal2

  # Global alignment
  alignment <- globalAlign(
    interScaledGlobal1$scaledData,
    interScaledGlobal2$scaledData,
    scores = list(query = interScaledGlobal1$traj, ref = interScaledGlobal2$traj),
    sigCalc = FALSE, numPerm = 20
  )

  # Alignment cost matrix heatmap
  ht <- Heatmap(
    matrix = alignment$localCostMatrix,
    col = circlize::colorRamp2(
      breaks = c(0, 0.1, 0.2, 0.4, 0.7, 1),
      colors = c(rev(brewer.pal(n = 6, name = "RdYlBu")))
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "Distance",
    width = unit(8, "inch"),
    height = unit(8, "inch"),
    column_title = paste0(specie1, " vs ", specie2),
    column_title_side = "top",
    column_title_gp = gpar(fontsize = 20),
    heatmap_legend_param = list(
      title = "Distance",
      at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0")
    )
  )

  size <- calc_ht_size(ht)
  pdf(paste0("./figure/", pair, ".heatmap.pdf"), width = size[1], height = size[2])
  print(ht)
  dev.off()
}

# =========================== CCA Integration & LISI Evaluation ===========================

load("./data/species.pseudo.RData")

# Convert all gene names to uppercase
rownames(human_pseu@assays$RNA@features)     <- toupper(rownames(human_pseu@assays$RNA@features))
rownames(macaque_pseu@assays$RNA@features)   <- toupper(rownames(macaque_pseu@assays$RNA@features))
rownames(mouse_pseu@assays$RNA@features)     <- toupper(rownames(mouse_pseu@assays$RNA@features))
rownames(treeshrew_pseu@assays$RNA@features) <- toupper(rownames(treeshrew_pseu@assays$RNA@features))

# Take intersection of all genes
gene_inter <- Reduce(intersect, lapply(
  ls(pattern = "_pseu$"), function(x) toupper(Features(get(x)))
))

human_pseu     <- human_pseu[Features(human_pseu) %in% gene_inter, ]
macaque_pseu   <- macaque_pseu[Features(macaque_pseu) %in% gene_inter, ]
mouse_pseu     <- mouse_pseu[Features(mouse_pseu) %in% gene_inter, ]
treeshrew_pseu <- treeshrew_pseu[Features(treeshrew_pseu) %in% gene_inter, ]

# Add species metadata
human_pseu$Species     <- "human"
macaque_pseu$Species   <- "macaque"
mouse_pseu$Species     <- "mouse"
treeshrew_pseu$Species <- "treeshrew"

# Merge all species
data.list <- list(human_pseu, macaque_pseu, mouse_pseu, treeshrew_pseu)
obj <- merge(
  x = data.list[[1]],
  y = do.call(c, data.list[-1]),
  add.cell.ids = c("human", "macaque", "mouse", "treeshrew"),
  project = "species"
)

# --- CCA Integration (Seurat v5) ---
obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:20, reduction.name = "umap")

obj <- IntegrateLayers(
  object = obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE,
  dims = 1:20
)

obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:20,
               reduction.name = "umap.cca")

# Set cell type factor levels
obj$NEW_CELL_TYPE <- factor(obj$NEW_CELL_TYPE,
                            levels = c("SG", "SC", "Early SD", "Late SD"))

# --- LISI Score Calculation ---
pca_coords <- Embeddings(obj, reduction = "pca")[, 1:20]
cca_coords <- Embeddings(obj, reduction = "integrated.cca")[, 1:20]
metadata <- obj@meta.data[, c("Species", "NEW_CELL_TYPE")]

lisi_pca <- compute_lisi(pca_coords, metadata,
                         label_colnames = c("Species", "NEW_CELL_TYPE"))
lisi_cca <- compute_lisi(cca_coords, metadata,
                         label_colnames = c("Species", "NEW_CELL_TYPE"))

# Print summary statistics
cat("\n--- iLISI (species mixing) ---\n")
cat("PCA iLISI mean:", mean(lisi_pca$Species), "\n")
cat("CCA iLISI mean:", mean(lisi_cca$Species), "\n")
cat("\n--- cLISI (cell type purity) ---\n")
cat("PCA cLISI mean:", mean(lisi_pca$NEW_CELL_TYPE), "\n")
cat("CCA cLISI mean:", mean(lisi_cca$NEW_CELL_TYPE), "\n")

# Calculate normalized iLISI
props <- table(obj$Species)
props <- props / sum(props)
max_iLISI <- 1 / sum(props^2)
cat(sprintf("Max theoretical iLISI: %.3f\n", max_iLISI))
cat(sprintf("CCA iLISI / Max iLISI: %.3f\n", median(lisi_cca$Species) / max_iLISI))

# LISI boxplots
plot_df <- data.frame(
  iLISI  = c(lisi_pca$Species, lisi_cca$Species),
  cLISI  = c(lisi_pca$NEW_CELL_TYPE, lisi_cca$NEW_CELL_TYPE),
  Method = factor(rep(c("PCA", "CCA"), each = ncol(obj)), levels = c("PCA", "CCA"))
)

# iLISI plot (species mixing)
p1 <- ggplot(plot_df, aes(x = Method, y = iLISI, fill = Method)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red") +
  theme_minimal() +
  labs(title = "iLISI Scores (Species Mixing)", y = "iLISI Score") +
  theme(legend.position = "none")

# cLISI plot (cell type purity)
p2 <- ggplot(plot_df, aes(x = Method, y = cLISI, fill = Method)) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red") +
  theme_minimal() +
  labs(title = "cLISI Scores (Cell Type Purity)", y = "cLISI Score") +
  theme(legend.position = "none")

p1 | p2
ggsave(filename = "./figure/LISI.pdf", plot = p1 | p2, width = 10, height = 5)

save(obj, file = "./data/species.merge.RData")

# =========================== Pairwise Trajectory Correlation ===========================

sharedMarkers <- Reduce(intersect, lapply(
  ls(pattern = "_pseu$"), function(x) toupper(Features(get(x)))
))

# Upper triangle for correlation heatmaps
compare_group <- outer(c("human", "macaque", "mouse", "treeshrew"),
                       c("human", "macaque", "mouse", "treeshrew"),
                       paste, sep = "-")
compare_group <- compare_group[upper.tri(compare_group)]

for (pair in compare_group) {
  specie1 <- strsplit(pair, "-")[[1]][1]
  specie2 <- strsplit(pair, "-")[[1]][2]

  interScaledGlobal1 <- get(paste0("interScaledGlobal_", specie1))
  interScaledGlobal2 <- get(paste0("interScaledGlobal_", specie2))

  mtx1 <- interScaledGlobal1$scaledData
  mtx2 <- interScaledGlobal2$scaledData
  rownames(mtx1) <- toupper(rownames(mtx1))
  rownames(mtx2) <- toupper(rownames(mtx2))

  # Correlation between pseudotime-binned expression profiles
  mtx <- cor(mtx1[sharedMarkers, ], mtx2[sharedMarkers, ],
             method = "pearson", use = "pairwise.complete.obs")

  ht <- Heatmap(
    matrix = mtx,
    name = "Correlation",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = paste0(specie1, " vs ", specie2)
  )

  size <- calc_ht_size(ht)
  pdf(paste0("./figure/", pair, ".heatmap.pdf"), width = size[1], height = size[2])
  print(ht)
  dev.off()
}

cat(">>> 03_cross_species_integration.R completed successfully!\n")