################################################################################
# 04_wgcna.R
# High-dimensional Weighted Gene Co-expression Network Analysis (hdWGCNA)
# across four species: human, macaque, tree shrew, mouse
#
# Functions:
#   1. Set up hdWGCNA object from integrated cross-species Seurat object
#   2. Construct metacells per cell type and species
#   3. Soft power threshold testing
#   4. Co-expression network construction and module detection
#   5. Module eigengene (hME) calculation and visualization
#   6. Module-trajectory correlation along pseudotime
#   7. Cell cycle score vs module correlation
#   8. Hub gene identification and dot plot visualization
#
# Cell type nomenclature (aligned with manuscript):
#   SG        - Spermatogonia
#   SC        - Spermatocytes
#   Early SD  - Early stage Spermatids (round spermatids, rSD in code)
#   Late SD   - Late stage Spermatids (elongating spermatids, eSD in code)
#   Sertoli   - Sertoli cells
#   OSC       - Other Somatic Cells
################################################################################

# =========================== Environment Setup ===========================

setwd("./cross_testis")

# ---------- Load R packages ----------
library(Seurat)            # Single-cell analysis framework
library(tidyverse)         # Data manipulation suite
library(WGCNA)             # Weighted gene co-expression network analysis
library(hdWGCNA)           # High-dimensional WGCNA for single-cell data
library(cowplot)           # Publication-ready plotting
library(patchwork)         # Plot composition
library(UCell)             # UCell gene signature scoring
library(reshape2)          # Data reshaping

# =========================== Load Integrated Data & Prepare hdWGCNA Object ===========================

load("./data/species.merge.RData")

# Create a fresh Seurat object from the integrated data
# Use raw counts and metadata from the CCA-integrated object
obj <- CreateSeuratObject(
  GetAssayData(obj, layer = "counts"),
  meta.data = obj@meta.data
)

# Transfer UMAP reduction from integrated object
obj[["umap.cca"]] <- CreateDimReducObject(
  embeddings = Embeddings(obj, reduction = "umap.cca") %>% as.matrix(),
  key = "umapcca_",
  assay = "RNA"
)
rm(obj)  # Free memory

# =========================== Setup for hdWGCNA ===========================

set.seed(2024)
enableWGCNAThreads(nThreads = 8)

# Select genes expressed in at least 5% of cells
obj <- SetupForWGCNA(
  obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "obj"
)

# =========================== Construct Metacells ===========================

# Aggregate cells into metacells grouped by cell type and individual to reduce noise
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("NEW_CELL_TYPE", "orig.ident"),
  reduction = 'umap.cca',
  k = 25,
  max_shared = 10,
  ident.group = 'NEW_CELL_TYPE'
)

# Normalize metacell expression matrix
obj <- NormalizeMetacells(obj)

# =========================== Set Expression Data for Network Construction ===========================

# Focus on germ cell types only
obj <- SetDatExpr(
  obj,
  group_name = c("SG", "SC", "Early SD", "Late SD"),
  group.by = 'NEW_CELL_TYPE',
  assay = 'RNA',
  slot = 'data'
)

# =========================== Test Soft Power Thresholds ===========================

obj <- TestSoftPowers(
  obj,
  networkType = 'signed'  # Signed network preserves directionality
)

# Plot soft power results
plot_list <- PlotSoftPowers(obj)
wrap_plots(plot_list, ncol = 2)

# View power table
power_table <- GetPowerTable(obj)
head(power_table)

# =========================== Construct Co-expression Network ===========================

obj <- ConstructNetwork(
  obj,
  tom_name = NULL,        # Do not write TOM to disk
  overwrite_tom = TRUE
)

# Visualize dendrogram
PlotDendrogram(obj, main = 'Cross-Species hdWGCNA Dendrogram')

# =========================== Calculate Module Eigengenes (hMEs) ===========================

# Normalize and scale full dataset for hME calculation
obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = VariableFeatures(obj))

# Compute module eigengenes
obj <- ModuleEigengenes(
  obj,
  group.by.vars = "orig.ident"
)

# Compute module connectivity
obj <- ModuleConnectivity(
  obj,
  harmonized = TRUE
)

# Rename modules to M1, M2, ..., M15
obj <- ResetModuleNames(
  obj,
  _name = "M"
)

# =========================== Module Expression Score (UCell) ===========================

obj <- ModuleExprScore(
  obj,
  n_genes = 25,
  method = 'UCell'
)

# =========================== Cell Cycle Score vs Module Correlation ===========================

# Calculate cell cycle scores
obj <- CellCycleScoring(
  obj,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)

df <- data.frame(
  S.Score   = obj@meta.data$S.Score,
  G2M.Score = obj@meta.data$G2M.Score
)

# Get module scores
modules <- obj@misc[["obj"]][["module_scores"]]

# Calculate correlation between cell cycle scores and each module
cor_matrix <- data.frame(
  Module        = names(modules),
  Cor_S_score   = sapply(modules, function(x) cor(df$S.Score, x, method = "pearson")),
  Cor_G2M_score = sapply(modules, function(x) cor(df$G2M.Score, x, method = "pearson")),
  Pval_S_score  = sapply(modules, function(x) cor.test(df$S.Score, x, method = "pearson")$p.value),
  Pval_G2M_score = sapply(modules, function(x) cor.test(df$G2M.Score, x, method = "pearson")$p.value)
)

# Prepare correlation heatmap data
cor_heatmap <- as.matrix(cor_matrix[, c("Cor_S_score", "Cor_G2M_score")])
rownames(cor_heatmap) <- cor_matrix$Module
colnames(cor_heatmap) <- c("S.Score", "G2M.Score")
cor_heatmap <- t(cor_heatmap)

p_matrix <- as.matrix(cor_matrix[, c("Pval_S_score", "Pval_G2M_score")])
rownames(p_matrix) <- cor_matrix$Module
colnames(p_matrix) <- c("S.Score", "G2M.Score")
p_matrix <- t(p_matrix)

# Helper function: format p-values as significance stars
format_p_stars <- function(p) {
  stars <- character(length(p))
  stars[p < 0.001] <- "***"
  stars[p < 0.01 & p >= 0.001] <- "**"
  stars[p < 0.05 & p >= 0.01] <- "*"
  stars[p >= 0.05] <- "ns"
  stars[p == 0] <- "***"
  return(stars)
}

# Helper function: format p-values in scientific notation
format_p_scientific <- function(p) {
  sapply(p, function(x) {
    if (x == 0) return("< 1e-324")
    if (x < 0.0001) return(sprintf("%.1e", x))
    return(sprintf("%.4f", x))
  })
}

# Reshape for ggplot2
plot_data <- melt(cor_heatmap, varnames = c("Score", "Module"), value.name = "Correlation")
pval_data <- melt(p_matrix, varnames = c("Score", "Module"), value.name = "Pvalue")
plot_data$Pvalue <- pval_data$Pvalue
plot_data$Significance <- format_p_stars(plot_data$Pvalue)
plot_data$Pvalue_label <- format_p_scientific(plot_data$Pvalue)

# Heatmap: Cell cycle scores vs module eigengenes
p <- ggplot(plot_data, aes(x = Module, y = Score, fill = Correlation)) +
  geom_tile(color = "grey50", size = 0.5) +
  geom_text(aes(label = paste0(sprintf("%.2f", Correlation), "\n", Significance)),
            size = 3.5, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limit = c(-0.5, 0.5),
    breaks = seq(-0.4, 0.4, by = 0.2),
    labels = c("-0.4", "-0.2", "0", "0.2", "0.4"),
    name = "Correlation"
  ) +
  theme_minimal() +
  labs(
    title = "Correlation between Cell Cycle Scores and Modules",
    subtitle = expression(italic("*** p < 0.001, ** p < 0.01, * p < 0.05, ns: not significant")),
    x = "Module", y = ""
  ) +
  coord_fixed(ratio = 0.8)

print(p)
ggsave(filename = "./figure/correlation_score_and_modules.pdf", width = 12, height = 6)

# =========================== Module hME Feature Plots ===========================

plot_list <- ModuleFeaturePlot(
  obj,
  reduction = "umap.cca",
  features = 'hMEs',
  order = TRUE
)
wrap_plots(plot_list, ncol = 5)

# =========================== Module Trajectory Plots ===========================

p <- PlotModuleTrajectory(
  obj,
  pseudotime_col = "pseudotime",
  ncol = 5
)
ggsave(p, file = "./figure/module.trajectory.pdf")

# =========================== Save hdWGCNA Object ===========================

write.csv(GetModules(obj), file = "./data/obj.Modules.csv")
save(obj, file = "./data/obj.modules.Rdata")

# =========================== Hub Gene Dot Plots ===========================

load("./data/obj.modules.Rdata")

# ---------- Helper function: Generate cross-species dot plot data ----------
# Creates scaled expression and percent expressed data for hub genes
# across cell types and species
generate_dotplot_data <- function(obj, hub_genes) {
  df <- FetchData(obj, vars = c(hub_genes, "Species", "NEW_CELL_TYPE"))

  data.plots <- list()
  for (species in unique(obj$Species)) {
    df_sub <- df[df$Species == species, ]

    data.plot <- lapply(
      X = unique(x = df_sub$NEW_CELL_TYPE),
      FUN = function(ident) {
        data.use <- df_sub[df_sub$NEW_CELL_TYPE == ident,
                           1:(ncol(x = df_sub) - 2), drop = FALSE]
        avg.exp <- apply(
          X = data.use, MARGIN = 2,
          FUN = function(x) { return(mean(x = expm1(x = x))) }
        )
        pct.exp <- apply(
          X = data.use, MARGIN = 2,
          FUN = PercentAbove, threshold = 0
        )
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
      }
    )
    names(x = data.plot) <- unique(x = df_sub$NEW_CELL_TYPE)

    data.plot <- lapply(
      X = names(x = data.plot),
      FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
      }
    )
    data.plot <- do.call(what = 'rbind', args = data.plot)

    # Scale average expression
    avg.exp.scaled <- sapply(
      X = unique(x = data.plot$features.plot),
      FUN = function(x) {
        data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
        data.use <- scale(x = log1p(data.use))
        return(data.use)
      }
    )

    # Scale percent expressed
    pct.exp.scaled <- sapply(
      X = unique(x = data.plot$features.plot),
      FUN = function(x) {
        data.use <- data.plot[data.plot$features.plot == x, 'pct.exp']
        data.use <- scales::rescale(x = data.use)
        return(data.use)
      }
    )

    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    pct.exp.scaled <- as.vector(x = t(x = pct.exp.scaled))

    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$pct.exp.scaled <- pct.exp.scaled
    data.plot$Species <- species
    data.plots[[species]] <- data.plot
  }

  data.plots <- do.call(what = 'rbind', args = data.plots)
  data.plots$pct.exp <- data.plots$pct.exp * 100
  data.plots$pct.exp.scaled <- data.plots$pct.exp.scaled * 100
  data.plots$Species <- factor(data.plots$Species,
                               levels = c("human", "macaque", "treeshrew", "mouse"))
  data.plots$features.plot <- factor(data.plots$features.plot, levels = hub_genes)

  return(data.plots)
}

# --- Dot plot for a specific module's hub genes (example: M15) ---
hub_df <- GetHubGenes(obj, n_hubs = 35)

module <- "M15"
hub_gene <- hub_df$gene_name[hub_df$module == module]

data.plots <- generate_dotplot_data(obj, hub_gene)

ggplot(data = data.plots, mapping = aes_string(x = "features.plot", y = "id")) +
  geom_point(mapping = aes_string(size = "pct.exp.scaled", color = "avg.exp.scaled")) +
  labs(
    x = "Features", y = "Identity",
    size = "Scale Percent Expressed",
    color = "Scale Mean Expression"
  ) +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(liidth = 1),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 0.5, vjust = 0.5)
  ) +
  facet_wrap(~Species, scales = "free_y", ncol = 1)

# =========================== Module-Specific Hub Gene Dot Plots ===========================

# Define representative hub genes for selected modules
# (Modules referenced in manuscript Figures 4C-D and Supplementary)
module_gene_lists <- list(
  "M1"  = c("BCAS3", "PDS5B", "PARG", "SMYD3", "PTGES3"),
  "M2"  = c("PPP2R2B", "HMGB4", "RETREG1", "HDAC11", "TSKS"),
  "M3"  = c("PABPC1", "H2AJ", "LRWD1", "SHCBP1L", "ZPBP2"),
  "M4"  = c("CNBD1", "CATSPER4", "CSNKA2IP", "RMDN2", "SPATA9"),
  "M5"  = c("PRKCQ", "MTUS2", "TMCC3", "RNF180", "FBXL17"),
  "M8"  = c("NPEPPS", "TTLL6", "SLC2A13", "EIF2B5", "LUC7L2"),
  "M9"  = c("FAM184A", "CDC20B", "CEP128", "ADAM2", "EXD1"),
  "M10" = c("JMJD1C", "DLG1", "NFAT5", "DPH6", "SPPL3"),
  "M15" = c("MLLT10", "SETX", "PTBP2", "NBEA", "WDR62")
)

# Additional modules
module_gene_lists_extra <- list(
  "M6"  = c("AIG1", "MVB12B", "SIKE1", "NAMPT", "P4HA2"),
  "M11" = c("ORMDL3", "NOL4", "AMN1", "AZIN2", "LRRTM3"),
  "M7"  = c("TEX38", "SYT16", "ADGRG1", "IGSF11", "ASGR1"),
  "M14" = c("LRBA", "ENOX1", "COL3A1", "PLD5", "PRKD1"),
  "M12" = c("DACH2", "SORCS1", "IQGAP1", "DOCK10", "PTPRD"),
  "M13" = c("STK31", "EML4", "BBS9", "TULP4", "DNAJC17")
)

# Combine all module gene lists
all_module_lists <- c(module_gene_lists, module_gene_lists_extra)

# Generate dot plots for each module
for (i in seq_along(all_module_lists)) {
  module_name <- names(all_module_lists)[i]
  hub_genes <- all_module_lists[[i]]

  cat(">>> Plotting hub gene dot plot for module:", module_name, "\n")

  data.plots <- generate_dotplot_data(obj, hub_genes)

  p <- ggplot(data = data.plots, mapping = aes_string(x = "features.plot", y = "id")) +
    geom_point(mapping = aes_string(size = "pct.exp.scaled", color = "avg.exp.scaled")) +
    labs(
      x = "Features", y = "Identity",
      size = "Scale Percent Expressed",
      color = "Scale Mean Expression",
      title = paste0("Module ", module_name, " Hub Genes")
    ) +
    scale_color_gradient(low = "lightgrey", high = "red") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(liidth = 1),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 0.5, vjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ) +
    facet_wrap(~Species, scales = "free_y", ncol = 1)

  ggsave(
    filename = paste0("./figure/",module_name, "-DotPlot.pdf"),
    plot = p, width = 4, height = 8
  )
}

cat(">>> 04_wgcna.R completed successfully!\n")