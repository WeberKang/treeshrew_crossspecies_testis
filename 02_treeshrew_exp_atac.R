################################################################################
# 02_treeshrew_exp_atac.R
# Downstream analysis of tree shrew single-nucleus multiome data
# (snRNA-seq + snATAC-seq)
#
# Functions:
#   1. Gene activity score calculation from ATAC data
#   2. Differential gene activity & differential peak analysis
#   3. Heatmap visualization of marker gene expression vs activity
#   4. ATAC peak genomic feature annotation (ChIPseeker)
#   5. Peak-to-TSS enrichment profiling across cell types
#   6. Transcriptional delay analysis (ATAC peak vs RNA peak timing)
#
# Species: Chinese tree shrew (Tupaia belangeri chinensis)
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
library(Signac)            # Single-cell chromatin accessibility analysis
library(ggplot2)           # Plotting
library(patchwork)         # Plot composition
library(clustree)          # Clustering tree visualization
library(reticulate)        # Python interface
library(openxlsx)          # Excel I/O
library(tidyverse)         # Data manipulation suite
library(GenomeInfoDb)      # Genome information utilities
library(GenomicFeatures)   # Genomic feature handling
library(rtracklayer)       # GTF/GFF file I/O
library(JASPAR2022)        # JASPAR motif database
library(TFBSTools)         # Transcription factor binding site tools
library(BSgenome.Treeshrew.yaolab.ts3)  # Tree shrew reference genome
library(pheatmap)           # Heatmap plotting
library(scales)             # Scale functions for ggplot2
library(ChIPseeker)         # ChIP peak annotation

# ---------- Configure Python environment ----------
use_condaenv(
  condaenv = "Renv",
  conda = "/home/weber/Software/miniconda3/condabin/conda",
  required = TRUE
)
py_config()
stopifnot(py_module_available('leidenalg'))
stopifnot(py_module_available('numpy'))
stopifnot(py_module_available('umap'))

# =========================== Gene Activity Score Calculation ===========================

load(file = "merge.s.Rdata")

# Calculate gene activity scores from ATAC counts (summarize peaks near each gene)
DefaultAssay(merge.s) <- "ATAC"
gene.activates <- GeneActivity(merge.s)

# Ensure unique row names
rownames(gene.activates) <- make.unique(rownames(gene.activates))

# Store gene activity as a new assay
merge.s[["activates"]] <- CreateAssayObject(counts = gene.activates)

# Normalize gene activity scores
merge.s <- NormalizeData(
  merge.s,
  assay = "activates",
  normalization.method = 'LogNormalize',
  scale.factor = 10000
)

save(merge.s, file = "merge.s.Rdata")

# =========================== Differential Gene Activity Analysis ===========================

load(file = "merge.s.Rdata")

DefaultAssay(merge.s) <- "activates"
Idents(merge.s) <- merge.s@meta.data$NEW_CELL_TYPE

# Find differentially active genes per cell type
de_peaks <- FindAllMarkers(
  merge.s,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
write.csv(de_peaks, file = "de_activates_gene_peak.csv")

# =========================== Helper Function: Save pheatmap to PDF ===========================

save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# =========================== Marker Gene Heatmap (Activity vs RNA) ===========================

# Define canonical marker genes for each germ cell type
plot_genes <- c(
  "FOXO1", "MKI67", "GCNA", "KIT", "NANOS3",          # SG (Spermatogonia)
  "SPO11", "SYCP1", "MSH4", "SYCE1", "OVOL2",          # SC (Spermatocytes)
  "FER1L5", "DNAH3", "C2CD6", "SPATA17", "SUN5",      # Early SD (round spermatids)
  "PRM3", "TNP2", "TNP1"                                # Late SD (elongating spermatids)
)

# --- Gene Activity heatmap ---
AverageExpression(merge.s, features = plot_genes, assay = "activates")$activates %>%
  as.data.frame() %>%
  dplyr::select(-c("OSC", "Sertoli")) %>%
  subset(select = c("SG", "SC", "Early SD", "Late SD")) -> avg_activates

avg_activates <- avg_activates[plot_genes, ]

p1 <- pheatmap(
  avg_activates,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Gene Activity",
  silent = TRUE
)
save_pheatmap_pdf(p1, "./figure/Activates.pdf")

# --- RNA expression heatmap ---
AverageExpression(merge.s, features = plot_genes, assay = "RNA")$RNA %>%
  as.data.frame() %>%
  dplyr::select(-c("OSC", "Sertoli")) %>%
  subset(select = c("SG", "SC", "Early SD", "Late SD")) -> avg_rna

avg_rna <- avg_rna[plot_genes[plot_genes %in% rownames(avg_rna)], ]

p2 <- pheatmap(
  avg_rna,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "RNA Expression",
  silent = TRUE
)
save_pheatmap_pdf(p2, "./figure/RNA.pdf")

# =========================== ATAC Count Distribution per Cell Type ===========================

df <- FetchData(merge.s, vars = c("nCount_ATAC", "NEW_CELL_TYPE"))

# Filter to germ cell types only
df <- df[df$NEW_CELL_TYPE %in% c("SG", "SC", "Early SD", "Late SD"), ]
df$NEW_CELL_TYPE <- factor(df$NEW_CELL_TYPE, levels = c("SG", "SC", "Early SD", "Late SD"))

ggplot(df, aes(x = NEW_CELL_TYPE, y = nCount_ATAC, fill = NEW_CELL_TYPE)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = hue_pal()(6)[1:4]) +
  theme_bw() +
  NoLegend()
ggsave("./figure/nCounts.ATAC.pdf", width = 5, height = 5)

# =========================== Differential Peak Analysis & Genomic Annotation ===========================

# Build TxDb from tree shrew GTF for peak annotation
TxDb.ts <- makeTxDbFromGFF(file = "TS.fix.gtf")
TxDb.ts <- keepStandardChromosomes(TxDb.ts, pruning.mode = "coarse")

# Define promoter region (±3000 bp from TSS)
promoter <- getPromoters(TxDb = TxDb.ts, upstream = 3000, downstream = 3000)

# --- Find differentially accessible peaks per cell type ---
DefaultAssay(merge.s) <- "ATAC"
Idents(merge.s) <- factor(
  merge.s@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
)

de_peaks <- FindAllMarkers(
  merge.s,
  only.pos = TRUE,
  logfc.threshold = 1,
  test.use = "wilcox"
)

# Filter significant peaks and sort
de_peaks <- de_peaks %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, -avg_log2FC)

write.csv(de_peaks, file = "de_peaks.csv")

# --- Per-cell-type peak annotation and TSS enrichment profiling ---
tagMatrixs <- list()
peakAnnos <- list()

for (celltype in levels(Idents(merge.s))) {
  # Extract peak coordinates for this cell type
  peak <- de_peaks[de_peaks$cluster == celltype, ]$gene
  peak <- do.call(rbind, strsplit(peak, "-"))
  peak <- GRanges(
    seqnames = peak[, 1],
    ranges = IRanges(start = as.numeric(peak[, 2]), end = as.numeric(peak[, 3])),
    strand = "*"
  )

  # Get tag matrix around promoter for TSS enrichment profile
  tagMatrix <- getTagMatrix(peak, windows = promoter)
  tagMatrixs[[celltype]] <- tagMatrix

  # Annotate peaks with genomic features
  peakAnno <- annotatePeak(
    peak = peak,
    TxDb = TxDb.ts,
    tssRegion = c(-3000, 3000),
    verbose = TRUE
  )
  peakAnnos[[celltype]] <- peakAnno
}

save(tagMatrixs, peakAnnos, file = "chipseeker.Rdata")

# --- Visualize peak annotations ---
load(file = "chipseeker.Rdata")

# Bar plot: genomic feature distribution of peaks per cell type
p <- plotAnnoBar(peakAnnos)
ggsave(plot = p, filename = "./figure/plotAnnoBar.pdf", width = 10, height = 10)

# Profile plot: TSS enrichment signal per cell type
p <- plotPeakProf(tagMatrixs, resample = 500, facet = "row", free_y = FALSE)
ggsave(plot = p, filename = "./figure/plotPeakProf.pdf", width = 10, height = 10)

# =========================== Transcriptional Delay Analysis ===========================
# Compare the timing of ATAC peak (chromatin opening) vs RNA peak (expression)
# for each gene to identify genes with delayed transcription

load(file = "merge.s.Rdata")

# Set cell type identities
Idents(merge.s) <- factor(
  merge.s$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
)

# --- Find marker genes/peaks for each assay ---
DefaultAssay(merge.s) <- "RNA"
df.RNA <- FindAllMarkers(
  merge.s,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

DefaultAssay(merge.s) <- "activates"
df.activates <- FindAllMarkers(
  merge.s,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Subset to germ cells only
merge.s <- subset(merge.s, subset = NEW_CELL_TYPE %in% c("SG", "SC", "Early SD", "Late SD"))

# --- Calculate average expression and activity per cell type ---
DefaultAssay(merge.s) <- "RNA"
avg.rna <- AverageExpression(merge.s, assays = "RNA")$RNA

DefaultAssay(merge.s) <- "activates"
avg.atac <- AverageExpression(merge.s, assays = "activates")$activates

# --- Take intersection of genes present in both assays ---
common.genes <- intersect(rownames(avg.rna), rownames(avg.atac))
avg.rna <- avg.rna[common.genes, ]
avg.atac <- avg.atac[common.genes, ]

# Define germ cell developmental order
cell_order <- c("SG", "SC", "Early SD", "Late SD")

# ----- Helper function: Identify peak stage, value, and fold change -----
# For each gene, find the cell type with the highest signal,
# the peak value, and the fold change relative to the second-highest
get_peak_info <- function(mat, cell_order) {
  res <- apply(mat[, cell_order], 1, function(x) {
    vals <- as.numeric(x)
    names(vals) <- cell_order
    o <- order(vals, decreasing = TRUE)
    peak_val <- vals[o[1]]
    peak_stage <- cell_order[o[1]]
    second_val <- vals[o[2]]
    fold_change <- ifelse(second_val == 0, NA, peak_val / second_val)
    data.frame(stage = peak_stage, value = peak_val, fold_change = fold_change)
  })
  res <- do.call(rbind, res)
  rownames(res) <- rownames(mat)
  return(res)
}

# --- Get peak timing for both ATAC (chromatin) and RNA (expression) ---
atac_peak <- get_peak_info(avg.atac, cell_order)
rna_peak  <- get_peak_info(avg.rna, cell_order)

# --- Quantify transcriptional delay ---
# Convert stage to numeric for delay calculation
stage_map <- c(SG = 1, SC = 2, `Early SD` = 3, `Late SD` = 4)

delay_df <- data.frame(
  gene             = common.genes,
  atac_stage       = atac_peak$stage,
  rna_stage        = rna_peak$stage,
  atac_peak        = atac_peak$value,
  rna_peak         = rna_peak$value,
  atac_fold_change = atac_peak$fold_change,
  rna_fold_change  = rna_peak$fold_change,
  atac_num         = stage_map[atac_peak$stage],
  rna_num          = stage_map[rna_peak$stage]
) %>%
  mutate(delay = rna_num - atac_num) %>%
  arrange(delay)

# --- Filter high-confidence delayed genes ---
# Criteria: sufficient signal in both assays, clear fold change, and RNA lags behind ATAC
high_confidence_genes <- delay_df %>%
  filter(
    atac_peak > 0.4,
    rna_peak > 0.2,
    atac_fold_change > 1.25,
    rna_fold_change > 1.25,
    delay > 0
  ) %>%
  arrange(desc(delay), desc(atac_fold_change), desc(rna_fold_change))

# --- Plot individual gene delay profiles ---
for (gene in high_confidence_genes$gene) {
  cat("Plotting delay profile:", gene, "\n")

  df <- data.frame(
    stage = factor(cell_order, levels = cell_order),
    RNA = as.numeric(avg.rna[gene, cell_order]),
    ATAC = as.numeric(avg.atac[gene, cell_order])
  )

  # Scale factor for dual-axis plotting
  max_rna <- max(df$RNA, na.rm = TRUE)
  max_atac <- max(df$ATAC, na.rm = TRUE)
  scale_factor <- max_atac / max_rna

  ggplot(df, aes(x = stage)) +
    # ATAC signal (left y-axis)
    geom_line(aes(y = ATAC, group = 1, color = "ATAC"), linewidth = 1.5) +
    geom_point(aes(y = ATAC, color = "ATAC"), size = 3) +
    # RNA expression (right y-axis, scaled)
    geom_line(aes(y = RNA * scale_factor, group = 1, color = "RNA"), linewidth = 1.5) +
    geom_point(aes(y = RNA * scale_factor, color = "RNA"), size = 3) +
    # Dual y-axis
    scale_y_continuous(
      name = "ATAC Signal",
      sec.axis = sec_axis(~ . / scale_factor, name = "RNA Expression")
    ) +
    scale_color_manual(values = c("ATAC" = "#E64B35", "RNA" = "#4DBBD5")) +
    ggtitle(paste0(gene, " - Transcriptional Delay")) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    )

  ggsave(
    filename = paste0("./figure/Delay/Delay_Plot_", gene, ".pdf"),
    width = 6, height = 4
  )
}

cat(">>> 02_treeshrew_exp_atac.R completed successfully!\n")