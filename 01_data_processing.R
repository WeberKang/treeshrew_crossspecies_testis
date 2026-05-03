################################################################################
# 01_data_processing.R
# Upstream processing of tree shrew testis single-nucleus multiome data &
# cross-species data preprocessing
#
# Functions:
#   1. Read and process tree shrew 10x Multiome (snRNA-seq + snATAC-seq) data
#   2. QC, normalization, dimensionality reduction, multiome integration (WNN)
#   3. Read published human, macaque, and mouse testis scRNA-seq data
#   4. Gene ID conversion (Ensembl ID -> Gene Symbol)
#   5. Unify cell type annotations: SG, SC, Early SD, Late SD, Sertoli, OSC
#
# Species abbreviations:
#   human    - Human (Homo sapiens)
#   macaque  - Rhesus macaque (Macaca mulatta)
#   mouse    - Mouse (Mus musculus)
#   treeshrew - Chinese tree shrew (Tupaia belangeri chinensis)
#
# Cell type nomenclature (aligned with manuscript):
#   SG        - Spermatogonia
#   SC        - Spermatocytes
#   Early SD  - Early stage Spermatids (round spermatids)
#   Late SD   - Late stage Spermatids (elongating spermatids)
#   Sertoli   - Sertoli cells
#   OSC       - Other Somatic Cells
################################################################################

# =========================== Environment Setup ===========================

# Set working directory
setwd("./cross_testis")

# ---------- Load R packages ----------
library(Seurat)          # Single-cell analysis framework
library(Signac)          # Single-cell chromatin accessibility analysis
library(ggplot2)         # Plotting
library(patchwork)       # Plot composition
library(clustree)        # Clustering tree visualization
library(reticulate)      # Python interface
library(openxlsx)        # Excel I/O
library(tidyverse)       # Data manipulation suite
library(GenomeInfoDb)    # Genome information utilities
library(GenomicFeatures) # Genomic feature handling
library(rtracklayer)     # GTF/GFF file I/O
library(data.table)      # Fast data reading (fread)
library(biomaRt)         # Ensembl gene annotation

# ---------- Configure Python environment for reticulate ----------
use_condaenv(
  condaenv = "Renv",
  conda = "/home/weber/Software/miniconda3/condabin/conda",
  required = TRUE
)
py_config()  # Verify Python configuration

# Check required Python modules
stopifnot(py_module_available('leidenalg'))
stopifnot(py_module_available('numpy'))
stopifnot(py_module_available('umap'))

# Clear workspace
rm(list = ls())

# =========================== Tree Shrew Genome Annotation ===========================

# Read tree shrew GTF annotation file
ts.Grange <- import("TS.fix.gtf")
ts.Grange <- keepStandardChromosomes(ts.Grange, pruning.mode = "coarse")

# Add transcript ID column for downstream annotation
mcols(ts.Grange)$tx_id <- mcols(ts.Grange)$transcript_id

# Unify chromosome naming style to UCSC
annotations <- ts.Grange
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "ts3"

# =========================== Read Tree Shrew 10x Multiome Data ===========================

# --- Sample 1 ---
s1 <- Read10X_h5(filename = "./data/sample1/filtered_feature_bc_matrix.h5")

# Separate RNA and ATAC counts
s1_rna_counts <- s1$`Gene Expression`
s1_atac_counts <- s1$Peaks

# Retain only tree shrew genes (TSDBGID prefix)
s1_rna_counts <- s1_rna_counts[
  contains(match = "TSDBGID", vars = s1_rna_counts@Dimnames[[1]]),
]

# --- Sample 2 ---
s2 <- Read10X_h5(filename = "./data/sample2/filtered_feature_bc_matrix.h5")

s2_rna_counts <- s2$`Gene Expression`
s2_atac_counts <- s2$Peaks

s2_rna_counts <- s2_rna_counts[
  contains(match = "TSDBGID", vars = s2_rna_counts@Dimnames[[1]]),
]

# --- Convert gene IDs (TSDBGID) to gene names ---
id2name <- read.csv(file = "gene.csv", header = TRUE, row.names = 2)

s1_rna_counts@Dimnames[[1]] <- id2name[s1_rna_counts@Dimnames[[1]], ]$gene_name_re
s2_rna_counts@Dimnames[[1]] <- id2name[s2_rna_counts@Dimnames[[1]], ]$gene_name_re

# Remove duplicated gene names
s1_rna_counts <- s1_rna_counts[unique(rownames(s1_rna_counts)), ]
s2_rna_counts <- s2_rna_counts[unique(rownames(s2_rna_counts)), ]

# =========================== Create Seurat Objects (RNA) ===========================

# Create RNA Seurat objects (retain genes expressed in at least 10 cells)
s1.RNA <- CreateSeuratObject(
  counts = s1_rna_counts,
  project = "s1.RNA",
  min.cells = 10,
  assay = "RNA"
)
s2.RNA <- CreateSeuratObject(
  counts = s2_rna_counts,
  project = "s2.RNA",
  min.cells = 10,
  assay = "RNA"
)

# =========================== Create ChromatinAssay Objects (ATAC) ===========================

# --- Filter peaks on non-standard chromosomes ---

# Sample 1
s1.grange.counts <- StringToGRanges(rownames(s1_atac_counts), sep = c(":", "-"))
s1.grange.use <- seqnames(s1.grange.counts) %in% standardChromosomes(s1.grange.counts)
s1_atac_counts <- s1_atac_counts[as.vector(s1.grange.use), ]

# Sample 2
s2.grange.counts <- StringToGRanges(rownames(s2_atac_counts), sep = c(":", "-"))
s2.grange.use <- seqnames(s2.grange.counts) %in% standardChromosomes(s2.grange.counts)
s2_atac_counts <- s2_atac_counts[as.vector(s2.grange.use), ]

# --- Create ATAC assays ---
s1.ATAC <- CreateChromatinAssay(
  counts = s1_atac_counts,
  sep = c(":", "-"),
  fragments = './data/sample1/atac_fragments.tsv.gz',
  min.cells = 10,
  annotation = annotations
)

s2.ATAC <- CreateChromatinAssay(
  counts = s2_atac_counts,
  sep = c(":", "-"),
  fragments = './data/sample2/atac_fragments.tsv.gz',
  min.cells = 10,
  annotation = annotations
)

# =========================== Merge RNA and ATAC Assays ===========================

# Add ATAC assay to RNA Seurat objects
s1 <- s1.RNA
s1[['ATAC']] <- s1.ATAC

s2 <- s2.RNA
s2[['ATAC']] <- s2.ATAC

# =========================== Quality Control (QC) ===========================

# Calculate mitochondrial gene percentage
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^MT-")
s2[["percent.mt"]] <- PercentageFeatureSet(s2, pattern = "^MT-")

# Visualize QC metric distributions
VlnPlot(s1,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
               "nCount_ATAC", "nFeature_ATAC"),
  ncol = 5, log = TRUE, pt.size = 0
) + NoLegend()

VlnPlot(s2,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
               "nCount_ATAC", "nFeature_ATAC"),
  ncol = 5, log = TRUE, pt.size = 0
) + NoLegend()

# --- Apply QC filters ---
s1 <- subset(
  x = s1,
  nFeature_RNA > 200 & nFeature_RNA < 6000 &
    percent.mt < 2 &
    nCount_RNA < 8000 &
    nCount_ATAC > 100 & nCount_ATAC < 10000
)

s2 <- subset(
  x = s2,
  nFeature_RNA > 200 & nFeature_RNA < 6000 &
    percent.mt < 2 &
    nCount_RNA < 8000 &
    nCount_ATAC > 100 & nCount_ATAC < 10000
)

# =========================== Merge Samples & RNA Analysis ===========================

# Merge the two biological replicates
merge.s <- merge(
  s1, s2,
  add.cell.ids = c("s1", "s2"),
  project = "Ts"
)

# --- RNA normalization and dimensionality reduction ---
DefaultAssay(merge.s) <- "RNA"
merge.s <- merge.s %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:7, min.dist = 0.1, reduction.key = 'rnaUMAP_') %>%
  RunTSNE(dims = 1:14, reduction.key = 'rnaTSNE_')

# Visualize RNA UMAP
DimPlot(merge.s,
  group.by = 'CELL_TYPE',
  label = TRUE, repel = FALSE,
  reduction = 'umap',
  pt.size = 1, label.size = 5
)

# =========================== ATAC Analysis ===========================

DefaultAssay(merge.s) <- "ATAC"

# Calculate nucleosome signal and TSS enrichment scores
merge.s <- NucleosomeSignal(object = merge.s)
merge.s <- TSSEnrichment(object = merge.s, fast = TRUE)

# Scatter plot: TSS enrichment vs ATAC counts
pdf(file = 'DensityScatter.pdf', width = 8, height = 8)
DensityScatter(merge.s,
  x = 'nCount_ATAC', y = 'TSS.enrichment',
  log_x = TRUE, quantiles = TRUE
)
dev.off()

# ATAC normalization and dimensionality reduction
merge.s <- RunTFIDF(merge.s)
merge.s <- FindTopFeatures(merge.s, min.cutoff = 'q5')
merge.s <- RunSVD(merge.s)

# Check correlation between SVD dimensions and sequencing depth
DepthCor(merge.s, n = 50)

# =========================== Save Merged Tree Shrew Object ===========================

save(merge.s, file = "merge.s.Rdata")

# =========================== WNN Multi-modal Integration ===========================

# NOTE: CELL_TYPE column will be populated after cell annotation step
merge.s <- FindMultiModalNeighbors(
  merge.s,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:7, 3:20)
)

merge.s <- RunUMAP(
  merge.s,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

# Visualize WNN UMAP
DimPlot(merge.s,
  reduction = "wnn.umap",
  group.by = "CELL_TYPE",
  label = TRUE, pt.size = 1, label.size = 5
)

# =========================== Read Published Cross-Species Data ===========================

# Read cell-level metadata (feature table with annotations and UMAP coordinates)
feature <- fread("./data/feature.csv")
rownames(feature) <- feature$CELL

# ----- Mouse data -----
cat(">>> Loading mouse data...\n")
counts <- fread("./data/Mouse.merged.raw.counts.txt")
counts <- as.matrix(counts, rownames = 1)

# Verify cell barcode consistency between counts and metadata
stopifnot(sum(rownames(feature) %in% colnames(counts)) == dim(counts)[2])

# Subset metadata to matching cells
meta.data <- as.data.frame(feature[rownames(feature) %in% colnames(counts), ])
rownames(meta.data) <- meta.data$CELL
meta.data <- meta.data[colnames(counts), ]

# Create Seurat object and embed published UMAP coordinates
mouse <- CreateSeuratObject(counts, project = "Mouse", meta.data = meta.data)
umap <- as.matrix(data.frame(
  umap_1 = meta.data[, "UMAP1"],
  umap_2 = meta.data[, "UMAP2"],
  row.names = rownames(meta.data)
))
mouse@reductions[["umap"]] <- CreateDimReducObject(
  embeddings = umap,
  key = "umap_",
  assay = DefaultAssay(mouse)
)

# Remove unannotated cells
mouse <- mouse[, !is.na(mouse@meta.data$CELL_TYPE)]
DimPlot(mouse, group.by = "CELL_TYPE", label = TRUE)
save(mouse, file = "./data/mouse.RData")

# ----- Macaque data -----
cat(">>> Loading macaque data...\n")
counts <- fread("./data/Macaque.merged.raw.counts.txt")
counts <- as.matrix(counts, rownames = 1)

stopifnot(sum(rownames(feature) %in% colnames(counts)) == dim(counts)[2])

meta.data <- as.data.frame(feature[rownames(feature) %in% colnames(counts), ])
rownames(meta.data) <- meta.data$CELL
meta.data <- meta.data[colnames(counts), ]

macaque <- CreateSeuratObject(counts, project = "Macaque", meta.data = meta.data)
umap <- as.matrix(data.frame(
  umap_1 = meta.data[, "UMAP1"],
  umap_2 = meta.data[, "UMAP2"],
  row.names = rownames(meta.data)
))
macaque@reductions[["umap"]] <- CreateDimReducObject(
  embeddings = umap,
  key = "umap_",
  assay = DefaultAssay(macaque)
)

macaque <- macaque[, !is.na(macaque@meta.data$CELL_TYPE)]
DimPlot(macaque, group.by = "CELL_TYPE", label = TRUE)
save(macaque, file = "./data/macaque.RData")

# ----- Human data -----
cat(">>> Loading human data...\n")
counts <- fread("./data/human.merged.raw.counts.txt")
counts <- as.matrix(counts, rownames = 1)

stopifnot(sum(rownames(feature) %in% colnames(counts)) == dim(counts)[2])

meta.data <- as.data.frame(feature[rownames(feature) %in% colnames(counts), ])
rownames(meta.data) <- meta.data$CELL
meta.data <- meta.data[colnames(counts), ]

human <- CreateSeuratObject(counts, project = "Human", meta.data = meta.data)
umap <- as.matrix(data.frame(
  umap_1 = meta.data[, "UMAP1"],
  umap_2 = meta.data[, "UMAP2"],
  row.names = rownames(meta.data)
))
human@reductions[["umap"]] <- CreateDimReducObject(
  embeddings = umap,
  key = "umap_",
  assay = DefaultAssay(human)
)

human <- human[, !is.na(human@meta.data$CELL_TYPE)]
DimPlot(human, group.by = "CELL_TYPE", label = TRUE)
save(human, file = "./data/human.RData")

# =========================== Gene ID Conversion: Ensembl -> Gene Symbol ===========================

load("./data/human.RData")
load("./data/macaque.RData")
load("./data/mouse.RData")

# Connect to Ensembl BioMart
ensembl <- useEnsembl(biomart = "genes", mirror = "asia")

# --- Human ---
ensembl_hsa <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene_names <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name"),
  values = Features(human),
  mart = ensembl_hsa
)
gene_names <- gene_names[nchar(gene_names$external_gene_name) > 0, ]
human <- human[Features(human) %in% gene_names$ensembl_gene_id, ]
gene_names <- gene_names[match(Features(human), gene_names$ensembl_gene_id), ]
gene_names$external_gene_name <- make.unique(gene_names$external_gene_name)
rownames(human@assays$RNA@features) <- gene_names$external_gene_name
save(human, file = "./data/human.id2name.RData")

# --- Mouse ---
ensembl_mmus <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
gene_names <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name"),
  values = Features(mouse),
  mart = ensembl_mmus
)
gene_names <- gene_names[nchar(gene_names$external_gene_name) > 0, ]
mouse <- mouse[Features(mouse) %in% gene_names$ensembl_gene_id, ]
gene_names <- gene_names[match(Features(mouse), gene_names$ensembl_gene_id), ]
gene_names$external_gene_name <- make.unique(gene_names$external_gene_name)
rownames(mouse@assays$RNA@features) <- gene_names$external_gene_name
save(mouse, file = "./data/mouse.id2name.RData")

# --- Macaque ---
ensembl_mmul <- useDataset(dataset = "mmulatta_gene_ensembl", mart = ensembl)
gene_names <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name"),
  values = Features(macaque),
  mart = ensembl_mmul
)
gene_names <- gene_names[nchar(gene_names$external_gene_name) > 0, ]
macaque <- macaque[Features(macaque) %in% gene_names$ensembl_gene_id, ]
gene_names <- gene_names[match(Features(macaque), gene_names$ensembl_gene_id), ]
gene_names$external_gene_name <- make.unique(gene_names$external_gene_name)
rownames(macaque@assays$RNA@features) <- gene_names$external_gene_name
save(macaque, file = "./data/macaque.id2name.RData")

# =========================== Tree Shrew Cell Type Annotation (Human Reference-Based) ===========================

load("./data/treeshrew.RData")
load("./data/human.id2name.RData")

# QC filtering for tree shrew
treeshrew <- subset(
  x = treeshrew,
  nFeature_RNA > 200 & nFeature_RNA < 6000 &
    percent.mt < 0.5 &
    nCount_RNA < 8000
)

# Normalization and dimensionality reduction
treeshrew <- treeshrew %>%
  JoinLayers() %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:25)

# Normalize human reference
human <- human %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:25)

# Transfer labels from human reference to tree shrew query
anchors <- FindTransferAnchors(
  reference = human,
  query = treeshrew,
  dims = 1:25,
  reference.reduction = 'pca'
)
pred <- TransferData(
  anchorset = anchors,
  refdata = human$CELL_TYPE,
  dims = 1:25
)
treeshrew <- AddMetaData(treeshrew, metadata = pred)

# Use predicted cell type as initial annotation
treeshrew@meta.data["CELL_TYPE"] <- treeshrew@meta.data["predicted.id"]
DimPlot(treeshrew, group.by = 'predicted.id', label = TRUE, repel = FALSE, reduction = 'umap')
save(treeshrew, file = "./data/treeshrew.id2name.RData")

# =========================== Unify Cell Type Nomenclature ===========================
# Map all species to unified cell types: SG, SC, Early SD, Late SD, Sertoli, OSC
# Original mapping: rSD -> Early SD; eSD -> Late SD; Other somatic -> OSC

load("./data/human.id2name.RData")
load("./data/treeshrew.id2name.RData")
load("./data/mouse.id2name.RData")
load("./data/macaque.id2name.RData")

# --- Human ---
Idents(human) <- human@meta.data$CELL_TYPE
NEW_CELL_TYPE <- c(
  "SG", "SG",               # First two clusters -> SG
  "OSC",                    # -> Other Somatic Cells
  "Early SD", "Early SD",   # rSD -> Early SD
  "SC", "SC", "SC",         # -> Spermatocytes
  "Sertoli",                # -> Sertoli cells
  "Late SD"                 # eSD -> Late SD
)
names(NEW_CELL_TYPE) <- levels(human)
human <- RenameIdents(human, NEW_CELL_TYPE)
human@meta.data$NEW_CELL_TYPE <- Idents(human)
human@meta.data$NEW_CELL_TYPE <- factor(
  human@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
)
save(human, file = "./data/human.newcelltype.RData")

# --- Tree shrew ---
Idents(treeshrew) <- treeshrew@meta.data$CELL_TYPE
NEW_CELL_TYPE <- c(
  "Late SD", "SC", "Early SD", "SG",
  "Early SD", "SC", "OSC", "SC",
  "SG", "Sertoli"
)
names(NEW_CELL_TYPE) <- levels(treeshrew)
treeshrew <- RenameIdents(treeshrew, NEW_CELL_TYPE)
treeshrew@meta.data$NEW_CELL_TYPE <- Idents(treeshrew)
treeshrew@meta.data$NEW_CELL_TYPE <- factor(
  treeshrew@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
)
save(treeshrew, file = "./data/treeshrew.newcelltype.RData")

# --- Macaque ---
Idents(macaque) <- macaque@meta.data$CELL_TYPE
NEW_CELL_TYPE <- c(
  "Early SD", "Early SD", "Late SD", "SC",
  "SG", "SC", "SC", "OSC",
  "SG", "Sertoli"
)
names(NEW_CELL_TYPE) <- levels(macaque)
macaque <- RenameIdents(macaque, NEW_CELL_TYPE)
macaque@meta.data$NEW_CELL_TYPE <- Idents(macaque)
macaque@meta.data$NEW_CELL_TYPE <- factor(
  macaque@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
)
save(macaque, file = "./data/macaque.newcelltype.RData")

# --- Mouse ---
Idents(mouse) <- mouse@meta.data$CELL_TYPE
NEW_CELL_TYPE <- c(
  "SC", "Early SD", "Sertoli", "Late SD", "SG", "OSC"
)
names(NEW_CELL_TYPE) <- levels(mouse)
mouse <- RenameIdents(mouse, NEW_CELL_TYPE)
mouse@meta.data$NEW_CELL_TYPE <- Idents(mouse)
mouse@meta.data$NEW_CELL_TYPE <- factor(
  mouse@meta.data$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
)
save(mouse, file = "./data/mouse.newcelltype.RData")

# =========================== Per-Species Clustering & Fine Annotation ===========================

# --- Human ---
species <- "human"
obj <- get(species)
obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:10, min.dist = 0.1) %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.5)
DimPlot(obj, group.by = c("NEW_CELL_TYPE", "seurat_clusters"),
        label = TRUE, pt.size = 1, label.size = 6)
assign(x = species, value = obj)

# --- Macaque ---
species <- "macaque"
obj <- get(species)
obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:11, min.dist = 0.3) %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.5)
DimPlot(obj, group.by = c("NEW_CELL_TYPE", "seurat_clusters"),
        label = TRUE, pt.size = 1, label.size = 6)
assign(x = species, value = obj)

# --- Mouse ---
species <- "mouse"
obj <- get(species)
obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:26, min.dist = 0.3) %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.5)
DimPlot(obj, group.by = c("NEW_CELL_TYPE", "seurat_clusters"),
        label = TRUE, pt.size = 1, label.size = 6)
assign(x = species, value = obj)

# --- Tree shrew (with additional cleaning steps) ---
species <- "treeshrew"
obj <- get(species)
obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:9, min.dist = 0.1) %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.5)

# Remove clusters with intronic contamination and low-quality cells
obj <- subset(obj, idents = c(6, 14), invert = TRUE)
obj <- subset(obj, prediction.score.max > 0.5)

# Re-run normalization and dimensionality reduction after cleaning
obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:7, min.dist = 0.1) %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.5)

# Fine-grained cell type annotation
new_cell_type <- c(
  "Late SD", "Late SD", "Early SD", "SC",
  "Early SD", "SC", "Early SD", "SC",
  "SG", "SC", "SC", "Early SD",
  "Early SD", "OSC", "Sertoli"
)
names(new_cell_type) <- levels(obj)
obj <- RenameIdents(obj, new_cell_type)
obj$NEW_CELL_TYPE <- Idents(obj)
obj$NEW_CELL_TYPE <- factor(
  obj$NEW_CELL_TYPE,
  levels = c("SG", "SC", "Early SD", "Late SD", "Sertoli", "OSC")
)
DimPlot(obj, group.by = c("NEW_CELL_TYPE", "seurat_clusters"),
        label = TRUE, pt.size = 1, label.size = 6)

# =========================== X Chromosome Gene Fraction Analysis (Sex Verification) ===========================

load("./data2/species.RData")
ensembl <- useEnsembl(biomart = "genes", mirror = "asia")

# --- Human X chromosome genes ---
ensembl_hsa <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
geneAttr <- getBM(
  filters = "external_gene_name",
  attributes = c("external_gene_name", "chromosome_name"),
  values = Features(human),
  mart = ensembl_hsa
)
Xgene <- geneAttr$external_gene_name[geneAttr$chromosome_name == "X"]
human[["chrX"]] <- PercentageFeatureSet(human, features = Xgene)

# --- Macaque X chromosome genes ---
ensembl_mmul <- useDataset(dataset = "mmulatta_gene_ensembl", mart = ensembl)
geneAttr <- getBM(
  filters = "external_gene_name",
  attributes = c("external_gene_name", "chromosome_name"),
  values = Features(macaque),
  mart = ensembl_mmul
)
Xgene <- geneAttr$external_gene_name[geneAttr$chromosome_name == "X"]
macaque[["chrX"]] <- PercentageFeatureSet(macaque, features = Xgene)

# --- Mouse X chromosome genes ---
ensembl_mmus <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
geneAttr <- getBM(
  filters = "external_gene_name",
  attributes = c("external_gene_name", "chromosome_name"),
  values = Features(mouse),
  mart = ensembl_mmus
)
Xgene <- geneAttr$external_gene_name[geneAttr$chromosome_name == "X"]
mouse[["chrX"]] <- PercentageFeatureSet(mouse, features = Xgene)

# --- Tree shrew X chromosome genes (extracted from GTF file) ---
gff <- readGFF("./data/TS.gtf")
mapid <- gff[gff$seqid == "chrX" & gff$type == "transcript",
             c("gene_id", "gene_name")]
mapid <- mapid[!(duplicated(mapid) | is.na(mapid$gene_name)), ]
Xgene <- unique(gsub(x = mapid$gene_name, pattern = "li[0-9]*", replacement = ""))
Xgene <- Xgene[Xgene %in% Features(treeshrew)]
treeshrew[["chrX"]] <- PercentageFeatureSet(treeshrew, features = Xgene)

# --- Generate violin plots of X chromosome gene fraction per species ---
for (species in c("human", "macaque", "mouse", "treeshrew")) {
  obj <- get(species)
  df <- obj@meta.data[, c("nFeature_RNA", "nCount_RNA", "chrX", "NEW_CELL_TYPE")]

  # Panel 1: Number of genes per cell type
  p1 <- ggplot(data = df, aes(x = NEW_CELL_TYPE, y = nFeature_RNA, fill = NEW_CELL_TYPE)) +
    geom_violin(trim = FALSE, color = "white") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(),
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    labs(x = "", y = "nGene", title = species) +
    NoLegend()

  # Panel 2: Number of UMIs per cell type
  p2 <- ggplot(data = df, aes(x = NEW_CELL_TYPE, y = nCount_RNA, fill = NEW_CELL_TYPE)) +
    geom_violin(trim = FALSE, color = "white") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line()
    ) +
    labs(x = "", y = "nUMI")

  # Panel 3: X chromosome gene percentage per cell type
  p3 <- ggplot(data = df, aes(x = NEW_CELL_TYPE, y = chrX, fill = NEW_CELL_TYPE)) +
    geom_violin(trim = FALSE, color = "white") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line()
    ) +
    labs(x = "", y = "%chrX") +
    NoLegend()

  # Combine panels and save
  p <- p1 / p2 / p3
  ggsave(
    filename = paste0("./figure/", species, ".pdf"),
    plot = p, width = 9, height = 12
  )
}

# =========================== Save Final Objects ===========================

save(human, macaque, mouse, treeshrew, file = "./data/all_species_processed.RData")

cat(">>> 01_data_processing.R completed successfully!\n")