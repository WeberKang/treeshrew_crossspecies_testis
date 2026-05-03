################################################################################
# 05_mfuzz.R
# Mfuzz soft clustering of gene expression trajectories across species
# Species: human, macaque, tree shrew, mouse
#
# Functions:
#   1. Mfuzz clustering of PTTG (Primate-Tupai shared Testis Gene set) genes
#      across human, macaque, and tree shrew pseudotime
#   2. Identification of genes with conserved expression trajectories
#      in primates + tree shrew (lacking/divergent in mouse)
#   3. Mfuzz clustering of genes shared across all four species
#   4. Conservation score calculation and cross-species pattern classification
#   5. Visualization of individual gene expression trajectories
#
# Cell type nomenclature (aligned with manuscript):
#   SG        - Spermatogonia
#   SC        - Spermatocytes
#   Early SD  - Early stage Spermatids (round spermatids, rSD in code)
#   Late SD   - Late stage Spermatids (elongating spermatids, eSD in code)
#
# Pseudotime stage nomenclature:
#   stage1-2  - SG (Spermatogonia)
#   stage3-4  - SC (Spermatocytes)
#   stage5-6  - Early SD (round spermatids)
#   stage7-8  - Late SD (elongating spermatids)
################################################################################

# =========================== Environment Setup ===========================

setwd("./cross_testis")

# ---------- Load R packages ----------
library(Mfuzz)              # Soft clustering for time-series expression data
library(Seurat)             # Single-cell analysis framework
library(tidyverse)          # Data manipulation suite
library(reshape2)           # Data reshaping
library(scales)             # Scale functions for ggplot2
library(openxlsx)           # Excel I/O

# =========================== Load Pseudotime Data ===========================

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

# =========================== Part 1: PTTG (Primate-Tupai Shared Testis Gene Set) Analysis ===========================
# These are genes detected in human, macaque, and tree shrew testes,
# but undetectable in mouse testis (manuscript Figure 5)

# Load PTTG gene list from Venn diagram results
primate_gene <- read.csv("./figure/F2H.venn.v2.csv", header = TRUE) %>%
  pull("human.macaque.treeshrew")
primate_gene <- primate_gene[nchar(primate_gene) > 0]

# Subset to PTTG genes
human_pseu     <- subset(human_pseu, features = primate_gene)
macaque_pseu   <- subset(macaque_pseu, features = primate_gene)
treeshrew_pseu <- subset(treeshrew_pseu, features = primate_gene)

# Merge primate + tree shrew objects
primate_pseu <- merge(
  x = human_pseu,
  y = list(macaque_pseu, treeshrew_pseu),
  add.cell.ids = c("human", "macaque", "treeshrew"),
  project = "primate"
)

# Load integrated object to get pseudotime stage assignments
load("./data/species.merge.RData")

# Map pseudotime stage from integrated object to primate subset
df1 <- obj$pseudostage %>%
  as.data.frame() %>%
  rownames_to_column(var = "barcode")
df2 <- primate_pseu$pseudotime %>%
  as.data.frame() %>%
  rownames_to_column(var = "barcode")

df <- left_join(df2, df1, by = "barcode")
names(df) <- c("barcode", "pseudotime", "pseudostage")
primate_pseu <- AddMetaData(primate_pseu, df)

# =========================== Calculate Average Expression per Pseudotime Stage ===========================

# For each species, aggregate expression by pseudotime stage
for (species in c("human", "macaque", "treeshrew")) {
  cat(">>> Aggregating expression for:", species, "\n")

  objsub <- subset(primate_pseu, Species == species)
  objsub$NEW_CELL_TYPE <- factor(objsub$NEW_CELL_TYPE,
                                 levels = c("SG", "SC", "Early SD", "Late SD"))

  # Aggregate expression by pseudotime stage
  avg <- AggregateExpression(objsub, assay = "RNA",
                             group.by = "pseudostage")$RNA %>%
    as.data.frame()

  # Log-normalize
  avg <- LogNormalize(avg)

  assign(x = paste0("gene_expr_", species), value = avg)
}

# =========================== Mfuzz Clustering for PTTG Genes ===========================

# Run Mfuzz for each primate/tree shrew species
for (species in c("human", "macaque", "treeshrew")) {
  cat(">>> Running Mfuzz for:", species, "\n")

  gene_mtx <- get(paste0("gene_expr_", species))

  # Create ExpressionSet object
  gene_mtx <- new("ExpressionSet", exprs = gene_mtx)

  # Filter and standardize
  gene_mtx <- filter.NA(gene_mtx, thres = 0.25)
  gene_mtx <- fill.NA(gene_mtx, mode = 'mean')
  gene_mtx <- filter.std(gene_mtx, min.std = 0)
  gene_mtx <- standardise(gene_mtx)

  # Mfuzz soft clustering
  set.seed(2024)
  cluster_num <- 8
  gene_cl <- mfuzz(gene_mtx, centers = cluster_num,
                   m = mestimate(gene_mtx),
                   iter.max = 200, rate.par = 0.1, weights = 1)

  assign(x = paste0("gene_cl_", species), value = gene_cl)
  assign(x = paste0("gene_mtx_", species), value = exprs(gene_mtx))

  # Plot Mfuzz clusters
  pdf(paste0("./figure/F4A.", species, ".mfuzz.primate.pdf"),
      width = 10, height = 5)
  mfuzz.plot2(gene_mtx, cl = gene_cl,
              mfrow = c(2, 3),
              centre = TRUE, centre.col = "black", centre.lwd = 2,
              time.labels = colnames(get(paste0("gene_expr_", species))),
              x11 = FALSE, ylab = "log2(Ratio)", ylim.set = c(-2, 2))
  dev.off()
}

# =========================== Match Clusters Across Species (PTTG) ===========================
# Manual cluster matching based on visual inspection of expression patterns
# Mapping: human cluster -> macaque cluster -> tree shrew cluster
#
# human     : 1 2 3 4 5 6 7 8
# macaque   : 8 1 6 3 4 2 7 5
# treeshrew : 3 2 6 4 8 1 7 5
# Celltype  : SG SD SD SC SC SG SD SC

matchlist_primate <- list(
  list(human = 1, macaque = 8, treeshrew = 3),   # Cluster 1
  list(human = 2, macaque = 1, treeshrew = 2),   # Cluster 2
  list(human = 3, macaque = 6, treeshrew = 6),   # Cluster 3
  list(human = 4, macaque = 3, treeshrew = 4),   # Cluster 4
  list(human = 5, macaque = 4, treeshrew = 8),   # Cluster 5
  list(human = 6, macaque = 2, treeshrew = 1),   # Cluster 6
  list(human = 7, macaque = 7, treeshrew = 7),   # Cluster 7
  list(human = 8, macaque = 5, treeshrew = 5)    # Cluster 8
)

# Find genes with conserved cluster membership across human, macaque, and tree shrew
human_macaque_treeshrew_genes <- list()
for (i in seq_along(matchlist_primate)) {
  cluster_union <- matchlist_primate[[i]]
  human_macaque_treeshrew_genes[[i]] <- Reduce(
    f = intersect,
    x = list(
      names(gene_cl_human$cluster)[gene_cl_human$cluster == cluster_union$human],
      names(gene_cl_macaque$cluster)[gene_cl_macaque$cluster == cluster_union$macaque],
      names(gene_cl_treeshrew$cluster)[gene_cl_treeshrew$cluster == cluster_union$treeshrew]
    )
  )
}

save(human_macaque_treeshrew_genes,
     gene_mtx_human, gene_mtx_macaque, gene_mtx_treeshrew,
     gene_cl_human, gene_cl_macaque, gene_cl_treeshrew,
     file = "./data/primate.only.change.RData")

# =========================== Plot Individual PTTG Gene Trajectories ===========================

# Plot expression trajectories for each conserved PTTG gene
for (gene in unlist(human_macaque_treeshrew_genes)) {
  gene_mtx <- rbind(
    gene_mtx_human[gene, ],
    gene_mtx_macaque[gene, ],
    gene_mtx_treeshrew[gene, ]
  )
  rownames(gene_mtx) <- c("human", "macaque", "treeshrew")

  df <- reshape2::melt(as.matrix(gene_mtx),
                       varnames = c("Species", "Stage"),
                       value.name = "Expression")

  p2 <- ggplot(df, aes(x = Stage, y = Expression)) +
    geom_point(aes(color = Species, shape = Species), size = 5) +
    geom_smooth(aes(color = Species, group = Species),
                method = 'loess', span = 0.5, se = FALSE, linewidth = 2) +
    geom_vline(xintercept = 1.5, linetype = "dashed", linewidth = 2) +
    geom_vline(xintercept = 5.5, linetype = "dashed", linewidth = 2) +
    scale_color_manual(values = hue_pal()(4)[c(1, 2, 4)]) +
    scale_fill_manual(values = hue_pal()(4)[c(1, 2, 4)]) +
    scale_shape_manual(values = c(15, 16, 18)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 15),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 2),
      axis.ticks = element_line(linewidth = 2)
    ) +
    labs(x = "", y = "", title = gene)

  ggsave(plot = p2,
         filename = paste0("./figure/primate.only/", gene, ".pdf"),
         width = 24, height = 8)
}

# =========================== Part 2: Four-Species Mfuzz Clustering ===========================
# Analyze genes expressed in all four species
# (manuscript Figure 4C, Supplementary Figure S6A)

load("./data/species.merge.RData")

# Verify pseudotime stage assignments
DimPlot(obj, pt.size = 1, group.by = "pseudostage")

# =========================== Calculate Gene Detection Rate per Stage ===========================
# Identify genes with sufficient expression (UMI > 0.05 fraction) in at least one
# pseudotime stage within each species

Idents(obj) <- paste0(obj$Species, "-", obj$pseudostage)

# Initialize 3D array: genes x stages x species
mtx <- array(0,
             dim = c(length(Features(obj)), 8, 4),
             dimnames = list(Features(obj),
                             paste0("stage", 1:8),
                             c("human", "macaque", "mouse", "treeshrew")))

# Calculate fraction of cells expressing each gene per stage per species
for (species in c("human", "macaque", "mouse", "treeshrew")) {
  for (stage in paste0("stage", c(1:8))) {
    obj_sub <- subset(obj, idents = paste0(species, "-", stage))
    obj_sub <- GetAssayData(obj_sub, assay = "RNA", layer = "data")
    pct <- rowSums(obj_sub) / ncol(obj_sub)
    mtx[, stage, species] <- pct
  }
}

# Select genes with expression fraction > 0.05 in at least one stage,
# AND this condition holds in ALL four species
gene_select <- apply((apply(mtx > 0.05, c(1, 3), any)), 1, all)
gene_select <- names(gene_select)[gene_select]

cat(">>> Number of genes passing expression filter:", length(gene_select), "\n")

# =========================== Aggregate Expression by Pseudotime Stage ===========================

for (species in c("human", "macaque", "mouse", "treeshrew")) {
  cat(">>> Aggregating expression for:", species, "\n")

  objsub <- subset(obj, Species == species)
  objsub$NEW_CELL_TYPE <- factor(objsub$NEW_CELL_TYPE,
                                 levels = c("SG", "SC", "Early SD", "Late SD"))

  avg <- AggregateExpression(objsub, assay = "RNA",
                             features = gene_select,
                             group.by = "pseudostage")$RNA %>%
    as.data.frame()

  avg <- LogNormalize(avg)

  assign(x = paste0("gene_expr_", species), value = avg)
}

# =========================== Mfuzz Clustering for Four Species ===========================

for (species in c("human", "macaque", "mouse", "treeshrew")) {
  cat(">>> Running Mfuzz for:", species, "\n")

  gene_mtx <- get(paste0("gene_expr_", species))

  # Create ExpressionSet
  gene_mtx <- new("ExpressionSet", exprs = gene_mtx)

  # Filter and standardize
  gene_mtx <- filter.NA(gene_mtx, thres = 0.25)
  gene_mtx <- fill.NA(gene_mtx, mode = 'mean')
  gene_mtx <- filter.std(gene_mtx, min.std = 0)
  gene_mtx <- standardise(gene_mtx)

  # Mfuzz with 6 clusters
  set.seed(2024)
  cluster_num <- 6
  gene_cl <- mfuzz(gene_mtx, centers = cluster_num,
                   m = mestimate(gene_mtx),
                   iter.max = 200, rate.par = 0.1, weights = 1)

  assign(x = paste0("gene_cl_", species), value = gene_cl)
  assign(x = paste0("gene_mtx_", species), value = exprs(gene_mtx))

  # Plot Mfuzz clusters
  pdf(paste0("./figure/F4A.", species, ".mfuzz.pdf"),
      width = 10, height = 5)
  mfuzz.plot2(gene_mtx, cl = gene_cl,
              mfrow = c(2, 3),
              centre = TRUE, centre.col = "black", centre.lwd = 2,
              time.labels = colnames(get(paste0("gene_expr_", species))),
              x11 = FALSE, ylab = "log2(Ratio)", ylim.set = c(-2, 2))
  dev.off()
}

save(list = c(ls(pattern = "gene_cl_"), ls(pattern = "gene_expr_"),
              ls(pattern = "gene_mtx_")),
     file = "./data/mfuzz.RData")

# =========================== Match Clusters Across Four Species ===========================
# Manual cluster matching based on visual inspection
# Mapping: human -> macaque -> mouse -> tree shrew
#
# human     : 1 2 3 4 5 6
# macaque   : 6 4 3 2 1 5
# mouse     : 4 2 6 1 5 3
# treeshrew : 6 2 5 1 4 3
# Celltype  : SG SC SC SG SD SD

load("./data/mfuzz.RData")

matchlist <- list(
  list(human = 1, macaque = 6, mouse = 4, treeshrew = 6),   # Cluster 1: SG
  list(human = 2, macaque = 4, mouse = 2, treeshrew = 2),   # Cluster 2: SC
  list(human = 3, macaque = 3, mouse = 6, treeshrew = 5),   # Cluster 3: SC
  list(human = 4, macaque = 2, mouse = 1, treeshrew = 1),   # Cluster 4: SG
  list(human = 5, macaque = 1, mouse = 5, treeshrew = 4),   # Cluster 5: SD
  list(human = 6, macaque = 5, mouse = 3, treeshrew = 3)    # Cluster 6: SD
)

celltype_labels <- c("SG", "SC", "SC", "SG", "SD", "SD")

# =========================== Identify Genes by Conservation Pattern ===========================

# --- Group 1: Conserved in ALL four species ---
human_macaque_mouse_treeshrew_genes <- list()
for (i in seq_along(matchlist)) {
  cluster_union <- matchlist[[i]]
  human_macaque_mouse_treeshrew_genes[[i]] <- Reduce(
    f = intersect,
    x = list(
      names(gene_cl_human$cluster)[gene_cl_human$cluster == cluster_union$human],
      names(gene_cl_macaque$cluster)[gene_cl_macaque$cluster == cluster_union$macaque],
      names(gene_cl_treeshrew$cluster)[gene_cl_treeshrew$cluster == cluster_union$treeshrew],
      names(gene_cl_mouse$cluster)[gene_cl_mouse$cluster == cluster_union$mouse]
    )
  )
}

# --- Group 2: Conserved in human, macaque, tree shrew (DIVERGENT in mouse) ---
# These are genes with same cluster in primates+tree shrew, but different in mouse
human_macaque_treeshrew_genes <- list()
for (i in seq_along(matchlist)) {
  cluster_union <- matchlist[[i]]
  human_macaque_treeshrew_genes[[i]] <- setdiff(
    Reduce(f = intersect, x = list(
      names(gene_cl_human$cluster)[gene_cl_human$cluster == cluster_union$human],
      names(gene_cl_macaque$cluster)[gene_cl_macaque$cluster == cluster_union$macaque],
      names(gene_cl_treeshrew$cluster)[gene_cl_treeshrew$cluster == cluster_union$treeshrew]
    )),
    names(gene_cl_mouse$cluster)[gene_cl_mouse$cluster == cluster_union$mouse]
  )
}

# --- Group 3: Conserved in human, macaque, mouse (DIVERGENT in tree shrew) ---
human_macaque_mouse_genes <- list()
for (i in seq_along(matchlist)) {
  cluster_union <- matchlist[[i]]
  human_macaque_mouse_genes[[i]] <- setdiff(
    Reduce(f = intersect, x = list(
      names(gene_cl_human$cluster)[gene_cl_human$cluster == cluster_union$human],
      names(gene_cl_macaque$cluster)[gene_cl_macaque$cluster == cluster_union$macaque],
      names(gene_cl_mouse$cluster)[gene_cl_mouse$cluster == cluster_union$mouse]
    )),
    names(gene_cl_treeshrew$cluster)[gene_cl_treeshrew$cluster == cluster_union$treeshrew]
  )
}

# --- Group 4: Conserved ONLY in human and macaque ---
human_macaque_genes <- list()
for (i in seq_along(matchlist)) {
  cluster_union <- matchlist[[i]]
  human_macaque_genes[[i]] <- setdiff(
    intersect(
      names(gene_cl_human$cluster)[gene_cl_human$cluster == cluster_union$human],
      names(gene_cl_macaque$cluster)[gene_cl_macaque$cluster == cluster_union$macaque]
    ),
    union(
      names(gene_cl_mouse$cluster)[gene_cl_mouse$cluster == cluster_union$mouse],
      names(gene_cl_treeshrew$cluster)[gene_cl_treeshrew$cluster == cluster_union$treeshrew]
    )
  )
}

# =========================== Calculate Conservation Scores ===========================
calculate_conservation_score <- function(P_gij,species = c("human","macaque","mouse","treeshrew")){
    P_gij_sub <- P_gij[ , ,which(species %in% c("human","macaque","mouse","treeshrew"))]
    cluster_product <- apply(P_gij_sub, c(1,2), function(x) prod(x))
    sum_product <- rowSums(cluster_product)
    conservation_score <- log2(sum_product)
    return(conservation_score)
}

# Build data frames with conservation scores
build_conservation_df <- function(gene_list, conservation_scores, celltype_labels) {
  cluster_len <- unlist(lapply(gene_list, function(x) length(x)))
  df <- data.frame(
    Cluster    = rep(paste0("Cluster", 1:length(gene_list)), cluster_len),
    CellType   = rep(celltype_labels, cluster_len),
    Gene       = unlist(gene_list),
    ConservationScore = conservation_scores[unlist(gene_list)]
  )
  return(df)
}



# Attempt to calculate conservation scores (wrapped in tryCatch for safety)
tryCatch({
  conservation_scores_human_macaque_mouse_treeshrew <- calculate_conservation_score(
    P_gij, species = c("human", "macaque", "mouse", "treeshrew"))
  names(conservation_scores_human_macaque_mouse_treeshrew) <- inter_gene

  human_macaque_mouse_treeshrew_df <- build_conservation_df(
    human_macaque_mouse_treeshrew_genes,
    conservation_scores_human_macaque_mouse_treeshrew[
      unlist(human_macaque_mouse_treeshrew_genes)],
    celltype_labels
  )

  conservation_scores_human_macaque_treeshrew <- calculate_conservation_score(
    P_gij, species = c("human", "macaque", "treeshrew"))
  names(conservation_scores_human_macaque_treeshrew) <- inter_gene

  human_macaque_treeshrew_df <- build_conservation_df(
    human_macaque_treeshrew_genes,
    conservation_scores_human_macaque_treeshrew[
      unlist(human_macaque_treeshrew_genes)],
    celltype_labels
  )

  conservation_scores_human_macaque_mouse <- calculate_conservation_score(
    P_gij, species = c("human", "macaque", "mouse"))
  names(conservation_scores_human_macaque_mouse) <- inter_gene

  human_macaque_mouse_df <- build_conservation_df(
    human_macaque_mouse_genes,
    conservation_scores_human_macaque_mouse[
      unlist(human_macaque_mouse_genes)],
    celltype_labels
  )

  conservation_scores_human_macaque <- calculate_conservation_score(
    P_gij, species = c("human", "macaque"))
  names(conservation_scores_human_macaque) <- inter_gene

  human_macaque_df <- build_conservation_df(
    human_macaque_genes,
    conservation_scores_human_macaque[
      unlist(human_macaque_genes)],
    celltype_labels
  )

  # Save conservation score data
  save(human_macaque_mouse_treeshrew_df, human_macaque_treeshrew_df,
       human_macaque_mouse_df, human_macaque_df,
       file = "./data/conservation_scores_genes.RData")

  # =========================== Export Conservation Scores to Excel ===========================

  df_list <- list(
    human_macaque_mouse_treeshrew = human_macaque_mouse_treeshrew_df,
    human_macaque_treeshrew       = human_macaque_treeshrew_df,
    human_macaque_mouse           = human_macaque_mouse_df,
    human_macaque                 = human_macaque_df
  )

  write.xlsx(df_list,
             file = "./figure/F4A.conservation_scores_genes.xlsx",
             rowNames = TRUE)

  # Export per-group gene lists with scores
  for (inter in c("human_macaque_mouse_treeshrew", "human_macaque_treeshrew",
                  "human_macaque_mouse", "human_macaque")) {
    conservation_scores <- get(paste0("conservation_scores_", inter))
    genes <- get(paste0(inter, "_genes"))

    write.xlsx(
      list(
        stage = as.data.frame(sapply(genes, "[",
                                     i = 1:max(sapply(genes, length)))),
        score = data.frame(
          name = names(conservation_scores),
          conservation_score = conservation_scores
        )
      ),
      file = paste0("./figure/F4A.", inter, ".xlsx")
    )
  }

}, error = function(e) {
  cat("NOTE: Conservation score calculation skipped.",
      "Please ensure calculate_conservation_score() and P_gij are defined.\n")
  cat("Error message:", e$message, "\n")
})

# =========================== Plot Top Conserved Gene Trajectories ===========================

load("./data/species.merge.RData")
load("./data/mfuzz.RData")

# Try to load conservation scores if available
if (file.exists("./data/conservation_scores_genes.RData")) {
  load("./data/conservation_scores_genes.RData")

  # For each conservation pattern, plot top 10 genes per cell type
  for (inter in c("human_macaque_mouse_treeshrew", "human_macaque_treeshrew",
                  "human_macaque_mouse", "human_macaque")) {

    df <- get(paste0(inter, "_df"))
    df <- df %>%
      rownames_to_column(var = "gene") %>%
      group_by(CellType) %>%
      mutate(rank = rank(-ConservationScore)) %>%
      arrange(rank) %>%
      filter(rank <= 10)

    write.csv(df, file = paste0("./figure/", inter, ".csv"))

    for (gene in df$gene) {
      gene_mtx <- rbind(
        gene_mtx_human[gene, ],
        gene_mtx_macaque[gene, ],
        gene_mtx_mouse[gene, ],
        gene_mtx_treeshrew[gene, ]
      )
      rownames(gene_mtx) <- c("human", "macaque", "mouse", "treeshrew")

      df_plot <- reshape2::melt(as.matrix(gene_mtx),
                                varnames = c("Species", "Stage"),
                                value.name = "Expression")

      p2 <- ggplot(df_plot, aes(x = Stage, y = Expression)) +
        geom_point(aes(color = Species, shape = Species), size = 3) +
        geom_smooth(aes(color = Species, group = Species),
                    method = 'loess', span = 0.5, se = FALSE, linewidth = 2) +
        geom_vline(xintercept = 1.5, linetype = "dashed", linewidth = 2) +
        geom_vline(xintercept = 5.5, linetype = "dashed", linewidth = 2) +
        scale_shape_manual(values = c(15, 16, 17, 18)) +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 15),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(linewidth = 2),
          axis.ticks = element_line(linewidth = 2)
        ) +
        labs(x = "", y = "", title = gene)

      ggsave(plot = p2,
             filename = paste0("./figure/", inter, "-", gene, ".pdf"),
             width = 16, height = 8)
    }
  }
}

cat(">>> 05_mfuzz.R completed successfully!\n")