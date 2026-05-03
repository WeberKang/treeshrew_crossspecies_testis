# TreeShrew Cross-Species Single-Cell Transcriptomics

Cross-species single-cell transcriptomic analysis of spermatogenesis in human, macaque, tree shrew, and mouse.

## Pipeline

| Script | Analysis |
|--------|----------|
| `01_data_processing.R` | Multiome data loading, QC, WNN integration |
| `02_treeshrew_exp_atac.R` | Gene activity scores, differential peaks, ChIPseeker |
| `03_cross_species_integration.R` | CCA integration, PCA, correlations, LISI |
| `04_wgcna.R` | hdWGCNA co-expression modules |
| `05_mfuzz.R` | Soft clustering, conservation scores |
| `06_sg.R` | Spermatogonia subtypes (SSC/Undiff/Diff) |
| `07_scsd.R` | Meiotic/post-meiotic subtypes |

## Species

- Human (*Homo sapiens*)
- Rhesus macaque (*Macaca mulatta*)
- Chinese tree shrew (*Tupaia belangeri chinensis*)
- Mouse (*Mus musculus*)

## Requirements

**R packages:**
```r
Seurat, Signac, tidyverse, monocle3, cellAlign, hdWGCNA, Mfuzz, ComplexHeatmap, ChIPseeker, biomaRt
```

**Python environment:**
```bash
conda create -n Renv python=3.9
conda activate Renv
pip install leidenalg numpy umap-learn
```

## Execution

```bash
Rscript 01_data_processing.R
Rscript 02_treeshrew_exp_atac.R
Rscript 03_cross_species_integration.R
Rscript 04_wgcna.R
Rscript 05_mfuzz.R
Rscript 06_sg.R
Rscript 07_scsd.R
```

## Data availability
The raw sequencing data reported in this study were deposited in the National Genomics Data Center (NGDC) under accession code [PRJCA047161](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA047161).

## Citation

```bibtex

```

## Contact

Wei-Bo Kang:kangweibo@mail.kiz.ac.cn

## License

GPL-3.0
