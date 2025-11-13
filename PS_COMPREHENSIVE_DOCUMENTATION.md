# PS (Perturbation-response Score) Analysis - Comprehensive Documentation

**Version**: Tutorial Repository
**Package Used**: scMAGeCK (R package)
**Method**: Perturbation-response Score (PS) Analysis
**Authors**: Bicna Song, Wei Li (Children's National Hospital)
**Date**: 2025-11-13

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Repository Overview](#repository-overview)
3. [Workflow 1: Demo 1 - Simple PS Example](#workflow-1-demo-1---simple-ps-example)
4. [Workflow 2: Demo 2 - BeeSTING-seq Analysis](#workflow-2-demo-2---beesting-seq-analysis)
5. [Workflow 3: HIV Perturb-seq Dataset](#workflow-3-hiv-perturb-seq-dataset)
6. [Workflow 4: Pancreatic Differentiation Dataset](#workflow-4-pancreatic-differentiation-dataset)
7. [scMAGeCK Functions Reference](#scmageck-functions-reference)
8. [Data Formats and Requirements](#data-formats-and-requirements)
9. [Complete Analysis Pipeline](#complete-analysis-pipeline)
10. [Cross-Reference Tables](#cross-reference-tables)

---

## Executive Summary

This repository contains **tutorials and analysis workflows** for calculating Perturbation-response Score (PS), a method to quantify diverse perturbation responses and discover novel biological insights in single-cell perturbation datasets.

### What is PS Analysis?

PS (Perturbation-response Score) quantifies how individual cells respond to perturbations (genetic, chemical, environmental, or mechanical). Unlike traditional differential expression analysis that averages responses across populations, PS captures **heterogeneous** responses at single-cell resolution.

### Key Capabilities

1. **Quantify perturbation effects** at single-cell resolution
2. **Identify heterogeneous responses** within perturbed cell populations
3. **Discover subpopulations** with distinct perturbation responses
4. **Visualize response patterns** across cell types and states

### Repository Statistics

| Metric | Count |
|--------|-------|
| Total Workflows | 4 (2 demos + 2 datasets) |
| Total Lines of Code | 1,085 lines |
| R Script Files | 1 (.R) |
| R Markdown Files | 3 (.Rmd) |
| Demo Datasets | 2 |
| Real Dataset Analyses | 2 |

### Key Technologies

- **R Package**: scMAGeCK (implements PS method)
- **Framework**: Seurat (single-cell RNA-seq analysis)
- **Data Types**: Perturb-seq, CRISPR screens, base editing screens
- **Visualization**: UMAP, feature plots, violin plots

---

## Repository Overview

### Purpose

This repository demonstrates how to:

1. Calculate PS scores from single-cell perturbation data
2. Process data from 10X Genomics platforms
3. Integrate PS analysis with Seurat workflows
4. Visualize and interpret perturbation heterogeneity
5. Perform downstream analyses on PS-scored cells

### Repository Structure

```
PS/
├── README.md                          # Main documentation
├── demo/
│   ├── demo1/
│   │   ├── ps_demo.R                 # Simple PS example (47 lines)
│   │   ├── barcode_rec.txt           # Barcode file
│   │   └── singles_dox_mki67_v3.RDS  # Mini Seurat object
│   └── demo2/
│       ├── PS_demo2.Rmd              # BeeSTING-seq workflow (175 lines)
│       ├── PS_demo2.nb.html          # Compiled notebook
│       └── counts/                    # 10X count matrices
└── datasets/
    ├── HIV_Perturb-seq/
    │   ├── hiv_perturbseq.Rmd        # HIV screen analysis (326 lines)
    │   ├── BARCODE_H13Ld2EGFP.txt    # Barcode file
    │   └── JKLAT_H13Ld2EGFP.combined.rds  # Seurat object
    └── Pancreatic_differentiation_scRNA-seq/
        ├── pancreatic_scrnaseq.Rmd   # Pancreatic analysis (537 lines)
        └── 10clones_seurat.rds        # Seurat object
```

### Workflows Summary

| Workflow | File | Lines | Purpose | Dataset Type |
|----------|------|-------|---------|--------------|
| Demo 1 | ps_demo.R | 47 | Quick start example | Mini Perturb-seq |
| Demo 2 | PS_demo2.Rmd | 175 | Full 10X pipeline | BeeSTING-seq (base editing) |
| HIV Perturb-seq | hiv_perturbseq.Rmd | 326 | Publication dataset | CRISPR Perturb-seq |
| Pancreatic Diff | pancreatic_scrnaseq.Rmd | 537 | Cell differentiation | Lineage tracing + CRISPR |

### When to Use Each Workflow

**Use Demo 1 if:**
- You want a quick 5-minute introduction
- You have pre-processed Seurat objects ready
- You just want to see PS calculation basics

**Use Demo 2 if:**
- You're starting from 10X count matrices
- You need the complete preprocessing pipeline
- You want to learn Seurat integration steps

**Use HIV Perturb-seq workflow if:**
- You're analyzing CRISPR Perturb-seq data
- You want to see publication-quality analysis
- You need examples of multiple perturbation comparisons

**Use Pancreatic workflow if:**
- You're studying cell differentiation contexts
- You have perturbations in multiple cell types
- You need cell-type-specific PS calculations

---

## Workflow 1: Demo 1 - Simple PS Example

**File**: `demo/demo1/ps_demo.R`
**Lines of Code**: 47
**Level**: Beginner
**Time to Run**: 2-5 minutes
**Data Type**: Mini Perturb-seq dataset (TP53 knockout)

### Purpose

This is the **quickest way** to learn PS calculation. It demonstrates the minimal code needed to:
1. Load a Seurat object and barcode file
2. Calculate PS scores for a single gene (TP53)
3. Visualize PS scores vs gene expression

### When to Use This Workflow

✅ **Use if:**
- You're new to PS analysis
- You have pre-processed Seurat objects
- You want a quick proof-of-concept
- You're teaching/learning the basics

❌ **Don't use if:**
- You need to start from raw count matrices
- You need comprehensive preprocessing steps
- You're analyzing multiple perturbations

### Algorithm - Step by Step

#### Step 1: Load Libraries (lines 1-2)
```r
library(scMAGeCK)  # PS calculation functions
library(Seurat)    # Single-cell analysis framework
```

#### Step 2: Load and Prepare Barcode File (lines 3-15)

**Lines 3-5**: Load barcode table
```r
BARCODE = "barcode_rec.txt"
bc_frame = read.table(BARCODE, header = T, as.is = T)
```

**Lines 7-11**: Alternative method (commented) - convert guide matrix to barcode format
```r
# If you have a guide expression matrix instead of a barcode file:
# bc_frame = guidematrix_to_triplet(rds_object[['CRISPR']]@counts, sobj)
# bc_frame[,'sgrna'] = bc_frame[,'barcode']
# bc_frame[,'gene'] = sub('-[0-9]+$', '', bc_frame[,'barcode'])
```

**Lines 14-15**: Clean cell identifiers
```r
bc_frame$cell = sub('-1', '', bc_frame$cell)  # Remove Seurat suffix
```

#### Step 3: Load Seurat Object (lines 18-20)

```r
RDS = "singles_dox_mki67_v3.RDS"
rds_object = readRDS(RDS)
```

**Input**: Pre-processed Seurat object containing:
- Normalized gene expression matrix
- Dimensionality reduction (PCA, t-SNE/UMAP)
- Cell metadata

#### Step 4: Calculate PS Scores (lines 22-27)

**Lines 26-27**: Core PS calculation
```r
eff_object <- scmageck_eff_estimate(
  rds_object,               # Seurat object with scRNA-seq data
  bc_frame,                 # Barcode table linking cells to guides
  perturb_gene = 'TP53',    # Gene to calculate PS for
  non_target_ctrl = 'NonTargetingControlGuideForHuman'  # Negative control
)
```

**What happens internally:**
1. Identifies cells with TP53-targeting guides
2. Identifies cells with non-targeting control guides
3. Finds differentially expressed genes between perturbation and control
4. Calculates PS score for each cell based on perturbation-response signature
5. Adds PS scores to Seurat metadata

#### Step 5: Alternative - Calculate PS for All Genes (lines 29-35, commented)

```r
# To analyze multiple genes:
# eff_object <- scmageck_eff_estimate(
#   rds_object, bc_frame,
#   perturb_gene = grep('NonTargetingControlGuideForHuman',
#                       unique(rds_object$gene),
#                       value = T, invert = T),  # All genes except NTC
#   non_target_ctrl = 'NonTargetingControlGuideForHuman'
# )
```

#### Step 6: Extract Results (lines 38-39)

```r
eff_estimat = eff_object$eff_matrix    # PS score matrix
rds_subset = eff_object$rds            # Updated Seurat object
```

**Return value structure:**
- `eff_matrix`: Matrix of PS scores (cells × genes)
- `rds`: Seurat object with new metadata column: `{gene}_eff` (e.g., `TP53_eff`)

#### Step 7: Visualize Results (lines 41-45)

**Line 42**: Visualize PS scores
```r
FeaturePlot(rds_subset, features = 'TP53_eff', reduction = 'tsne')
```
- Shows spatial distribution of TP53 perturbation effects
- High scores = strong TP53 perturbation response
- Reveals response heterogeneity

**Line 45**: Compare to gene expression
```r
FeaturePlot(rds_subset, features = 'TP53', reduction = 'tsne')
```
- Shows TP53 mRNA expression
- **Key insight**: PS scores often reveal patterns NOT visible in gene expression

### Input Requirements

#### 1. Barcode File (`barcode_rec.txt`)

**Format**: Tab-separated text file with required columns:

| Column | Type | Description |
|--------|------|-------------|
| cell | character | Cell barcode (e.g., "AAATCAACGGGTGA-1") |
| barcode | character | Guide barcode identifier |
| sgrna | character | sgRNA sequence or identifier |
| gene | character | Target gene name |
| read_count | integer | Number of reads for this guide in this cell |
| umi_count | integer | Number of UMIs (unique molecular identifiers) |

**Example:**
```
cell                barcode         sgrna                   gene    read_count  umi_count
AAATCAACGGGTGA-1   NF1_sg_118      AGTCAGTACTGAGCACAACA    NF1     1187        30
AAATCAACGGGTGA-1   CDKN2A_sg_70    TCTTGGTGACCCTCCGGATT    CDKN2A  1           1
```

**Important notes:**
- One row per guide detected in a cell
- Cells with multiple guides have multiple rows
- Cell barcodes must match Seurat object exactly

#### 2. Seurat Object (`singles_dox_mki67_v3.RDS`)

**Required components:**
- Normalized RNA expression matrix
- Dimensionality reduction (PCA, t-SNE, or UMAP)
- Cell metadata (optional but recommended)

**Preprocessing already done:**
- QC filtering
- Normalization
- Variable feature selection
- Scaling
- Dimensionality reduction

### Output

#### 1. PS Score Matrix (`eff_estimat`)
- **Type**: Matrix (cells × genes)
- **Range**: Typically 0 to 1+ (higher = stronger response)
- **Interpretation**:
  - 0-0.2: Weak/no response
  - 0.2-0.5: Moderate response
  - 0.5+: Strong response

#### 2. Updated Seurat Object (`rds_subset`)
- **New metadata column**: `TP53_eff` (PS scores)
- **Subset**: Only cells with guides (perturbation or control)
- **Use for**: Visualization, downstream analysis, differential expression

### Expected Results

**From line 42 visualization:**
- TP53 PS scores show clear clustering patterns
- High PS score regions indicate strong TP53 knockout effects
- Spatial patterns reveal functional heterogeneity

**Key Discovery (comparing lines 42 vs 45):**
- PS scores reveal perturbation effects that gene expression cannot
- TP53 expression alone doesn't show response heterogeneity
- PS captures transcriptome-wide perturbation signatures

### Usage Example - Modified

```r
# Calculate PS for a different gene
eff_object <- scmageck_eff_estimate(
  rds_object, bc_frame,
  perturb_gene = 'KRAS',  # Different target
  non_target_ctrl = 'NonTargetingControlGuideForHuman'
)

# Extract and visualize
rds_subset = eff_object$rds
FeaturePlot(rds_subset, features = 'KRAS_eff', reduction = 'tsne')
```

### Dependencies

#### External Packages
- **scMAGeCK**: Core PS calculation (`scmageck_eff_estimate`)
- **Seurat**: Single-cell framework (`FeaturePlot`, `readRDS`)

#### Data Files
- `barcode_rec.txt`: Guide assignment table (546 KB)
- `singles_dox_mki67_v3.RDS`: Seurat object (3 MB)

### Common Issues and Solutions

**Issue 1**: Cell barcode mismatch
```r
# Problem: Seurat cell names have "-1" suffix, barcode file doesn't
# Solution (line 15):
bc_frame$cell = sub('-1', '', bc_frame$cell)
```

**Issue 2**: Converting from guide matrix
```r
# If your data is in matrix format (guides × cells):
bc_frame = guidematrix_to_triplet(guide_matrix, seurat_obj)
bc_frame[,'sgrna'] = bc_frame[,'barcode']
bc_frame[,'gene'] = sub('-[0-9]+$', '', bc_frame[,'barcode'])
```

**Issue 3**: Multiple perturbations
```r
# Calculate PS for multiple genes at once:
gene_list = c('TP53', 'KRAS', 'EGFR')
eff_object <- scmageck_eff_estimate(
  rds_object, bc_frame,
  perturb_gene = gene_list,
  non_target_ctrl = 'NonTargetingControlGuideForHuman'
)
```

---

## Workflow 2: Demo 2 - BeeSTING-seq Analysis

**File**: `demo/demo2/PS_demo2.Rmd`
**Lines of Code**: 175
**Level**: Intermediate
**Time to Run**: 15-30 minutes
**Data Type**: BeeSTING-seq (base editing Perturb-seq from 10X Genomics 5' protocol)
**Dataset Source**: Morris et al. Science 2023, GEO: GSE171452

### Purpose

This workflow demonstrates the **complete end-to-end pipeline** from raw 10X count matrices to PS score visualization. It includes:

1. Reading 10X count matrices (gene expression + guide counts)
2. Creating Seurat objects with both RNA and CRISPR assays
3. Standard Seurat preprocessing (QC, normalization, clustering, UMAP)
4. Converting guide matrix to barcode format
5. Preprocessing guide expression matrix
6. Calculating PS scores
7. Downstream analysis and visualization

### When to Use This Workflow

✅ **Use if:**
- Starting from 10X Genomics count matrices
- Need complete preprocessing pipeline
- Working with dual-assay data (RNA + CRISPR)
- Want to learn full Seurat integration
- Analyzing base editing screens

❌ **Don't use if:**
- You already have preprocessed Seurat objects (use Demo 1)
- You're not using 10X data format
- You need specialized QC thresholds

### Algorithm - Step by Step

#### Step 1: Load Libraries (lines 15-22, 36-37)

```r
library(Seurat)      # Single-cell analysis framework
library(ggplot2)     # Plotting
library(patchwork)   # Combine plots
library(scales)      # Scale functions
library(dplyr)       # Data manipulation
library(scMAGeCK)    # PS calculation
```

**Optional installation** (lines 28-32):
```r
library(devtools)
install_github('weililab/scMAGeCK')
```

#### Step 2: Read 10X Count Matrices (lines 43-51)

**Lines 45-47**: Read 10X matrix format
```r
exp_mat_GDO = ReadMtx(
  mtx = 'counts/GSM7108136_BeeSTINGseq_GDO-A_matrix.mtx.gz',
  cells = 'counts/GSM7108136_BeeSTINGseq_GDO-A_barcodes.tsv.gz',
  features = 'counts/GSM7108136_BeeSTINGseq_GDO-A_features.tsv.gz'
)
```

**What this does:**
- Reads sparse matrix format (efficient for large datasets)
- `mtx`: Count matrix (features × cells)
- `cells`: Cell barcodes
- `features`: Gene/guide names

**Lines 49-50**: Separate RNA and CRISPR data
```r
sobj = CreateSeuratObject(counts = exp_mat_GDO[1:36601,])    # Genes
sobj[['CRISPR']] = CreateAssayObject(counts = exp_mat_GDO[36602:36939,])  # Guides
```

**Key insight**: 10X combines both assays in one matrix
- Rows 1-36601: Gene expression features
- Rows 36602-36939: CRISPR guide features (338 guides)

#### Step 3: Quality Control (lines 58-65)

**Line 59**: Calculate mitochondrial percentage
```r
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
```

**Lines 60-64**: Visualize QC metrics
```r
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

**Metrics to check:**
- `nFeature_RNA`: Number of genes detected per cell
- `nCount_RNA`: Total UMI counts per cell
- `percent.mt`: Mitochondrial gene percentage (high = dying cells)

**Line 66** (commented): Optional filtering
```r
# sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 10)
```

#### Step 4: Seurat Preprocessing Pipeline (lines 68-77)

**Line 68**: Normalize, find variable features, scale (chained operations)
```r
sobj <- NormalizeData(object = sobj) %>% 
        FindVariableFeatures() %>% 
        ScaleData()
```

**Details:**
- `NormalizeData`: Log-normalization (default: scale factor 10,000)
- `FindVariableFeatures`: Identify most variable genes (default: 2,000)
- `ScaleData`: Z-score transformation

**Line 70**: Principal component analysis
```r
sobj <- RunPCA(object = sobj)
```
- Reduces dimensionality for downstream analysis
- Uses variable features by default

**Line 72**: UMAP dimensionality reduction
```r
sobj <- RunUMAP(object = sobj, dims = 1:10)
```
- Uses first 10 PCs
- Creates 2D representation for visualization

**Lines 75-77**: Clustering
```r
sobj <- FindNeighbors(sobj, dims = 1:10)
sobj <- FindClusters(sobj, resolution = 0.5)
DimPlot(sobj, reduction = "umap")
```
- Constructs k-nearest neighbor graph
- Identifies clusters using Louvain algorithm
- Resolution 0.5 = moderate cluster granularity

#### Step 5: Prepare Barcode File (lines 90-100)

**Line 93**: Convert guide matrix to barcode format
```r
bc_frame = guidematrix_to_triplet(sobj[['CRISPR']]@counts, sobj)
```

**What this does:**
- Converts sparse matrix (guides × cells) to long format
- Each row = one guide detected in one cell
- Includes read counts and UMI counts

**Lines 94-95**: Add required columns
```r
bc_frame[,'sgrna'] = bc_frame[,'barcode']  # Use barcode as sgRNA ID
bc_frame[,'gene'] = sub('-[0-9]+$', '', bc_frame[,'barcode'])  # Extract gene name
```

**Example transformation:**
- Barcode: `CD55-Q86X-1`
- Gene extracted: `CD55-Q86X`

**Line 99**: Filter low-quality guides
```r
bc_frame <- bc_frame[bc_frame$read_count > 1,]
```
- Removes guides with only 1 read (likely noise)
- Improves data quality

#### Step 6: Pre-process Guide Expression (lines 107-114)

**Line 109**: Switch to CRISPR assay
```r
DefaultAssay(sobj) = 'CRISPR'
```

**Line 111**: Pre-process RDS with guide information
```r
sobj <- pre_processRDS(bc_frame, sobj)
```

**What this does:**
- Adds guide metadata to Seurat object
- Identifies cells with single guides (singlets)
- Calculates guide expression metrics
- Prepares object for PS calculation

**Line 113**: Visualize guide distribution
```r
FeaturePlot(sobj, features = c('NTC-GV2','CD55-Q86X'))
```
- `NTC-GV2`: Non-targeting control
- `CD55-Q86X`: CD55 stop codon mutation

#### Step 7: Calculate PS Scores (lines 120-131)

**Line 121**: Switch back to RNA assay
```r
DefaultAssay(sobj) <- 'RNA'
```

**Lines 125-126**: Run PS calculation
```r
eff_object <- scmageck_eff_estimate(
  sobj, bc_frame,
  perturb_gene = c('SNP-36'),       # Target perturbation
  non_target_ctrl = 'NTC-GV2',      # Negative control
  subset_rds = T,                    # Subset to guide-containing cells only
  lambda = 0,                        # No background correction
  target_gene_max = 100              # Maximum target genes to use
)
```

**Parameters explained:**
- `perturb_gene`: Which perturbation to calculate PS for
- `non_target_ctrl`: Control cells for comparison
- `subset_rds = T`: Return only cells with guides
- `lambda = 0`: Disable background correlation correction
- `target_gene_max = 100`: Use top 100 differentially expressed genes

**Lines 128-129**: Extract results
```r
eff_estimat = eff_object$eff_matrix
rds_subset = eff_object$rds
```

#### Step 8: Visualize PS Scores (lines 136-140)

**Line 137**: Visualize PS scores
```r
FeaturePlot(rds_subset, features = 'SNP.36_eff')
```
- Note: Hyphens converted to periods in feature names
- Shows spatial distribution of perturbation response

**Line 138**: Compare to target gene expression
```r
FeaturePlot(rds_subset, features = 'APPBP2')
```

#### Step 9: Downstream Analysis (lines 144-161)

**Line 146**: Extract data for analysis
```r
data_f = FetchData(rds_subset, vars = c('sgrna','gene','APPBP2','SNP.36_eff'))
```

**Line 150**: Categorize cells by PS score
```r
data_f$PS = ifelse(data_f$gene == 'NTC-GV2', 'NTC',
                   ifelse(data_f$SNP.36_eff > 0.5, 'PS_High', 'PS_Low'))
```

**Categorization logic:**
- `NTC`: Non-targeting control cells
- `PS_High`: Strong responders (PS > 0.5)
- `PS_Low`: Weak responders (PS ≤ 0.5)

**Lines 153-157**: Violin plot comparison
```r
p <- ggplot(data_f, aes(x = PS, y = APPBP2)) + 
  geom_violin() + 
  theme_bw()
print(p)
```

**Line 167-172**: ECDF plot (empirical cumulative distribution)
```r
plot(ecdf(data_f$APPBP2[data_f$gene == 'NTC-GV2']),
     ylim = c(0.9,1), xlab = 'APPBP2 expression', ylab = 'Fraction')
plot(ecdf(data_f$APPBP2[data_f$SNP.36_eff < 0.2]), add = T, col = 'blue')
plot(ecdf(data_f$APPBP2[data_f$SNP.36_eff > 0.2]), add = T, col = 'red')
legend('bottomright', c('NT','PS low','PS high'), col = c('black','blue','red'), lwd = 1, pch = 20)
```

**What this shows:**
- Distribution shift in target gene expression
- Separation between high/low PS score groups
- Effect size visualization

### Input Requirements

#### 1. 10X Count Matrices (in `counts/` directory)

Three files per sample:
- `*_matrix.mtx.gz`: Sparse count matrix (Market Matrix format)
- `*_barcodes.tsv.gz`: Cell barcodes (one per line)
- `*_features.tsv.gz`: Feature names (genes + guides)

**Feature file structure:**
```
ENSEMBL_ID    GENE_NAME    feature_type
...
Guide-1       Guide-1      CRISPR Guide
Guide-2       Guide-2      CRISPR Guide
```

#### 2. Data Characteristics

**From this dataset:**
- Total features: 36,939 (36,601 genes + 338 guides)
- Feature types: Mixed (Gene Expression + CRISPR Guide Capture)
- Platform: 10X Genomics 5' chemistry
- Protocol: BeeSTING-seq (base editing screen)

### Output

#### 1. Processed Seurat Object
- **Assays**: RNA (gene expression) + CRISPR (guide counts)
- **Reductions**: PCA, UMAP
- **Metadata**: Clusters, QC metrics, guide assignments, PS scores

#### 2. PS Scores
- **Format**: Metadata column `SNP.36_eff`
- **Range**: Continuous scores (typically 0-1+)
- **Subset**: Only cells with guide assignments

#### 3. Visualizations
- QC plots (violin, scatter)
- UMAP with clusters
- Guide distribution maps
- PS score feature plots
- Violin plots (PS groups vs target gene)
- ECDF plots (cumulative distributions)

### Key Differences from Demo 1

| Aspect | Demo 1 | Demo 2 |
|--------|--------|--------|
| Starting point | Preprocessed RDS | Raw 10X matrices |
| Data format | Single assay | Dual assay (RNA + CRISPR) |
| QC steps | None | Full pipeline |
| Normalization | Pre-done | Shown explicitly |
| Guide format | Barcode table | Matrix (converted) |
| Complexity | 47 lines | 175 lines |
| Time to run | 2-5 min | 15-30 min |

### Usage Example - Adapt to Your Data

```r
# Read your 10X data
exp_mat = ReadMtx(
  mtx = 'your_data/matrix.mtx.gz',
  cells = 'your_data/barcodes.tsv.gz',
  features = 'your_data/features.tsv.gz'
)

# Determine split point (check your feature file)
n_genes = 20000  # Adjust based on your data
n_guides = 100   # Adjust based on your data

sobj = CreateSeuratObject(counts = exp_mat[1:n_genes,])
sobj[['CRISPR']] = CreateAssayObject(counts = exp_mat[(n_genes+1):(n_genes+n_guides),])

# Continue with standard workflow...
```

### Dependencies

#### R Packages
- **Seurat** (≥4.0): Single-cell framework
- **scMAGeCK**: PS calculation
- **ggplot2**: Visualization
- **patchwork**: Combine plots
- **dplyr**: Data manipulation
- **scales**: Scaling functions

#### System Requirements
- R version ≥4.0
- Sufficient memory for count matrix (~2-4 GB for this dataset)

### Common Issues and Solutions

**Issue 1**: Can't determine split between genes and guides
```r
# Solution: Check feature file
features = read.table('counts/*_features.tsv.gz', sep = '\t')
table(features$V3)  # Shows feature types
```

**Issue 2**: Out of memory
```r
# Solution: Use sparse matrices, don't convert to dense
library(Matrix)  # Already loaded by Seurat
# Avoid: as.matrix(large_sparse_matrix)
```

**Issue 3**: PS scores all zeros or NAs
```r
# Check guide assignments
table(sobj$gene)  # Should show guide distribution
table(is.na(sobj$gene))  # Check for NAs

# Ensure you have enough cells per condition
table(bc_frame$gene)  # Need >10 cells per perturbation
```

**Issue 4**: Hyphen vs period in feature names
```r
# Seurat converts hyphens to periods in metadata
# 'SNP-36' becomes 'SNP.36_eff'
# Use: make.names() to predict the conversion
make.names('SNP-36')  # Returns 'SNP.36'
```

### Best Practices

1. **QC Filtering**: Always inspect QC metrics before filtering
2. **Guide Quality**: Filter guides with very low read counts
3. **Normalization**: Use standard Seurat pipeline for consistency
4. **Control Cells**: Ensure sufficient non-targeting control cells (>50)
5. **Target Genes**: Adjust `target_gene_max` based on effect size
6. **Visualization**: Always compare PS scores to gene expression

---

