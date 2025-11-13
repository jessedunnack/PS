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

## Workflow 3: HIV Perturb-seq Dataset

**File**: `datasets/HIV_Perturb-seq/hiv_perturbseq.Rmd`
**Lines of Code**: 326
**Level**: Advanced
**Time to Run**: 1-2 hours
**Data Type**: CRISPR Perturb-seq (HIV latency screen)
**Publication**: HIV reactivation screen studying BRD4, CCNT1, and other targets
**Authors**: Bicna Song, Wei Li (Children's National Hospital)

### Purpose

This workflow demonstrates **publication-quality analysis** of a real Perturb-seq dataset studying HIV latency. It showcases:

1. Working with pre-integrated multi-sample Seurat objects
2. Comprehensive QC and clustering analysis
3. Multi-perturbation PS score calculation
4. Advanced downstream analysis (differential expression, signature scoring)
5. Identifying perturbation-responsive subpopulations
6. Publication-quality figure generation

### When to Use This Workflow

✅ **Use if:**
- Analyzing multi-gene CRISPR screens
- Working with pre-integrated Seurat objects
- Need publication-quality analysis examples
- Studying biological phenotypes (e.g., HIV reactivation)
- Want advanced downstream analysis patterns

❌ **Don't use if:**
- You need basic preprocessing (see Demo 2)
- You're new to PS analysis (start with Demo 1)
- You need cell-type-specific analysis (see Workflow 4)

### Scientific Context

**Research Question**: Which genes regulate HIV latency reversal?
**Approach**: CRISPR knockout screen with HIV-GFP reporter
**Key Findings**:
- BRD4 knockout induces HIV reactivation (high GFP) in subset of cells
- CCNT1 knockout shows similar effect
- PS scores identify specific responder subpopulation (cluster 8)

### Algorithm - Step by Step

#### Step 1: Setup and Load Data (lines 18-27)

```r
library(Seurat)
library(scMAGeCK)
library(ggplot2)
library(hdf5r)

feat_c = readRDS(file = "JKLAT_H13Ld2EGFP.combined.rds")
BARCODE = 'BARCODE_H13Ld2EGFP.txt'
```

**Input**: Pre-integrated Seurat object
- Multiple samples combined
- Already normalized and integrated
- Contains HIV-GFP reporter expression

#### Step 2: Visualize Guide Distribution (lines 34-39)

**Line 34**: Visualize guide distribution across cells
```r
featurePlot(RDS = feat_c, BARCODE = BARCODE, TYPE = "Dis")
```

**Line 38**: Add guide metadata to Seurat object
```r
feat_c = pre_processRDS(BARCODE = BARCODE, RDS = feat_c)
```

**What this does:**
- Parses barcode file
- Adds guide identity to cell metadata
- Calculates guide UMI counts
- Identifies singlets (cells with 1 guide) vs multiplets

#### Step 3: Filter for Singlet Cells (lines 41-43)

```r
feat_c_singlet = subset(x = feat_c, subset = nFeature_sgRNA_guides == 1)
```

**Why singlets only?**
- Cells with 1 guide = unambiguous perturbation identity
- Cells with multiple guides = confounded effects
- Standard practice in Perturb-seq analysis

#### Step 4: Quality Control (lines 46-60)

**Lines 47-48**: Visualize QC metrics by sample
```r
feat_c_singlet <- SetIdent(feat_c_singlet, value = "orig.ident")
VlnPlot(feat_c_singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

**Line 53**: Apply QC filters
```r
feat_c_singlet <- subset(feat_c_singlet, subset = 
  nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
```

**Filters explained:**
- `nFeature_RNA > 200`: Remove low-quality cells (too few genes)
- `nFeature_RNA < 7500`: Remove doublets (too many genes)
- `percent.mt < 15`: Remove dying cells (high mitochondrial %)

**Lines 58-60**: Recheck QC after filtering
```r
VlnPlot(feat_c_singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

#### Step 5: Standard Seurat Workflow (lines 73-133)

**Lines 75-76**: Normalization and variable feature selection
```r
feat_c_singlet <- NormalizeData(feat_c_singlet, normalization.method = "LogNormalize", scale.factor = 10000)
feat_c_singlet <- FindVariableFeatures(feat_c_singlet, selection.method = "vst", nfeatures = 2000)
```

**Lines 78-87**: Identify and plot top variable genes
```r
top10 <- head(VariableFeatures(feat_c_singlet), 10)
plot1 <- VariableFeaturePlot(feat_c_singlet)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
```

**Line 92**: Scale data
```r
feat_c_singlet <- ScaleData(feat_c_singlet)
```

**Lines 95-96**: PCA
```r
feat_c_singlet <- RunPCA(feat_c_singlet, features = VariableFeatures(object = feat_c_singlet))
print(feat_c_singlet[["pca"]], dims = 1:5, nfeatures = 5)
```

**Lines 111-112**: Determine dimensionality with JackStraw
```r
feat_c_singlet <- JackStraw(feat_c_singlet, num.replicate = 100)
feat_c_singlet <- ScoreJackStraw(feat_c_singlet, dims = 1:20)
```

**JackStraw analysis**: Statistical test to determine significant PCs
- Resamples data 100 times
- Tests significance of each PC
- Guides choice of dimensions for clustering

**Lines 123-125**: Clustering
```r
feat_c_singlet <- FindNeighbors(feat_c_singlet, dims = 1:10)
feat_c_singlet <- FindClusters(feat_c_singlet, resolution = 0.5)
head(Idents(feat_c_singlet), 5)
```

**Lines 127-132**: UMAP
```r
feat_c_singlet <- RunUMAP(feat_c_singlet, dims = 1:10)
DimPlot(feat_c_singlet, reduction = "umap")
DimPlot(feat_c_singlet, reduction = "umap", group.by = "orig.ident")
```

#### Step 6: Guide Expression Visualization (lines 149-176)

**Lines 149-155**: Get list of perturbed genes
```r
DefaultAssay(feat_c_singlet) <- "sgRNA"
perturbed_gene_list <- feat_c_singlet@assays[["sgRNA"]]@counts@Dimnames[[1]]
perturbed_gene_list <- perturbed_gene_list[-4]  # Remove one guide
```

**Lines 159-163**: Plot all guide distributions
```r
for (i in perturbed_gene_list) {
  print(i)
  FeaturePlot(feat_c_singlet, features = i)
  ggsave(paste0(i,"_sgrnaplot.jpg"), plot = last_plot(), device = "jpg", width = 7, height = 5)
}
```

**Lines 167-176**: Plot target gene expression
```r
DefaultAssay(feat_c_singlet) <- "RNA"
for (i in perturbed_gene_list) {
  FeaturePlot(feat_c_singlet, features = i)
  ggsave(paste0(i,"_featureplot.jpg"), plot = last_plot(), device = "jpg", width = 7, height = 5)
}
```

**Key comparison**: Guide presence vs gene expression
- Shows knockout efficiency
- Validates guide-gene linkage

#### Step 7: Prepare for PS Calculation (lines 189-195)

**Line 185**: Check guide distribution
```r
table(feat_c_singlet@meta.data$gene)
```

**Line 190**: Load barcode frame
```r
bc_frame = read.table(BARCODE, header = T)
```

**Line 194**: Define control
```r
non_target_ctrl = "Non-Targeting"
```

#### Step 8: Calculate PS Scores for All Perturbations (lines 198-204)

**Lines 198-199**: Multi-gene PS calculation
```r
eff_obj <- scmageck_eff_estimate(
  feat_c_singlet, bc_frame,
  perturb_gene = perturbed_gene_list,  # All genes at once
  non_target_ctrl,
  scale_factor = 3  # Adjust for effect size
)
```

**Parameters:**
- `perturb_gene = perturbed_gene_list`: Calculate PS for multiple genes
- `scale_factor = 3`: Amplifies PS scores for visualization
- Higher scale_factor = more sensitive to small effects

**Lines 202-203**: Extract results
```r
eff_estimat = eff_obj$eff_matrix
rds_subset = eff_obj$rds
```

#### Step 9: Visualize PS Scores and Phenotype (lines 207-242)

**Line 208**: Check updated object
```r
DefaultAssay(rds_subset) <- "RNA"
DimPlot(rds_subset, group.by = "orig.ident")
```

**Line 212**: Color by perturbation
```r
DimPlot(rds_subset, group.by = "gene")
```

**Lines 218-227**: Plot all PS scores
```r
for(pb in perturbed_gene_list){
  p = FeaturePlot(rds_subset, features = paste(pb,'eff',sep='_'))
  print(p)
  p = FeaturePlot(rds_subset, features = pb)
  print(p)
}
```

**Line 218**: **Figure 4d** - HIV-GFP expression (phenotype)
```r
FeaturePlot(rds_subset, features = "H13Ld2EGFP", order = TRUE)
```
- `H13Ld2EGFP`: HIV reporter gene
- `order = TRUE`: Plot high values on top
- Shows HIV reactivation pattern

**Lines 233-234**: **Figure 4c** - BRD4 PS score
```r
pb = "BRD4"
FeaturePlot(rds_subset, features = paste(pb,'eff',sep='_'), order = TRUE) + ggtitle("BRD4 PS score")
```

**Lines 240-241**: **Figure 4f** - CCNT1 PS score
```r
pb = "CCNT1"
FeaturePlot(rds_subset, features = paste(pb,'eff',sep='_'), order = TRUE) + ggtitle("CCNT1 PS score")
```

#### Step 10: Identify Marker Genes for Responsive Cluster (lines 246-262)

**Line 247**: Find markers for cluster 8
```r
cluster8.markers <- FindMarkers(feat_c_singlet, ident.1 = 8, min.pct = 0.25)
```

**Why cluster 8?**
- Visual inspection revealed BRD4 PS-high cells cluster here
- Distinct transcriptional state
- Hypothesis: This cluster represents HIV-reactivatable cells

**Lines 251-255**: Export marker genes
```r
head(cluster8.markers, n = 10)
write.table(cluster8.markers, file = 'cluster8_markers.txt', sep = '\t', quote = F, row.names = T)
```

**Lines 259-261**: Alternative marker search (no filter)
```r
cluster8.markers.2 <- FindMarkers(feat_c_singlet, ident.1 = 8)
write.table(cluster8.markers.2, file = 'cluster8_markers.v2.txt', sep = '\t', quote = F, row.names = T)
```

#### Step 11: BRD4 Target Gene Signature Analysis (lines 268-285)

**Line 269**: Subset to BRD4-perturbed cells only
```r
rds_brd4 <- subset(feat_c_singlet, cells = rownames(feat_c_singlet@meta.data)[feat_c_singlet@meta.data$gene == 'BRD4'])
```

**Lines 271-272**: Load BRD4 target genes
```r
brd4_target <- read.table('BRD4_targets.txt', header = T)
brd4_target <- brd4_target[,1]
```

**Lines 273-276**: Calculate signature score
```r
z <- GetAssayData(rds_brd4)
z <- z[rownames(z) %in% brd4_target,]
zmean <- colMeans(z)  # Average expression of BRD4 targets
rds_brd4 <- AddMetaData(rds_brd4, zmean, col.name = 'BRD4_targets')
```

**Lines 277-284**: **Figure S6e** - Compare signature across clusters
```r
z <- rds_brd4@meta.data
z$cluster <- ifelse(z$seurat_clusters == 8, 'Cluster 8', 'Others')
wxs <- wilcox.test(z$BRD4_targets[z$seurat_clusters == 8],
                   z$BRD4_targets[z$seurat_clusters != 8])
ggplot(z, aes(x = cluster, y = BRD4_targets, fill = cluster)) + 
  geom_violin() +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  theme_classic() +
  ggtitle('BRD4 targets', subtitle = paste('p=', wxs$p.value))
```

**Key finding**: Cluster 8 has significantly different BRD4 target expression

#### Step 12: Differential Expression in BRD4 PS+ Cells (lines 289-318)

**Lines 290-292**: Add BRD4 guide presence metadata
```r
gene_e <- GetAssayData(feat_c_singlet, assay = 'sgRNA')
pfr <- ifelse(gene_e['BRD4',] > 0, 1, 0)
feat_c_singlet <- AddMetaData(feat_c_singlet, pfr, col.name = 'BRD4_guides')
```

**Lines 296-300**: **Figure S6f** - DE analysis within cluster 8
```r
rds_cluster8 <- subset(feat_c_singlet, idents = '8')
de_brd4 = FindMarkers(rds_cluster8, ident.1 = 1, ident.2 = 0, group.by = 'BRD4_guides')
```
- Compares cluster 8 cells WITH vs WITHOUT BRD4 guides
- Identifies BRD4-specific effects in responsive cluster

**Lines 305-308**: Alternative - Compare to ALL other cells
```r
pfr2 <- ifelse(pfr == 1 & feat_c_singlet$seurat_clusters == '8', 1, 0)
feat_c_singlet <- AddMetaData(feat_c_singlet, pfr2, col.name = 'cluster8_BRD4_guides')
de_brd4_2 <- FindMarkers(feat_c_singlet, ident.1 = 1, ident.2 = 0, group.by = 'cluster8_BRD4_guides')
```

**Lines 312-317**: Volcano plot with HIV-GFP highlighted
```r
de_brd4_2$colors <- ifelse(rownames(de_brd4_2) == "H13Ld2EGFP", "red", "black")
plot(de_brd4_2$avg_log2FC, -log10(de_brd4_2$p_val_adj),
     pch = 20, xlab = 'log2 Fold Change', ylab = '-log10 (adj p value)',
     col = de_brd4_2$colors) + title(main = "Cells with strong BRD4 perturbation vs. other cells")
text(2.2, 11, labels = "GFP", col = "red")
```

**Key result**: HIV-GFP is top upregulated gene in BRD4 PS-high cells

### Input Requirements

#### 1. Seurat Object (`JKLAT_H13Ld2EGFP.combined.rds`)
- **Type**: Combined/integrated Seurat object
- **Assays**: RNA, sgRNA
- **Status**: Pre-normalized, pre-integrated
- **Contains**: HIV-GFP reporter expression (`H13Ld2EGFP` gene)

#### 2. Barcode File (`BARCODE_H13Ld2EGFP.txt`)
- **Size**: 1.9 MB (large, many cells)
- **Format**: Standard barcode format (cell, barcode, sgrna, gene, read_count, umi_count)

#### 3. BRD4 Target Genes (`BRD4_targets.txt`)
- **Source**: Literature/databases (genes regulated by BRD4)
- **Format**: One gene per line
- **Purpose**: Gene signature scoring

### Output

#### 1. PS Scores for Multiple Perturbations
- BRD4_eff, CCNT1_eff, and others
- Added to Seurat metadata
- Visualized on UMAP

#### 2. Publication Figures
- Figure 4b: UMAP by sample
- Figure 4c: BRD4 PS score distribution
- Figure 4d: HIV-GFP expression
- Figure 4f: CCNT1 PS score distribution
- Figure S6a: QC violin plots
- Figure S6b: Clustering UMAP
- Figure S6d: BRD4 guide distribution
- Figure S6e: BRD4 target signature comparison
- Figure S6f: Volcano plot

#### 3. Marker Gene Tables
- `cluster8_markers.txt`: Cluster 8 defining genes
- `cluster8_markers.v2.txt`: Unfiltered markers
- Differential expression results

### Key Biological Insights

**Discovery 1: PS reveals responsive subpopulation**
- Not all BRD4-knockout cells respond equally
- Cluster 8 enriched for high PS scores
- Represents HIV-reactivatable cell state

**Discovery 2: BRD4 PS correlates with HIV reactivation**
- High BRD4 PS → High HIV-GFP expression
- PS score is better predictor than BRD4 expression alone
- Validates PS method for phenotypic studies

**Discovery 3: Cluster 8 has distinct BRD4 target signature**
- BRD4 target genes expressed differently
- Pre-existing transcriptional state
- Permissive for HIV reactivation

### Advanced Analysis Patterns Demonstrated

1. **Multi-perturbation analysis**: Calculate PS for many genes at once
2. **Signature scoring**: Aggregate expression of gene sets
3. **Cluster-specific DE**: Find markers defining responsive populations
4. **Conditional DE**: Compare perturbed vs control within clusters
5. **Phenotype correlation**: Link PS scores to biological readout (GFP)
6. **Statistical testing**: Wilcoxon tests for group comparisons

### Usage Example - Adapt to Your Screen

```r
# For a different CRISPR screen with phenotype

# 1. Load your data
sobj = readRDS("your_perturbseq.rds")
bc_frame = read.table("your_barcodes.txt", header = T)

# 2. Get list of perturbed genes
gene_list = unique(bc_frame$gene)
gene_list = gene_list[gene_list != "Non-Targeting"]  # Remove control

# 3. Calculate PS for all genes
eff_obj <- scmageck_eff_estimate(
  sobj, bc_frame,
  perturb_gene = gene_list,
  non_target_ctrl = "Non-Targeting",
  scale_factor = 3
)

# 4. Identify phenotype-associated perturbations
# (correlate PS scores with your phenotype of interest)
phenotype = FetchData(eff_obj$rds, vars = "YOUR_PHENOTYPE_GENE")
ps_scores = FetchData(eff_obj$rds, vars = grep("_eff", colnames(eff_obj$rds@meta.data), value = T))

cors = cor(phenotype, ps_scores, method = "spearman")
print(sort(cors[1,], decreasing = TRUE))  # Top hits
```

### Dependencies

#### R Packages
- **Seurat**: Single-cell framework
- **scMAGeCK**: PS calculation
- **ggplot2**: Plotting
- **hdf5r**: Reading H5 files (if needed)

#### External Files
- BRD4 target gene list (optional, for signature analysis)

### Common Issues and Solutions

**Issue 1**: Cluster identification differs
```r
# Solution: Clustering is stochastic; set seed for reproducibility
set.seed(42)
feat_c_singlet <- FindClusters(feat_c_singlet, resolution = 0.5)
```

**Issue 2**: Scale factor too high/low
```r
# PS scores too high (>2): reduce scale_factor
# PS scores too low (<0.2): increase scale_factor
# Typical range: 1-5
```

**Issue 3**: Not enough control cells
```r
# Check control cell count
table(bc_frame$gene)["Non-Targeting"]
# Need at least 50-100 for robust statistics
```

**Issue 4**: Memory issues with JackStraw
```r
# JackStraw is memory-intensive; can skip for large datasets
# Use ElbowPlot instead to choose dimensions
ElbowPlot(sobj)
```

---

## Workflow 4: Pancreatic Differentiation Dataset

**File**: `datasets/Pancreatic_differentiation_scRNA-seq/pancreatic_scrnaseq.Rmd`
**Lines of Code**: 537
**Level**: Expert
**Time to Run**: 2-4 hours
**Data Type**: Lineage tracing + CRISPR knockout in pancreatic differentiation
**Scientific Context**: Cell-type-specific perturbation responses during organ development
**Authors**: Bicna Song, Wei Li (Children's National Hospital)

### Purpose

This workflow demonstrates the **most advanced PS analysis** - calculating cell-type-specific and context-dependent perturbation responses. It showcases:

1. Working with developmental/differentiation datasets
2. Cell-type-specific PS score calculation
3. Multiple response patterns for single perturbation
4. Custom target gene selection per cell type
5. Comparing PS scores across differentiation stages
6. Advanced parameter tuning (background correction, lambda, target gene filtering)

### When to Use This Workflow

✅ **Use if:**
- Analyzing perturbations in developmental contexts
- Studying cell-type-specific perturbation responses
- Have heterogeneous cell populations (multiple cell types/states)
- Need to discover context-dependent effects
- Want to see advanced scMAGeCK parameter usage

❌ **Don't use if:**
- Your dataset is homogeneous (one cell type)
- You're new to PS analysis (start with Demo 1)
- You don't need cell-type-specific analysis

### Scientific Context

**Research Question**: How do HHEX, FOXA1, and CCDC6 knockouts affect pancreatic differentiation?

**Key Challenge**: Same perturbation causes different responses in different cell types:
- **HHEX knockout** in DE (definitive endoderm) cells → Pattern 1
- **HHEX knockout** in PP (pancreatic progenitor) cells → Pattern 2
- **CCDC6 knockout** in LV/DUO cells → Pattern 1
- **CCDC6 knockout** in DE cells → Pattern 2

**Solution**: Calculate PS scores using cell-type-specific target genes

### Cell Type Legend

Abbreviations used throughout:
- **DE**: Definitive Endoderm (clusters 0,1,5,13)
- **PP**: Pancreatic Progenitors (cluster 6)
- **PP-transition**: PP in transition (clusters 3,6,8,12)
- **LV/DUO**: Liver/Duodenum (clusters 4,7,10)
- **DE-transition**: DE in transition (cluster 2)
- **WT**: Wild-type control (clone 47-WT)

### Algorithm - Step by Step

#### Part 1: Initial Setup and Global PS Calculation (lines 1-108)

##### Step 1: Load Libraries and Data (lines 18-27)

```r
library(Seurat)
library(scMAGeCK)
library(hdf5r)
library(ggplot2)

rds <- readRDS("10clones_seurat.rds")
```

**Input**: Pre-processed Seurat object with:
- Multiple clones (different CRISPR knockouts)
- Already clustered and annotated
- Multiple differentiation states

##### Step 2: Visualize Starting Data (lines 29-31)

```r
DimPlot(rds, reduction = "umap")
```

##### Step 3: Prepare Barcode File from Metadata (lines 33-48)

**Lines 34-39**: Extract guide info from Seurat metadata
```r
bc_frame <- rds@meta.data[,c("sgrna","gene","nCount_sgRNA")]
bc_frame[,"cell"] <- rownames(bc_frame)
bc_frame[,"barcode"] <- bc_frame[,"sgrna"]
colnames(bc_frame) <- c("sgrna","gene","read_count","cell","barcode")
bc_frame <- bc_frame[,c("cell","barcode","sgrna","gene","read_count")]
bc_frame[,"umi_count"] <- bc_frame[,"read_count"]
```

**Key difference from other workflows**: Guide info already in metadata, not separate file

**Line 44**: Clean up missing values
```r
bc_frame <- bc_frame[!is.na(bc_frame$gene),]
```

**Line 46**: Fix gene name format
```r
bc_frame$gene <- sub("_","-",bc_frame$gene)  # Convert underscores to hyphens
```

##### Step 4: Preprocess RDS (lines 54-60)

**Line 55**: Add guide metadata
```r
rds <- pre_processRDS(bc_frame, rds)
```

**Line 59**: Check clone distribution
```r
table(rds$gene)
```

##### Step 5: Define Target Genes and Control (lines 62-68)

**Lines 63-65**: List of perturbed clones
```r
targetgenelist <- c("43-HHEX", "44-HHEX", "45-HHEX", "51-FOXA1", 
                    "53-FOXA1/2", "57-OTUD5", "58-OTUD5", 
                    "61-CCDC6", "62-CCDC6")
targetgenelist_geneid <- c("HHEX","HHEX","HHEX","FOXA1","FOXA2",
                           "OTUD5","OTUD5","CCDC6","CCDC6")
```

**Line 67**: Negative control
```r
negative_ctrl_gene <- "47-WT"  # Wild-type control clone
```

##### Step 6: Initial PS Calculation (lines 71-79)

**Lines 71-73**: Calculate PS with automatic target gene discovery
```r
eff_obj <- scmageck_eff_estimate(
  rds, bc_frame, targetgenelist, negative_ctrl_gene,
  scale_factor = 6,           # High amplification for developmental effects
  assay_for_cor = "RNA",      # Use RNA assay for correlation
  perturb_gene_exp_id_list = targetgenelist_geneid  # Gene expression IDs
)
```

**Parameters explained:**
- `scale_factor = 6`: Higher than previous workflows (developmental effects subtle)
- `assay_for_cor = "RNA"`: Specifies which assay to use
- `perturb_gene_exp_id_list`: Maps clone IDs to gene expression IDs

##### Step 7: Visualize Initial Results (lines 82-92)

**Lines 82-91**: Plot PS scores and gene expression for all clones
```r
for(pb_i in 1:length(targetgenelist)){
  pb <- targetgenelist[pb_i]
  p <- FeaturePlot(rds_subset, features = paste(pb, "eff", sep = "_"))
  print(p)
  pb_gene <- targetgenelist_geneid[pb_i]
  p <- FeaturePlot(rds_subset, features = pb_gene)
  print(p)
}
```

#### Part 2: CCDC6 Analysis - Multiple Response Patterns (lines 110-269)

**Key Discovery**: CCDC6 knockout shows different effects in different cell types

##### Step 8: Identify Cell-Type-Specific Differentially Expressed Genes (lines 113-150)

**Lines 114-115**: DE analysis in LV/DUO cells
```r
mk1_47 <- FindMarkers(rds_c47, group.by = "gene", ident.1 = "61-CCDC6", 
                      ident.2 = "47-WT", logfc.threshold = 0.1)
```

**Lines 117-118**: DE analysis in DE cells
```r
mk1_01513 <- FindMarkers(rds_c01513, group.by = "gene", ident.1 = "61-CCDC6",
                         ident.2 = "47-WT", logfc.threshold = 0.1)
```

**Lines 120-121**: DE analysis in PP cells
```r
mk1_6_PP <- FindMarkers(rds_c6_PP, group.by = "gene", ident.1 = "61-CCDC6",
                        ident.2 = "47-WT", logfc.threshold = 0.1)
```

**Lines 123-124**: DE analysis in PP-transition cells
```r
mk1_63812_PP <- FindMarkers(rds_c63812_PP, group.by = "gene", ident.1 = "61-CCDC6",
                            ident.2 = "47-WT", logfc.threshold = 0.1)
```

**Why separate DE analyses?**
- Different cell types have different baseline expression
- CCDC6 affects different genes in different contexts
- Need cell-type-specific target genes for PS calculation

##### Step 9: Re-calculate CCDC6 PS with Global Target Genes (lines 154-189)

**Lines 162-170**: Calculate PS using all differentially expressed genes
```r
eff_obj2 <- scmageck_eff_estimate(
  rds_subset2, bc_frame, targetgenelist, negative_ctrl_gene,
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = targetgenelist_geneid,
  background_correction = T  # NEW: Enable background correction
)
```

**New parameter:**
- `background_correction = T`: Removes correlations present in control cells

##### Step 10: Cell-Type-Specific PS - PP/PP-transition (lines 192-209)

**Lines 194-196**: Select target genes specific to PP/PP-transition
```r
target_gene_63812 <- rownames(mk1_63812_PP)[
  abs(mk1_63812_PP$avg_log2FC) > 0.25 & mk1_63812_PP$p_val_adj < 0.05]
```

**Filtering criteria:**
- `|log2FC| > 0.25`: Significant fold change
- `p_val_adj < 0.05`: Statistically significant

**Lines 196-201**: Calculate PS using PP-specific target genes
```r
eff_obj5 <- scmageck_eff_estimate(
  rds_subset2, bc_frame, targetgenelist, negative_ctrl_gene,
  perturb_target_gene = target_gene_63812,  # Cell-type-specific targets!
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = targetgenelist_geneid,
  background_correction = T
)
```

**NEW parameter:**
- `perturb_target_gene`: Manually specify which genes to use for PS calculation

**Lines 207-208**: **Figure 5f** - PP/PP-transition pattern
```r
FeaturePlot(eff_obj5$rds, features = grep("eff", colnames(eff_obj5$rds@meta.data), value = T)) + 
  ggtitle("Pattern 1: PP/PP in transition")
```

##### Step 11: Cell-Type-Specific PS - DE cells (lines 212-229)

**Lines 214-221**: Calculate PS using DE-specific target genes
```r
target_gene_01513 <- rownames(mk1_01513)[
  abs(mk1_01513$avg_log2FC) > 0.25 & mk1_01513$p_val_adj < 0.05]

eff_obj3 <- scmageck_eff_estimate(
  rds_subset2, bc_frame, targetgenelist, negative_ctrl_gene,
  perturb_target_gene = target_gene_01513,  # DE-specific targets
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = targetgenelist_geneid,
  background_correction = T
)
```

**Lines 227-228**: **Figure 5f** - DE pattern
```r
FeaturePlot(eff_obj3$rds, features = grep("eff", colnames(eff_obj3$rds@meta.data), value = T)) + 
  ggtitle("Pattern 2: DE")
```

**Key insight**: Same perturbation (CCDC6), different patterns in different cell types!

##### Step 12: Cell-Type-Specific PS - LV/DUO cells (lines 232-249)

**Lines 234-241**: LV/DUO-specific PS
```r
target_gene_47 <- rownames(mk1_47)[
  abs(mk1_47$avg_log2FC) > 0.25 & mk1_47$p_val_adj < 0.05]

eff_obj4 <- scmageck_eff_estimate(
  rds_subset2, bc_frame, targetgenelist, negative_ctrl_gene,
  perturb_target_gene = target_gene_47,  # LV/DUO-specific
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = targetgenelist_geneid,
  background_correction = T
)
```

**Line 247**: **Figure S10a** - LV/DUO pattern
```r
FeaturePlot(eff_obj4$rds, features = grep("eff", colnames(eff_obj4$rds@meta.data), value = T)) + 
  ggtitle("Pattern 1: CCD6 PS score from LV/DUO")
```

##### Step 13: Cell-Type-Specific PS - DE-transition (lines 252-269)

**Lines 254-261**: DE-transition-specific PS
```r
target_gene_c2 <- rownames(mk1_2_DE)[
  abs(mk1_2_DE$avg_log2FC) > 0.25 & mk1_2_DE$p_val_adj < 0.05]

eff_obj6 <- scmageck_eff_estimate(
  rds_subset2, bc_frame, targetgenelist, negative_ctrl_gene,
  perturb_target_gene = target_gene_c2,  # DE-transition-specific
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = targetgenelist_geneid,
  background_correction = T
)
```

**Line 267**: **Figure S10b** - DE-transition pattern
```r
FeaturePlot(eff_obj6$rds, features = grep("eff", colnames(eff_obj6$rds@meta.data), value = T)) + 
  ggtitle("Pattern 2: CCD6 PS score from DE in transition")
```


#### Part 3: HHEX Analysis - Multiple Clones, Multiple Patterns (lines 271-484)

**Scientific Context**: Three HHEX knockout clones (43, 44, 45) analyzed

##### Step 14: Initial HHEX DE Analysis (lines 313-335)

**Lines 314-317**: DE in different cell types
```r
mkhhex_1 <- FindMarkers(rds_subset3_rerun, group.by = "gene", ident.1 = "45-HHEX",
                        ident.2 = "47-WT", subset.ident = c("0","1","5","13"), 
                        logfc.threshold = 0.1)  # DE cells

mkhhex_2 <- FindMarkers(rds_subset3_rerun, group.by = "gene", ident.1 = "45-HHEX",
                        ident.2 = "47-WT", subset.ident = c("4","7","10"), 
                        logfc.threshold = 0.1)  # LV/DUO cells

mkhhex_3 <- FindMarkers(rds_subset3_rerun, group.by = "gene", ident.1 = "45-HHEX",
                        ident.2 = "47-WT", subset.ident = c("6"), 
                        logfc.threshold = 0.1)  # PP cells

mkhhex_36812 <- FindMarkers(rds_subset3_rerun, group.by = "gene", ident.1 = "45-HHEX",
                            ident.2 = "47-WT", subset.ident = c("3","6","8","12"), 
                            logfc.threshold = 0.1)  # PP/PP-transition
```

**Lines 318-319**: DE using ALL three HHEX clones combined
```r
mkhhex_3_all <- FindMarkers(rds_subset3, group.by = "gene", 
                            ident.1 = c("43-HHEX", "44-HHEX", "45-HHEX"),
                            ident.2 = "47-WT", subset.ident = c("6"), 
                            logfc.threshold = 0.1)

mkhhex_36812_all <- FindMarkers(rds_subset3, group.by = "gene",
                                ident.1 = c("43-HHEX", "44-HHEX", "45-HHEX"),
                                ident.2 = "47-WT", subset.ident = c("3","6","8","12"),
                                logfc.threshold = 0.1)
```

**Strategy**: Combine multiple clones for increased power

##### Step 15: Define Cell-Type-Specific Target Genes (lines 339-344)

**Lines 339-343**: Select targets with stringent thresholds
```r
hhex_cluster01513_target_gene <- mkhhex_1[
  mkhhex_1$p_val_adj < 0.05 & abs(mkhhex_1$avg_log2FC) > 0.5, 1]

hhex_cluster4710_target_gene <- mkhhex_2[
  mkhhex_2$p_val_adj < 0.25, 1]  # Relaxed threshold (very few diff exp genes)

hhex_cluster6_target_gene <- mkhhex_3[
  mkhhex_3$p_val_adj < 0.05 & abs(mkhhex_3$avg_log2FC) > 0.75, 1]

hhex_cluster36812_target_gene <- mkhhex_36812[
  mkhhex_36812$p_val_adj < 0.05 & abs(mkhhex_36812$avg_log2FC) > 0.7, 1]
```

**Threshold tuning**: Different cell types need different stringency
- LV/DUO: Relaxed (p < 0.25) - few DE genes
- DE: Moderate (p < 0.05, |FC| > 0.5)
- PP: Stringent (p < 0.05, |FC| > 0.75)

##### Step 16: HHEX PS - DE Cells (lines 356-383)

**Lines 361-367**: Calculate PS for DE cells
```r
eff_obj8 <- scmageck_eff_estimate(
  rds_subset4, bc_frame, targetgenelist, negative_ctrl_gene,
  perturb_target_gene = hhex_cluster01513_target_gene,  # DE-specific
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = targetgenelist_geneid,
  lambda = 0.0,              # Disable regularization
  background_correction = T
)
```

**New parameter:**
- `lambda = 0.0`: Disable L2 regularization penalty
- When to use: Clean data, strong effects, want maximum sensitivity

**Lines 371-382**: Visualize all three HHEX clones
```r
for(pb_i in 1:length(targetgenelist)){
  pb <- targetgenelist[pb_i]
  p <- FeaturePlot(rds_subset4_rerun, features = paste(pb, "eff", sep="_"))
  print(p)
  pb_gene <- targetgenelist_geneid[pb_i]
  p <- FeaturePlot(rds_subset4_rerun, features = pb_gene)
  print(p)
}
```

##### Step 17: HHEX PS - PP/PP-transition (lines 387-410)

**Lines 388-394**: PP/PP-transition PS
```r
eff_obj9 <- scmageck_eff_estimate(
  rds_subset4, bc_frame, targetgenelist, negative_ctrl_gene,
  perturb_target_gene = hhex_cluster36812_target_gene,  # PP-specific
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = targetgenelist_geneid,
  lambda = 0.0,
  background_correction = T
)
```

##### Step 18: HHEX PS - LV/DUO Cells (lines 414-437)

**Lines 415-421**: LV/DUO PS
```r
eff_obj10 <- scmageck_eff_estimate(
  rds_subset4, bc_frame, targetgenelist, negative_ctrl_gene,
  perturb_target_gene = hhex_cluster4710_target_gene,  # LV/DUO-specific
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = targetgenelist_geneid,
  lambda = 0.0,
  background_correction = T
)
```

##### Step 19: Merge HHEX Clones for Combined Analysis (lines 440-472)

**Lines 445-456**: Merge clone labels
```r
z <- rds_subset5$gene
z[z == "43-HHEX"] <- "HHEX"
z[z == "44-HHEX"] <- "HHEX"
z[z == "45-HHEX"] <- "HHEX"
rds_subset5 <- AddMetaData(rds_subset5, z, col.name = "gene")

bc_frame2 <- bc_frame
z <- bc_frame2$gene
z[z == "43-HHEX"] <- "HHEX"
z[z == "44-HHEX"] <- "HHEX"
z[z == "45-HHEX"] <- "HHEX"
bc_frame2$gene <- z
```

**Why merge?**: Increases statistical power by pooling all HHEX knockout cells

##### Step 20: Calculate Combined HHEX PS (lines 461-477)

**Lines 465-471**: PS with merged clones
```r
eff_obj11 <- scmageck_eff_estimate(
  rds_subset5, bc_frame2, 
  targetgenelist = c("HHEX"),                    # Single merged label
  negative_ctrl_gene,
  perturb_target_gene = hhex_cluster36812_target_gene,
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = c("HHEX"),          # Single gene ID
  background_correction = T,
  lambda = 0
)
```

**Lines 482-483**: **Figure 5c** - Combined HHEX PS
```r
FeaturePlot(subset(rds_subset5_rerun, gene == "HHEX"), features = "HHEX_eff") + 
  ggtitle("HHEX PS score")
```

#### Part 4: FOXA1 Analysis (lines 486-532)

##### Step 21: FOXA1 PS Calculation (lines 497-511)

**Lines 501-505**: Calculate PS for FOXA1 clones
```r
eff_obj12 <- scmageck_eff_estimate(
  rds_subset6, bc_frame, 
  targetgenelist = c("51-FOXA1", "53-FOXA1/2"),
  negative_ctrl_gene,
  scale_factor = 6,
  assay_for_cor = "RNA",
  perturb_gene_exp_id_list = c("FOXA1", "FOXA1"),  # Both target FOXA1
  background_correction = T
)
```

**Note**: Clone 53 targets both FOXA1 and FOXA2, but mapped to FOXA1 expression

##### Step 22: Visualize FOXA1 Results (lines 515-532)

**Line 516**: **Figure S9a** - FOXA1 PS (clone 51)
```r
FeaturePlot(eff_obj12$rds, features = "51-FOXA1_eff") + 
  ggtitle("FOXA1 PS score (clone 51)")
```

**Line 522**: **Figure S9b** - FOXA1 PS (clone 53)
```r
FeaturePlot(eff_obj12$rds, features = "53-FOXA1/2_eff") + 
  ggtitle("FOXA1 PS score (clone 53)")
```

**Line 530**: **Figure S9c** - FOXA1 expression
```r
FeaturePlot(rds_subset6_rerun, features = "FOXA1") + 
  ggtitle("FOXA1 expression (all clones)")
```

### Input Requirements

#### 1. Seurat Object (`10clones_seurat.rds`)

**Characteristics:**
- **10 clones**: 9 CRISPR knockouts + 1 wild-type control
- **Multiple cell types**: DE, PP, LV/DUO, transition states
- **Pre-clustered**: Seurat clusters already assigned
- **Metadata includes**: Clone identity, guide counts, cluster assignments

**Unique feature**: Guide assignments already in metadata (not separate file)

#### 2. Clone Identity Mapping

| Clone ID | Target Gene | Number of Clones |
|----------|-------------|------------------|
| 43-HHEX | HHEX | 1 of 3 |
| 44-HHEX | HHEX | 2 of 3 |
| 45-HHEX | HHEX | 3 of 3 |
| 51-FOXA1 | FOXA1 | 1 of 2 |
| 53-FOXA1/2 | FOXA1 + FOXA2 | 1 of 2 |
| 57-OTUD5 | OTUD5 | 1 of 2 |
| 58-OTUD5 | OTUD5 | 2 of 2 |
| 61-CCDC6 | CCDC6 | 1 of 2 |
| 62-CCDC6 | CCDC6 | 2 of 2 |
| 47-WT | None (control) | 1 |

### Output

#### 1. Publication Figures

**Main Figures:**
- **Figure 5c**: HHEX PS score (combined clones)
- **Figure 5f**: CCDC6 PS patterns (PP vs DE)

**Supplementary Figures:**
- **Figure S9a**: FOXA1 PS (clone 51)
- **Figure S9b**: FOXA1 PS (clone 53)
- **Figure S9c**: FOXA1 expression
- **Figure S10a**: CCDC6 PS in LV/DUO
- **Figure S10b**: CCDC6 PS in DE-transition

#### 2. Multiple PS Score Versions per Gene

For each perturbation, multiple PS scores calculated:
- **Global PS**: Using all differentially expressed genes
- **Cell-type-specific PS**: Using genes DE in specific cell types
- **Different patterns**: Same gene, different cell types

Example for CCDC6:
- `61-CCDC6_eff` (global)
- `61-CCDC6_eff_PP` (PP-specific targets)
- `61-CCDC6_eff_DE` (DE-specific targets)
- `61-CCDC6_eff_LV` (LV/DUO-specific targets)

#### 3. Differential Expression Tables

Saved to `table/` directory:
- Cell-type-specific DE results
- Multiple comparisons per perturbation
- Merged comparisons across cell types

### Key Biological Insights

**Discovery 1: Context-dependent perturbation responses**
- Same genetic perturbation → different transcriptional responses in different cell types
- Cannot use single PS calculation for heterogeneous datasets

**Discovery 2: Cell-type-specific target genes improve PS accuracy**
- Custom target gene selection per cell type reveals patterns missed by global analysis
- Different thresholds needed for different cell types (some have subtle effects)

**Discovery 3: Multiple clones can be combined**
- Merging clones with same target increases statistical power
- Useful when individual clones have low cell counts

**Discovery 4: Developmental stage matters**
- Transition states (DE-transition, PP-transition) show intermediate patterns
- Perturbation effects depend on differentiation stage

### Advanced Techniques Demonstrated

1. **Cell-type-specific DE analysis** with `subset.ident`
2. **Custom target gene selection** with `perturb_target_gene`
3. **Background correction** to remove control correlations
4. **Lambda tuning** for regularization control
5. **Scale factor optimization** for subtle developmental effects
6. **Clone merging** for increased statistical power
7. **Threshold adaptation** per cell type
8. **Multi-pattern discovery** from single perturbation

### scMAGeCK Parameters - Complete Reference

Based on usage in this workflow:

```r
scmageck_eff_estimate(
  RDS,                           # Seurat object
  BARCODE,                       # Barcode data frame
  perturb_gene,                  # Gene(s) to calculate PS for
  non_target_ctrl,               # Control label
  
  # Optional - commonly used:
  scale_factor = 6,              # Amplification factor (1-10, default 3)
  assay_for_cor = "RNA",         # Which assay to use
  perturb_gene_exp_id_list,      # Map clone IDs to gene expression IDs
  
  # Advanced:
  perturb_target_gene,           # Custom target genes (override auto-discovery)
  target_gene_min = 50,          # Minimum target genes (default 30)
  target_gene_max = 100,         # Maximum target genes (default 200)
  lambda = 0.0,                  # Regularization penalty (default 0.01)
  background_correction = TRUE,  # Remove control correlations (default FALSE)
  subset_rds = TRUE              # Return only guide-containing cells (default FALSE)
)
```

### Parameter Tuning Guide

**scale_factor**: How much to amplify PS scores
- Default: 3
- Weak effects (developmental): 5-10
- Strong effects (essential genes): 1-3
- This workflow: 6 (subtle developmental effects)

**target_gene_max/min**: Number of target genes to use
- Default: 30-200
- Subtle effects: 50-100 (more genes = better signal)
- Strong effects: 20-50 (fewer genes sufficient)
- This workflow: Varies by cell type

**lambda**: L2 regularization penalty
- Default: 0.01
- Clean data: 0 (no regularization needed)
- Noisy data: 0.01-0.1 (penalize overfitting)
- This workflow: 0 (high-quality data)

**background_correction**: Remove control cell correlations
- Default: FALSE
- Use when: Control cells show structure unrelated to perturbation
- This workflow: TRUE (removes developmental gradients from controls)

### Usage Example - Adapt to Your Dataset

```r
# Scenario: Multi-cell-type perturbation screen

# 1. Identify cell types in your data
DimPlot(sobj, group.by = "seurat_clusters")
# Annotate clusters manually

# 2. For each perturbation, find DE genes per cell type
de_celltype1 <- FindMarkers(sobj, ident.1 = "Perturbation", ident.2 = "Control",
                            subset.ident = clusters_celltype1)
de_celltype2 <- FindMarkers(sobj, ident.1 = "Perturbation", ident.2 = "Control",
                            subset.ident = clusters_celltype2)

# 3. Select target genes with appropriate thresholds
targets_celltype1 <- rownames(de_celltype1)[
  abs(de_celltype1$avg_log2FC) > 0.25 & de_celltype1$p_val_adj < 0.05]
targets_celltype2 <- rownames(de_celltype2)[
  abs(de_celltype2$avg_log2FC) > 0.25 & de_celltype2$p_val_adj < 0.05]

# 4. Calculate cell-type-specific PS
ps_celltype1 <- scmageck_eff_estimate(
  sobj, bc_frame, "Perturbation", "Control",
  perturb_target_gene = targets_celltype1,
  scale_factor = 6,
  background_correction = TRUE
)

ps_celltype2 <- scmageck_eff_estimate(
  sobj, bc_frame, "Perturbation", "Control",
  perturb_target_gene = targets_celltype2,
  scale_factor = 6,
  background_correction = TRUE
)

# 5. Compare patterns
FeaturePlot(ps_celltype1$rds, features = "Perturbation_eff") + ggtitle("Cell Type 1")
FeaturePlot(ps_celltype2$rds, features = "Perturbation_eff") + ggtitle("Cell Type 2")
```

### Dependencies

#### R Packages
- **Seurat**: Single-cell framework
- **scMAGeCK**: PS calculation
- **ggplot2**: Visualization
- **hdf5r**: File I/O
- **dplyr**: Data manipulation

### Common Issues and Solutions

**Issue 1**: Different patterns not visible
```r
# Solution: Use cell-type-specific target genes, not global
# Calculate DE per cell type, then use perturb_target_gene
```

**Issue 2**: PS scores too weak
```r
# Solution 1: Increase scale_factor (try 6-10)
# Solution 2: Relax target gene thresholds (more genes)
# Solution 3: Enable background_correction
```

**Issue 3**: Clone merging causes errors
```r
# Solution: Ensure both metadata and barcode frame updated
# Check: table(rds$gene) and table(bc_frame$gene) should match
```

**Issue 4**: Not enough target genes found
```r
# Check DE results
table(de_results$p_val_adj < 0.05)  # How many significant?

# Solution: Relax thresholds
targets <- rownames(de_results)[
  abs(de_results$avg_log2FC) > 0.1 &  # Lower from 0.25
  de_results$p_val_adj < 0.1]          # Lower from 0.05
```

**Issue 5**: Control cells cluster with perturbed cells
```r
# Solution: Enable background_correction = TRUE
# This removes correlations present in control population
```

### Best Practices for Heterogeneous Datasets

1. **Cluster first**, annotate cell types
2. **Run DE per cell type** separately
3. **Inspect DE results** to choose appropriate thresholds
4. **Calculate cell-type-specific PS** using `perturb_target_gene`
5. **Compare patterns** across cell types
6. **Merge clones** if needed for power
7. **Tune parameters** based on effect size
8. **Validate** with known marker genes

---

## scMAGeCK Functions Reference

This section documents all scMAGeCK functions used across the four workflows. These are the core functions that enable PS analysis.

### Function Inventory

| Function | Type | Purpose | Used In |
|----------|------|---------|---------|
| scmageck_eff_estimate | Core | Calculate PS scores | All workflows |
| guidematrix_to_triplet | Utility | Convert guide matrix to barcode format | Demo 2, HIV, Pancreatic |
| pre_processRDS | Utility | Add guide metadata to Seurat | Demo 2, HIV, Pancreatic |
| featurePlot | Visualization | Plot guide distribution | HIV |

---

### scmageck_eff_estimate()

**Purpose**: Core function to calculate Perturbation-response Score (PS) for single cells

**Category**: Analysis / PS Calculation

**File References**:
- demo/demo1/ps_demo.R:26-27
- demo/demo2/PS_demo2.Rmd:125-126
- datasets/HIV_Perturb-seq/hiv_perturbseq.Rmd:198-199
- datasets/Pancreatic_differentiation_scRNA-seq/pancreatic_scrnaseq.Rmd:71-73 (and many more)

#### Algorithm

**What it does** (conceptual overview):

1. **Identify cell populations** (lines internal to function):
   - Extract cells with perturbation guides
   - Extract cells with control guides
   - Validate sufficient cells in each group (typically >10 per group)

2. **Differential expression analysis**:
   - Compare perturbation vs control cells
   - Find significantly differentially expressed genes
   - Rank genes by effect size and significance

3. **Target gene selection**:
   - If `perturb_target_gene` provided: use specified genes
   - Otherwise: auto-select top DE genes
   - Filter by `target_gene_min` and `target_gene_max` parameters
   - Default: 30-200 genes

4. **Expression signature calculation**:
   - For each cell, calculate signature score based on target genes
   - Weighted by gene importance (fold change, p-value)
   - Normalize across cells

5. **Background correction** (if `background_correction = TRUE`):
   - Calculate correlation structure in control cells
   - Remove background correlations from perturbation signature
   - Reduces false positives from cell type/state effects

6. **PS score calculation**:
   - Scale signature scores by `scale_factor`
   - Clip negative values to zero (PS ≥ 0)
   - Higher score = stronger perturbation response

7. **Return results**:
   - Add PS scores to Seurat metadata as `{gene}_eff`
   - Return both PS matrix and updated Seurat object

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| **RDS** | Seurat object | required | Seurat object with normalized RNA expression |
| **BARCODE** | data.frame | required | Barcode table (cell, barcode, sgrna, gene, read_count, umi_count) |
| **perturb_gene** | character vector | required | Gene(s) to calculate PS for |
| **non_target_ctrl** | character | required | Label for non-targeting control cells |
| scale_factor | numeric | 3 | Amplification factor for PS scores (1-10) |
| assay_for_cor | character | "RNA" | Which assay to use for correlation |
| perturb_gene_exp_id_list | character vector | NULL | Map perturbation labels to gene expression IDs |
| perturb_target_gene | character vector | NULL | Custom target genes (overrides auto-discovery) |
| target_gene_min | integer | 30 | Minimum number of target genes to use |
| target_gene_max | integer | 200 | Maximum number of target genes to use |
| lambda | numeric | 0.01 | L2 regularization penalty (0 = no regularization) |
| background_correction | logical | FALSE | Remove control cell correlations |
| subset_rds | logical | FALSE | Return only cells with guide assignments |

#### Return Value

**Type**: List with two components

**Structure**:
```r
list(
  eff_matrix = matrix,      # PS score matrix (cells × genes)
  rds = Seurat object       # Updated Seurat object
)
```

**eff_matrix**: Matrix of PS scores
- Rows: Cell barcodes
- Columns: Perturbation genes
- Values: PS scores (typically 0 to 1+)
- Interpretation: Higher = stronger perturbation response

**rds**: Seurat object with added metadata
- New columns: `{gene}_eff` for each perturbation
- Example: If `perturb_gene = "TP53"`, adds column `TP53_eff`
- If `subset_rds = TRUE`: Only returns cells with guides
- If `subset_rds = FALSE`: All cells (non-guide cells have NA)

#### Usage Examples

**Example 1: Basic usage** (from Demo 1)
```r
eff_object <- scmageck_eff_estimate(
  rds_object,
  bc_frame,
  perturb_gene = 'TP53',
  non_target_ctrl = 'NonTargetingControlGuideForHuman'
)

# Extract results
ps_scores = eff_object$eff_matrix
rds_updated = eff_object$rds

# Visualize
FeaturePlot(rds_updated, features = 'TP53_eff')
```

**Example 2: Multiple genes** (from HIV workflow)
```r
eff_obj <- scmageck_eff_estimate(
  sobj,
  bc_frame,
  perturb_gene = c('BRD4', 'CCNT1', 'CDK9'),
  non_target_ctrl = 'Non-Targeting',
  scale_factor = 3
)

# Now have: BRD4_eff, CCNT1_eff, CDK9_eff in metadata
```

**Example 3: Custom target genes** (from Pancreatic workflow)
```r
# First, find DE genes
de_results <- FindMarkers(sobj, ident.1 = "Perturbation", ident.2 = "Control")
target_genes <- rownames(de_results)[
  abs(de_results$avg_log2FC) > 0.25 & de_results$p_val_adj < 0.05]

# Then, use custom targets
eff_obj <- scmageck_eff_estimate(
  sobj,
  bc_frame,
  perturb_gene = 'CCDC6',
  non_target_ctrl = '47-WT',
  perturb_target_gene = target_genes,  # Custom targets
  scale_factor = 6,
  background_correction = TRUE
)
```

**Example 4: Advanced parameters** (from Pancreatic workflow)
```r
eff_obj <- scmageck_eff_estimate(
  sobj,
  bc_frame,
  perturb_gene = 'HHEX',
  non_target_ctrl = '47-WT',
  perturb_target_gene = custom_targets,
  scale_factor = 6,
  assay_for_cor = "RNA",
  lambda = 0.0,                      # Disable regularization
  background_correction = TRUE,      # Remove control structure
  target_gene_min = 50,
  target_gene_max = 100
)
```

#### Dependencies

**Required packages:**
- Seurat: For Seurat object handling
- Matrix: For sparse matrix operations (usually loaded by Seurat)

**Calls internally:**
- Seurat differential expression functions
- Statistical tests (Wilcoxon, t-test)
- Correlation calculations

#### Common Issues

**Issue 1**: "Not enough cells in perturbation group"
- **Cause**: Fewer than ~10 cells with the perturbation guide
- **Solution**: Check guide distribution with `table(bc_frame$gene)`
- **Workaround**: Merge similar perturbations or use lower filtering thresholds

**Issue 2**: PS scores all very low (<0.1)
- **Cause**: Weak perturbation effect or scale_factor too low
- **Solution**: Increase `scale_factor` (try 5-10)
- **Alternative**: Relax `target_gene_max` to include more genes

**Issue 3**: PS scores all very high (>5)
- **Cause**: Strong perturbation effect or scale_factor too high
- **Solution**: Decrease `scale_factor` (try 1-2)

**Issue 4**: PS scores don't match biological expectation
- **Cause**: Background structure (cell type, cell cycle) dominates signal
- **Solution**: Enable `background_correction = TRUE`
- **Alternative**: Regress out confounders before PS calculation

**Issue 5**: Gene name mismatch error
- **Cause**: `perturb_gene_exp_id_list` doesn't match actual gene names
- **Solution**: Check gene names with `rownames(sobj)` and ensure exact match
- **Example**: "TP53" in expression but labeled "p53" in barcode file

**Issue 6**: Function runs very slowly
- **Cause**: Large dataset or many perturbations
- **Solution**: Run perturbations separately, or subset to cells of interest
- **Optimization**: Use `subset_rds = TRUE` to reduce output size

---

### guidematrix_to_triplet()

**Purpose**: Convert guide count matrix (guides × cells) to long-format barcode table

**Category**: Utility / Data Transformation

**File References**:
- demo/demo1/ps_demo.R:9 (commented example)
- demo/demo2/PS_demo2.Rmd:93

#### Algorithm

1. **Input validation**:
   - Takes sparse matrix of guide counts (guides × cells)
   - Takes Seurat object for cell barcode reference

2. **Matrix to triplet conversion**:
   - For each non-zero entry in matrix:
     - Extract cell barcode (column name)
     - Extract guide barcode (row name)
     - Extract count value

3. **Calculate statistics**:
   - read_count: Direct count value
   - umi_count: Typically same as read_count for UMI data

4. **Format output**:
   - Create data frame with required columns
   - One row per guide detected per cell

#### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| **guide_matrix** | Matrix | Sparse matrix of guide counts (guides × cells) |
| **seurat_obj** | Seurat | Seurat object (for cell barcode reference) |

#### Return Value

**Type**: data.frame

**Columns**:
- `cell`: Cell barcode
- `barcode`: Guide barcode identifier
- `read_count`: Number of reads
- `umi_count`: Number of UMIs

**Note**: Still needs `sgrna` and `gene` columns to be added manually

#### Usage Example

```r
# Extract guide count matrix from Seurat CRISPR assay
guide_matrix = sobj[['CRISPR']]@counts

# Convert to triplet format
bc_frame = guidematrix_to_triplet(guide_matrix, sobj)

# Add required columns
bc_frame[,'sgrna'] = bc_frame[,'barcode']
bc_frame[,'gene'] = sub('-[0-9]+$', '', bc_frame[,'barcode'])  # Extract gene from barcode

# Now ready for scmageck_eff_estimate()
```

#### When to Use

✅ **Use when**:
- Your guide counts are in matrix format (common with 10X Feature Barcoding)
- You have a CRISPR assay in your Seurat object
- Need to convert from wide to long format

❌ **Don't use when**:
- You already have a barcode table file
- Your data is already in triplet/long format

---

### pre_processRDS()

**Purpose**: Add guide metadata to Seurat object and prepare for PS calculation

**Category**: Utility / Data Preprocessing

**File References**:
- demo/demo2/PS_demo2.Rmd:111
- datasets/HIV_Perturb-seq/hiv_perturbseq.Rmd:38
- datasets/Pancreatic_differentiation_scRNA-seq/pancreatic_scrnaseq.Rmd:55

#### Algorithm

1. **Parse barcode file**:
   - Read barcode table
   - Extract guide assignments per cell

2. **Identify singlets vs multiplets**:
   - Count guides per cell
   - Cells with 1 guide = singlets
   - Cells with 2+ guides = multiplets

3. **Add metadata to Seurat**:
   - `gene`: Target gene for each cell
   - `sgrna`: Guide sequence for each cell  
   - `nFeature_sgRNA_guides`: Number of guides per cell
   - Additional guide-related metrics

4. **Quality metrics**:
   - UMI counts for guides
   - Read counts for guides
   - Guide expression levels

#### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| **BARCODE** | data.frame or path | Barcode table or path to barcode file |
| **RDS** | Seurat object | Seurat object to annotate |

#### Return Value

**Type**: Seurat object

**Added metadata columns**:
- `gene`: Target gene identity (character)
- `sgrna`: Guide RNA identity (character)
- `nFeature_sgRNA_guides`: Number of guides per cell (integer)
- `umi_count_sgRNA`: Total guide UMIs (numeric)
- Additional guide statistics

#### Usage Example

```r
# From barcode file
bc_frame = read.table("BARCODE.txt", header = T)
sobj = pre_processRDS(bc_frame, sobj)

# Check results
table(sobj$gene)  # Distribution of guides
table(sobj$nFeature_sgRNA_guides)  # Singlets vs multiplets

# Filter for singlets
sobj_singlet = subset(sobj, nFeature_sgRNA_guides == 1)
```

#### When to Use

✅ **Use before**:
- Calculating PS scores
- Any guide-based analysis
- Filtering for singlet cells

✅ **Use after**:
- Creating Seurat object
- QC filtering
- Normalization

---

### featurePlot()

**Purpose**: Visualize guide distribution across cells

**Category**: Visualization

**File References**:
- datasets/HIV_Perturb-seq/hiv_perturbseq.Rmd:34

**Note**: This appears to be a custom function from scMAGeCK, not extensively used in the workflows. Most visualization uses Seurat's `FeaturePlot()` instead.

#### Usage Example

```r
featurePlot(RDS = sobj, BARCODE = bc_frame, TYPE = "Dis")
```

---

