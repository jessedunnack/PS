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

