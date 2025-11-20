# PS Method - Comprehensive Knowledge Base

**Repository**: PS (Perturbation-response Score Analysis)
**Purpose**: Complete knowledge base for comparative evaluation of Perturb-seq analysis methods
**Date**: 2025-11-13
**Version**: 1.0

---

## Overview

This folder contains a comprehensive knowledge base about the **PS (Perturbation-response Score)** method for analyzing single-cell CRISPR screen data. It was created to enable systematic comparison and evaluation of different Perturb-seq analysis approaches.

---

## Document Structure

### 📄 Core Documents

#### 1. **METHOD_OVERVIEW.md**
- **What it covers**: High-level summary of PS method
- **Key sections**:
  - Executive summary and core innovation
  - Problem PS solves vs limitations of existing methods
  - How PS works (algorithmic overview)
  - Performance benchmarks
  - Validated biological applications
  - When to use PS vs when not to
- **Audience**: Anyone needing quick understanding of PS
- **Length**: ~25 pages

#### 2. **MATHEMATICAL_FOUNDATION.md**
- **What it covers**: Detailed mathematical formulation
- **Key sections**:
  - Problem formulation and model equations
  - Constrained quadratic optimization
  - Algorithm steps with equations
  - Statistical framework
  - Computational complexity analysis
  - Parameter selection guidance
  - Model assumptions and violations
- **Audience**: Methodologists, statisticians, algorithm developers
- **Length**: ~30 pages

#### 3. **STRENGTHS_WEAKNESSES_COMPARISON.md**
- **What it covers**: Critical analysis and method comparison
- **Key sections**:
  - Detailed strengths (7 major advantages)
  - Detailed weaknesses (7 limitations)
  - Head-to-head comparisons vs Mixscape, SCEPTRE, scMAGeCK-LR, MAST, GSFA, Mixscale
  - When to use PS vs alternatives (decision trees)
  - Complementary vs competitive methods
  - Recommended workflows
  - Future improvements needed
- **Audience**: Users selecting analysis methods, comparative evaluations
- **Length**: ~35 pages

---

## Quick Access Guide

### If you want to know...

**"What is PS?"**
→ Read: `1_METHOD_OVERVIEW.md` - Executive Summary

**"How does PS work mathematically?"**
→ Read: `2_MATHEMATICAL_FOUNDATION.md` - Sections 1-3

**"Should I use PS for my data?"**
→ Read: `3_STRENGTHS_WEAKNESSES_COMPARISON.md` - Decision Tree section

**"How does PS compare to Mixscape?"**
→ Read: `3_STRENGTHS_WEAKNESSES_COMPARISON.md` - "vs Mixscape" section

**"What are PS's limitations?"**
→ Read: `3_STRENGTHS_WEAKNESSES_COMPARISON.md` - "Weaknesses" section (7 major limitations)

**"What biological discoveries were made with PS?"**
→ Read: `1_METHOD_OVERVIEW.md` - "Validated Biological Applications"

**"What parameters do I need to tune?"**
→ Read: `2_MATHEMATICAL_FOUNDATION.md` - "Parameter Selection" section

**"Can PS handle my experimental design?"**
→ Read: `1_METHOD_OVERVIEW.md` - "Supported Experimental Platforms"

---

## Key Findings Summary

### What Makes PS Unique

1. **First and only method** to quantify partial perturbations at single-cell level
2. **Continuous dose-response scores** (0-1 scale) rather than binary classification
3. **Discovers heterogeneity** in perturbation responses (cell-state, cell-type dependencies)
4. **Published in Nature Cell Biology** (2025) with extensive validation

### PS vs Alternatives - Quick Comparison

| Capability | PS | Mixscape | scMAGeCK-LR | SCEPTRE | MAST |
|------------|-----|----------|-------------|---------|------|
| Cell-level scores | ✅ | ❌ | ❌ | ❌ | ❌ |
| Partial perturbations | ✅ | ❌ | ❌ | ❌ | ❌ |
| Dose-response | ✅ | ❌ | ❌ | ❌ | ❌ |
| Heterogeneity | ✅ | ❌ | ❌ | ❌ | ❌ |
| Statistical testing | ❌ | ✅ | ✅ | ✅ | ✅ |
| FDR control | ❌ | ❌ | ❌ | ✅ | ❌ |
| Speed | Medium | Fast | Fast | Medium | Slow |

### When to Use PS

✅ **Use PS when**:
- Studying response heterogeneity
- Need dose-response quantification
- CRISPR efficiency varies
- Cell-state dependencies expected
- Sufficient cells per condition (>50)

❌ **Don't use PS when**:
- Discovery screen needing FDR control (use SCEPTRE)
- Only need average effects (use scMAGeCK-LR)
- Binary classification sufficient (use Mixscape)
- Very low cell counts (<10)
- No transcriptional signature

### Major Strengths

1. ⭐ **Partial perturbation quantification** - Validated on synthetic & real data
2. ⭐ **Cell-state-dependent effects** - Discovered in HIV latency, essential genes
3. ⭐ **Outperforms existing methods** - 40% vs 0-5% expression correlation (vs Mixscape)
4. ⭐ **Extensive validation** - 8+ datasets, multiple platforms, Nature Cell Biology
5. ⭐ **Novel discoveries** - Buffered vs sensitive genes, context-dependent functions

### Major Limitations

1. ⚠️ **Survival bias** - Cannot measure cells killed by perturbation
2. ⚠️ **Requires DEGs** - Fails with weak transcriptional responses
3. ⚠️ **No statistical testing** - No p-values or FDR control
4. ⚠️ **Interpretability** - PS magnitude context-dependent
5. ⚠️ **Computational cost** - Slower than simple DE tests
6. ⚠️ **Assumes linearity** - May fail for highly nonlinear responses
7. ⚠️ **No uncertainty quantification** - No confidence intervals

---

## Comprehensive Method Details

### Algorithm Overview

```
┌─────────────────────────────────────────────────────────────┐
│  INPUT: scRNA-seq data + guide assignments                   │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 1: Identify Target Genes                               │
│  • Wilcoxon test: perturbed vs control                      │
│  • Select 30-200 significant DEGs                            │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 2: Estimate Average Effects (scMAGeCK-LR)             │
│  • Ridge regression: B = (D^T D + λI)^-1 D^T Y              │
│  • β scores = average effect per gene                        │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 3: Optimize Cell-Specific Scores                       │
│  • Constrained quadratic optimization                        │
│  • min Σ[observed - predicted]^2 + λΣψ                      │
│  • Constraints: 0 ≤ ψ ≤ U, ψ=0 for controls                 │
│  • Normalize: PS = ψ/U                                       │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  OUTPUT: PS scores [0,1] for each cell × perturbation       │
└─────────────────────────────────────────────────────────────┘
```

### Mathematical Core

**Expression Model**:
```
Y = Y₀ + Ψ × B + ε
```

**Optimization Objective**:
```
minimize: Σⱼᵢ [yⱼᵢ - yⱼᵢ⁰ - Σₖ ψⱼₖβₖᵢ]² + λΣⱼₖ ψⱼₖ

subject to:
  0 ≤ ψⱼₖ ≤ U  (for perturbed cells)
  ψⱼₖ = 0      (for control cells)
```

**Properties**:
- Convex objective → global optimum
- Non-negativity → biologically interpretable
- Regularized → prevents overfitting

### Performance Metrics

**Synthetic Data**:
- 50% knockdown → PS = 0.32-0.34 (accurate)
- 100% knockdown → PS > 0.8

**vs Mixscape (CRISPRi)**:
- Expression correlation: 40-41% (PS) vs 0-5% (Mixscape)
- AUC (T cell screen): 0.73 (PS) vs 0.65 (Mixscape)

**ECCITE-seq (protein prediction)**:
- PS outperforms: 19/25 genes (76%)

**Scale**: Successfully tested on 586,000 cells, 18,595 perturbations

---

## Biological Applications Demonstrated

### 1. Genome-Scale CRISPR Screens

**Dataset**: 18,595 genes in Jurkat T cells (586K+ cells)

**Discoveries**:
- Identified TCR complex components
- Ranked signaling molecules by essentiality
- Validated against published hits (AUC = 0.73)

### 2. HIV Latency Reactivation

**Question**: Which genes regulate HIV latency reversal?

**Discoveries**:
- BRD4 and CCNT1 show cell-state-dependent effects
- Identified responsive subpopulation (cluster 8)
- Nonlinear relationship between PS and HIV-GFP expression

**Impact**: Mechanistic understanding of latency heterogeneity

### 3. Essential Gene Dosage Responses

**Dataset**: 2,285 essential genes in K562 cells

**Discoveries**:
- **Buffered genes**: Proteasome subunits (require strong perturbation)
- **Sensitive genes**: Immediate effects at moderate perturbation
- Compensation mechanisms: Upregulation of complex members

**Impact**: Predicts drug target vulnerabilities

### 4. Pancreatic Differentiation

**Dataset**: 10 clones across differentiation stages

**Discoveries**:
- CCDC6 has cell-type-specific functions
- Different PS patterns in DE vs PP vs LV/DUO
- Novel role in liver vs pancreatic lineage decision

**Impact**: Context-dependent developmental gene functions

---

## Software & Implementation

### Package Information

**Name**: scMAGeCK (PS module)

**Language**: R

**Installation**:
```r
library(devtools)
install_github('weililab/scMAGeCK')
```

**Dependencies**:
- Seurat (primary interface)
- Matrix (sparse matrices)
- quadprog (optimization)

### Usage Example

```r
library(scMAGeCK)
library(Seurat)

# Load data
rds <- readRDS("seurat_object.rds")
bc_frame <- read.table("barcode_file.txt", header = TRUE)

# Calculate PS scores
eff_obj <- scmageck_eff_estimate(
  RDS = rds,
  BARCODE = bc_frame,
  perturb_gene = c("TP53", "KRAS"),
  non_target_ctrl = "NTC",
  scale_factor = 3
)

# Extract results
ps_scores <- eff_obj$eff_matrix  # PS score matrix
rds_updated <- eff_obj$rds       # Seurat object with PS in metadata

# Visualize
FeaturePlot(rds_updated, features = "TP53_eff")
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| scale_factor | 3 | Amplification factor (1-10) |
| lambda | 0.01 | Regularization strength |
| target_gene_min | 30 | Minimum target genes |
| target_gene_max | 200 | Maximum target genes |
| background_correction | FALSE | Remove control correlations |

---

## Validation & Benchmarking

### Datasets Used

**Synthetic**: 8 simulated datasets (varying DEGs, efficiency)

**Real Perturb-seq**:
1. K562 CROP-seq (low MOI, 23 genes)
2. K562 CROP-seq (high MOI, 342 genes)
3. CRISPRi with sgRNA mismatches (efficiency validation)
4. Jurkat genome-scale CRISPRi (18,595 genes, 586K cells)
5. ECCITE-seq (25 genes, protein validation)
6. HIV latency Perturb-seq (10 genes, 7-8K cells per condition)
7. Essential gene CRISPRi (2,285 genes)
8. Pancreatic differentiation (10 clones, 20K cells)

**External Validation**:
- Published T cell CRISPR screen (385 positive, 1,297 negative hits)
- Flow cytometry (CCDC6 pancreatic effects)
- TNF-α stimulation (CCNT1 validation)

### Performance Summary

✅ **Accuracy**: Correctly quantifies 50% and 100% knockdowns
✅ **Sensitivity**: Detects partial perturbations missed by Mixscape
✅ **Specificity**: Strong correlation with perturbed gene expression
✅ **Scalability**: Tested on 500K+ cells
✅ **Reproducibility**: Published in Nature Cell Biology (peer-reviewed)

---

## Citations & Resources

### Primary Publication

> Song B, Liu D, Dai W, et al. **Decoding heterogeneous single-cell perturbation responses**. *Nature Cell Biology*. 2025;27(3):493-504. doi:10.1038/s41556-025-01626-9

### Preprint

> Song B, Liu D, Dai W, et al. **Decoding Heterogenous Single-cell Perturbation Responses**. *bioRxiv*. 2023. doi:10.1101/2023.10.30.564796

### Code & Data

- **GitHub**: https://github.com/davidliwei/PS
- **scMAGeCK**: https://github.com/weililab/scMAGeCK
- **GEO**: GSE247601

### Related Methods

- **Mixscape**: Papalexi et al., Genome Biology, 2021
- **SCEPTRE**: Barry et al., Genome Biology, 2024
- **scMAGeCK-LR**: Li et al., Genome Biology, 2018
- **GSFA**: Yao et al., Nature Methods, 2023
- **MAST**: Finak et al., Genome Biology, 2015

---

## Use This Knowledge Base To...

### For Method Selection

1. Read overview of PS capabilities
2. Check decision tree for your use case
3. Review limitations that might affect you
4. Compare to alternatives for your needs

### For Implementation

1. Review algorithm details
2. Check parameter guidance
3. Follow usage examples
4. Consult troubleshooting

### For Comparative Evaluation

1. Use standardized metrics (AUC, correlation, accuracy)
2. Compare across same datasets
3. Consider complementarity (not just competition)
4. Match method to research question

### For Paper Writing

1. Cite primary publication (Nature Cell Biology, 2025)
2. Reference specific capabilities demonstrated
3. Acknowledge limitations
4. Compare to relevant alternatives

---

## Key Takeaways

### What PS Adds to the Field

**Before PS**:
- Binary perturbation classifications
- Average effect estimates only
- Heterogeneity not quantified
- Dose-response required titration experiments

**After PS**:
- Continuous cell-level scores
- Individual response heterogeneity
- Dose-response from single experiment
- Cell-state and cell-type specificity

### Research Enabled by PS

1. **Mechanistic understanding** of perturbation heterogeneity
2. **Drug target prioritization** (buffered vs sensitive genes)
3. **CRISPR QC** (efficiency validation)
4. **Subpopulation discovery** (responder identification)
5. **Context-dependent functions** (cell-type-specific roles)

### Impact on Perturb-seq Field

- **Paradigm shift**: From average effects to single-cell heterogeneity
- **New questions**: Why do cells differ? What determines response?
- **Integration**: Combines with existing methods (not replacement)
- **Validation**: Raises bar for method development (Nature Cell Biology)

---

## Document Maintenance

**Version**: 1.0
**Last Updated**: 2025-11-13
**Created By**: Comprehensive documentation system for PS tutorial repository
**Purpose**: Enable systematic comparative evaluation of Perturb-seq analysis methods

**Recommended Updates**:
- When new PS applications published
- When compared to emerging methods
- When limitations addressed in updates
- When new biological insights discovered

---

## Contact & Support

**Package Maintainer**: Wei Li Lab (Children's National Hospital)
**Code Issues**: https://github.com/davidliwei/PS/issues
**scMAGeCK Issues**: https://github.com/weililab/scMAGeCK/issues

**For This Knowledge Base**: Created as part of PS tutorial repository documentation

---

**END OF INDEX**

This knowledge base provides comprehensive information for evaluating PS relative to other Perturb-seq analysis methods. All documents are designed for systematic comparison and method selection.
