# PS (Perturbation-response Score) Method - Overview

**Method Name**: Perturbation-response Score (PS)
**Package**: scMAGeCK-PS
**Publication**: Song et al., Nature Cell Biology, 2025; 27(3):493-504
**DOI**: 10.1038/s41556-025-01626-9
**Preprint**: bioRxiv 2023.10.30.564796
**Code**: https://github.com/davidliwei/PS (integrated into scMAGeCK)
**Data**: GEO GSE247601

---

## Executive Summary

PS (Perturbation-response Score) is a computational method for **quantifying heterogeneous cellular responses to genetic perturbations at single-cell resolution**. Unlike existing methods that provide binary classifications (perturbed vs unperturbed), PS quantifies the **dose-dependent strength** of perturbation effects in individual cells on a continuous scale (0-1).

### The Core Innovation

PS uses **constrained quadratic optimization** to estimate how strongly each cell responds to a perturbation, enabling discovery of:
1. **Partial perturbations** - Cells with incomplete gene knockdown
2. **Dose-response relationships** - Buffered vs sensitive genes
3. **Cell-state dependencies** - Perturbation effects that vary by cell context
4. **Response heterogeneity** - Subpopulations with different sensitivities

---

## What Problem Does PS Solve?

### The Heterogeneity Challenge in Perturb-seq

Traditional Perturb-seq analysis assumes:
- All cells with a guide are "fully perturbed"
- All perturbation effects are uniform across cells
- Binary classification: perturbed vs control

**Reality**:
- CRISPR efficiency varies widely (20-100% knockdown)
- Cell state affects perturbation penetrance
- Compensation mechanisms create dose-responses
- Subpopulations show different sensitivities

### Limitations of Existing Methods

**Mixscape** (Seurat):
- Binary classification: perturbed vs non-perturbed
- Cannot quantify partial effects
- Assigns posterior probability = 1.0 to all perturbed cells
- Misses dose-response relationships

**scMAGeCK-LR** (Linear Regression):
- Estimates average perturbation effects
- Cannot capture cell-to-cell heterogeneity
- No individual cell scores

**SCEPTRE**:
- Association testing for low-MOI data
- Focuses on calibration and power
- Does not quantify individual cell responses

**MAST**:
- General differential expression tool
- Not designed for perturbation heterogeneity
- No partial effect quantification

---

## How PS Works (High-Level)

### Three-Step Algorithm

```
┌─────────────────────────────────────────────────────────────┐
│  Step 1: IDENTIFY TARGET GENES                               │
│  • Find differentially expressed genes (DEGs)                │
│  • Wilcoxon rank-sum test: perturbed vs control             │
│  • These genes form the "perturbation signature"             │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  Step 2: ESTIMATE AVERAGE EFFECTS                            │
│  • Use scMAGeCK-LR linear regression                         │
│  • Calculate β scores (average effect per gene)              │
│  • B = (D^T D + λI)^-1 D^T Y                                │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  Step 3: OPTIMIZE CELL-SPECIFIC RESPONSES                    │
│  • Solve constrained quadratic optimization                  │
│  • min Σ [observed - predicted]^2 + λΣψ                     │
│  • Constraints: 0 ≤ ψ ≤ U (for perturbed cells)             │
│  •              ψ = 0 (for control cells)                    │
│  • Normalize: PS = ψ/U to range [0,1]                       │
└─────────────────────────────────────────────────────────────┘
```

### Mathematical Foundation

**Core Equations**:

1. **Expression Model**:
   ```
   Y = Y₀ + Ψ × B + ε
   ```
   - Y = observed expression (log-transformed)
   - Y₀ = baseline expression
   - Ψ = PS score matrix (cells × perturbations)
   - B = effect coefficients (perturbations × genes)
   - ε = Gaussian noise

2. **Optimization Objective**:
   ```
   minimize: Σⱼᵢ [yⱼᵢ - yⱼᵢ⁰ - Σₖ ψⱼₖβₖᵢ]² + λΣⱼₖ ψⱼₖ
   ```

3. **Constraints**:
   - **Non-negativity**: ψⱼₖ ≥ 0 (perturbation can't reduce effects)
   - **Presence**: ψⱼₖ = 0 if cell j doesn't have perturbation k
   - **Upper bound**: ψⱼₖ ≤ U (normalization factor)

4. **Final PS Score**:
   ```
   PSⱼₖ = ψⱼₖ / U
   ```
   - Range: [0, 1]
   - 0 = no detectable perturbation effect
   - 1 = maximum perturbation effect

---

## Key Capabilities

### 1. Partial Perturbation Quantification

**Validation** (Synthetic Data):
- 50% knockdown → median PS = 0.32-0.34
- 100% knockdown → median PS > 0.8

**Application**:
- Identifies cells with incomplete CRISPR editing
- Reveals compensation mechanisms
- Enables dose-response analysis without manual titration

### 2. Cell-State-Dependent Effects

**Example** (HIV Latency):
- BRD4 knockout effect varies by T cell activation state
- CCNT1 shows nonlinear response to HIV-GFP in stimulated cells
- PS identifies responsive subpopulation (cluster 8)

### 3. Essential Gene Dosage Responses

**Discovery** (2,285 Essential Genes):
- **Buffered genes**: Require strong perturbation for effect (e.g., proteasome subunits)
- **Sensitive genes**: Show effects at moderate perturbation
- Reveals compensatory upregulation in protein complexes

### 4. Cell-Type-Specific Functions

**Example** (Pancreatic Differentiation):
- CCDC6 knockout shows different responses in:
  - Definitive endoderm (DE)
  - Pancreatic progenitors (PP)
  - Liver/duodenum (LV/DUO)
- Different target genes per cell type
- Context-dependent developmental roles

---

## Performance Benchmarks

### vs Mixscape (CRISPRi Data)

| Metric | PS | Mixscape |
|--------|-----|----------|
| Genes with expression correlation | 40-41% | 0-5% |
| AUC (T cell screen) | 0.73 | 0.65 |
| Partial effect detection | ✅ Yes | ❌ No |
| Continuous scoring | ✅ Yes | ❌ Binary only |

### ECCITE-seq (Protein Prediction)

- PS outperforms mixscape: 19/25 genes (76%)
- Superior on moderate-effect genes
- AUC > 0.8 for 21/25 genes with strong transcriptomic changes

### Genome-Scale Validation

**T Cell Activation Screen**:
- 18,595 genes tested
- 586,000+ cells analyzed
- Successfully identified TCR complex components and signaling molecules
- Correlated with published positive hits (385) vs negative (1,297)

---

## Validated Biological Applications

### 1. Genome-Scale Screens

**Platform**: CRISPRi/a Perturb-seq
**Scale**: 18,000+ genes
**Capability**: Identify essential genes, signaling pathways, dose-responses

### 2. Phenotype Association

**Example**: HIV reactivation, T cell activation
**Capability**: Link PS scores to functional readouts (GFP, cell surface markers, proteomics)

### 3. Developmental Biology

**Example**: Pancreatic differentiation, lineage decisions
**Capability**: Cell-type-specific PS scores reveal context-dependent gene functions

### 4. Drug Target Identification

**Potential**: Identify genes with buffered vs sensitive responses
**Application**: Prioritize targets requiring complete inhibition vs partial modulation

---

## Supported Experimental Platforms

✅ **Perturb-seq** (CRISPR + scRNA-seq)
✅ **CROP-seq** (Combinatorial perturbations)
✅ **ECCITE-seq** (Transcriptome + proteomics)
✅ **sci-Plex** (Multiplexed chemical perturbations)
✅ **General multiplexed scRNA-seq** (Any perturbation + single-cell readout)

**Requirements**:
- Single-cell RNA-seq data
- Perturbation identity per cell (guide barcode or condition label)
- Control cells for comparison
- Sufficient cells per perturbation (typically >10, ideally >50)

---

## Integration with Existing Workflows

PS is **orthogonal** to existing scRNA-seq analysis pipelines:

1. **Seurat/Scanpy**: Can use Seurat/Scanpy objects as input
2. **Confounding correction**: Compatible with methods that remove batch effects, cell cycle, etc.
3. **scMAGeCK-LR**: PS extends scMAGeCK-LR with cell-level scores
4. **Mixscape**: Can be used alongside mixscape for validation

**Typical Workflow Integration**:
```
Raw Data → QC → Normalization → Seurat/Scanpy → PS Calculation → Downstream Analysis
                                       ↓
                              (optional: confounding correction)
```

---

## Software Implementation

### Package: scMAGeCK

**Language**: R
**Dependencies**: Seurat (primary interface)
**Installation**:
```r
library(devtools)
install_github('weililab/scMAGeCK')
```

**Core Function**:
```r
scmageck_eff_estimate(
  RDS,                    # Seurat object
  BARCODE,               # Guide assignments
  perturb_gene,          # Gene(s) to analyze
  non_target_ctrl,       # Control label
  scale_factor = 3,      # Amplification factor
  ...                    # Advanced parameters
)
```

**Output**:
- PS score matrix (cells × perturbations)
- Updated Seurat object with PS scores as metadata
- Target gene lists
- Effect coefficients (β scores)

---

## When to Use PS

### Ideal Use Cases

✅ **When you need to**:
- Quantify dose-response relationships
- Identify partial perturbations
- Discover cell-state dependencies
- Find responsive subpopulations
- Validate CRISPR efficiency
- Correlate PS with phenotypes

✅ **Dataset characteristics**:
- Heterogeneous perturbation efficiency
- Sufficient cells per condition (>50 recommended)
- Clear transcriptional response to perturbation
- Need for continuous scoring (not just binary)

### When NOT to Use PS

❌ **Inappropriate scenarios**:
- Very weak perturbations (few DEGs)
- Insufficient cell counts (<10 per condition)
- Survival bias dominates (lethal perturbations kill all strongly affected cells)
- Only need average effect estimates (use scMAGeCK-LR)
- Need statistical testing framework (use SCEPTRE)

---

## Key Strengths

1. **First method for partial perturbation quantification**
2. **Dose-response discovery without titration experiments**
3. **Cell-state and cell-type heterogeneity analysis**
4. **Validated on multiple platforms and scales**
5. **Published in Nature Cell Biology (high-impact)**
6. **Open-source and well-documented**
7. **Integrates with standard single-cell workflows**

---

## Key Limitations

1. **Survival bias**: Cannot assess effects in dead cells (lethal perturbations)
2. **Power constraints**: May only capture partial spectrum in essential gene contexts
3. **Requires signature genes**: Performance depends on DEG quality and quantity
4. **Computational cost**: Optimization step slower than simple differential expression
5. **Interpretability**: PS magnitude meaning depends on context (not absolute measure)

---

## Quick Comparison Table

| Feature | PS | Mixscape | scMAGeCK-LR | SCEPTRE |
|---------|-----|----------|-------------|---------|
| **Continuous scoring** | ✅ | ❌ | ❌ | ❌ |
| **Partial perturbations** | ✅ | ❌ | ❌ | ❌ |
| **Cell-level scores** | ✅ | ❌ | ❌ | ❌ |
| **Dose-response** | ✅ | ❌ | ❌ | ❌ |
| **Average effects** | ✅ | ✅ | ✅ | ❌ |
| **Statistical testing** | ❌ | ✅ | ✅ | ✅ |
| **Binary classification** | ❌ | ✅ | ❌ | ❌ |
| **Low-MOI calibration** | ❌ | ❌ | ❌ | ✅ |

---

## Citations

**Primary Publication**:
> Song B, Liu D, Dai W, et al. Decoding heterogeneous single-cell perturbation responses. Nature Cell Biology. 2025;27(3):493-504. doi:10.1038/s41556-025-01626-9

**Preprint**:
> Song B, Liu D, Dai W, et al. Decoding Heterogenous Single-cell Perturbation Responses. bioRxiv. 2023. doi:10.1101/2023.10.30.564796

---

## Authors

**Lead Authors**: Bicna Song, Wei Li (Children's National Hospital)

**Contributors**: Dingyu Liu, Weiwei Dai, Natalie McMyn, Qingyang Wang, Dapeng Yang, Adam Krejci, Anatoly Vasilyev, Nicole Untermoser, Anke Loregger, Dongyuan Song, Breanna Williams, Bess Rosen, Xiaolong Cheng, Lumen Chao, Hanuman T. Kale, Hao Zhang, Yarui Diao, Tilmann Bürckstümmer, Jenet M. Siliciano, Jingyi Jessica Li, Robert Siliciano, Danwei Huangfu

---

**Document Version**: 1.0
**Last Updated**: 2025-11-13
**Purpose**: Comprehensive method overview for comparative evaluation of Perturb-seq analysis approaches
