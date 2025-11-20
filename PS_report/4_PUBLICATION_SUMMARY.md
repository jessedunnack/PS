# PS Method - Publication Summary

**Full Title**: Decoding heterogeneous single-cell perturbation responses

**Authors**: Bicna Song, Dingyu Liu, Weiwei Dai, Natalie McMyn, Qingyang Wang, Dapeng Yang, Adam Krejci, Anatoly Vasilyev, Nicole Untermoser, Anke Loregger, Dongyuan Song, Breanna Williams, Bess Rosen, Xiaolong Cheng, Lumen Chao, Hanuman T. Kale, Hao Zhang, Yarui Diao, Tilmann Bürckstümmer, Jenet M. Siliciano, Jingyi Jessica Li, Robert Siliciano, Danwei Huangfu, and Wei Li

**Journal**: Nature Cell Biology
**Year**: 2025
**Volume**: 27
**Issue**: 3
**Pages**: 493-504
**DOI**: 10.1038/s41556-025-01626-9

**Preprint**: bioRxiv 2023.10.30.564796 (posted November 29, 2023)

---

## Abstract (Paraphrased)

Genetic perturbations produce heterogeneous responses in single cells, but existing methods treat cells as uniformly perturbed or unperturbed. This paper introduces the **perturbation-response score (PS)**, which quantifies the strength of perturbation effects in individual cells using constrained quadratic optimization.

**Key Applications**:
1. **Genome-scale CRISPRi screen**: Analyzed 18,595 genes in 586,000+ Jurkat T cells, identifying TCR signaling components and dose-response relationships
2. **HIV latency**: Discovered state-dependent functions of BRD4 and CCNT1 in viral reactivation
3. **Essential genes**: Classified 2,285 genes as "buffered" (requiring strong perturbation) or "sensitive" (showing moderate effects), revealing compensatory mechanisms in protein complexes
4. **Pancreatic differentiation**: Uncovered novel CCDC6 role in directing endoderm toward liver rather than pancreatic lineage

**Performance**: PS outperforms mixscape in detecting partial gene perturbations and predicting protein expression from transcriptomic data.

---

## Main Findings

### Finding 1: PS Quantifies Partial Perturbations

**Validation**: Synthetic data with 50% vs 100% knockdown efficiency
- 50% perturbation → median PS = 0.32-0.34
- 100% perturbation → median PS > 0.8

**Real Data**: CRISPRi with intentional sgRNA mismatches
- Perfect match > 1bp mismatch > 2bp mismatch > control
- Mixscape cannot differentiate these gradations

**Significance**: First method to quantify dosage without experimental titration

### Finding 2: Genome-Scale Screen Reveals Dose-Response Diversity

**Dataset**: 18,595 genes in Jurkat T cells (CRISPRi Perturb-seq)

**Discovery**: Essential genes show two response patterns
1. **Buffered genes** (e.g., proteasome subunits): 
   - Weak PS at moderate perturbation
   - Compensation mechanisms upregulate related genes
   - Drug targets requiring complete inhibition

2. **Sensitive genes**: 
   - Strong PS even at partial perturbation
   - Immediate functional impact
   - Amenable to partial inhibition

**Validation**: Correlated with published T cell CRISPR screen (AUC = 0.73 vs 0.65 for mixscape)

### Finding 3: Cell-State-Dependent HIV Reactivation

**Question**: Why does BRD4/CCNT1 knockout variably reactivate latent HIV?

**Discovery**: 
- BRD4 effects depend on T cell activation state
- CCNT1 shows nonlinear relationship with HIV-GFP in stimulated cells
- PS identifies responsive subpopulation (cluster 8)
- Cluster 8 has distinct BRD4 target gene signature

**Validation**: 
- TNF-α stimulation confirms CCNT1 state-dependency
- Signature genes enriched in responsive cells

**Significance**: Mechanistic understanding of latency reversal heterogeneity

### Finding 4: Context-Dependent CCDC6 Function in Development

**Question**: Does CCDC6 knockout affect pancreatic differentiation?

**Discovery**:
- Different PS patterns in DE (definitive endoderm) vs PP (pancreatic progenitors) vs LV/DUO (liver/duodenum)
- Cell-type-specific transcriptional responses
- CCDC6 knockout increases liver fate (HNF4A+) and decreases pancreatic fate (PDX1+)

**Validation**: Flow cytometry confirms increased HNF4A+, decreased PDX1+ populations

**Significance**: Novel developmental function, reveals limitations of average-effect analysis

### Finding 5: PS Outperforms Mixscape

**Comparison Metrics**:

| Metric | PS | Mixscape |
|--------|-----|----------|
| Expression correlation (CRISPRi low MOI) | 40% of genes | 0% of genes |
| Expression correlation (CRISPRi high MOI) | 41% of genes | <5% of genes |
| AUC (T cell screen validation) | 0.73 | 0.65 |
| ECCITE-seq protein prediction | 19/25 genes | 6/25 genes |
| Partial perturbation detection | ✅ | ❌ |

**Conclusion**: PS provides more biologically meaningful scores

---

## Methods Overview

### Step 1: Target Gene Identification
- Wilcoxon rank-sum test (perturbed vs control)
- Select 30-200 significant DEGs
- Forms "perturbation signature"

### Step 2: Average Effect Estimation
- scMAGeCK-LR ridge regression
- Estimate β scores (average effect per gene)
- Formula: B = (D^T D + λI)^(-1) D^T Y

### Step 3: Cell-Specific PS Optimization
- Constrained quadratic programming
- Objective: minimize ||Y - Y₀ - Ψ×B||² + λΣψ
- Constraints: 0 ≤ ψ ≤ U for perturbed, ψ = 0 for control
- Normalize: PS = ψ/U

### Optimization Properties
- Convex objective (global optimum guaranteed)
- Non-negativity constraint (biological interpretability)
- L₁ regularization (sparsity, prevent overfitting)

---

## Datasets

### Synthetic Data
- 8 simulated datasets (scDesign3)
- Varying perturbation efficiency (50%, 100%)
- Varying number of DEGs (10-500)

### Real Data

1. **Published K562 CROP-seq** (low MOI, 23 genes)
2. **Published K562 CROP-seq** (high MOI, 342 genes)
3. **Published CRISPRi** with sgRNA mismatches
4. **New: Genome-scale Jurkat CRISPRi** (18,595 genes, 586K+ cells)
5. **Published ECCITE-seq** (25 perturbations, protein validation)
6. **New: HIV latency Perturb-seq** (10 genes, 2 conditions, 7-8K cells each)
7. **New: Essential gene CRISPRi** (2,285 genes in K562)
8. **New: Pancreatic differentiation** (10 clones, 2 stages, 20,678 cells)

**Data Availability**: GEO GSE247601

---

## Key Innovations

### Methodological

1. **Constrained optimization framework** for single-cell scores
2. **Non-negativity constraint** ensures biological interpretability
3. **Presence constraint** prevents spurious scoring
4. **Two-step approach** (average then individual) ensures identifiability

### Biological

1. **Dose-response without titration** experiments
2. **Cell-state dependency** discovery framework
3. **Intrinsic vs extrinsic** determinants of response heterogeneity
4. **Compensation mechanism** revelation in protein complexes

### Computational

1. **Scalable** to 500K+ cells
2. **Open-source** R package (scMAGeCK)
3. **Integrates** with Seurat/standard workflows
4. **Validated** across multiple platforms

---

## Limitations Acknowledged

### In Paper

1. **Survival bias**: "PS probably only reflects perturbation responses in a fraction of cells, rather than the full spectrum" when strong perturbations kill cells

2. **Requires sufficient DEGs**: Performance depends on clear transcriptional signature

3. **Suggested combination**: "Combining PS with prediction methods that model uneven perturbed/control cell distributions" for lethal perturbations

### Implied but Not Explicitly Stated

1. No uncertainty quantification (confidence intervals)
2. Assumes linear dose-response relationship
3. No hypothesis testing framework (p-values)
4. Computational cost higher than simple DE tests

---

## Impact and Significance

### Scientific Impact

**Paradigm Shift**: From binary (perturbed/not) to continuous dosage quantification

**New Biology Discovered**:
- Buffered vs sensitive essential genes
- Cell-state-dependent HIV latency mechanisms
- Context-dependent CCDC6 developmental function

**Method Advancement**: Outperforms existing standard (Mixscape)

### Publication Impact

**Journal**: Nature Cell Biology (Impact Factor ~20)
**Acceptance**: Indicates high-quality validation and peer review
**Timing**: 2023 preprint → 2025 publication (thorough review process)

### Field Impact

**Enables New Questions**:
- What causes response heterogeneity?
- Which genes are buffered vs sensitive?
- How does cell state affect perturbation penetrance?

**Practical Applications**:
- Drug target prioritization
- CRISPR efficiency QC
- Responder identification for therapies

---

## Comparison to Related Work

### vs Mixscape (Papalexi et al. 2021)

**Mixscape Approach**: LDA-based binary classification

**PS Advantages**:
- Continuous scores vs binary
- Better expression correlation (40% vs 0-5%)
- Higher AUC (0.73 vs 0.65)
- Detects partial perturbations

**Use Case**: PS supersedes Mixscape for heterogeneity questions

### vs scMAGeCK-LR (Li et al. 2018)

**scMAGeCK-LR**: Average perturbation effects (β scores)

**PS Relationship**: Extends scMAGeCK-LR with cell-level resolution

**Complementarity**: PS uses scMAGeCK-LR for Step 2, then adds cell scores

### vs SCEPTRE (Barry et al. 2024)

**SCEPTRE**: Calibrated association testing for low-MOI data

**Different Goals**:
- SCEPTRE: Discovery (which genes significant?)
- PS: Characterization (how heterogeneous?)

**Complementarity**: Use SCEPTRE for discovery, PS for follow-up

### vs GSFA (Yao et al. 2023)

**GSFA**: Bayesian factor analysis with uncertainty quantification

**Trade-offs**:
- GSFA: Uncertainty, factors, slower
- PS: Deterministic, direct scores, faster

**Use Case**: GSFA when uncertainty critical, PS when speed/interpretability prioritized

---

## Author Contributions (Implied from Authorship)

**Lead Authors**: Bicna Song, Wei Li (method development, paper writing)

**Key Collaborators**:
- **HIV studies**: J.M. Siliciano, R. Siliciano (Johns Hopkins)
- **Pancreatic differentiation**: D. Huangfu (MSKCC)
- **Statistics**: J.J. Li (UCLA)
- **Experimental validation**: Multiple groups

**Multi-Institutional**: Collaboration across 4+ institutions

---

## Reproducibility

### Code Availability
- ✅ GitHub: https://github.com/davidliwei/PS
- ✅ Integrated into scMAGeCK package
- ✅ Well-documented with examples

### Data Availability
- ✅ GEO: GSE247601
- ✅ Processed data for reproduction
- ✅ Raw data for validation

### Documentation
- ✅ Detailed methods in paper
- ✅ Supplementary materials
- ✅ README in code repository
- ✅ Tutorial repository (this repo)

---

## Reception and Adoption

### Publication Timeline
- **Nov 2023**: bioRxiv preprint posted
- **Mar 2025**: Nature Cell Biology publication
- **14 months**: Review and revision period (thorough validation)

### Citations (as of 2025-01-01)
- Preprint citations: Multiple (exact number in publication databases)
- Early adoption by Perturb-seq community

### Community Response
- Inclusion in PerturBase database (2025)
- Compared in recent method benchmarking papers
- Tutorial materials (this repository)

---

## Future Directions (Implied)

### Method Extensions Needed

1. **Uncertainty quantification**: Bayesian variant or bootstrap intervals
2. **Nonlinear models**: Sigmoid dose-responses, thresholds
3. **Epistasis modeling**: Combinatorial perturbation interactions
4. **Temporal dynamics**: PS(t) for differentiation/time-series
5. **Low-MOI calibration**: SCEPTRE-style validation for sparse data

### Biological Applications

1. **Drug screens**: Dose-response heterogeneity in compound libraries
2. **Cancer therapy**: Identify resistant subpopulations
3. **Immunology**: T cell response heterogeneity
4. **Development**: Cell fate decisions and lineage commitment
5. **Synthetic biology**: Circuit characterization at single-cell level

### Computational Improvements

1. **GPU acceleration**: Faster optimization for mega-scale screens
2. **Sparse implementations**: Memory efficiency
3. **Approximation methods**: Stochastic gradient for ultra-large datasets
4. **Integration tools**: Automated workflows, GUI

---

## Key Quotes (Paraphrased Concepts)

"Perturbations produce heterogeneous responses, not uniform effects"

"PS quantifies the strength of cellular response, not just presence/absence"

"Essential genes show buffering - compensation maintains function despite partial perturbation"

"Cell state determines perturbation penetrance - same knockout, different effects"

"PS enables discovery of intrinsic and extrinsic determinants of response heterogeneity"

---

## Bottom Line

**What this paper did**: Introduced first method for continuous, cell-level quantification of perturbation effects

**Why it matters**: Enables studying response heterogeneity systematically

**Key achievement**: Validated across 8+ datasets, outperforms existing methods, published in top-tier journal

**Impact**: Paradigm shift from average effects to single-cell dosage responses

**Recommendation**: Should be standard tool for Perturb-seq heterogeneity studies

---

**Document Version**: 1.0
**Last Updated**: 2025-11-13
**Purpose**: Publication summary for PS method evaluation

