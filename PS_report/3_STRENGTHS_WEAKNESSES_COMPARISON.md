# PS Method - Critical Analysis and Comparison

**Document Purpose**: Comprehensive evaluation of strengths, weaknesses, and comparison to alternative Perturb-seq analysis methods

---

## Strengths of PS Method

### 1. **Unique Capability: Quantifies Partial Perturbations**

**Innovation**: First and only method to provide continuous dosage scores at single-cell level

**Validation**:
- Synthetic data: 50% knockdown → PS ≈ 0.32-0.34
- CRISPRi mismatch: Correctly ranks guide efficiency
- ECCITE-seq: PS correlates with protein knockdown efficiency

**Impact**: Enables dose-response analysis without experimental titration

**Competitors Cannot Do This**:
- Mixscape: Binary classification only
- scMAGeCK-LR: Average effects only
- SCEPTRE: Statistical testing only

### 2. **Discovers Biological Heterogeneity**

**Cell-State Dependencies**:
- HIV latency: BRD4 effects depend on T cell activation state
- Essential genes: Identifies buffered vs sensitive dosage responses
- Subpopulation discovery: Cluster 8 in HIV screen

**Cell-Type Specificities**:
- Pancreatic differentiation: CCDC6 shows different responses in DE vs PP vs LV/DUO
- Enables context-dependent gene function discovery

**Clinical Relevance**: Drug responses likely heterogeneous - PS can identify responder subsets

### 3. **Strong Validation Across Multiple Platforms**

**Platforms Tested**:
- ✅ Perturb-seq (CRISPR + scRNA-seq)
- ✅ CROP-seq (low and high MOI)
- ✅ ECCITE-seq (transcriptome + protein)
- ✅ CRISPRi with sgRNA mismatches
- ✅ Genome-scale screens (18,595 genes, 586K cells)

**Datasets**:
- 8 synthetic datasets
- 8+ real experimental datasets
- Multiple species (human, mouse)
- Multiple cell types (T cells, K562, iPSCs)

**Performance**: Published in Nature Cell Biology (high-impact validation)

### 4. **Integrates with Standard Workflows**

**Compatibility**:
- Uses Seurat objects (most common scRNA-seq format)
- Orthogonal to confounding correction methods
- Can combine with mixscape, scMAGeCK-LR
- Works with any scRNA-seq normalization

**Ease of Use**:
- Single function call: `scmageck_eff_estimate()`
- Automatic parameter selection (reasonable defaults)
- Well-documented R package

### 5. **Mathematical Rigor**

**Optimization**:
- Convex objective → guaranteed global optimum
- Constrained formulation → biologically interpretable
- Regularization → prevents overfitting

**Reproducibility**:
- Deterministic solution (not stochastic)
- Stable across initializations
- Robust to parameter choices within reasonable ranges

### 6. **Outperforms Existing Methods**

**vs Mixscape**:
| Metric | PS | Mixscape | Winner |
|--------|-----|----------|---------|
| Expression correlation | 40-41% | 0-5% | **PS** |
| AUC (T cell screen) | 0.73 | 0.65 | **PS** |
| ECCITE-seq protein prediction | 19/25 genes | 6/25 genes | **PS** |
| Partial effect detection | ✅ | ❌ | **PS** |

**Novel Discoveries**:
- 2,285 essential genes classified as buffered/sensitive
- Cell-state-dependent HIV reactivation mechanisms
- CCDC6 novel role in pancreatic vs liver fate decisions

---

## Weaknesses and Limitations

### 1. **Survival Bias Cannot Be Addressed**

**Problem**: Lethal perturbations kill strongly-affected cells before measurement

**Example**: Essential genes
- Cells with strong perturbation die
- Only measure cells with compensation/escape
- PS scores reflect survivors, not full spectrum

**Stated in Paper**: "PS probably only reflects perturbation responses in a fraction of cells, rather than the full spectrum"

**Mitigation**: Combine with prediction methods modeling uneven distributions

**Severity**: **HIGH** - Fundamental limitation for essential gene studies

### 2. **Requires Sufficient DEGs (Target Genes)**

**Problem**: PS performance degrades with few target genes

**Requirements**:
- Minimum: ~10 DEGs
- Recommended: 30-200 DEGs
- Optimal: Clear transcriptional signature

**Failure Modes**:
- Weak perturbations → few DEGs → noisy PS
- Post-transcriptional effects → no DEGs → PS fails
- Non-transcriptional phenotypes → no signature → cannot use PS

**Mitigation**: 
- Relax statistical thresholds to get more DEGs
- Use ECCITE-seq for protein readouts

**Severity**: **MEDIUM** - Excludes certain perturbation types

### 3. **Interpretability Challenges**

**PS Score Meaning**:
- Not absolute measure of "biological effect"
- Relative to target gene signature chosen
- Context-dependent (cell type, state)
- Scale depends on `scale_factor` parameter

**Questions Users Face**:
- What does PS = 0.5 mean biologically?
- How to compare PS across experiments?
- When is difference in PS meaningful?

**No Confidence Intervals**: Unlike SCEPTRE, no statistical uncertainty quantification

**Severity**: **MEDIUM** - Users need guidance on interpretation

### 4. **Computational Cost**

**Optimization Step**:
- Quadratic programming slower than simple DE tests
- Scales as O(NMG × iterations)
- 586K cells took hours (not instant)

**Memory Requirements**:
- Dense matrix operations
- Cannot easily use sparse matrices in optimization

**Comparison**:
- Mixscape: Minutes
- PS: Hours (for same dataset)

**Mitigation**: Parallelization across perturbations possible

**Severity**: **LOW-MEDIUM** - Acceptable for most screens, issue for ultra-large datasets

### 5. **Assumes Linearity and Additivity**

**Linearity Assumption**: y = y₀ + ψ×β

**Violations**:
- Highly nonlinear dose-responses (thresholds, saturation)
- Feedback loops
- Bistable systems

**Additivity Assumption**: Multiple perturbations add linearly

**Violations**:
- Epistasis (genetic interactions)
- Synergy/antagonism in combinatorial perturbations

**Impact**: PS may mis-estimate for complex regulatory networks

**Severity**: **MEDIUM** - Limits applicability to highly nonlinear systems

### 6. **Limited Statistical Framework**

**No P-Values**: PS provides scores, not hypothesis tests

**Cannot Directly Answer**:
- "Is this perturbation significant?"
- "What is my false discovery rate?"
- "How many cells are truly perturbed?"

**Comparison to SCEPTRE**:
- SCEPTRE: Calibrated p-values, FDR control
- PS: Continuous scores (requires thresholding)

**Use Case Mismatch**:
- Discovery screen → SCEPTRE better (controlled FDR)
- Heterogeneity characterization → PS better

**Severity**: **MEDIUM** - Depends on research question

### 7. **Not Validated for Low-MOI Edge Cases**

**SCEPTRE's Niche**: Low-MOI data with sparse perturbations

**PS Tested Mostly On**: High-quality, well-powered datasets

**Unknown Performance**:
- Very low cell counts (<10 per condition)
- Ultra-sparse perturbation matrices
- Single-cell-resolution interactions

**Severity**: **LOW** - Most Perturb-seq experiments have sufficient cells

---

## Comparison to Alternative Methods

### vs Mixscape

**Mixscape** (Seurat, Papalexi et al. 2021):

**What It Does**:
- Classifies cells as perturbed vs non-perturbed
- Uses LDA to find perturbation signature
- Posterior probability per cell

**Strengths**:
- Fast computation
- Built into Seurat (widely used)
- Good for binary classification

**Weaknesses**:
- No partial effect quantification (assigns probability ≈ 1.0 uniformly)
- Poor correlation with expression (0-5% of genes)
- Lower AUC scores vs PS

**When to Use**:
- Need binary classification only
- Very fast results required
- Simple visualization of perturbed vs control

**When PS Better**:
- Need dose-response info
- Studying incomplete perturbations
- Quantifying heterogeneity

### vs scMAGeCK-LR

**scMAGeCK-LR** (Li et al. 2018):

**What It Does**:
- Ridge regression to estimate average perturbation effects
- Gene-level β scores
- Ranking of perturbation impacts

**Strengths**:
- Average effect estimation (β scores)
- Gene set enrichment analysis
- Handles negative selection

**Weaknesses**:
- No cell-level scores
- Cannot detect heterogeneity
- Averages mask subpopulations

**Relationship to PS**:
- PS uses scMAGeCK-LR as Step 2 (estimate β)
- PS extends scMAGeCK-LR with cell scores

**When to Use scMAGeCK-LR**:
- Only need average effects
- Gene ranking sufficient
- No heterogeneity of interest

**When PS Better**:
- Need cell-level resolution
- Studying dose-responses
- Subpopulation discovery

### vs SCEPTRE

**SCEPTRE** (Barry et al. 2024):

**What It Does**:
- Association testing for low-MOI Perturb-seq
- Calibrated p-values
- FDR control

**Strengths**:
- Statistical rigor (calibrated tests)
- Low-MOI power
- FDR control for discovery

**Weaknesses**:
- No cell-level scores
- No dose-response quantification
- Binary (significant vs not)

**Complementarity**:
- SCEPTRE: Discovery (which genes affected?)
- PS: Characterization (how heterogeneous? what dosage?)

**When to Use SCEPTRE**:
- Discovery screen with FDR control
- Low-MOI data
- Need statistical significance

**When PS Better**:
- Characterizing known perturbations
- Dose-response analysis
- Heterogeneity is primary question

**Ideal Workflow**: SCEPTRE for discovery → PS for characterization

### vs MAST

**MAST** (Finak et al. 2015):

**What It Does**:
- Hurdle model for scRNA-seq DE
- Handles zero-inflation
- Flexible covariate adjustment

**Strengths**:
- General DE testing framework
- Handles zero-inflation explicitly
- Covariate adjustment

**Weaknesses**:
- Not specific to perturbations
- No cell-level scores
- Slower than simple tests

**Use Case Difference**:
- MAST: General DE testing
- PS: Perturbation-specific heterogeneity

**When to Use MAST**:
- DE testing with covariates
- Zero-inflation is major concern
- Not focused on perturbation dosage

### vs GSFA

**GSFA** (Guided Factor Analysis, Yao et al. 2023):

**What It Does**:
- Bayesian factor analysis
- Posterior inference
- Gene set enrichment

**Strengths**:
- Uncertainty quantification
- Handles gene correlations
- Interpretable factors

**Weaknesses**:
- Computationally expensive
- Complex model
- Many tuning parameters

**Comparison**:
- GSFA: Bayesian, uncertainty, factors
- PS: Optimization, deterministic, direct scores

**When to Use GSFA**:
- Need uncertainty quantification
- Factor interpretation important
- Computational time not limiting

**When PS Better**:
- Need fast results
- Interpretable cell scores priority
- Large-scale screens

### vs Mixscale

**Mixscale** (Extension of Mixscape, 2024):

**What It Does**:
- Continuous perturbation efficiency scores
- Extends Mixscape with scalar values
- Models downstream variation

**Strengths**:
- Improvement over Mixscape
- Continuous scoring (like PS)

**Weaknesses**:
- Not yet widely adopted/validated
- Unclear performance vs PS

**Comparison**:
- Both provide continuous scores
- PS has more extensive validation
- Mathematical framework differs

**Status**: Too new for definitive comparison (2024 preprint)

---

## When to Use PS vs Alternatives

### Decision Tree

```
START: What is your primary question?

├─ "Which genes are significantly affected?"
│   → Use: SCEPTRE (for discovery with FDR control)
│
├─ "What are the average effects of perturbations?"
│   → Use: scMAGeCK-LR (gene-level ranking)
│
├─ "Which cells are perturbed?" (binary classification)
│   → Use: Mixscape (fast, built into Seurat)
│
├─ "How strongly is each cell affected?" (dose-response)
│   → Use: PS ✓ (unique capability)
│
├─ "What causes heterogeneity in responses?"
│   ├─ Cell state?
│   ├─ CRISPR efficiency?
│   ├─ Compensation mechanisms?
│   └─ → Use: PS ✓ (characterization)
│
└─ "Which factors drive variation?"
    → Use: GSFA (factor analysis)
```

### Use Case Matrix

| Research Question | Best Method | Alternative | PS Applicable? |
|-------------------|-------------|-------------|----------------|
| Gene discovery (FDR control) | SCEPTRE | scMAGeCK-LR | No |
| Average perturbation effects | scMAGeCK-LR | MAST | No |
| Binary cell classification | Mixscape | threshold on PS | Yes (not optimal) |
| **Dose-response quantification** | **PS** | None | **Yes (unique)** |
| **Response heterogeneity** | **PS** | Mixscale? | **Yes (best)** |
| Partial perturbation detection | **PS** | None | **Yes (unique)** |
| Cell-state dependencies | **PS** | stratified DE | **Yes** |
| Subpopulation discovery | **PS** | clustering + DE | **Yes** |
| CRISPR QC (efficiency) | **PS** | sgRNA UMI | **Yes** |
| Epistasis detection | GSFA | linear models | Partly |
| Protein-level validation | **PS** | correlation | **Yes** |

---

## Complementary vs Competitive Methods

### Complementary (Use Together)

**PS + SCEPTRE**:
1. SCEPTRE: Identify significant perturbations
2. PS: Characterize heterogeneity of hits

**PS + scMAGeCK-LR**:
- Already integrated (PS uses LR for β estimation)
- LR provides gene rankings, PS adds cell scores

**PS + Seurat/Scanpy**:
- Standard preprocessing
- PS adds perturbation-specific analysis

### Competitive (Choose One)

**PS vs Mixscape**:
- Same input, different outputs
- PS supersedes Mixscape for heterogeneity questions

**PS vs Mixscale**:
- Both attempt continuous scoring
- Head-to-head comparison not yet published

---

## Recommended Workflows

### Workflow 1: Discovery Screen

**Goal**: Find genes affecting phenotype X

```
1. SCEPTRE: Test all genes, FDR < 0.05
2. scMAGeCK-LR: Rank hits by effect size
3. PS: Characterize top hits (dose-response, heterogeneity)
4. Validate: Protein, functional assays
```

### Workflow 2: Heterogeneity Study

**Goal**: Understand why same perturbation has different effects

```
1. PS: Calculate cell-level scores
2. Clustering: Group cells by PS
3. DE analysis: Find genes differing between PS-high vs PS-low
4. Functional analysis: Pathway enrichment, GO terms
```

### Workflow 3: CRISPR QC

**Goal**: Validate guide efficiency

```
1. PS: Score all cells
2. Correlation: PS vs perturbed gene expression
3. QC flag: Low correlation = poor guide
4. Filter: Remove low-efficiency guides
```

### Workflow 4: Drug Target Prioritization

**Goal**: Find genes with buffered vs sensitive responses

```
1. PS: Calculate scores for essential genes
2. Classify: PS distribution shape (bimodal vs uniform)
3. Buffered: Require high doses (challenging drug targets)
4. Sensitive: Partial inhibition sufficient (better targets)
```

---

## Future Improvements Needed

### 1. Uncertainty Quantification

**Current Gap**: No confidence intervals or credible regions

**Needed**: Bootstrap or Bayesian variant providing:
- Standard errors for PS scores
- Credible intervals
- Significance testing

### 2. Nonlinear Extensions

**Current Limitation**: Assumes linear dose-response

**Needed**: Nonparametric or piecewise-linear formulation:
- Sigmoid dose-responses
- Threshold effects
- Saturation

### 3. Epistasis Modeling

**Current Limitation**: Additive perturbation effects

**Needed**: Interaction terms:
- ψₖ₁ₖ₂ for combinatorial perturbations
- Network-based priors

### 4. Integration with Time-Series

**Current Limitation**: Single time point

**Needed**: Temporal PS:
- PS(t) over differentiation
- Infer dynamic responses

### 5. Low-MOI Calibration

**Current Gap**: Not explicitly validated for sparse data

**Needed**: SCEPTRE-style calibration analysis:
- Permutation tests
- Empirical null distributions
- FDR control

### 6. Computational Efficiency

**Current Limitation**: Hours for 500K+ cells

**Needed**: 
- Sparse matrix optimizations
- GPU acceleration
- Approximation methods (e.g., stochastic gradient)

---

## Critical Assessment Summary

### What PS Does Uniquely Well

1. ✅ Partial perturbation quantification
2. ✅ Dose-response characterization
3. ✅ Heterogeneity discovery (cell-state, cell-type)
4. ✅ Continuous cell-level scores

### What PS Cannot Do

1. ❌ Address survival bias
2. ❌ Statistical hypothesis testing (no p-values)
3. ❌ Model nonlinear dose-responses
4. ❌ Directly handle epistasis

### Overall Assessment

**Innovation Level**: ⭐⭐⭐⭐⭐ (5/5)
- First method for partial perturbation quantification
- Novel constrained optimization approach

**Validation Quality**: ⭐⭐⭐⭐⭐ (5/5)
- Published in Nature Cell Biology
- Multiple datasets, platforms
- Outperforms existing methods

**Practical Utility**: ⭐⭐⭐⭐ (4/5)
- Very useful for heterogeneity questions
- Moderate computational cost
- Interpretation requires expertise

**Statistical Rigor**: ⭐⭐⭐ (3/5)
- Strong optimization framework
- Limited uncertainty quantification
- No hypothesis testing

**Usability**: ⭐⭐⭐⭐ (4/5)
- Well-documented R package
- Reasonable defaults
- Integrates with Seurat

**Recommendation**: **Highly Recommended** for studies where:
- Response heterogeneity is the question
- Dose-response relationships are of interest
- Sufficient cells per condition (>50)
- Strong transcriptional perturbation signature

**Not Recommended** for:
- Discovery screens (use SCEPTRE)
- Weak/post-transcriptional perturbations
- Simple binary classification (use Mixscape)
- Ultra-low cell counts (<10)

---

**Document Version**: 1.0
**Last Updated**: 2025-11-13
**Purpose**: Critical evaluation for method comparison and selection

