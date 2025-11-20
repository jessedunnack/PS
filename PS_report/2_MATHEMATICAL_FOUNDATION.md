# PS Method - Mathematical and Statistical Foundation

**Document Purpose**: Detailed mathematical formulation, optimization algorithm, and statistical framework

---

## Table of Contents

1. [Problem Formulation](#problem-formulation)
2. [Mathematical Model](#mathematical-model)
3. [Constrained Quadratic Optimization](#constrained-quadratic-optimization)
4. [Algorithm Steps](#algorithm-steps)
5. [Statistical Framework](#statistical-framework)
6. [Computational Complexity](#computational-complexity)
7. [Parameter Selection](#parameter-selection)
8. [Mathematical Properties](#mathematical-properties)

---

## Problem Formulation

### Biological Setup

Given:
- **N cells** indexed by j = 1, ..., N
- **M perturbations** (genes) indexed by k = 1, ..., M
- **G target genes** indexed by i = 1, ..., G
- **Expression matrix** Y ∈ ℝ^(N×G)
- **Perturbation assignments** D ∈ {0,1}^(N×M)

Goal:
Estimate **ψⱼₖ** ∈ [0,1] = strength of perturbation k's effect in cell j

### Conceptual Framework

**Assumption**: Cell's expression profile is a mixture of:
1. Baseline expression (Y₀)
2. Perturbation-induced changes (scaled by PS score ψ)
3. Random noise (ε)

**Key Insight**: Different cells with same perturbation may show different response strengths due to:
- CRISPR efficiency variation
- Cell state differences
- Stochastic gene expression
- Compensation mechanisms

---

## Mathematical Model

### Step 1: Base Linear Model (scMAGeCK-LR)

Before PS calculation, estimate average perturbation effects:

```
Y = Y₀ + D × B + ε
```

Where:
- **Y** ∈ ℝ^(N×G): Log-transformed expression matrix
- **Y₀** ∈ ℝ^(N×G): Baseline expression (inferred from controls)
- **D** ∈ {0,1}^(N×M): Perturbation indicator matrix
  - Dⱼₖ = 1 if cell j has perturbation k, else 0
- **B** ∈ ℝ^(M×G): Effect coefficient matrix
  - βₖᵢ = average effect of perturbation k on gene i
- **ε** ~ N(0, σ²): Gaussian noise

**Solving for B** (Ridge Regression):
```
B = (D^T D + λI)^(-1) D^T Y
```

Where λ is regularization parameter.

### Step 2: PS-Enhanced Model

Replace binary perturbation matrix D with continuous PS scores Ψ:

```
Y = Y₀ + Ψ × B + ε
```

Where:
- **Ψ** ∈ ℝ^(N×M): PS score matrix (to be estimated)
  - ψⱼₖ = strength of perturbation k in cell j
  - Range: [0, U] where U is upper bound

**Key Difference from Step 1**:
- D is binary (0 or 1) - assumes uniform perturbation
- Ψ is continuous (0 to U) - allows heterogeneous responses

### Step 3: Normalization

Final PS scores normalized to [0,1]:

```
PSⱼₖ = ψⱼₖ / U
```

---

## Constrained Quadratic Optimization

### Objective Function

**Minimize**:
```
L(Ψ) = Σⱼ₌₁ᴺ Σᵢ₌₁ᴳ [yⱼᵢ - yⱼᵢ⁰ - Σₖ₌₁ᴹ ψⱼₖ βₖᵢ]² + λ Σⱼ₌₁ᴺ Σₖ₌₁ᴹ ψⱼₖ
```

**Term 1** (Reconstruction Error):
```
Σⱼᵢ [yⱼᵢ - yⱼᵢ⁰ - Σₖ ψⱼₖ βₖᵢ]²
```
- Minimizes difference between observed and predicted expression
- Penalizes poor fit to data

**Term 2** (L₁ Regularization):
```
λ Σⱼₖ ψⱼₖ
```
- Encourages sparsity (smaller PS scores when not strongly supported)
- Prevents overfitting
- λ > 0 is regularization strength

**Note on L₁ Formulation**:
Typically L₁ uses |ψⱼₖ|, but non-negativity constraint (ψⱼₖ ≥ 0) makes absolute value unnecessary:
```
|ψⱼₖ| = ψⱼₖ when ψⱼₖ ≥ 0
```

### Constraints

**1. Non-negativity**:
```
ψⱼₖ ≥ 0  for all j, k
```
Biological meaning: Perturbation cannot reduce response below baseline

**2. Presence Constraint**:
```
ψⱼₖ = 0  if Dⱼₖ = 0
```
Biological meaning: Cells without perturbation k cannot have non-zero PSₖ

**3. Upper Bound**:
```
ψⱼₖ ≤ U  if Dⱼₖ = 1
```
Where U is normalization constant (typically chosen based on data scale)

### Matrix Form

Rewrite objective as quadratic program:

Let **ψ** = vec(Ψ) be vectorized form (stacking columns)

**Minimize**:
```
½ ψ^T Q ψ + c^T ψ
```

Subject to:
```
A_eq ψ = b_eq    (equality constraints for ψⱼₖ = 0)
0 ≤ ψ ≤ u        (box constraints)
```

Where:
- **Q**: Hessian matrix (positive semi-definite)
- **c**: Linear term coefficients
- **A_eq**, **b_eq**: Encode presence constraints
- **u**: Vector of upper bounds

### Optimization Algorithm

**Method**: Newton's Method (or variant)

**Iterative Update**:
```
ψ^(t+1) = ψ^(t) - [∇²L(ψ^(t))]^(-1) ∇L(ψ^(t))
```

With projection onto constraint set.

**Gradient**:
```
∇L(ψⱼₖ) = -2 Σᵢ [yⱼᵢ - yⱼᵢ⁰ - Σₖ' ψⱼₖ' βₖ'ᵢ] βₖᵢ + λ
```

**Hessian**:
```
∇²L(ψⱼₖ, ψⱼₖ') = 2 Σᵢ βₖᵢ βₖ'ᵢ
```

---

## Algorithm Steps (Detailed)

### Input

- **Seurat object** (RDS): Normalized expression data
- **Barcode table**: Guide assignments per cell
- **Target perturbations**: List of genes to analyze
- **Control label**: Non-targeting control cells

### Step 1: Target Gene Identification

**Goal**: Find genes affected by perturbation

**Method**: Wilcoxon Rank-Sum Test

For each perturbation k:
1. Define cell sets:
   - S₊ = {cells with perturbation k}
   - S₋ = {control cells}

2. For each gene i:
   - Test H₀: Expression distribution same in S₊ vs S₋
   - Compute p-value using Wilcoxon test
   - Adjust for multiple testing (FDR)

3. Select significant genes:
   - FDR < 0.05 (typical threshold)
   - |log₂ FC| > threshold (optional)
   - Rank by significance and effect size

4. Filter target gene list:
   - Keep top N genes (N ∈ [target_gene_min, target_gene_max])
   - Default: 30 ≤ N ≤ 200

**Output**: Gene set Gₖ for perturbation k

### Step 2: Average Effect Estimation

**Goal**: Estimate average perturbation effects (β scores)

**Method**: Ridge Regression (scMAGeCK-LR)

1. Construct perturbation matrix D:
   - Dⱼₖ = 1 if cell j has guide for gene k
   - Dⱼₖ = 0 otherwise

2. Extract target gene expression:
   - Y_subset = Y[:, Gₖ] (only target genes)

3. Estimate baseline:
   - Y₀ = mean expression in control cells

4. Solve ridge regression:
   ```
   B = argmin_B ||Y - Y₀ - D×B||² + λ||B||²
   ```
   Closed form:
   ```
   B = (D^T D + λI)^(-1) D^T (Y - Y₀)
   ```

5. **Output**: Effect matrix B ∈ ℝ^(M×|Gₖ|)

### Step 3: PS Optimization

**Goal**: Estimate cell-specific response strengths Ψ

**Method**: Constrained Quadratic Programming

1. Initialize:
   - ψⱼₖ = Dⱼₖ (start with binary indicator)
   - Or ψⱼₖ ~ U(0, U) (random initialization)

2. Set up optimization:
   - Objective: L(Ψ) as defined above
   - Constraints:
     - ψⱼₖ = 0 if Dⱼₖ = 0
     - 0 ≤ ψⱼₖ ≤ U if Dⱼₖ = 1

3. Solve using quadratic programming solver:
   - E.g., CVXOPT, quadprog (R), scipy.optimize (Python)
   - Typically converges in 10-100 iterations

4. Normalize:
   ```
   PSⱼₖ = ψⱼₖ / U
   ```

**Output**: PS score matrix PS ∈ [0,1]^(N×M)

### Step 4: Post-processing

1. **Scaling** (optional):
   ```
   PS_scaled = scale_factor × PS
   ```
   - Default scale_factor = 3
   - Amplifies scores for visualization

2. **Add to Seurat metadata**:
   - Column name: "{gene}_eff"
   - Example: "TP53_eff", "BRD4_eff"

3. **Return**:
   - PS matrix
   - Updated Seurat object

---

## Statistical Framework

### Hypothesis Testing Analog

While PS provides continuous scores (not p-values), can interpret:

- **Null hypothesis**: ψⱼₖ = 0 (no perturbation effect)
- **Alternative**: ψⱼₖ > 0 (detectable effect)

**Implicit threshold**:
- Regularization term λ acts as penalty for non-zero scores
- Effectively sets minimum evidence required

### Confidence in PS Scores

**Higher confidence when**:
1. More target genes (larger Gₖ)
2. Stronger effect sizes (larger |βₖᵢ|)
3. Better model fit (lower residuals)
4. More cells per condition

**Uncertainty quantification**:
- Not directly provided by method
- Could bootstrap cell sampling for intervals
- Posterior credible intervals not available (not Bayesian)

### Relationship to Expression Correlation

**Expected**: PS ∼ -corr(perturbed gene expression, PS score)

**Validation approach**:
- For knockout: expect negative correlation
- If gene X knocked out, cells with high PS should have low X expression
- PS paper shows 40% of genes meet this (vs 0-5% for mixscape)

---

## Computational Complexity

### Time Complexity

**Step 1** (Target gene identification):
- Wilcoxon test per gene: O(n log n) where n = cells
- Total: O(G × n log n) for G genes

**Step 2** (Ridge regression):
- Matrix multiplication: O(M²N + MNG)
- Inversion: O(M³)
- Total: O(M³ + M²N + MNG)
- Typically fast (M, G << N usually)

**Step 3** (PS optimization):
- Per iteration: O(NMG) for gradient
- Convergence: typically 10-100 iterations
- Total: O(iterations × NMG)
- **Bottleneck** for large datasets

**Overall**: O(NMG × iterations)

### Space Complexity

- Expression matrix Y: O(NG)
- PS matrix Ψ: O(NM)
- Effect matrix B: O(MG)
- **Total**: O(NG + NM + MG) ≈ O(NG) typically

### Scalability

**Tested on**:
- 586,000 cells (genome-scale CRISPRi)
- 18,595 perturbations
- Successfully completed

**Limitations**:
- Memory: Dense matrix operations require RAM ~ O(NG)
- Time: Quadratic optimization slower than simple DE tests
- Parallelization: Can process perturbations independently

---

## Parameter Selection

### Key Parameters

#### 1. λ (Regularization Strength)

**Default**: 0.01

**Effect**:
- λ → 0: No regularization, risk overfitting
- λ → ∞: All PS scores → 0, underfitting

**Tuning**:
- Cross-validation on synthetic data
- Elbow plot: reconstruction error vs sparsity
- Typical range: [0, 0.1]

#### 2. U (Upper Bound)

**Purpose**: Normalization factor

**Typical value**: Determined by data scale

**Effect**:
- Final PS = ψ/U, so U sets the scale
- Larger U → smaller final PS values
- Usually set automatically by algorithm

#### 3. target_gene_min, target_gene_max

**Default**: min=30, max=200

**Effect**:
- Too few genes → noisy estimates
- Too many genes → diluted signal, slower computation

**Tuning**:
- Depends on effect size
- Strong perturbations: 20-50 genes sufficient
- Weak perturbations: 50-300 genes may be needed

#### 4. scale_factor

**Default**: 3

**Purpose**: Amplify PS scores for visualization/analysis

**Effect**:
- Only affects final output scale, not optimization
- Applied after normalization: PS_final = scale_factor × PS

**Tuning**:
- Weak effects (developmental): 5-10
- Moderate effects: 2-4
- Strong effects: 1-2

---

## Mathematical Properties

### Convexity

**Objective function** L(Ψ) is:
- **Convex** in Ψ (sum of convex terms)
- **Quadratic** form ensures unique global minimum

**Implications**:
- Optimization always converges
- No local minima issues
- Solution independent of initialization (in theory)

### Identifiability

**Question**: Is Ψ uniquely determined?

**Answer**: Partially

- Given B, Ψ is uniquely determined (convex optimization)
- Given Ψ, B is uniquely determined (ridge regression)
- Joint (B, Ψ) may have multiple solutions (trade-off)

**Resolution**: Two-step approach fixes this
- Step 2: Estimate B assuming uniform perturbation
- Step 3: Estimate Ψ given fixed B

### Sensitivity to Target Gene Selection

**Key question**: How does PS change with different target genes?

**Expectation**:
- More target genes → more stable estimates
- Including non-target genes → noise increases
- Optimal set: true causal downstream genes

**Robustness**:
- PS relatively stable with ≥30 genes
- Large changes (>0.2) suggest unstable signature

### Comparison to PCA/Matrix Factorization

**Similarity**:
- Both decompose expression into components
- PCA: Y ≈ U × V^T (unconstrained)
- PS: Y ≈ Ψ × B (constrained)

**Key Differences**:
- PS incorporates biological constraints (presence, non-negativity)
- PS uses pre-estimated B (informed by perturbation labels)
- PS scores directly interpretable as perturbation strength
- PCA components not biologically interpretable without rotation

---

## Model Assumptions

### Explicit Assumptions

1. **Linearity**: Expression changes linear in PS score
   - y ≈ y₀ + ψ×β
   - May not hold for highly nonlinear responses

2. **Additivity**: Multiple perturbations combine additively
   - For cell with k₁ and k₂: y ≈ y₀ + ψₖ₁β₁ + ψₖ₂β₂
   - Ignores epistasis/interactions

3. **Gaussian noise**: ε ~ N(0, σ²)
   - Reasonable for log-transformed counts
   - Violated for zero-inflated genes

4. **Shared effects**: βₖᵢ same across all cells with perturbation k
   - Allows PS to vary, but average effect fixed
   - Cell-type-specific effects require separate PS calculations

### Violations and Robustness

**When assumptions violated**:
- **Nonlinearity**: PS may not capture full dose-response curve
- **Non-additivity**: Combinatorial perturbations may be mis-scored
- **Heavy-tailed noise**: Outliers may inflate scores
- **Heterogeneous effects**: May need cell-type-stratified analysis

**Robustness strategies**:
- Use robust regression (Huber loss) instead of squared error
- Stratify by cell type before PS calculation
- Validate PS scores against independent readouts

---

## Comparison to Alternative Formulations

### Alternative 1: Bayesian Hierarchical Model

Could model:
```
yⱼᵢ ~ N(yⱼᵢ⁰ + Σₖ ψⱼₖ βₖᵢ, σ²)
ψⱼₖ ~ Beta(α, β) if Dⱼₖ = 1, else ψⱼₖ = 0
βₖᵢ ~ N(0, τ²)
```

**Advantages**:
- Uncertainty quantification (credible intervals)
- Hierarchical shrinkage

**Disadvantages**:
- Computational cost (MCMC or variational inference)
- More tuning parameters
- Slower for large datasets

### Alternative 2: Non-negative Matrix Factorization

Could use NMF:
```
Y ≈ W × H
```
Where W ≥ 0, H ≥ 0

Then interpret W as PS scores, H as effects.

**Advantages**:
- Non-negativity built in
- Unsupervised (no labels needed)

**Disadvantages**:
- No incorporation of perturbation labels
- Components not directly interpretable
- No presence constraints

### Why Current Formulation?

**Trade-offs**:
- Convex optimization → fast, reliable
- Constrained → biologically interpretable
- Regularized → prevents overfitting
- Two-step → identifiable solution

**Chosen for**: Balance of speed, interpretability, and accuracy

---

## Validation Against Ground Truth

### Synthetic Data Validation

**Setup**:
- Simulate scRNA-seq with known PS values
- Generate data: y = y₀ + ψ_true × β + noise
- Run PS method, compare ψ_estimated vs ψ_true

**Results**:
- Correlation: r > 0.8 for moderate noise
- RMSE: <0.15 for 100+ cells per condition
- Calibration: 50% knockdown → PS ≈ 0.33 (accurate)

### Experimental Validation

**ECCITE-seq** (protein as ground truth):
- PS predicts protein abundance
- Outperforms mixscape: 19/25 genes
- Strong correspondence for high-effect genes

**CRISPRi titration** (sgRNA mismatch):
- Mismatches reduce CRISPR efficiency
- PS correctly ranks: perfect > 1bp mismatch > 2bp > control
- Mixscape fails to differentiate

---

## Software Implementation Notes

### Optimization Solver

**R package**: scMAGeCK uses R's quadprog or similar

**Key functions**:
- `quadprog::solve.QP()` for convex QP
- Handles box constraints and equality constraints
- Efficient for moderate-scale problems (N×M < 10^6)

### Numerical Stability

**Issues**:
- Matrix inversion: (D^T D + λI)^(-1) may be ill-conditioned
- Very small λ → numerical errors

**Solutions**:
- Always use regularization (λ ≥ 0.001)
- Check condition number of D^T D
- Scale expression data to similar ranges

---

**Document Version**: 1.0
**Last Updated**: 2025-11-13
**Purpose**: Mathematical foundation for PS method evaluation
