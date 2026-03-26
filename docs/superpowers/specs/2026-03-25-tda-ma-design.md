# TDA-MA: Topological Data Analysis for Meta-Analysis — Design Spec

**Date:** 2026-03-25
**Author:** Mahmood Ahmad
**Target:** Single-file HTML app, ~3K lines, `C:\Models\TDA_MA\tda_ma.html`
**Paper target:** Research Synthesis Methods or JRSS Series A
**Phase:** 1 (Lean MVP — subgroup discovery + persistence barcode)

---

## 1. Purpose

TDA-MA is the world's first browser-based tool for applying topological data analysis to meta-analysis. It discovers hidden subgroups in study-level data using the Mapper algorithm and validates topological structure via persistent homology — without requiring the user to pre-specify subgroup variables. It then re-pools effects within each discovered subgroup and tests between-subgroup heterogeneity.

**Problem it solves:** Current subgroup analysis requires the analyst to specify which variable to stratify by (age, dose, risk of bias). If the true structure is multivariate or unexpected, it is missed. TDA finds structure the analyst didn't know to look for.

---

## 2. Architecture

Single-file HTML app following the project's standard pattern (embedded JS + CSS, no server, no external dependencies). All TDA computations in pure JavaScript.

### Pipeline

```
User data (CSV/paste/demo)
  → Feature selection (checkboxes for which dimensions to include)
  → Standardization (z-scores per dimension)
  → Distance matrix (Euclidean or Mahalanobis)
  → ├─ Mapper graph (lens → overlapping bins → cluster → nerve)
    ├─ Persistent homology (Rips filtration → Union-Find H0 + boundary reduction H1)
    └─ Cluster re-pooling (DL meta-analysis per subgroup + between-subgroup Q-test)
```

### 5 Tabs

1. **Data** — demo datasets, CSV import, column mapping, dimension selection, distance metric
2. **Mapper** — interactive network graph (hero visualization), force-directed layout
3. **Persistence** — barcode (default) + persistence diagram (Advanced toggle)
4. **Subgroups** — forest plots per cluster, Q-test, covariate profiles
5. **Export** — R code, CSV cluster assignments, PNG figures, JSON results

---

## 3. TDA Engine

### 3.1 Distance Matrix

- User selects dimensions via checkboxes (effect, SE, log(n), year, RoB, etc.)
- Each dimension standardized to z-scores (subtract mean, divide by SD)
- Euclidean distance (default): `d(i,j) = sqrt(sum((z_i - z_j)^2))`
- Mahalanobis distance (option): `d(i,j) = sqrt((z_i - z_j)^T S^{-1} (z_i - z_j))` where S is the covariance matrix of the selected dimensions
- **Fallback:** If k ≤ number of selected dimensions (singular covariance matrix), automatically fall back to Euclidean with a console warning
- Output: symmetric k×k matrix, k = number of studies

### 3.2 Persistent Homology (Vietoris-Rips Filtration)

**H0 (connected components):**
- Sort all pairwise distances ascending
- Process edges in order using Union-Find
- Each study is born at epsilon=0
- A component "dies" when it merges into an older component
- Track (birth, death) pairs

**H1 (loops/cycles):**
- Build boundary matrix for 2-simplices (triangles)
- Reduce via column operations (standard persistence algorithm)
- A 1-cycle is born when a triangle's boundary creates a new cycle
- It dies when a higher simplex fills it
- Track (birth, death) pairs

**Significance:** Features with persistence > 2× median persistence are flagged as significant (simplified bootstrap threshold per Fasy et al. 2014).

**Complexity:** O(k² log k) for H0 via Union-Find. H1 uses boundary matrix reduction — O(k³) worst case, but meta-analyses rarely exceed k=200 studies, so this is fine in-browser.

### 3.3 Mapper Algorithm

**Lens function** (user-selectable):
- First principal component of selected dimensions (default)
- Effect size
- Standard error
- Publication year
- Any selected numeric column

**Cover construction:**
- Divide lens range into N intervals (default N=10, slider: 5-20)
- Overlap between adjacent intervals (default 50%, slider: 25-75%)

**Clustering within intervals:**
- Single-linkage hierarchical clustering on the distance matrix restricted to studies in each interval
- Automatic cutoff: first gap > 1.5× median nearest-neighbor distance within the interval
- Each cluster within each interval becomes a Mapper node

**Nerve construction:**
- Create one node per cluster
- Edge between two nodes if they share ≥1 study (a study can appear in overlapping intervals)
- Node size proportional to number of studies
- Node color = mean effect size in that cluster (diverging colormap: blue = protective, white = null, red = harmful)

### 3.4 Cluster Extraction and Re-Pooling

- Connected components of the Mapper graph = topological subgroups
- If only 1 component: spectral bisection (cut at longest edge)
- Each subgroup: DerSimonian-Laird random-effects pooled estimate + 95% CI + I² + tau²
- Between-subgroup test: Q_between = Q_total - sum(Q_within), df = n_clusters - 1, p from chi-squared
- Covariate profile table: mean/SD of each input dimension per cluster

---

## 4. Demo Datasets

### 4.1 BCG Vaccine (13 studies)

Known moderator: latitude (tropical vs temperate climates). TDA should discover 2 clusters without being told about latitude.

Fields: `study, yi (log RR), sei, latitude, year, alloc_method`

Expected output: Mapper splits into tropical (weak/no effect) and temperate (strong protective effect) clusters. Persistence barcode shows one dominant H0 feature.

Validation: TDA clusters vs known latitude split should agree >80%.

### 4.2 Magnesium for Acute MI (16 studies)

Known issue: small-study effects + one mega-trial (ISIS-4) that reversed the pooled conclusion.

Fields: `study, yi (log OR), sei, n, year, blinding`

Expected output: ISIS-4 isolated as topologically distant from the cluster of small positive trials. Persistence diagram shows ISIS-4 as a long-lived H0 feature.

### 4.3 Simulated Clusters (20 studies)

Planted structure: 10 studies with effect ~ N(-0.5, 0.1), 10 with effect ~ N(0.1, 0.1), plus 1 bridge study.

Expected output: TDA recovers 2 clusters with >90% accuracy. Between-cluster Q-test significant.

---

## 5. UI Design

### Tab 1 — Data

- Dropdown: "BCG Vaccine", "Magnesium MI", "Simulated Clusters", "Custom"
- Custom: textarea for CSV paste + file upload button
- Auto-detect columns or manual assignment via dropdowns
- Dimension checkboxes: all numeric columns listed, effect + SE checked by default
- Distance metric toggle: Euclidean / Mahalanobis
- "Run TDA" button (primary action)

### Tab 2 — Mapper Graph (hero tab)

- Force-directed network layout rendered on HTML5 Canvas
- Node size = study count, node color = mean effect (diverging blue-white-red)
- Edge thickness = shared study count
- Hover tooltip: study IDs, mean effect [95% CI], covariate means
- Click node: highlights cluster studies, scrolls to forest plot in Subgroups tab
- Control panel:
  - Lens function dropdown (PC1, effect, SE, year)
  - Bins slider (5-20, default 10)
  - Overlap slider (25-75%, default 50%)
- Live re-computation on parameter change

### Tab 3 — Persistence

- Default: barcode visualization (horizontal bars)
  - H0 bars in blue, H1 bars in orange
  - Significant features (persistence > 2× median) drawn thicker with label
- Toggle: persistence diagram (birth vs death scatter with diagonal)
- Annotation text: "K significant connected components detected" / "M significant loops detected"
- Explanation tooltip: "A long bar means a robust topological feature. A short bar means noise."

### Tab 4 — Subgroups

- One forest plot per Mapper-discovered cluster (standard forest plot: study labels, effect + CI lines, diamond for pooled)
- Summary table:

  | Cluster | k | Pooled Effect | 95% CI | I² | tau² |
  |---------|---|--------------|--------|-----|------|

- Between-cluster Q-test: Q = X, df = K-1, p = Y
- Covariate profile table: mean ± SD of each input dimension per cluster
- Text summary: "TDA identified K topologically distinct subgroups. Between-cluster heterogeneity: Q = X, df = K-1, p = Y"
- **Exploratory warning:** "Note: TDA subgroup discovery is exploratory (hypothesis-generating). Clusters are data-driven and should be validated in independent datasets before informing clinical decisions."

### Tab 5 — Export

- **R code**: generates equivalent analysis using `TDAstats` (persistent homology) + `TDAmapper` (Mapper) + `metafor` (re-pooling)
- **CSV**: study ID, cluster assignment, effect, SE, all covariates
- **PNG**: Mapper graph, persistence barcode, forest plots (via canvas-to-blob)
- **JSON**: full results object (distance matrix, persistence pairs, Mapper graph, cluster assignments, pooled estimates)

### Standard Features

- Dark mode toggle (CSS custom properties)
- Responsive layout (works on tablet)
- Keyboard accessible (tab navigation, Enter to activate)
- Seeded PRNG: xoshiro128** for reproducible clustering (user can set seed)
- `escapeHtml()` on all user input rendered to DOM

---

## 6. Testing Strategy

### Unit Tests (in-browser test suite, run via button)

1. Distance matrix: symmetry, triangle inequality, zero diagonal
2. Persistence: birth ≤ death for all features
3. Union-Find: correctness on known 5-point example
4. Mapper: nerve is valid (no duplicate edges, node IDs match)
5. DL pooling: matches metafor on BCG dataset (tolerance 1e-4)
6. Q-test: matches manual computation
7. Simulated dataset: recovers planted clusters (≥90% accuracy)
8. BCG: discovers latitude-correlated clusters (≥80% agreement)

### Validation

- R cross-validation script using `TDAstats::calculate_homology()` for persistence
- Compare persistence diagrams: same (birth, death) pairs within tolerance
- Compare Mapper output against `TDAmapper::mapper()` on BCG dataset

---

## 7. Phase 2 Roadmap (future, not in MVP)

- Full persistence diagram visualization with interactive filtration slider
- Multiple filtration types: Alpha complex, Witness complex
- UMAP/t-SNE dimensionality reduction for 3D point cloud visualization
- Sheaf-theoretic NMA consistency analysis (companion module)
- Conformal prediction intervals per subgroup
- WebR validation tier
- Integration with TruthCert/MetaSprint import
- Bottleneck distance between persistence diagrams of different meta-analyses

---

## 8. Success Criteria

1. BCG demo correctly separates tropical vs temperate studies without latitude as input
2. Simulated demo recovers planted clusters with ≥90% accuracy
3. Mapper graph renders interactively with <500ms update on parameter change
4. All unit tests pass
5. R cross-validation matches persistence pairs within tolerance
6. Re-pooled estimates match metafor within 1e-4
7. Single HTML file, no external dependencies, works offline
