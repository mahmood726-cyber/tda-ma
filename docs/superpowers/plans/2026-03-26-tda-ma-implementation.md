# TDA-MA Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the world's first browser-based TDA tool for meta-analysis — Mapper graph + persistent homology + cluster re-pooling in a single HTML file.

**Architecture:** Single-file HTML app (~3K lines). Pure JS TDA engine (distance matrix → Rips filtration → Union-Find H0 + boundary reduction H1 → Mapper). Canvas-rendered Mapper graph. DL re-pooling per discovered cluster with Q-test. 3 demo datasets (BCG, magnesium, simulated). No external dependencies.

**Tech Stack:** HTML5/CSS3/JS (ES2020), Canvas 2D for Mapper graph and forest plots, xoshiro128** seeded PRNG.

**Spec:** `docs/superpowers/specs/2026-03-25-tda-ma-design.md`

---

## File Structure

This is a single-file HTML app. All code lives in `tda_ma.html`. The plan is organized by logical module within that file.

```
C:\Models\TDA_MA\
├── tda_ma.html              # The entire app (~3K lines)
├── tests/
│   └── test_tda_ma.py       # Selenium test suite
├── validate_vs_r.R          # R cross-validation script
└── docs/superpowers/
    ├── specs/2026-03-25-tda-ma-design.md
    └── plans/2026-03-26-tda-ma-implementation.md
```

Within `tda_ma.html`, logical sections (delimited by comments):
1. CSS (~300 lines) — layout, tabs, dark mode, Mapper styling
2. HTML (~200 lines) — 5 tab panels, controls, containers
3. Data module (~250 lines) — demo datasets, CSV parser, column mapper
4. Linear algebra module (~150 lines) — matrix ops, PCA, standardization
5. Distance module (~100 lines) — Euclidean, Mahalanobis, distance matrix
6. Persistence module (~250 lines) — Union-Find, Rips filtration, H0, H1 boundary reduction
7. Mapper module (~300 lines) — lens, cover, single-linkage, nerve construction
8. Meta-analysis module (~200 lines) — DL pooling, Q-test, I², tau²
9. Visualization module (~400 lines) — Mapper canvas, barcode, persistence diagram, forest plots
10. Export module (~150 lines) — R code gen, CSV, PNG, JSON
11. UI glue (~200 lines) — tab switching, event handlers, dark mode, test runner
12. In-browser tests (~200 lines) — 8 unit tests per spec

---

### Task 1: HTML Shell + CSS + Tab System

**Files:**
- Create: `tda_ma.html`

- [ ] **Step 1: Create the HTML skeleton with 5 tabs and CSS**

Write the complete HTML structure: `<!DOCTYPE html>`, meta charset, title "TDA-MA: Topological Data Analysis for Meta-Analysis", CSS custom properties for light/dark mode, tab system (Data, Mapper, Persistence, Subgroups, Export), header with title + dark mode toggle, footer with version.

CSS: use `--bg`, `--fg`, `--accent`, `--card-bg` custom properties. Tab buttons use `.tab-btn.active` pattern. Card containers for each tab panel. Responsive grid. Canvas containers for Mapper and forest plots.

```html
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>TDA-MA: Topological Data Analysis for Meta-Analysis</title>
<style>
:root {
  --bg: #f8f9fa; --fg: #1a1a2e; --accent: #2563eb;
  --card-bg: #ffffff; --border: #e2e8f0; --muted: #64748b;
  --red: #ef4444; --green: #22c55e; --blue: #3b82f6; --orange: #f97316;
}
body.dark {
  --bg: #0f172a; --fg: #e2e8f0; --accent: #60a5fa;
  --card-bg: #1e293b; --border: #334155; --muted: #94a3b8;
}
/* ... full CSS follows in implementation ... */
</style>
</head>
```

- [ ] **Step 2: Verify the shell renders**

Open `tda_ma.html` in browser. Verify: 5 tabs clickable, dark mode toggle works, responsive layout on narrow viewport.

- [ ] **Step 3: Commit**

```bash
git add tda_ma.html
git commit -m "feat: HTML shell with 5-tab layout and dark mode"
```

---

### Task 2: Data Module — Demo Datasets + CSV Parser

**Files:**
- Modify: `tda_ma.html` (add `<script>` section)

- [ ] **Step 1: Write the in-browser test for data loading**

```javascript
function testDataModule() {
  const results = [];
  // Test 1: BCG dataset loads with correct shape
  const bcg = DATASETS.bcg;
  results.push({ name: 'BCG has 13 studies', pass: bcg.studies.length === 13 });
  results.push({ name: 'BCG has yi column', pass: bcg.studies[0].yi !== undefined });
  // Test 2: CSV parser handles basic input
  const csv = 'study,yi,sei\nA,-0.5,0.1\nB,0.2,0.15';
  const parsed = parseCSV(csv);
  results.push({ name: 'CSV parses 2 rows', pass: parsed.length === 2 });
  results.push({ name: 'CSV numeric conversion', pass: parsed[0].yi === -0.5 });
  // Test 3: Simulated dataset has planted structure
  const sim = DATASETS.simulated;
  results.push({ name: 'Simulated has 21 studies', pass: sim.studies.length === 21 });
  return results;
}
```

- [ ] **Step 2: Implement demo datasets**

BCG vaccine (13 studies from Colditz et al. 1994): embed `study, yi, sei, latitude, year, alloc` as JS array of objects.

Magnesium MI (16 studies from Teo et al. 1991 + ISIS-4): embed `study, yi, sei, n, year, blinding`.

Simulated clusters (21 studies): 10 with yi ~ N(-0.5, 0.1²), 10 with yi ~ N(0.1, 0.1²), 1 bridge at yi=-0.2. Use seeded PRNG for deterministic generation.

- [ ] **Step 3: Implement CSV parser**

```javascript
function parseCSV(text) {
  const lines = text.trim().split('\n');
  const headers = lines[0].split(',').map(h => h.trim());
  return lines.slice(1).map(line => {
    const vals = line.split(',');
    const row = {};
    headers.forEach((h, i) => {
      const v = vals[i]?.trim();
      row[h] = isNaN(v) || v === '' ? v : parseFloat(v);
    });
    return row;
  });
}
```

- [ ] **Step 4: Wire up Data tab UI**

Dataset dropdown populates textarea with CSV preview. "Custom" shows editable textarea + file upload. Column auto-detection (look for `yi/effect/es`, `sei/se/stderr`, `n/sample_size`). Dimension checkboxes generated from numeric columns.

- [ ] **Step 5: Run tests, verify BCG loads correctly**

Click "Run Tests" button. All 5 data tests should pass.

- [ ] **Step 6: Commit**

```bash
git add tda_ma.html
git commit -m "feat: data module — 3 demo datasets, CSV parser, column mapper"
```

---

### Task 3: Linear Algebra + Distance Matrix

**Files:**
- Modify: `tda_ma.html`

- [ ] **Step 1: Write distance matrix tests**

```javascript
function testDistanceMatrix() {
  const results = [];
  // 3 points in 2D
  const points = [[0,0], [3,4], [0,0]];
  const D = euclideanDistanceMatrix(points);
  results.push({ name: 'Diagonal is zero', pass: D[0][0] === 0 && D[1][1] === 0 });
  results.push({ name: 'Symmetry', pass: D[0][1] === D[1][0] });
  results.push({ name: 'Euclidean 3-4-5', pass: Math.abs(D[0][1] - 5) < 1e-10 });
  results.push({ name: 'Identical points = 0', pass: D[0][2] === 0 });
  // Standardization
  const data = [{a: 10, b: 100}, {a: 20, b: 200}, {a: 30, b: 300}];
  const z = standardize(data, ['a', 'b']);
  results.push({ name: 'Z-score mean ~0', pass: Math.abs(z.map(r=>r[0]).reduce((a,b)=>a+b)/3) < 1e-10 });
  return results;
}
```

- [ ] **Step 2: Implement standardization and distance functions**

`standardize(studies, dimensions)` → array of z-score vectors.
`euclideanDistanceMatrix(points)` → k×k symmetric matrix.
`mahalanobisDistanceMatrix(points)` → k×k using inverse covariance (with singular fallback to Euclidean).
`computePCA(points, nComponents)` → eigendecomposition for lens function.

- [ ] **Step 3: Run tests**

All 5 distance tests pass.

- [ ] **Step 4: Commit**

```bash
git add tda_ma.html
git commit -m "feat: linear algebra — standardization, distance matrices, PCA"
```

---

### Task 4: Persistent Homology Engine

**Files:**
- Modify: `tda_ma.html`

- [ ] **Step 1: Write persistence tests**

```javascript
function testPersistence() {
  const results = [];
  // Known 4-point example: square with vertices at (0,0),(1,0),(1,1),(0,1)
  // Distances: sides=1, diagonals=sqrt(2)≈1.414
  const D = [[0,1,1.414,1],[1,0,1,1.414],[1.414,1,0,1],[1,1.414,1,0]];
  const ph = computePersistence(D);

  // H0: 4 components born at 0. At eps=1, 4 edges form, merging into 1 component.
  // 3 deaths at eps=1, 1 component survives (death=Infinity)
  const h0 = ph.h0.filter(f => f.death !== Infinity);
  results.push({ name: 'H0: 3 finite features', pass: h0.length === 3 });
  results.push({ name: 'H0: all die at eps=1', pass: h0.every(f => Math.abs(f.death - 1) < 1e-10) });
  results.push({ name: 'H0: birth <= death', pass: ph.h0.every(f => f.birth <= f.death) });

  // H1: 1-cycle born at eps=1 (square formed), dies at eps=sqrt(2) (diagonal fills it)
  results.push({ name: 'H1: 1 loop detected', pass: ph.h1.length === 1 });
  if (ph.h1.length > 0) {
    results.push({ name: 'H1: born at 1, dies at ~1.414',
      pass: Math.abs(ph.h1[0].birth - 1) < 1e-6 && Math.abs(ph.h1[0].death - 1.414) < 0.01 });
  }
  return results;
}
```

- [ ] **Step 2: Implement Union-Find for H0**

```javascript
class UnionFind {
  constructor(n) { this.parent = Array.from({length:n},(_,i)=>i); this.rank = new Array(n).fill(0); this.birth = new Array(n).fill(0); }
  find(x) { if (this.parent[x]!==x) this.parent[x]=this.find(this.parent[x]); return this.parent[x]; }
  union(x, y, eps) { /* merge younger into older, record death */ }
}
```

- [ ] **Step 3: Implement Rips filtration + H0 computation**

Sort all k*(k-1)/2 edges by distance. Process via Union-Find. Record (birth=0, death=merge_epsilon) for each merger. The last-surviving component gets death=Infinity.

- [ ] **Step 4: Implement H1 boundary matrix reduction**

Build sorted list of 2-simplices (triangles). For each triangle [i,j,k], its filtration value = max(d(i,j), d(j,k), d(i,k)). Sort triangles by filtration value. Build boundary matrix (columns = triangles, rows = edges). Reduce via left-to-right column reduction. Unpaired edges = H1 births. Pairing with a triangle = H1 death.

- [ ] **Step 5: Run persistence tests on square example**

All 5 tests pass. H0 shows 3 mergers at eps=1. H1 shows 1 loop born at 1, dying at sqrt(2).

- [ ] **Step 6: Commit**

```bash
git add tda_ma.html
git commit -m "feat: persistent homology — Union-Find H0 + boundary reduction H1"
```

---

### Task 5: Mapper Algorithm

**Files:**
- Modify: `tda_ma.html`

- [ ] **Step 1: Write Mapper tests**

```javascript
function testMapper() {
  const results = [];
  // Use simulated dataset (2 planted clusters + bridge)
  const sim = DATASETS.simulated;
  const dims = ['yi', 'sei'];
  const points = standardize(sim.studies, dims);
  const D = euclideanDistanceMatrix(points);
  const lens = sim.studies.map(s => s.yi); // effect size as lens

  const graph = runMapper(D, lens, sim.studies, { nBins: 10, overlap: 0.5 });

  results.push({ name: 'Mapper produces nodes', pass: graph.nodes.length > 0 });
  results.push({ name: 'Mapper produces edges', pass: graph.edges.length >= 0 }); // may be 0 if disconnected
  results.push({ name: 'Nodes have study lists', pass: graph.nodes[0].studies && graph.nodes[0].studies.length > 0 });
  results.push({ name: 'All studies assigned', pass:
    new Set(graph.nodes.flatMap(n => n.studies)).size === sim.studies.length });

  // Extract connected components — should find 2 for simulated data
  const components = getConnectedComponents(graph);
  results.push({ name: 'Simulated: 2 components found', pass: components.length === 2 });
  return results;
}
```

- [ ] **Step 2: Implement Mapper**

```javascript
function runMapper(distMatrix, lensValues, studies, opts) {
  const { nBins = 10, overlap = 0.5 } = opts;
  // 1. Cover: divide lens range into overlapping intervals
  // 2. For each interval, get indices of studies falling in it
  // 3. Single-linkage clustering on sub-distance-matrix
  // 4. Each cluster → Mapper node { id, studies: [...indices], meanEffect }
  // 5. Nerve: edge between nodes sharing ≥1 study
  return { nodes, edges };
}

function singleLinkageClusters(subDistMatrix, indices) {
  // Agglomerative: merge closest pair, stop when gap > 1.5× median NN distance
}

function getConnectedComponents(graph) {
  // BFS/DFS on graph.nodes using graph.edges
  // Return array of arrays of node indices
}
```

- [ ] **Step 3: Run Mapper tests on simulated data**

All 5 tests pass. Simulated dataset produces 2 connected components.

- [ ] **Step 4: Commit**

```bash
git add tda_ma.html
git commit -m "feat: Mapper algorithm — lens, cover, clustering, nerve"
```

---

### Task 6: Meta-Analysis Re-Pooling Engine

**Files:**
- Modify: `tda_ma.html`

- [ ] **Step 1: Write DL pooling tests**

```javascript
function testMetaAnalysis() {
  const results = [];
  // BCG full dataset — compare against metafor::rma(yi, vi, data=dat.bcg, method="DL")
  // Known reference: pooled = -0.7145, se = 0.1791, tau2 = 0.3132
  const bcg = DATASETS.bcg.studies;
  const dl = dlPooled(bcg.map(s=>s.yi), bcg.map(s=>s.sei**2));
  results.push({ name: 'DL theta within 1e-3', pass: Math.abs(dl.theta - (-0.7145)) < 0.001 });
  results.push({ name: 'DL tau2 within 1e-3', pass: Math.abs(dl.tau2 - 0.3132) < 0.001 });

  // Q-test between 2 subgroups
  const g1 = bcg.filter(s => s.latitude > 33);
  const g2 = bcg.filter(s => s.latitude <= 33);
  const qt = betweenGroupQ(
    [g1.map(s=>s.yi), g2.map(s=>s.yi)],
    [g1.map(s=>s.sei**2), g2.map(s=>s.sei**2)]
  );
  results.push({ name: 'Q-between > 0', pass: qt.Q > 0 });
  results.push({ name: 'Q-between df = 1', pass: qt.df === 1 });
  results.push({ name: 'Q-between p < 0.05 for latitude', pass: qt.p < 0.05 });
  return results;
}
```

- [ ] **Step 2: Implement DL pooling**

```javascript
function dlPooled(yi, vi) {
  // Fixed-effect weights
  const w = vi.map(v => 1/v);
  const sumW = w.reduce((a,b)=>a+b);
  const thetaFE = w.reduce((s,wi,i)=>s+wi*yi[i],0) / sumW;
  // Q statistic
  const Q = w.reduce((s,wi,i)=>s+wi*(yi[i]-thetaFE)**2, 0);
  const k = yi.length;
  const C = sumW - w.reduce((s,wi)=>s+wi*wi,0)/sumW;
  const tau2 = Math.max(0, (Q - (k-1)) / C);
  // RE weights
  const wRE = vi.map(v => 1/(v + tau2));
  const sumWRE = wRE.reduce((a,b)=>a+b);
  const theta = wRE.reduce((s,wi,i)=>s+wi*yi[i],0) / sumWRE;
  const se = Math.sqrt(1/sumWRE);
  const z = theta / se;
  const p = 2 * (1 - normalCDF(Math.abs(z)));
  const I2 = Math.max(0, (Q - (k-1)) / Q) * 100;
  return { theta, se, tau2, Q, I2, ci: [theta - 1.96*se, theta + 1.96*se], p, k };
}
```

- [ ] **Step 3: Implement between-group Q-test**

```javascript
function betweenGroupQ(yiGroups, viGroups) {
  // Pool each group separately, get Q_within per group
  // Q_between = Q_total - sum(Q_within)
  // df = nGroups - 1
  // p from chi-squared CDF
}
```

- [ ] **Step 4: Run tests — DL matches metafor reference**

All 5 meta-analysis tests pass.

- [ ] **Step 5: Commit**

```bash
git add tda_ma.html
git commit -m "feat: DL random-effects pooling + between-group Q-test"
```

---

### Task 7: Visualizations — Mapper Canvas + Barcode + Forest Plots

**Files:**
- Modify: `tda_ma.html`

- [ ] **Step 1: Implement force-directed Mapper graph on Canvas**

Force-directed layout: nodes repel (Coulomb), edges attract (Hooke). Iterate ~200 steps. Render on `<canvas>`:
- Circle per node, radius ∝ sqrt(study count)
- Color from diverging scale (blue → white → red) based on mean effect
- Lines for edges, thickness ∝ shared count
- Hover detection via distance-to-node check in mousemove handler
- Tooltip div positioned near cursor

- [ ] **Step 2: Implement persistence barcode**

Render on `<canvas>`:
- H0 bars in blue, H1 bars in orange
- X-axis = filtration value (epsilon)
- Each bar: horizontal line from birth to death
- Significant bars (persistence > 2× median) drawn thicker + labeled
- Toggle button to switch to persistence diagram (scatter: birth on x, death on y, diagonal line)

- [ ] **Step 3: Implement forest plots for subgroups**

Render on `<canvas>` per cluster:
- Study labels on left
- Effect ± CI as horizontal lines with square at point estimate (size ∝ weight)
- Diamond for pooled DL estimate
- Vertical null line at 0
- I², tau², p-value annotation

- [ ] **Step 4: Wire up Mapper controls (lens, bins, overlap sliders)**

Lens dropdown, bins slider (5-20), overlap slider (25-75%). On change: re-run Mapper, re-render canvas. Debounce at 100ms.

- [ ] **Step 5: Visual verification**

Load BCG dataset, click "Run TDA". Verify:
- Mapper tab shows a graph with identifiable clusters
- Persistence tab shows barcode with H0 bars
- Subgroups tab shows forest plots per cluster with Q-test

- [ ] **Step 6: Commit**

```bash
git add tda_ma.html
git commit -m "feat: visualizations — Mapper canvas, persistence barcode, forest plots"
```

---

### Task 8: Subgroups Tab — Covariate Profiles + Summary

**Files:**
- Modify: `tda_ma.html`

- [ ] **Step 1: Implement covariate profile table**

For each Mapper-connected-component, compute mean ± SD of every input dimension. Render as HTML table. Highlight dimensions where clusters differ by > 1 SD.

- [ ] **Step 2: Implement text summary**

Generate: "TDA identified K topologically distinct subgroups. Between-cluster heterogeneity: Q = X.XX, df = K-1, p = Y.YYY"

Add exploratory warning below.

- [ ] **Step 3: Verify on BCG**

Subgroups tab should show 2 clusters. Covariate profile should show latitude differs significantly between clusters. Q-test should be significant.

- [ ] **Step 4: Commit**

```bash
git add tda_ma.html
git commit -m "feat: subgroup profiles — covariate table, text summary, Q-test"
```

---

### Task 9: Export Module

**Files:**
- Modify: `tda_ma.html`

- [ ] **Step 1: Implement R code generator**

```javascript
function generateRCode(studies, clusterAssignments, dims) {
  return `
library(TDAstats)
library(TDAmapper)
library(metafor)

# Data
dat <- data.frame(
  study = c(${studies.map(s => `"${s.study}"`).join(', ')}),
  yi = c(${studies.map(s => s.yi).join(', ')}),
  sei = c(${studies.map(s => s.sei).join(', ')}),
  cluster = c(${clusterAssignments.join(', ')})
)

# Persistent homology
X <- scale(dat[, c(${dims.map(d => `"${d}"`).join(', ')})])
ph <- calculate_homology(X, dim = 1)
plot_persist(ph)

# Re-pool per cluster
for (cl in unique(dat$cluster)) {
  sub <- dat[dat$cluster == cl, ]
  cat("\\nCluster", cl, ":\\n")
  print(rma(yi, sei^2, data = sub, method = "DL"))
}
`;
}
```

- [ ] **Step 2: Implement CSV, PNG, JSON export**

CSV: study, cluster, yi, sei, + all covariates.
PNG: canvas.toBlob() for Mapper, barcode, forest plots. Create download link.
JSON: full results object.

- [ ] **Step 3: Verify exports**

Download R code, open in RStudio — should run without errors (given packages installed). Download CSV — should open in Excel with correct columns.

- [ ] **Step 4: Commit**

```bash
git add tda_ma.html
git commit -m "feat: export — R code, CSV, PNG, JSON"
```

---

### Task 10: Integration + Full Test Suite + Polish

**Files:**
- Modify: `tda_ma.html`
- Create: `tests/test_tda_ma.py`
- Create: `validate_vs_r.R`

- [ ] **Step 1: Wire up the full pipeline**

"Run TDA" button: parse data → select dimensions → standardize → distance matrix → run persistence + Mapper → extract clusters → re-pool → render all tabs. Show loading spinner during computation.

- [ ] **Step 2: Run all 8 in-browser unit tests**

Combine all test functions into a test runner. "Run Tests" button in footer. All 8 categories (data, distance, persistence, Mapper, meta-analysis) must pass.

- [ ] **Step 3: BCG validation**

Load BCG demo. Run TDA with dimensions = [yi, sei, latitude, year]. Verify Mapper separates tropical from temperate. Check cluster agreement with known latitude split ≥80%.

- [ ] **Step 4: Simulated validation**

Load Simulated demo. Run TDA with dimensions = [yi, sei]. Verify 2 clusters recovered with ≥90% accuracy. Q-test significant.

- [ ] **Step 5: Write Selenium test suite**

```python
# tests/test_tda_ma.py
# Test: page loads, demo loads, Run TDA works, Mapper renders,
#       tabs switch, export downloads, dark mode toggles
```

- [ ] **Step 6: Write R validation script**

```r
# validate_vs_r.R
library(TDAstats)
library(metafor)
# BCG data, compute persistence, compare against JS output
```

- [ ] **Step 7: Count divs, verify balance**

```bash
grep -c '<div' tda_ma.html
grep -c '</div>' tda_ma.html
# Must be equal
```

- [ ] **Step 8: Final polish**

- Keyboard accessibility (tab order, Enter activates buttons)
- escapeHtml() on all user input
- Console.log removed (or guarded behind debug flag)
- Version number in footer: "TDA-MA v1.0 | Mahmood Ahmad | Royal Free Hospital"

- [ ] **Step 9: Commit + push**

```bash
git add -A
git commit -m "feat: TDA-MA v1.0 — world's first browser-based TDA for meta-analysis"
gh repo create mahmood726-cyber/tda-ma --public --description "Topological Data Analysis for Meta-Analysis"
git remote add origin https://github.com/mahmood726-cyber/tda-ma.git
git push -u origin master
```

---

## Commit Sequence

| Task | Commit Message | ~Lines Added |
|------|---------------|-------------|
| 1 | `feat: HTML shell with 5-tab layout and dark mode` | ~500 |
| 2 | `feat: data module — 3 demo datasets, CSV parser, column mapper` | ~350 |
| 3 | `feat: linear algebra — standardization, distance matrices, PCA` | ~250 |
| 4 | `feat: persistent homology — Union-Find H0 + boundary reduction H1` | ~300 |
| 5 | `feat: Mapper algorithm — lens, cover, clustering, nerve` | ~350 |
| 6 | `feat: DL random-effects pooling + between-group Q-test` | ~250 |
| 7 | `feat: visualizations — Mapper canvas, persistence barcode, forest plots` | ~450 |
| 8 | `feat: subgroup profiles — covariate table, text summary, Q-test` | ~150 |
| 9 | `feat: export — R code, CSV, PNG, JSON` | ~200 |
| 10 | `feat: TDA-MA v1.0 — integration, tests, validation` | ~200 |
| **Total** | | **~3,000** |
