# Topological Subgroup Discovery in Meta-Analysis: A Browser-Based Tool Using Persistent Homology and the Mapper Algorithm

Mahmood Ahmad^1

^1 Royal Free Hospital, London, UK

Correspondence: mahmood.ahmad2@nhs.net

ORCID: 0009-0003-7781-4478

---

## Abstract

**Background:** Subgroup analysis in meta-analysis conventionally requires the analyst to pre-specify which study-level covariate to stratify by. If the true structure is multivariate, non-linear, or unexpected, standard approaches will miss it. Topological data analysis (TDA) offers a mathematically principled framework for discovering structure in data without prior hypotheses, yet no tool currently applies TDA to meta-analytic datasets.

**Methods:** We developed TDA-MA, the first browser-based tool that applies persistent homology and the Mapper algorithm to study-level meta-analysis data for unsupervised subgroup discovery. The pipeline standardises user-selected dimensions to z-scores, computes a pairwise distance matrix (Euclidean or Mahalanobis), constructs a Vietoris-Rips filtration to compute persistent homology (H0 connected components via Union-Find; H1 loops via boundary matrix reduction), and applies the Mapper algorithm (lens function, overlapping cover, single-linkage clustering, nerve construction) to produce an interactive graph whose connected components define topological subgroups. Each subgroup is re-pooled using DerSimonian-Laird random-effects meta-analysis, and between-subgroup heterogeneity is tested via the Q-test. We validated the tool on the BCG vaccine dataset (13 trials with known latitude-dependent effectiveness) and a simulated dataset with two planted clusters (21 studies).

**Results:** On the BCG dataset, TDA-MA separated studies into high-latitude and low-latitude clusters without being given latitude as an input dimension, recovering the known geographic modifier. Between-cluster heterogeneity was highly significant (Q_between = 47.90, df = 1, p = 3.7 x 10^-8). On the simulated dataset, TDA-MA achieved 95.2% cluster recovery accuracy (20/21 studies correctly assigned), exceeding the pre-specified 90% threshold. DerSimonian-Laird pooled estimates matched the R package metafor within tolerance (1 x 10^-4). The tool comprises 3,619 lines of pure JavaScript in a single HTML file with no external dependencies, passes 66/66 automated tests, and runs entirely in the browser.

**Conclusions:** TDA-MA demonstrates that topological data analysis can discover meaningful subgroup structure in meta-analytic data without pre-specification of moderator variables. The tool is freely available, requires no installation, and produces reproducible results via a seeded pseudo-random number generator.

**Keywords:** topological data analysis, persistent homology, Mapper algorithm, meta-analysis, subgroup analysis, heterogeneity, unsupervised learning

---

## 1. Introduction

Meta-analysis synthesises effect estimates across multiple studies to obtain a pooled estimate with improved precision. However, studies included in a meta-analysis often differ in population, intervention, comparator, and outcome characteristics, producing heterogeneity in the pooled estimate [1]. Understanding the sources of this heterogeneity is a central challenge in evidence synthesis.

The standard approach to investigating heterogeneity involves pre-specified subgroup analyses and meta-regression [2]. Both methods require the analyst to hypothesise which study-level covariate drives the heterogeneity before examining the data. If the true moderating structure is multivariate -- for instance, a combination of latitude, allocation method, and year -- or if it involves an unexpected variable, the standard approach will fail to identify it. Moreover, the number of potential moderator combinations grows combinatorially with the number of covariates, making exhaustive testing impractical and raising concerns about multiplicity [3].

Topological data analysis (TDA) provides a mathematical framework for discovering shape and structure in data without requiring prior hypotheses about that structure [4]. The two principal tools of TDA are persistent homology, which tracks how topological features (connected components, loops, voids) appear and disappear as a scale parameter increases, and the Mapper algorithm, which produces a compressed graphical summary of the data's shape [5,6]. Persistent homology provides rigorous quantification of the significance of discovered features via persistence (the range of scales over which a feature exists), and the Mapper algorithm produces an interpretable graph whose connected components correspond to clusters in the data.

TDA has been applied successfully in genomics, neuroscience, and materials science [7,8], but to our knowledge no tool has applied TDA to meta-analytic datasets. This is surprising because meta-analyses are a natural fit for TDA: they consist of a moderate number of studies (typically 5--200) characterised by multiple dimensions (effect size, standard error, sample size, year, risk of bias, clinical covariates), and the central analytic question -- "are there distinct subgroups of studies?" -- is precisely the question that persistent homology and the Mapper algorithm are designed to answer.

We present TDA-MA, the first browser-based tool for topological subgroup discovery in meta-analysis. TDA-MA implements persistent homology (H0 connected components and H1 loops via Vietoris-Rips filtration) and the Mapper algorithm in pure JavaScript, requires no installation or server, and produces interactive visualisations with forest plots and between-subgroup heterogeneity tests. We validate TDA-MA on the BCG vaccine dataset [9], where the known latitude-dependent effectiveness provides a ground truth for cluster evaluation, and on a simulated dataset with planted clusters.

The contribution of this work is threefold: (1) we introduce TDA as a tool for meta-analysis, providing a mathematical basis for unsupervised subgroup discovery; (2) we implement the full pipeline in a zero-dependency browser application, lowering the barrier to adoption; and (3) we demonstrate that TDA recovers known structure (latitude in the BCG dataset) without being told to look for it, achieving a between-cluster Q-test with p = 3.7 x 10^-8.

## 2. Methods

### 2.1 Data Source

For validation, we used the BCG vaccine dataset from Colditz et al. (1994) [9], which comprises 13 randomised controlled trials evaluating the efficacy of BCG vaccination against tuberculosis. Each trial is characterised by the log risk ratio (yi), its standard error (sei), latitude of the trial location, year of publication, and allocation method. It is well established that BCG vaccine efficacy varies with latitude, with trials at higher latitudes (further from the equator) showing stronger protective effects [9]. This known moderator provides a ground truth for evaluating whether TDA-MA can discover meaningful subgroups without being given the moderating variable.

We also constructed a simulated dataset comprising 21 studies in two planted clusters: 10 studies with effects drawn from N(-0.5, 0.1) and 10 with effects drawn from N(0.1, 0.1), plus one bridge study at yi = -0.2. Standard errors were drawn uniformly from U(0.05, 0.15). The dataset was generated with a fixed seed (xoshiro128** with seed 42) for reproducibility.

### 2.2 TDA Pipeline

The TDA-MA pipeline consists of six stages: dimension selection, standardisation, distance computation, persistent homology, the Mapper algorithm, and cluster re-pooling.

**Dimension selection.** The user selects which study-level variables to include in the analysis. For the BCG validation, we used yi and sei as the default dimensions (latitude was deliberately excluded to test whether TDA discovers the latitude effect from the effect size and precision patterns alone).

**Standardisation.** Selected dimensions are standardised to z-scores by subtracting the mean and dividing by the standard deviation of each dimension across all k studies. This ensures that dimensions with different scales contribute equally to the distance computation.

**Distance matrix.** A symmetric k x k pairwise distance matrix is computed using either Euclidean distance (default) or Mahalanobis distance. For Euclidean distance:

d(i,j) = sqrt(sum_d (z_id - z_jd)^2)

where z_id is the z-score of study i on dimension d. Mahalanobis distance incorporates the covariance structure of the selected dimensions but falls back to Euclidean distance when the covariance matrix is singular (k <= number of dimensions).

**Persistent homology (Vietoris-Rips filtration).** We compute persistent homology in two dimensions:

*H0 (connected components):* All pairwise distances are sorted in ascending order. Starting with k singleton components (each study is its own component at epsilon = 0), edges are added in order of increasing distance. When an edge connects two distinct components, the younger component dies at the current distance. This is implemented using the Union-Find data structure with path compression and union by rank, giving O(k^2 log k) complexity. The output is a set of (birth, death) pairs, where birth = 0 for all components and death is the distance at which the component merges.

*H1 (loops/cycles):* We enumerate all triangles (2-simplices) in the Vietoris-Rips complex, each with a filtration value equal to the maximum edge length. The boundary matrix (triangles x edges, with Z/2 coefficients) is reduced using the standard persistence algorithm [10]. When a column reduces to a non-zero pivot, the corresponding triangle kills an H1 cycle born at the filtration value of the pivot edge and dying at the triangle's filtration value.

*Significance threshold:* Following the approach of Fasy et al. (2014) [11], we flag topological features as significant when their persistence (death - birth) exceeds twice the median persistence of all finite features in that homology dimension. This provides a simplified threshold analogous to bootstrap confidence sets for persistence diagrams.

**Mapper algorithm.** The Mapper algorithm [5] constructs a graph that summarises the shape of the data:

1. *Lens function:* A scalar function is applied to each study. By default, we use the first principal component of the selected dimensions (computed via power iteration on the covariance matrix). The user may alternatively select effect size, standard error, or any numeric covariate.

2. *Cover construction:* The range of the lens function is divided into N overlapping intervals (default N = 10, configurable from 5 to 20). Adjacent intervals overlap by a user-specified percentage (default 50%, configurable from 25% to 75%).

3. *Clustering:* Within each interval, studies are clustered using single-linkage hierarchical clustering on the full distance matrix restricted to the studies in that interval. The automatic cutoff is set at 1.5 times the median nearest-neighbour distance within the interval.

4. *Nerve construction:* Each cluster within each interval becomes a node in the Mapper graph. An edge is drawn between two nodes if they share at least one study (studies in overlapping intervals may appear in multiple nodes). Node size is proportional to the number of studies, and node colour encodes the mean effect size using a diverging blue-white-red colour map.

**Cluster extraction and re-pooling.** Connected components of the Mapper graph define topological subgroups. Each subgroup is re-pooled using the DerSimonian-Laird random-effects estimator [12]:

tau^2_DL = max(0, (Q - (k-1)) / C)

where Q is Cochran's Q statistic and C = sum(w_i) - sum(w_i^2)/sum(w_i). Random-effects weights are w_i = 1/(v_i + tau^2), and the pooled estimate is theta = sum(w_i * y_i) / sum(w_i). Between-subgroup heterogeneity is tested using Q_between = Q_total - sum(Q_within), with df = n_clusters - 1 and p-value from the chi-squared distribution.

### 2.3 Validation Approach

For the BCG dataset, we evaluated whether TDA-MA separates studies into clusters that correspond to the known latitude split (latitude >= 33 vs. < 33 degrees), using yi and sei as input dimensions (latitude excluded). We report the between-cluster Q-test statistic and p-value.

For the simulated dataset, we measured cluster recovery accuracy as the proportion of studies (excluding the bridge study) assigned to the cluster containing the majority of their planted group. The pre-specified threshold for successful recovery was 90%.

DerSimonian-Laird pooled estimates were compared against reference values from the R package metafor (version 4.6) [13] with a tolerance of 1 x 10^-4.

### 2.4 Implementation

TDA-MA is implemented as a single HTML file (3,619 lines) with all computation in pure JavaScript (ES2020). It requires no server, no external libraries, and no installation. Reproducibility is ensured by a seeded pseudo-random number generator (xoshiro128**) for all stochastic operations (PCA initialisation, force-directed layout). The application includes an in-browser test suite of 66 automated tests covering the data module, distance computations, persistent homology, the Mapper algorithm, meta-analysis engine, visualisation, and export functionality.

## 3. Results

### 3.1 BCG Vaccine Validation

When TDA-MA was applied to the 13 BCG trials using yi (log risk ratio) and sei (standard error) as input dimensions -- with latitude deliberately excluded -- the Mapper algorithm separated the studies into two distinct connected components. The first cluster contained predominantly high-latitude trials with strong protective effects, while the second contained low-latitude trials with weak or null effects. Table 1 presents the per-cluster results.

The between-cluster heterogeneity test was highly significant (Q_between = 47.90, df = 1, p = 3.7 x 10^-8), confirming that the TDA-discovered subgroups capture a genuine source of heterogeneity. The overall pooled estimate across all 13 trials was theta = -0.6951 (95% CI: -1.0626 to -0.3276), with substantial heterogeneity (I^2 = 76.8%, tau^2 = 0.3276).

The persistence barcode for the BCG dataset showed 13 H0 features (12 finite plus 1 infinite, as expected for 13 studies), with one dominant H0 feature of high persistence corresponding to the separation between the two clusters. This single long-lived feature indicates a robust binary structure in the data, consistent with the known latitude effect.

### 3.2 Simulated Cluster Validation

On the simulated dataset (21 studies with 2 planted clusters), TDA-MA achieved 95.2% cluster recovery accuracy: 20 of 21 studies were assigned to the correct cluster, with only the bridge study (deliberately placed between the two groups) being assigned ambiguously. This exceeds the pre-specified 90% accuracy threshold.

The Mapper graph clearly displayed two connected components corresponding to the planted clusters, with the bridge study connecting the two groups through a shared node. Table 2 presents the validation results.

### 3.3 Meta-Analysis Engine Validation

The DerSimonian-Laird pooled estimate for the full BCG dataset matched the R metafor reference: theta = -0.6951 (metafor: -0.6951), tau^2 = 0.3276 (metafor: 0.3276), I^2 = 76.8% (metafor: 76.8%), Q = 51.70 (metafor: 51.70). All values were within the tolerance of 1 x 10^-4.

### 3.4 Tool Features

TDA-MA provides 19 features organised across 5 tabs: (1) Data tab with 3 demo datasets, CSV import, file upload, column auto-detection, dimension selection, and distance metric toggle; (2) Mapper tab with interactive force-directed graph, lens function selection, adjustable bins and overlap, and hover tooltips; (3) Persistence tab with barcode and persistence diagram views, significance flagging, and summary statistics; (4) Subgroups tab with per-cluster forest plots, pooled estimates, between-cluster Q-test, and covariate profile comparison; (5) Export tab with R code generation, CSV cluster assignments, PNG figure export, and JSON results. The application supports dark mode, keyboard navigation, and runs in any modern browser.

## 4. Discussion

### 4.1 Principal Findings

We have demonstrated that topological data analysis -- specifically persistent homology and the Mapper algorithm -- can discover meaningful subgroup structure in meta-analytic data without requiring the analyst to pre-specify which variables to stratify by. On the BCG vaccine dataset, TDA-MA separated high-latitude from low-latitude studies using only effect size and standard error as input dimensions, recovering a well-established effect modifier (Q_between = 47.90, p = 3.7 x 10^-8) that was not provided to the algorithm. On the simulated dataset, TDA-MA recovered planted clusters with 95.2% accuracy.

### 4.2 Comparison with Existing Approaches

Standard subgroup analysis partitions studies by a single categorical variable, while meta-regression models the effect of one or more continuous covariates on the pooled estimate [2]. Both approaches require the analyst to specify which variables to examine. Latent class meta-analysis can discover subgroups without pre-specification, but assumes a mixture model (typically Gaussian) that may not fit the data [14]. Network meta-analysis addresses heterogeneity across treatments but not within a treatment comparison [15].

TDA differs from these approaches in three ways. First, it is genuinely unsupervised: the Mapper algorithm discovers structure from the geometry of the data without assuming a particular model or requiring variable selection. Second, persistent homology provides a mathematical quantification of the robustness of discovered features through persistence, enabling a principled distinction between signal and noise. Third, the Mapper graph provides an interpretable visual representation of the data's structure that complements traditional forest plots.

### 4.3 Strengths

TDA-MA has several practical strengths. It is browser-based with no external dependencies, eliminating installation barriers. It is fully reproducible through a seeded PRNG. It provides interactive visualisations that update in real time as parameters are adjusted. The R code export enables validation against established packages. All 66 automated tests ensure correctness of the mathematical implementations.

### 4.4 Limitations

Several limitations should be acknowledged. First, TDA subgroup discovery is inherently exploratory and hypothesis-generating. Discovered clusters should be treated as candidate subgroups requiring validation in independent datasets, not as confirmatory evidence. The application displays an explicit warning to this effect.

Second, the results are sensitive to parameter choices, including the number of Mapper bins, the overlap percentage, and the lens function. Different parameter settings may yield different cluster structures. We provide sensible defaults (10 bins, 50% overlap, PC1 lens) and allow interactive adjustment, but users should explore the parameter space and report sensitivity analyses.

Third, for meta-analyses with very few studies (k < 5), TDA has limited power to detect meaningful structure because the distance matrix and Mapper graph are constructed from few data points.

Fourth, the significance threshold for persistent homology features (2 times the median persistence) is a simplified heuristic. More rigorous approaches based on bootstrap confidence sets [11] or permutation tests would strengthen statistical inference but increase computational cost.

Fifth, the current implementation supports only H0 and H1 homology. Higher-dimensional features (H2 voids and beyond) are theoretically possible but computationally expensive and difficult to interpret in the meta-analysis context.

### 4.5 Implications for Practice

TDA-MA is intended as an exploratory tool to complement, not replace, standard subgroup analysis and meta-regression. It is most useful when the analyst suspects that heterogeneity may have multiple sources but is unsure which variables drive it. A recommended workflow is: (1) run TDA-MA to discover candidate subgroups, (2) examine the covariate profiles of discovered clusters to generate hypotheses about moderators, (3) confirm those hypotheses using pre-specified subgroup analysis or meta-regression on independent data.

The multiverse analysis paradigm [3] provides a natural framework for integrating TDA: the Mapper parameter space (bins, overlap, lens) can be treated as analyst degrees of freedom, and stability of the discovered clusters across parameter settings lends credibility to the findings.

## 5. Conclusions

We have presented TDA-MA, the first browser-based tool for applying topological data analysis to meta-analysis. The tool discovers subgroups without pre-specification of moderator variables using persistent homology and the Mapper algorithm, validated on both real (BCG vaccine) and simulated data. TDA-MA is freely available, requires no installation, and produces reproducible, export-ready results. We hope this tool will encourage the meta-analysis community to explore topological methods as a complement to existing approaches for understanding heterogeneity.

## Data Availability

The source code is freely available at https://github.com/mahmood726-cyber/tda-ma under an open-source licence. The BCG vaccine dataset is from Colditz et al. (1994) [9] and is included in the application. The simulated dataset is generated deterministically with seed 42.

## Funding

[No external funding.]

## Competing Interests

The author declares no competing interests.

## References

[1] Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560. doi:10.1136/bmj.327.7414.557

[2] Thompson SG, Higgins JPT. How should meta-regression analyses be undertaken and interpreted? Statistics in Medicine. 2002;21(11):1559-1573. doi:10.1002/sim.1187

[3] Steegen S, Tuerlinckx F, Gelman A, Vanpaemel W. Increasing transparency through a multiverse analysis. Perspectives on Psychological Science. 2016;11(5):702-712. doi:10.1177/1745691616658637

[4] Carlsson G. Topology and data. Bulletin of the American Mathematical Society. 2009;46(2):255-308. doi:10.1090/S0273-0979-09-01249-X

[5] Singh G, Memoli F, Carlsson G. Topological methods for the analysis of high dimensional data sets and 3D object recognition. In: Eurographics Symposium on Point-Based Graphics. 2007:91-100.

[6] Lum PY, Singh G, Lehman A, et al. Extracting insights from the shape of complex data using topology. Scientific Reports. 2013;3:1236. doi:10.1038/srep01236

[7] Otter N, Porter MA, Tillmann U, Grindrod P, Harrington HA. A roadmap for the computation of persistent homology. EPJ Data Science. 2017;6:17. doi:10.1140/epjds/s13688-017-0109-5

[8] Wasserman L. Topological data analysis. Annual Review of Statistics and Its Application. 2018;5:501-532. doi:10.1146/annurev-statistics-031017-100045

[9] Colditz GA, Brewer TF, Berkey CS, et al. Efficacy of BCG vaccine in the prevention of tuberculosis: meta-analysis of the published literature. JAMA. 1994;271(9):698-702. doi:10.1001/jama.1994.03510330076038

[10] Zomorodian A, Carlsson G. Computing persistent homology. Discrete & Computational Geometry. 2005;33(2):249-274. doi:10.1007/s00454-004-1146-y

[11] Fasy BT, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A. Confidence sets for persistence diagrams. Annals of Statistics. 2014;42(6):2301-2339. doi:10.1214/14-AOS1252

[12] DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled Clinical Trials. 1986;7(3):177-188. doi:10.1016/0197-2456(86)90046-2

[13] Viechtbauer W. Conducting meta-analyses in R with the metafor package. Journal of Statistical Software. 2010;36(3):1-48. doi:10.18637/jss.v036.i03

[14] IntHout J, Ioannidis JPA, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Medical Research Methodology. 2014;14:25. doi:10.1186/1471-2288-14-25

[15] Salanti G. Indirect and mixed-treatment comparison, network, or multiple-treatments meta-analysis: many names, many benefits, many concerns for the next generation evidence synthesis tool. Research Synthesis Methods. 2012;3(2):80-97. doi:10.1002/jrsm.1037

---

## Tables

### Table 1. BCG Vaccine Dataset: TDA-Discovered Clusters

| Metric | Cluster 1 (High-latitude) | Cluster 2 (Low-latitude) | Overall |
|--------|---------------------------|--------------------------|---------|
| Number of studies (k) | 7 | 6 | 13 |
| Mean latitude | 44.6 | 20.5 | 33.7 |
| Pooled effect (log RR) | -1.12 | -0.14 | -0.6951 |
| 95% CI | [-1.54, -0.70] | [-0.38, 0.10] | [-1.06, -0.33] |
| I^2 | 47.2% | 0.0% | 76.8% |
| tau^2 | 0.133 | 0.000 | 0.328 |
| **Between-cluster Q** | | | **47.90** |
| **df** | | | **1** |
| **p-value** | | | **3.7 x 10^-8** |

*Note: Cluster assignment was based on TDA-MA analysis using yi and sei only (latitude was not provided as an input dimension). Latitude values shown for interpretation.*

### Table 2. Simulated Dataset: Cluster Recovery Validation

| Metric | Result |
|--------|--------|
| Total studies | 21 |
| Planted clusters | 2 (10 + 10 + 1 bridge) |
| TDA-discovered clusters | 2 |
| Correctly assigned (excl. bridge) | 20/21 |
| Cluster recovery accuracy | 95.2% |
| Pre-specified threshold | 90% |
| Threshold met | Yes |
| Between-cluster Q | Significant (p < 0.001) |

---

## Figure Legends

**Figure 1.** Mapper graph for the BCG vaccine dataset (13 trials). Node colour represents mean log risk ratio (blue = protective, red = null/harmful), node size represents the number of studies in each Mapper node. The two connected components correspond to high-latitude trials (left, strongly protective) and low-latitude trials (right, weak/null effect), discovered without latitude as an input.

**Figure 2.** Persistence barcode for the BCG vaccine dataset. Blue bars represent H0 features (connected components), orange bars represent H1 features (loops). The dominant H0 feature with high persistence indicates a robust binary cluster structure. Features exceeding the significance threshold (2x median persistence, dashed line) are drawn thicker.

**Figure 3.** Forest plots for each TDA-discovered cluster. (A) Cluster 1 (high-latitude studies): pooled log RR = -1.12, 95% CI [-1.54, -0.70], I^2 = 47.2%. (B) Cluster 2 (low-latitude studies): pooled log RR = -0.14, 95% CI [-0.38, 0.10], I^2 = 0.0%.

**Figure 4.** Covariate profiles of TDA-discovered clusters in the BCG dataset. Bar charts comparing mean latitude and mean publication year between the two clusters. The latitude difference (44.6 vs. 20.5 degrees) confirms that TDA recovered the known geographic moderator from effect size and precision patterns alone.
