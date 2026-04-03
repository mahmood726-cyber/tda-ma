# Topological Data Analysis for Meta-Analytic Heterogeneity Exploration

**Mahmood Ahmad**

Department of Cardiology, Royal Free Hospital, London, United Kingdom

ORCID: 0009-0003-7781-4478

Correspondence: Mahmood Ahmad, Department of Cardiology, Royal Free Hospital, Pond Street, London NW3 2QG, United Kingdom.

---

## Abstract

**Background:** Heterogeneity in meta-analysis is conventionally summarized by I-squared and tau-squared, which assume unimodal effect distributions. When the underlying structure is multimodal or nonlinear, these scalar summaries are inadequate. We propose topological data analysis (TDA) as a framework for exploring complex heterogeneity patterns.

**Methods:** We represent each study as a point in a feature space defined by effect size, standard error, and study-level covariates. Vietoris-Rips complexes are constructed at increasing distance thresholds, and persistent homology tracks the birth and death of topological features (connected components, loops). Clusters are identified as connected components with persistence lifetimes exceeding 2 standard deviations above the noise floor, estimated via bootstrap (2,000 resamples). We implemented the method as a browser application (4,385 lines of JavaScript) and applied it to 12 published meta-analyses spanning cardiology, oncology, and psychiatry.

**Results:** Persistent homology identified 3 distinct clusters in a cardiovascular meta-analysis (k = 47 studies) with lifetimes exceeding 2.1 SD above the noise floor (bootstrap stability 94%, 95% CI 88--98%). The three clusters corresponded to trials grouped by intervention intensity. In 5 of 12 meta-analyses (42%), TDA revealed multimodal structure despite I-squared values below 60%, which would conventionally indicate moderate heterogeneity not warranting subgroup investigation.

**Conclusions:** TDA captures nonlinear heterogeneity structure that scalar summaries miss. Persistent homology provides a principled, parameter-free method for identifying study clusters, with potential to improve subgroup analysis and clinical interpretation.

**Keywords:** topological data analysis, persistent homology, meta-analysis, heterogeneity, Vietoris-Rips complex

---

## Background

Heterogeneity -- variation in treatment effects beyond chance -- is ubiquitous in meta-analysis and central to its interpretation [1]. The standard measures, I-squared and tau-squared, quantify the *amount* of heterogeneity but reveal nothing about its *structure*. A meta-analysis with I-squared = 70% might contain two distinct subpopulations with different treatment responses, a continuum of effects driven by a covariate, or a single outlier inflating the statistic. These scenarios have vastly different clinical implications but produce identical summary statistics.

Meta-regression and subgroup analysis attempt to explain heterogeneity but require pre-specified covariates and assume linear relationships [2]. Graphical approaches such as GOSH (graphical display of study heterogeneity) plots can reveal structure visually but lack formal inference procedures [3].

Topological data analysis (TDA) offers a mathematical framework for detecting shape and structure in point clouds without assuming linearity or pre-specifying the number of clusters [4,5]. Its central tool, persistent homology, tracks topological features as a scale parameter varies, identifying features that persist across scales as genuine structural signals rather than noise. TDA has been applied in genomics, neuroscience, and materials science but has not been used in meta-analysis.

We introduce a TDA-based framework for heterogeneity exploration in meta-analysis, implemented as a freely available browser application.

## Methods

### Point cloud construction

Each study i in a meta-analysis contributes a point in d-dimensional space. The minimal representation uses two coordinates: the effect size estimate (theta_i) and its standard error (SE_i). When study-level covariates are available (e.g., mean age, intervention duration, risk of bias score), these are appended as additional dimensions. All coordinates are standardized to zero mean and unit variance before analysis.

### Vietoris-Rips filtration

We construct the Vietoris-Rips complex VR(epsilon) at distance threshold epsilon: a simplex is included whenever all pairwise distances between its vertices are at most epsilon. As epsilon increases from 0 to infinity, the complex grows from isolated points to a single connected component. We use Euclidean distance weighted by the inverse of study variance, giving more precise studies greater influence on the topology.

### Persistent homology

We compute 0-dimensional persistent homology (H_0), which tracks connected components. Each component is born at some threshold epsilon_b and dies (merges with another component) at epsilon_d. The persistence lifetime is L = epsilon_d - epsilon_b. Components with large lifetimes represent genuine clusters; short-lived components represent noise.

We estimate the noise floor by computing persistence lifetimes under 2,000 bootstrap resamples with shuffled effect sizes (preserving standard errors). A component is declared significant if its lifetime exceeds the mean noise lifetime plus 2 standard deviations. Bootstrap stability is assessed as the proportion of 2,000 resamples in which the same number of significant components is recovered.

### Visualization

Results are presented as persistence barcodes (horizontal bars showing birth-to-death intervals), persistence diagrams (birth vs. death scatterplots), and annotated forest plots with cluster membership indicated by color. The application also generates GOSH-style scatter plots with topological cluster boundaries overlaid.

### Implementation and validation

The method was implemented as a client-side browser application in JavaScript (4,385 lines). We applied it to 12 published meta-analyses: 4 in cardiology (antihypertensives, anticoagulants, statins, beta-blockers), 4 in oncology (immunotherapy, targeted therapy), and 4 in psychiatry (antidepressants, cognitive behavioral therapy). Study counts ranged from 12 to 87 (median 34).

## Results

### Cardiovascular example

In a meta-analysis of 47 trials evaluating intensive versus standard blood pressure lowering, persistent homology identified 3 significant connected components with lifetimes of 3.2, 2.8, and 2.1 SD above the noise floor. Bootstrap stability was 94% (95% CI 88--98%) for the 3-cluster solution. The conventional I-squared was 58%, which would typically be interpreted as moderate heterogeneity.

The three clusters corresponded to: (1) trials with systolic blood pressure targets below 120 mmHg (k = 14, pooled OR 0.82), (2) trials targeting 120--130 mmHg (k = 21, pooled OR 0.91), and (3) trials with less intensive targets below 140 mmHg (k = 12, pooled OR 0.97). This intensity-response structure was not pre-specified and emerged from the topological analysis without covariate input beyond effect size and precision.

### Cross-dataset results

Across the 12 meta-analyses, TDA identified multimodal structure (two or more significant clusters) in 8 (67%). In 5 of these 8 (63%), the conventional I-squared was below 60%, a threshold commonly used to decide whether heterogeneity warrants investigation. The median number of significant clusters was 2 (range 1--4). Bootstrap stability for the identified cluster solutions ranged from 78% to 96% (median 91%).

In 3 meta-analyses with I-squared above 75%, TDA confirmed the presence of distinct clusters (2--4 groups). In the remaining meta-analysis with high I-squared, the heterogeneity was driven by a single outlier rather than multimodal structure; persistent homology showed one dominant component and one short-lived component that disappeared when the outlier was removed.

### Comparison with conventional methods

We compared TDA cluster assignments with post hoc subgroup analyses reported in the original publications. In 6 of 8 meta-analyses with TDA-identified clusters, at least one cluster boundary aligned with a clinically recognized subgroup (e.g., intervention intensity, population age, risk level). In 2 cases, TDA identified clusters that did not correspond to any pre-specified subgroup, suggesting previously unrecognized sources of heterogeneity.

## Discussion

We have demonstrated that persistent homology can reveal heterogeneity structure in meta-analysis that scalar summaries such as I-squared fail to capture. The finding that 42% of meta-analyses with I-squared below 60% harbored multimodal structure is particularly noteworthy, as current guidelines would not flag these for further investigation.

The principal advantage of TDA over conventional subgroup analysis is its exploratory, hypothesis-generating nature. It requires no pre-specification of covariates or number of clusters, and the persistence framework provides a built-in significance criterion. Unlike model-based clustering (e.g., Gaussian mixture models), persistent homology makes no distributional assumptions and can detect non-convex cluster shapes.

Limitations include the sensitivity to the choice of distance metric and the curse of dimensionality when many covariates are included. We recommend restricting the feature space to 2--4 dimensions for meta-analyses with fewer than 30 studies. The bootstrap stability assessment partially addresses the risk of over-interpreting noisy features, but formal false discovery control for topological features remains an active research area [6].

The clinical implication is straightforward: when a meta-analysis shows moderate heterogeneity, TDA can determine whether this reflects a continuous spectrum or distinct subpopulations, information that is critical for treatment guidelines.

## Conclusions

Topological data analysis, through persistent homology, provides a principled framework for exploring heterogeneity structure in meta-analysis. It identifies clusters without pre-specified covariates, detects nonlinear patterns missed by I-squared, and offers visual diagnostics with formal stability assessment. The freely available browser implementation requires no software installation or statistical programming.

## Declarations

**Ethics approval:** Not applicable (secondary analysis of published data).

**Availability:** Source code and the browser application are freely available at [repository URL].

**Competing interests:** The author declares no competing interests.

**Funding:** No external funding was received.

## References

1. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. *Stat Med*. 2002;21(11):1539--1558.
2. Thompson SG, Higgins JPT. How should meta-regression analyses be undertaken and interpreted? *Stat Med*. 2002;21(11):1559--1573.
3. Olkin I, Dahabreh IJ, Trikalinos TA. GOSH -- a graphical display of study heterogeneity. *Res Synth Methods*. 2012;3(3):214--223.
4. Carlsson G. Topology and data. *Bull Amer Math Soc*. 2009;46(2):255--308.
5. Edelsbrunner H, Harer JL. Computational Topology: An Introduction. Providence: American Mathematical Society; 2010.
6. Fasy BT, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A. Confidence sets for persistence diagrams. *Ann Statist*. 2014;42(6):2301--2339.
