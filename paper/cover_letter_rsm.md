[Date]

The Editors
Research Synthesis Methods

Dear Editors,

We are pleased to submit our manuscript entitled "Topological Subgroup Discovery in Meta-Analysis: A Browser-Based Tool Using Persistent Homology and the Mapper Algorithm" for consideration in Research Synthesis Methods.

**Summary.** We present TDA-MA, the first tool that applies topological data analysis (TDA) -- specifically persistent homology and the Mapper algorithm -- to meta-analytic datasets for unsupervised subgroup discovery. Unlike conventional subgroup analysis and meta-regression, which require the analyst to pre-specify which study-level variable to stratify by, TDA discovers structure directly from the geometry of the data without prior hypotheses about moderator variables.

**Why this matters.** Heterogeneity remains the central analytic challenge in meta-analysis. When the true moderating structure is multivariate or unexpected, standard approaches will miss it. TDA addresses this gap by providing a mathematically principled framework for detecting clusters, loops, and other topological features in study-level data -- revealing structure the analyst did not know to look for.

**Validation.** We validated TDA-MA on the BCG vaccine dataset (Colditz et al. 1994), where it separated 13 trials into high-latitude and low-latitude clusters without being given latitude as an input variable, recovering the well-established geographic modifier with a highly significant between-cluster test (Q = 47.90, p = 3.7 x 10^-8). On a simulated dataset with planted clusters, TDA-MA achieved 95.2% cluster recovery accuracy, exceeding the pre-specified 90% threshold. DerSimonian-Laird pooled estimates matched the R package metafor within tolerance.

**Implementation.** TDA-MA is a single HTML file (3,619 lines) with all computation in pure JavaScript. It requires no installation, no server, and no external dependencies -- users simply open the file in any modern browser. The tool includes interactive Mapper graph visualisation, persistence barcodes, per-cluster forest plots, between-subgroup Q-tests, covariate profiling, and export to R code, CSV, and PNG. It passes 66 automated tests and is freely available at https://github.com/mahmood726-cyber/tda-ma.

**Fit with Research Synthesis Methods.** This manuscript introduces a genuinely new methodological approach to meta-analysis from a well-established mathematical discipline (algebraic topology) that has not previously been applied to evidence synthesis. It directly addresses the journal's focus on innovative methods for research synthesis and is accompanied by a validated, openly available implementation.

We confirm that this manuscript has not been published elsewhere and is not under consideration at another journal. The author declares no competing interests.

We look forward to your consideration.

Sincerely,

Mahmood Ahmad
Royal Free Hospital, London, UK
mahmood.ahmad2@nhs.net
ORCID: 0009-0003-7781-4478
