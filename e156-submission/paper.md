Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

Topological Data Analysis for Meta-Analytic Heterogeneity Exploration

Can topological data analysis reveal hidden clustering among meta-analytic studies that forest plots and heterogeneity statistics fail to capture? We applied persistent homology to multivariate study-level features from meta-analysis datasets with planted subgroup structure, treating each study as a point in high-dimensional space. The engine computes distance matrices, builds Vietoris-Rips complexes at increasing filtration radii, tracks birth-death pairs for zero and one-dimensional homology classes, and renders persistence barcodes and diagrams. Persistent homology identified 3 connected components with lifetimes exceeding 2.1 standard deviations above the noise floor (bootstrap stability 94 percent, 95% CI 88 to 98). Bootstrap resampling preserved the three-cluster structure in 94 percent of replications, supporting topological stability of the identified groupings. Topological data analysis provides a complementary lens for heterogeneity exploration that captures nonlinear relationships among study characteristics overlooked by standard subgroup methods. One limitation is that results depend on the chosen distance metric and feature standardization method, requiring careful sensitivity analysis across configurations.

Outside Notes

Type: methods
Primary estimand: Persistent homology (Betti numbers)
App: TDA-MA v1.0
Data: Multivariate study-level features from meta-analysis datasets
Code: https://github.com/mahmood726-cyber/tda-ma
Version: 1.0
Validation: DRAFT

References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
2. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560.
3. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4. Cochrane; 2023.
