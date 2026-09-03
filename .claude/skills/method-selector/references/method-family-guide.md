# Method-Family Routing Guide

Use this guide to form a shortlist, not as a menu that must be exhausted.

## Evaluation and ranking

- Credible baselines: equal-weight normalized score or an existing operational rule, when they complete the real task.
- Main families: entropy/TOPSIS, PCA-assisted evaluation, grey/fuzzy evaluation when their assumptions fit.
- Key risks: redundant indicators, unjustified directions or weights, weight dominance, concentrated scores, unstable top-k ranks.

## Regression and relationship modeling

- Credible baselines: pooled OLS, simple correlation analysis.
- Main families (linear): multiple regression, mixed-effects models (LME/GLMM) for clustered/panel data.
- Main families (nonlinear): GAM/GAMM (spline-based, data-driven nonlinearity), polynomial regression, piecewise/segmented regression, quantile regression.
- Main families (high-dimensional): LASSO/elastic net, random forest, gradient boosting regression.
- Prefer mixed-effects models (LME/GLMM/GAMM) when data has repeated measures, clustering, or hierarchical structure — ICC > 0.1 is a strong signal.
- Prefer GAM/GAMM over linear models when nonlinearity is plausible — let the data decide the functional form rather than assuming linearity.
- Key risks: multicollinearity (check VIF), residual non-normality, overdispersion, boundary effects for bounded responses, ignoring clustering inflates significance.

## Prediction

- Credible baselines: seasonal naive, last value, moving average, or simple regression selected to match the time structure.
- Main families: exponential smoothing, ARIMA/SARIMA, regression, tree boosting, small-data grey models.
- Key risks: leakage, invalid split, short series, nonstationarity, over-capacity, poor interval coverage.

## Optimization

- Credible baselines: current policy, a feasible greedy rule, or a relaxed exact formulation.
- Main families: LP/MILP, network flow, dynamic or nonlinear programming, justified metaheuristics.
- Key risks: missing constraints, infeasibility, meaningless objectives, nonimplementable solutions, excessive runtime.

## Classification and clustering

- Credible baselines: rule-based or majority/stratified reference when meaningful; Z-score threshold rules for anomaly detection.
- Main families (supervised): logistic/tree models, SVM, random forest, XGBoost, neural networks.
- Main families (unsupervised/anomaly): Mahalanobis distance, Isolation Forest, LOF, one-class SVM, DBSCAN, Gaussian mixture models.
- Main families (clustering): k-means, hierarchical clustering, spectral clustering, DBSCAN.
- When labeled positive samples are absent or extremely rare, prefer unsupervised anomaly detection over supervised classification.
- For multivariate anomaly detection with correlated features, prefer Mahalanobis distance or Isolation Forest over univariate thresholds.
- Key risks: label absence, class imbalance, arbitrary cluster count, instability, accuracy-only evaluation, threshold sensitivity in anomaly detection.

## Mechanism and simulation

- Credible baselines: simplified algebraic or deterministic scenario model that preserves the core mechanism.
- Main families: differential/difference equations, compartment models, Monte Carlo, discrete-event or agent-based simulation.
- Key risks: unidentified parameters, unit errors, invalid boundary conditions, hidden distribution assumptions, too few replications.

## Graph and routing

- Credible baselines: feasible direct rule or nearest-neighbor heuristic.
- Main families: shortest path, flow, matching, spanning tree, TSP/VRP formulations.
- Key risks: unrealistic edges or weights, omitted operational constraints, mathematically valid but unusable routes.
