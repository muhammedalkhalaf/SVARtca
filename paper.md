---
title: 'SVARtca: Transmission Channel Analysis for Structural VAR Models in R'
tags:
  - R
  - econometrics
  - structural VAR
  - impulse response functions
  - transmission channels
  - directed acyclic graphs
  - macroeconomics
authors:
  - name: Muhammad Abdullah Alkhalaf
    orcid: 0009-0002-2677-9246
    corresponding: true
    email: muhammedalkhalaf@gmail.com
    affiliation: 1
affiliations:
  - name: Rufyq Elngeh for Academic and Business Services, Riyadh, Saudi Arabia
    index: 1
date: 18 March 2026
bibliography: paper.bib
---

# Summary

`SVARtca` is an R package that implements Transmission Channel Analysis (TCA) for structural vector autoregressive (SVAR) models, following the methodology of @Wegner2025. TCA decomposes impulse response functions (IRFs) into contributions from distinct causal transmission channels using a systems form representation and directed acyclic graph (DAG) path analysis. The package supports overlapping channels, exhaustive 3-way and 4-way decompositions via the inclusion-exclusion principle, and provides publication-ready visualizations. `SVARtca` is a parallel R implementation of the original MATLAB toolbox by Wegner, Lieb, and Smeekes, making this methodology accessible to the large R-based macroeconometrics community.

# Statement of Need

Structural VAR models are the workhorse of empirical macroeconomics for analyzing the dynamic effects of structural shocks. Standard impulse response analysis quantifies the total effect of a shock on each variable over time, but does not reveal *through which channels* these effects propagate. Understanding transmission channels is crucial for policy analysis—for example, whether monetary policy affects output primarily through the interest rate channel, the credit channel, or the exchange rate channel.

@Wegner2025 introduced TCA as a formal framework for decomposing IRFs into channel-specific contributions. The methodology represents the SVAR in systems form and uses DAG-based path analysis to isolate direct, indirect, and channel-specific effects. When channels overlap (share intermediate variables), the inclusion-exclusion principle ensures proper attribution without double-counting.

The original implementation is available as a MATLAB toolbox (`tca-matlab-toolbox`), which limits accessibility for researchers who work primarily in R. The `vars` package [@Pfaff2008vars] provides standard SVAR estimation and IRF computation in R, but does not implement channel decomposition. `SVARtca` bridges this gap by providing a native R implementation that integrates naturally with existing R-based VAR workflows.

# Features

## Channel Decomposition

The core function decomposes IRFs into contributions from user-specified transmission channels. Each channel is defined by a set of intermediate variables through which the structural shock propagates from the impulse variable to the response variable.

```r
library(SVARtca)

# Define a 4-variable SVAR: monetary policy -> output
# Variables: interest_rate, credit, exchange_rate, output
# Channels: credit channel and exchange rate channel

# Estimate SVAR (using vars package or custom estimation)
# svar_model <- ...

# Define channels
channels <- list(
  credit = c("credit"),
  exchange_rate = c("exchange_rate")
)

# Perform TCA
tca_result <- tca(svar_model,
                  impulse = "interest_rate",
                  response = "output",
                  channels = channels,
                  horizon = 20)
summary(tca_result)
plot(tca_result)
```

## Overlapping Channels

When transmission channels share intermediate variables, `SVARtca` applies the inclusion-exclusion principle to produce an exhaustive decomposition. For $k$ potentially overlapping channels, the total IRF is decomposed into $2^k - 1$ non-overlapping components plus a direct effect.

```r
# Channels that share the "expectations" variable
channels <- list(
  interest_rate = c("real_rate", "expectations"),
  wealth = c("asset_prices", "expectations")
)

tca_result <- tca(svar_model,
                  impulse = "monetary_policy",
                  response = "consumption",
                  channels = channels,
                  horizon = 20,
                  overlap = "inclusion_exclusion")
```

## Exhaustive Decompositions

For systems with 3 or 4 channels, the package provides exhaustive decomposition functions that compute all possible channel combinations and their interactions, ensuring that the sum of all components equals the total IRF at each horizon.

```r
# 3-way exhaustive decomposition
decomp <- tca_exhaustive(svar_model,
                         impulse = "oil_shock",
                         response = "gdp",
                         channels = list(
                           supply = c("production"),
                           demand = c("consumption"),
                           financial = c("stock_market")
                         ),
                         horizon = 24)
plot(decomp, type = "stacked")
```

## Visualization

The package provides `ggplot2`-based visualization functions for channel decomposition results, including stacked area plots, grouped bar charts, and individual channel IRF plots with confidence bands.

# Implementation

`SVARtca` is implemented in pure R with dependencies on `Matrix` for efficient sparse matrix operations, `ggplot2` for visualization, and `rlang` for non-standard evaluation. The package follows the systems form representation of the SVAR:

$$\mathbf{y}_t = \sum_{j=1}^{p} \mathbf{\Phi}_j \mathbf{y}_{t-j} + \mathbf{B} \boldsymbol{\varepsilon}_t$$

where the structural moving average representation $\mathbf{y}_t = \sum_{s=0}^{\infty} \mathbf{\Theta}_s \boldsymbol{\varepsilon}_{t-s}$ is decomposed via the DAG adjacency structure implied by the contemporaneous impact matrix $\mathbf{B}$.

The package is compatible with SVAR objects from the `vars` package and accepts custom coefficient matrices for maximum flexibility.

# Acknowledgements

The author acknowledges Emanuel Wegner, Lenard Lieb, and Stephan Smeekes for developing the TCA methodology and the original MATLAB toolbox. This R implementation was developed independently as a parallel port to make TCA accessible to the R community.

# References
