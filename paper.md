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

`SVARtca` is an R package implementing Transmission Channel Analysis (TCA) for structural VAR models following @Wegner2025. TCA decomposes impulse response functions into contributions from distinct causal transmission channels using a systems form representation and DAG path analysis. The package supports overlapping channels, exhaustive 3-way and 4-way decompositions via the inclusion-exclusion principle, and publication-ready `ggplot2` visualizations. It is a parallel R implementation of the original MATLAB toolbox by Wegner, Lieb, and Smeekes.

# Statement of Need

Standard SVAR impulse response analysis quantifies total shock effects but does not reveal *through which channels* effects propagate. @Wegner2025 introduced TCA to decompose IRFs into channel-specific contributions via DAG-based path analysis. The original MATLAB implementation limits accessibility for the large R-based macroeconometrics community. The `vars` package [@Pfaff2008vars] provides SVAR estimation but not channel decomposition. `SVARtca` bridges this gap with a native R implementation that interfaces directly with `vars::VAR` objects.

# Usage

## From a VAR Object

The `tca_from_var` function extracts coefficient matrices from a fitted `vars::VAR` object and performs TCA directly:

```r
library(SVARtca)
library(vars)

var_est <- VAR(macro_data, p = 2, type = "const")

result <- tca_from_var(var_est,
                       from = "interest_rate",
                       intermediates = c("credit", "exchange_rate"),
                       h = 20, mode = "overlapping")
print(result, target = 3)
plot_tca(result, target = "output")
```

## Low-Level Systems Form Interface

For custom SVAR specifications, construct the systems form directly:

```r
sf <- tca_systems_form(Phi0, As, h = 20, order = 1:K)
result <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                      intermediates = c(2, 3), K = K, h = 20,
                      order = 1:K, mode = "exhaustive_4way",
                      var_names = c("y1", "y2", "y3", "y4"))
```

## Decomposition Modes

Three modes are available. *Overlapping* computes per-channel contributions that may sum to more than the total when channels share intermediate variables. *Exhaustive 3-way* decomposes into an inclusive channel, an exclusive channel, and the direct effect. *Exhaustive 4-way* applies full inclusion-exclusion, yielding four non-overlapping components: each channel alone, both jointly, and the direct effect.

# Implementation

The SVAR in systems form:

$$\mathbf{y}_t = \sum_{j=1}^{p} \mathbf{\Phi}_j \mathbf{y}_{t-j} + \mathbf{B} \boldsymbol{\varepsilon}_t$$

The structural MA representation $\mathbf{y}_t = \sum_{s=0}^{\infty} \mathbf{\Theta}_s \boldsymbol{\varepsilon}_{t-s}$ is decomposed via the DAG adjacency structure implied by the contemporaneous impact matrix $\mathbf{B}$.

`SVARtca` is implemented in pure R with dependencies on `Matrix`, `ggplot2`, and `rlang`. It accepts Cholesky or user-supplied identification and is compatible with SVAR objects from `vars`.

# Acknowledgements

The author acknowledges Emanuel Wegner, Lenard Lieb, and Stephan Smeekes for developing TCA and the original MATLAB toolbox.

# References
