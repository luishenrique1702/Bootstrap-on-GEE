# Bootstrap on GEE

This repository contains R scripts for statistical applications using Generalized Estimating Equations (GEE). These codes were developed to implement Bootstrap resampling with BCa (Bias-Corrected and accelerated) intervals and simulated envelopes for randomized quantile residuals in GEE models.

## Overview
The project focuses on two main applications for longitudinal data analysis:
- **Application 1 (Bernoulli):** Analysis of the `polypharm` dataset, focusing on polypharmacy among the elderly.
- **Application 2 (Gamma):** Analysis of the `Milk` dataset, modeling protein content in milk over time.
- 
## Academic Reference
The codes in this repository are extracted from my Master's dissertation:
> **Toledo, L. H. C. P.** (2025). *Técnicas de reamostragem e detecção de influência mascarada em equações de estimação generalizadas*. Dissertation (Master's Degree) - Institute of Mathematics and Statistics, University of São Paulo (IME-USP).

The full document is available at the **[USP Digital Library of Theses and Dissertations](https://www.teses.usp.br/teses/disponiveis/45/45133/tde-08092025-150349/publico/Dissertacao_Luis_H_Toledo_Versao_Corrigida.pdf)**.

## Data Sources
The scripts use datasets available in specific R libraries. Ensure these are installed to load the data:
- **`polypharm`**: Found in the **`aplore3`** package. It contains data from the Well-Being of the Elderly (WBE) study.
- **`Milk`**: Found in the **`gamlss.data`** (or loaded via **`gamlss`**) package. It consists of longitudinal data on the protein content of cows' milk.

## Methodology
The implementations include:
1. **GEE Modeling:** Parameter estimation using the `glmtoolbox` package, considering AR(2) correlation structures.
2. **Bootstrap BCa:** Resampling method for obtaining Bias-Corrected and accelerated confidence intervals, providing robust inference for regression, dispersion and correlation parameters.
3. **Simulated Envelopes:** Diagnostic tool for randomized quantile residuals. The simulation process employs **Gaussian copulas** to account for the intra-cluster dependence structure, ensuring the envelopes correctly reflect the GEE framework.

## Required Packages
To run these scripts, install the following dependencies:
```R
install.packages(c("glmtoolbox", "dplyr", "ggplot2", "aplore3", "mvnfast", "gamlss", "tidyverse"))
