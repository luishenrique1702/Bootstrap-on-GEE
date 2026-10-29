# Bootstrap on GEE

This repository contains the R code for statistical applications using Generalized Estimating Equations (GEE). The code implements cluster bootstrap resampling with BCa (Bias-Corrected and accelerated) intervals for all the GEE parameters and simulated envelopes for quantile residuals in GEE models, together with the simulation studies that assess these procedures.

## Overview
The project focuses on two applications for longitudinal data analysis:
- **Application 1 (Bernoulli):** Analysis of the `polypharm` dataset, relating the use of polypharmacy (more than three medications) by 500 patients observed from 2002 to 2008 to their characteristics, with a GEE logistic model.
- **Application 2 (Gamma):** Analysis of the `Milk` dataset, modeling the protein content in the milk of 79 cows over the weeks after calving, with a GEE-gamma log-linear model.

All the results of the paper *On bootstrap in generalized estimating equations* (Toledo and Paula) and of its supplementary materials are produced by the code in this repository.

## Repository Structure
| Path | Content |
|---|---|
| `functions.R` | Functions shared by both applications: GEE fitting, cluster bootstrap and BCa intervals, quantile residuals, heat plots, Gaussian copula and NORTA calibration, simulated envelopes and global test, simulation tools. |
| `Application 1 - Bernoulli/Bernoulli.Rmd` | Application 1: data and model selection, bootstrap, envelope, heat plots and worm plots, envelope diagnostics, type I error and power of the envelope test, simulation study. |
| `Application 2 - Gamma/Gamma.Rmd` | Application 2: data, bootstrap, envelope, heat plots and worm plots, sensitivity analysis, type I error of the envelope test, simulation study. |
| `*.html` (in each folder) | Rendered reports, with the code and all the outputs (tables and figures). Download them and open them in a web browser. |

The rendered reports can be viewed online:
- [Application 1 (Bernoulli)](https://luishenrique1702.github.io/Bootstrap-on-GEE/Application%201%20-%20Bernoulli/Bernoulli.html)
- [Application 2 (Gamma)](https://luishenrique1702.github.io/Bootstrap-on-GEE/Application%202%20-%20Gamma/Gamma.html)
  
## Academic Reference
The paper is based on the Master's dissertation:
> **Toledo, L. H. C. P.** (2025). *Técnicas de reamostragem e detecção de influência mascarada em equações de estimação generalizadas*. Dissertation (Master's Degree) - Institute of Mathematics and Statistics, University of São Paulo (IME-USP).

The full document is available at the **[USP Digital Library of Theses and Dissertations](https://www.teses.usp.br/teses/disponiveis/45/45133/tde-08092025-150349/publico/Dissertacao_Luis_H_Toledo_Versao_Corrigida.pdf)**.

## Data Sources
The datasets are available in R packages:
- **`polypharm`**: Found in the **`aplore3`** package.
- **`Milk`**: Found in the **`nlme`** package. It consists of longitudinal data on the protein content of cows' milk.

## Methodology
The implementation includes:
1. **GEE Modeling:** Parameter estimation with the `glmtoolbox` package, with AR(2) working correlation structures and lags measured in time units.
2. **Bootstrap BCa:** Cluster bootstrap with percentile, basic and BCa intervals, providing inference for the regression coefficients and also for the dispersion and correlation parameters.
3. **Heat Plots:** Quantile residuals against the fitted values and worm plots drawn as density plots of many randomizations of the residuals (function `heat_plot` in `functions.R`), with smooth trend curves that indicate where the model fails.
4. **Simulated Envelopes:** Diagnostic tool for quantile residuals. The responses are simulated through Gaussian copulas that reproduce the intra-cluster dependence of the fitted model (with NORTA calibration for binary responses), the model is refitted to every simulated sample, and a global Monte Carlo test complements the graph.
5. **Simulation Studies:** Bias of the estimators, coverage of the Wald and bootstrap intervals, and type I error and power of the envelope test.

## Required Packages
To run the code, install the following dependencies:
```R
install.packages(c("glmtoolbox", "aplore3", "nlme", "gamlss", "pbivnorm", "Matrix", "ggplot2", "pbapply", "viridisLite", "rmarkdown"))
```

## Running the Code
Keep the folder structure (each `.Rmd` loads `../functions.R`). Open an `.Rmd` file in RStudio and click **Knit**, or run, from the repository folder:
```R
rmarkdown::render("Application 1 - Bernoulli/Bernoulli.Rmd")
rmarkdown::render("Application 2 - Gamma/Gamma.Rmd")
```
The two applications are independent. The computations run in parallel; the complete run, with the Monte Carlo sizes of the paper, takes about 15 hours with 10 processes. The results are saved in the `results/` subfolder of each application and reused in later runs. The simulations of the paper were run with R 4.4.3 and `glmtoolbox` 0.1.12, and the report was rendered from their stored results.
