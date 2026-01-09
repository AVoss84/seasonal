# Seasonal Time Series Forecasting with Bayesian Model Averaging

R code used in: [Vosseler, A., Weber, E. (2018): Forecasting seasonal time series data: a Bayesian model averaging approach](https://www.researchgate.net/publication/312372580_Forecasting_seasonal_time_series_data_A_Bayesian_model_averaging_approach)

## Overview

This repository provides R implementations for Bayesian estimation and forecasting of **Periodic Autoregressive (PAR)** models with potential structural breaks. The methodology uses **Gibbs sampling** for posterior inference and **Bayesian Model Averaging (BMA)** to account for model uncertainty, including uncertainty about the presence and timing of structural breaks.

## Key Features

- Bayesian estimation of PAR(p) models for quarterly (S=4) and monthly (S=12) time series
- Structural break detection with unknown break date
- Gibbs sampling for posterior inference
- Bayesian Model Averaging across models with/without breaks
- Marginal posterior computation for break date estimation
- Forecast generation with predictive density estimation
- Statistical tests for forecast evaluation

## Repository Structure

### Main Estimation Scripts

| File | Description |
|------|-------------|
| [par4_estimation_gibbs_B_with_BMA_marg.R](utils/par4_estimation_gibbs_B_with_BMA_marg.R) | Complete estimation and forecasting pipeline for **quarterly** PAR(p) models with BMA |
| [par12_estimation_gibbs_with_BMA_marg.R](utils/par12_estimation_gibbs_with_BMA_marg.R) | Complete estimation and forecasting pipeline for **monthly** PAR(p) models with BMA |

### Core Functions

| File | Description |
|------|-------------|
| [bayes.par.R](utils/bayes.par.R) | Bayesian PAR model setup for quarterly data (S=4). Creates design matrices, computes OLS estimates, and information criteria (AIC, BIC, HQ, FPE) |
| [bayes.par12.R](utils/bayes.par12.R) | Bayesian PAR model setup for monthly data (S=12). Analogous to `bayes.par.R` |
| [gibbs4.R](utils/gibbs4.R) | Gibbs sampler for quarterly PAR(p) models. Draws from posterior distributions of coefficients, variance, and break dates |
| [gibbs12.R](utils/gibbs12.R) | Gibbs sampler for monthly PAR(p) models. Analogous to `gibbs4.R` |

### Marginal Posterior for Break Date

| File | Description |
|------|-------------|
| [margpostTb.R](utils/margpostTb.R) | Computes marginal posterior density of break date $T_b$ for quarterly models. Includes HPD interval computation |
| [margpostTb12.R](utils/margpostTb12.R) | Computes marginal posterior density of break date $T_b$ for monthly models |

### Design Matrix Construction

| File | Description |
|------|-------------|
| [design.R](utils/design.R) | Constructs seasonal dummy design matrices with mean or effect coding |
| [design.break.R](utils/design.break.R) | Constructs design matrices accommodating structural breaks at time $T_b$ |
| [strbr.comp.R](utils/strbr.comp.R) | Comprehensive structural break component builder for various deterministic specifications |
| [mats01.R](utils/mats01.R) | Constructs $\Omega_0$ and $\Omega_1$ matrices for the companion form representation of PAR models |

### Simulation Functions

| File | Description |
|------|-------------|
| [par_simulation4.R](utils/par_simulation4.R) | Simulates quarterly PAR(1), PAR(2), PAR(3), PAR(4), and broken PAR processes. Includes unit root restrictions |
| [par_simulation12.R](utils/par_simulation12.R) | Simulates monthly PAR(1), PAR(2), and broken PAR processes |

### Forecasting Functions

| File | Description |
|------|-------------|
| [forec.sarima.R](utils/forec.sarima.R) | Forecasting using SARIMA models with seasonal dummies (benchmark comparison) |
| [forec.sdummy.R](utils/forec.sdummy.R) | Forecasting using seasonal dummy regression models |

### Forecast Evaluation

| File | Description |
|------|-------------|
| [predictive_tests.R](utils/predictive_tests.R) | Statistical tests for comparing predictive accuracy: Sign test, Wilcoxon test, and Bayesian posterior probability of equal predictive ability |
| [signtest.R](utils/signtest.R) | Implementation of the sign test for paired comparisons with confidence intervals |

## Dependencies

The code requires the following R packages:

```r
install.packages(c("MASS", "Bolstad", "TSA", "pear"))
```

## Usage Example

```r
# Source simulation functions
source("utils/par_simulation4.R")

# Simulate a quarterly PAR(1) process
y <- par1(N = 200, mu = c(1, 1, 1, 1), beta = 0, 
          phi = c(0.87, 1.2, 0.95, 0.85), start = c(1976, 2))$y

# Source Gibbs sampler
gibbs4 <- dget("utils/gibbs4.R")

# Run Gibbs sampling with structural break detection
results <- gibbs4(y = y, p = 1, n.ahead = 8, M = 5000, burnin = 1000,
                  sbreak = TRUE, in.samp.fc = TRUE, 
                  sdeter = "const", tau = 3)
```

## Model Specification

The PAR(p) model is specified as:

$$y_t = \sum_{s=1}^{S} D_{s,t} \left( \mu_s + \sum_{j=1}^{p} \phi_{s,j} y_{t-j} \right) + \varepsilon_t$$

where:
- $S$ = number of seasons (4 for quarterly, 12 for monthly)
- $D_{s,t}$ = seasonal dummy for season $s$
- $\mu_s$ = season-specific intercept
- $\phi_{s,j}$ = season-specific autoregressive coefficients
- $\varepsilon_t \sim N(0, \sigma^2)$

## Citation

If you use this code, please cite:

```bibtex
@article{vosseler2018forecasting,
  title={Forecasting seasonal time series data: A Bayesian model averaging approach},
  author={Vosseler, Alexander and Weber, Enzo},
  journal={Computational Statistics},
  volume={33},
  pages={1733--1765},
  year={2018},
  publisher={Springer}
}
```

## License

This project is provided for academic and research purposes.
