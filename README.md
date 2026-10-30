# wLC-DLNM: A penalized distributed lag non-linear Lee–Carter model for weekly mortality

This repository contains the R code for the paper

> **A penalized distributed lag non-linear Lee–Carter framework for regional weekly mortality forecasting**
> Jens Robben and Karim Barigou
> [arXiv:2509.24087](https://arxiv.org/abs/2509.24087) · [PDF](https://arxiv.org/pdf/2509.24087) · [Online Appendix](https://jensrobben.github.io/articles/wLC-DLNM/Online%20Appendix.pdf)

The framework forecasts **weekly, regional, age-specific mortality**. It separates the usual seasonal pattern from short-term excess deaths caused by **heat, cold and influenza**.

---

## What the model does

The model extends the Lee–Carter model to weekly data for several regions. The log force of mortality for age group $`x`$, year $`t`$, week $`w`$ and region $`r`$ is:

```math
\log \mu_{x,t,w,r} =
\underbrace{\alpha_{x,r} + \beta_x \kappa_{t,r} + \gamma_x \lambda_{w,r}}_{\text{seasonal baseline}} +
\underbrace{\delta_x f^{(1)}_r(\text{Tavg}_{t,w,r}) + \epsilon_x f^{(2)}_r(\text{ILI}_{t,w,r})}_{\text{temperature and influenza effects}}
```

Death counts follow a negative binomial distribution, which allows for overdispersion.

Main features:

- **Seasonal baseline.** A yearly mortality index $`\kappa_{t,r}`$ and a separate week-of-year effect $`\lambda_{w,r}`$, each with its own age sensitivity.
- **Distributed lag non-linear models (DLNMs).** These capture the delayed, non-linear effects of weekly average temperature and influenza-like illness (ILI) incidence. The cross-bases are centred on their week-of-year means, so the DLNM terms describe departures from typical conditions only.
- **Spatial smoothing.** A Laplacian penalty pulls the temperature and influenza effects of neighbouring regions towards each other. The smoothing parameters are chosen by an effective-degrees-of-freedom BIC.
- **Forecasting.** SARIMA models for the mortality index, temperature and ILI. Dependence between regions is captured with Gaussian or t-copulas fitted on shrunk correlation matrices.
- **Uncertainty.** Parametric bootstrap confidence intervals for the parameters, and simulated prediction intervals for future deaths.

## Case study

- **Mortality:** individual death records from INSEE, aggregated to weekly death counts for women by 5-year age group (50–54 to 95+).
- **Regions:** the 12 administrative regions of metropolitan France (excluding Corsica).
- **Period:** 1990–2019. The model is fitted on 1990–2014 and tested out of sample on 2015–2019. It is also used to estimate excess mortality over 2020–2024.
- **Benchmarks:** nested versions of the model, plus the Li–Lee model, Guibert et al. (2026), Min et al. (2026) and Robben et al. (2026).

---

## Repository structure

The scripts fall into four groups, marked by the prefix in their file names.

### Main pipeline (run in order)

| Script | What it does |
|---|---|
| `0. Read E-OBS - NUTS 2.R` | Reads the gridded E-OBS daily temperatures, intersects them with the regional boundaries and computes population-weighted daily regional averages (`FR_NUTS1_tg_daily.txt`). |
| `1. Read data.R` | Cleans the INSEE death records, aggregates them by ISO week, 5-year age group, sex and region, and computes weekly exposures. It then merges in the influenza and temperature data to build the final data set (`Df_var.txt`). |
| `2. Fit model.R` | Builds the cross-basis matrices, selects the smoothing parameters, fits the model and computes the bootstrap and Hessian-based confidence intervals. It also produces the relative risk plots, residual diagnostics, the Poisson comparison and the 2020–2024 excess-mortality forecasts (Sections 5.1–5.3 and 5.5; Suppl. Mat. B, C, D and F). |

### `B.` Benchmarks and out-of-sample evaluation (Section 5.4)

| Script | Model |
|---|---|
| `B. Benchmark model 1-3.R` | Models 1–3: annual regional, weekly national and weekly regional Lee–Carter models without covariates |
| `B. Benchmark model 4-6.R` | Models 4–6: with temperature, influenza and both DLNMs, under observed and simulated covariates. Also contains the uncertainty decomposition (Table 6) and the (S)ARIMAX comparison (Suppl. Mat. E). |
| `B. Benchmark model 6 (no penal.).R` | Model 6 without spatial penalty (Suppl. Mat. C) |
| `B. Benchmark model 7 (Li-Lee).R` | Li–Lee model (N. Li & Lee, 2005) |
| `B. Benchmark model 8 (Guibert et. al).R` | Two-stage DLNM Li–Lee model (Guibert et al., 2026) |
| `B. Benchmark model 9 (Min et. al).R` | Weekly DLNM Li–Lee model (Min et al., 2026) |
| `B. Benchmark model 10 (Robben et. al).R` | Three-state regime-switching model (Robben et al., 2026) |

### `S.` Sensitivity analyses (Supplementary Material)

| Script | Analysis |
|---|---|
| `S. Sensitivity to DLNM specification.R` | Maximum lags, spline bases and knot positions (Suppl. Mat. A) |
| `S. Sensitivity to temperature specification.R` | Weekly mean, maximum, minimum and maximum 3-day average temperature (Suppl. Mat. G) |
| `S. Sensitivity to dispersion specification.R` | Week-specific dispersion parameter (Section 3.1) |

### `F.` Function files

These are sourced by the scripts above and are not run on their own:

- `F. Weekly Lee-Carter (NB).R`, `F. Weekly DLNM (NB).R`, `F. Weekly Lee-Carter + DLNM (NB).R`: likelihoods, gradients, Hessians and the blockwise Newton fitting routine for the proposed model.
- `F. ... (POI).R`: the same under a Poisson assumption.
- `F. ... + time disp (NB).R`: versions with a week-specific dispersion parameter.
- `F. Annual, regional Lee-Carter (NB).R`, `F. Annual, regional Li-Lee (NB).R`: annual benchmark models.
- `F. Guibert et al. approach.R`, `F. Min et. al approach.R`, `F. Robben et al. approach.R`: the literature benchmarks.
- `F. Bootstrap CI.R`: the parametric bootstrap.

---

## Data

All data are publicly available, but they are **not included** in this repository. The scripts expect the data in folders **next to** the repository folder, as follows:

```
parent folder/
├── wLC-DLNM/                 # this repository
├── Deaths/                   # deces-YYYY.csv: INSEE death records, one file per year
├── Population/
│   ├── estim-pop-dep-sexe-aq-1975-2026.xlsx   # INSEE population estimates
│   └── PopulationCount.nc                      # gridded population (GPW, CIESIN)
├── Temperature/
│   └── tg_ens_mean_0.1deg_reg_v31.0e.nc        # E-OBS daily mean temperature
├── Sentinel/
│   └── inc-3-RDD.csv                           # Sentinelles ILI incidence by region
├── NUTS/
│   └── NUTS_RG_03M_2024_3035.shp (+ .dbf, .prj, .shx)   # Eurostat NUTS boundaries
├── v_departement_2023.csv    # INSEE geographic code tables
├── v_region_2023.csv
├── Results updated/          # intermediate results (written by the scripts)
└── Figures updated/          # figures (written by the scripts)
```

| Data | Source |
|---|---|
| Individual death records | [INSEE – fichiers des personnes décédées](https://www.insee.fr/fr/information/4769950) |
| Population estimates | [INSEE – population estimates](https://www.insee.fr/fr/statistiques/8331297) |
| Department and region code tables | [INSEE – Code officiel géographique 2023](https://www.insee.fr/fr/information/7766585) |
| Daily mean temperature (E-OBS v31.0e, 0.1°) | [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets/insitu-gridded-observations-europe) |
| Gridded population | NASA SEDAC – Gridded Population of the World (CIESIN) |
| Regional boundaries | [Eurostat GISCO – NUTS 2024](https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics) |
| Influenza-like illness incidence | [Réseau Sentinelles](https://www.sentiweb.fr/) |

---

## Reproducing the results

### Requirements

The code is written in R and uses the following packages:

```r
install.packages(c(
  "dplyr", "tidyr", "tidyverse", "tibble", "purrr", "data.table", "collapse",
  "lubridate", "ISOweek", "readxl", "imputeTS", "stringr",
  "sf", "spdep", "ncdf4", "metR",
  "dlnm", "tsModel", "splines", "mgcv", "Matrix", "MASS", "abind", "matrixStats",
  "forecast", "copula", "corpcor", "mvnfast", "scoringRules",
  "parallel", "doParallel", "snow", "doSNOW", "pbapply", "progress", "cli",
  "ggplot2", "ggpubr", "Polychrome", "viridis", "scales", "xtable",
  "plyr", "rstudioapi"
))
```

### Steps

1. Clone the repository and put the data in the folder layout shown above:
   ```bash
   git clone https://github.com/jensrobben/wLC-DLNM.git
   ```
2. Open the scripts in **RStudio**. Each script sets its working directory to its own folder with `rstudioapi`, so it must be run from RStudio. If you run it elsewhere, replace the `setwd(...)` line with the path to the repository.
3. Run `0. Read E-OBS - NUTS 2.R`, then `1. Read data.R`, then `2. Fit model.R`.
4. Run the `B.` scripts for the out-of-sample comparison and the `S.` scripts for the sensitivity analyses. Run `B. Benchmark model 4-6.R` before `B. Benchmark model 1-3.R`, because the latter reads the scores that the former saves.

**Note on saved results.** Some steps are slow, such as the smoothing-parameter grid search, the bootstrap with 5 000 replications and the forecast simulations. Their results are saved to `../Results updated/` with `saveRDS()` and then read back with `readRDS()`. In the scripts as provided, the `saveRDS()` lines are commented out. **On a first run, uncomment them** so the intermediate files are created. The same applies to the `fwrite()` lines in `1. Read data.R` that write `Df_final.txt` and `Df_var.txt`. Some steps run in parallel; set the number of cores in `makeCluster()` to suit your machine.

---

## Paper and supplementary material

- **Paper:** [arXiv:2509.24087](https://arxiv.org/abs/2509.24087). The Supplementary Material (sensitivity analyses, model comparisons and benchmark specifications, referred to as "Suppl. Mat." in the scripts above) is at the end of the paper PDF.
- **Online Appendix:** [Online Appendix.pdf](https://jensrobben.github.io/articles/wLC-DLNM/Online%20Appendix.pdf), with the in-sample fits, relative risk curves and relative risk surfaces for all regions and age groups.

## Citation

If you use this code, please cite the paper:

```bibtex
@misc{robben_barigou_2025_wlcdlnm,
  author        = {Robben, Jens and Barigou, Karim},
  title         = {A penalized distributed lag non-linear {Lee--Carter} framework for regional weekly mortality forecasting},
  year          = {2025},
  eprint        = {2509.24087},
  archivePrefix = {arXiv},
  url           = {https://arxiv.org/abs/2509.24087}
}
```

## Contact

Jens Robben, Research Centre for Longevity Risk, University of Amsterdam: j.robben@uva.nl
