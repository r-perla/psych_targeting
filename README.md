# Psychological Targeting Meta-Analysis
This repository hosts materials for the paper [The (In)Effectiveness of Psychological Targeting: A Meta-Analytic Review](https://doi.org/10.1002/mar.70073). The project comprises:

- **Stan Model**: A meta-analytic model that combines ML prediction accuracies with advertising targeting effects across personality dimensions. See [model.stan](model.stan).
- **Simulation Study**: A data-generating simulation demonstrating four experimental designs for estimating a personalization (targeting) effect. See [design_simulation.Rmd](design_simulation.Rmd).
- **Closed-Form Analysis**: An analytic (non-simulated) derivation showing why alternative designs fail to identify the personalization effect and by how much they deviate. See [design_failures_closed_form.R](design_failures_closed_form.R).

## Overview
The meta-analysis evaluates advertising targeting effects using a Bayesian model in Stan. The empirical design question is addressed two ways:

1. **Simulation (RMarkdown):** Generates synthetic data under  
   $y=b_0+b_1T+b_2C+b_3T\cdot C+\varepsilon$ with randomized treatment $C$ and trait $T\sim\mathcal N(0,1)$. It compares four designs:
   - *Trait–Treatment Interaction* (identifies $b_3$)
   - *Selected Population* (restrict $T>q$; mixes $b_2$ and $b_3$)
   - *Matched–Mismatched* (extremes of $T$; scales $b_3$ by a truncated-normal mean)
   - *Trait-only* within $C=1$ (confounds $b_1$ with $b_3$)

2. **Closed-Form (R script):** Derives the limiting estimands of each design without simulation. It highlights **estimand mismatch**: designs 2–4 converge to parameters that are not $b_3$. For upper-tail threshold $p\in(0.5,1)$ with $q=\mathrm{qnorm}(p)$, the truncated-normal mean is $\lambda_p=\phi(q)/\{1-\Phi(q)\}$. The designs converge to:
   - Trait×Treatment: $b_3$
   - Selected Population: $b_2 + b_3\,\lambda_p$
   - Matched–Mismatched: $b_3\,\lambda_p$
   - Trait-only (within $C=1$): $b_1 + b_3$

## Repository Structure
- **design_simulation.Rmd**: Simulation code: DGP, grid search, and 3D bias visualization.
- **design_simulation.html**: Rendered simulation (download to view).
- **design_failures_closed_form.R**: Closed-form derivations and 2×2 heatmaps of deviation from $b_3$.
- **model.stan**: Stan model for the meta-analysis.
- **README.md**: This file.

## Getting Started
### 1) Closed-Form Analysis (recommended for speed and lack of sampling error)
- Run from the command line:
  ```Rscript design_failures_closed_form.R```

- The script prints a tidy table of design-specific limits and their **deviation from (b_3)**, plus a small 2×2 heatmap figure.
- To change the threshold or grid:
  - Edit `p` (e.g., `p <- 0.75` for top quartile, `0.5` for median split, `0.9` for top decile).
  - Edit `b1_vals` and `b2_vals` to explore main-effect contamination.

### 2) Simulation Study
- Open [design_simulation.Rmd](design_simulation.Rmd) in RStudio or run via:

  ```r
  rmarkdown::render("design_simulation.Rmd")
  ```
- The document sets global parameters, generates data with `gen_data()`, computes design estimates via `get_results()`, runs a grid with `grid_search()`, and visualizes surfaces using `plot_3d_bias()`.

## Dependencies
**Closed-Form Script**
`dplyr`, `tidyr`, `ggplot2`

**Simulation Study**
`dplyr`, `ggplot2`, `plotly`, `gridExtra`, `foreach`, `doParallel`, `readr`

## Reproducibility
- **Closed-Form**: Fully deterministic; no random simulation involved.
- **Simulation**: Seeded in the Rmd.
