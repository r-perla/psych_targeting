# Psychological Targeting Meta-Analysis

This repository hosts the materials for a meta-analysis on psychological targeting studies (specifics will be added after peer-review). The project includes:

- **Stan Model**: Implements a meta-analysis that combines machine learning prediction accuracies with advertising targeting effects across personality dimensions. See the [Stan model](model.stan) for details.
- **Simulation Study**: Demonstrates four experimental designs to estimate the targeting effect. The simulation (in [design_simulation.Rmd](design_simulation.Rmd)) shows that only one design yields an unbiased estimate. The designs are:
  - *Trait-Treatment Interaction*
  - *Selected Population Design*
  - *Matched-mismatched Design*
  - *Trait-only Design*

## Overview

The meta-analysis uses a Stan model to assess advertising targeting effects, while the simulation study evaluates the performance of different experimental designs. The simulation generates synthetic data, applies various regression models, and visualizes bias in the estimates, highlighting the strengths and weaknesses of each design.

## Repository Structure

- **design_simulation.html**: Rendered HTML output for the simulation study. Unfortunetaly, GitHub does not support rendering files of this size. You may download the file and open it in any browser.
- **design_simulation.Rmd**: R Markdown file containing code to generate simulation data, run the designs, conduct a grid search over parameter combinations, and visualize estimation bias.
- **model.stan**: Stan file for the meta-analysis model.
- **README.md**: This file.

## Getting Started

1. **Simulation Study**  
   - Open [design_simulation.Rmd](design_simulation.Rmd) in RStudio or run it using your preferred R environment.
   - The script sets global parameters, generates data via `gen_data()`, executes the designs using `get_results()`, performs a grid search with `grid_search()`, and creates a 3D visualization of bias using `plot_3d_bias()`.

2. **Meta-Analysis Model**  
   - Use [model.stan](model.stan) with your Stan-compatible interface to conduct the meta-analysis on advertising targeting effects.

## Reproducibility

Ensure that your R environment has the required packages:
- dplyr
- ggplot2
- plotly
- gridExtra
- foreach
- doParallel
- readr

Seed for reproducibility is set in the simulation script.

## Author
Raphael Perla

## License
[Specify your project's license here]