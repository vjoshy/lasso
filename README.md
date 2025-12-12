# LASSO Regression: Computational Methods

This repository contains R code for the STAT*6801 final project implementing Coordinate Descent and Proximal Gradient Descent algorithms for LASSO regression, with Strong Sequential Screening Rules.

## Repository Structure
```
├── functions/
│   ├── lasso.R                    # Core algorithms (CD, PGD, screening rules)
│   ├── cv.R                       # Cross-validation and regularization path
│   ├── generate_datasets.R        # Synthetic data generation
│   ├── run_simulation.R           # Sequential simulation script
│   ├── run_simulation_parallel.R  # Parallel simulation script
│   └── simulation_results.R       # Results analysis and plotting
├── data/                          # Generated datasets (not tracked)
└── test.R                         # Basic testing script (not updated)
```

## Core Implementations

The main algorithms are in `functions/lasso.R`:

- `cd_lasso()` — Coordinate Descent with soft-thresholding
- `pgd_lasso()` — Proximal Gradient Descent (ISTA)
- `screen_variables()` — Basic Strong Rule screening
- `screen_variables_sequential()` — Sequential Strong Rule screening

## Running the Simulation

1. Generate datasets: `source("functions/generate_datasets.R")`
2. Run simulation: `source("functions/run_simulation_parallel.R")`
3. Analyze results: `source("functions/simulation_results.R")`

## Dependencies

- `MASS` , `glmnet` , `dplyr`, `tidyr`, `ggplot2`  `foreach`, `doParallel` 

## Author

Vinay Joshy  
University of Guelph, Fall 2025