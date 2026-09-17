# dynamic-survival-delisting
This repository contains the core R code accompanying the manuscript
"Dynamic Survival Analysis for Predicting US Stock Delisting (Yoon and Fan)".

## Contents

- `setup.R`: set up packages and source functions (mostly modified functions of LTRCforests and pcoxtime packages)
- `min_depth_sim_linear.R`: simulation study for minimal depth variable importance measure
- `train_pcox.R`: estimation of the regularized Cox models (LASSO, ridge, elastic net)
- `train_forest.R`: estimation of the LTRC-CIF model
- `eval_TVAUC.R`: incident/dynamic time-dependent AUC
- `eval_BS_TV.R`: time-dependent Brier Scores and Integrated Brier Scores
- the folder `functions` contains helper functions for running the simulations, training the models and evaluating model performance

## Data

The empirical analysis uses CRSP and S&P Global Market Intelligence
Compustat data obtained through WRDS. These proprietary data cannot be
redistributed through this repository.

## Software

Analyses were conducted in R. Required packages include:
`LTRCforests`, `pcoxtime`, `partykit`, `survival`, `risksetROC`, `mice`, ...

## Reproducibility

The simulation study can be reproduced using the code provided in this
repository. Reproduction of the empirical analysis requires access to
the underlying CRSP and Compustat data.
