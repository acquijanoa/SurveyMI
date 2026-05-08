# SurveyMI

`SurveyMI` is an R package for **multiple imputation under complex survey designs**.  
Its central goal is to produce imputations and model estimates that respect survey features such as unequal weights, stratification, and clustering, instead of treating the data as if it came from simple random sampling.

## Objective

Missing data is common in health, social, and economic surveys. Standard multiple imputation approaches are often built for i.i.d. data and may under-represent uncertainty when applied directly to complex survey samples.

`SurveyMI` is designed to address this gap by combining:

- sequential/chained imputation of incomplete variables,
- survey-aware weighting through multinomial bootstrap replicate weights,
- repeated model fitting across weighted imputations.

In short, the package objective is to support **valid, practical inference with incomplete complex survey data**.

## What the package does

At a high level, `survey_mi()`:

1. Takes a `survey::svydesign()` object.
2. Generates `B` multinomial bootstrap replicate weights.
3. Imputes missing values using weighted conditional models.
4. Fits the requested analysis model for each replicate.
5. Returns replicate-level coefficients plus an empirical summary.

Current model families for imputation/analysis include:

- `gaussian`
- `binomial`
- `poisson`

## Installation

You can install from source in the repository:

```r
# install.packages("devtools")
devtools::install(".")
```

Or during development:

```r
devtools::load_all()
```

## Minimal example

```r
library(survey)
library(SurveyMI)

# Example survey design
dsgn <- svydesign(
  data = data,
  strata = ~strat,
  ids = ~BGid,
  weights = ~bghhsub_s2,
  nest = TRUE
)

fit <- survey_mi(
  design = dsgn,
  analysis = y_gfr ~ x17 + y_bmi + x12 + x14 + x18,
  variables = c("y_gfr", "y_bmi", "x12", "x14", "x17", "x18"),
  methods = list(
    gaussian = "x17",
    binomial = "x18"
  ),
  analysis_family = "gaussian",
  B = 100,
  n_rounds = 10,
  seed = 123
)

fit$summary
```

## Main exported functions

- `survey_mi()`: end-to-end MI + analysis for complex surveys.
- `multinom_weights()`: generates multinomial bootstrap replicate weights.
- `has_rcpp_fitters()`: checks whether Rcpp accelerated fitters are available.

## Status

This repository is currently in early development (`0.0.0.9000`).  
Interfaces and defaults may evolve as methods are refined and validated.
