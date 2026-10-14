# DynamicDOOR

Estimate visit-specific and weighted-average win probabilities for repeated
ordinal outcomes using AIPW or TMLE, with optional cross-fitting and variance
correction. Larger outcome categories are more desirable; ties receive half credit.

## Installation

Install from GitHub:

```r
install.packages("remotes")  # if needed
remotes::install_github("DynamicDOOR/DynamicDOOR")
```

This is an anonymous review copy. The maintainer name and email in
`DESCRIPTION` are placeholders for review purposes.

## Example

```r
library(DynamicDOOR)

obs_data <- gen_toy_data(n = 1000, K = 4, seed = 12)$observed_data
times <- c(5, 10, 15, 20, 25, 30)
predictors <- paste0("X", 1:5)

fit <- fit_dynamic_door(
  obs_data,
  treatment = "Treatment",
  outcome_columns = paste0("PrimaryDOOROutcome_collapsed_", times),
  trt_predictors = predictors,
  outcome_predictors = predictors,
  outcome_levels = as.character(1:4),
  times = times,
  method = "TMLE-P", cf = TRUE, V = 5, seed = 1,
  outcome = "polr", varc = TRUE, correction = "proposition",
  scale = FALSE, ci_trajectory = "both", plot = TRUE
)

fit$results      # Visit-specific estimates and weighted summary
fit$joint_tests  # Joint trajectory tests
```

## Data and estimators

Use one row per participant, treatment coded 0/1, baseline covariates, and one
outcome column per visit. Analysis columns must be complete. Specify at least
three `outcome_levels`, ordered from least to most desirable. For long-format
data, use `long_to_wide()`.

| Estimator | `method` | `cf` |
| --- | --- | --- |
| AIPW | `"AIPW"` | `FALSE` |
| AIPW-CF | `"AIPW"` | `TRUE` |
| TMLE | `"TMLE"` | `FALSE` |
| TMLE-CF, fold-specific targeting | `"TMLE"` | `TRUE` |
| TMLE-P, pooled targeting | `"TMLE-P"` | `TRUE` |

Outcome models are fitted by treatment arm and visit. Choose `outcome = "polr"`
for proportional-probit regression, `"multinom"` for multinomial logistic
regression, or `"tree"` for classification trees.

## Inference

Specify `outcome`, `varc`, and `correction` explicitly in every call.

| `varc` | `correction` | Procedure |
| --- | --- | --- |
| `FALSE` | `NULL` | Ordinary fitted influence-function covariance |
| `TRUE` | `"proposition"` | Parametric correction for `polr` or `multinom` |
| `TRUE` | `"corollary"` | Propensity-score correction under a correctly specified fitted propensity model |
| `TRUE` | `"bootstrap"` | Covariance from refitted bootstrap estimates |

Use `ci_trajectory = "pointwise"`, `"simul"`, or `"both"`. Simultaneous intervals
are coordinate projections of a joint Wald region. Corrected output columns
end in `_corrected`. The weighted summary uses equal visit weights unless
`weights` is supplied.

## Documentation

See `?fit_dynamic_door` for all arguments and the
tutorial in `vignettes/Tutorial_dynamic_DOOR.Rmd` for long-format data,
all five estimators, and trajectory plots.
