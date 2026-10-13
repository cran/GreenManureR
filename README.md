# GreenManureR

Green manure decomposition, nutrient release and subsequent crop response.

GreenManureR provides tools for experiments with Mucuna, sunhemp, dhaincha and
other green manure and cover crops.

| Function | Purpose |
|---|---|
| `decomp_fit()` | Fit single-pool, double-pool and asymptotic decomposition models |
| `decomp_t50()` | Time to 50% (or any %) mass loss, with Monte-Carlo intervals |
| `nutrient_release()` | N, P, K and other nutrient release over time |
| `greenmanure_n_balance()` | N inputs, crop uptake, apparent balance, ANR, NFRV |
| `carryover_effect()` | Residual effects on kharif and rabi crops |
| `gm_joint_analysis()` | Link decomposition, release and crop response |
| `greenmanure_report()` | Standardized tables, plots and a Word report |

## Installation

```r
# once on CRAN
install.packages("GreenManureR")
```

## Example

```r
library(GreenManureR)
fit <- decomp_fit(gm_litterbag, group = "treatment")
decomp_t50(fit)
```

See `vignette("greenmanure-workflow")` for a full walk-through.
