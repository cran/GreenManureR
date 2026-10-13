#' GreenManureR: Green Manure Decomposition, Nutrient Release and Crop Response
#'
#' Tools for green manure and cover crop experiments (Mucuna, sunhemp,
#' dhaincha and others): litterbag decomposition models, nutrient release,
#' nitrogen balance, residual (carry-over) effects on succeeding crops, a
#' joint analysis linking the three, and standardized reporting.
#'
#' @section Main functions:
#' * [decomp_fit()] fits single-pool, double-pool and asymptotic decay models.
#' * [decomp_t50()] estimates the time to 50% (or any other %) mass loss.
#' * [nutrient_release()] quantifies N, P, K and other nutrient release.
#' * [greenmanure_n_balance()] computes N inputs, crop uptake and balances.
#' * [carryover_effect()] quantifies residual effects on kharif and rabi crops.
#' * [gm_joint_analysis()] links decomposition, release and crop response.
#' * [greenmanure_report()] writes standardized tables, figures and a report.
#'
#' Example data sets (simulated) are provided as [gm_litterbag] and
#' [gm_agronomy]; see also [simulate_greenmanure()].
#'
#' @keywords internal
#' @importFrom stats AIC aggregate anova coef cor deviance fitted lm nls
#' @importFrom stats nls.control p.adjust predict qt quantile reshape resid
#' @importFrom stats rnorm sd setNames uniroot vcov
#' @importFrom utils write.csv
"_PACKAGE"

# Column names used in ggplot2 aesthetics and nls formulas
utils::globalVariables(c(
  ".t", ".m", "grp", "time", "pred", "released", "rel", "fitted", "resid",
  "ANR_mean", "ANR_se", "fert", "trt", "pct_change", "pct_lwr", "pct_upr",
  "treatment", "strat", "k", "A", "k1", "r", "C", "m"
))
