## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5, fig.height = 4)
has_gg <- requireNamespace("ggplot2", quietly = TRUE)
library(GreenManureR)

## ----data---------------------------------------------------------------------
head(gm_litterbag, 3)
head(gm_agronomy, 3)

## ----decomp-------------------------------------------------------------------
fit <- decomp_fit(gm_litterbag, time = "time", mass = "mass", group = "treatment")
fit$table[fit$table$selected, c("treatment", "model", "k", "k_slow", "pool_A", "R2", "AICc")]

## ----t50----------------------------------------------------------------------
decomp_t50(fit)

## ----decomp-plot, eval = has_gg-----------------------------------------------
plot(fit)

## ----nutrient-----------------------------------------------------------------
nr <- nutrient_release(gm_litterbag, nutrients = c("N", "P", "K"),
                       group = "treatment", biomass = "biomass_tha",
                       at = c(14, 30, 60))
nr$t50[, c("treatment", "nutrient", "model", "t50")]
subset(nr$release_at, nutrient == "N" & time == 30)

## ----nutrient-plot, eval = has_gg---------------------------------------------
plot(nr)

## ----nbal---------------------------------------------------------------------
nb <- greenmanure_n_balance(gm_agronomy, by = "season")
s <- nb$summary
s[s$fert_N == 0 & s$treatment != "Control",
  c("season", "treatment", "N_gm_input_mean", "N_uptake_mean", "ANR_mean", "NFRV_mean")]

## ----carry--------------------------------------------------------------------
co <- carryover_effect(gm_agronomy, response = "grain_yield", strata = "fert_N",
                       season_levels = c("Kharif", "Rabi"))
co$effects[co$effects$fert_N == 0,
           c("season", "treatment", "diff", "pct_change", "p_adj", "retention_vs_first")]
co$pooled

## ----carry-plot, eval = has_gg------------------------------------------------
plot(co)

## ----joint--------------------------------------------------------------------
ja <- suppressWarnings(
  gm_joint_analysis(fit, nr, nb, co, nutrient_name = "N", at_day = 60, fert_N_level = 0))
ja$data

## ----report, eval = FALSE-----------------------------------------------------
#  greenmanure_report(fit, nr, nb, co, ja, outdir = "gm_report")

