#' Simulated litterbag decomposition data
#'
#' Simulated litterbag data for three green manure species, with mass and
#' nutrient concentrations over 112 days. Generated with
#' [simulate_greenmanure()] (seed 42); the values are synthetic.
#'
#' @format A data frame with 108 rows and 8 variables:
#' \describe{
#'   \item{treatment}{Green manure species (Mucuna, Sunhemp, Dhaincha).}
#'   \item{rep}{Litterbag replicate (1-4).}
#'   \item{time}{Days after incorporation (0-112).}
#'   \item{mass}{Dry mass remaining (% of initial).}
#'   \item{N, P, K}{Nutrient concentrations (% dry weight).}
#'   \item{biomass_tha}{Green manure dry matter incorporated (t/ha).}
#' }
#' @source Simulated; see `data-raw/simulate-data.R`.
#' @seealso [simulate_greenmanure()], [decomp_fit()], [nutrient_release()]
"gm_litterbag"

#' Simulated two-season agronomic data
#'
#' Simulated field data for a control and three green manures at three
#' fertilizer N levels, in four blocks, for a kharif crop followed by a rabi
#' crop. Generated with [simulate_greenmanure()] (seed 42); the values are
#' synthetic.
#'
#' @format A data frame with 96 rows and 10 variables:
#' \describe{
#'   \item{treatment}{Control or green manure species.}
#'   \item{fert_N}{Fertilizer N applied (kg/ha).}
#'   \item{block}{Block (1-4).}
#'   \item{season}{Crop season (Kharif, Rabi).}
#'   \item{gm_biomass}{Green manure dry matter (t/ha).}
#'   \item{gm_N}{Green manure N concentration (%).}
#'   \item{grain_yield, straw_yield}{Yields (t/ha).}
#'   \item{grain_N, straw_N}{N concentrations (%).}
#' }
#' @source Simulated; see `data-raw/simulate-data.R`.
#' @seealso [simulate_greenmanure()], [greenmanure_n_balance()],
#'   [carryover_effect()]
"gm_agronomy"
