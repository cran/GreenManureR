#' Nutrient release from decomposing green manure residues
#'
#' Calculates the percentage of the initial nutrient content remaining in
#' litterbags, fits decay models to each nutrient and reports the time to 50%
#' release and model-based cumulative release at chosen times.
#'
#' Nutrient remaining (% of initial content) is
#' `mass * concentration / initial concentration`, where `mass` is the
#' percentage of the initial dry mass remaining and the initial concentration
#' is the group mean at the earliest sampling time (which must therefore be
#' present in `data`). Values above 100% indicate net immobilization.
#' Decay models are fitted with [decomp_fit()] to the nutrient-remaining
#' percentages.
#'
#' @param data Long data frame of litterbag observations including the
#'   initial (time zero) sampling.
#' @param time,mass Names of the time and mass-remaining (%) columns.
#' @param nutrients Names of the columns holding nutrient concentrations
#'   (% dry weight), e.g. `c("N", "P", "K")`.
#' @param group Character vector of grouping columns, or `NULL`.
#' @param biomass Optional. Name of a column, or a single number, giving the
#'   dry-matter input in t/ha. If supplied, nutrient input and release are
#'   also reported in kg/ha.
#' @param at Times at which model-based cumulative release is reported.
#' @param models Models passed to [decomp_fit()].
#'
#' @return An object of class `"gm_nutrient"`: a list with `data` (long table
#'   with `remaining`, `released` and optionally `input_kg_ha`,
#'   `released_kg_ha`), `fits` (a `gm_decomp` object per nutrient), `t50`,
#'   `release_at`, `gcols` and `at`.
#' @seealso [decomp_fit()], [decomp_t50()]
#' @examples
#' nr <- nutrient_release(gm_litterbag, nutrients = c("N", "P", "K"),
#'                        group = "treatment", biomass = "biomass_tha",
#'                        at = c(14, 30, 60))
#' nr
#' subset(nr$release_at, nutrient == "N")
#' @export
nutrient_release <- function(data, time = "time", mass = "mass",
                             nutrients = c("N", "P", "K"), group = NULL,
                             biomass = NULL, at = c(7, 14, 30, 60, 90),
                             models = c("single", "double", "asymptotic")) {
  d <- as.data.frame(data)
  .check_cols(d, c(time, mass, nutrients, group))
  if (is.character(biomass)) .check_cols(d, biomass)
  if (is.null(group)) d$group <- "all"
  gcols <- if (is.null(group)) "group" else group
  d$.t <- as.numeric(d[[time]]); d$.m <- as.numeric(d[[mass]])
  key <- .gm_key(d, gcols); t0 <- min(d$.t, na.rm = TRUE)
  bm <- if (is.null(biomass)) NULL else if (is.character(biomass)) {
    as.numeric(d[[biomass]])
  } else rep(as.numeric(biomass), nrow(d))

  long <- list()
  for (nu in nutrients) {
    conc <- as.numeric(d[[nu]]); at0 <- d$.t == t0
    C0 <- unname(tapply(conc[at0], key[at0], mean, na.rm = TRUE)[key])
    o <- d[, gcols, drop = FALSE]
    o$nutrient <- nu; o$time <- d$.t; o$mass <- d$.m
    o$conc <- conc; o$conc0 <- C0
    o$remaining <- d$.m * conc / C0
    o$released  <- 100 - o$remaining
    if (!is.null(bm)) {
      o$input_kg_ha <- bm * 1000 * C0 / 100
      o$released_kg_ha <- o$input_kg_ha * o$released / 100
    }
    long[[nu]] <- o
  }
  long <- do.call(rbind, long); rownames(long) <- NULL

  fits <- list(); t50 <- list(); rel <- list()
  for (nu in nutrients) {
    ln <- long[long$nutrient == nu, , drop = FALSE]
    fit <- decomp_fit(ln, time = "time", mass = "remaining", group = gcols,
                      models = models)
    fits[[nu]] <- fit
    t <- decomp_t50(fit); t$nutrient <- nu; t50[[nu]] <- t
    pr <- predict(fit, times = at); pr$nutrient <- nu
    pr$remaining_pct <- pr$pred; pr$released_pct <- 100 - pr$pred; pr$pred <- NULL
    if (!is.null(bm)) {
      inp <- tapply(ln$input_kg_ha, .gm_key(ln, gcols), mean, na.rm = TRUE)
      pr$released_kg_ha <- unname(inp[.gm_key(pr, gcols)]) * pr$released_pct / 100
    }
    rel[[nu]] <- pr
  }
  structure(list(data = long, fits = fits, t50 = do.call(rbind, t50),
                 release_at = do.call(rbind, rel), gcols = gcols, at = at),
            class = "gm_nutrient")
}


#' Nitrogen inputs, crop uptake and apparent balance
#'
#' Plot-level nitrogen accounting for green manure experiments, with
#' treatment summaries.
#'
#' All quantities are in kg N/ha unless stated otherwise:
#'
#' * `N_gm_input`: green manure dry matter (t/ha) x 1000 x N percent / 100.
#' * `N_uptake`: grain plus straw N removal.
#' * `N_balance`: (green manure N + fertilizer N) minus crop N uptake, an
#'   apparent balance (losses and soil N are not measured).
#' * `ANR`: apparent N recovery from green manure (percent): (uptake - mean
#'   control uptake at the same stratum and fertilizer N) / green manure N
#'   input x 100.
#' * `AE_gmN`: additional grain (kg) per kg of green manure N.
#' * `NFRV`: fertilizer N replacement value, the fertilizer N needed in the
#'   control plots (quadratic response fitted per `by` stratum, requiring at
#'   least three fertilizer levels) to match the green manure plot yield,
#'   minus the fertilizer N applied. `NA` when the yield lies outside the
#'   fitted response range.
#'
#' @param data Plot-level data frame.
#' @param treatment Column with the treatment label.
#' @param control Label of the treatment without green manure.
#' @param rep Column identifying replicates or blocks (not used in the
#'   calculations; kept for documentation and future use).
#' @param by Optional columns defining strata (e.g. season, site, year).
#' @param gm_biomass,gm_N Columns with green manure dry matter (t/ha) and N
#'   concentration (%).
#' @param fert_N Column with fertilizer N applied (kg/ha). If absent it is
#'   taken as zero.
#' @param grain_yield,grain_N,straw_yield,straw_N Columns with yields (t/ha)
#'   and N concentrations (%). Straw columns may be absent.
#' @param nfrv Logical; compute the fertilizer N replacement value?
#'
#' @return An object of class `"gm_nbalance"`: a list with `plot` (plot-level
#'   results), `summary` (means, standard deviations and standard errors per
#'   stratum, treatment and fertilizer level), `by` and `treatment`.
#' @seealso [carryover_effect()]
#' @examples
#' nb <- greenmanure_n_balance(gm_agronomy, by = "season")
#' nb
#' head(nb$summary[, c("season", "treatment", "fert_N", "ANR_mean")])
#' @export
greenmanure_n_balance <- function(data, treatment = "treatment",
                                  control = "Control", rep = "block",
                                  by = NULL, gm_biomass = "gm_biomass",
                                  gm_N = "gm_N", fert_N = "fert_N",
                                  grain_yield = "grain_yield",
                                  grain_N = "grain_N",
                                  straw_yield = "straw_yield",
                                  straw_N = "straw_N", nfrv = TRUE) {
  d <- as.data.frame(data)
  .check_cols(d, c(treatment, gm_biomass, gm_N, grain_yield, grain_N, by))
  d$.trt   <- as.character(d[[treatment]])
  if (!control %in% d$.trt) {
    stop("Control treatment '", control, "' not found in column '", treatment, "'.",
         call. = FALSE)
  }
  d$fert_N <- .get(d, fert_N)
  d$grain_yield <- .get(d, grain_yield, NA)
  d$N_gm_input  <- .get(d, gm_biomass) * 1000 * .get(d, gm_N) / 100
  d$N_fert      <- d$fert_N
  d$N_input_total <- d$N_gm_input + d$N_fert
  d$N_uptake_grain <- d$grain_yield * 1000 * .get(d, grain_N) / 100
  d$N_uptake_straw <- .get(d, straw_yield) * 1000 * .get(d, straw_N) / 100
  d$N_uptake <- d$N_uptake_grain + d$N_uptake_straw
  d$N_balance <- d$N_input_total - d$N_uptake
  d$NHI <- ifelse(d$N_uptake > 0, d$N_uptake_grain / d$N_uptake, NA)

  mk <- .gm_key(d, c(by, "fert_N")); isc <- d$.trt == control
  cu <- tapply(d$N_uptake[isc], mk[isc], mean, na.rm = TRUE)
  cy <- tapply(d$grain_yield[isc], mk[isc], mean, na.rm = TRUE)
  d$ctl_N_uptake <- unname(cu[mk]); d$ctl_yield <- unname(cy[mk])
  hasgm <- d$N_gm_input > 0
  d$ANR    <- ifelse(hasgm, (d$N_uptake - d$ctl_N_uptake) / d$N_gm_input * 100, NA)
  d$AE_gmN <- ifelse(hasgm, (d$grain_yield - d$ctl_yield) * 1000 / d$N_gm_input, NA)

  d$NFRV <- NA_real_
  if (nfrv) {
    sk <- .gm_key(d, by)
    for (s in unique(sk)) {
      ctl <- d[sk == s & isc, , drop = FALSE]
      if (length(unique(ctl$fert_N)) < 3 || sum(!is.na(ctl$grain_yield)) < 6) next
      m <- lm(grain_yield ~ fert_N + I(fert_N^2), data = ctl)
      nmax <- 2 * max(ctl$fert_N)
      pr <- function(n) predict(m, data.frame(fert_N = n))
      for (i in which(sk == s & hasgm & !is.na(d$grain_yield))) {
        y <- d$grain_yield[i]
        neq <- if (y <= pr(0)) 0 else if (pr(nmax) < y) NA_real_ else
          uniroot(function(n) pr(n) - y, c(0, nmax))$root
        d$NFRV[i] <- neq - d$fert_N[i]
      }
    }
  }
  vars <- c("N_gm_input", "N_fert", "N_uptake", "N_balance", "NHI",
            "grain_yield", "ANR", "AE_gmN", if (nfrv) "NFRV")
  structure(list(plot = d[, setdiff(names(d), ".trt"), drop = FALSE],
                 summary = .summ(d, vars, c(by, treatment, "fert_N")),
                 by = by, treatment = treatment),
            class = "gm_nbalance")
}
