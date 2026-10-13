#' Simulate litterbag and two-season agronomic data
#'
#' Generates synthetic data with the structure expected by the package
#' functions: a litterbag decomposition experiment with three green manure
#' species (Mucuna, sunhemp, dhaincha) and a replicated field experiment
#' with a control, three fertilizer N levels and two successive crops
#' (kharif and rabi). The values are simulated for testing and teaching and
#' do not represent real measurements.
#'
#' @param seed Integer seed. The global random number state is left
#'   unchanged. Use `NULL` to draw without fixing a seed.
#'
#' @return A list with data frames `litterbag` (see [gm_litterbag]) and
#'   `agronomy` (see [gm_agronomy]).
#' @examples
#' sim <- simulate_greenmanure(seed = 1)
#' str(sim$litterbag)
#' @export
simulate_greenmanure <- function(seed = 42) {
  .with_seed(seed, function() {
    sp <- c("Mucuna", "Sunhemp", "Dhaincha")
    P <- list(
      Mucuna   = list(A = .50, k1 = .055, r = .08, N0 = 2.9, P0 = .22, K0 = 1.8,
                      kN = .035, kP = .025, kK = .20, bm = 4.5),
      Sunhemp  = list(A = .40, k1 = .045, r = .06, N0 = 2.3, P0 = .19, K0 = 1.5,
                      kN = .028, kP = .020, kK = .18, bm = 6.0),
      Dhaincha = list(A = .55, k1 = .065, r = .07, N0 = 2.8, P0 = .24, K0 = 2.0,
                      kN = .045, kP = .030, kK = .22, bm = 5.5))
    times <- c(0, 7, 14, 21, 28, 42, 56, 84, 112)
    lit <- expand.grid(treatment = sp, rep = 1:4, time = times,
                       stringsAsFactors = FALSE)
    attr(lit, "out.attrs") <- NULL
    nz <- function(n, s) 1 + rnorm(n, 0, s)
    n <- nrow(lit)
    mass <- N <- P_ <- K <- bm <- numeric(n)
    for (i in seq_len(n)) {
      p <- P[[lit$treatment[i]]]; t <- lit$time[i]
      M <- 100 * (p$A * exp(-p$k1 * t) + (1 - p$A) * exp(-p$k1 * p$r * t)) *
        (if (t == 0) 1 else nz(1, .025))
      conc <- function(C0, k) C0 * (100 * exp(-k * t) *
                                      (if (t == 0) nz(1, .02) else nz(1, .03))) / M
      mass[i] <- M
      N[i] <- conc(p$N0, p$kN); P_[i] <- conc(p$P0, p$kP); K[i] <- conc(p$K0, p$kK)
      bm[i] <- p$bm
    }
    lit$mass <- mass; lit$N <- N; lit$P <- P_; lit$K <- K; lit$biomass_tha <- bm

    trts <- c("Control", sp); fert <- c(0, 60, 120)
    agr <- expand.grid(treatment = trts, fert_N = fert, block = 1:4,
                       season = c("Kharif", "Rabi"), stringsAsFactors = FALSE)
    attr(agr, "out.attrs") <- NULL
    na <- nrow(agr)
    bef <- setNames(rnorm(8, 0, .12), paste(rep(c("Kharif", "Rabi"), each = 4), 1:4))
    gme <- c(Control = 0, Mucuna = .85, Sunhemp = .70, Dhaincha = .80)
    bmv <- vapply(agr$treatment, function(s) if (s == "Control") 0 else P[[s]]$bm, 1)
    nv  <- vapply(agr$treatment, function(s) if (s == "Control") 0 else P[[s]]$N0, 1)
    agr$gm_biomass <- bmv * nz(na, .05)
    agr$gm_N <- nv
    Nf <- agr$fert_N; isk <- agr$season == "Kharif"
    base <- ifelse(isk, 3.0 + .020 * Nf - .00006 * Nf^2, 2.6 + .014 * Nf - .00004 * Nf^2)
    eff  <- gme[agr$treatment] * (1 - Nf / 240) * ifelse(isk, 1, .45)
    agr$grain_yield <- unname(base + eff + bef[paste(agr$season, agr$block)] +
                                rnorm(na, 0, .12))
    agr$straw_yield <- agr$grain_yield * 1.25 * nz(na, .04)
    agr$grain_N <- ifelse(isk, 1.10, 1.90) * nz(na, .04)
    agr$straw_N <- ifelse(isk, .45, .55) * nz(na, .05)
    list(litterbag = lit, agronomy = agr)
  })
}
