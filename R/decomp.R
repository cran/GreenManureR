#' Fit litterbag decomposition models
#'
#' Fits single-pool, double-pool and asymptotic decay models to litterbag
#' mass-loss data, separately for each group (for example each green manure
#' species), and ranks the models by AICc.
#'
#' With `M` the percentage of the initial mass remaining and `t` the time
#' since incorporation, the models are:
#'
#' * single: `M = 100 * exp(-k * t)`
#' * double: `M = 100 * (A * exp(-kf * t) + (1 - A) * exp(-ks * t))`, fitted
#'   with `ks = r * kf` and `0 < r <= 1`, so the fast pool can never decay
#'   more slowly than the slow pool.
#' * asymptotic: `M = C + (100 - C) * exp(-k * t)`, where `C` is the
#'   recalcitrant residue (percent of initial mass).
#'
#' The initial mass is fixed at 100%, so time-zero observations should be
#' included (or the data expressed relative to the initial mass). Each model
#' is fitted by nonlinear least squares (port algorithm with bounds) from a
#' small grid of starting values, and the best-fitting converged solution is
#' kept.
#'
#' @param data A data frame in long format with one row per litterbag and
#'   sampling date.
#' @param time Name of the column holding time since incorporation
#'   (e.g. days).
#' @param mass Name of the column holding the mass remaining.
#' @param group Character vector of grouping columns (species, site, ...), or
#'   `NULL` to fit a single group.
#' @param models Any of `"single"`, `"double"`, `"asymptotic"`.
#' @param mass_scale Either `"percent"` (default) or `"fraction"`, the scale
#'   of `mass`.
#'
#' @return An object of class `"gm_decomp"`, a list with components
#'   `table` (one row per group and model with parameters `k`, `k_slow`,
#'   `pool_A`, `C_resistant`, fit statistics `RSS`, `RMSE`, `R2`, `AIC`,
#'   `AICc`, `dAICc`, Akaike `weight`, and `selected` marking the lowest
#'   AICc), `fits` (the `nls` objects), `data`, `gcols` and `groups`. For the
#'   double-pool model `k` is the fast-pool rate and `k_slow` the slow-pool
#'   rate.
#'
#' @seealso [decomp_t50()], [nutrient_release()]
#' @examples
#' fit <- decomp_fit(gm_litterbag, time = "time", mass = "mass",
#'                   group = "treatment")
#' fit
#' fit$table[fit$table$selected, c("treatment", "model", "k", "R2")]
#' @export
decomp_fit <- function(data, time = "time", mass = "mass", group = NULL,
                       models = c("single", "double", "asymptotic"),
                       mass_scale = c("percent", "fraction")) {
  mass_scale <- match.arg(mass_scale)
  models <- match.arg(models, c("single", "double", "asymptotic"),
                      several.ok = TRUE)
  d <- as.data.frame(data)
  .check_cols(d, c(time, mass, group))
  if (is.null(group)) d$group <- "all"
  gcols <- if (is.null(group)) "group" else group
  d$.t <- as.numeric(d[[time]])
  d$.m <- as.numeric(d[[mass]]) * if (mass_scale == "fraction") 100 else 1
  d <- d[is.finite(d$.t) & is.finite(d$.m), , drop = FALSE]
  d$.key <- .gm_key(d, gcols)
  npar <- c(single = 1, double = 3, asymptotic = 2)

  rows <- list(); fits <- list(); groups <- list()
  for (key in unique(d$.key)) {
    dd <- d[d$.key == key, , drop = FALSE]
    gv <- dd[1, gcols, drop = FALSE]; rownames(gv) <- NULL
    groups[[key]] <- cbind(gv, .key = key)
    n <- nrow(dd); tss <- sum((dd$.m - mean(dd$.m))^2)
    for (mod in models) {
      if (n < npar[[mod]] + 2) next
      f <- .fit_one(mod, dd$.t, dd$.m)
      if (is.null(f)) next
      p <- coef(f); rss <- deviance(f); kk <- npar[[mod]] + 1
      row <- data.frame(gv, model = mod, n = n, npar = npar[[mod]],
        k = NA_real_, k_slow = NA_real_, pool_A = NA_real_,
        C_resistant = NA_real_, RSS = rss, RMSE = sqrt(rss / n),
        R2 = 1 - rss / tss, AIC = AIC(f),
        AICc = if (n - kk - 1 > 0) AIC(f) + 2 * kk * (kk + 1) / (n - kk - 1) else NA_real_,
        .key = key, stringsAsFactors = FALSE, row.names = NULL)
      if (mod == "single") row$k <- p[["k"]]
      if (mod == "asymptotic") { row$k <- p[["k"]]; row$C_resistant <- p[["C"]] }
      if (mod == "double") {
        row$k <- p[["k1"]]; row$k_slow <- p[["k1"]] * p[["r"]]
        row$pool_A <- p[["A"]]
      }
      rows[[length(rows) + 1]] <- row
      fits[[paste(key, mod, sep = "::")]] <- f
    }
  }
  if (!length(rows)) {
    stop("No model could be fitted; check the data (at least 4 points per group are needed).",
         call. = FALSE)
  }
  tb <- do.call(rbind, rows)
  tb$selected <- FALSE; tb$dAICc <- NA_real_; tb$weight <- NA_real_
  for (key in unique(tb$.key)) {
    i <- which(tb$.key == key)
    crit <- if (all(is.na(tb$AICc[i]))) tb$AIC[i] else tb$AICc[i]
    crit[is.na(crit)] <- Inf
    tb$dAICc[i] <- crit - min(crit)
    w <- exp(-0.5 * tb$dAICc[i]); tb$weight[i] <- w / sum(w)
    tb$selected[i[which.min(crit)]] <- TRUE
  }
  structure(list(table = tb, fits = fits, data = d, gcols = gcols,
                 groups = do.call(rbind, groups)), class = "gm_decomp")
}


#' Time to a given percentage of mass (or nutrient) loss
#'
#' Solves the fitted decay model for the time at which a proportion `p` of
#' the initial mass has been lost (`p = 0.5` gives t50), with an optional
#' Monte-Carlo confidence interval.
#'
#' The interval is obtained by drawing `nsim` parameter vectors from a
#' multivariate normal distribution centred on the estimates with the
#' covariance matrix of the `nls` fit, discarding draws outside the
#' admissible parameter space, and taking quantiles of the resulting times.
#' If the model never reaches the target (for example an asymptotic model
#' with a residue above the target), `NA` is returned.
#'
#' @param fit An object returned by [decomp_fit()].
#' @param p Proportion lost, strictly between 0 and 1 (0.5 gives t50, 0.95
#'   gives t95).
#' @param model `"best"` (lowest AICc, default) or the name of a model.
#' @param ci Logical; compute the Monte-Carlo interval?
#' @param nsim Number of Monte-Carlo draws.
#' @param level Confidence level.
#' @param seed Integer seed used for the draws. The global random number
#'   state is left unchanged. Use `NULL` to draw without fixing a seed.
#'
#' @return A data frame with the grouping columns, the model used, `p`, the
#'   time (column named e.g. `t50`) and the interval limits `lwr` and `upr`.
#' @seealso [decomp_fit()]
#' @examples
#' fit <- decomp_fit(gm_litterbag, group = "treatment")
#' decomp_t50(fit)
#' decomp_t50(fit, p = 0.95, ci = FALSE)
#' @export
decomp_t50 <- function(fit, p = 0.5, model = "best", ci = TRUE,
                       nsim = 2000, level = 0.95, seed = 123) {
  if (!inherits(fit, "gm_decomp")) {
    stop("`fit` must be an object returned by decomp_fit().", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
    stop("`p` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  tb <- fit$table; tmax <- max(fit$data$.t); target <- 100 * (1 - p)
  lab <- paste0("t", round(p * 100))

  run <- function() {
    out <- list()
    for (key in unique(tb$.key)) {
      tk <- tb[tb$.key == key, , drop = FALSE]
      mod <- if (model == "best") tk$model[tk$selected][1] else model
      f <- fit$fits[[paste(key, mod, sep = "::")]]
      if (is.null(f)) next
      est <- coef(f)
      tp <- .t_at(mod, est, target, tmax)
      lo <- hi <- NA_real_
      if (ci) {
        vc <- tryCatch(vcov(f), error = function(e) NULL)
        R <- if (!is.null(vc) && all(is.finite(vc))) {
          tryCatch(chol(vc + diag(1e-12, nrow(vc))), error = function(e) NULL)
        } else NULL
        if (!is.null(R)) {
          Z <- matrix(rnorm(nsim * length(est)), nsim) %*% R
          P <- sweep(Z, 2, est, "+"); colnames(P) <- names(est)
          P <- P[.valid_par(mod, P), , drop = FALSE]
          if (nrow(P) >= 50) {
            tv <- apply(P, 1, function(r) .t_at(mod, r, target, tmax))
            q <- quantile(tv, c((1 - level) / 2, 1 - (1 - level) / 2), na.rm = TRUE)
            lo <- unname(q[1]); hi <- unname(q[2])
          }
        }
      }
      row <- data.frame(tk[1, fit$gcols, drop = FALSE], model = mod, p = p,
                        value = tp, lwr = lo, upr = hi, row.names = NULL)
      names(row)[names(row) == "value"] <- lab
      out[[length(out) + 1]] <- row
    }
    do.call(rbind, out)
  }
  .with_seed(seed, run)
}
