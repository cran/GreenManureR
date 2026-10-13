#' Residual (carry-over) effects of green manure on succeeding crops
#'
#' Quantifies the effect of each green manure treatment relative to a
#' control on the yield (or any response) of successive crops, for example a
#' kharif crop followed by a rabi crop.
#'
#' For every stratum and season the model `response ~ block + treatment` is
#' fitted (randomized complete block design; without a block term if
#' `block = NULL`). Treatment-versus-control contrasts use the model
#' estimates and standard errors, with p-values adjusted across treatments
#' by `adjust`. The percentage change is the contrast divided by the control
#' mean; the log response ratio uses the pooled residual variance in a
#' delta-method standard error. `retention_vs_first` expresses each seasonal
#' effect as a percentage of the effect in the first season.
#'
#' The pooled analysis fits `response ~ season:block + treatment +
#' treatment:season` across seasons; a significant treatment x season term
#' indicates that the effect changes between crops. Seasonal variances are
#' assumed homogeneous, which should be checked. The cumulative effect is the
#' sum of seasonal contrasts, with a standard error that assumes independence
#' between seasons and is therefore approximate.
#'
#' @param data Plot-level data frame.
#' @param response Column with the response (e.g. grain yield).
#' @param treatment Column with treatment labels.
#' @param control Label of the control treatment.
#' @param season Column identifying the crop or season.
#' @param block Column identifying blocks, or `NULL`.
#' @param strata Optional columns to analyse separately (e.g. `"fert_N"`,
#'   site).
#' @param season_levels Order of seasons in the crop sequence (defaults to
#'   order of appearance).
#' @param adjust Method passed to [stats::p.adjust()].
#' @param level Confidence level for intervals.
#'
#' @return An object of class `"gm_carryover"`: a list with `effects`
#'   (contrasts per stratum, season and treatment), `anova` (treatment
#'   F-tests within season), `pooled` (pooled analysis, or `NULL` with one
#'   season), `cumulative`, `response` and `control`.
#' @seealso [greenmanure_n_balance()], [gm_joint_analysis()]
#' @examples
#' co <- carryover_effect(gm_agronomy, response = "grain_yield",
#'                        strata = "fert_N",
#'                        season_levels = c("Kharif", "Rabi"))
#' co
#' @export
carryover_effect <- function(data, response = "grain_yield",
                             treatment = "treatment", control = "Control",
                             season = "season", block = "block",
                             strata = NULL, season_levels = NULL,
                             adjust = "holm", level = 0.95) {
  d <- as.data.frame(data)
  .check_cols(d, c(response, treatment, season, block, strata))
  d$.y <- as.numeric(d[[response]])
  d$.trt <- as.character(d[[treatment]])
  d$.season <- as.character(d[[season]])
  d$.blk <- if (is.null(block)) "1" else as.character(d[[block]])
  d <- d[is.finite(d$.y), , drop = FALSE]
  if (!control %in% d$.trt) {
    stop("Control treatment '", control, "' not found in column '", treatment, "'.",
         call. = FALSE)
  }
  if (is.null(season_levels)) season_levels <- unique(d$.season)
  d$.sk <- .gm_key(d, strata)
  a <- 1 - level
  eff <- list(); aov_l <- list(); pool <- list()

  for (sk in unique(d$.sk)) {
    ds <- d[d$.sk == sk, , drop = FALSE]
    sv <- if (is.null(strata)) data.frame(stratum = "all") else ds[1, strata, drop = FALSE]
    rownames(sv) <- NULL
    for (se in season_levels) {
      dd <- ds[ds$.season == se, , drop = FALSE]
      if (!control %in% dd$.trt || length(unique(dd$.trt)) < 2) next
      dd$trt <- factor(dd$.trt, levels = c(control, setdiff(sort(unique(dd$.trt)), control)))
      dd$blk <- factor(dd$.blk)
      f <- if (nlevels(dd$blk) > 1) lm(.y ~ blk + trt, dd) else lm(.y ~ trt, dd)
      dfr <- f$df.residual; mse <- sum(resid(f)^2) / dfr
      an <- anova(f)
      aov_l[[length(aov_l) + 1]] <- data.frame(sv, season = se,
        F = an["trt", "F value"], df1 = an["trt", "Df"], df2 = dfr,
        p = an["trt", "Pr(>F)"], row.names = NULL)
      cf <- summary(f)$coefficients
      cf <- cf[grepl("^trt", rownames(cf)), , drop = FALSE]
      mc <- mean(dd$.y[dd$.trt == control]); nc <- sum(dd$.trt == control)
      tq <- qt(1 - a / 2, dfr)
      e <- data.frame(sv, season = se, treatment = sub("^trt", "", rownames(cf)),
                      row.names = NULL, stringsAsFactors = FALSE)
      e$n_t <- as.integer(table(dd$.trt)[e$treatment]); e$n_c <- nc
      e$mean_t <- as.numeric(tapply(dd$.y, dd$.trt, mean)[e$treatment]); e$mean_c <- mc
      e$diff <- cf[, "Estimate"]; e$se <- cf[, "Std. Error"]
      e$lwr <- e$diff - tq * e$se; e$upr <- e$diff + tq * e$se
      e$t <- cf[, "t value"]; e$df <- dfr; e$p <- cf[, "Pr(>|t|)"]
      e$p_adj <- p.adjust(e$p, adjust)
      e$pct_change <- e$diff / mc * 100
      e$pct_lwr <- e$lwr / mc * 100; e$pct_upr <- e$upr / mc * 100
      e$lnRR <- log(e$mean_t / mc)
      e$lnRR_se <- sqrt(mse / (e$n_t * e$mean_t^2) + mse / (nc * mc^2))
      e$sig <- cut(e$p_adj, c(-Inf, .001, .01, .05, .1, Inf),
                   c("***", "**", "*", ".", "ns"))
      eff[[length(eff) + 1]] <- e
    }
    dp <- ds[ds$.season %in% season_levels, , drop = FALSE]
    if (length(unique(dp$.season)) > 1 && control %in% dp$.trt) {
      dp$trt <- factor(dp$.trt, levels = c(control, setdiff(sort(unique(dp$.trt)), control)))
      dp$sb <- interaction(dp$.season, dp$.blk, drop = TRUE)
      dp$s <- factor(dp$.season)
      an <- anova(lm(.y ~ sb + trt + trt:s, dp))
      pool[[length(pool) + 1]] <- data.frame(sv,
        term = c("treatment", "treatment x season"),
        F = an[c("trt", "trt:s"), "F value"], df1 = an[c("trt", "trt:s"), "Df"],
        df2 = an["Residuals", "Df"], p = an[c("trt", "trt:s"), "Pr(>F)"],
        row.names = NULL)
    }
  }
  if (!length(eff)) stop("No season contained both the control and another treatment.", call. = FALSE)
  eff <- do.call(rbind, eff); rownames(eff) <- NULL
  gk <- c(if (is.null(strata)) "stratum" else strata, "treatment")
  eff$season <- factor(eff$season, levels = season_levels)
  ek <- .gm_key(eff, gk); first <- season_levels[1]
  d1 <- eff$diff[eff$season == first]; names(d1) <- ek[eff$season == first]
  eff$retention_vs_first <- 100 * eff$diff / unname(d1[ek])

  cum <- do.call(rbind, lapply(split(eff, ek), function(x) {
    data.frame(x[1, gk, drop = FALSE], seasons = paste(x$season, collapse = " + "),
               cum_diff = sum(x$diff), cum_se = sqrt(sum(x$se^2)), row.names = NULL)
  }))
  rownames(cum) <- NULL
  structure(list(effects = eff, anova = do.call(rbind, aov_l),
                 pooled = if (length(pool)) do.call(rbind, pool) else NULL,
                 cumulative = cum, response = response, control = control),
            class = "gm_carryover")
}


#' Joint analysis of decomposition, nutrient release and crop response
#'
#' Builds one row per treatment that joins the decomposition rate and t50,
#' the cumulative release of one nutrient by a chosen day, apparent N
#' recovery and the percentage yield response in each season, and reports
#' Spearman rank correlations among these variables.
#'
#' With only a few treatments (for example three green manure species) the
#' correlations are descriptive only, and a warning is issued when fewer than
#' five rows are available. To draw inference, pool several sites or years.
#'
#' @param decomp A `gm_decomp` object from [decomp_fit()] (one grouping
#'   column named as `key`).
#' @param nutrient A `gm_nutrient` object from [nutrient_release()].
#' @param nbal A `gm_nbalance` object from [greenmanure_n_balance()].
#' @param carry A `gm_carryover` object from [carryover_effect()].
#' @param key Name of the treatment column shared by all inputs.
#' @param nutrient_name Nutrient (as named in `nutrient_release()`) whose
#'   release is linked to crop response.
#' @param at_day Day at which model-based cumulative release is evaluated.
#' @param anr_season,fert_N_level Optional filters restricting the plots used
#'   to average apparent N recovery: the first `by` stratum of `nbal` and the
#'   fertilizer N level.
#'
#' @return An object of class `"gm_joint"`: a list with `data` (the linkage
#'   table) and `cor` (Spearman correlation matrix, or `NULL`).
#' @seealso [decomp_fit()], [nutrient_release()], [greenmanure_n_balance()],
#'   [carryover_effect()]
#' @examples
#' fit <- decomp_fit(gm_litterbag, group = "treatment")
#' nr  <- nutrient_release(gm_litterbag, group = "treatment",
#'                         biomass = "biomass_tha")
#' nb  <- greenmanure_n_balance(gm_agronomy, by = "season")
#' co  <- carryover_effect(gm_agronomy, strata = "fert_N")
#' ja  <- suppressWarnings(
#'   gm_joint_analysis(fit, nr, nb, co, at_day = 60, fert_N_level = 0))
#' ja$data
#' @export
gm_joint_analysis <- function(decomp, nutrient, nbal, carry,
                              key = "treatment", nutrient_name = "N",
                              at_day = 60, anr_season = NULL,
                              fert_N_level = NULL) {
  stopifnot(inherits(decomp, "gm_decomp"), inherits(nutrient, "gm_nutrient"),
            inherits(nbal, "gm_nbalance"), inherits(carry, "gm_carryover"))
  if (length(decomp$gcols) != 1 || length(nutrient$gcols) != 1) {
    stop("decomp_fit() and nutrient_release() must use a single grouping column.",
         call. = FALSE)
  }
  if (!nutrient_name %in% names(nutrient$fits)) {
    stop("Nutrient '", nutrient_name, "' not found in `nutrient`.", call. = FALSE)
  }
  tb <- decomp$table[decomp$table$selected, ]
  t50 <- decomp_t50(decomp, ci = FALSE)
  J <- merge(tb[, c(decomp$gcols, "model", "k")], t50[, c(decomp$gcols, "t50")],
             by = decomp$gcols)
  names(J)[names(J) == "k"] <- "k_mass"; names(J)[1] <- key

  pr <- predict(nutrient$fits[[nutrient_name]], times = at_day)
  pr$released <- 100 - pr$pred
  rel <- data.frame(pr[[nutrient$gcols]], pr$released)
  names(rel) <- c(key, paste0(nutrient_name, "_release_d", at_day))
  J <- merge(J, rel, by = key)

  nb <- nbal$plot
  if (!is.null(fert_N_level)) nb <- nb[nb$fert_N == fert_N_level, ]
  if (!is.null(anr_season) && !is.null(nbal$by)) nb <- nb[nb[[nbal$by[1]]] == anr_season, ]
  nb <- nb[nb[[nbal$treatment]] != carry$control & !is.na(nb$ANR), ]
  A <- aggregate(nb[, c("N_gm_input", "ANR", "NFRV")], list(nb[[nbal$treatment]]),
                 mean, na.rm = TRUE)
  names(A)[1] <- key; J <- merge(J, A, by = key, all.x = TRUE)

  e <- carry$effects
  R <- aggregate(e$pct_change, list(e$treatment, e$season), mean)
  names(R) <- c(key, "season", "pct")
  Rw <- reshape(R, idvar = key, timevar = "season", direction = "wide")
  names(Rw) <- sub("^pct\\.", "yield_gain_pct_", names(Rw))
  J <- merge(J, Rw, by = key, all.x = TRUE)

  num <- names(J)[vapply(J, is.numeric, logical(1))]
  cor_m <- if (nrow(J) >= 3) {
    suppressWarnings(cor(J[, num], method = "spearman", use = "pairwise.complete.obs"))
  } else NULL
  if (nrow(J) < 5) {
    warning("Only ", nrow(J), " treatment rows: correlations are descriptive. ",
            "Pool sites or years for inference.", call. = FALSE)
  }
  structure(list(data = J, cor = cor_m), class = "gm_joint")
}
