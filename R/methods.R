#' Methods for GreenManureR objects
#'
#' Print, prediction and plotting methods for the objects returned by
#' [decomp_fit()], [nutrient_release()], [greenmanure_n_balance()],
#' [carryover_effect()] and [gm_joint_analysis()]. Plot methods return
#' `ggplot` objects and require the 'ggplot2' package.
#'
#' @param x,object An object returned by one of the package functions.
#' @param y Unused; present for consistency with the [plot()] generic.
#' @param digits Number of significant digits to print.
#' @param times Times (e.g. days) at which to predict mass remaining.
#' @param model `"best"` (lowest AICc) or the name of a model.
#' @param ... Further arguments (ignored).
#'
#' @return `print` methods return `x` invisibly. `predict.gm_decomp()`
#'   returns a data frame with the grouping columns, `model`, `time` and
#'   `pred` (predicted mass remaining, %). `plot` methods return a `ggplot`
#'   object.
#' @examples
#' fit <- decomp_fit(gm_litterbag, group = "treatment")
#' predict(fit, times = c(30, 60, 90))
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(fit)
#' @name gm_methods
NULL

.print_df <- function(df, digits) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], signif, digits)
  print(df, row.names = FALSE)
}

#' @rdname gm_methods
#' @export
print.gm_decomp <- function(x, digits = 4, ...) {
  cat("<gm_decomp> ", length(unique(x$table$.key)), " group(s); models: ",
      paste(unique(x$table$model), collapse = ", "), "\n\n", sep = "")
  .print_df(x$table[, setdiff(names(x$table), c(".key", "AIC", "RSS")), drop = FALSE], digits)
  invisible(x)
}

#' @rdname gm_methods
#' @export
predict.gm_decomp <- function(object, times, model = "best", ...) {
  tb <- object$table; out <- list()
  for (key in unique(tb$.key)) {
    tk <- tb[tb$.key == key, , drop = FALSE]
    mod <- if (model == "best") tk$model[tk$selected][1] else model
    f <- object$fits[[paste(key, mod, sep = "::")]]
    if (is.null(f)) next
    out[[length(out) + 1]] <- data.frame(tk[1, object$gcols, drop = FALSE],
      model = mod, time = times, pred = .mass_fun(mod, coef(f), times),
      row.names = NULL)
  }
  do.call(rbind, out)
}

#' @rdname gm_methods
#' @export
print.gm_nutrient <- function(x, digits = 3, ...) {
  cat("<gm_nutrient> nutrients:", paste(names(x$fits), collapse = ", "), "\n\n")
  cat("Time to 50% release (best model):\n")
  .print_df(x$t50, digits)
  invisible(x)
}

#' @rdname gm_methods
#' @export
print.gm_nbalance <- function(x, digits = 3, ...) {
  cat("<gm_nbalance> treatment means (kg N/ha unless stated)\n\n")
  s <- x$summary
  .print_df(s[, c(x$by, x$treatment, "fert_N", grep("_mean$", names(s), value = TRUE))], digits)
  invisible(x)
}

#' @rdname gm_methods
#' @export
print.gm_carryover <- function(x, digits = 3, ...) {
  cat("<gm_carryover> response:", x$response, "| control:", x$control, "\n\n")
  e <- x$effects
  cols <- intersect(c("stratum", "fert_N", "season", "treatment", "mean_t",
                      "mean_c", "diff", "pct_change", "p_adj", "sig",
                      "retention_vs_first"), names(e))
  .print_df(e[, cols, drop = FALSE], digits)
  if (!is.null(x$pooled)) {
    cat("\nPooled analysis across seasons:\n")
    .print_df(x$pooled, digits)
  }
  invisible(x)
}

#' @rdname gm_methods
#' @export
print.gm_joint <- function(x, digits = 3, ...) {
  cat("<gm_joint> treatment-level linkage table\n\n")
  .print_df(x$data, digits)
  if (!is.null(x$cor)) {
    cat("\nSpearman correlations:\n")
    print(round(x$cor, 2))
  }
  invisible(x)
}

#' @rdname gm_methods
#' @export
plot.gm_decomp <- function(x, y, model = "best", ...) {
  .need("ggplot2")
  tt <- seq(0, max(x$data$.t), length.out = 200)
  pr <- predict(x, tt, model); pr$grp <- .gm_key(pr, x$gcols)
  ob <- x$data; ob$grp <- ob$.key
  ggplot2::ggplot() +
    ggplot2::geom_point(data = ob, ggplot2::aes(.t, .m, colour = grp), alpha = 0.6) +
    ggplot2::geom_line(data = pr, ggplot2::aes(time, pred, colour = grp), linewidth = 0.9) +
    ggplot2::labs(x = "Days after incorporation",
                  y = "Residue mass remaining (%)", colour = NULL) +
    ggplot2::theme_bw()
}

#' @rdname gm_methods
#' @export
plot.gm_nutrient <- function(x, y, ...) {
  .need("ggplot2")
  tt <- seq(0, max(x$data$time), length.out = 200); pr <- list()
  for (nu in names(x$fits)) { p <- predict(x$fits[[nu]], tt); p$nutrient <- nu; pr[[nu]] <- p }
  pr <- do.call(rbind, pr); pr$grp <- .gm_key(pr, x$gcols); pr$rel <- 100 - pr$pred
  ob <- x$data; ob$grp <- .gm_key(ob, x$gcols)
  ggplot2::ggplot() +
    ggplot2::geom_point(data = ob, ggplot2::aes(time, released, colour = grp), alpha = 0.5) +
    ggplot2::geom_line(data = pr, ggplot2::aes(time, rel, colour = grp), linewidth = 0.9) +
    ggplot2::facet_wrap(~nutrient) +
    ggplot2::labs(x = "Days after incorporation",
                  y = "Cumulative nutrient release (%)", colour = NULL) +
    ggplot2::theme_bw()
}

#' @rdname gm_methods
#' @export
plot.gm_nbalance <- function(x, y, ...) {
  .need("ggplot2")
  s <- x$summary; s <- s[!is.na(s$ANR_mean), ]
  s$fert <- factor(s$fert_N); s$trt <- s[[x$treatment]]
  s$season <- if (length(x$by)) s[[x$by[1]]] else "all"
  ggplot2::ggplot(s, ggplot2::aes(trt, ANR_mean, fill = fert)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.75) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ANR_mean - ANR_se, ymax = ANR_mean + ANR_se),
                           position = ggplot2::position_dodge(0.8), width = 0.25) +
    ggplot2::facet_wrap(~season) +
    ggplot2::labs(x = NULL, y = "Apparent N recovery from green manure (%)",
                  fill = "Fert. N (kg/ha)") +
    ggplot2::theme_bw()
}

#' @rdname gm_methods
#' @export
plot.gm_carryover <- function(x, y, ...) {
  .need("ggplot2")
  e <- x$effects
  e$strat <- if ("fert_N" %in% names(e)) paste0("Fert N ", e$fert_N) else "all"
  ggplot2::ggplot(e, ggplot2::aes(pct_change, treatment, colour = strat)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = pct_lwr, xmax = pct_upr),
                             position = ggplot2::position_dodge(0.6)) +
    ggplot2::facet_wrap(~season) +
    ggplot2::labs(x = paste0("Change in ", x$response, " vs ", x$control, " (%, 95% CI)"),
                  y = NULL, colour = NULL) +
    ggplot2::theme_bw()
}

.plot_resid <- function(fit) {
  .need("ggplot2")
  r <- list()
  for (key in unique(fit$table$.key)) {
    tk <- fit$table[fit$table$.key == key, ]; mod <- tk$model[tk$selected][1]
    f <- fit$fits[[paste(key, mod, sep = "::")]]
    r[[key]] <- data.frame(grp = key, fitted = fitted(f), resid = resid(f))
  }
  r <- do.call(rbind, r)
  ggplot2::ggplot(r, ggplot2::aes(fitted, resid, colour = grp)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(x = "Fitted mass remaining (%)", y = "Residual", colour = NULL,
                  title = "Residuals of selected decay models") +
    ggplot2::theme_bw()
}
