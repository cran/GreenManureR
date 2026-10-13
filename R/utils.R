# Internal helpers (not exported) ---------------------------------------------

.gm_key <- function(d, cols) {
  if (is.null(cols) || length(cols) == 0) return(rep("all", nrow(d)))
  do.call(paste, c(lapply(cols, function(g) as.character(d[[g]])), sep = " | "))
}

.need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this function. ",
         "Install it with install.packages(\"", pkg, "\").", call. = FALSE)
  }
}

.get <- function(d, name, default = 0) {
  if (is.null(name) || !name %in% names(d)) rep(default, nrow(d)) else as.numeric(d[[name]])
}

.check_cols <- function(d, cols, what = "data") {
  miss <- setdiff(cols, names(d))
  if (length(miss)) {
    stop("Column(s) not found in ", what, ": ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.with_seed <- function(seed, expr_fun) {
  if (is.null(seed)) expr_fun() else withr::with_seed(seed, expr_fun())
}

# Mass remaining (% of initial) for the three decay models
.mass_fun <- function(model, p, t) {
  switch(model,
    single     = 100 * exp(-p[["k"]] * t),
    double     = 100 * (p[["A"]] * exp(-p[["k1"]] * t) +
                        (1 - p[["A"]]) * exp(-p[["k1"]] * p[["r"]] * t)),
    asymptotic = p[["C"]] + (100 - p[["C"]]) * exp(-p[["k"]] * t),
    stop("Unknown model: ", model, call. = FALSE))
}

.valid_par <- function(model, P) {
  switch(model,
    single     = P[, "k"] > 0,
    double     = P[, "A"] > 0 & P[, "A"] < 1 & P[, "k1"] > 0 &
                 P[, "r"] > 0 & P[, "r"] <= 1,
    asymptotic = P[, "C"] >= 0 & P[, "C"] < 100 & P[, "k"] > 0)
}

# Time at which remaining mass reaches `target` (% of initial); NA if never
.t_at <- function(model, p, target, tmax) {
  f <- function(t) .mass_fun(model, p, t) - target
  if (!is.finite(f(0)) || f(0) <= 0) return(0)
  hi <- max(tmax, 1); n <- 0
  while (f(hi) > 0 && n < 30) { hi <- hi * 2; n <- n + 1 }
  if (f(hi) > 0) return(NA_real_)
  uniroot(f, c(0, hi), tol = 1e-8)$root
}

.k_start <- function(t, m) {
  ok <- t > 0 & m > 0
  if (sum(ok) < 2) return(0.02)
  k <- -unname(coef(lm(log(m[ok] / 100) ~ 0 + t[ok]))[1])
  if (!is.finite(k) || k <= 0) 0.02 else k
}

# Fit one model with a small multi-start grid; returns an nls object or NULL
.fit_one <- function(model, t, m) {
  k0  <- .k_start(t, m)
  ctl <- nls.control(maxiter = 500, warnOnly = TRUE, minFactor = 1e-12)
  dat <- list(t = t, m = m)
  best <- NULL
  try_fit <- function(formula, start, lower, upper) {
    f <- tryCatch(suppressWarnings(
      nls(formula, data = dat, start = start, lower = lower, upper = upper,
          algorithm = "port", control = ctl)), error = function(e) NULL)
    if (!is.null(f) && isTRUE(f$convInfo$isConv) &&
        (is.null(best) || deviance(f) < deviance(best))) best <<- f
  }
  if (model == "single") {
    for (k in k0 * c(0.5, 1, 2))
      try_fit(m ~ 100 * exp(-k * t), list(k = k), c(k = 1e-8), c(k = 10))
  } else if (model == "double") {
    for (A in c(0.3, 0.5, 0.7)) for (k1 in k0 * c(2, 4)) for (r in c(0.05, 0.2))
      try_fit(m ~ 100 * (A * exp(-k1 * t) + (1 - A) * exp(-k1 * r * t)),
              list(A = A, k1 = k1, r = r),
              c(A = 0.001, k1 = 1e-8, r = 1e-6), c(A = 0.999, k1 = 10, r = 1))
  } else if (model == "asymptotic") {
    for (C in c(5, 20, 40)) for (k in k0 * c(1, 3))
      try_fit(m ~ C + (100 - C) * exp(-k * t), list(C = C, k = k),
              c(C = 0, k = 1e-8), c(C = 99.9, k = 10))
  }
  best
}

# Mean / sd / se of `vars` within `groups`
.summ <- function(d, vars, groups) {
  key <- .gm_key(d, groups)
  u <- !duplicated(key)
  out <- d[u, groups, drop = FALSE]
  rownames(out) <- NULL
  for (v in vars) {
    x <- d[[v]]
    out[[paste0(v, "_mean")]] <- unname(tapply(x, key, mean, na.rm = TRUE)[key[u]])
    out[[paste0(v, "_sd")]]   <- unname(tapply(x, key, sd,   na.rm = TRUE)[key[u]])
    n <- unname(tapply(!is.na(x), key, sum)[key[u]])
    out[[paste0(v, "_se")]]   <- out[[paste0(v, "_sd")]] / sqrt(n)
  }
  out
}
