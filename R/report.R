#' Standardized tables, diagnostic plots and Word report
#'
#' Writes tables (CSV), figures (PNG) and, if the 'officer' package is
#' installed, a Word report containing all tables and figures for any
#' combination of analysis objects.
#'
#' Nothing is written unless `outdir` is supplied; the directory is created
#' if needed. Figures require 'ggplot2' and the Word report requires
#' 'officer'; when either is missing the corresponding output is skipped
#' with a message and the remaining outputs are still produced.
#'
#' @param decomp,nutrient,nbal,carry,joint Objects returned by
#'   [decomp_fit()], [nutrient_release()], [greenmanure_n_balance()],
#'   [carryover_effect()] and [gm_joint_analysis()]. Any may be `NULL`.
#' @param outdir Directory in which to write the files (required).
#' @param docx Logical; also write `GreenManureR_report.docx`?
#' @param dpi Resolution of the PNG figures.
#' @param digits Significant digits used in the exported tables.
#'
#' @return Invisibly, a list with `tables` (data frames) and `figures`
#'   (file paths).
#' @examples
#' \donttest{
#' fit <- decomp_fit(gm_litterbag, group = "treatment")
#' out <- file.path(tempdir(), "gm_report")
#' greenmanure_report(decomp = fit, outdir = out, docx = FALSE)
#' list.files(out)
#' }
#' @export
greenmanure_report <- function(decomp = NULL, nutrient = NULL, nbal = NULL,
                               carry = NULL, joint = NULL, outdir,
                               docx = TRUE, dpi = 300, digits = 3) {
  if (missing(outdir) || !is.character(outdir) || length(outdir) != 1) {
    stop("Please supply `outdir`, the directory in which to write the report files.",
         call. = FALSE)
  }
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  have_gg <- requireNamespace("ggplot2", quietly = TRUE)
  if (!have_gg) message("Package 'ggplot2' not installed: figures skipped.")
  tabs <- list(); figs <- list()
  rnd <- function(x) {
    n <- vapply(x, is.numeric, logical(1)); x[n] <- lapply(x[n], signif, digits); x
  }
  addtab <- function(nm, x, caption) {
    x <- rnd(x); tabs[[nm]] <<- list(x = x, cap = caption)
    write.csv(x, file.path(outdir, paste0(nm, ".csv")), row.names = FALSE)
  }
  addfig <- function(nm, p, caption, w = 7, h = 4.5) {
    if (!have_gg) return(invisible(NULL))
    fp <- file.path(outdir, paste0(nm, ".png"))
    ggplot2::ggsave(fp, p, width = w, height = h, dpi = dpi)
    figs[[nm]] <<- list(path = fp, cap = caption, w = w, h = h)
  }

  if (!is.null(decomp)) {
    addtab("T1_decomposition_models", decomp$table[, setdiff(names(decomp$table), ".key")],
           "Table 1. Decomposition model parameters, fit statistics and AICc ranking.")
    addtab("T2_t50", decomp_t50(decomp), "Table 2. Time to 50% mass loss (days) with 95% CI.")
    if (have_gg) {
      addfig("F1_decomposition_curves", plot(decomp), "Figure 1. Residue decomposition curves.")
      addfig("F2_decomposition_residuals", .plot_resid(decomp), "Figure 2. Residual diagnostics.")
    }
  }
  if (!is.null(nutrient)) {
    addtab("T3_nutrient_t50", nutrient$t50, "Table 3. Time to 50% nutrient release (days).")
    addtab("T4_nutrient_release_at", nutrient$release_at,
           "Table 4. Model-based cumulative nutrient release at selected times.")
    if (have_gg) addfig("F3_nutrient_release", plot(nutrient), "Figure 3. Cumulative nutrient release.", 8, 4.5)
  }
  if (!is.null(nbal)) {
    addtab("T5_N_balance", nbal$summary,
           "Table 5. N inputs, crop uptake and apparent balance (kg N/ha).")
    if (have_gg) addfig("F4_N_recovery", plot(nbal), "Figure 4. Apparent N recovery from green manure.")
  }
  if (!is.null(carry)) {
    addtab("T6_carryover_effects", carry$effects[, setdiff(names(carry$effects), "sig")],
           "Table 6. Residual effects of green manure on succeeding crops vs control.")
    addtab("T7_carryover_anova", carry$anova, "Table 7. Treatment F-tests within season.")
    if (!is.null(carry$pooled)) {
      addtab("T8_carryover_pooled", carry$pooled,
             "Table 8. Pooled analysis across seasons (treatment x season).")
    }
    if (have_gg) addfig("F5_carryover", plot(carry), "Figure 5. Yield response to green manure by season.", 8, 4.5)
  }
  if (!is.null(joint)) {
    addtab("T9_joint_table", joint$data,
           "Table 9. Integrated decomposition, release and response table.")
    if (!is.null(joint$cor)) {
      addtab("T10_joint_correlations",
             data.frame(variable = rownames(joint$cor), joint$cor, check.names = FALSE),
             "Table 10. Spearman rank correlations among integrated variables.")
    }
  }

  if (docx) {
    if (!requireNamespace("officer", quietly = TRUE)) {
      message("Package 'officer' not installed: Word report skipped (CSV and PNG files written).")
    } else {
      doc <- officer::read_docx()
      doc <- officer::body_add_par(doc, "GreenManureR analysis report", style = "heading 1")
      doc <- officer::body_add_par(doc, paste("Generated:", format(Sys.time(), "%d %B %Y")))
      for (t in tabs) {
        doc <- officer::body_add_par(doc, t$cap, style = "Normal")
        doc <- officer::body_add_table(doc, t$x)
        doc <- officer::body_add_par(doc, "")
      }
      for (f in figs) {
        doc <- officer::body_add_img(doc, f$path, width = min(f$w, 6.3),
                                     height = min(f$w, 6.3) * f$h / f$w)
        doc <- officer::body_add_par(doc, f$cap, style = "Normal")
      }
      print(doc, target = file.path(outdir, "GreenManureR_report.docx"))
    }
  }
  message("Report files written to: ", normalizePath(outdir))
  invisible(list(tables = lapply(tabs, `[[`, "x"), figures = lapply(figs, `[[`, "path")))
}
