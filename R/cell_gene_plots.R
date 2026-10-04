#' Plot tumor/reference donor expression
#' @param comparison Output of compare_cell_gene_expression.
#' @return ggplot object; only donors used in the stated test are plotted.
#' @export
plot_cell_gene_expression <- function(comparison) {
  d <- comparison$values[comparison$values$used_in_test, , drop = FALSE]
  if (!nrow(d)) stop("No eligible expression comparisons to plot.", call. = FALSE)
  d$tissue_status <- factor(d$tissue_status, c("healthy", "adjacent_normal", "tumor"))
  ggplot2::ggplot(d, ggplot2::aes(x = cohort_label, y = analysis_score, fill = tissue_status)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, position = ggplot2::position_dodge(0.75)) +
    ggplot2::geom_point(ggplot2::aes(color = tissue_status), alpha = 0.25, size = 0.6,
      position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75, seed = 42)) +
    ggplot2::scale_fill_manual(values = c(tumor = "#D95F59", healthy = "#57A89A", adjacent_normal = "#5F8BB4"), drop = TRUE) +
    ggplot2::scale_color_manual(values = c(tumor = "#D95F59", healthy = "#57A89A", adjacent_normal = "#5F8BB4"), guide = "none") +
    ggplot2::labs(x = "Cancer cohort", y = if (comparison$score_mode == "joint") paste0(comparison$target_gene, "/cell-context joint score") else "Normalized log expression", fill = "Specimen",
      title = paste(comparison$target_gene, "tumor/reference comparison"),
      subtitle = comparison$interpretation,
      caption = "Biological donors; source tables report n, paired/unpaired test, raw P and BH FDR.") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")
}

#' Plot a cancer-by-endpoint survival association grid
#' @param result Output of run_cell_gene_survival.
#' @param model median or continuous.
#' @return ggplot. Fill is log2 HR, size is capped -log10(global BH FDR),
#'   black outline indicates raw Wald P<0.05; crosses are not estimable.
#' @export
plot_cell_gene_survival_grid <- function(result, model = c("median", "continuous")) {
  model <- match.arg(model)
  d <- result$results[result$results$model == model, , drop = FALSE]
  d$cohort_label <- factor(d$cohort_label, rev(sort(unique(d$cohort_label))))
  d$endpoint <- factor(d$endpoint, result$settings$endpoints)
  d$log2_hr <- log2(d$HR)
  d$fdr_size <- pmin(10, -log10(pmax(d$FDR, 1e-10)))
  d$raw_p <- ifelse(d$wald_p < 0.05, "P < 0.05", "P >= 0.05")
  ok <- d[d$status == "ok", , drop = FALSE]
  missing <- d[d$status != "ok", , drop = FALSE]
  ggplot2::ggplot(d, ggplot2::aes(endpoint, cohort_label)) +
    ggplot2::geom_point(data = ok, ggplot2::aes(fill = log2_hr, size = fdr_size, color = raw_p), shape = 21, stroke = 0.7) +
    ggplot2::geom_point(data = missing, shape = 4, size = 1.8, color = "#777777") +
    ggplot2::scale_fill_gradient2(low = "#287C8E", mid = "white", high = "#C44941", midpoint = 0, name = "log2 HR") +
    ggplot2::scale_color_manual(values = c("P < 0.05" = "black", "P >= 0.05" = "#B5B5B5"), name = "Wald P") +
    ggplot2::scale_size_continuous(range = c(1.5, 5), limits = c(0, 10), name = "-log10 FDR") +
    ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = "Endpoint", y = "Cancer cohort", title = paste(result$settings$target_gene, model, "survival association"),
      subtitle = paste(result$interpretation, "|", result$settings$score_mode), caption = "HR > 1: higher hazard in the high-score group. Cross: not estimable, not nonsignificant.") +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Plot GSEA leading-edge statistic ridgelines
#' @param result Output of run_cell_gene_enrichment.
#' @param cohort One analyzed cohort.
#' @param top_n Top pathways by adjusted P, with alphabetical tie breaking.
#' @param significant_only Whether to show only raw P<0.05 pathways.
#' @return ggplot ridgeline of leading-edge moderated t statistics. This is
#'   not a distribution of GSEA enrichment scores or a Sankey/alluvial plot.
#' @export
plot_cell_gene_enrichment <- function(result, cohort, top_n = 20L, significant_only = FALSE) {
  tcgasig_need("ggridges")
  tcgasig_validate_threshold(top_n, "top_n")
  e <- result$enrichment[result$enrichment$cohort_label == cohort & result$enrichment$status == "ok", , drop = FALSE]
  if (significant_only) e <- e[is.finite(e$pval) & e$pval < 0.05, , drop = FALSE]
  e <- e[order(e$padj, e$pathway), , drop = FALSE]
  e <- utils::head(e, top_n)
  d <- result$leading_edge[result$leading_edge$cohort_label == cohort & result$leading_edge$pathway %in% e$pathway, , drop = FALSE]
  if (!nrow(d)) stop("No eligible leading-edge results.", call. = FALSE)
  labels <- paste0(sub("^HALLMARK_", "", e$pathway), "  (NES=", sprintf("%.2f", e$NES), ")")
  d$label <- factor(labels[match(d$pathway, e$pathway)], rev(labels))
  ggplot2::ggplot(d, ggplot2::aes(x = rank_statistic, y = label, fill = padj)) +
    ggridges::geom_density_ridges(scale = 0.8, rel_min_height = 0.01, color = "#333333", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
    ggplot2::scale_fill_gradient(low = "#59B5AD", high = "#E4C76B", name = "BH FDR") +
    ggplot2::labs(x = "Leading-edge moderated t (high minus low)", y = NULL,
      title = paste(cohort, "GSEA"), subtitle = result$interpretation,
      caption = "All tested genes were ranked. Ridge density is descriptive, not patient-level uncertainty.") +
    ggplot2::theme_classic(base_size = 10)
}

#' Plot endpoint-specific unadjusted Kaplan-Meier curves
#' @param result Output of run_cell_gene_survival.
#' @param cohort,endpoint One analyzed cohort and endpoint.
#' @return ggplot KM curves for the exact endpoint model population, showing
#'   unadjusted log-rank P. Cox adjustment does not make KM adjusted.
#' @export
plot_cell_gene_km <- function(result, cohort, endpoint = "OS") {
  d <- result$inputs[result$inputs$cohort_label == cohort & result$inputs$endpoint == endpoint &
    result$inputs$exclusion == "included", , drop = FALSE]
  r <- result$results[result$results$cohort_label == cohort & result$results$endpoint == endpoint &
    result$results$model == "median", , drop = FALSE]
  if (!nrow(d) || nrow(r) != 1L || r$status != "ok") stop("KM requested for a nonestimable model.", call. = FALSE)
  d$months <- d$endpoint_time / 30.4375
  fit <- survival::survfit(survival::Surv(months, endpoint_status) ~ group, data = d)
  lr <- survival::survdiff(survival::Surv(months, endpoint_status) ~ group, d)
  p <- stats::pchisq(lr$chisq, 1, lower.tail = FALSE)
  pplot <- survminer::ggsurvplot(fit, data = d, palette = c("#287C8E", "#C44941"),
    conf.int = TRUE, risk.table = FALSE,
    legend.title = if (result$settings$score_mode == "joint") paste(result$settings$target_gene, result$settings$cell_type, "joint score") else result$settings$target_gene,
    legend.labs = c(paste0("Low (n=", r$low_n, ")"), paste0("High (n=", r$high_n, ")")),
    xlab = "Time (months; days / 30.4375)", ylab = paste(endpoint, "probability"),
    ggtheme = ggplot2::theme_classic(base_size = 10))$plot
  pplot + ggplot2::labs(title = paste(cohort, endpoint), subtitle = result$interpretation,
    caption = paste0("Unadjusted log-rank P=", signif(p, 3), "; groups fixed before endpoint exclusions. Cox results are separate."))
}
