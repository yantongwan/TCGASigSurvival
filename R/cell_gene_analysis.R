#' Run the cell-gene expression, enrichment and survival workflow
#' @param data Prepared patient-level cell expression data.
#' @param target_gene,cell_type Exact gene and cell label.
#' @param pathways Named gene-set list, or NULL to skip enrichment.
#' @param enrichment_cohorts Cohorts for enrichment, NULL for all.
#' @param reference,paired Tumor/reference comparison settings.
#' @param endpoints Clinical outcomes. PFS requires its own measured columns.
#' @param covariates Clinical adjustment columns for DE and Cox.
#' @param context_subset all or high (an explicitly labelled bulk sensitivity).
#' @param allow_proxy Explicit permission for a bulk_proxy example.
#' @param out_dir,prefix Output location.
#' @param min_donors,min_patients,min_events,min_group Eligibility thresholds.
#' @param make_plots Write separate source-backed PDF/PNG panels.
#' @param score_mode gene or joint, the latter retaining the legacy scoring formula.
#' @return List of module results, module_status and provenance. Unsupported
#'   reference or outcome comparisons are reported, never relabelled as controls.
#' @export
run_cell_gene_analysis <- function(data, target_gene = "ACLY", cell_type = "Treg",
    pathways = NULL, enrichment_cohorts = NULL, reference = "healthy", paired = FALSE,
    endpoints = c("OS", "PFI", "DSS", "DFI"), covariates = character(),
    context_subset = "all", allow_proxy = FALSE,
    out_dir = file.path(getwd(), "results", "cell_gene"), prefix = "cell_gene",
    min_donors = 3L, min_patients = 40L, min_events = 10L, min_group = 10L,
    make_plots = TRUE, score_mode = c("gene", "joint")) {
  tcgasig_cell_check(data, allow_proxy)
  score_mode <- match.arg(score_mode)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  comparison <- compare_cell_gene_expression(data, target_gene, cell_type, reference,
    paired, min_donors, allow_proxy, out_dir, paste0(prefix, "_expression"), score_mode = score_mode)
  survival <- run_cell_gene_survival(data, target_gene, cell_type, endpoints, covariates,
    context_subset, min_patients, min_events, min_group, allow_proxy,
    out_dir, paste0(prefix, "_survival"), score_mode = score_mode)
  enrichment <- if (!is.null(pathways)) run_cell_gene_enrichment(data, target_gene,
    cell_type, pathways, cohorts = enrichment_cohorts, covariates = covariates,
    context_subset = context_subset, min_group = min_group, allow_proxy = allow_proxy,
    out_dir = out_dir, prefix = paste0(prefix, "_enrichment"), score_mode = score_mode) else NULL
  status <- data.frame(module = c("expression", "survival", "enrichment"),
    status = c(if (any(comparison$tests$status == "ok")) "completed_with_qc" else "reference_or_sample_blocked",
      if (any(survival$results$status == "ok")) "completed_with_qc" else "not_estimable",
      if (is.null(enrichment)) "not_requested" else if (nrow(enrichment$enrichment)) "completed_with_qc" else "not_estimable"),
    measurement = data$measurement)
  result <- list(comparison = comparison, survival = survival, enrichment = enrichment,
    module_status = status, provenance = data$provenance)
  if (make_plots) {
    save_plot <- function(p, name, width, height) {
      for (ext in c("pdf", "png")) ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_", name, ".", ext)),
        p, width = width, height = height, dpi = 300)
    }
    if (any(comparison$values$used_in_test)) save_plot(plot_cell_gene_expression(comparison), "expression", 10, 4)
    save_plot(plot_cell_gene_survival_grid(survival), "survival_median", 7, max(5, length(unique(survival$results$cohort_label)) * 0.2 + 2))
    save_plot(plot_cell_gene_survival_grid(survival, "continuous"), "survival_continuous", 7, max(5, length(unique(survival$results$cohort_label)) * 0.2 + 2))
    if (!is.null(enrichment) && nrow(enrichment$leading_edge)) {
      for (cohort in unique(enrichment$leading_edge$cohort_label)) {
        save_plot(plot_cell_gene_enrichment(enrichment, cohort), paste0(cohort, "_GSEA_ridges"), 10, 7)
      }
    }
  }
  tcgasig_write_cell_tables(result, out_dir, prefix)
  saveRDS(result, file.path(out_dir, paste0(prefix, "_analysis.rds")))
  writeLines(c(data$provenance, paste("Measurement:", data$measurement), utils::capture.output(utils::sessionInfo())),
    file.path(out_dir, paste0(prefix, "_provenance_session.txt")))
  result
}
