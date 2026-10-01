#' Run multiple signatures or individual genes from one prepared matrix
#' @param prepared Prepared-data list or RDS path containing all required genes.
#' @param signature_sets Named list of gene-symbol vectors. A one-gene vector
#'   with min_signature_genes=1 gives single-gene survival analysis.
#' @param target_gene Optional target for the historical target-plus-signature model.
#' @param out_dir Output directory.
#' @param output_prefix Basename for combined results.
#' @param make_summary_plots Write per-signature forest plots and heatmaps.
#' @param ... Further parameters to run_signature_only_survival or
#'   run_signature_survival, such as min_signature_genes and make_km.
#' @return A list with combined results, per-signature analyses and failures.
#'   batch_FDR adjusts across all successful signature/cohort tests within
#'   cutoff_method and default_cohort; the historical per-signature FDR is retained.
#' @export
run_signature_batch <- function(prepared, signature_sets, target_gene = NULL,
    out_dir = file.path(getwd(), "results"), output_prefix = "panTCGA_batch",
    make_summary_plots = TRUE, ...) {
  if (is.character(prepared)) prepared <- readRDS(prepared)
  if (!is.list(signature_sets) || !length(signature_sets) || is.null(names(signature_sets)) ||
      anyNA(names(signature_sets)) || any(!nzchar(names(signature_sets))) ||
      anyDuplicated(vapply(names(signature_sets), tcgasig_safe_name, character(1)))) {
    stop("signature_sets requires unique, nonempty, safe names.", call. = FALSE)
  }
  dots <- list(...)
  if (any(names(dots) %in% c("prepared", "target_gene", "signature_genes", "signature_name", "out_dir", "output_prefix"))) {
    stop("Reserved batch arguments supplied in ...", call. = FALSE)
  }
  analyses <- list()
  failures <- list()
  rows <- list()
  for (name in names(signature_sets)) {
    args <- c(list(prepared = prepared, signature_genes = signature_sets[[name]],
      signature_name = name, out_dir = out_dir,
      output_prefix = paste0(tcgasig_safe_name(output_prefix), "_", tcgasig_safe_name(name))), dots)
    fn <- run_signature_only_survival
    if (!is.null(target_gene)) { fn <- run_signature_survival; args$target_gene <- target_gene }
    analysis <- tryCatch(do.call(fn, args), error = function(e) e)
    if (inherits(analysis, "error")) {
      failures[[length(failures) + 1L]] <- data.frame(signature_name = name, reason = conditionMessage(analysis))
      next
    }
    analyses[[name]] <- analysis
    if (nrow(analysis$results)) {
      row <- data.table::copy(analysis$results)
      row$batch_signature <- name
      rows[[name]] <- row
      if (make_summary_plots) {
        fn <- if (is.null(target_gene)) plot_signature_only_results else plot_signature_results
        plot_args <- list(results = row, signature_name = name, out_dir = out_dir, output_prefix = args$output_prefix)
        if (!is.null(target_gene)) plot_args$target_gene <- target_gene
        analyses[[name]]$plots <- do.call(fn, plot_args)
      }
    }
  }
  result <- if (length(rows)) data.table::rbindlist(rows, fill = TRUE) else data.table::data.table(batch_signature = character())
  if (nrow(result)) result[, batch_FDR := stats::p.adjust(wald_p, method = "BH"), by = .(cutoff_method, default_cohort)]
  failure_dt <- if (length(failures)) data.table::rbindlist(failures) else data.table::data.table(signature_name = character(), reason = character())
  if (!identical(dots$write_files, FALSE)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(result, file.path(out_dir, paste0(tcgasig_safe_name(output_prefix), "_results.tsv")), sep = "\t")
    data.table::fwrite(failure_dt, file.path(out_dir, paste0(tcgasig_safe_name(output_prefix), "_failures.tsv")), sep = "\t")
  }
  list(results = result, analyses = analyses, failures = failure_dt)
}

#' Prepare one TCGA matrix and run a batch of signatures
#' @inheritParams run_signature_batch
#' @param project_dir Project directory.
#' @param data_dir Public TCGA input directory.
#' @return The batch result list with a prepared element.
#' @export
run_pan_tcga_signature_batch <- function(signature_sets, target_gene = NULL,
    project_dir = getwd(), data_dir = file.path(project_dir, "data", "xena"),
    out_dir = file.path(project_dir, "results"), output_prefix = "panTCGA_batch", ...) {
  prepared <- prepare_signature_only_data(unique(c(target_gene, unlist(signature_sets, use.names = FALSE))),
    signature_name = "Batch", project_dir = project_dir, data_dir = data_dir,
    out_dir = out_dir, output_prefix = paste0(tcgasig_safe_name(output_prefix), "_input"), write_files = FALSE)
  result <- run_signature_batch(prepared, signature_sets, target_gene, out_dir, output_prefix, ...)
  result$prepared <- prepared
  result
}
