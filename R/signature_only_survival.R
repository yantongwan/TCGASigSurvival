tcgasig_signature_only_output_prefix <- function(signature_name, output_prefix = NULL) {
  if (!is.null(output_prefix) && nzchar(output_prefix)) return(output_prefix)
  paste0("panTCGA_", tcgasig_safe_name(signature_name), "_signature_only")
}

prepare_signature_only_data <- function(
    signature_genes,
    signature_name = "Signature",
    project_dir = getwd(),
    data_dir = file.path(project_dir, "data", "xena"),
    out_dir = file.path(project_dir, "results"),
    expression_file = file.path(data_dir, "tcga_RSEM_gene_tpm.gz"),
    survival_file = file.path(data_dir, "Survival_SupplementalTable_S1_20171025_xena_sp.tsv"),
    phenotype_file = file.path(data_dir, "TcgaTargetGTEX_phenotype.txt.gz"),
    gencode_file = file.path(data_dir, "gencode.v23.annotation.gtf.gz"),
    output_prefix = NULL,
    write_files = TRUE,
    verbose = TRUE) {
  signature_genes <- unique(tcgasig_clean_missing(signature_genes))
  signature_genes <- signature_genes[!is.na(signature_genes)]
  if (length(signature_genes) == 0) tcgasig_stop("signature_genes must contain at least one gene symbol.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  output_prefix <- tcgasig_signature_only_output_prefix(signature_name, output_prefix)

  for (path in c(expression_file, survival_file, phenotype_file, gencode_file)) {
    if (!file.exists(path)) {
      tcgasig_stop("Missing input file: ", path, ". Download UCSC Xena files first.")
    }
  }

  tcgasig_message(verbose, "Preparing pan-TCGA ", signature_name, " signature-only data.")
  tcgasig_message(verbose, "Interpretation: bulk RNA-seq ", signature_name,
                  " signature high/low survival; no target gene is added.")

  gene_map <- tcgasig_parse_gencode_gene_map(gencode_file, signature_genes, verbose = verbose)
  missing_from_gencode <- setdiff(signature_genes, gene_map$gene_name)

  selected_expr <- tcgasig_read_selected_expression(expression_file, gene_map$gene_id, verbose = verbose)
  selected_expr <- merge(gene_map[, .(gene_id, gene_name, gene_type)], selected_expr,
                         by = "gene_id", all.y = TRUE, sort = FALSE)
  sample_cols <- setdiff(names(selected_expr), c("gene_id", "gene_name", "gene_type"))
  if (length(sample_cols) == 0) tcgasig_stop("Selected expression matrix contains no sample columns.")

  transformed <- tcgasig_detect_and_transform_expression(selected_expr, sample_cols)
  collapsed_expr <- tcgasig_collapse_duplicate_symbols(transformed$expr, sample_cols)

  present_symbols <- collapsed_expr$gene_name
  missing_symbols <- setdiff(signature_genes, present_symbols)
  if (length(intersect(signature_genes, present_symbols)) == 0) {
    tcgasig_stop("None of the requested signature genes were found after expression extraction.")
  }

  survival <- tcgasig_normalize_survival(survival_file)
  phenotype <- tcgasig_normalize_phenotype(phenotype_file)
  sample_metadata <- tcgasig_build_sample_metadata(sample_cols, survival, phenotype)
  cohort_samples <- tcgasig_build_cohort_samples(sample_metadata)

  prepared <- list(
    signature_only = TRUE,
    signature_name = signature_name,
    signature_genes = signature_genes,
    expression = collapsed_expr,
    sample_metadata = sample_metadata,
    cohort_samples = cohort_samples,
    expression_scale = transformed$scale,
    expression_scale_note = transformed$note,
    expression_summary = data.table::data.table(
      metric = c("sampled_min", "sampled_max", "sampled_q99"),
      value = c(transformed$min, transformed$max, transformed$q99)
    ),
    missing_from_gencode = missing_from_gencode,
    missing_from_expression = missing_symbols,
    output_prefix = output_prefix
  )

  if (isTRUE(write_files)) {
    data.table::fwrite(gene_map, file.path(out_dir, paste0(output_prefix, "_gene_annotation.tsv")), sep = "\t")
    data.table::fwrite(collapsed_expr, file.path(out_dir, paste0(output_prefix, "_selected_log_expression.tsv")), sep = "\t")
    data.table::fwrite(sample_metadata, file.path(out_dir, paste0(output_prefix, "_sample_metadata.tsv")), sep = "\t")
    data.table::fwrite(cohort_samples, file.path(out_dir, paste0(output_prefix, "_cohort_samples.tsv")), sep = "\t")
    saveRDS(prepared, file.path(out_dir, paste0(output_prefix, "_prepared_data.rds")))

    qc <- data.table::data.table(
      metric = c(
        "signature_name",
        "requested_signature_genes",
        "signature_genes_present_after_expression_collapse",
        "signature_genes_missing_after_expression_collapse",
        "expression_scale",
        "expression_scale_note",
        "expression_rows_after_symbol_collapse",
        "expression_samples",
        "tcga_samples_with_os",
        "cohort_sample_rows",
        "default_cohorts",
        "skcm_modes_created"
      ),
      value = c(
        signature_name,
        paste(signature_genes, collapse = ","),
        paste(intersect(signature_genes, present_symbols), collapse = ","),
        ifelse(length(setdiff(signature_genes, present_symbols)) == 0, "none", paste(setdiff(signature_genes, present_symbols), collapse = ",")),
        transformed$scale,
        transformed$note,
        as.character(nrow(collapsed_expr)),
        as.character(length(sample_cols)),
        as.character(nrow(sample_metadata[has_os == TRUE & is_tcga_barcode == TRUE])),
        as.character(nrow(cohort_samples)),
        paste(sort(unique(cohort_samples[default_cohort == TRUE, cohort_label])), collapse = ","),
        paste(sort(unique(cohort_samples[cancer_type == "SKCM", sample_mode])), collapse = ",")
      )
    )
    data.table::fwrite(qc, file.path(out_dir, paste0(output_prefix, "_preparation_QC.tsv")), sep = "\t")
  }

  tcgasig_message(verbose, "Prepared signature expression rows: ", nrow(collapsed_expr),
                  "; cohort sample rows: ", nrow(cohort_samples))
  tcgasig_message(verbose, "Expression scale: ", transformed$scale, " | ", transformed$note)
  prepared
}

tcgasig_plot_signature_only_km <- function(dt, cancer_label, sample_mode, cutoff_method,
                                           cutoff_value, cox_row, out_file,
                                           signature_name) {
  plot_dt <- data.table::copy(dt)
  group_counts <- plot_dt[, .(
    n = .N,
    events = sum(OS.status == 1L, na.rm = TRUE)
  ), by = group]
  group_count_value <- function(group_name, column_name) {
    value <- group_counts[group == group_name, get(column_name)]
    if (length(value) == 0 || is.na(value[1])) 0L else as.integer(value[1])
  }
  low_label <- sprintf("Low n=%d, events=%d", group_count_value("low", "n"), group_count_value("low", "events"))
  high_label <- sprintf("High n=%d, events=%d", group_count_value("high", "n"), group_count_value("high", "events"))
  plot_dt[, group_label := factor(group, levels = c("low", "high"), labels = c(low_label, high_label))]
  fit <- survival::survfit(survival::Surv(OS.time, OS.status) ~ group_label, data = plot_dt)
  annotation <- sprintf(
    "Cox Wald P = %s\nHR high vs low = %.2f (95%% CI %.2f-%.2f)\nCutoff = %.3f",
    tcgasig_format_p(cox_row$wald_p),
    cox_row$HR,
    cox_row$CI_lower,
    cox_row$CI_upper,
    cutoff_value
  )
  title <- sprintf(
    "TCGA %s %s signature high/low (%s, %s)",
    cancer_label,
    signature_name,
    cutoff_method,
    sample_mode
  )
  g <- survminer::ggsurvplot(
    fit,
    data = plot_dt,
    conf.int = FALSE,
    risk.table = FALSE,
    pval = annotation,
    pval.coord = c(0, 0.13),
    legend.title = "",
    legend.labs = c(low_label, high_label),
    xlab = "Overall survival time (days)",
    ylab = "Survival probability",
    title = title,
    palette = c("#2C7FB8", "#D95F0E"),
    ggtheme = ggplot2::theme_classic(base_size = 12)
  )
  g$plot <- g$plot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "top"
    )
  tcgasig_open_pdf(out_file, width = 7, height = 5.5)
  print(g$plot)
  grDevices::dev.off()
  invisible(out_file)
}

run_signature_only_survival <- function(
    prepared,
    signature_genes = NULL,
    signature_name = NULL,
    out_dir = file.path(getwd(), "results"),
    output_prefix = NULL,
    min_signature_genes = 6,
    min_patients = 80,
    min_events = 20,
    include_skcm_primary_only = TRUE,
    cutoff_methods = c("optimized", "median"),
    optimized_minprop = 0.10,
    make_km = TRUE,
    write_files = TRUE,
    verbose = TRUE) {
  if (is.character(prepared) && length(prepared) == 1L) prepared <- readRDS(prepared)
  if (!is.list(prepared)) tcgasig_stop("prepared must be a prepared-data list or an RDS path.")

  signature_name <- if (is.null(signature_name)) tcgasig_prepared_value(prepared, "signature_name", "Signature") else signature_name
  signature_genes <- if (is.null(signature_genes)) {
    tcgasig_prepared_value(prepared, "signature_genes")
  } else {
    signature_genes
  }
  signature_genes <- unique(tcgasig_clean_missing(signature_genes))
  signature_genes <- signature_genes[!is.na(signature_genes)]
  if (length(signature_genes) == 0) tcgasig_stop("signature_genes must contain at least one gene symbol.")

  output_prefix <- tcgasig_signature_only_output_prefix(signature_name, output_prefix)
  cutoff_methods <- match.arg(cutoff_methods, choices = c("optimized", "median"), several.ok = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  km_dir <- file.path(out_dir, "KM_plots")
  if (isTRUE(make_km)) dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)

  expr <- data.table::as.data.table(prepared$expression)
  cohort_samples <- data.table::as.data.table(prepared$cohort_samples)
  if (!all(c("gene_name", "gene_id") %in% names(expr))) {
    tcgasig_stop("Prepared expression table must contain gene_name and gene_id.")
  }
  if (!all(c("cohort_label", "patient_id", "sample_barcode", "cancer_type", "sample_mode", "OS.time", "OS.status") %in% names(cohort_samples))) {
    tcgasig_stop("Prepared cohort sample table lacks required columns.")
  }
  if (!isTRUE(include_skcm_primary_only)) {
    cohort_samples <- cohort_samples[!(cancer_type == "SKCM" & sample_mode == "primary_only" & default_cohort == FALSE)]
  }

  tcgasig_message(verbose, "Running pan-TCGA ", signature_name, " signature-only survival analysis.")
  tcgasig_message(verbose, "Minimum patients: ", min_patients, "; minimum OS death events: ", min_events)
  tcgasig_message(verbose, "Interpretation: bulk RNA-seq ", signature_name,
                  " signature high/low survival; no target gene is added.")

  cohort_defs <- unique(cohort_samples[, .(cohort_label, cancer_type, sample_mode, default_cohort)])
  data.table::setorder(cohort_defs, cancer_type, sample_mode)

  all_results <- list()
  all_inputs <- list()
  skip_rows <- list()
  signature_safe <- tcgasig_safe_name(signature_name)

  for (i in seq_len(nrow(cohort_defs))) {
    cohort_label <- cohort_defs$cohort_label[i]
    cancer_type <- cohort_defs$cancer_type[i]
    sample_mode <- cohort_defs$sample_mode[i]
    default_cohort <- cohort_defs$default_cohort[i]
    tcgasig_message(verbose, "Cohort ", cohort_label, " (", sample_mode, ")")

    cs <- data.table::copy(cohort_samples[cohort_label == cohort_defs$cohort_label[i]])
    available_samples <- intersect(cs$sample_barcode, names(expr))
    cs <- cs[sample_barcode %in% available_samples]
    if (nrow(cs) == 0) {
      skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
        cohort_label, cancer_type, sample_mode, default_cohort,
        reason = "no_expression_samples_after_matching",
        n_patients = 0L,
        n_events = 0L
      )
      next
    }

    data.table::setorder(cs, patient_id, sample_barcode)
    cs <- cs[, .SD[1], by = patient_id]
    sample_cols <- cs$sample_barcode

    signature_available_global <- intersect(signature_genes, expr$gene_name)
    signature_z_list <- list()
    zero_var_signature <- character()
    for (gene in signature_available_global) {
      vals <- tcgasig_extract_gene_vector(expr, gene, sample_cols)
      z <- tcgasig_zscore_or_na(vals)
      if (all(is.na(z))) {
        zero_var_signature <- c(zero_var_signature, gene)
      } else {
        signature_z_list[[gene]] <- z
      }
    }
    usable_signature <- names(signature_z_list)
    missing_signature <- sort(unique(c(setdiff(signature_genes, signature_available_global), zero_var_signature)))
    if (length(usable_signature) < min_signature_genes) {
      skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
        cohort_label, cancer_type, sample_mode, default_cohort,
        reason = paste0("fewer_than_", min_signature_genes, "_usable_signature_genes:", paste(usable_signature, collapse = ",")),
        n_patients = length(unique(cs$patient_id)),
        n_events = sum(cs$OS.status == 1L, na.rm = TRUE)
      )
      next
    }

    signature_mat <- do.call(cbind, signature_z_list)
    signature_score <- rowMeans(signature_mat, na.rm = TRUE)

    analysis <- data.table::data.table(
      cohort_label = cohort_label,
      cancer_type = cancer_type,
      sample_mode = sample_mode,
      default_cohort = default_cohort,
      patient_id = cs$patient_id,
      sample_barcode = cs$sample_barcode,
      OS.time = as.numeric(cs$OS.time),
      OS.status = as.integer(cs$OS.status),
      signature_name = signature_name,
      signature_score = as.numeric(signature_score),
      signature_genes_present = paste(usable_signature, collapse = ","),
      signature_genes_missing = ifelse(length(missing_signature) == 0, "none", paste(missing_signature, collapse = ","))
    )
    analysis[, signature_z := tcgasig_zscore_or_na(signature_score)]
    analysis[, signature_only_score := signature_score]
    analysis[, target_signature_score := signature_score]
    analysis <- analysis[
      is.finite(OS.time) &
        OS.time > 0 &
        OS.status %in% c(0L, 1L) &
        is.finite(signature_only_score)
    ]

    n_patients <- nrow(analysis)
    n_events <- sum(analysis$OS.status == 1L, na.rm = TRUE)
    if (n_patients < min_patients || n_events < min_events) {
      skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
        cohort_label, cancer_type, sample_mode, default_cohort,
        reason = sprintf("below_threshold:n=%d_events=%d", n_patients, n_events),
        n_patients = n_patients,
        n_events = n_events
      )
      next
    }

    all_inputs[[length(all_inputs) + 1L]] <- analysis

    for (method in cutoff_methods) {
      grouped <- tcgasig_make_groups(analysis, method, minprop = optimized_minprop)
      if (!is.na(grouped$error)) {
        skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
          cohort_label, cancer_type, sample_mode, default_cohort,
          reason = paste0(method, "_cutoff_failed:", grouped$error),
          n_patients = n_patients,
          n_events = n_events
        )
        next
      }

      model_dt <- grouped$data
      if (length(unique(model_dt$group)) < 2) {
        skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
          cohort_label, cancer_type, sample_mode, default_cohort,
          reason = paste0(method, "_cutoff_single_group"),
          n_patients = n_patients,
          n_events = n_events
        )
        next
      }

      cox_row <- tryCatch(tcgasig_run_cox(model_dt, method, grouped$cutoff), error = function(e) e)
      if (inherits(cox_row, "error")) {
        skip_rows[[length(skip_rows) + 1L]] <- data.table::data.table(
          cohort_label, cancer_type, sample_mode, default_cohort,
          reason = paste0(method, "_cox_failed:", conditionMessage(cox_row)),
          n_patients = n_patients,
          n_events = n_events
        )
        next
      }

      result_row <- data.table::data.table(
        cancer_type = cancer_type,
        sample_mode = sample_mode,
        cohort_label = cohort_label,
        default_cohort = default_cohort,
        n_patients = n_patients,
        n_events = n_events,
        signature_name = signature_name,
        signature_score_mean = mean(analysis$signature_score, na.rm = TRUE),
        signature_genes_present = paste(usable_signature, collapse = ","),
        signature_genes_missing = ifelse(length(missing_signature) == 0, "none", paste(missing_signature, collapse = ",")),
        cutoff_method = method
      )
      result_row[, (paste0(signature_safe, "_score_mean")) := signature_score_mean]
      result_row[, (paste0(signature_safe, "_genes_present")) := signature_genes_present]
      result_row[, (paste0(signature_safe, "_genes_missing")) := signature_genes_missing]
      result_row <- cbind(result_row, cox_row[, setdiff(names(cox_row), "cutoff_method"), with = FALSE])
      result_row[, direction := data.table::fifelse(
        is.na(HR), NA_character_,
        data.table::fifelse(HR > 1, "high_worse", data.table::fifelse(HR < 1, "high_better", "neutral"))
      )]
      all_results[[length(all_results) + 1L]] <- result_row

      if (isTRUE(make_km)) {
        km_prefix <- if (cohort_label == cancer_type) cancer_type else cohort_label
        km_file <- file.path(
          km_dir,
          paste0(tcgasig_safe_name(km_prefix), "_", signature_safe, "_signature_only_KM_", method, ".pdf")
        )
        tcgasig_plot_signature_only_km(model_dt, km_prefix, sample_mode, method,
                                       grouped$cutoff, cox_row, km_file,
                                       signature_name)
      }
    }
  }

  survival_input <- if (length(all_inputs) > 0) data.table::rbindlist(all_inputs, fill = TRUE) else data.table::data.table()
  skips <- if (length(skip_rows) > 0) {
    data.table::rbindlist(skip_rows, fill = TRUE)
  } else {
    data.table::data.table(
      cohort_label = character(),
      cancer_type = character(),
      sample_mode = character(),
      default_cohort = logical(),
      reason = character(),
      n_patients = integer(),
      n_events = integer()
    )
  }
  if (length(all_results) == 0) {
    if (isTRUE(write_files)) data.table::fwrite(skips, file.path(out_dir, paste0(output_prefix, "_skipped_cohorts.tsv")), sep = "\t")
    tcgasig_stop("No cancer type passed analysis thresholds.")
  }

  results <- data.table::rbindlist(all_results, fill = TRUE)
  results[, FDR := stats::p.adjust(wald_p, method = "BH"), by = cutoff_method]
  desired <- c(
    "cancer_type", "sample_mode", "cohort_label", "default_cohort",
    "n_patients", "n_events", "signature_name",
    "signature_score_mean", paste0(signature_safe, "_score_mean"),
    "signature_genes_present", "signature_genes_missing",
    paste0(signature_safe, "_genes_present"), paste0(signature_safe, "_genes_missing"),
    "cutoff_method", "cutoff_value", "high_n", "low_n",
    "high_events", "low_events", "HR", "CI_lower", "CI_upper",
    "wald_p", "FDR", "direction"
  )
  data.table::setcolorder(results, c(intersect(desired, names(results)), setdiff(names(results), desired)))
  data.table::setorder(results, cutoff_method, FDR, cancer_type, sample_mode)

  if (isTRUE(write_files)) {
    data.table::fwrite(survival_input, file.path(out_dir, paste0(output_prefix, "_survival_input.tsv")), sep = "\t")
    data.table::fwrite(skips, file.path(out_dir, paste0(output_prefix, "_skipped_cohorts.tsv")), sep = "\t")
    data.table::fwrite(results, file.path(out_dir, paste0(output_prefix, "_survival_results.tsv")), sep = "\t")

    default_results <- results[default_cohort == TRUE]
    qc_lines <- c(
      paste0("pan-TCGA ", signature_name, " signature-only survival QC report"),
      paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "",
      paste0("Interpretation: bulk RNA-seq ", signature_name,
             " signature high/low survival; no target gene was added."),
      paste0("Signature name: ", signature_name),
      paste0("Signature genes: ", paste(signature_genes, collapse = ", ")),
      paste0("Minimum usable signature genes per cohort: ", min_signature_genes),
      paste0("Minimum patients per cohort: ", min_patients),
      paste0("Minimum death events per cohort: ", min_events),
      paste0("Expression scale: ", tcgasig_prepared_value(prepared, "expression_scale", "unknown")),
      paste0("Expression scale note: ", tcgasig_prepared_value(prepared, "expression_scale_note", "unknown")),
      "",
      paste0("Analyzed cohort-method rows: ", nrow(results)),
      paste0("Analyzed default cancer types: ", length(unique(default_results$cancer_type))),
      paste0("Skipped cohort rows: ", nrow(skips)),
      paste0("SKCM modes present in prepared data: ", paste(sort(unique(cohort_samples[cancer_type == "SKCM", sample_mode])), collapse = ", ")),
      paste0("SKCM primary-only included in survival run: ", include_skcm_primary_only),
      "",
      "Skipped cohorts:",
      if (nrow(skips) == 0) "none" else paste(
        skips$cohort_label,
        skips$sample_mode,
        skips$reason,
        paste0("n=", skips$n_patients),
        paste0("events=", skips$n_events),
        sep = "\t"
      )
    )
    writeLines(qc_lines, file.path(out_dir, paste0(output_prefix, "_QC_report.txt")))
  }

  tcgasig_message(verbose, "Wrote signature-only survival results for ", nrow(results), " cohort-method rows.")
  list(
    results = results,
    survival_input = survival_input,
    skipped_cohorts = skips,
    output_prefix = output_prefix,
    out_dir = out_dir,
    km_dir = km_dir
  )
}

plot_signature_only_forest <- function(
    results = NULL,
    results_file = NULL,
    signature_name = NULL,
    out_dir = file.path(getwd(), "results"),
    output_prefix = NULL,
    width = 10.5,
    height = NULL) {
  results <- tcgasig_results_table(results, results_file)
  if (is.null(signature_name)) signature_name <- if ("signature_name" %in% names(results)) unique(results$signature_name)[1] else "Signature"
  output_prefix <- tcgasig_signature_only_output_prefix(signature_name, output_prefix)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_dt <- tcgasig_prepare_plot_dt(results)
  forest <- ggplot2::ggplot(plot_dt, ggplot2::aes(y = cancer_plot)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    ggplot2::geom_segment(
      ggplot2::aes(x = CI_lower, xend = CI_upper, yend = cancer_plot, color = direction),
      linewidth = 0.6,
      alpha = 0.9
    ) +
    ggplot2::geom_point(ggplot2::aes(x = HR, color = direction), size = 2.4) +
    ggplot2::facet_wrap(~method_label, nrow = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(
      values = c(high_worse = "#B2182B", high_better = "#2166AC", neutral = "grey35"),
      na.value = "grey60"
    ) +
    ggplot2::labs(
      title = paste0("Pan-TCGA ", signature_name, " Signature-Only Survival"),
      subtitle = "Default cohorts use primary tumors except SKCM, where tumor_all is code 01 + 06",
      x = paste0("Hazard ratio, ", signature_name, " signature high vs low (log scale)"),
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (is.null(height)) height <- max(5, 0.22 * length(unique(plot_dt$cancer_type)) + 2.5)
  forest_pdf <- file.path(out_dir, paste0(output_prefix, "_forest_plot.pdf"))
  forest_png <- file.path(out_dir, paste0(output_prefix, "_forest_plot.png"))
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(forest_pdf, forest, width = width, height = height, device = pdf_device)
  ggplot2::ggsave(forest_png, forest, width = width, height = height, dpi = 300)
  list(plot = forest, pdf = forest_pdf, png = forest_png)
}

plot_signature_only_heatmap <- function(
    results = NULL,
    results_file = NULL,
    signature_name = NULL,
    out_dir = file.path(getwd(), "results"),
    output_prefix = NULL,
    width = 7.5,
    height = NULL) {
  results <- tcgasig_results_table(results, results_file)
  if (is.null(signature_name)) signature_name <- if ("signature_name" %in% names(results)) unique(results$signature_name)[1] else "Signature"
  output_prefix <- tcgasig_signature_only_output_prefix(signature_name, output_prefix)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  heat_dt <- tcgasig_prepare_plot_dt(results)
  heat_dt[, log2_HR := log2(HR)]
  heat_dt[, significance := data.table::fifelse(FDR < 0.05, "*", "")]
  heatmap <- ggplot2::ggplot(heat_dt, ggplot2::aes(x = method_label, y = cancer_plot, fill = log2_HR)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = significance), size = 4, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "log2(HR)"
    ) +
    ggplot2::labs(
      title = paste0("Pan-TCGA ", signature_name, " Signature-Only Cox Direction"),
      subtitle = "* marks FDR < 0.05; red means high group has higher OS hazard",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  if (is.null(height)) height <- max(5, 0.22 * length(unique(heat_dt$cancer_type)) + 2.2)
  heatmap_pdf <- file.path(out_dir, paste0(output_prefix, "_heatmap.pdf"))
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(heatmap_pdf, heatmap, width = width, height = height, device = pdf_device)
  list(plot = heatmap, pdf = heatmap_pdf)
}

plot_signature_only_results <- function(
    results = NULL,
    results_file = NULL,
    signature_name = NULL,
    out_dir = file.path(getwd(), "results"),
    output_prefix = NULL) {
  forest <- plot_signature_only_forest(
    results = results,
    results_file = results_file,
    signature_name = signature_name,
    out_dir = out_dir,
    output_prefix = output_prefix
  )
  heatmap <- plot_signature_only_heatmap(
    results = results,
    results_file = results_file,
    signature_name = signature_name,
    out_dir = out_dir,
    output_prefix = output_prefix
  )
  list(forest = forest, heatmap = heatmap)
}

run_pan_tcga_signature_only_survival <- function(
    signature_genes,
    signature_name = "Signature",
    project_dir = getwd(),
    data_dir = file.path(project_dir, "data", "xena"),
    out_dir = file.path(project_dir, "results"),
    output_prefix = NULL,
    min_signature_genes = 6,
    min_patients = 80,
    min_events = 20,
    include_skcm_primary_only = TRUE,
    cutoff_methods = c("optimized", "median"),
    optimized_minprop = 0.10,
    make_km = TRUE,
    make_summary_plots = TRUE,
    write_files = TRUE,
    verbose = TRUE) {
  output_prefix <- tcgasig_signature_only_output_prefix(signature_name, output_prefix)
  prepared <- prepare_signature_only_data(
    signature_genes = signature_genes,
    signature_name = signature_name,
    project_dir = project_dir,
    data_dir = data_dir,
    out_dir = out_dir,
    output_prefix = output_prefix,
    write_files = write_files,
    verbose = verbose
  )
  survival <- run_signature_only_survival(
    prepared = prepared,
    signature_genes = signature_genes,
    signature_name = signature_name,
    out_dir = out_dir,
    output_prefix = output_prefix,
    min_signature_genes = min_signature_genes,
    min_patients = min_patients,
    min_events = min_events,
    include_skcm_primary_only = include_skcm_primary_only,
    cutoff_methods = cutoff_methods,
    optimized_minprop = optimized_minprop,
    make_km = make_km,
    write_files = write_files,
    verbose = verbose
  )
  plots <- NULL
  if (isTRUE(make_summary_plots)) {
    plots <- plot_signature_only_results(
      results = survival$results,
      signature_name = signature_name,
      out_dir = out_dir,
      output_prefix = output_prefix
    )
    qc_file <- file.path(out_dir, paste0(output_prefix, "_QC_report.txt"))
    if (isTRUE(write_files) && file.exists(qc_file)) {
      cat(
        "\nPlot outputs:\n",
        paste0("- ", plots$forest$pdf, "\n"),
        paste0("- ", plots$forest$png, "\n"),
        paste0("- ", plots$heatmap$pdf, "\n"),
        file = qc_file,
        append = TRUE
      )
    }
  }
  list(prepared = prepared, survival = survival, plots = plots)
}
