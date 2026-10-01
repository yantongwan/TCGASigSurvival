default_treg_signature <- function() {
  c(
    "FOXP3", "IL2RA", "CTLA4", "CCR8", "TIGIT", "TNFRSF18", "TNFRSF4",
    "IKZF2", "LAYN", "ENTPD1", "LRRC32", "BATF", "IL2RB"
  )
}

tcgasig_safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

tcgasig_output_prefix <- function(target_gene, signature_name, output_prefix = NULL) {
  if (!is.null(output_prefix) && nzchar(output_prefix)) return(output_prefix)
  paste0("panTCGA_", tcgasig_safe_name(target_gene), "_", tcgasig_safe_name(signature_name))
}

tcgasig_message <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  }
}

tcgasig_stop <- function(...) {
  stop(paste0(...), call. = FALSE)
}

tcgasig_clean_missing <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "N/A", "na", "null", "NULL", "--", "'--",
             "not reported", "Not Reported", "unknown", "Unknown")] <- NA_character_
  x
}

tcgasig_to_numeric <- function(x) {
  suppressWarnings(as.numeric(tcgasig_clean_missing(x)))
}

tcgasig_first_non_missing <- function(x) {
  x <- tcgasig_clean_missing(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_character_ else x[1]
}

tcgasig_extract_attr <- function(attributes, key) {
  pattern <- paste0(key, " \"([^\"]+)\"")
  out <- stringr::str_match(attributes, pattern)[, 2]
  tcgasig_clean_missing(out)
}

tcgasig_sample_type_code <- function(barcode) {
  barcode <- tcgasig_clean_missing(barcode)
  out <- rep(NA_character_, length(barcode))
  ok <- !is.na(barcode) & nchar(barcode) >= 15
  out[ok] <- substr(barcode[ok], 14, 15)
  out
}

tcgasig_patient_from_barcode <- function(barcode) {
  barcode <- tcgasig_clean_missing(barcode)
  out <- rep(NA_character_, length(barcode))
  ok <- !is.na(barcode) & nchar(barcode) >= 12
  out[ok] <- substr(barcode[ok], 1, 12)
  out
}

tcgasig_read_tsv_auto <- function(path) {
  if (!file.exists(path)) tcgasig_stop("Missing input file: ", path)
  if (grepl("\\.gz$", path)) {
    tcgasig_read_delimited(path)
  } else {
    data.table::fread(path, sep = "\t", data.table = TRUE)
  }
}

tcgasig_parse_gencode_gene_map <- function(gtf_file, needed_symbols, verbose = TRUE) {
  tcgasig_message(verbose, "Parsing GENCODE annotation: ", gtf_file)
  cmd <- sprintf("gzip -cd %s | awk '$3==\"gene\"'", shQuote(gtf_file))
  gtf <- if (!tcgasig_has_unix_streaming()) {
    tcgasig_stream_rows(gtf_file, gtf = TRUE)
  } else tryCatch(
    data.table::fread(cmd = cmd, sep = "\t", header = FALSE, quote = "", fill = TRUE, data.table = TRUE),
    error = function(e) {
      tcgasig_message(verbose, "Falling back to full GTF read because shell filtering failed.")
      tcgasig_stream_rows(gtf_file, gtf = TRUE)
    }
  )
  if (ncol(gtf) < 9) tcgasig_stop("GENCODE GTF parsing failed; expected at least 9 columns: ", gtf_file)
  if ("V3" %in% names(gtf)) gtf <- gtf[gtf[["V3"]] == "gene"]
  gene_map <- data.table::data.table(
    gene_id = tcgasig_extract_attr(gtf[[9]], "gene_id"),
    gene_name = tcgasig_extract_attr(gtf[[9]], "gene_name"),
    gene_type = tcgasig_extract_attr(gtf[[9]], "gene_type")
  )
  gene_map <- unique(gene_map[!is.na(gene_id) & !is.na(gene_name)])
  gene_map <- gene_map[gene_name %in% needed_symbols]
  if (nrow(gene_map) == 0) tcgasig_stop("No requested marker symbols were found in GENCODE annotation.")
  gene_map[]
}

tcgasig_read_selected_expression <- function(expr_file, gene_ids, verbose = TRUE) {
  gene_ids <- unique(tcgasig_clean_missing(gene_ids))
  gene_ids <- gene_ids[!is.na(gene_ids)]
  if (length(gene_ids) == 0) tcgasig_stop("No Ensembl gene IDs are available for expression extraction.")

  ids_file <- tempfile("tcgasig_gene_ids_")
  writeLines(gene_ids, ids_file)
  on.exit(unlink(ids_file), add = TRUE)

  tcgasig_message(verbose, "Streaming selected target/signature rows from Xena expression matrix.")
  cmd <- sprintf(
    "gzip -cd %s | awk -v ids=%s 'BEGIN{while ((getline id < ids) > 0) keep[id]=1} NR==1 || ($1 in keep)'",
    shQuote(expr_file),
    shQuote(ids_file)
  )
  expr <- if (tcgasig_has_unix_streaming()) {
    data.table::fread(cmd = cmd, sep = "\t", data.table = TRUE)
  } else tcgasig_stream_rows(expr_file, keys = gene_ids)
  if (nrow(expr) == 0 || ncol(expr) < 2) {
    tcgasig_stop("No expression rows were extracted from: ", expr_file)
  }
  data.table::setnames(expr, 1, "gene_id")
  expr[]
}

tcgasig_detect_and_transform_expression <- function(expr_dt, sample_cols) {
  expr_dt[, (sample_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = sample_cols]
  probe_cols <- sample_cols[seq_len(min(length(sample_cols), 400L))]
  vals <- as.numeric(unlist(expr_dt[, ..probe_cols], use.names = FALSE))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) tcgasig_stop("Expression matrix has no finite values in selected genes.")

  min_val <- min(vals, na.rm = TRUE)
  max_val <- max(vals, na.rm = TRUE)
  q99 <- as.numeric(stats::quantile(vals, 0.99, na.rm = TRUE))

  if (min_val < 0 || max_val <= 35) {
    scale <- "already_log_transformed"
    note <- "Expression values were kept as Xena log-scale values; no second log transform was applied."
  } else if (q99 > 50 || max_val > 100) {
    scale <- "raw_tpm_converted_to_log2_tpm_plus_1"
    note <- "Expression values looked like raw TPM and were converted to log2(TPM + 1)."
    expr_dt[, (sample_cols) := lapply(.SD, function(x) log2(pmax(as.numeric(x), 0) + 1)), .SDcols = sample_cols]
  } else {
    scale <- "ambiguous_but_treated_as_log_transformed"
    note <- "Expression values were in a low numeric range and were treated as already log transformed."
  }

  list(expr = expr_dt, scale = scale, note = note, min = min_val, max = max_val, q99 = q99)
}

tcgasig_collapse_duplicate_symbols <- function(expr_dt, sample_cols) {
  expr_dt[, row_mean_expression := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = sample_cols]
  data.table::setorder(expr_dt, gene_name, -row_mean_expression, gene_id)
  collapsed <- expr_dt[, .SD[1], by = gene_name]
  collapsed[, row_mean_expression := NULL]
  data.table::setcolorder(collapsed, c("gene_name", "gene_id", setdiff(names(collapsed), c("gene_name", "gene_id"))))
  collapsed[]
}

tcgasig_normalize_survival <- function(survival_file) {
  survival <- tcgasig_read_tsv_auto(survival_file)
  required <- c("sample", "_PATIENT", "cancer type abbreviation", "OS", "OS.time")
  missing <- setdiff(required, names(survival))
  if (length(missing) > 0) {
    tcgasig_stop("Survival table lacks required column(s): ", paste(missing, collapse = ", "))
  }
  survival[, sample_barcode := tcgasig_clean_missing(sample)]
  survival[, patient_id := tcgasig_clean_missing(`_PATIENT`)]
  survival[is.na(patient_id), patient_id := tcgasig_patient_from_barcode(sample_barcode)]
  survival[, cancer_type := tcgasig_clean_missing(`cancer type abbreviation`)]
  survival[, OS.status := as.integer(tcgasig_to_numeric(OS))]
  survival[, OS.time := tcgasig_to_numeric(OS.time)]
  survival <- survival[!is.na(sample_barcode) & !is.na(patient_id) & !is.na(cancer_type)]
  survival[]
}

tcgasig_normalize_phenotype <- function(phenotype_file) {
  phenotype <- tcgasig_read_tsv_auto(phenotype_file)
  if (!"sample" %in% names(phenotype)) tcgasig_stop("Phenotype table lacks a sample column: ", phenotype_file)
  phenotype[, sample_barcode := tcgasig_clean_missing(sample)]
  phenotype[, phenotype_sample_type := if ("_sample_type" %in% names(phenotype)) tcgasig_clean_missing(`_sample_type`) else NA_character_]
  phenotype[, phenotype_study := if ("_study" %in% names(phenotype)) tcgasig_clean_missing(`_study`) else NA_character_]
  phenotype[, phenotype_primary_site := if ("_primary_site" %in% names(phenotype)) tcgasig_clean_missing(`_primary_site`) else NA_character_]
  phenotype[, phenotype_detailed_category := if ("detailed_category" %in% names(phenotype)) tcgasig_clean_missing(detailed_category) else NA_character_]
  phenotype[, .(
    sample_barcode,
    phenotype_sample_type,
    phenotype_study,
    phenotype_primary_site,
    phenotype_detailed_category
  )]
}

tcgasig_build_sample_metadata <- function(sample_cols, survival, phenotype) {
  samples <- data.table::data.table(sample_barcode = sample_cols)
  samples[, patient_id := tcgasig_patient_from_barcode(sample_barcode)]
  samples[, sample_type_code := tcgasig_sample_type_code(sample_barcode)]
  samples[, is_tcga_barcode := grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", sample_barcode)]

  surv_sample <- survival[, .(
    sample_barcode,
    survival_patient_id = patient_id,
    cancer_type,
    OS.status,
    OS.time
  )]
  samples <- merge(samples, surv_sample, by = "sample_barcode", all.x = TRUE, sort = FALSE)
  samples[is.na(patient_id), patient_id := survival_patient_id]
  samples[, survival_patient_id := NULL]

  patient_surv <- survival[
    !is.na(patient_id),
    .(
      cancer_type_patient = tcgasig_first_non_missing(cancer_type),
      OS.status_patient = tcgasig_to_numeric(tcgasig_first_non_missing(OS.status)),
      OS.time_patient = tcgasig_to_numeric(tcgasig_first_non_missing(OS.time))
    ),
    by = patient_id
  ]
  samples <- merge(samples, patient_surv, by = "patient_id", all.x = TRUE, sort = FALSE)
  samples[is.na(cancer_type), cancer_type := cancer_type_patient]
  samples[is.na(OS.status), OS.status := as.integer(OS.status_patient)]
  samples[is.na(OS.time), OS.time := OS.time_patient]
  samples[, c("cancer_type_patient", "OS.status_patient", "OS.time_patient") := NULL]

  samples <- merge(samples, phenotype, by = "sample_barcode", all.x = TRUE, sort = FALSE)
  samples[, sample_type := phenotype_sample_type]
  samples[is.na(sample_type) & sample_type_code == "01", sample_type := "Primary Tumor"]
  samples[is.na(sample_type) & sample_type_code == "06", sample_type := "Metastatic"]
  samples[, has_os := !is.na(OS.status) & !is.na(OS.time) & OS.time > 0]
  samples[]
}

tcgasig_build_cohort_samples <- function(sample_metadata) {
  meta <- sample_metadata[
    is_tcga_barcode == TRUE &
      !is.na(cancer_type) &
      has_os == TRUE &
      sample_type_code %in% c("01", "06")
  ]

  non_skcm <- meta[cancer_type != "SKCM" & sample_type_code == "01"]
  non_skcm[, sample_mode := "primary_only"]
  non_skcm[, default_cohort := TRUE]

  skcm_all <- meta[cancer_type == "SKCM" & sample_type_code %in% c("01", "06")]
  skcm_all[, sample_mode := "tumor_all"]
  skcm_all[, default_cohort := TRUE]

  skcm_primary <- meta[cancer_type == "SKCM" & sample_type_code == "01"]
  skcm_primary[, sample_mode := "primary_only"]
  skcm_primary[, default_cohort := FALSE]

  cohorts <- data.table::rbindlist(list(non_skcm, skcm_all, skcm_primary), fill = TRUE)
  cohorts[, cohort_label := data.table::fifelse(cancer_type == "SKCM" & sample_mode == "primary_only", "SKCM_primary_only", cancer_type)]
  cohorts[, km_label := cohort_label]
  data.table::setorder(cohorts, cohort_label, patient_id, sample_barcode)
  cohorts <- cohorts[, .SD[1], by = .(cohort_label, patient_id)]
  cohorts[]
}

prepare_signature_data <- function(
    target_gene,
    signature_genes = default_treg_signature(),
    signature_name = "Treg",
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
  if (missing(target_gene) || !nzchar(target_gene)) tcgasig_stop("target_gene must be provided.")
  signature_genes <- unique(tcgasig_clean_missing(signature_genes))
  signature_genes <- signature_genes[!is.na(signature_genes)]
  if (length(signature_genes) == 0) tcgasig_stop("signature_genes must contain at least one gene symbol.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  output_prefix <- tcgasig_output_prefix(target_gene, signature_name, output_prefix)

  for (path in c(expression_file, survival_file, phenotype_file, gencode_file)) {
    if (!file.exists(path)) {
      tcgasig_stop("Missing input file: ", path, ". Download UCSC Xena files first.")
    }
  }

  tcgasig_message(verbose, "Preparing pan-TCGA ", target_gene, "-", signature_name, " data.")
  tcgasig_message(verbose, "Interpretation: bulk RNA-seq ", target_gene, "-", signature_name,
                  " signature high/low; not direct cell-specific ", target_gene, " expression.")

  needed_symbols <- unique(c(target_gene, signature_genes))
  gene_map <- tcgasig_parse_gencode_gene_map(gencode_file, needed_symbols, verbose = verbose)

  missing_from_gencode <- setdiff(needed_symbols, gene_map$gene_name)
  if (target_gene %in% missing_from_gencode) {
    tcgasig_stop("Target gene ", target_gene, " was not found in GENCODE annotation.")
  }

  selected_expr <- tcgasig_read_selected_expression(expression_file, gene_map$gene_id, verbose = verbose)
  selected_expr <- merge(gene_map[, .(gene_id, gene_name, gene_type)], selected_expr,
                         by = "gene_id", all.y = TRUE, sort = FALSE)
  sample_cols <- setdiff(names(selected_expr), c("gene_id", "gene_name", "gene_type"))
  if (length(sample_cols) == 0) tcgasig_stop("Selected expression matrix contains no sample columns.")

  transformed <- tcgasig_detect_and_transform_expression(selected_expr, sample_cols)
  collapsed_expr <- tcgasig_collapse_duplicate_symbols(transformed$expr, sample_cols)

  present_symbols <- collapsed_expr$gene_name
  missing_symbols <- setdiff(needed_symbols, present_symbols)
  if (target_gene %in% missing_symbols) {
    tcgasig_stop("Target gene ", target_gene, " was not found after expression extraction.")
  }

  survival <- tcgasig_normalize_survival(survival_file)
  phenotype <- tcgasig_normalize_phenotype(phenotype_file)
  sample_metadata <- tcgasig_build_sample_metadata(sample_cols, survival, phenotype)
  cohort_samples <- tcgasig_build_cohort_samples(sample_metadata)

  prepared <- list(
    target_gene = target_gene,
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
        "target_gene",
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
        target_gene,
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

  tcgasig_message(verbose, "Prepared expression rows: ", nrow(collapsed_expr),
                  "; cohort sample rows: ", nrow(cohort_samples))
  tcgasig_message(verbose, "Expression scale: ", transformed$scale, " | ", transformed$note)
  prepared
}

tcgasig_zscore_or_na <- function(x) {
  x <- as.numeric(x)
  sx <- stats::sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

tcgasig_extract_gene_vector <- function(expr, gene, samples) {
  row <- expr[gene_name == gene]
  if (nrow(row) == 0) return(rep(NA_real_, length(samples)))
  vals <- as.numeric(unlist(row[1, ..samples], use.names = FALSE))
  names(vals) <- samples
  vals
}

tcgasig_format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

tcgasig_open_pdf <- function(path, width, height) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(path, width = width, height = height)
  } else {
    grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  }
}

tcgasig_run_cox <- function(dt, cutoff_method, cutoff_value) {
  model_dt <- data.table::copy(dt)
  model_dt[, group := factor(group, levels = c("low", "high"))]
  if (length(unique(model_dt$group)) < 2) tcgasig_stop("Grouping produced fewer than two groups.")
  checked <- tcgasig_fit_checked(model_dt, "group", "grouphigh")
  if (!is.null(checked$error)) tcgasig_stop(checked$error)
  group_counts <- model_dt[, .(
    group_n = .N,
    group_events = sum(OS.status == 1L, na.rm = TRUE)
  ), by = group]
  data.table::data.table(
    cutoff_method = cutoff_method,
    cutoff_value = cutoff_value,
    high_n = group_counts[group == "high", group_n],
    low_n = group_counts[group == "low", group_n],
    high_events = group_counts[group == "high", group_events],
    low_events = group_counts[group == "low", group_events],
    HR = checked$row$HR,
    CI_lower = checked$row$CI_lower,
    CI_upper = checked$row$CI_upper,
    wald_p = checked$row$wald_p
  )
}

tcgasig_make_groups <- function(analysis_dt, method, minprop = 0.10) {
  dt <- data.table::copy(analysis_dt)
  if (method == "median") {
    cutoff <- stats::median(dt$target_signature_score, na.rm = TRUE)
    dt[, group := ifelse(target_signature_score >= cutoff, "high", "low")]
    return(list(data = dt, cutoff = cutoff, error = NA_character_))
  }

  if (method == "optimized") {
    cut_obj <- tryCatch(
      survminer::surv_cutpoint(
        dt,
        time = "OS.time",
        event = "OS.status",
        variables = "target_signature_score",
        minprop = minprop
      ),
      error = function(e) e
    )
    if (inherits(cut_obj, "error")) return(list(data = NULL, cutoff = NA_real_, error = conditionMessage(cut_obj)))
    cutoff <- as.numeric(cut_obj$cutpoint["target_signature_score", "cutpoint"])
    dt[, group := ifelse(target_signature_score > cutoff, "high", "low")]
    return(list(data = dt, cutoff = cutoff, error = NA_character_))
  }

  tcgasig_stop("Unknown cutoff method: ", method)
}

tcgasig_plot_km <- function(dt, cancer_label, sample_mode, cutoff_method, cutoff_value,
                            cox_row, out_file, target_gene, signature_name) {
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
    "TCGA %s %s-%s signature high/low (%s, %s)",
    cancer_label,
    target_gene,
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

tcgasig_prepared_value <- function(prepared, name, default = NULL) {
  if (!is.null(prepared[[name]])) prepared[[name]] else default
}

run_signature_survival <- function(
    prepared,
    target_gene = NULL,
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

  target_gene <- if (is.null(target_gene)) tcgasig_prepared_value(prepared, "target_gene") else target_gene
  signature_name <- if (is.null(signature_name)) tcgasig_prepared_value(prepared, "signature_name", "Signature") else signature_name
  signature_genes <- if (is.null(signature_genes)) {
    tcgasig_prepared_value(prepared, "signature_genes", tcgasig_prepared_value(prepared, "treg_signature", default_treg_signature()))
  } else {
    signature_genes
  }
  signature_genes <- unique(tcgasig_clean_missing(signature_genes))
  signature_genes <- signature_genes[!is.na(signature_genes)]
  output_prefix <- tcgasig_output_prefix(target_gene, signature_name, output_prefix)

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

  tcgasig_message(verbose, "Running pan-TCGA ", target_gene, "-", signature_name, " survival analysis.")
  tcgasig_message(verbose, "Minimum patients: ", min_patients, "; minimum OS death events: ", min_events)
  tcgasig_message(verbose, "Interpretation: bulk RNA-seq ", target_gene, "-", signature_name,
                  " signature high/low; not direct cell-specific ", target_gene, " expression.")

  cohort_defs <- unique(cohort_samples[, .(cohort_label, cancer_type, sample_mode, default_cohort)])
  data.table::setorder(cohort_defs, cancer_type, sample_mode)

  all_results <- list()
  all_inputs <- list()
  skip_rows <- list()
  target_safe <- tcgasig_safe_name(target_gene)
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

    target_expr <- tcgasig_extract_gene_vector(expr, target_gene, sample_cols)
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
      target_gene = target_gene,
      signature_name = signature_name,
      target_gene_expr = as.numeric(target_expr[cs$sample_barcode]),
      signature_score = as.numeric(signature_score),
      signature_genes_present = paste(usable_signature, collapse = ","),
      signature_genes_missing = ifelse(length(missing_signature) == 0, "none", paste(missing_signature, collapse = ","))
    )
    analysis[, target_z := tcgasig_zscore_or_na(target_gene_expr)]
    analysis[, signature_z := tcgasig_zscore_or_na(signature_score)]
    analysis[, target_signature_score := target_z + signature_z]
    analysis <- analysis[
      is.finite(OS.time) &
        OS.time > 0 &
        OS.status %in% c(0L, 1L) &
        is.finite(target_signature_score)
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
        target_gene = target_gene,
        signature_name = signature_name,
        target_gene_mean = mean(analysis$target_gene_expr, na.rm = TRUE),
        signature_genes_present = paste(usable_signature, collapse = ","),
        signature_genes_missing = ifelse(length(missing_signature) == 0, "none", paste(missing_signature, collapse = ",")),
        cutoff_method = method
      )
      result_row[, (paste0(target_safe, "_mean")) := target_gene_mean]
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
          paste0(tcgasig_safe_name(km_prefix), "_", target_safe, "_", signature_safe, "_KM_", method, ".pdf")
        )
        tcgasig_plot_km(model_dt, km_prefix, sample_mode, method, grouped$cutoff,
                        cox_row, km_file, target_gene, signature_name)
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
    "n_patients", "n_events", "target_gene", "signature_name",
    "target_gene_mean", paste0(target_safe, "_mean"),
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
      paste0("pan-TCGA ", target_gene, "-", signature_name, " survival QC report"),
      paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "",
      paste0("Interpretation: bulk RNA-seq ", target_gene, "-", signature_name,
             " signature high/low; not direct cell-specific ", target_gene, " expression."),
      paste0("Target gene: ", target_gene),
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

  tcgasig_message(verbose, "Wrote survival results for ", nrow(results), " cohort-method rows.")
  list(
    results = results,
    survival_input = survival_input,
    skipped_cohorts = skips,
    output_prefix = output_prefix,
    out_dir = out_dir,
    km_dir = km_dir
  )
}

tcgasig_format_fdr_label <- function(fdr) {
  ifelse(
    is.na(fdr), "",
    ifelse(fdr < 0.001, "FDR<0.001", paste0("FDR=", formatC(fdr, format = "f", digits = 3)))
  )
}

tcgasig_results_table <- function(results = NULL, results_file = NULL) {
  if (is.list(results) && !data.table::is.data.table(results) && !is.data.frame(results) && !is.null(results$results)) {
    results <- results$results
  }
  if (is.null(results)) {
    if (is.null(results_file)) tcgasig_stop("Provide either results or results_file.")
    results <- data.table::fread(results_file, sep = "\t", data.table = TRUE)
  }
  data.table::as.data.table(results)
}

tcgasig_prepare_plot_dt <- function(results) {
  needed <- c("cancer_type", "sample_mode", "default_cohort", "cutoff_method", "HR", "CI_lower", "CI_upper", "FDR", "wald_p", "direction")
  missing <- setdiff(needed, names(results))
  if (length(missing) > 0) tcgasig_stop("Results table lacks required column(s): ", paste(missing, collapse = ", "))
  plot_dt <- results[default_cohort == TRUE & is.finite(HR) & is.finite(CI_lower) & is.finite(CI_upper)]
  if (nrow(plot_dt) == 0) tcgasig_stop("No default cohort rows with finite HR/CI are available for plotting.")

  order_dt <- plot_dt[
    cutoff_method == "optimized",
    .(order_hr = HR[which.min(FDR)]),
    by = cancer_type
  ]
  if (nrow(order_dt) == 0) order_dt <- plot_dt[, .(order_hr = HR[1]), by = cancer_type]
  order_dt[, order_value := log2(order_hr)]
  data.table::setorder(order_dt, order_value)
  plot_dt[, cancer_plot := factor(cancer_type, levels = order_dt$cancer_type)]
  plot_dt[, method_label := factor(cutoff_method, levels = c("optimized", "median"), labels = c("Optimized cutoff", "Median cutoff"))]
  plot_dt[, fdr_label := tcgasig_format_fdr_label(FDR)]
  plot_dt[]
}

plot_signature_forest <- function(
    results = NULL,
    results_file = NULL,
    target_gene = NULL,
    signature_name = NULL,
    out_dir = file.path(getwd(), "results"),
    output_prefix = NULL,
    width = 10.5,
    height = NULL) {
  results <- tcgasig_results_table(results, results_file)
  if (is.null(target_gene)) target_gene <- unique(results$target_gene)[1]
  if (is.null(signature_name)) signature_name <- if ("signature_name" %in% names(results)) unique(results$signature_name)[1] else "Signature"
  output_prefix <- tcgasig_output_prefix(target_gene, signature_name, output_prefix)
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
      title = paste0("Pan-TCGA ", target_gene, "-", signature_name, " Signature Survival"),
      subtitle = paste0("Default cohorts use primary tumors except SKCM, where tumor_all is code 01 + 06; target gene: ", target_gene),
      x = paste0("Hazard ratio, ", target_gene, "-", signature_name, " signature high vs low (log scale)"),
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

plot_signature_heatmap <- function(
    results = NULL,
    results_file = NULL,
    target_gene = NULL,
    signature_name = NULL,
    out_dir = file.path(getwd(), "results"),
    output_prefix = NULL,
    width = 7.5,
    height = NULL) {
  results <- tcgasig_results_table(results, results_file)
  if (is.null(target_gene)) target_gene <- unique(results$target_gene)[1]
  if (is.null(signature_name)) signature_name <- if ("signature_name" %in% names(results)) unique(results$signature_name)[1] else "Signature"
  output_prefix <- tcgasig_output_prefix(target_gene, signature_name, output_prefix)
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
      title = paste0("Pan-TCGA ", target_gene, "-", signature_name, " Signature Cox Direction"),
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

plot_signature_results <- function(
    results = NULL,
    results_file = NULL,
    target_gene = NULL,
    signature_name = NULL,
    out_dir = file.path(getwd(), "results"),
    output_prefix = NULL) {
  forest <- plot_signature_forest(
    results = results,
    results_file = results_file,
    target_gene = target_gene,
    signature_name = signature_name,
    out_dir = out_dir,
    output_prefix = output_prefix
  )
  heatmap <- plot_signature_heatmap(
    results = results,
    results_file = results_file,
    target_gene = target_gene,
    signature_name = signature_name,
    out_dir = out_dir,
    output_prefix = output_prefix
  )
  list(forest = forest, heatmap = heatmap)
}

tcgasig_survival_result_tables <- function(survival = NULL, results = NULL,
                                           survival_input = NULL,
                                           results_file = NULL,
                                           survival_input_file = NULL) {
  if (is.list(survival) && !is.null(survival[["survival"]])) survival <- survival[["survival"]]
  if (is.list(survival) && !is.null(survival[["results"]])) {
    if (is.null(results)) results <- survival[["results"]]
    if (is.null(survival_input)) survival_input <- survival[["survival_input"]]
  }
  if (is.null(results)) {
    if (is.null(results_file)) tcgasig_stop("Provide results, survival, or results_file.")
    results <- data.table::fread(results_file, sep = "\t", data.table = TRUE)
  }
  if (is.null(survival_input)) {
    if (is.null(survival_input_file)) tcgasig_stop("Provide survival_input, survival, or survival_input_file.")
    survival_input <- data.table::fread(survival_input_file, sep = "\t", data.table = TRUE)
  }
  list(
    results = data.table::as.data.table(results),
    survival_input = data.table::as.data.table(survival_input)
  )
}

tcgasig_score_column <- function(survival_input) {
  if ("target_signature_score" %in% names(survival_input)) return("target_signature_score")
  if ("target_Treg_score" %in% names(survival_input)) return("target_Treg_score")
  legacy_candidates <- grep("_Treg_score$", names(survival_input), value = TRUE)
  if (length(legacy_candidates) == 1L) return(legacy_candidates)
  tcgasig_stop("survival_input lacks target_signature_score or target_Treg_score.")
}

tcgasig_expression_summary <- function(x, value_col = "target_gene_expr") {
  data.table::data.table(
    n = sum(is.finite(x[[value_col]])),
    mean = mean(x[[value_col]], na.rm = TRUE),
    median = stats::median(x[[value_col]], na.rm = TRUE),
    q25 = as.numeric(stats::quantile(x[[value_col]], 0.25, na.rm = TRUE)),
    q75 = as.numeric(stats::quantile(x[[value_col]], 0.75, na.rm = TRUE))
  )
}

tcgasig_wilcox_p <- function(dt) {
  if (length(unique(dt$group)) < 2) return(NA_real_)
  out <- tryCatch(
    stats::wilcox.test(target_gene_expr ~ group, data = dt)$p.value,
    error = function(e) NA_real_
  )
  as.numeric(out)
}

plot_target_expression_by_signature <- function(
    survival = NULL,
    results = NULL,
    survival_input = NULL,
    results_file = NULL,
    survival_input_file = NULL,
    target_gene = NULL,
    signature_name = NULL,
    cutoff_method = "optimized",
    default_only = TRUE,
    out_dir = file.path(getwd(), "results"),
    output_prefix = NULL,
    write_files = TRUE,
    width = 9,
    delta_height = NULL,
    boxplot_width = 11,
    boxplot_height = NULL) {
  tables <- tcgasig_survival_result_tables(
    survival = survival,
    results = results,
    survival_input = survival_input,
    results_file = results_file,
    survival_input_file = survival_input_file
  )
  results <- tables$results
  survival_input <- tables$survival_input
  if (is.null(target_gene)) {
    target_gene <- if ("target_gene" %in% names(results)) unique(results$target_gene)[1] else unique(survival_input$target_gene)[1]
  }
  if (is.null(signature_name)) {
    signature_name <- if ("signature_name" %in% names(results)) unique(results$signature_name)[1] else "Signature"
  }
  if (!"target_gene_expr" %in% names(survival_input)) {
    legacy_expr_col <- paste0(target_gene, "_expr")
    if (legacy_expr_col %in% names(survival_input)) {
      survival_input[, target_gene_expr := get(legacy_expr_col)]
    } else {
      tcgasig_stop("survival_input lacks target_gene_expr.")
    }
  }
  score_col <- tcgasig_score_column(survival_input)
  output_prefix <- tcgasig_output_prefix(target_gene, signature_name, output_prefix)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  needed_results <- c("cohort_label", "cancer_type", "sample_mode", "default_cohort", "cutoff_method", "cutoff_value")
  missing_results <- setdiff(needed_results, names(results))
  if (length(missing_results) > 0) {
    tcgasig_stop("results lacks required column(s): ", paste(missing_results, collapse = ", "))
  }
  needed_input <- c("cohort_label", "cancer_type", "sample_mode", "default_cohort", "patient_id", "sample_barcode", "target_gene_expr", score_col)
  missing_input <- setdiff(needed_input, names(survival_input))
  if (length(missing_input) > 0) {
    tcgasig_stop("survival_input lacks required column(s): ", paste(missing_input, collapse = ", "))
  }

  cutoff_method_value <- cutoff_method
  cutoff_dt <- results[cutoff_method == cutoff_method_value]
  if (isTRUE(default_only)) cutoff_dt <- cutoff_dt[default_cohort == TRUE]
  cutoff_dt <- unique(cutoff_dt[, .(cohort_label, cancer_type, sample_mode, default_cohort, cutoff_method, cutoff_value)])
  if (nrow(cutoff_dt) == 0) tcgasig_stop("No cutoff rows found for method: ", cutoff_method)

  plot_dt <- merge(
    survival_input,
    cutoff_dt,
    by = c("cohort_label", "cancer_type", "sample_mode", "default_cohort"),
    all = FALSE,
    sort = FALSE
  )
  plot_dt <- plot_dt[is.finite(target_gene_expr) & is.finite(get(score_col))]
  if (nrow(plot_dt) == 0) tcgasig_stop("No samples remain after matching survival input to cutoff rows.")
  if (cutoff_method == "median") {
    plot_dt[, group := ifelse(get(score_col) >= cutoff_value, "high", "low")]
  } else {
    plot_dt[, group := ifelse(get(score_col) > cutoff_value, "high", "low")]
  }
  plot_dt[, group := factor(group, levels = c("low", "high"))]

  baseline <- plot_dt[, tcgasig_expression_summary(.SD), by = .(cohort_label, cancer_type, sample_mode, default_cohort)]
  data.table::setnames(
    baseline,
    c("n", "mean", "median", "q25", "q75"),
    c("overall_n", "overall_mean", "overall_median", "overall_q25", "overall_q75")
  )
  group_summary <- plot_dt[, tcgasig_expression_summary(.SD), by = .(cohort_label, cancer_type, sample_mode, default_cohort, group)]
  data.table::setnames(
    group_summary,
    c("n", "mean", "median", "q25", "q75"),
    c("group_n", "group_mean", "group_median", "group_q25", "group_q75")
  )
  summary_dt <- merge(group_summary, baseline, by = c("cohort_label", "cancer_type", "sample_mode", "default_cohort"), all.x = TRUE, sort = FALSE)
  p_dt <- plot_dt[, .(wilcox_p = tcgasig_wilcox_p(.SD)), by = .(cohort_label, cancer_type, sample_mode, default_cohort)]
  summary_dt <- merge(summary_dt, p_dt, by = c("cohort_label", "cancer_type", "sample_mode", "default_cohort"), all.x = TRUE, sort = FALSE)
  summary_dt[, target_gene := target_gene]
  summary_dt[, signature_name := signature_name]
  summary_dt[, cutoff_method := cutoff_method]
  summary_dt[, delta_mean_vs_overall := group_mean - overall_mean]
  summary_dt[, delta_median_vs_overall := group_median - overall_median]
  data.table::setcolorder(summary_dt, c(
    "target_gene", "signature_name", "cancer_type", "sample_mode", "cohort_label", "default_cohort",
    "cutoff_method", "group", "group_n", "group_mean", "group_median", "group_q25", "group_q75",
    "overall_n", "overall_mean", "overall_median", "overall_q25", "overall_q75",
    "delta_mean_vs_overall", "delta_median_vs_overall", "wilcox_p"
  ))

  order_dt <- summary_dt[group == "high", .(order_value = delta_mean_vs_overall[1]), by = cancer_type]
  data.table::setorder(order_dt, order_value)
  summary_dt[, cancer_plot := factor(cancer_type, levels = order_dt$cancer_type)]
  plot_dt[, cancer_plot := factor(cancer_type, levels = order_dt$cancer_type)]

  delta_plot <- ggplot2::ggplot(summary_dt, ggplot2::aes(x = delta_mean_vs_overall, y = cancer_plot, color = group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = group_n), alpha = 0.9, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_color_manual(values = c(low = "#2C7FB8", high = "#D95F0E")) +
    ggplot2::scale_size_continuous(range = c(1.8, 4.2)) +
    ggplot2::labs(
      title = paste0("Bulk ", target_gene, " Expression by ", signature_name, " Signature Group"),
      subtitle = paste0("Groups use ", cutoff_method, " target-signature cutpoints; x-axis is group mean minus all-sample mean within each cancer"),
      x = paste0("Mean ", target_gene, " expression minus overall mean"),
      y = NULL,
      color = NULL,
      size = "n"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "bottom"
    )

  baseline_for_plot <- unique(plot_dt[, .(cohort_label, overall_median = stats::median(target_gene_expr, na.rm = TRUE))])
  box_dt <- merge(plot_dt, baseline_for_plot, by = "cohort_label", all.x = TRUE, sort = FALSE)
  box_plot <- ggplot2::ggplot(box_dt, ggplot2::aes(x = group, y = target_gene_expr, fill = group)) +
    ggplot2::geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.75) +
    ggplot2::geom_jitter(width = 0.12, size = 0.35, alpha = 0.35, color = "grey25") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = overall_median), linetype = "dashed", color = "grey35", linewidth = 0.35) +
    ggplot2::facet_wrap(~cancer_plot, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(low = "#2C7FB8", high = "#D95F0E")) +
    ggplot2::labs(
      title = paste0("Bulk ", target_gene, " Expression in ", signature_name, " Signature-High and Signature-Low Tumors"),
      subtitle = "Dashed line is the overall median target-gene expression within each cancer",
      x = NULL,
      y = paste0(target_gene, " expression"),
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 8),
      legend.position = "bottom"
    )

  if (is.null(delta_height)) delta_height <- max(5, 0.22 * length(unique(summary_dt$cancer_type)) + 2.5)
  if (is.null(boxplot_height)) boxplot_height <- max(7, 0.55 * ceiling(length(unique(box_dt$cancer_type)) / 4) + 5)
  delta_pdf <- file.path(out_dir, paste0(output_prefix, "_target_expression_delta_", cutoff_method, ".pdf"))
  delta_png <- file.path(out_dir, paste0(output_prefix, "_target_expression_delta_", cutoff_method, ".png"))
  boxplot_pdf <- file.path(out_dir, paste0(output_prefix, "_target_expression_boxplot_", cutoff_method, ".pdf"))
  summary_tsv <- file.path(out_dir, paste0(output_prefix, "_target_expression_summary_", cutoff_method, ".tsv"))

  if (isTRUE(write_files)) {
    summary_write <- data.table::copy(summary_dt)
    summary_write[, cancer_plot := NULL]
    data.table::fwrite(summary_write, summary_tsv, sep = "\t")
    pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
    ggplot2::ggsave(delta_pdf, delta_plot, width = width, height = delta_height, device = pdf_device)
    ggplot2::ggsave(delta_png, delta_plot, width = width, height = delta_height, dpi = 300)
    ggplot2::ggsave(boxplot_pdf, box_plot, width = boxplot_width, height = boxplot_height, device = pdf_device)
  }

  list(
    summary = summary_dt,
    plot_data = plot_dt,
    delta_plot = delta_plot,
    box_plot = box_plot,
    summary_tsv = summary_tsv,
    delta_pdf = delta_pdf,
    delta_png = delta_png,
    boxplot_pdf = boxplot_pdf
  )
}

run_multi_signature_target_expression <- function(
    target_gene,
    signature_sets,
    project_dir = getwd(),
    data_dir = file.path(project_dir, "data", "xena"),
    out_dir = file.path(project_dir, "results"),
    output_prefix = NULL,
    min_signature_genes = 6,
    min_patients = 80,
    min_events = 20,
    include_skcm_primary_only = TRUE,
    cutoff_method = "optimized",
    optimized_minprop = 0.10,
    make_km = FALSE,
    write_files = TRUE,
    verbose = TRUE) {
  if (!is.list(signature_sets) || length(signature_sets) == 0) {
    tcgasig_stop("signature_sets must be a named list of gene-symbol vectors.")
  }
  if (is.null(names(signature_sets)) || any(!nzchar(names(signature_sets)))) {
    names(signature_sets) <- paste0("Signature", seq_along(signature_sets))
  }
  all_signature_genes <- unique(unlist(signature_sets, use.names = FALSE))
  multi_prefix <- if (is.null(output_prefix)) {
    paste0("panTCGA_", tcgasig_safe_name(target_gene), "_multi_signature")
  } else {
    output_prefix
  }
  prepared <- prepare_signature_data(
    target_gene = target_gene,
    signature_genes = all_signature_genes,
    signature_name = "MultiSignature",
    project_dir = project_dir,
    data_dir = data_dir,
    out_dir = out_dir,
    output_prefix = multi_prefix,
    write_files = write_files,
    verbose = verbose
  )

  survival_runs <- list()
  expression_runs <- list()
  for (signature_name in names(signature_sets)) {
    prefix <- paste0("panTCGA_", tcgasig_safe_name(target_gene), "_", tcgasig_safe_name(signature_name))
    survival_runs[[signature_name]] <- run_signature_survival(
      prepared = prepared,
      target_gene = target_gene,
      signature_genes = signature_sets[[signature_name]],
      signature_name = signature_name,
      out_dir = out_dir,
      output_prefix = prefix,
      min_signature_genes = min_signature_genes,
      min_patients = min_patients,
      min_events = min_events,
      include_skcm_primary_only = include_skcm_primary_only,
      cutoff_methods = cutoff_method,
      optimized_minprop = optimized_minprop,
      make_km = make_km,
      write_files = write_files,
      verbose = verbose
    )
    expression_runs[[signature_name]] <- plot_target_expression_by_signature(
      survival = survival_runs[[signature_name]],
      target_gene = target_gene,
      signature_name = signature_name,
      cutoff_method = cutoff_method,
      out_dir = out_dir,
      output_prefix = prefix,
      write_files = write_files
    )
  }

  combined_summary <- data.table::rbindlist(lapply(expression_runs, `[[`, "summary"), fill = TRUE)
  combined_summary[, signature_name := factor(as.character(signature_name), levels = names(signature_sets))]
  combined_plot <- ggplot2::ggplot(combined_summary, ggplot2::aes(x = delta_mean_vs_overall, y = cancer_plot, color = group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(size = group_n), alpha = 0.9, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::facet_wrap(~signature_name, nrow = 1) +
    ggplot2::scale_color_manual(values = c(low = "#2C7FB8", high = "#D95F0E")) +
    ggplot2::scale_size_continuous(range = c(1.5, 3.8)) +
    ggplot2::labs(
      title = paste0("Bulk ", target_gene, " Expression Across Cell Signature Groups"),
      subtitle = paste0("Groups use ", cutoff_method, " target-signature cutpoints; baseline is all samples within each cancer"),
      x = paste0("Mean ", target_gene, " expression minus overall mean"),
      y = NULL,
      color = NULL,
      size = "n"
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
  combined_tsv <- file.path(out_dir, paste0(multi_prefix, "_target_expression_summary_", cutoff_method, ".tsv"))
  combined_pdf <- file.path(out_dir, paste0(multi_prefix, "_target_expression_delta_", cutoff_method, ".pdf"))
  combined_png <- file.path(out_dir, paste0(multi_prefix, "_target_expression_delta_", cutoff_method, ".png"))
  if (isTRUE(write_files)) {
    combined_write <- data.table::copy(combined_summary)
    combined_write[, cancer_plot := NULL]
    data.table::fwrite(combined_write, combined_tsv, sep = "\t")
    pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
    plot_width <- max(9, 4.5 * length(signature_sets))
    plot_height <- max(5, 0.22 * length(unique(combined_summary$cancer_type)) + 2.5)
    ggplot2::ggsave(combined_pdf, combined_plot, width = plot_width, height = plot_height, device = pdf_device)
    ggplot2::ggsave(combined_png, combined_plot, width = plot_width, height = plot_height, dpi = 300)
  }
  list(
    prepared = prepared,
    survival_runs = survival_runs,
    expression_runs = expression_runs,
    combined_summary = combined_summary,
    combined_plot = combined_plot,
    combined_tsv = combined_tsv,
    combined_pdf = combined_pdf,
    combined_png = combined_png
  )
}

run_pan_tcga_signature_survival <- function(
    target_gene,
    signature_genes = default_treg_signature(),
    signature_name = "Treg",
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
    make_expression_plot = TRUE,
    write_files = TRUE,
    verbose = TRUE) {
  output_prefix <- tcgasig_output_prefix(target_gene, signature_name, output_prefix)
  prepared <- prepare_signature_data(
    target_gene = target_gene,
    signature_genes = signature_genes,
    signature_name = signature_name,
    project_dir = project_dir,
    data_dir = data_dir,
    out_dir = out_dir,
    output_prefix = output_prefix,
    write_files = write_files,
    verbose = verbose
  )
  survival <- run_signature_survival(
    prepared = prepared,
    target_gene = target_gene,
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
    plots <- plot_signature_results(
      results = survival$results,
      target_gene = target_gene,
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
  target_expression <- NULL
  if (isTRUE(make_expression_plot) && "optimized" %in% cutoff_methods) {
    target_expression <- plot_target_expression_by_signature(
      survival = survival,
      target_gene = target_gene,
      signature_name = signature_name,
      cutoff_method = "optimized",
      out_dir = out_dir,
      output_prefix = output_prefix,
      write_files = write_files
    )
    qc_file <- file.path(out_dir, paste0(output_prefix, "_QC_report.txt"))
    if (isTRUE(write_files) && file.exists(qc_file)) {
      cat(
        "\nTarget expression outputs:\n",
        paste0("- ", target_expression$summary_tsv, "\n"),
        paste0("- ", target_expression$delta_pdf, "\n"),
        paste0("- ", target_expression$delta_png, "\n"),
        paste0("- ", target_expression$boxplot_pdf, "\n"),
        file = qc_file,
        append = TRUE
      )
    }
  }
  list(prepared = prepared, survival = survival, plots = plots, target_expression = target_expression)
}
