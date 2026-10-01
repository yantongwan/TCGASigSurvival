default_cptac_tcga_map <- function() {
  data.table::data.table(
    cptac_cancer = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"),
    tcga_cancer = c("BRCA", "KIRC", "COAD", "GBM", "HNSC", "LUSC", "LUAD", "OV", "PAAD", "UCEC"),
    note = c(
      "CPTAC breast cancer maps to TCGA BRCA.",
      "CPTAC clear cell renal cell carcinoma maps to TCGA KIRC.",
      "CPTAC colon adenocarcinoma maps to TCGA COAD.",
      "CPTAC glioblastoma maps to TCGA GBM.",
      "CPTAC head and neck squamous cell carcinoma maps to TCGA HNSC.",
      "CPTAC lung squamous cell carcinoma maps to TCGA LUSC.",
      "CPTAC lung adenocarcinoma maps to TCGA LUAD.",
      "CPTAC ovarian cancer maps to TCGA OV.",
      "CPTAC pancreatic ductal adenocarcinoma maps to TCGA PAAD.",
      "CPTAC uterine corpus endometrial carcinoma maps to TCGA UCEC."
    )
  )
}

default_cd8_signature <- function() {
  c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GZMA", "IFNG")
}

tcgasig_cptac_signature_genes <- function(signature_name, signature_genes = NULL) {
  if (!is.null(signature_genes)) return(unique(toupper(signature_genes)))
  if (toupper(signature_name) == "CD8") return(default_cd8_signature())
  default_treg_signature()
}

tcgasig_cptac_default_dir <- function(project_dir) {
  file.path(project_dir, "PanCancer_STRAP_Proteogenomics")
}

tcgasig_existing_file_info <- function(paths) {
  data.table::rbindlist(lapply(names(paths), function(key) {
    path <- paths[[key]]
    exists <- file.exists(path)
    info <- if (exists) file.info(path) else NULL
    data.table::data.table(
      file_key = key,
      path = path,
      exists = exists,
      size_bytes = if (exists) as.numeric(info$size) else NA_real_,
      mtime = if (exists) as.character(info$mtime) else NA_character_
    )
  }), fill = TRUE)
}

tcgasig_read_required_table <- function(path, label) {
  if (!file.exists(path)) tcgasig_stop("Missing CPTAC ", label, " file: ", path)
  tcgasig_read_tsv_auto(path)
}

tcgasig_read_optional_table <- function(path) {
  if (!file.exists(path)) return(data.table::data.table())
  tcgasig_read_tsv_auto(path)
}

tcgasig_prefer_existing_path <- function(primary, fallback = NULL) {
  if (file.exists(primary) || is.null(fallback)) return(primary)
  fallback
}

tcgasig_load_cptac_joint_tables <- function(proteogenomics_dir, target_gene, signature_name = "Treg") {
  joint_dir <- file.path(proteogenomics_dir, "results", "joint_tables")
  if (!dir.exists(joint_dir)) return(data.table::data.table())
  pattern <- paste0("_", tcgasig_safe_name(target_gene), "_", tcgasig_safe_name(signature_name), "_joint\\.tsv\\.gz$")
  files <- list.files(joint_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0 && identical(signature_name, "Treg")) {
    pattern <- paste0("_", tcgasig_safe_name(target_gene), "_joint\\.tsv\\.gz$")
    files <- list.files(joint_dir, pattern = pattern, full.names = TRUE)
  }
  if (length(files) == 0) return(data.table::data.table())
  tables <- lapply(files, function(path) {
    dt <- tcgasig_read_tsv_auto(path)
    dt[, source_file := path]
    dt
  })
  data.table::rbindlist(tables, fill = TRUE)
}

load_cptac_proteogenomic_results <- function(
    project_dir = getwd(),
    proteogenomics_dir = tcgasig_cptac_default_dir(project_dir),
    target_gene = "STRAP",
    signature_name = "Treg",
    signature_genes = NULL,
    include_joint_tables = FALSE,
    verbose = TRUE) {
  if (!dir.exists(proteogenomics_dir)) {
    tcgasig_stop("Missing proteogenomics_dir: ", proteogenomics_dir)
  }

  results_dir <- file.path(proteogenomics_dir, "results")
  qc_dir <- file.path(results_dir, "qc")
  tables_dir <- file.path(results_dir, "tables")
  reports_dir <- file.path(proteogenomics_dir, "reports")
  target_gene <- toupper(target_gene)
  signature_genes <- tcgasig_cptac_signature_genes(signature_name, signature_genes)

  old_strap <- identical(target_gene, "STRAP")
  old_treg <- identical(signature_name, "Treg")
  prefix <- paste(target_gene, signature_name, sep = "_")
  paths <- c(
    audit = tcgasig_prefer_existing_path(file.path(qc_dir, paste0("cptac_", prefix, "_data_audit.tsv")), file.path(qc_dir, "cptac_data_audit.tsv")),
    joint_build_qc = tcgasig_prefer_existing_path(file.path(qc_dir, paste0("cptac_", prefix, "_joint_table_build_qc.tsv")), file.path(qc_dir, "cptac_joint_table_build_qc.tsv")),
    missingness = tcgasig_prefer_existing_path(file.path(qc_dir, paste0("cptac_", prefix, "_missingness.tsv")), file.path(qc_dir, "cptac_missingness.tsv")),
    rna_protein_correlations = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("cptac_", target_gene, "_rna_protein_correlations.tsv"))),
    rna_protein_meta_analysis = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("cptac_", target_gene, "_rna_protein_meta_analysis.tsv"))),
    treg_associations = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("cptac_", prefix, "_associations.tsv")), if (old_treg) file.path(tables_dir, paste0("cptac_", target_gene, "_Treg_associations.tsv")) else NULL),
    decoupling_treg_associations = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("cptac_", target_gene, "_decoupling_", signature_name, "_associations.tsv")), if (old_treg && old_strap) file.path(tables_dir, "cptac_decoupling_Treg_associations.tsv") else NULL),
    survival_univariable = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("cptac_", prefix, "_survival_univariable.tsv")), if (old_treg && old_strap) file.path(tables_dir, "cptac_survival_univariable.tsv") else NULL),
    survival_multivariable = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("cptac_", prefix, "_survival_multivariable.tsv")), if (old_treg && old_strap) file.path(tables_dir, "cptac_survival_multivariable.tsv") else NULL),
    survival_skipped = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("cptac_", prefix, "_survival_skipped.tsv")), if (old_treg && old_strap) file.path(tables_dir, "cptac_survival_skipped.tsv") else NULL),
    survival_model_comparison = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("cptac_", prefix, "_survival_model_comparison.tsv")), if (old_treg && old_strap) file.path(tables_dir, "cptac_survival_model_comparison.tsv") else NULL),
    integrated_summary = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("tcga_cptac_", prefix, "_integrated_summary.tsv")), if (old_treg && old_strap) file.path(tables_dir, "tcga_cptac_STRAP_integrated_summary.tsv") else NULL),
    replication_status = tcgasig_prefer_existing_path(file.path(tables_dir, paste0("tcga_cptac_", prefix, "_replication_status.tsv")), if (old_treg && old_strap) file.path(tables_dir, "tcga_cptac_replication_status.tsv") else NULL),
    final_report = tcgasig_prefer_existing_path(file.path(reports_dir, paste0("final_analysis_report_", prefix, ".md")), if (old_treg && old_strap) file.path(reports_dir, "final_analysis_report.md") else NULL),
    methods_report = file.path(reports_dir, "methods_and_interpretation.md")
  )

  tcgasig_message(verbose, "Loading CPTAC proteogenomic validation results from: ", proteogenomics_dir)
  out <- list(
    project_dir = project_dir,
    proteogenomics_dir = proteogenomics_dir,
    target_gene = target_gene,
    signature_name = signature_name,
    signature_genes = signature_genes,
    tcga_cptac_map = default_cptac_tcga_map(),
    source_files = tcgasig_existing_file_info(paths),
    audit = tcgasig_read_required_table(paths[["audit"]], "audit"),
    joint_build_qc = tcgasig_read_optional_table(paths[["joint_build_qc"]]),
    missingness = tcgasig_read_optional_table(paths[["missingness"]]),
    rna_protein_correlations = tcgasig_read_required_table(paths[["rna_protein_correlations"]], "RNA-protein correlation"),
    rna_protein_meta_analysis = tcgasig_read_optional_table(paths[["rna_protein_meta_analysis"]]),
    treg_associations = tcgasig_read_optional_table(paths[["treg_associations"]]),
    decoupling_treg_associations = tcgasig_read_optional_table(paths[["decoupling_treg_associations"]]),
    survival_univariable = tcgasig_read_optional_table(paths[["survival_univariable"]]),
    survival_multivariable = tcgasig_read_optional_table(paths[["survival_multivariable"]]),
    survival_skipped = tcgasig_read_optional_table(paths[["survival_skipped"]]),
    survival_model_comparison = tcgasig_read_optional_table(paths[["survival_model_comparison"]]),
    integrated_summary = tcgasig_read_optional_table(paths[["integrated_summary"]]),
    replication_status = tcgasig_read_optional_table(paths[["replication_status"]]),
    joint_tables = if (isTRUE(include_joint_tables)) tcgasig_load_cptac_joint_tables(proteogenomics_dir, target_gene, signature_name) else data.table::data.table(),
    interpretation_notes = c(
      "CPTAC is an independent proteogenomic validation cohort, not patient-level TCGA data.",
      "RNA-protein matching is performed only within the same CPTAC patient/cohort.",
      "TCGA-CPTAC integration is at the cancer/effect-size level only.",
      paste0(signature_name, "_score is a bulk RNA marker-module signal, not an immune-cell fraction or cell-intrinsic ", target_gene, " measurement."),
      "Optimized cutoffs are exploratory; continuous Cox and median cutoffs are the more stable primary views."
    )
  )
  class(out) <- c("tcgasig_cptac_proteogenomics", class(out))
  out
}

summarize_cptac_proteogenomic_results <- function(results, fdr_cutoff = 0.05) {
  if (is.character(results) && length(results) == 1L) {
    results <- load_cptac_proteogenomic_results(proteogenomics_dir = results, verbose = FALSE)
  }
  if (!inherits(results, "tcgasig_cptac_proteogenomics") && !is.list(results)) {
    tcgasig_stop("results must be a load_cptac_proteogenomic_results() object or a proteogenomics directory path.")
  }

  audit <- data.table::as.data.table(results$audit)
  corr <- data.table::as.data.table(results$rna_protein_correlations)
  surv <- data.table::as.data.table(results$survival_univariable)
  skipped <- data.table::as.data.table(results$survival_skipped)
  integ <- data.table::as.data.table(results$integrated_summary)

  qc_metrics <- data.table::data.table(
    metric = c(
      "audited_cptac_cancers",
      "audit_ok_or_warning",
      "rna_protein_correlation_cohorts",
      "rna_protein_fdr_significant_cohorts",
      "survival_result_rows",
      "survival_skipped_rows",
      "integrated_tcga_cptac_cancers"
    ),
    value = c(
      nrow(audit),
      sum(audit$status %in% c("ok", "warning"), na.rm = TRUE),
      nrow(corr),
      if ("spearman_FDR" %in% names(corr)) sum(corr$spearman_FDR < fdr_cutoff, na.rm = TRUE) else NA_integer_,
      nrow(surv),
      nrow(skipped),
      nrow(integ)
    )
  )

  cohort_inventory <- if (nrow(audit) > 0) {
    keep <- intersect(
      c("cptac_cancer", "tcga_cancer", "target_gene", "status", "reason", "n_matched_ids",
        "n_target_rna_protein_complete", "n_strap_rna_protein_complete", "n_signature_score_complete", "n_treg_score_complete", "n_survival_complete"),
      names(audit)
    )
    audit[, ..keep]
  } else {
    data.table::data.table()
  }

  concordance <- if (nrow(corr) > 0) {
    keep <- intersect(
      c("cptac_cancer", "tcga_cancer", "n", "spearman_rho", "spearman_ci_lower",
        "spearman_ci_upper", "spearman_p", "spearman_FDR"),
      names(corr)
    )
    corr[, ..keep][order(spearman_FDR)]
  } else {
    data.table::data.table()
  }

  survival_status <- if (nrow(skipped) > 0) {
    skipped[, .(
      skipped_models = .N,
      skipped_reasons = paste(unique(reason), collapse = " | ")
    ), by = .(cptac_cancer, tcga_cancer)]
  } else {
    data.table::data.table()
  }

  replication <- if (nrow(integ) > 0) {
    target_rho_col <- paste0("CPTAC_", results$target_gene, "_RNA_protein_rho")
    tcga_hr_col <- paste0("TCGA_RNA_", results$signature_name, "_HR_median")
    cptac_rna_hr_col <- paste0("CPTAC_RNA_", results$signature_name, "_HR")
    cptac_protein_hr_col <- paste0("CPTAC_protein_", results$signature_name, "_HR")
    keep <- intersect(
      c("tcga_cancer", "cptac_cancer", tcga_hr_col, cptac_rna_hr_col, cptac_protein_hr_col,
        "TCGA_RNA_Treg_HR_median", "CPTAC_RNA_Treg_HR", "CPTAC_protein_Treg_HR",
        target_rho_col, "CPTAC_STRAP_RNA_protein_rho", "replication_status", "not_testable_reason"),
      names(integ)
    )
    integ[, ..keep]
  } else {
    data.table::data.table()
  }

  list(
    qc_metrics = qc_metrics,
    cohort_inventory = cohort_inventory,
    rna_protein_concordance = concordance,
    survival_status = survival_status,
    tcga_cptac_replication = replication,
    interpretation_notes = results$interpretation_notes
  )
}

run_cptac_proteogenomic_validation <- function(
    project_dir = getwd(),
    proteogenomics_dir = tcgasig_cptac_default_dir(project_dir),
    target_gene = "STRAP",
    signature_name = "Treg",
    signature_genes = NULL,
    rerun = FALSE,
    runner = file.path(proteogenomics_dir, "scripts", "run_all.sh"),
    runner_args = "--all",
    include_joint_tables = FALSE,
    verbose = TRUE) {
  if (isTRUE(rerun)) {
    if (!file.exists(runner)) tcgasig_stop("Missing CPTAC proteogenomics runner: ", runner)
    runner <- normalizePath(runner, mustWork = TRUE)
    proteogenomics_dir <- normalizePath(proteogenomics_dir, mustWork = TRUE)
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(proteogenomics_dir)
    tcgasig_message(verbose, "Running CPTAC proteogenomic validation pipeline: bash ", runner, " ", paste(runner_args, collapse = " "))
    signature_genes <- tcgasig_cptac_signature_genes(signature_name, signature_genes)
    env <- c(
      paste0("TARGET_GENE=", shQuote(toupper(target_gene))),
      paste0("SIGNATURE_NAME=", shQuote(signature_name)),
      paste0("SIGNATURE_GENES=", shQuote(paste(signature_genes, collapse = ",")))
    )
    status <- system2("bash", c(shQuote(runner), vapply(runner_args, shQuote, character(1))), env = env)
    if (!identical(status, 0L)) {
      tcgasig_stop("CPTAC proteogenomic validation runner failed with exit status: ", status)
    }
  }

  results <- load_cptac_proteogenomic_results(
    project_dir = project_dir,
    proteogenomics_dir = proteogenomics_dir,
    target_gene = target_gene,
    signature_name = signature_name,
    signature_genes = signature_genes,
    include_joint_tables = include_joint_tables,
    verbose = verbose
  )
  results$summary <- summarize_cptac_proteogenomic_results(results)
  results
}

print.tcgasig_cptac_proteogenomics <- function(x, ...) {
  audit <- data.table::as.data.table(x$audit)
  corr <- data.table::as.data.table(x$rna_protein_correlations)
  cat("TCGASigSurvival CPTAC proteogenomic validation results\n")
  cat("  proteogenomics_dir:", x$proteogenomics_dir, "\n")
  cat("  target_gene:", x$target_gene, "\n")
  cat("  signature_name:", x$signature_name, "\n")
  if (nrow(audit) > 0) {
    cat("  audited CPTAC cancers:", nrow(audit), "\n")
    cat("  ok/warning cancers:", sum(audit$status %in% c("ok", "warning"), na.rm = TRUE), "\n")
  }
  if (nrow(corr) > 0 && "spearman_FDR" %in% names(corr)) {
    cat("  RNA-protein correlation cohorts:", nrow(corr), "\n")
    cat("  FDR < 0.05 cohorts:", sum(corr$spearman_FDR < 0.05, na.rm = TRUE), "\n")
  }
  invisible(x)
}
