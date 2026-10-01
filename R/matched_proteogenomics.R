#' Join RNA, protein and clinical tables within matching cancer-patient keys
#' @param rna Data frame with cancer_type, patient_id and target_rna;
#'   optional signature_score.
#' @param protein Data frame with cancer_type, patient_id and target_protein.
#' @param clinical Data frame with cancer_type, patient_id, OS.time and OS.status.
#' @return A list with data (strict RNA-protein intersection, clinical left join)
#'   and qc input/intersection counts. No identifier rewriting or cross-cohort
#'   patient matching is performed.
#' @export
match_proteogenomic_tables <- function(rna, protein, clinical) {
  tables <- list(rna = as.data.frame(rna), protein = as.data.frame(protein), clinical = as.data.frame(clinical))
  required <- list(rna = c("target_rna"), protein = "target_protein", clinical = c("OS.time", "OS.status"))
  for (name in names(tables)) {
    dt <- tables[[name]]
    if (!all(c("cancer_type", "patient_id", required[[name]]) %in% names(dt))) {
      stop(name, " lacks required columns.", call. = FALSE)
    }
    if (anyNA(dt[, c("cancer_type", "patient_id")]) ||
        any(!nzchar(as.character(dt$cancer_type))) || any(!nzchar(as.character(dt$patient_id))) ||
        anyDuplicated(dt[, c("cancer_type", "patient_id")])) {
      stop(name, " must contain unique, nonempty cancer-patient keys.", call. = FALSE)
    }
  }
  keys <- c("cancer_type", "patient_id")
  if (length(intersect(setdiff(names(tables$rna), keys), setdiff(names(tables$protein), keys))) ||
      length(intersect(setdiff(c(names(tables$rna), names(tables$protein)), keys), setdiff(names(tables$clinical), keys)))) {
    stop("Non-key columns overlap between input tables.", call. = FALSE)
  }
  joined <- merge(tables$rna, tables$protein, by = keys, all = FALSE)
  joined <- merge(joined, tables$clinical, by = keys, all.x = TRUE)
  list(data = joined, qc = data.frame(source = c(names(tables), "rna_protein_matched"),
    n = c(vapply(tables, nrow, integer(1)), nrow(joined))))
}

#' Analyze an independently matched proteogenomic cohort in R
#' @param data A matched data frame, delimited file path or the output of
#'   match_proteogenomic_tables. Columns: cancer_type, patient_id, target_rna,
#'   target_protein, OS.time (days), OS.status (0/1); optional signature_score.
#' @param target_gene Human gene symbol used for labeling.
#' @param signature_name Name of the optional RNA marker module.
#' @param out_dir Output directory.
#' @param output_prefix Output basename.
#' @param min_patients,min_events Minimum patients/deaths for survival models.
#' @param covariates Extra clinical columns to adjust for.
#' @param write_files Write tables and QC.
#' @param make_plots Write an RNA-protein scatter plot per cohort.
#' @return RNA-protein Spearman correlations, continuous/median Cox results,
#'   skipped models and missingness. Signature adjustment is reported as a
#'   separate model. All inputs must already represent the same patient within
#'   the same study; cancer labels alone cannot establish patient matching.
#' @export
run_matched_proteogenomic_analysis <- function(data, target_gene = "STRAP",
    signature_name = "Signature", out_dir = file.path(getwd(), "results"),
    output_prefix = NULL, min_patients = 80L, min_events = 20L,
    covariates = character(), write_files = TRUE, make_plots = TRUE) {
  if (is.character(data)) data <- tcgasig_read_delimited(data)
  if (is.list(data) && !is.data.frame(data) && !is.null(data$data)) data <- data$data
  data <- as.data.frame(data)
  required <- c("cancer_type", "patient_id", "target_rna", "target_protein", "OS.time", "OS.status")
  if (!all(required %in% names(data))) stop("Matched input requires: ", paste(required, collapse = ", "), call. = FALSE)
  if (anyNA(data[, c("cancer_type", "patient_id")]) ||
      any(!nzchar(as.character(data$cancer_type))) || any(!nzchar(as.character(data$patient_id))) ||
      anyDuplicated(data[, c("cancer_type", "patient_id")])) stop("Matched cancer-patient keys must be unique and nonempty.", call. = FALSE)
  tcgasig_validate_threshold(min_patients, "min_patients")
  tcgasig_validate_threshold(min_events, "min_events")
  if (!all(covariates %in% names(data)) || any(make.names(covariates) != covariates) ||
      any(covariates %in% c(required, "signature_score", "value_z", "signature_z", "group"))) {
    stop("Invalid covariates.", call. = FALSE)
  }
  for (name in required[3:6]) data[[name]] <- tcgasig_to_numeric(data[[name]])
  if ("signature_score" %in% names(data)) data$signature_score <- tcgasig_to_numeric(data$signature_score)
  prefix <- if (is.null(output_prefix)) paste0("CPTAC_", tcgasig_safe_name(target_gene), "_", tcgasig_safe_name(signature_name)) else tcgasig_safe_name(output_prefix)
  correlations <- list(); models <- list(); skipped <- list(); qc <- list()
  for (cancer in unique(data$cancer_type)) {
    dt <- data[data$cancer_type == cancer, , drop = FALSE]
    paired <- is.finite(dt$target_rna) & is.finite(dt$target_protein)
    valid_os <- is.finite(dt$OS.time) & dt$OS.time > 0 & dt$OS.status %in% c(0, 1)
    qc[[cancer]] <- data.frame(cancer_type = cancer, n_matched = nrow(dt),
      n_rna_protein_complete = sum(paired), n_survival_complete = sum(valid_os),
      n_rna_missing = sum(!is.finite(dt$target_rna)), n_protein_missing = sum(!is.finite(dt$target_protein)))
    pair_dt <- dt[paired, ]
    if (nrow(pair_dt) >= 4L && stats::sd(pair_dt$target_rna) > 0 && stats::sd(pair_dt$target_protein) > 0) {
      corr <- stats::cor.test(pair_dt$target_rna, pair_dt$target_protein, method = "spearman", exact = FALSE)
      correlations[[cancer]] <- data.frame(cancer_type = cancer, n = nrow(pair_dt),
        spearman_rho = unname(corr$estimate), spearman_p = corr$p.value)
    } else skipped[[length(skipped) + 1L]] <- data.frame(cancer_type = cancer,
      predictor = "RNA-protein", model = "correlation", reason = "Insufficient pairs or constant values.")
    for (predictor in c("target_rna", "target_protein")) {
      for (adjustment in c("unadjusted", if ("signature_score" %in% names(dt)) "signature_adjusted")) {
        needed <- c(predictor, covariates, if (adjustment == "signature_adjusted") "signature_score")
        valid <- valid_os & stats::complete.cases(dt[, needed, drop = FALSE])
        for (column in needed) if (is.numeric(dt[[column]])) valid <- valid & is.finite(dt[[column]])
        d <- dt[valid, , drop = FALSE]
        n <- nrow(d); events <- sum(d$OS.status == 1)
        d$value_z <- tcgasig_zscore_or_na(d[[predictor]])
        if (adjustment == "signature_adjusted") d$signature_z <- tcgasig_zscore_or_na(d$signature_score)
        for (method in c("continuous", "median")) {
          reason <- NULL
          if (n < min_patients || events < min_events) reason <- "Below patient/event thresholds."
          else if (any(!is.finite(d$value_z)) || (adjustment == "signature_adjusted" && any(!is.finite(d$signature_z)))) reason <- "Constant predictor or signature."
          if (!is.null(reason)) {
            skipped[[length(skipped) + 1L]] <- data.frame(cancer_type = cancer, predictor = predictor,
              model = paste(method, adjustment), reason = reason); next
          }
          if (method == "median") d$group <- factor(ifelse(d$value_z >= stats::median(d$value_z), "high", "low"), levels = c("low", "high"))
          terms <- c(if (method == "continuous") "value_z" else "group",
            if (adjustment == "signature_adjusted") "signature_z", covariates)
          fit <- tcgasig_fit_checked(d, paste(terms, collapse = " + "), if (method == "continuous") "value_z" else "grouphigh")
          if (!is.null(fit$error)) skipped[[length(skipped) + 1L]] <- data.frame(cancer_type = cancer,
            predictor = predictor, model = paste(method, adjustment), reason = fit$error)
          else models[[length(models) + 1L]] <- cbind(data.frame(cancer_type = cancer, predictor = predictor,
            method = method, adjustment = adjustment, n_patients = n, n_events = events), fit$row)
        }
      }
    }
  }
  combine <- function(x, template) if (length(x)) data.table::rbindlist(x, fill = TRUE) else template
  result <- list(
    correlations = combine(correlations, data.table::data.table(cancer_type = character(), n = integer(), spearman_rho = numeric(), spearman_p = numeric())),
    survival = combine(models, data.table::data.table(cancer_type = character(), predictor = character(), method = character(), adjustment = character(),
      n_patients = integer(), n_events = integer(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), wald_p = numeric(), PH_global_p = numeric())),
    skipped = combine(skipped, data.table::data.table(cancer_type = character(), predictor = character(), model = character(), reason = character())),
    qc = combine(qc, data.table::data.table(cancer_type = character(), n_matched = integer(), n_rna_protein_complete = integer(), n_survival_complete = integer(),
      n_rna_missing = integer(), n_protein_missing = integer())), target_gene = target_gene, signature_name = signature_name)
  if (nrow(result$correlations)) result$correlations$spearman_FDR <- stats::p.adjust(result$correlations$spearman_p, "BH")
  if (nrow(result$survival)) result$survival[, FDR := stats::p.adjust(wald_p, "BH"), by = .(predictor, method, adjustment)]
  if (write_files || make_plots) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (write_files) for (name in c("correlations", "survival", "skipped", "qc")) {
    data.table::fwrite(result[[name]], file.path(out_dir, paste0(prefix, "_", name, ".tsv")), sep = "\t")
  }
  if (make_plots) {
    d <- data[is.finite(data$target_rna) & is.finite(data$target_protein), ]
    if (nrow(d)) {
      p <- ggplot2::ggplot(d, ggplot2::aes(x = target_rna, y = target_protein)) +
        ggplot2::geom_point(alpha = 0.6, size = 1.3, colour = "#287B8E") +
        ggplot2::facet_wrap(~cancer_type, scales = "free") + ggplot2::theme_bw() +
        ggplot2::labs(x = paste(target_gene, "normalized RNA"), y = paste(target_gene, "normalized protein"))
      ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_rna_protein.pdf")), p, width = 8, height = 6)
    }
  }
  result
}
