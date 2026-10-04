#' Analyze cell-type gene expression against multiple survival endpoints
#'
#' Groups are defined once within each cohort from all eligible tumor profiles,
#' before endpoint-specific missingness or covariate exclusions. The high group
#' is expression >= median; ties are never randomly split. PFS is accepted only
#' from explicit PFS columns and is never substituted with TCGA PFI.
#' @param data Prepared cell expression data.
#' @param target_gene,cell_type Exact gene and cell label.
#' @param endpoints Endpoint names with endpoint.time and endpoint.status columns.
#' @param covariates Additional clinical covariates. Constant covariates are
#'   recorded and omitted within a cohort; rank-deficient models are blocked.
#' @param context_subset all or high. high is a bulk sensitivity subset defined
#'   by within-cohort median cell_score, not isolated cells.
#' @param min_patients,min_events,min_group Minimum model eligibility thresholds.
#' @param allow_proxy Explicit permission to analyze a bulk proxy.
#' @param out_dir,prefix,write_files Output settings.
#' @param score_mode gene or legacy-compatible joint target/cell-context score.
#' @return Results including every requested cohort/endpoint/model, fixed
#'   groups and per-endpoint patient-level exclusions. Models are median high
#'   versus low and continuous expression per SD. Wald P and log-rank P are
#'   separate; BH correction includes all cancer/endpoint tests per model.
#' @export
run_cell_gene_survival <- function(data, target_gene = "ACLY", cell_type = "Treg",
    endpoints = c("OS", "PFI", "DSS", "DFI"), covariates = character(),
    context_subset = c("all", "high"), min_patients = 40L, min_events = 10L,
    min_group = 10L, allow_proxy = FALSE, out_dir = file.path(getwd(), "results"),
    prefix = "cell_gene", write_files = TRUE, score_mode = c("gene", "joint")) {
  tcgasig_cell_check(data, allow_proxy)
  context_subset <- match.arg(context_subset)
  score_mode <- match.arg(score_mode)
  for (n in c("min_patients", "min_events", "min_group")) tcgasig_validate_threshold(get(n), n)
  if (!length(endpoints) || anyNA(endpoints) || any(make.names(endpoints) != endpoints)) stop("Invalid endpoints.", call. = FALSE)
  endpoints <- unique(endpoints)
  covariates <- unique(covariates)
  groups <- tcgasig_cell_groups(data, target_gene, cell_type, context_subset, score_mode)
  results <- inputs <- list()
  for (cohort in unique(groups$cohort_label)) {
    initial <- groups[groups$cohort_label == cohort, , drop = FALSE]
    for (endpoint in endpoints) {
      d <- initial
      columns <- paste0(endpoint, c(".time", ".status"))
      available <- all(columns %in% names(d))
      d$endpoint_time <- if (available) tcgasig_to_numeric(d[[columns[1]]]) else NA_real_
      d$endpoint_status <- if (available) tcgasig_to_numeric(d[[columns[2]]]) else NA_real_
      d$exclusion <- ifelse(!d$context_selected, "outside_context_subset",
        ifelse(is.na(d$group), "gene_or_group_not_estimable",
          ifelse(!is.finite(d$endpoint_time) | d$endpoint_time <= 0 | !d$endpoint_status %in% c(0, 1),
            "endpoint_missing_or_invalid", "included")))
      if (!available) d$exclusion[d$exclusion == "endpoint_missing_or_invalid"] <- "endpoint_not_measured"
      candidate <- d[d$exclusion == "included", , drop = FALSE]
      cc <- tcgasig_covariate_rows(candidate, covariates)
      model_dt <- cc$data
      d$exclusion[d$exclusion == "included" & !d$profile_id %in% model_dt$profile_id] <- "covariate_missing"
      d$endpoint <- endpoint
      inputs[[paste(cohort, endpoint)]] <- d
      n <- nrow(model_dt); events <- sum(model_dt$endpoint_status == 1)
      high <- sum(model_dt$group == "high"); low <- sum(model_dt$group == "low")
      for (model in c("median", "continuous")) {
        row <- data.frame(cohort_label = cohort, cancer_type = initial$cancer_type[1],
          cell_type = cell_type, target_gene = target_gene, measurement = data$measurement,
          context_subset = context_subset, score_mode = score_mode, endpoint = endpoint, model = model,
          total_tumor_n = nrow(initial), context_n = sum(initial$context_selected),
          n_patients = n, n_events = events, high_n = high, low_n = low,
          gene_cutoff = initial$gene_cutoff[1], score_cutoff = initial$score_cutoff[1], context_cutoff = initial$context_cutoff[1],
          covariates_used = paste(cc$active, collapse = ";"),
          covariates_constant = paste(cc$constant, collapse = ";"),
          HR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, wald_p = NA_real_,
          logrank_p = NA_real_, PH_global_p = NA_real_, status = "ok", reason = "")
        reason <- if (!available) "endpoint_not_measured" else if (n < min_patients) "insufficient_patients" else
          if (events < min_events) "insufficient_events" else if (min(high, low) < min_group) "insufficient_or_tied_groups" else ""
        model_dt$gene_z <- tcgasig_zscore_or_na(model_dt$analysis_score)
        if (!nzchar(reason) && any(!is.finite(model_dt$gene_z))) reason <- "constant_target_expression"
        predictor <- if (model == "median") "group" else "gene_z"
        term <- if (model == "median") "grouphigh" else "gene_z"
        rhs <- paste(c(predictor, cc$active), collapse = " + ")
        if (!nzchar(reason)) {
          design <- stats::model.matrix(stats::as.formula(paste("~", rhs)), model_dt)
          if (qr(design)$rank < ncol(design)) reason <- "rank_deficient_design"
        }
        if (!nzchar(reason)) {
          # Reuse the package's convergence-checked Cox helper with endpoint-local aliases.
          model_dt$OS.time <- model_dt$endpoint_time
          model_dt$OS.status <- model_dt$endpoint_status
          fit <- tcgasig_fit_checked(model_dt, rhs, term)
          if (is.null(fit$error)) {
            for (v in names(fit$row)) row[[v]] <- fit$row[[v]]
            if (model == "median") {
              lr <- survival::survdiff(survival::Surv(endpoint_time, endpoint_status) ~ group, model_dt)
              row$logrank_p <- stats::pchisq(lr$chisq, df = 1, lower.tail = FALSE)
            }
          } else reason <- fit$error
        }
        if (nzchar(reason)) { row$status <- "not_estimable"; row$reason <- reason }
        results[[paste(cohort, endpoint, model)]] <- row
      }
    }
  }
  res <- as.data.frame(data.table::rbindlist(results, fill = TRUE))
  res$FDR <- NA_real_
  res$FDR_within_endpoint <- NA_real_
  for (model in unique(res$model)) {
    ii <- which(res$model == model)
    res$FDR[ii] <- stats::p.adjust(res$wald_p[ii], "BH")
    for (endpoint in endpoints) {
      jj <- ii[res$endpoint[ii] == endpoint]
      res$FDR_within_endpoint[jj] <- stats::p.adjust(res$wald_p[jj], "BH")
    }
  }
  result <- list(results = res, groups = groups,
    inputs = as.data.frame(data.table::rbindlist(inputs, fill = TRUE)),
    measurement = data$measurement, interpretation = tcgasig_cell_label(data),
    settings = list(target_gene = target_gene, cell_type = cell_type, endpoints = endpoints,
      covariates = covariates, context_subset = context_subset, score_mode = score_mode,
      direction = "high minus low; HR > 1 = higher hazard in the high-score group",
      cutoff = "cohort-specific median fixed before endpoint exclusions", time_unit = "days"))
  result$files <- if (write_files) tcgasig_write_cell_tables(result, out_dir, prefix) else character()
  result
}
