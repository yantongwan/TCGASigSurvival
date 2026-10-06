#' Test cell marker scores or estimated fractions against survival
#'
#' The predictor is abundance alone, not a gene plus cell-context joint score.
#' Cutoffs and standard deviations are fixed before endpoint/clinical exclusions.
#' High means strictly above median, with all ties assigned low. Constant/zero
#' predictors, tied groups and failed Cox fits remain explicit output rows.
#' @param abundance Long table with sample_barcode, patient_id, cohort_label,
#'   cell_type, value, measurement (marker_score or estimated_fraction), method,
#'   unit, denominator, and role (primary or auxiliary).
#' @param clinical Patient/sample/cohort table and endpoint.time/status columns.
#' @param endpoints Distinct endpoints; absent columns are not measured.
#' @param covariates Clinical adjustment columns.
#' @param adjustment Label for the adjustment set.
#' @param min_patients,min_events,min_group Eligibility thresholds.
#' @param min_events_per_parameter Exploratory events/design-parameter floor.
#' @param out_dir,write_files Output settings.
#' @return Results, fixed groups and explicit patient-level exclusion records.
#'   Continuous HR is per reference-cohort SD; estimated fractions additionally
#'   have HR per 10 percentage points. FDR spans both methods and five primary
#'   cells within cohort/endpoint/model/adjustment, separately from auxiliary cells.
#' @export
run_cell_abundance_survival <- function(abundance, clinical, endpoints = "OS",
    covariates = character(), adjustment = "unadjusted", min_patients = 40L,
    min_events = 10L, min_group = 10L, min_events_per_parameter = 5,
    out_dir = file.path(getwd(), "results"), write_files = TRUE) {
  a <- as.data.frame(abundance, stringsAsFactors = FALSE)
  required <- c("sample_barcode", "patient_id", "cohort_label", "cell_type", "value",
    "measurement", "method", "unit", "denominator", "role")
  if (!all(required %in% names(a)) || !nrow(a) || !is.numeric(a$value) ||
      anyNA(a[, setdiff(required, "value"), drop = FALSE]) ||
      any(!a$measurement %in% c("marker_score", "estimated_fraction")) ||
      any(!a$role %in% c("primary", "auxiliary"))) stop("Invalid abundance table/schema.", call. = FALSE)
  keys <- c("cohort_label", "measurement", "method", "cell_type", "patient_id")
  if (anyDuplicated(a[, keys]) || anyDuplicated(a[, setdiff(keys, "patient_id") , drop = FALSE] |> cbind(sample_barcode = a$sample_barcode))) {
    stop("Duplicate abundance rows/patient observations.", call. = FALSE)
  }
  fraction <- a$measurement == "estimated_fraction"
  if (any(is.finite(a$value[fraction]) & (a$value[fraction] < 0 | a$value[fraction] > 1))) stop("Fraction values must be in [0,1].", call. = FALSE)
  if (any(a$unit[fraction] != "fraction") || any(a$denominator[fraction] == "not_a_fraction")) stop("Declare fraction unit and denominator.", call. = FALSE)
  clinical <- tcgasig_abundance_clinical(clinical)
  ci <- match(paste(a$cohort_label, a$sample_barcode, sep = "\r"),
    paste(clinical$cohort_label, clinical$sample_barcode, sep = "\r"))
  if (anyNA(ci) || any(a$patient_id != clinical$patient_id[ci])) stop("Abundance and clinical patient/sample identities must match exactly.", call. = FALSE)
  endpoints <- unique(endpoints)
  if (!length(endpoints) || anyNA(endpoints) || any(!nzchar(endpoints)) || any(make.names(endpoints) != endpoints)) stop("Invalid endpoints.", call. = FALSE)
  covariates <- unique(covariates)
  if (!all(covariates %in% names(clinical)) || any(make.names(covariates) != covariates) ||
      any(covariates %in% c(required, "abundance_z", "group", "OS.time", "OS.status", "endpoint_time", "endpoint_status"))) stop("Use existing nonreserved clinical covariates.", call. = FALSE)
  for (v in c("min_patients", "min_events", "min_group")) tcgasig_validate_threshold(get(v), v)
  if (length(min_events_per_parameter) != 1 || !is.finite(min_events_per_parameter) || min_events_per_parameter < 0) stop("Invalid event/parameter threshold.", call. = FALSE)
  extra <- setdiff(names(clinical), names(a))
  d <- cbind(a, clinical[ci, extra, drop = FALSE])
  strata <- unique(a[, c("cohort_label", "measurement", "method", "cell_type")])
  rows <- inputs <- groups <- list()
  for (i in seq_len(nrow(strata))) {
    keep <- rep(TRUE, nrow(d))
    for (v in names(strata)) keep <- keep & d[[v]] == strata[[v]][i]
    base <- d[keep, , drop = FALSE]
    if (any(vapply(base[, c("unit", "denominator", "role"), drop = FALSE],
        function(v) length(unique(v)) != 1L, logical(1)))) {
      stop("Unit, denominator and role must be constant within each cell/method/cohort.", call. = FALSE)
    }
    valid_value <- is.finite(base$value)
    median <- if (any(valid_value)) stats::median(base$value[valid_value]) else NA_real_
    sd <- stats::sd(base$value[valid_value]); mean <- mean(base$value[valid_value])
    base$cutoff <- median; base$reference_sd <- sd
    base$group <- factor(ifelse(valid_value, ifelse(base$value > median, "high", "low"), NA), levels = c("low", "high"))
    base$abundance_z <- if (is.finite(sd) && sd > 0) (base$value - mean) / sd else NA_real_
    groups[[as.character(i)]] <- base[, c(required, "cutoff", "reference_sd", "group", "abundance_z")]
    for (endpoint in endpoints) {
      e <- base; cols <- paste0(endpoint, c(".time", ".status")); available <- all(cols %in% names(e))
      e$endpoint_time <- if (available) tcgasig_to_numeric(e[[cols[1]]]) else NA_real_
      e$endpoint_status <- if (available) tcgasig_to_numeric(e[[cols[2]]]) else NA_real_
      outcome_reason <- if (!available) rep("endpoint_not_measured", nrow(e)) else
        ifelse(!is.finite(e$endpoint_time) | e$endpoint_time <= 0 | !e$endpoint_status %in% c(0, 1),
          "endpoint_missing_or_invalid", "included")
      e$exclusion <- ifelse(!valid_value, "predictor_not_estimable", outcome_reason)
      candidate <- e[e$exclusion == "included", , drop = FALSE]
      cc <- tcgasig_covariate_rows(candidate, covariates)
      m <- cc$data
      e$exclusion[e$exclusion == "included" & !e$sample_barcode %in% m$sample_barcode] <- "covariate_missing"
      e$endpoint <- endpoint; e$adjustment <- adjustment
      inputs[[paste(i, endpoint)]] <- e
      n <- nrow(m); events <- sum(m$endpoint_status == 1); high <- sum(m$group == "high"); low <- sum(m$group == "low")
      for (model in c("continuous", "median")) {
        row <- data.frame(cohort_label = base$cohort_label[1], measurement = base$measurement[1],
          method = base$method[1], cell_type = base$cell_type[1], role = base$role[1],
          unit = base$unit[1], denominator = base$denominator[1], endpoint = endpoint,
          model = model, adjustment = adjustment, initial_n = nrow(base), n_patients = n,
          n_events = events, high_n = high, low_n = low, cutoff = median, reference_sd = sd,
          covariates_used = paste(cc$active, collapse = ";"), covariates_constant = paste(cc$constant, collapse = ";"),
          design_parameters = NA_integer_, HR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
          HR_per_10pp = NA_real_, CI_lower_10pp = NA_real_, CI_upper_10pp = NA_real_,
          wald_p = NA_real_, logrank_p = NA_real_, PH_global_p = NA_real_, PH_predictor_p = NA_real_,
          status = "ok", reason = "")
        reason <- if (!available) "endpoint_not_measured" else if (n < min_patients) "insufficient_patients" else
          if (events < min_events) "insufficient_events" else if (!is.finite(sd) || sd <= 0) "constant_or_unestimable_predictor" else
          if (model == "median" && min(high, low) < min_group) "insufficient_or_tied_groups" else ""
        predictor <- if (model == "continuous") "abundance_z" else "group"
        term <- if (model == "continuous") "abundance_z" else "grouphigh"
        rhs <- paste(c(predictor, cc$active), collapse = " + ")
        if (!nzchar(reason)) {
          mm <- tryCatch(stats::model.matrix(stats::as.formula(paste("~", rhs)), m), error = function(e) e)
          if (inherits(mm, "error")) reason <- conditionMessage(mm) else {
            row$design_parameters <- ncol(mm) - 1L
            if (qr(mm)$rank < ncol(mm)) reason <- "rank_deficient_design" else
              if (events < min_events_per_parameter * row$design_parameters) reason <- "insufficient_events_per_parameter"
          }
        }
        if (!nzchar(reason)) {
          m$OS.time <- m$endpoint_time; m$OS.status <- m$endpoint_status
          fit <- tcgasig_fit_checked(m, rhs, term)
          if (is.null(fit$error)) {
            for (v in names(fit$row)) row[[v]] <- fit$row[[v]]
            ph <- tryCatch({
              ph_table <- survival::cox.zph(fit$fit)$table
              ph_table[if (term %in% rownames(ph_table)) term else predictor, "p"]
            }, error = function(e) NA_real_)
            row$PH_predictor_p <- ph
            if (model == "continuous" && row$measurement == "estimated_fraction") {
              row$HR_per_10pp <- exp(log(row$HR) * 0.1 / sd)
              row$CI_lower_10pp <- exp(log(row$CI_lower) * 0.1 / sd)
              row$CI_upper_10pp <- exp(log(row$CI_upper) * 0.1 / sd)
              if (any(!is.finite(unlist(row[, c("HR_per_10pp", "CI_lower_10pp", "CI_upper_10pp")])))) {
                row$HR_per_10pp <- row$CI_lower_10pp <- row$CI_upper_10pp <- NA_real_
              }
            }
            if (model == "median") {
              lr <- survival::survdiff(survival::Surv(endpoint_time, endpoint_status) ~ group, data = m)
              row$logrank_p <- stats::pchisq(lr$chisq, 1, lower.tail = FALSE)
            }
          } else reason <- fit$error
        }
        if (nzchar(reason)) { row$status <- "not_estimable"; row$reason <- reason }
        rows[[paste(i, endpoint, model)]] <- row
      }
    }
  }
  result <- as.data.frame(data.table::rbindlist(rows))
  result$FDR <- result$FDR_within_method <- NA_real_
  family <- interaction(result$cohort_label, result$endpoint, result$model, result$adjustment, result$role, drop = TRUE)
  for (f in unique(family)) {
    ii <- which(family == f); result$FDR[ii] <- stats::p.adjust(result$wald_p[ii], "BH", n = length(ii))
    for (method in unique(result$measurement[ii])) {
      jj <- ii[result$measurement[ii] == method]
      result$FDR_within_method[jj] <- stats::p.adjust(result$wald_p[jj], "BH", n = length(jj))
    }
  }
  out <- list(results = result, inputs = as.data.frame(data.table::rbindlist(inputs, fill = TRUE)),
    groups = as.data.frame(data.table::rbindlist(groups)), settings = list(
      cutoff = "strictly above reference-cohort median is high; ties assigned low; no endpoint-guided cutoffs",
      continuous_effect = "per 1 reference-cohort SD; fractions additionally per +10 percentage points",
      FDR = "both methods within cohort/endpoint/model/adjustment/primary or auxiliary family",
      covariates = covariates, min_events_per_parameter = min_events_per_parameter))
  if (write_files) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    for (v in c("results", "inputs", "groups")) data.table::fwrite(out[[v]], file.path(out_dir, paste0("abundance_", v, ".tsv")), sep = "\t", na = "NA")
  }
  out
}
