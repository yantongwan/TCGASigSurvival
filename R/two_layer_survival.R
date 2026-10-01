tcgasig_validate_threshold <- function(x, label, minimum = 1) {
  if (length(x) != 1L || !is.finite(x) || x < minimum || x != as.integer(x)) {
    stop(label, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
}

tcgasig_score_module <- function(expr, genes, samples, minimum) {
  genes <- unique(tcgasig_clean_missing(genes))
  genes <- genes[!is.na(genes)]
  used <- character()
  z <- list()
  for (gene in intersect(genes, expr$gene_name)) {
    values <- tcgasig_extract_gene_vector(expr, gene, samples)
    values[!is.finite(values)] <- NA_real_
    scored <- tcgasig_zscore_or_na(values)
    if (any(is.finite(scored))) {
      z[[gene]] <- scored
      used <- c(used, gene)
    }
  }
  score <- rep(NA_real_, length(samples))
  counts <- integer(length(samples))
  if (length(z)) {
    mat <- do.call(cbind, z)
    counts <- rowSums(is.finite(mat))
    score <- rowMeans(mat, na.rm = TRUE)
    score[counts < minimum | !is.finite(score)] <- NA_real_
  }
  list(score = score, counts = counts, used = used, missing = setdiff(genes, used))
}

tcgasig_fit_checked <- function(dt, rhs, term) {
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    survival::coxph(stats::as.formula(paste("survival::Surv(OS.time, OS.status) ~", rhs)),
                    data = dt, x = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) return(list(error = conditionMessage(fit)))
  if (length(warnings)) return(list(error = paste(unique(warnings), collapse = " | ")))
  s <- summary(fit)
  if (!term %in% rownames(s$coefficients)) return(list(error = "Requested coefficient absent."))
  ci <- s$conf.int[term, ]
  p <- s$coefficients[term, "Pr(>|z|)"]
  values <- unname(c(ci["exp(coef)"], ci["lower .95"], ci["upper .95"], p))
  if (any(!is.finite(values)) || any(values[1:3] <= 0)) {
    return(list(error = "Nonfinite or nonpositive hazard ratio/confidence interval."))
  }
  ph <- tryCatch(survival::cox.zph(fit)$table["GLOBAL", "p"], error = function(e) NA_real_)
  list(error = NULL, row = data.frame(HR = values[1], CI_lower = values[2],
       CI_upper = values[3], wald_p = values[4], PH_global_p = unname(ph)), fit = fit)
}

#' Run cell-context and state-signature survival models
#'
#' Per-cohort, per-gene z scores are averaged independently for the two
#' modules. Model A compares state high/low while adjusting for continuous cell
#' score; model B compares state high/low within the cell-score-high subset;
#' model C tests continuous state score per SD adjusted for cell score.
#' Models with convergence warnings or nonfinite estimates are logged and
#' excluded from inference. Optimized cutoffs are exploratory and their Wald
#' p values do not correct for cutoff searching.
#' @param prepared Prepared-data list or RDS path.
#' @param state_signature_genes Human symbols for the state module.
#' @param cell_signature_genes Human symbols for the cell-context module.
#' @param state_name Name of the state signature.
#' @param cell_name Name of the cell context.
#' @param cell_scores Optional data frame with sample_barcode and cell_score
#'   from an independently computed cell abundance/deconvolution method.
#'   IDs must be unique and match exactly; missing scores are excluded.
#' @param cell_score_method Label documenting the external method when
#'   cell_scores is supplied, otherwise marker_module.
#' @param covariates Additional clinical columns to adjust for in all models.
#' @param out_dir Output directory.
#' @param output_prefix Safe basename for output files.
#' @param min_state_genes,min_cell_genes Minimum usable genes per sample.
#' @param min_patients,min_events Minimum full-cohort patients and deaths.
#' @param min_subset_patients,min_subset_events Minimum cell-high patients and deaths.
#' @param cutoff_methods Median and/or optimized cutoffs.
#' @param optimized_minprop Minimum group fraction for optimized cutoffs.
#' @param include_skcm_primary_only Include the SKCM primary-only sensitivity cohort.
#' @param make_km Write unadjusted Kaplan-Meier plots for the cell-high subset.
#' @param make_summary_plots Write forest plots for group comparisons.
#' @param write_files Write tables, QC and session information.
#' @param verbose Emit progress messages.
#' @return A list containing scores, adjusted, continuous, cell_high, skipped,
#'   gene_qc, cohort_qc, files and settings. FDR families are separate by
#'   model, cutoff method and default_cohort; sensitivity cohorts are separate.
#' @examples
#' demo <- tcgasig_demo_data(n = 80)
#' prep <- prepare_expression_survival_data(demo$expression, demo$clinical)
#' result <- run_two_layer_survival(prep, demo$state_genes,
#'   cell_signature_genes = demo$cell_genes, cutoff_methods = "median",
#'   min_patients = 40, min_events = 10, min_subset_patients = 20,
#'   min_subset_events = 5, make_km = FALSE, make_summary_plots = FALSE,
#'   write_files = FALSE, verbose = FALSE)
#' result$continuous
#' @export
run_two_layer_survival <- function(prepared, state_signature_genes,
    cell_signature_genes = default_naive_cd4_signature(),
    state_name = "State", cell_name = "NaiveCD4", cell_scores = NULL,
    cell_score_method = "marker_module", covariates = character(),
    out_dir = file.path(getwd(), "results"), output_prefix = NULL,
    min_state_genes = 6L, min_cell_genes = 6L, min_patients = 80L,
    min_events = 20L, min_subset_patients = 40L, min_subset_events = 10L,
    cutoff_methods = c("median", "optimized"), optimized_minprop = 0.10,
    include_skcm_primary_only = TRUE, make_km = TRUE,
    make_summary_plots = TRUE, write_files = TRUE, verbose = TRUE) {
  if (is.character(prepared) && length(prepared) == 1L) prepared <- readRDS(prepared)
  if (!is.list(prepared) || is.null(prepared$expression) || is.null(prepared$cohort_samples)) {
    stop("prepared requires expression and cohort_samples.", call. = FALSE)
  }
  state_signature_genes <- unique(tcgasig_clean_missing(state_signature_genes))
  state_signature_genes <- state_signature_genes[!is.na(state_signature_genes)]
  cell_signature_genes <- unique(tcgasig_clean_missing(cell_signature_genes))
  cell_signature_genes <- cell_signature_genes[!is.na(cell_signature_genes)]
  if (!length(state_signature_genes)) stop("The state gene set is empty.", call. = FALSE)
  limits <- list(min_state_genes = min_state_genes, min_cell_genes = min_cell_genes,
    min_patients = min_patients, min_events = min_events,
    min_subset_patients = min_subset_patients, min_subset_events = min_subset_events)
  for (name in names(limits)) tcgasig_validate_threshold(limits[[name]], name)
  cutoff_methods <- match.arg(cutoff_methods, c("median", "optimized"), several.ok = TRUE)
  if (length(optimized_minprop) != 1L || !is.finite(optimized_minprop) ||
      optimized_minprop <= 0 || optimized_minprop >= 0.5) stop("optimized_minprop must be between 0 and 0.5.", call. = FALSE)
  expr <- data.table::copy(data.table::as.data.table(prepared$expression))
  cs <- data.table::copy(data.table::as.data.table(prepared$cohort_samples))
  required <- c("cohort_label", "cancer_type", "sample_mode", "default_cohort",
                "patient_id", "sample_barcode", "OS.time", "OS.status")
  if (!all(required %in% names(cs)) || !all(c("gene_name", "gene_id") %in% names(expr))) {
    stop("Prepared tables lack required identifiers/columns.", call. = FALSE)
  }
  if (anyNA(cs[, c("cohort_label", "patient_id", "sample_barcode", "default_cohort"), with = FALSE]) ||
      !is.logical(cs$default_cohort) || anyDuplicated(expr$gene_name)) {
    stop("Prepared identifiers must be nonmissing, genes unique, default_cohort logical.", call. = FALSE)
  }
  if (anyDuplicated(cs[, c("cohort_label", "patient_id"), with = FALSE])) {
    stop("Prepared cohorts contain duplicate patients.", call. = FALSE)
  }
  covariates <- unique(covariates)
  if (!all(covariates %in% names(cs)) || any(make.names(covariates) != covariates) ||
      any(covariates %in% c(required, "state_score", "cell_score", "state_score_z", "cell_score_z", "group"))) {
    stop("covariates must name existing, syntactic, nonreserved clinical columns.", call. = FALSE)
  }
  if (!include_skcm_primary_only) cs <- cs[!(cs$cancer_type == "SKCM" & !cs$default_cohort), ]
  if (!is.null(cell_scores)) {
    if (!all(c("sample_barcode", "cell_score") %in% names(cell_scores)) ||
        anyDuplicated(cell_scores$sample_barcode) || anyNA(cell_scores$sample_barcode) ||
        !is.numeric(cell_scores$cell_score)) {
      stop("cell_scores requires unique sample_barcode and numeric cell_score.", call. = FALSE)
    }
    if (identical(cell_score_method, "marker_module") || !nzchar(cell_score_method)) {
      stop("Set cell_score_method to document external cell_scores.", call. = FALSE)
    }
  } else if (!length(cell_signature_genes)) stop("The cell gene set is empty.", call. = FALSE)
  output_prefix <- if (is.null(output_prefix)) paste0("panTCGA_", tcgasig_safe_name(cell_name),
    "_", tcgasig_safe_name(state_name), "_two_layer") else tcgasig_safe_name(output_prefix)
  if (write_files || make_km || make_summary_plots) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  km_dir <- file.path(out_dir, "KM_plots_two_layer")
  if (make_km) dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list(scores = list(), adjusted = list(), continuous = list(), cell_high = list(),
               skipped = list(), gene_qc = list(), cohort_qc = list())
  add <- function(name, row) rows[[name]][[length(rows[[name]]) + 1L]] <<- row
  skip <- function(label, model, method, reason, n, events) {
    add("skipped", data.frame(cohort_label = label, model = model, cutoff_method = method,
      reason = reason, n_patients = n, n_events = events))
  }
  samples <- setdiff(names(expr), c("gene_name", "gene_id", "gene_type"))
  for (label in unique(cs$cohort_label)) {
    dt <- as.data.frame(cs[cs$cohort_label == label & cs$sample_barcode %in% samples, ])
    if (!nrow(dt)) { skip(label, "all", NA_character_, "No overlapping expression samples.", 0, 0); next }
    if (nrow(unique(dt[, c("cancer_type", "sample_mode", "default_cohort")])) != 1L) {
      stop("Inconsistent cohort metadata: ", label, call. = FALSE)
    }
    tcgasig_message(verbose, "Two-layer analysis: ", label)
    state <- tcgasig_score_module(expr, state_signature_genes, dt$sample_barcode, min_state_genes)
    cell <- if (is.null(cell_scores)) tcgasig_score_module(expr, cell_signature_genes,
       dt$sample_barcode, min_cell_genes) else list(
       score = cell_scores$cell_score[match(dt$sample_barcode, cell_scores$sample_barcode)],
       counts = rep(NA_integer_, nrow(dt)), used = character(), missing = character())
    for (layer in c("state", "cell")) {
      module <- if (layer == "state") state else cell
      requested <- if (layer == "state") state_signature_genes else if (is.null(cell_scores)) cell_signature_genes else character()
      if (length(requested)) add("gene_qc", data.frame(cohort_label = label, layer = layer,
        gene = requested, usable = requested %in% module$used))
    }
    dt$state_score <- state$score
    dt$cell_score <- cell$score
    dt$state_genes_n <- state$counts
    dt$cell_genes_n <- cell$counts
    valid <- is.finite(dt$OS.time) & dt$OS.time > 0 & dt$OS.status %in% c(0, 1) &
      is.finite(dt$state_score) & is.finite(dt$cell_score)
    if (length(covariates)) {
      valid <- valid & stats::complete.cases(dt[, covariates, drop = FALSE])
      for (column in covariates) if (is.numeric(dt[[column]])) valid <- valid & is.finite(dt[[column]])
    }
    input_n <- nrow(dt)
    dt <- dt[valid, , drop = FALSE]
    n <- nrow(dt)
    events <- sum(dt$OS.status == 1)
    dt$state_score_z <- tcgasig_zscore_or_na(dt$state_score)
    dt$cell_score_z <- tcgasig_zscore_or_na(dt$cell_score)
    cell_cut <- if (n) stats::median(dt$cell_score) else NA_real_
    dt$cell_group <- ifelse(dt$cell_score >= cell_cut, "high", "low")
    add("scores", dt)
    correlation <- if (n > 2 && all(is.finite(dt$state_score_z)) && all(is.finite(dt$cell_score_z)))
      stats::cor(dt$state_score, dt$cell_score, method = "spearman") else NA_real_
    add("cohort_qc", data.frame(cohort_label = label, input_n = input_n, analyzed_n = n,
      excluded_n = input_n - n, n_events = events, cell_cutoff = cell_cut,
      state_cell_spearman = correlation, cell_score_method = cell_score_method,
      overlapping_genes = paste(intersect(cell_signature_genes, state_signature_genes), collapse = ",")))
    if (n < min_patients || events < min_events) {
      skip(label, "all", NA_character_, "Below full-cohort patient/event thresholds.", n, events); next
    }
    if (any(!is.finite(dt$state_score_z)) || any(!is.finite(dt$cell_score_z))) {
      skip(label, "all", NA_character_, "Constant state or cell score.", n, events); next
    }
    meta <- dt[1, c("cohort_label", "cancer_type", "sample_mode", "default_cohort"), drop = FALSE]
    extra <- if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""
    fitted <- tcgasig_fit_checked(dt, paste0("state_score_z + cell_score_z", extra), "state_score_z")
    if (is.null(fitted$error)) {
      cell_s <- summary(fitted$fit)
      add("continuous", cbind(meta, n_patients = n, n_events = events, fitted$row,
        cell_HR_per_1sd = cell_s$conf.int["cell_score_z", "exp(coef)"],
        cell_wald_p = cell_s$coefficients["cell_score_z", "Pr(>|z|)"]))
    } else skip(label, "continuous", NA_character_, fitted$error, n, events)
    for (method in cutoff_methods) {
      for (model in c("adjusted", "cell_high")) {
        model_dt <- if (model == "adjusted") dt else dt[dt$cell_group == "high", , drop = FALSE]
        model_n <- nrow(model_dt)
        model_events <- sum(model_dt$OS.status == 1)
        if (model == "cell_high" && (model_n < min_subset_patients || model_events < min_subset_events)) {
          skip(label, model, method, "Below cell-high patient/event thresholds.", model_n, model_events); next
        }
        model_dt$target_signature_score <- model_dt$state_score
        grouped <- tcgasig_make_groups(data.table::as.data.table(model_dt), method, optimized_minprop)
        if (!is.na(grouped$error)) { skip(label, model, method, grouped$error, model_n, model_events); next }
        model_dt <- as.data.frame(grouped$data)
        model_dt$group <- factor(model_dt$group, levels = c("low", "high"))
        if (any(table(model_dt$group) < 2L)) {
          skip(label, model, method, "Fewer than two patients in a state group.", model_n, model_events); next
        }
        rhs <- paste0(if (model == "adjusted") "group + cell_score_z" else "group", extra)
        fitted <- tcgasig_fit_checked(model_dt, rhs, "grouphigh")
        if (!is.null(fitted$error)) { skip(label, model, method, fitted$error, model_n, model_events); next }
        row <- cbind(meta, cutoff_method = method, cutoff_value = grouped$cutoff,
          cell_cutoff = cell_cut, n_patients = model_n, n_events = model_events,
          high_n = sum(model_dt$group == "high"), low_n = sum(model_dt$group == "low"), fitted$row)
        add(model, row)
        if (make_km && model == "cell_high") {
          tcgasig_plot_signature_only_km(data.table::as.data.table(model_dt), label,
            paste0(cell_name, "-high"), method, grouped$cutoff, row,
            file.path(km_dir, paste0(output_prefix, "_", tcgasig_safe_name(label), "_", method, ".pdf")),
            state_name)
        }
      }
    }
  }
  meta_empty <- data.table::data.table(cohort_label = character(), cancer_type = character(),
    sample_mode = character(), default_cohort = logical())
  model_empty <- cbind(meta_empty, data.table::data.table(n_patients = integer(), n_events = integer(),
    HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), wald_p = numeric(), PH_global_p = numeric()))
  group_empty <- cbind(model_empty, data.table::data.table(cutoff_method = character(),
    cutoff_value = numeric(), cell_cutoff = numeric(), high_n = integer(), low_n = integer()))
  score_empty <- data.table::copy(cs[0L])
  for (column in c("state_score", "cell_score", "state_genes_n", "cell_genes_n", "state_score_z", "cell_score_z")) score_empty[[column]] <- numeric()
  score_empty$cell_group <- character()
  templates <- list(scores = score_empty, adjusted = group_empty, continuous = model_empty,
    cell_high = group_empty, skipped = data.table::data.table(cohort_label = character(), model = character(),
      cutoff_method = character(), reason = character(), n_patients = integer(), n_events = integer()),
    gene_qc = data.table::data.table(cohort_label = character(), layer = character(), gene = character(), usable = logical()),
    cohort_qc = data.table::data.table(cohort_label = character(), input_n = integer(), analyzed_n = integer(),
      excluded_n = integer(), n_events = integer(), cell_cutoff = numeric(), state_cell_spearman = numeric(),
      cell_score_method = character(), overlapping_genes = character()))
  result <- lapply(names(rows), function(name) if (length(rows[[name]]))
    data.table::rbindlist(rows[[name]], fill = TRUE) else data.table::copy(templates[[name]]))
  names(result) <- names(rows)
  for (name in c("adjusted", "continuous", "cell_high")) {
    dt <- result[[name]]
    if (!nrow(dt)) next
    families <- if (name == "continuous") "default_cohort" else c("cutoff_method", "default_cohort")
    dt[, FDR := stats::p.adjust(wald_p, method = "BH"), by = families]
    dt[, direction := ifelse(HR > 1, "state_high_worse", "state_high_better")]
    if (name == "continuous") dt[, cell_FDR := stats::p.adjust(cell_wald_p, method = "BH"), by = families]
    data.table::setorderv(dt, c("FDR", "wald_p"))
    result[[name]] <- dt
  }
  result$settings <- c(list(cell_name = cell_name, state_name = state_name,
    cell_signature_genes = cell_signature_genes, state_signature_genes = state_signature_genes,
    cell_score_method = cell_score_method, covariates = covariates,
    cutoff_methods = cutoff_methods, optimized_minprop = optimized_minprop,
    FDR_families = "model x cutoff_method x default_cohort"), limits)
  result$files <- list()
  if (write_files) {
    for (name in names(rows)) {
      path <- file.path(out_dir, paste0(output_prefix, "_", name, ".tsv"))
      data.table::fwrite(result[[name]], path, sep = "\t")
      result$files[[name]] <- path
    }
    result$files$qc <- file.path(out_dir, paste0(output_prefix, "_QC_report.txt"))
    writeLines(c("TCGASigSurvival two-layer analysis", paste("Cell context:", cell_name),
      paste("State:", state_name), paste("Cell score method:", cell_score_method),
      "Scores are cohort-relative; marker scores do not identify individual cells.",
      "Median and continuous models are preferred; optimized cutoffs are exploratory.",
      "BH families: model x cutoff_method x default_cohort.",
      "Convergence failures/nonfinite estimates are recorded in the skipped table.",
      paste("Additional covariates:", paste(covariates, collapse = ",")),
      utils::capture.output(utils::str(result$settings))), result$files$qc)
    result$files$session <- file.path(out_dir, paste0(output_prefix, "_sessionInfo.txt"))
    writeLines(utils::capture.output(utils::sessionInfo()), result$files$session)
  }
  if (make_summary_plots) {
    result$plots <- plot_two_layer_results(result, out_dir = out_dir, output_prefix = output_prefix)
  }
  result
}

#' Prepare public TCGA inputs and run two-layer survival analysis
#' @inheritParams run_two_layer_survival
#' @param project_dir Project directory.
#' @param data_dir Directory containing the four public TCGA inputs.
#' @param ... Further arguments to run_two_layer_survival, excluding prepared.
#' @return The two-layer result list with an additional prepared element.
#' @export
run_pan_tcga_two_layer_survival <- function(state_signature_genes,
    cell_signature_genes = default_naive_cd4_signature(), state_name = "State",
    cell_name = "NaiveCD4", project_dir = getwd(),
    data_dir = file.path(project_dir, "data", "xena"),
    out_dir = file.path(project_dir, "results"), output_prefix = NULL, ...) {
  prefix <- if (is.null(output_prefix)) paste0("panTCGA_", tcgasig_safe_name(cell_name),
    "_", tcgasig_safe_name(state_name), "_two_layer") else tcgasig_safe_name(output_prefix)
  prepared <- prepare_signature_only_data(unique(c(cell_signature_genes, state_signature_genes)),
    signature_name = paste(cell_name, state_name, sep = "_"), project_dir = project_dir,
    data_dir = data_dir, out_dir = out_dir, output_prefix = paste0(prefix, "_input"),
    write_files = FALSE)
  result <- run_two_layer_survival(prepared, state_signature_genes, cell_signature_genes,
    state_name, cell_name, out_dir = out_dir, output_prefix = prefix, ...)
  result$prepared <- prepared
  result
}

#' Plot two-layer group-comparison forest plots
#' @param result Result from run_two_layer_survival.
#' @param out_dir Output directory.
#' @param output_prefix Output basename.
#' @return A list of PDF and PNG paths.
#' @export
plot_two_layer_results <- function(result, out_dir = file.path(getwd(), "results"),
    output_prefix = "two_layer") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list()
  for (model in c("adjusted", "cell_high")) {
    dt <- as.data.frame(result[[model]])
    if (!nrow(dt)) next
    dt <- dt[dt$default_cohort & is.finite(dt$HR) & dt$CI_lower > 0, , drop = FALSE]
    for (method in unique(dt$cutoff_method)) {
      d <- dt[dt$cutoff_method == method, , drop = FALSE]
      d$cancer_plot <- factor(d$cancer_type, levels = d$cancer_type[order(d$HR)])
      d$significant <- d$FDR < 0.05
      p <- ggplot2::ggplot(d, ggplot2::aes(x = cancer_plot, y = HR, ymin = CI_lower, ymax = CI_upper, colour = significant)) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_pointrange() + ggplot2::coord_flip() + ggplot2::scale_y_log10() +
        ggplot2::scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "#B8323C"), guide = "none") +
        ggplot2::labs(x = NULL, y = "HR: state high vs low", title = paste(model, method, sep = " / ")) +
        ggplot2::theme_bw(base_size = 10)
      prefix <- file.path(out_dir, paste0(tcgasig_safe_name(output_prefix), "_", model, "_forest_", method))
      height <- max(4, nrow(d) * 0.23 + 1.5)
      ggplot2::ggsave(paste0(prefix, ".pdf"), p, width = 7, height = height)
      ggplot2::ggsave(paste0(prefix, ".png"), p, width = 7, height = height, dpi = 300)
      paths[[paste(model, method, sep = "_")]] <- c(pdf = paste0(prefix, ".pdf"), png = paste0(prefix, ".png"))
    }
  }
  paths
}
