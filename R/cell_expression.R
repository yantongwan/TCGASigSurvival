tcgasig_need <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install the optional package '", package, "' for this operation.", call. = FALSE)
  }
}

tcgasig_cell_check <- function(data, allow_proxy = FALSE) {
  if (!inherits(data, "tcgasig_cell_data")) stop("Use prepare_cell_expression_data first.", call. = FALSE)
  if (identical(data$measurement, "bulk_proxy") && !isTRUE(allow_proxy)) {
    stop("bulk_proxy is not cell-specific expression; set allow_proxy = TRUE explicitly.", call. = FALSE)
  }
  invisible(data)
}

tcgasig_cell_label <- function(data) {
  switch(data$measurement, bulk_proxy = "Bulk RNA with a cell-context proxy",
    deconvolved = "Estimated cell-specific expression", sorted = "Sorted-cell expression",
    single_cell_pseudobulk = "Patient-level cell-type pseudobulk")
}

tcgasig_cell_log <- function(data) {
  if (data$expression_scale == "counts") {
    tcgasig_need("edgeR")
    x <- data$expression
    keys <- interaction(data$metadata$cohort_label, data$metadata$cell_type, drop = TRUE)
    for (ii in split(seq_len(ncol(x)), keys)) {
      y <- edgeR::calcNormFactors(edgeR::DGEList(data$expression[, ii, drop = FALSE]))
      x[, ii] <- edgeR::cpm(y, log = TRUE, prior.count = 0.5)
    }
    return(x)
  }
  data$expression
}

tcgasig_write_cell_tables <- function(x, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- character()
  for (name in names(x)) if (is.data.frame(x[[name]])) {
    path <- file.path(out_dir, paste0(tcgasig_safe_name(prefix), "_", name, ".tsv"))
    data.table::fwrite(x[[name]], path, sep = "\t", na = "NA")
    files[name] <- path
  }
  files
}

#' Prepare patient-level cell-type expression data
#'
#' One expression column is one biological donor, cell type and tissue status,
#' not one cell. Clinical columns are retained without filtering on OS, so
#' missing OS does not exclude a patient from PFI, PFS, DSS or DFI analyses.
#' Deconvolved expression is an estimate, not a measured cell-specific value.
#' Bulk data require explicit opt-in at each analysis step.
#' @param expression Named numeric gene-by-profile matrix. Counts must be raw
#'   nonnegative integer counts; log values must already be normalized.
#' @param metadata Data frame with profile_id, patient_id, cancer_type,
#'   cell_type and tissue_status (tumor, healthy or adjacent_normal).
#'   Optional cohort_label defaults to cancer_type. Clinical outcome columns
#'   use endpoint.time and endpoint.status, in days and 0/1 respectively.
#' @param measurement sorted, single_cell_pseudobulk, deconvolved or bulk_proxy.
#' @param expression_scale log, tpm (converted to log2(TPM+1)) or counts.
#' @param provenance Required source/annotation/estimation description.
#' @return A tcgasig_cell_data list. Input rows are never silently aggregated.
#' @export
prepare_cell_expression_data <- function(expression, metadata,
    measurement = c("sorted", "single_cell_pseudobulk", "deconvolved", "bulk_proxy"),
    expression_scale = c("log", "tpm", "counts"), provenance) {
  measurement <- match.arg(measurement)
  expression_scale <- match.arg(expression_scale)
  if (missing(provenance) || length(provenance) != 1L || is.na(provenance) || !nzchar(provenance)) {
    stop("Document provenance, cell annotation and/or estimation method.", call. = FALSE)
  }
  expression <- as.matrix(expression)
  if (!is.numeric(expression) || !nrow(expression) || !ncol(expression) ||
      is.null(rownames(expression)) || is.null(colnames(expression)) ||
      anyNA(expression) || any(!is.finite(expression))) {
    stop("expression must be a finite named numeric matrix.", call. = FALSE)
  }
  for (ids in dimnames(expression)) if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Expression gene/profile identifiers must be nonempty and unique.", call. = FALSE)
  }
  if (expression_scale != "log" && any(expression < 0)) stop("TPM/counts cannot be negative.", call. = FALSE)
  if (expression_scale == "counts" && any(abs(expression - round(expression)) > 1e-7)) {
    stop("Counts must be raw integers, not TPM or normalized expression.", call. = FALSE)
  }
  if (expression_scale == "counts" && any(colSums(expression) == 0)) stop("Empty count libraries.", call. = FALSE)
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  required <- c("profile_id", "patient_id", "cancer_type", "cell_type", "tissue_status")
  if (!all(required %in% names(meta))) stop("metadata requires: ", paste(required, collapse = ", "), call. = FALSE)
  for (name in required) {
    meta[[name]] <- tcgasig_clean_missing(as.character(meta[[name]]))
    if (anyNA(meta[[name]])) stop("Missing identifier/status: ", name, call. = FALSE)
  }
  if (anyDuplicated(meta$profile_id) || !setequal(meta$profile_id, colnames(expression))) {
    stop("Metadata and expression must match exactly by unique profile_id.", call. = FALSE)
  }
  meta <- meta[match(colnames(expression), meta$profile_id), , drop = FALSE]
  if (!all(meta$tissue_status %in% c("tumor", "healthy", "adjacent_normal"))) {
    stop("Distinguish tumor, healthy donor and adjacent_normal explicitly.", call. = FALSE)
  }
  if (!"cohort_label" %in% names(meta)) meta$cohort_label <- meta$cancer_type
  if (anyNA(meta$cohort_label) || any(!nzchar(meta$cohort_label))) stop("Invalid cohort_label.", call. = FALSE)
  by_cohort <- split(meta$cancer_type, meta$cohort_label)
  if (any(vapply(by_cohort, function(x) length(unique(x)) != 1L, logical(1)))) {
    stop("Each cohort_label must identify one cancer type.", call. = FALSE)
  }
  key <- meta[, c("cohort_label", "patient_id", "cell_type", "tissue_status")]
  if (anyDuplicated(key)) stop("Multiple profiles per donor/state: aggregate technical samples explicitly.", call. = FALSE)
  if (expression_scale == "tpm") expression <- log2(expression + 1)
  structure(list(expression = expression, metadata = meta, measurement = measurement,
    expression_scale = if (expression_scale == "tpm") "log" else expression_scale,
    provenance = provenance), class = "tcgasig_cell_data")
}

#' Aggregate annotated single cells into donor-level pseudobulk
#' @param counts Named raw gene-by-cell counts, dense or a Matrix object.
#' @param cells Cell metadata with cell_id, patient_id, cancer_type, cell_type,
#'   tissue_status and optional clinical/covariate columns. Retained clinical
#'   metadata must be identical across a donor's cells.
#' @param cell_type Exact annotated cell-type label to retain.
#' @param min_cells Minimum cells per donor/tissue profile; exclusions are recorded.
#' @param provenance Description of study, annotation, count layer and QC.
#' @return Prepared count data with cell_qc table, or an error if all fail QC.
#' @export
aggregate_cell_pseudobulk <- function(counts, cells, cell_type = "Treg",
    min_cells = 20L, provenance) {
  tcgasig_validate_threshold(min_cells, "min_cells")
  required <- c("cell_id", "patient_id", "cancer_type", "cell_type", "tissue_status")
  cells <- as.data.frame(cells, stringsAsFactors = FALSE)
  if (!all(required %in% names(cells)) || anyNA(cells[, required]) ||
      anyDuplicated(cells$cell_id) || !setequal(cells$cell_id, colnames(counts))) {
    stop("Cells must match raw count columns exactly, with unique IDs and complete donor metadata.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) stop("Unique gene symbols required.", call. = FALSE)
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0) || any(abs(values - round(values)) > 1e-7)) stop("Raw integer counts required.", call. = FALSE)
  if (!"cohort_label" %in% names(cells)) cells$cohort_label <- cells$cancer_type
  selected <- cells[cells$cell_type == cell_type, , drop = FALSE]
  if (!nrow(selected)) stop("No annotated cells of type ", cell_type, call. = FALSE)
  keys <- interaction(selected$cohort_label, selected$patient_id, selected$tissue_status, drop = TRUE)
  groups <- split(seq_len(nrow(selected)), keys)
  meta_rows <- count_rows <- qc_rows <- list()
  for (i in seq_along(groups)) {
    d <- selected[groups[[i]], , drop = FALSE]
    cols <- setdiff(names(d), c("cell_id", "profile_id"))
    if (any(vapply(d[, cols, drop = FALSE], function(v) length(unique(v)) != 1L, logical(1)))) {
      stop("Conflicting donor clinical metadata in pseudobulk group.", call. = FALSE)
    }
    id <- sprintf("PB_%05d", i)
    m <- d[1L, cols, drop = FALSE]
    m$profile_id <- id
    m$n_cells <- nrow(d)
    qc_rows[[i]] <- cbind(m, retained = nrow(d) >= min_cells)
    if (nrow(d) < min_cells) next
    x <- counts[, match(d$cell_id, colnames(counts)), drop = FALSE]
    if (inherits(x, "Matrix")) {
      tcgasig_need("Matrix")
      total <- Matrix::rowSums(x)
    } else total <- rowSums(x)
    count_rows[[id]] <- total
    meta_rows[[id]] <- m
  }
  if (!length(count_rows)) stop("No donor profiles pass min_cells.", call. = FALSE)
  expr <- do.call(cbind, count_rows)
  rownames(expr) <- rownames(counts)
  result <- prepare_cell_expression_data(expr, data.table::rbindlist(meta_rows, fill = TRUE),
    "single_cell_pseudobulk", "counts", provenance)
  result$cell_qc <- as.data.frame(data.table::rbindlist(qc_rows, fill = TRUE))
  result
}

#' Compare a cell-type gene between tumor and reference donors
#' @param data Prepared cell expression data.
#' @param target_gene Exact gene symbol, e.g. human ACLY (mouse Acly is not mapped).
#' @param cell_type Exact cell label.
#' @param reference healthy or adjacent_normal. Never treated interchangeably.
#' @param paired Whether to use only tumor/reference pairs from the same donor.
#' @param min_donors Minimum donors in each group or minimum complete pairs.
#' @param allow_proxy Explicit permission to analyze bulk proxies.
#' @param out_dir,prefix Output directory and basename.
#' @param write_files Write source tables.
#' @param score_mode gene (normalized expression) or joint (target z plus
#'   cell-score z, compatible with the legacy target/signature method).
#' @return List with per-donor values, per-cohort tests, measurement and files.
#'   Tests are two-sided Wilcoxon tests; effect is tumor minus reference median
#'   (median within-donor difference for paired tests), on normalized log scale.
#' @export
compare_cell_gene_expression <- function(data, target_gene = "ACLY", cell_type = "Treg",
    reference = c("healthy", "adjacent_normal"), paired = FALSE, min_donors = 3L,
    allow_proxy = FALSE, out_dir = file.path(getwd(), "results"), prefix = "cell_gene",
    write_files = TRUE, score_mode = c("gene", "joint")) {
  tcgasig_cell_check(data, allow_proxy)
  reference <- match.arg(reference)
  score_mode <- match.arg(score_mode)
  tcgasig_validate_threshold(min_donors, "min_donors")
  if (!target_gene %in% rownames(data$expression)) stop("Target gene not measured.", call. = FALSE)
  meta <- data$metadata[data$metadata$cell_type == cell_type, , drop = FALSE]
  if (!nrow(meta)) stop("Cell type not present.", call. = FALSE)
  meta <- tcgasig_cell_scores(data, target_gene, meta, score_mode)
  meta$measurement <- data$measurement
  meta$used_in_test <- FALSE
  rows <- list()
  for (cohort in unique(meta$cohort_label)) {
    ii <- which(meta$cohort_label == cohort & meta$tissue_status %in% c("tumor", reference))
    d <- meta[ii, , drop = FALSE]
    a <- d[d$tissue_status == "tumor" & is.finite(d$analysis_score), , drop = FALSE]
    b <- d[d$tissue_status == reference & is.finite(d$analysis_score), , drop = FALSE]
    if (!paired && length(intersect(a$patient_id, b$patient_id))) {
      stop("Shared tumor/reference donors require paired = TRUE; do not treat them as independent.", call. = FALSE)
    }
    if (paired) {
      ids <- intersect(a$patient_id, b$patient_id)
      a <- a[match(ids, a$patient_id), , drop = FALSE]
      b <- b[match(ids, b$patient_id), , drop = FALSE]
    }
    n_a <- nrow(a); n_b <- nrow(b)
    status <- if (!n_b) "reference_not_available" else if (min(n_a, n_b) < min_donors) "insufficient_donors" else "ok"
    effect <- p <- NA_real_
    if (status == "ok") {
      effect <- if (paired) stats::median(a$analysis_score - b$analysis_score) else
        stats::median(a$analysis_score) - stats::median(b$analysis_score)
      p <- suppressWarnings(stats::wilcox.test(a$analysis_score, b$analysis_score,
        paired = paired, exact = FALSE)$p.value)
      if (!is.finite(p)) { p <- NA_real_; status <- "not_estimable" }
      meta$used_in_test[meta$profile_id %in% c(a$profile_id, b$profile_id)] <- TRUE
    }
    rows[[cohort]] <- data.frame(cohort_label = cohort, cancer_type = d$cancer_type[1],
      target_gene = target_gene, cell_type = cell_type, measurement = data$measurement, score_mode = score_mode,
      reference = reference, paired = paired, tumor_n = n_a, reference_n = n_b,
      tumor_minus_reference = effect, p_value = p, status = status)
  }
  tests <- as.data.frame(data.table::rbindlist(rows, fill = TRUE))
  tests$FDR <- stats::p.adjust(tests$p_value, "BH")
  result <- list(values = meta, tests = tests, measurement = data$measurement,
    interpretation = tcgasig_cell_label(data), target_gene = target_gene, reference = reference,
    score_mode = score_mode)
  result$files <- if (write_files) tcgasig_write_cell_tables(result, out_dir, prefix) else character()
  result
}

tcgasig_cell_scores <- function(data, target_gene, meta, score_mode) {
  x <- tcgasig_cell_log(data)
  meta$gene_expression <- as.numeric(x[target_gene, match(meta$profile_id, colnames(x))])
  meta$analysis_score <- meta$gene_expression
  if (score_mode == "joint") {
    if (!"cell_score" %in% names(meta) || !is.numeric(meta$cell_score)) {
      stop("joint mode requires a documented numeric cell_score.", call. = FALSE)
    }
    meta$analysis_score <- NA_real_
    for (cohort in unique(meta$cohort_label)) {
      ii <- which(meta$cohort_label == cohort)
      tt <- ii[meta$tissue_status[ii] == "tumor"]
      # Apply tumor reference moments to normal specimens; never standardize controls separately.
      gene_sd <- stats::sd(meta$gene_expression[tt], na.rm = TRUE)
      cell_sd <- stats::sd(meta$cell_score[tt], na.rm = TRUE)
      if (is.finite(gene_sd) && gene_sd > 0 && is.finite(cell_sd) && cell_sd > 0) {
        meta$analysis_score[ii] <- (meta$gene_expression[ii] - mean(meta$gene_expression[tt], na.rm = TRUE)) / gene_sd +
          (meta$cell_score[ii] - mean(meta$cell_score[tt], na.rm = TRUE)) / cell_sd
      }
    }
  }
  meta
}

tcgasig_cell_groups <- function(data, target_gene, cell_type, context_subset, score_mode = "gene") {
  if (!target_gene %in% rownames(data$expression)) stop("Target gene not measured.", call. = FALSE)
  meta <- data$metadata[data$metadata$cell_type == cell_type & data$metadata$tissue_status == "tumor", , drop = FALSE]
  if (!nrow(meta)) stop("No tumor profiles for cell type.", call. = FALSE)
  meta <- tcgasig_cell_scores(data, target_gene, meta, score_mode)
  meta$context_selected <- TRUE
  meta$context_cutoff <- meta$gene_cutoff <- NA_real_
  meta$score_cutoff <- NA_real_
  meta$group <- NA_character_
  for (cohort in unique(meta$cohort_label)) {
    ii <- which(meta$cohort_label == cohort)
    if (context_subset == "high") {
      if (!"cell_score" %in% names(meta)) stop("context_subset high requires numeric cell_score.", call. = FALSE)
      score <- meta$cell_score[ii]
      cut <- if (any(is.finite(score))) stats::median(score[is.finite(score)]) else NA_real_
      meta$context_cutoff[ii] <- cut
      meta$context_selected[ii] <- is.finite(score) & score >= cut
    }
    jj <- ii[meta$context_selected[ii] & is.finite(meta$analysis_score[ii])]
    if (!length(jj)) next
    cut <- stats::median(meta$analysis_score[jj])
    meta$score_cutoff[ii] <- cut
    if (score_mode == "gene") meta$gene_cutoff[ii] <- cut
    meta$group[jj] <- ifelse(meta$analysis_score[jj] >= cut, "high", "low")
  }
  meta$group <- factor(meta$group, levels = c("low", "high"))
  meta
}

tcgasig_covariate_rows <- function(d, covariates) {
  if (!all(covariates %in% names(d)) || any(make.names(covariates) != covariates) ||
      any(covariates %in% c("group", "gene_expression", "gene_z", "endpoint_time", "endpoint_status"))) {
    stop("Use existing, syntactic, nonreserved covariate names.", call. = FALSE)
  }
  keep <- if (length(covariates)) stats::complete.cases(d[, covariates, drop = FALSE]) else rep(TRUE, nrow(d))
  for (v in covariates) if (is.numeric(d[[v]])) keep <- keep & is.finite(d[[v]])
  d <- d[keep, , drop = FALSE]
  active <- covariates[vapply(covariates, function(v) length(unique(d[[v]])) > 1L, logical(1))]
  for (v in active) if (is.character(d[[v]]) || is.factor(d[[v]])) d[[v]] <- droplevels(factor(d[[v]]))
  list(data = d, active = active, constant = setdiff(covariates, active))
}
