#' Exploratory marker panels for five immune-cell contexts
#'
#' Treg and CD8 retain the package's historical panels. The other panels are
#' pragmatic human marker modules, not validated cell counters. CD8 includes
#' cytotoxic/NK-associated genes; CD4 includes shared T-cell genes; myeloid
#' modules overlap with inflammatory and polarization states. Users may supply
#' independently justified panels. No score is a cell percentage.
#' @return Named list of human gene-symbol vectors.
#' @export
default_immune_signatures <- function() {
  list(CD8 = default_cd8_signature(), Treg = default_treg_signature(),
    Neutrophils = c("CEACAM8", "FCGR3B", "CSF3R", "CXCR2", "FPR1", "S100A8", "S100A9"),
    Macrophages = c("C1QA", "C1QB", "C1QC", "CSF1R", "CD68", "CD163", "MSR1", "MRC1"),
    CD4 = c("CD4", "CD3D", "CD3E", "TRAC", "CD2", "IL7R", "LTB"))
}

tcgasig_abundance_matrix <- function(expression) {
  x <- as.matrix(expression)
  if (!is.numeric(x) || !nrow(x) || !ncol(x) || is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) || anyDuplicated(rownames(x)) ||
      anyDuplicated(colnames(x)) || any(!nzchar(rownames(x))) || any(!nzchar(colnames(x)))) {
    stop("Expression must be a numeric matrix with unique nonempty gene and sample names.", call. = FALSE)
  }
  x
}

tcgasig_abundance_clinical <- function(clinical, samples = NULL) {
  d <- as.data.frame(clinical, stringsAsFactors = FALSE)
  required <- c("sample_barcode", "patient_id", "cohort_label")
  if (!all(required %in% names(d)) || anyNA(d[, required]) ||
      any(vapply(d[, required, drop = FALSE], function(v) any(!nzchar(as.character(v))), logical(1)))) {
    stop("Clinical data require nonmissing sample_barcode, patient_id and cohort_label.", call. = FALSE)
  }
  if (anyDuplicated(d[, c("cohort_label", "patient_id")]) ||
      anyDuplicated(d[, c("cohort_label", "sample_barcode")])) {
    stop("Use one specimen per independent patient within each cohort.", call. = FALSE)
  }
  if (!is.null(samples) && !all(d$sample_barcode %in% samples)) {
    stop("Clinical sample IDs must match expression exactly; unmatched IDs are not silently dropped.", call. = FALSE)
  }
  d
}

#' Compute marker scores without a target gene
#' @param expression Named gene-by-sample matrix.
#' @param clinical Patient/sample/cohort table, one sample per patient per cohort.
#' @param signatures Named gene lists; default exploratory panels are not fractions.
#' @param expression_scale log or tpm; tpm is transformed with log2(TPM+1).
#' @param min_coverage Minimum fraction of requested markers available per sample.
#' @return List with long-form values, gene_qc, sample_qc and overlap audit.
#' @export
score_immune_signatures <- function(expression, clinical,
    signatures = default_immune_signatures(), expression_scale = c("log", "tpm"),
    min_coverage = 0.7) {
  x <- tcgasig_abundance_matrix(expression)
  d <- tcgasig_abundance_clinical(clinical, colnames(x))
  expression_scale <- match.arg(expression_scale)
  if (expression_scale == "tpm") {
    if (any(x < 0, na.rm = TRUE)) stop("TPM must be nonnegative.", call. = FALSE)
    x <- log2(x + 1)
  }
  if (!is.list(signatures) || !length(signatures) || is.null(names(signatures)) ||
      anyNA(names(signatures)) || any(!nzchar(names(signatures))) || anyDuplicated(names(signatures))) {
    stop("Use uniquely named nonempty signature lists.", call. = FALSE)
  }
  if (length(min_coverage) != 1L || !is.finite(min_coverage) || min_coverage <= 0 || min_coverage > 1) {
    stop("min_coverage must be in (0,1].", call. = FALSE)
  }
  signatures <- lapply(signatures, function(g) {
    g <- unique(tcgasig_clean_missing(as.character(g)))
    g <- g[!is.na(g)]
    if (!length(g)) stop("Empty signature.", call. = FALSE)
    g
  })
  values <- gene_qc <- sample_qc <- list()
  for (cohort in unique(d$cohort_label)) {
    m <- d[d$cohort_label == cohort, , drop = FALSE]
    for (cell in names(signatures)) {
      genes <- signatures[[cell]]; minimum <- ceiling(length(genes) * min_coverage)
      z <- list(); reasons <- rep("absent", length(genes)); names(reasons) <- genes
      for (gene in intersect(genes, rownames(x))) {
        v <- x[gene, m$sample_barcode]; v[!is.finite(v)] <- NA_real_
        s <- stats::sd(v, na.rm = TRUE)
        if (is.finite(s) && s > 0) {
          z[[gene]] <- (v - mean(v, na.rm = TRUE)) / s
          reasons[gene] <- "usable"
        } else reasons[gene] <- "constant_or_insufficient_observations"
      }
      score <- rep(NA_real_, nrow(m)); counts <- integer(nrow(m))
      if (length(z)) {
        zz <- do.call(cbind, z); counts <- rowSums(is.finite(zz)); score <- rowMeans(zz, na.rm = TRUE)
        score[counts < minimum | !is.finite(score)] <- NA_real_
      }
      key <- paste(cohort, cell)
      values[[key]] <- data.frame(sample_barcode = m$sample_barcode, patient_id = m$patient_id,
        cohort_label = cohort, cell_type = cell, value = score, measurement = "marker_score",
        method = "mean_gene_z", unit = "mean_gene_z", denominator = "not_a_fraction", role = "primary")
      gene_qc[[key]] <- data.frame(cohort_label = cohort, cell_type = cell, gene = genes,
        status = unname(reasons), requested_genes = length(genes), minimum_genes = minimum)
      sample_qc[[key]] <- data.frame(cohort_label = cohort, cell_type = cell,
        sample_barcode = m$sample_barcode, genes_used = counts, genes_requested = length(genes),
        status = ifelse(is.finite(score), "ok", "insufficient_marker_coverage"))
    }
  }
  pairs <- expand.grid(cell_a = names(signatures), cell_b = names(signatures), stringsAsFactors = FALSE)
  pairs$shared_genes <- vapply(seq_len(nrow(pairs)), function(i)
    paste(intersect(signatures[[pairs$cell_a[i]]], signatures[[pairs$cell_b[i]]]), collapse = ";"), character(1))
  list(values = as.data.frame(data.table::rbindlist(values)),
    gene_qc = as.data.frame(data.table::rbindlist(gene_qc)),
    sample_qc = as.data.frame(data.table::rbindlist(sample_qc)), overlap = pairs, signatures = signatures,
    interpretation = "Exploratory bulk marker/state signal, not a cell fraction or cell-specific expression.")
}

#' Estimate immune fractions with the external quanTIseq implementation
#'
#' Requires full-transcriptome linear TPM, not log expression or counts. The
#' TIL10 reference is supplied by quantiseqr, not redistributed here. CD4 is
#' total CD4 (non-regulatory CD4 plus Treg); CD4_nonTreg and model-defined M1/M2
#' are auxiliary outputs. Other cells are retained in the denominator.
#' @param tpm Named human gene-symbol by sample matrix of linear TPM.
#' @param clinical Patient/sample/cohort table.
#' @param min_signature_coverage Minimum effective TIL10 gene coverage.
#' @return Long-form values, full fractions, signature/sample QC and settings.
#' @export
estimate_immune_fractions <- function(tpm, clinical, min_signature_coverage = 0.8) {
  if (!requireNamespace("quantiseqr", quietly = TRUE)) {
    stop("Install optional dependency with BiocManager::install('quantiseqr').", call. = FALSE)
  }
  x <- tcgasig_abundance_matrix(tpm)
  d <- tcgasig_abundance_clinical(clinical, colnames(x))
  if (any(!is.finite(x)) || any(x < 0) || any(colSums(x) <= 0) || nrow(x) < 1000L) {
    stop("quanTIseq requires finite nonnegative full-transcriptome TPM (at least 1000 genes).", call. = FALSE)
  }
  if (max(x) < 50) stop("TPM range is unsafe for quantiseqr's automatic log detection; do not pass marker-only/log data.", call. = FALSE)
  if (length(min_signature_coverage) != 1L || !is.finite(min_signature_coverage) ||
      min_signature_coverage <= 0 || min_signature_coverage > 1) stop("Invalid signature coverage threshold.", call. = FALSE)
  ref_path <- system.file("extdata/TIL10_signature.txt", package = "quantiseqr", mustWork = TRUE)
  ref <- utils::read.delim(ref_path, row.names = 1, check.names = FALSE)
  noisy <- scan(system.file("extdata/TIL10_rmgenes.txt", package = "quantiseqr"), what = character(), quiet = TRUE)
  aberrant <- scan(system.file("extdata/TIL10_TCGA_aberrant_immune_genes.txt", package = "quantiseqr"), what = character(), quiet = TRUE)
  effective <- setdiff(rownames(ref), union(noisy, aberrant))
  # Exact-symbol coverage is conservative; the upstream method additionally maps aliases.
  coverage <- sum(effective %in% rownames(x)) / length(effective)
  gene_qc <- data.frame(gene = rownames(ref), in_input = rownames(ref) %in% rownames(x),
    removed_noisy = rownames(ref) %in% noisy, removed_tumor_aberrant = rownames(ref) %in% aberrant)
  if (coverage < min_signature_coverage) stop("Insufficient effective TIL10 signature coverage: ", round(coverage, 3), call. = FALSE)
  selected <- unique(d$sample_barcode)
  fractions <- as.data.frame(quantiseqr::run_quantiseq(x[, selected, drop = FALSE],
    signature_matrix = "TIL10", is_arraydata = FALSE, is_tumordata = TRUE,
    scale_mRNA = TRUE, method = "lsei", rm_genes = "default", return_se = FALSE))
  required <- c("T.cells.CD8", "Tregs", "Neutrophils", "Macrophages.M1", "Macrophages.M2", "T.cells.CD4", "Other")
  if (!all(c("Sample", required) %in% names(fractions)) || anyDuplicated(fractions$Sample) ||
      !setequal(fractions$Sample, selected)) stop("Unexpected quanTIseq output schema/sample IDs.", call. = FALSE)
  f <- as.matrix(fractions[, setdiff(names(fractions), "Sample"), drop = FALSE])
  if (any(!is.finite(f)) || any(f < -1e-8 | f > 1 + 1e-8) || any(abs(rowSums(f) - 1) > 1e-6)) {
    stop("Invalid estimated fractions or composition sum.", call. = FALSE)
  }
  mapping <- list(CD8 = "T.cells.CD8", Treg = "Tregs", Neutrophils = "Neutrophils",
    Macrophages = c("Macrophages.M1", "Macrophages.M2"), CD4 = c("T.cells.CD4", "Tregs"),
    CD4_nonTreg = "T.cells.CD4", Macrophages_M1 = "Macrophages.M1", Macrophages_M2 = "Macrophages.M2")
  j <- match(d$sample_barcode, fractions$Sample)
  values <- lapply(names(mapping), function(cell) data.frame(sample_barcode = d$sample_barcode,
    patient_id = d$patient_id, cohort_label = d$cohort_label, cell_type = cell,
    value = rowSums(fractions[j, mapping[[cell]], drop = FALSE]), measurement = "estimated_fraction",
    method = "quanTIseq_TIL10", unit = "fraction", denominator = "all_modeled_cells_including_Other",
    role = if (cell %in% names(default_immune_signatures())) "primary" else "auxiliary"))
  list(values = as.data.frame(data.table::rbindlist(values)), fractions = fractions, gene_qc = gene_qc,
    sample_qc = data.frame(sample_barcode = fractions$Sample, sum_fractions = rowSums(f),
      immune_fraction = 1 - fractions$Other, status = "ok"),
    settings = list(package = "quantiseqr", version = as.character(utils::packageVersion("quantiseqr")),
      reference_md5 = unname(tools::md5sum(ref_path)), effective_signature_coverage = coverage,
      is_tumordata = TRUE, scale_mRNA = TRUE, method = "lsei", rm_genes = "default",
      CD4 = "non-regulatory CD4 + Treg", Macrophages = "model-defined M1 + M2",
      denominator = "all modeled cells including Other", input = "full-transcriptome linear TPM"))
}

#' Run both marker-score and estimated-fraction survival analyses
#' @param expression Named full-transcriptome human expression matrix.
#' @param clinical One independent tumor specimen per patient per cohort.
#' @param expression_scale tpm or xena_log2_tpm_0.001; explicit, never guessed.
#' @param signatures Named exploratory marker panels.
#' @param endpoints Distinct outcome names; PFI is never substituted for PFS.
#' @param covariate_sets Named list of adjustment sets, including unadjusted.
#' @param run_deconvolution Whether to run optional quanTIseq analysis.
#' @param out_dir Output directory.
#' @param min_patients,min_events,min_group Eligibility thresholds.
#' @return Scores, fractions, combined survival results and module status.
#' @export
run_immune_abundance_analysis <- function(expression, clinical,
    expression_scale = c("tpm", "xena_log2_tpm_0.001"), signatures = default_immune_signatures(),
    endpoints = "OS", covariate_sets = list(unadjusted = character()), run_deconvolution = TRUE,
    out_dir = file.path(getwd(), "results", "immune_abundance"),
    min_patients = 40L, min_events = 10L, min_group = 10L) {
  expression_scale <- match.arg(expression_scale)
  x <- tcgasig_abundance_matrix(expression)
  if (any(!is.finite(x))) stop("Resolve nonfinite expression before unified analysis.", call. = FALSE)
  if (expression_scale == "tpm") {
    if (any(x < 0)) stop("TPM must be nonnegative.", call. = FALSE)
    tpm <- x; logx <- log2(x + 1)
  } else {
    tpm <- 2^x - 0.001
    if (any(tpm < -1e-6)) stop("Values inconsistent with declared Xena log2(TPM+0.001) scale.", call. = FALSE)
    tpm[tpm < 0] <- 0; logx <- x
  }
  if (!is.list(covariate_sets) || !length(covariate_sets) || is.null(names(covariate_sets)) ||
      anyDuplicated(names(covariate_sets)) || any(!nzchar(names(covariate_sets)))) stop("Use named covariate_sets.", call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- score_immune_signatures(logx, clinical, signatures)
  fraction <- if (run_deconvolution) tryCatch(estimate_immune_fractions(tpm, clinical), error = function(e) e) else NULL
  status <- data.frame(module = c("marker_scores", "estimated_fractions"),
    status = c("completed", if (!run_deconvolution) "not_requested" else if (inherits(fraction, "error")) "blocked" else "completed"),
    reason = c("", if (inherits(fraction, "error")) conditionMessage(fraction) else ""))
  values <- marker$values
  if (is.list(fraction)) values <- rbind(values, fraction$values)
  analyses <- lapply(names(covariate_sets), function(name)
    run_cell_abundance_survival(values, clinical, endpoints = endpoints,
      covariates = covariate_sets[[name]], adjustment = name, min_patients = min_patients,
      min_events = min_events, min_group = min_group, write_files = FALSE))
  result <- list(marker = marker, fraction = if (inherits(fraction, "error")) NULL else fraction,
    abundance = values, survival = list(results = as.data.frame(data.table::rbindlist(lapply(analyses, `[[`, "results"))),
      inputs = as.data.frame(data.table::rbindlist(lapply(analyses, `[[`, "inputs"))),
      groups = analyses[[1]]$groups), module_status = status,
    settings = list(expression_scale = expression_scale, endpoints = endpoints, covariate_sets = covariate_sets,
      marker_interpretation = marker$interpretation, direction = "higher score/fraction or high minus low"))
  tables <- list(abundance = values, marker_gene_qc = marker$gene_qc, marker_sample_qc = marker$sample_qc,
    marker_overlap = marker$overlap, survival_results = result$survival$results,
    survival_inputs = result$survival$inputs, fixed_groups = result$survival$groups, module_status = status)
  if (is.list(fraction)) tables <- c(tables, list(full_fractions = fraction$fractions,
    fraction_gene_qc = fraction$gene_qc, fraction_sample_qc = fraction$sample_qc))
  for (name in names(tables)) data.table::fwrite(tables[[name]], file.path(out_dir, paste0(name, ".tsv")), sep = "\t", na = "NA")
  saveRDS(result, file.path(out_dir, "immune_abundance_analysis.rds"))
  result
}
