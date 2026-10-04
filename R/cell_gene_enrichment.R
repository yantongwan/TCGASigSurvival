#' Differential expression and preranked enrichment of cell-gene groups
#'
#' Tests use independent donor profiles, never individual cells. Raw counts use
#' edgeR filtering/TMM and limma-voom; normalized log expression uses limma-trend.
#' Gene groups are median based and defined before covariate missingness. The
#' grouping gene is excluded from DE/enrichment by default to avoid circular
#' evidence. Bulk-proxy DE is whole-tissue DE, not Treg-intrinsic DE.
#' @param data Prepared donor-level cell expression data.
#' @param target_gene,cell_type Exact gene and cell label.
#' @param pathways Named list of gene sets with symbols matching expression.
#'   No internet downloads are performed by this function.
#' @param cohorts Cohort labels to analyze, default all available.
#' @param covariates Clinical adjustment columns.
#' @param context_subset all or high, the latter requiring cell_score.
#' @param exclude_genes Genes to exclude from both testing and ranked universe.
#' @param min_group Minimum biological donors per expression group.
#' @param min_size,max_size Pathway size bounds after matching tested genes.
#' @param seed Deterministic fgsea seed; caller RNG state is restored.
#' @param allow_proxy Explicit permission to analyze a bulk proxy.
#' @param out_dir,prefix,write_files Output settings.
#' @param score_mode gene or joint target/cell-context score. Joint-mode
#'   marker genes, when documented in data$cell_genes, are also excluded.
#' @return Differential tables, all tested enrichment results, leading-edge
#'   gene statistics, donor inputs, QC and ranks. Positive logFC/NES favors
#'   high expression. GSEA uses all tested genes, not only significant DEGs.
#' @export
run_cell_gene_enrichment <- function(data, target_gene = "ACLY", cell_type = "Treg",
    pathways, cohorts = NULL, covariates = character(), context_subset = c("all", "high"),
    exclude_genes = target_gene, min_group = 10L, min_size = 15L, max_size = 500L,
    seed = 42L, allow_proxy = FALSE, out_dir = file.path(getwd(), "results"),
    prefix = "cell_gene", write_files = TRUE, score_mode = c("gene", "joint")) {
  tcgasig_cell_check(data, allow_proxy)
  tcgasig_need("limma"); tcgasig_need("fgsea"); tcgasig_need("BiocParallel")
  if (data$expression_scale == "counts") tcgasig_need("edgeR")
  context_subset <- match.arg(context_subset)
  score_mode <- match.arg(score_mode)
  if (score_mode == "joint") exclude_genes <- unique(c(exclude_genes, data$cell_genes))
  tcgasig_validate_threshold(min_group, "min_group", 2L)
  tcgasig_validate_threshold(min_size, "min_size")
  tcgasig_validate_threshold(max_size, "max_size", min_size)
  if (!is.list(pathways) || !length(pathways) || is.null(names(pathways)) ||
      anyNA(names(pathways)) || any(!nzchar(names(pathways))) || anyDuplicated(names(pathways))) {
    stop("pathways must be a nonempty, uniquely named gene-set list.", call. = FALSE)
  }
  pathways <- lapply(pathways, function(x) unique(as.character(x[!is.na(x)])))
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  groups <- tcgasig_cell_groups(data, target_gene, cell_type, context_subset, score_mode)
  if (is.null(cohorts)) cohorts <- unique(groups$cohort_label)
  if (!all(cohorts %in% groups$cohort_label) || !length(cohorts)) stop("Unknown/empty cohorts.", call. = FALSE)
  rows <- enrich <- edges <- qc <- inputs <- ranks <- list()
  for (cohort in sort(unique(cohorts))) {
    initial <- groups[groups$cohort_label == cohort, , drop = FALSE]
    candidate <- initial[!is.na(initial$group), , drop = FALSE]
    cc <- tcgasig_covariate_rows(candidate, unique(covariates))
    d <- cc$data
    initial$exclusion <- ifelse(!initial$context_selected, "outside_context_subset",
      ifelse(is.na(initial$group), "group_not_estimable",
        ifelse(!initial$profile_id %in% d$profile_id, "covariate_missing", "included")))
    inputs[[cohort]] <- initial
    high <- sum(d$group == "high"); low <- sum(d$group == "low")
    q <- data.frame(cohort_label = cohort, input_n = nrow(initial), analyzed_n = nrow(d),
      high_n = high, low_n = low, measurement = data$measurement,
      method = if (data$expression_scale == "counts") "TMM_limma_voom" else "limma_trend", score_mode = score_mode,
      covariates_used = paste(cc$active, collapse = ";"),
      covariates_constant = paste(cc$constant, collapse = ";"),
      excluded_genes = paste(exclude_genes, collapse = ";"), genes_tested = 0L,
      pathways_tested = 0L, rank_ties = NA_integer_, gsea_warnings = "", status = "ok", reason = "")
    if (min(high, low) < min_group) { q$status <- "not_estimable"; q$reason <- "insufficient_or_tied_groups"; qc[[cohort]] <- q; next }
    design <- stats::model.matrix(stats::reformulate(c("group", cc$active)), d)
    if (qr(design)$rank < ncol(design) || nrow(design) - ncol(design) < 2L) {
      q$status <- "not_estimable"; q$reason <- "rank_deficient_or_no_residual_degrees"; qc[[cohort]] <- q; next
    }
    x <- data$expression[, match(d$profile_id, colnames(data$expression)), drop = FALSE]
    if (data$expression_scale == "counts") {
      y <- edgeR::DGEList(x)
      keep <- edgeR::filterByExpr(y, design = design) & !rownames(x) %in% exclude_genes
      y <- y[keep, , keep.lib.sizes = FALSE]
      if (nrow(y) < min_size) { q$status <- "not_estimable"; q$reason <- "too_few_tested_genes"; qc[[cohort]] <- q; next }
      y <- edgeR::calcNormFactors(y)
      fit <- limma::eBayes(limma::lmFit(limma::voom(y, design, plot = FALSE), design))
    } else {
      keep <- apply(x, 1, stats::sd) > 0 & !rownames(x) %in% exclude_genes
      x <- x[keep, , drop = FALSE]
      if (nrow(x) < min_size) { q$status <- "not_estimable"; q$reason <- "too_few_tested_genes"; qc[[cohort]] <- q; next }
      fit <- limma::eBayes(limma::lmFit(x, design), trend = TRUE)
    }
    de <- limma::topTable(fit, coef = "grouphigh", number = Inf, sort.by = "none")
    de$gene <- rownames(de)
    de$cohort_label <- cohort
    de$measurement <- data$measurement
    de$exploratory_p_lt_0_05 <- de$P.Value < 0.05
    rows[[cohort]] <- de
    rank <- de$t[is.finite(de$t)]
    names(rank) <- de$gene[is.finite(de$t)]
    # Alphabetical secondary order is stable; do not jitter ties to chase enrichment.
    rank <- rank[order(-rank, names(rank))]
    ranks[[cohort]] <- rank
    q$genes_tested <- length(rank)
    q$rank_ties <- sum(duplicated(rank))
    eligible <- lapply(pathways, intersect, y = names(rank))
    eligible <- eligible[lengths(eligible) >= min_size & lengths(eligible) <= max_size]
    if (!length(eligible)) { q$status <- "enrichment_unavailable"; q$reason <- "no_pathways_match_universe"; qc[[cohort]] <- q; next }
    q$pathways_tested <- length(eligible)
    gsea_warnings <- character()
    gsea <- as.data.frame(withCallingHandlers(
      fgsea::fgseaMultilevel(eligible, rank, minSize = min_size,
        maxSize = max_size, eps = 0, BPPARAM = BiocParallel::SerialParam()),
      warning = function(w) {
        gsea_warnings <<- c(gsea_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }))
    q$gsea_warnings <- paste(unique(gsea_warnings), collapse = " | ")
    if (nrow(gsea)) {
      for (i in seq_len(nrow(gsea))) if (length(gsea$leadingEdge[[i]])) {
        genes <- gsea$leadingEdge[[i]]
        edges[[paste(cohort, i)]] <- data.frame(cohort_label = cohort,
          pathway = gsea$pathway[i], gene = genes, rank_statistic = unname(rank[genes]),
          NES = gsea$NES[i], pval = gsea$pval[i], padj = gsea$padj[i])
      }
      gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, character(1), collapse = ";")
      gsea$cohort_label <- cohort
      gsea$measurement <- data$measurement
      gsea$status <- ifelse(is.finite(gsea$pval) & is.finite(gsea$NES), "ok", "not_estimable")
      gsea$reason <- ifelse(gsea$status == "ok", "", "fgsea_returned_nonfinite_statistics")
      enrich[[cohort]] <- gsea
    }
    qc[[cohort]] <- q
  }
  empty_de <- data.frame(gene = character(), cohort_label = character(), logFC = numeric(),
    t = numeric(), P.Value = numeric(), adj.P.Val = numeric())
  empty_enrich <- data.frame(pathway = character(), cohort_label = character(), NES = numeric(),
    pval = numeric(), padj = numeric(), leadingEdge = character(), status = character())
  empty_edges <- data.frame(cohort_label = character(), pathway = character(), gene = character(),
    rank_statistic = numeric(), NES = numeric(), pval = numeric(), padj = numeric())
  result <- list(differential = if (length(rows)) as.data.frame(data.table::rbindlist(rows, fill = TRUE)) else empty_de,
    enrichment = if (length(enrich)) as.data.frame(data.table::rbindlist(enrich, fill = TRUE)) else empty_enrich,
    leading_edge = if (length(edges)) as.data.frame(data.table::rbindlist(edges, fill = TRUE)) else empty_edges,
    inputs = as.data.frame(data.table::rbindlist(inputs, fill = TRUE)),
    qc = as.data.frame(data.table::rbindlist(qc, fill = TRUE)), ranks = ranks,
    measurement = data$measurement, interpretation = tcgasig_cell_label(data),
    settings = list(direction = "high minus low", rank = "limma moderated t, all tested genes",
      excluded_genes = exclude_genes, seed = seed, score_mode = score_mode, pathways_names = names(pathways)))
  result$files <- if (write_files) tcgasig_write_cell_tables(result, out_dir, prefix) else character()
  result
}

#' Read user-supplied GMT gene sets
#' @param file GMT file (optionally gzip compressed).
#' @return Named pathway list; identifiers are not remapped between species.
#' @export
read_gmt_pathways <- function(file) {
  con <- if (grepl("[.]gz$", file)) gzfile(file, "rt") else file(file, "rt")
  on.exit(close(con))
  rows <- strsplit(readLines(con, warn = FALSE), "\t", fixed = TRUE)
  if (!length(rows) || any(lengths(rows) < 3L)) stop("GMT needs name, description and genes per row.", call. = FALSE)
  names <- vapply(rows, `[`, character(1), 1L)
  if (anyDuplicated(names) || any(!nzchar(names))) stop("Unique nonempty GMT names required.", call. = FALSE)
  stats::setNames(lapply(rows, function(x) unique(x[-c(1, 2)][nzchar(x[-c(1, 2)])])), names)
}
