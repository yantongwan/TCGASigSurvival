#' Build explicitly labelled TCGA bulk cell-context data
#'
#' Reuses a prepared target/marker cache, but selects specimens without an OS
#' eligibility filter. Primary tumors are preferred; SKCM permits metastases.
#' Adjacent normal (barcode 11) is not a healthy donor. Technical specimens are
#' selected deterministically and every exclusion is returned in specimen_qc.
#' @param prepared Prepared target/marker data list or RDS path.
#' @param clinical_file TCGA-CDR Xena table containing all outcome endpoints.
#' @param cell_type Cell-context label.
#' @param cell_genes Marker symbols. The target gene should not be in this set.
#' @param min_cell_genes Minimum usable markers.
#' @return tcgasig_cell_data with bulk_proxy measurement and selection QC.
#' @export
prepare_tcga_cell_context <- function(prepared, clinical_file, cell_type = "Treg",
    cell_genes = default_treg_signature(), min_cell_genes = 6L) {
  if (is.character(prepared)) prepared <- readRDS(prepared)
  tcgasig_validate_threshold(min_cell_genes, "min_cell_genes")
  expr <- data.table::as.data.table(prepared$expression)
  meta <- as.data.frame(prepared$sample_metadata, stringsAsFactors = FALSE)
  required <- c("sample_barcode", "patient_id", "sample_type_code", "cancer_type")
  if (!all(required %in% names(meta)) || !"gene_name" %in% names(expr)) stop("Cache lacks expression/specimen metadata.", call. = FALSE)
  samples <- intersect(meta$sample_barcode, names(expr))
  meta <- meta[meta$sample_barcode %in% samples, , drop = FALSE]
  meta$selection_reason <- ifelse(is.na(meta$cancer_type), "cancer_type_missing",
    ifelse(meta$sample_type_code %in% c("01", "11") |
      (meta$cancer_type == "SKCM" & meta$sample_type_code == "06"), "eligible", "specimen_type_out_of_scope"))
  meta$tissue_status <- ifelse(meta$sample_type_code == "11", "adjacent_normal", "tumor")
  meta <- meta[order(meta$cancer_type, meta$patient_id, meta$tissue_status,
    meta$sample_type_code, meta$sample_barcode, na.last = TRUE), , drop = FALSE]
  key <- paste(meta$cancer_type, meta$patient_id, meta$tissue_status, sep = "\r")
  eligible <- meta$selection_reason == "eligible"
  duplicate <- rep(FALSE, nrow(meta)); duplicate[eligible] <- duplicated(key[eligible])
  meta$selection_reason[duplicate] <- "duplicate_specimen_not_selected"
  specimen_qc <- meta[, c(required, "tissue_status", "selection_reason")]
  meta <- meta[meta$selection_reason == "eligible", required, drop = FALSE]
  meta$tissue_status <- ifelse(meta$sample_type_code == "11", "adjacent_normal", "tumor")
  clin <- as.data.frame(tcgasig_read_delimited(clinical_file), stringsAsFactors = FALSE)
  needed <- c("_PATIENT", "cancer type abbreviation")
  if (!all(needed %in% names(clin))) stop("Clinical table lacks patient/cancer keys.", call. = FALSE)
  clin$key <- paste(clin[["_PATIENT"]], clin[["cancer type abbreviation"]], sep = "\r")
  clinical_columns <- intersect(c("key", "age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage",
    unlist(lapply(c("OS", "PFI", "PFS", "DSS", "DFI"), function(e) c(e, paste0(e, ".time"))))), names(clin))
  consensus <- unique(clin[, clinical_columns, drop = FALSE])
  if (anyDuplicated(consensus$key)) stop("Conflicting clinical donor outcomes/covariates; resolve before joining.", call. = FALSE)
  clinical_qc <- as.data.frame(table(clin$key), stringsAsFactors = FALSE)
  names(clinical_qc) <- c("patient_cancer_key", "source_specimen_rows")
  clin <- clin[!duplicated(clin$key), , drop = FALSE]
  jj <- match(paste(meta$patient_id, meta$cancer_type, sep = "\r"), clin$key)
  for (endpoint in c("OS", "PFI", "PFS", "DSS", "DFI")) {
    if (all(c(endpoint, paste0(endpoint, ".time")) %in% names(clin))) {
      meta[[paste0(endpoint, ".status")]] <- tcgasig_to_numeric(clin[[endpoint]][jj])
      meta[[paste0(endpoint, ".time")]] <- tcgasig_to_numeric(clin[[paste0(endpoint, ".time")]][jj])
    }
  }
  if ("age_at_initial_pathologic_diagnosis" %in% names(clin)) meta$age <- tcgasig_to_numeric(clin$age_at_initial_pathologic_diagnosis[jj])
  if ("gender" %in% names(clin)) meta$sex <- tcgasig_clean_missing(clin$gender[jj])
  if ("ajcc_pathologic_tumor_stage" %in% names(clin)) meta$stage <- tcgasig_clean_missing(clin$ajcc_pathologic_tumor_stage[jj])
  meta$profile_id <- meta$sample_barcode
  meta$cell_type <- cell_type
  meta$cell_score <- NA_real_
  gene_qc <- list()
  for (cohort in unique(meta$cancer_type)) {
    ii <- which(meta$cancer_type == cohort)
    tumor <- ii[meta$tissue_status[ii] == "tumor"]
    markers <- list()
    used <- character()
    for (gene in intersect(cell_genes, expr$gene_name)) {
      v <- tcgasig_extract_gene_vector(expr, gene, meta$sample_barcode[ii])
      ref <- tcgasig_extract_gene_vector(expr, gene, meta$sample_barcode[tumor])
      sd <- stats::sd(ref, na.rm = TRUE)
      if (length(ref) > 1L && is.finite(sd) && sd > 0) {
        markers[[gene]] <- (v - mean(ref, na.rm = TRUE)) / sd
        used <- c(used, gene)
      }
    }
    if (length(markers) >= min_cell_genes) {
      z <- do.call(cbind, markers)
      score <- rowMeans(z, na.rm = TRUE)
      score[rowSums(is.finite(z)) < min_cell_genes | !is.finite(score)] <- NA_real_
      meta$cell_score[ii] <- score
    }
    gene_qc[[cohort]] <- data.frame(cohort_label = cohort, gene = cell_genes, usable = cell_genes %in% used)
  }
  x <- as.matrix(expr[, meta$profile_id, with = FALSE])
  rownames(x) <- expr$gene_name
  data <- prepare_cell_expression_data(x, meta, "bulk_proxy", "log",
    "TCGA Toil bulk RNA; marker z scores referenced to all tumor donors per cancer; TCGA-CDR outcomes. Not cell-specific expression or true deconvolution.")
  data$specimen_qc <- specimen_qc
  data$marker_qc <- as.data.frame(data.table::rbindlist(gene_qc))
  data$cell_genes <- cell_genes
  data$clinical_qc <- clinical_qc
  data
}

#' Read a sample subset of the TCGA Xena expression matrix
#' @param expression_file Local gzipped TCGA Toil expression file.
#' @param annotation_file Matching local GENCODE v23 GTF gzip file.
#' @param sample_ids Exact selected expression column names.
#' @param gene_type GENCODE gene biotype to retain, default protein_coding.
#' @return List with named log-expression matrix and gene mapping/QC. Duplicate
#'   symbols retain the row with highest mean across selected samples, recorded
#'   in gene_qc. This loader assumes the documented Xena log scale, never counts.
#' @export
read_tcga_expression_subset <- function(expression_file, annotation_file, sample_ids,
    gene_type = "protein_coding") {
  if (!length(sample_ids) || anyNA(sample_ids) || anyDuplicated(sample_ids)) stop("Unique sample_ids required.", call. = FALSE)
  gtf <- tcgasig_stream_rows(annotation_file, gtf = TRUE)
  map <- data.frame(gene_id = tcgasig_extract_attr(gtf[[9]], "gene_id"),
    gene_name = tcgasig_extract_attr(gtf[[9]], "gene_name"),
    gene_type = tcgasig_extract_attr(gtf[[9]], "gene_type"), stringsAsFactors = FALSE)
  map <- unique(map[map$gene_type %in% gene_type & !is.na(map$gene_name), ])
  con <- if (grepl("[.]gz$", expression_file)) gzfile(expression_file, "rt") else file(expression_file, "rt")
  on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
  if (!all(sample_ids %in% header)) stop("Some requested samples are absent from Xena expression.", call. = FALSE)
  jj <- match(sample_ids, header)
  rows <- list()
  ids <- character()
  repeat {
    lines <- readLines(con, n = 20L, warn = FALSE)
    if (!length(lines)) break
    keys <- sub("\t.*$", "", lines)
    keep <- keys %in% map$gene_id
    for (i in which(keep)) {
      fields <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
      if (length(fields) != length(header)) stop("Malformed expression row: ", keys[i], call. = FALSE)
      rows[[length(rows) + 1L]] <- as.numeric(fields[jj])
      ids <- c(ids, keys[i])
    }
  }
  if (!length(rows)) stop("No expression genes matched annotation.", call. = FALSE)
  x <- do.call(rbind, rows)
  colnames(x) <- sample_ids
  gene_qc <- map[match(ids, map$gene_id), , drop = FALSE]
  gene_qc$mean_selected_expression <- rowMeans(x)
  order <- order(gene_qc$gene_name, -gene_qc$mean_selected_expression, gene_qc$gene_id)
  retained <- !duplicated(gene_qc$gene_name[order])
  gene_qc$retained <- FALSE; gene_qc$retained[order[retained]] <- TRUE
  x <- x[order[retained], , drop = FALSE]
  rownames(x) <- gene_qc$gene_name[order[retained]]
  list(expression = x, gene_qc = gene_qc, expression_scale = "Xena_provided_log")
}
