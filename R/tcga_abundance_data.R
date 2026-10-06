#' Prepare full-transcriptome primary-tumor expression and patient outcomes
#'
#' Uses existing local Xena files; never downloads sequencing data. Clinical
#' conflicts block joining rather than silently selecting a specimen row.
#' Selection is independent of survival availability and of cell estimates.
#' @param project_dir Existing TCGA project containing data/xena.
#' @param cancer_type TCGA abbreviation, default PAAD.
#' @param cache_file Optional new output RDS cache. Fingerprints, sample IDs,
#'   loader version and expression-scale declaration are checked before reuse.
#' @return Expression matrix, clinical table, specimen and gene QC, fingerprints.
#' @export
prepare_tcga_abundance_data <- function(project_dir = getwd(), cancer_type = "PAAD", cache_file = NULL) {
  paths <- c(expression = file.path(project_dir, "data/xena/tcga_RSEM_gene_tpm.gz"),
    annotation = file.path(project_dir, "data/xena/gencode.v23.annotation.gtf.gz"),
    clinical = file.path(project_dir, "data/xena/Survival_SupplementalTable_S1_20171025_xena_sp.tsv"))
  if (!all(file.exists(paths))) stop("Missing local Xena inputs: ", paste(paths[!file.exists(paths)], collapse = ";"), call. = FALSE)
  manifest <- data.frame(role = names(paths), path = normalizePath(paths),
    bytes = file.info(paths)$size, md5 = unname(tools::md5sum(paths)))
  if (!is.null(cache_file) && file.exists(cache_file)) {
    previous <- readRDS(cache_file)
    if (!identical(previous$manifest, manifest) || !identical(previous$cancer_type, cancer_type) ||
        !identical(previous$loader_version, "immune_full_transcriptome_v1")) {
      stop("Existing abundance cache fingerprint/settings mismatch; use a new versioned cache path.", call. = FALSE)
    }
    return(previous)
  }
  source <- as.data.frame(tcgasig_read_delimited(paths["clinical"]), stringsAsFactors = FALSE)
  source <- source[source[["cancer type abbreviation"]] == cancer_type, , drop = FALSE]
  if (!nrow(source)) stop("Cancer absent from clinical source.", call. = FALSE)
  keep <- c("_PATIENT", "age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage", "histological_type")
  keep <- intersect(c(keep, unlist(lapply(c("OS", "PFI", "DSS", "DFI", "PFS"), function(e) c(e, paste0(e, ".time"))))), names(source))
  donor <- unique(source[, keep, drop = FALSE])
  if (anyNA(donor[["_PATIENT"]]) || anyDuplicated(donor[["_PATIENT"]])) stop("Conflicting or missing patient clinical data.", call. = FALSE)
  con <- gzfile(paths["expression"], "rt")
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]; close(con)
  samples <- header[-1L]
  ids <- substr(samples, 1, 12)
  samples <- samples[ids %in% donor[["_PATIENT"]]]
  specimens <- data.frame(sample_barcode = samples, patient_id = substr(samples, 1, 12),
    sample_type_code = substr(samples, 14, 15), stringsAsFactors = FALSE)
  specimens <- specimens[order(specimens$patient_id, specimens$sample_barcode), , drop = FALSE]
  specimens$selection_reason <- ifelse(specimens$sample_type_code == "01", "selected_primary_tumor", "non_primary_tumor_specimen")
  selected <- specimens$selection_reason == "selected_primary_tumor"
  dup <- rep(FALSE, nrow(specimens)); dup[selected] <- duplicated(specimens$patient_id[selected])
  specimens$selection_reason[dup] <- "duplicate_primary_specimen_not_selected"
  metadata <- specimens[specimens$selection_reason == "selected_primary_tumor", , drop = FALSE]
  if (!nrow(metadata)) stop("No primary tumor expression specimens.", call. = FALSE)
  j <- match(metadata$patient_id, donor[["_PATIENT"]])
  metadata$cancer_type <- cancer_type; metadata$cohort_label <- cancer_type
  metadata$age <- tcgasig_to_numeric(donor$age_at_initial_pathologic_diagnosis[j])
  metadata$sex <- tcgasig_clean_missing(donor$gender[j])
  metadata$stage_raw <- tcgasig_clean_missing(donor$ajcc_pathologic_tumor_stage[j])
  metadata$histological_type <- tcgasig_clean_missing(donor$histological_type[j])
  stage <- toupper(trimws(metadata$stage_raw))
  metadata$stage_group <- ifelse(grepl("^STAGE (III|IV)", stage), "III_IV",
    ifelse(grepl("^STAGE (I|II)([AB]|$)", stage), "I_II", NA_character_))
  for (endpoint in c("OS", "PFI", "DSS", "DFI", "PFS")) {
    if (all(c(endpoint, paste0(endpoint, ".time")) %in% names(donor))) {
      metadata[[paste0(endpoint, ".status")]] <- tcgasig_to_numeric(donor[[endpoint]][j])
      metadata[[paste0(endpoint, ".time")]] <- tcgasig_to_numeric(donor[[paste0(endpoint, ".time")]][j])
    }
  }
  message("Streaming full transcriptome for ", nrow(metadata), " independent ", cancer_type, " patients.")
  full <- read_tcga_expression_subset(paths["expression"], paths["annotation"], metadata$sample_barcode, gene_type = NULL)
  minimum <- min(full$expression)
  if (abs(minimum - log2(0.001)) > 0.01) stop("Xena minimum differs from declared log2(TPM+0.001) baseline; verify input scale.", call. = FALSE)
  tpm <- 2^full$expression - 0.001
  if (any(tpm < -1e-6)) stop("Negative back-transformed TPM beyond rounding tolerance.", call. = FALSE)
  out <- list(expression = full$expression, clinical = metadata, specimen_qc = specimens,
    gene_qc = full$gene_qc, manifest = manifest, cancer_type = cancer_type,
    expression_scale = "xena_log2_tpm_0.001", loader_version = "immune_full_transcriptome_v1",
    transform_qc = data.frame(raw_min = minimum, raw_max = max(full$expression),
      negative_rounding_values = sum(tpm < 0), min_backtransformed = min(tpm),
      gene_symbols = nrow(full$expression), patients = ncol(full$expression)),
    source_urls = c(expression = "https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz",
      clinical = "https://pancanatlas.xenahubs.net/Survival_SupplementalTable_S1_20171025_xena_sp"))
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cache_file)
  }
  out
}
