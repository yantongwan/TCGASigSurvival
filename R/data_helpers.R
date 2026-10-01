#' Default Naive CD4 cell-context markers
#'
#' A pragmatic bulk RNA marker module for T-cell and naive/central-memory
#' context. These markers overlap with other lymphocyte states; the score is
#' not a validated cell fraction or a uniquely Naive CD4 measurement.
#' @return A character vector of human gene symbols.
#' @export
default_naive_cd4_signature <- function() {
  c("CD3D", "CD3E", "TRAC", "CD2", "CD4", "IL7R", "CCR7", "TCF7",
    "LEF1", "SELL", "LTB", "MAL")
}

#' Read a gene signature from a delimited file
#' @param file A CSV, TSV or plain text file, optionally gzip compressed.
#' @param gene_column Column name or index in the input table.
#' @param uppercase Whether to uppercase symbols. This is case normalization,
#'   not orthologue mapping between species.
#' @param header Whether the input has a header row.
#' @return Unique nonempty gene symbols.
#' @export
read_signature_genes <- function(file, gene_column = 1L, uppercase = FALSE, header = TRUE) {
  x <- tcgasig_read_delimited(file, header = header)
  if (is.character(gene_column) && !gene_column %in% names(x)) {
    stop("Missing gene column: ", gene_column, call. = FALSE)
  }
  if (is.numeric(gene_column) && (length(gene_column) != 1L ||
      !is.finite(gene_column) || gene_column < 1L || gene_column > ncol(x))) {
    stop("gene_column is outside the input table.", call. = FALSE)
  }
  genes <- tcgasig_clean_missing(as.character(x[[gene_column]]))
  if (isTRUE(uppercase)) genes <- toupper(genes)
  genes <- unique(genes[!is.na(genes)])
  if (!length(genes)) stop("No nonempty genes were found.", call. = FALSE)
  genes
}

tcgasig_read_delimited <- function(file, header = TRUE) {
  if (!file.exists(file)) stop("Missing input file: ", file, call. = FALSE)
  if (!grepl("\\.gz$", file, ignore.case = TRUE)) {
    return(data.table::fread(file, header = header))
  }
  con <- gzfile(file, "rt")
  on.exit(close(con))
  data.table::fread(text = paste(readLines(con, warn = FALSE), collapse = "\n"), header = header)
}

tcgasig_has_unix_streaming <- function() {
  .Platform$OS.type != "windows" && all(nzchar(Sys.which(c("gzip", "awk"))))
}

tcgasig_stream_rows <- function(path, keys = NULL, gtf = FALSE) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  selected <- list()
  if (!gtf) selected[[1L]] <- readLines(con, n = 1L, warn = FALSE)
  repeat {
    lines <- readLines(con, n = if (gtf) 10000L else 100L, warn = FALSE)
    if (!length(lines)) break
    if (gtf) {
      keep <- !startsWith(lines, "#") & grepl("^[^\t]*\t[^\t]*\tgene\t", lines)
    } else keep <- sub("\t.*$", "", lines) %in% keys
    if (any(keep)) selected[[length(selected) + 1L]] <- lines[keep]
  }
  text <- paste(unlist(selected, use.names = FALSE), collapse = "\n")
  if (!nzchar(text)) stop("No rows extracted from ", path, call. = FALSE)
  data.table::fread(text = text, sep = "\t", header = !gtf, quote = "", fill = gtf)
}

#' Prepare a custom expression and survival cohort
#' @param expression Numeric gene-by-sample matrix with gene symbols as row
#'   names and sample identifiers as column names.
#' @param clinical Data frame with sample_barcode, patient_id, cancer_type,
#'   OS.time (days) and OS.status (0=censored, 1=dead). Additional covariates
#'   are retained. Optional cohort_label, sample_mode and default_cohort are
#'   retained; otherwise they are filled from cancer_type, custom and TRUE.
#' @param expression_scale Either log (already normalized log expression) or
#'   tpm (nonnegative TPM, transformed with log2(TPM+1)).
#' @param duplicate_patients Either error or first_sample, within each cohort.
#' @return A prepared-data list accepted by the survival analysis functions.
#' @export
prepare_expression_survival_data <- function(expression, clinical,
    expression_scale = c("log", "tpm"), duplicate_patients = c("error", "first_sample")) {
  expression_scale <- match.arg(expression_scale)
  duplicate_patients <- match.arg(duplicate_patients)
  expression <- as.matrix(expression)
  if (!is.numeric(expression) || !nrow(expression) || !ncol(expression) ||
      is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression must be a named numeric gene-by-sample matrix.", call. = FALSE)
  }
  if (anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression)) ||
      anyNA(rownames(expression)) || anyNA(colnames(expression)) ||
      any(!nzchar(rownames(expression))) || any(!nzchar(colnames(expression)))) {
    stop("Expression gene and sample names must be nonempty and unique.", call. = FALSE)
  }
  required <- c("sample_barcode", "patient_id", "cancer_type", "OS.time", "OS.status")
  cs <- data.table::copy(data.table::as.data.table(clinical))
  if (!all(required %in% names(cs))) stop("clinical requires: ", paste(required, collapse = ", "), call. = FALSE)
  for (column in required[1:3]) {
    cs[[column]] <- tcgasig_clean_missing(as.character(cs[[column]]))
    if (anyNA(cs[[column]])) stop("Missing clinical identifiers in ", column, call. = FALSE)
  }
  if (anyDuplicated(cs$sample_barcode)) stop("Clinical sample identifiers must be unique.", call. = FALSE)
  if (!"cohort_label" %in% names(cs)) cs$cohort_label <- cs$cancer_type
  if (!"sample_mode" %in% names(cs)) cs$sample_mode <- "custom"
  if (!"default_cohort" %in% names(cs)) cs$default_cohort <- TRUE
  if (anyNA(cs$cohort_label) || anyNA(cs$default_cohort) || !is.logical(cs$default_cohort)) {
    stop("cohort_label and logical default_cohort must be nonmissing.", call. = FALSE)
  }
  cs$OS.time <- tcgasig_to_numeric(cs$OS.time)
  cs$OS.status <- tcgasig_to_numeric(cs$OS.status)
  cs <- cs[cs$sample_barcode %in% colnames(expression), ]
  invalid <- !is.finite(cs$OS.time) | cs$OS.time <= 0 | !cs$OS.status %in% c(0, 1)
  invalid_n <- sum(invalid)
  cs <- cs[!invalid, ]
  if (!nrow(cs)) stop("No expression samples have valid survival information.", call. = FALSE)
  keys <- paste(cs$cohort_label, cs$patient_id, sep = "\r")
  duplicated_n <- sum(duplicated(keys))
  if (duplicated_n && duplicate_patients == "error") {
    stop("Multiple samples per patient within a cohort; choose first_sample explicitly.", call. = FALSE)
  }
  data.table::setorderv(cs, c("cohort_label", "patient_id", "sample_barcode"))
  cs <- cs[!duplicated(cs[, c("cohort_label", "patient_id"), with = FALSE]), ]
  if (expression_scale == "tpm") {
    if (any(expression < 0, na.rm = TRUE)) stop("TPM values cannot be negative.", call. = FALSE)
    expression <- log2(expression + 1)
  }
  expr <- data.table::as.data.table(expression)
  expr$gene_name <- rownames(expression)
  expr$gene_id <- rownames(expression)
  expr$gene_type <- "user_supplied"
  data.table::setcolorder(expr, c("gene_name", "gene_id", "gene_type", colnames(expression)))
  list(expression = expr, cohort_samples = cs, sample_metadata = cs,
       signature_genes = rownames(expression), signature_name = "Custom",
       expression_scale = if (expression_scale == "tpm") "log2_tpm_plus_1" else "user_supplied_log",
       input_qc = list(invalid_survival_rows = invalid_n, duplicate_patient_rows = duplicated_n))
}

#' Inventory patients and survival events by cohort
#' @param prepared A prepared-data list or its RDS path.
#' @return A data frame with cancer, cohort, sample mode, patients and deaths.
#' @export
summarize_tcga_cohorts <- function(prepared) {
  if (is.character(prepared)) prepared <- readRDS(prepared)
  cs <- data.table::as.data.table(prepared$cohort_samples)
  if (!nrow(cs)) return(data.frame())
  cs[, .(n_samples = .N, n_patients = data.table::uniqueN(patient_id),
         n_events = sum(OS.status == 1, na.rm = TRUE)),
     by = .(cancer_type, cohort_label, sample_mode, default_cohort)]
}

#' Public TCGA data source manifest
#' @return A data frame of file roles, names and public source URLs.
#' @export
tcga_data_sources <- function() {
  data.frame(
    role = c("expression", "survival", "phenotype", "annotation"),
    filename = c("tcga_RSEM_gene_tpm.gz", "Survival_SupplementalTable_S1_20171025_xena_sp.tsv",
                 "TcgaTargetGTEX_phenotype.txt.gz", "gencode.v23.annotation.gtf.gz"),
    url = c("https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz",
            "https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp",
            "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz",
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.annotation.gtf.gz"),
    stringsAsFactors = FALSE)
}

#' Download public TCGA inputs explicitly
#' Downloads are never triggered by package installation or loading. Partial
#' files are kept under a .part suffix; only a successful transfer is promoted.
#' Existing files are reused and recorded with a local MD5 fingerprint. That
#' fingerprint is not a comparison against a publisher checksum.
#' @param project_dir Directory containing data/xena.
#' @param data_dir Destination directory.
#' @param roles Source roles to download.
#' @param overwrite Whether to replace existing completed files.
#' @param timeout Download timeout in seconds.
#' @param quiet Suppress download progress.
#' @return A manifest with paths, sizes, local MD5 and retrieval time.
#' @export
download_tcga_data <- function(project_dir = getwd(),
    data_dir = file.path(project_dir, "data", "xena"),
    roles = c("expression", "survival", "phenotype", "annotation"),
    overwrite = FALSE, timeout = 3600, quiet = FALSE) {
  sources <- tcga_data_sources()
  if (!length(roles) || !all(roles %in% sources$role)) stop("Unknown or empty roles.", call. = FALSE)
  if (length(timeout) != 1 || !is.finite(timeout) || timeout <= 0) stop("timeout must be positive.", call. = FALSE)
  sources <- sources[sources$role %in% roles, ]
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  old <- options(timeout = max(timeout, getOption("timeout")))
  on.exit(options(old))
  sources$path <- file.path(data_dir, sources$filename)
  sources$status <- "reused_existing"
  for (i in seq_len(nrow(sources))) {
    dest <- sources$path[i]
    if (!file.exists(dest) || overwrite) {
      part <- paste0(dest, ".part")
      status <- utils::download.file(sources$url[i], part, method = "libcurl", mode = "wb", quiet = quiet)
      if (status != 0 || !file.exists(part) || file.info(part)$size == 0) {
        stop("Download failed; partial file retained: ", part, call. = FALSE)
      }
      if (grepl("\\.gz$", dest)) {
        con <- file(part, "rb")
        magic <- readBin(con, "raw", 2L)
        close(con)
        if (!identical(magic, as.raw(c(31, 139)))) stop("Downloaded file lacks a gzip header: ", part, call. = FALSE)
      }
      if (!file.rename(part, dest)) stop("Cannot promote completed download: ", part, call. = FALSE)
      sources$status[i] <- "downloaded"
    }
  }
  sources$bytes <- file.info(sources$path)$size
  sources$local_md5 <- unname(tools::md5sum(sources$path))
  sources$recorded_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  data.table::fwrite(sources, file.path(data_dir, "TCGASigSurvival_download_manifest.tsv"), sep = "\t")
  sources
}

#' Generate a reproducible synthetic example
#' No patient or unpublished signature data are included. Random state is
#' restored after generation.
#' @param n Patients per synthetic cohort, at least 40.
#' @param seed Random seed.
#' @return A list with expression, clinical, cell_genes and state_genes.
#' @examples
#' demo <- tcgasig_demo_data(n = 80)
#' prepared <- prepare_expression_survival_data(demo$expression, demo$clinical)
#' summarize_tcga_cohorts(prepared)
#' @export
tcgasig_demo_data <- function(n = 120L, seed = 42L) {
  if (length(n) != 1L || !is.finite(n) || n < 40L) stop("n must be at least 40.", call. = FALSE)
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  n <- as.integer(n)
  cell_genes <- default_naive_cd4_signature()
  state_genes <- c("ATF4", "DUSP1", "JUN", "TNFAIP3", "GADD45A", "SGK1")
  genes <- c(cell_genes, state_genes, "ZC3H12C")
  ids <- sprintf("SIM_%04d", seq_len(2L * n))
  cell <- stats::rnorm(2L * n)
  state <- 0.3 * cell + stats::rnorm(2L * n)
  expr <- matrix(stats::rnorm(length(genes) * length(ids), sd = 0.5),
                 nrow = length(genes), dimnames = list(genes, ids)) + 5
  expr[cell_genes, ] <- sweep(expr[cell_genes, ], 2L, cell, "+")
  expr[state_genes, ] <- sweep(expr[state_genes, ], 2L, state, "+")
  death <- stats::rexp(length(ids), rate = exp(0.35 * state - 0.2 * cell) / 1000)
  censor <- stats::rexp(length(ids), rate = 1 / 1700)
  clinical <- data.frame(sample_barcode = ids, patient_id = ids,
    cancer_type = rep(c("SIM_A", "SIM_B"), each = n),
    OS.time = pmin(death, censor), OS.status = as.integer(death <= censor),
    age = round(stats::runif(length(ids), 35, 80)))
  list(expression = expr, clinical = clinical, cell_genes = cell_genes, state_genes = state_genes)
}
