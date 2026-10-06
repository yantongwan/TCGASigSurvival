# Copy verified PDF examples only. No patient tables or expression data are copied.
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args)) args[1] else getwd(), mustWork = TRUE)
package <- if (file.exists(file.path(project, "TCGASigSurvival", "DESCRIPTION")))
  file.path(project, "TCGASigSurvival") else if (file.exists("DESCRIPTION")) normalizePath(".") else
  stop("Cannot locate TCGASigSurvival source directory.")
legacy <- "results/ACLY_Treg_Figure_features_20261003/figures"
immune <- "results/PAAD_Immune_Abundance_20261004/figures"
assets <- data.frame(
  figure = 1:7,
  file = c("tcga_acly_tumor_adjacent.pdf", "brca_acly_treg_hallmark.pdf",
    "pan_tcga_acly_treg_survival.pdf", "paad_immune_abundance_forest.pdf",
    "paad_cd4_marker_score_km.pdf", "paad_cd4_estimated_fraction_km.pdf",
    "paad_ductal_immune_abundance_forest.pdf"),
  source_file = c(file.path(legacy, c("A_bulk_ACLY_paired_adjacent.pdf",
    "B_BRCA_joint_Hallmark_ridges.pdf", "C_joint_survival_grid.pdf")),
    file.path(immune, c("PAAD_OS_continuous_unadjusted.pdf", "PAAD_OS_KM_CD4_marker_score.pdf",
      "PAAD_OS_KM_CD4_estimated_fraction.pdf", "PAAD_Ductal_OS_continuous_unadjusted.pdf"))),
  source_table = c("bulk_ACLY_paired_tests.tsv;bulk_ACLY_paired_values.tsv",
    "BRCA_joint_enrichment.tsv;BRCA_joint_leading_edge.tsv", "ACLY_Treg_joint_survival_results.tsv",
    "analysis/survival_results.tsv", "analysis/survival_results.tsv;analysis/survival_inputs.tsv",
    "analysis/survival_results.tsv;analysis/survival_inputs.tsv", "analysis/survival_results.tsv"),
  stringsAsFactors = FALSE)
sources <- file.path(project, assets$source_file)
if (!all(file.exists(sources))) stop("Missing completed example PDFs: ", paste(sources[!file.exists(sources)], collapse = ";"))
assets$bytes <- file.info(sources)$size
assets$md5 <- unname(tools::md5sum(sources))
dest <- file.path(package, "docs", "figures"); dir.create(dest, recursive = TRUE, showWarnings = FALSE)
targets <- file.path(dest, assets$file)
existing <- file.exists(targets)
if (any(existing) && any(unname(tools::md5sum(targets[existing])) != assets$md5[existing]))
  stop("Existing README PDFs differ; inspect them instead of overwriting.")
if (any(!existing) && !all(file.copy(sources[!existing], targets[!existing], overwrite = FALSE)))
  stop("PDF copy failed.")
stopifnot(identical(unname(tools::md5sum(targets)), assets$md5),
  identical(unname(tools::md5sum(sources)), assets$md5))
utils::write.table(assets, file.path(dest, "manifest.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
cat(nrow(assets), "PDF examples copied and MD5-verified; original analyses unchanged.\n")
