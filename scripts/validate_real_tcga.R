library(TCGASigSurvival)
args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) args[1] else stop("Pass your TCGA project directory.")
out <- file.path(project, "release_validation", "real_tcga")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# Exercise the public wrapper with canonical markers, without redistributing
# the project's unpublished state signature.
result <- run_pan_tcga_two_layer_survival(
  project_dir = project, cell_signature_genes = default_naive_cd4_signature(),
  state_signature_genes = default_cd8_signature(),
  cell_name = "NaiveCD4", state_name = "CD8ModuleSmokeTest",
  min_state_genes = 4, min_cell_genes = 6,
  cutoff_methods = "median", make_km = FALSE, make_summary_plots = FALSE,
  out_dir = out
)
stopifnot(nrow(result$continuous) >= 10L, nrow(result$adjusted) >= 10L)
stopifnot(all(is.finite(result$continuous$HR)), all(is.finite(result$adjusted$HR)))
inventory <- summarize_tcga_cohorts(result$prepared)
data.table::fwrite(inventory, file.path(out, "cohort_inventory.tsv"), sep = "\t")
writeLines(c("Public TCGA two-layer wrapper smoke test passed.",
  paste("Default cancer cohorts:", length(unique(inventory$cancer_type[inventory$default_cohort]))),
  paste("Continuous models:", nrow(result$continuous)),
  paste("Adjusted median models:", nrow(result$adjusted)),
  paste("Subset median models:", nrow(result$cell_high)),
  "This is a software validation run, not a trained-immunity analysis."), file.path(out, "VALIDATION.txt"))
cat("Real TCGA smoke test passed: ", out, "\n", sep = "")
