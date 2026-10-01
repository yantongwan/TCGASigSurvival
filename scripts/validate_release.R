library(TCGASigSurvival)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else file.path(tempdir(), "tcgasig_release_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
demo <- tcgasig_demo_data(120)
prepared <- prepare_expression_survival_data(demo$expression, demo$clinical)
two <- run_two_layer_survival(prepared, demo$state_genes, demo$cell_genes,
  state_name = "DemoState", cutoff_methods = c("median", "optimized"),
  out_dir = file.path(out_dir, "two_layer"))
stopifnot(nrow(two$continuous) == 2L, nrow(two$adjusted) == 4L, nrow(two$cell_high) == 4L)
joint <- run_signature_survival(prepared, target_gene = "ZC3H12C",
  signature_genes = demo$cell_genes, signature_name = "NaiveCD4",
  out_dir = file.path(out_dir, "joint"), cutoff_methods = "median")
stopifnot(nrow(joint$results) == 2L)
plots <- plot_signature_results(results = joint$results, target_gene = "ZC3H12C",
  signature_name = "NaiveCD4", out_dir = file.path(out_dir, "joint"))
expr <- plot_target_expression_by_signature(survival = joint, target_gene = "ZC3H12C",
  signature_name = "NaiveCD4", cutoff_method = "median", out_dir = file.path(out_dir, "joint"))
batch <- run_signature_batch(prepared, list(JUN = "JUN", ATF4 = "ATF4"),
  min_signature_genes = 1, cutoff_methods = "median", out_dir = file.path(out_dir, "batch"))
stopifnot(nrow(batch$results) == 4L)
matched <- demo$clinical
matched$target_rna <- as.numeric(demo$expression["ZC3H12C", ])
matched$target_protein <- matched$target_rna + as.numeric(demo$expression["ATF4", ]) * 0.2
matched$signature_score <- colMeans(demo$expression[demo$cell_genes, ])
protein <- run_matched_proteogenomic_analysis(matched, target_gene = "ZC3H12C",
  out_dir = file.path(out_dir, "protein"))
stopifnot(nrow(protein$correlations) == 2L, nrow(protein$survival) == 16L)
writeLines(c("All synthetic end-to-end checks passed.", utils::capture.output(utils::sessionInfo())),
  file.path(out_dir, "VALIDATION.txt"))
cat("Validation completed: ", out_dir, "\n", sep = "")
