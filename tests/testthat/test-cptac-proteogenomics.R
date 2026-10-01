test_that("default CPTAC to TCGA mapping is available", {
  mapping <- default_cptac_tcga_map()
  expect_true(all(c("cptac_cancer", "tcga_cancer", "note") %in% names(mapping)))
  expect_equal(mapping[cptac_cancer == "CCRCC", tcga_cancer], "KIRC")
  expect_equal(mapping[cptac_cancer == "COAD", tcga_cancer], "COAD")
  expect_equal(nrow(mapping), 10)
})

test_that("CPTAC reader works with a portable result directory", {
  proteogenomics_dir <- tempfile("cptac_fixture_")
  dir.create(file.path(proteogenomics_dir, "results", "qc"), recursive = TRUE)
  dir.create(file.path(proteogenomics_dir, "results", "tables"), recursive = TRUE)
  on.exit(unlink(proteogenomics_dir, recursive = TRUE))
  data.table::fwrite(data.frame(cptac_cancer = "COAD", tcga_cancer = "COAD",
    status = "ok", reason = "synthetic fixture", n_target_rna_protein_complete = 90),
    file.path(proteogenomics_dir, "results", "qc", "cptac_data_audit.tsv"), sep = "\t")
  data.table::fwrite(data.frame(cptac_cancer = "COAD", tcga_cancer = "COAD",
    n = 90, spearman_rho = 0.5, spearman_p = 0.001, spearman_FDR = 0.001),
    file.path(proteogenomics_dir, "results", "tables", "cptac_STRAP_rna_protein_correlations.tsv"), sep = "\t")
  cptac <- load_cptac_proteogenomic_results(proteogenomics_dir = proteogenomics_dir, verbose = FALSE)
  expect_s3_class(cptac, "tcgasig_cptac_proteogenomics")
  expect_true("COAD" %in% cptac$audit$cptac_cancer)
  expect_equal(cptac$rna_protein_correlations$n, 90)
  expect_true("COAD" %in% cptac$rna_protein_correlations$cptac_cancer)

  summary <- summarize_cptac_proteogenomic_results(cptac)
  expect_true("qc_metrics" %in% names(summary))
  expect_equal(summary$qc_metrics[metric == "audited_cptac_cancers", value], 1)
})
