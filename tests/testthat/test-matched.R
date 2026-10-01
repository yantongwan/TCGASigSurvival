test_that("matched tables never join the same identifier across different cancers", {
  rna <- data.frame(cancer_type = c("A", "B"), patient_id = c("P1", "P1"), target_rna = c(1, 2))
  protein <- data.frame(cancer_type = "A", patient_id = "P1", target_protein = 3)
  clinical <- data.frame(cancer_type = "A", patient_id = "P1", OS.time = 10, OS.status = 1)
  result <- match_proteogenomic_tables(rna, protein, clinical)
  expect_equal(nrow(result$data), 1L)
  expect_equal(result$data$cancer_type, "A")
  expect_error(match_proteogenomic_tables(rbind(rna, rna[1, ]), protein, clinical), "unique")
})

test_that("portable proteogenomic statistics use the actual matched pairs", {
  demo <- tcgasig_demo_data(100)
  dt <- demo$clinical
  dt$target_rna <- as.numeric(demo$expression["ZC3H12C", ])
  dt$target_protein <- dt$target_rna + as.numeric(demo$expression["ATF4", ]) * 0.2
  dt$signature_score <- colMeans(demo$expression[demo$cell_genes, ])
  result <- run_matched_proteogenomic_analysis(dt, min_patients = 40, min_events = 10,
    write_files = FALSE, make_plots = FALSE)
  expect_equal(nrow(result$correlations), 2L)
  expect_equal(nrow(result$survival), 16L)
  a <- dt[dt$cancer_type == "SIM_A", ]
  expect_equal(result$correlations[result$correlations$cancer_type == "SIM_A", ]$spearman_rho,
    unname(stats::cor(a$target_rna, a$target_protein, method = "spearman")))
})
