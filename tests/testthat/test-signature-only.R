test_that("signature-only survival API is available", {
  expect_true(is.function(run_pan_tcga_signature_only_survival))
  expect_true(is.function(run_signature_only_survival))
  expect_true(is.function(prepare_signature_only_data))
})

test_that("signature-only output prefix is stable", {
  expect_equal(
    TCGASigSurvival:::tcgasig_signature_only_output_prefix("Naive CD4 Trained Immunity"),
    "panTCGA_Naive_CD4_Trained_Immunity_signature_only"
  )
})
