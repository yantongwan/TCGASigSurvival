test_that("custom inputs exclude invalid OS and reject repeated patients", {
  demo <- tcgasig_demo_data(60)
  demo$clinical$OS.status[1] <- 2
  demo$clinical$OS.time[2] <- 0
  prepared <- prepare_expression_survival_data(demo$expression, demo$clinical)
  expect_equal(nrow(prepared$cohort_samples), 118L)
  expect_equal(prepared$input_qc$invalid_survival_rows, 2L)
  demo$clinical$patient_id[4] <- demo$clinical$patient_id[3]
  expect_error(prepare_expression_survival_data(demo$expression, demo$clinical), "Multiple samples")
  prepared <- prepare_expression_survival_data(demo$expression, demo$clinical, duplicate_patients = "first_sample")
  expect_equal(nrow(prepared$cohort_samples), 117L)
  expect_true(demo$clinical$sample_barcode[3] %in% prepared$cohort_samples$sample_barcode)
})

test_that("synthetic examples preserve caller RNG", {
  set.seed(99)
  before <- .Random.seed
  x <- tcgasig_demo_data(60)
  expect_identical(.Random.seed, before)
  expect_identical(x, tcgasig_demo_data(60))
})

test_that("gene files normalize case only when requested", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("Gene", "Jun", "JUN", "", "NA", "Dusp1"), path)
  expect_equal(read_signature_genes(path), c("Jun", "JUN", "Dusp1"))
  expect_equal(read_signature_genes(path, uppercase = TRUE), c("JUN", "DUSP1"))
  expect_error(read_signature_genes(path, "absent"), "Missing gene column")
  unlink(path)
})

test_that("R streaming extracts only selected rows from gzip", {
  path <- tempfile(fileext = ".gz")
  con <- gzfile(path, "wt")
  writeLines(c("gene\tS1\tS2", "ENSG1\t1\t2", "ENSG2\t3\t4", "ENSG3\t5\t6"), con)
  close(con)
  x <- TCGASigSurvival:::tcgasig_stream_rows(path, keys = "ENSG2")
  expect_equal(x[[1]], "ENSG2")
  expect_equal(as.numeric(x[1, 2:3]), c(3, 4))
  unlink(path)
})
