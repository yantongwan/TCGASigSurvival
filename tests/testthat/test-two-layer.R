two_layer_demo <- function() {
  demo <- tcgasig_demo_data(100)
  list(demo = demo, prepared = prepare_expression_survival_data(demo$expression, demo$clinical))
}

run_test_layers <- function(x, ...) {
  run_two_layer_survival(x$prepared, x$demo$state_genes, x$demo$cell_genes,
    cutoff_methods = "median", min_patients = 40, min_events = 10,
    min_subset_patients = 20, min_subset_events = 5,
    make_km = FALSE, make_summary_plots = FALSE, write_files = FALSE, verbose = FALSE, ...)
}

test_that("two-layer estimates agree with an independent Cox fit", {
  x <- two_layer_demo()
  result <- run_test_layers(x, covariates = "age")
  expect_equal(nrow(result$adjusted), 2L)
  expect_equal(nrow(result$continuous), 2L)
  expect_equal(nrow(result$cell_high), 2L)
  scores <- as.data.frame(result$scores[result$scores$cohort_label == "SIM_A", ])
  state_expected <- colMeans(t(scale(t(x$demo$expression[x$demo$state_genes, scores$sample_barcode]))))
  expect_equal(scores$state_score, unname(state_expected), tolerance = 1e-10)
  direct <- survival::coxph(survival::Surv(OS.time, OS.status) ~ state_score_z + cell_score_z + age, data = scores)
  got <- result$continuous[result$continuous$cohort_label == "SIM_A", ]
  expect_equal(got$HR, unname(exp(stats::coef(direct)["state_score_z"])), tolerance = 1e-10)
  expect_equal(result$continuous$FDR, stats::p.adjust(result$continuous$wald_p, "BH"))
  expect_equal(sum(scores$cell_group == "high"), 50L)
})

test_that("marker coverage is enforced per patient", {
  x <- two_layer_demo()
  gene_rows <- match(x$demo$state_genes, x$prepared$expression$gene_name)
  sample <- x$demo$clinical$sample_barcode[1]
  data.table::set(x$prepared$expression, i = gene_rows, j = sample, value = NA_real_)
  result <- run_test_layers(x)
  expect_false(sample %in% result$scores$sample_barcode)
  expect_equal(result$cohort_qc[result$cohort_qc$cohort_label == "SIM_A", ]$excluded_n, 1L)
})

test_that("external cell scores require unique IDs and method provenance", {
  x <- two_layer_demo()
  external <- data.frame(sample_barcode = x$demo$clinical$sample_barcode,
    cell_score = as.numeric(x$demo$expression["CCR7", ]))
  expect_error(run_test_layers(x, cell_scores = external), "cell_score_method")
  external <- rbind(external, external[1, ])
  expect_error(run_test_layers(x, cell_scores = external, cell_score_method = "external_test"), "unique")
  external <- external[-nrow(external), ]
  result <- run_test_layers(x, cell_scores = external, cell_score_method = "external_test")
  expect_equal(nrow(result$continuous), 2L)
  expect_true(all(result$cohort_qc$cell_score_method == "external_test"))
})

test_that("perfect separation is logged instead of reported as a finite effect", {
  d <- data.frame(OS.time = seq_len(40), OS.status = c(rep(1, 20), rep(0, 20)),
    group = factor(rep(c("high", "low"), each = 20), levels = c("low", "high")))
  fit <- TCGASigSurvival:::tcgasig_fit_checked(d, "group", "grouphigh")
  expect_false(is.null(fit$error))
})

test_that("batch single-gene analyses run from custom inputs", {
  x <- two_layer_demo()
  result <- run_signature_batch(x$prepared, list(JUN = "JUN", ATF4 = "ATF4"),
    min_signature_genes = 1, min_patients = 40, min_events = 10,
    cutoff_methods = "median", make_km = FALSE, make_summary_plots = FALSE,
    write_files = FALSE, verbose = FALSE)
  expect_equal(nrow(result$results), 4L)
  expect_equal(result$results$batch_FDR, stats::p.adjust(result$results$wald_p, "BH"))
})

test_that("empty result files keep readable headers without spurious warnings", {
  x <- two_layer_demo()
  out <- tempfile("two_layer_empty_")
  on.exit(unlink(out, recursive = TRUE))
  expect_no_warning(result <- run_two_layer_survival(x$prepared, x$demo$state_genes,
    x$demo$cell_genes, min_patients = 1000, cutoff_methods = "median", out_dir = out,
    make_km = FALSE, make_summary_plots = FALSE, verbose = FALSE))
  expect_equal(nrow(result$adjusted), 0L)
  expect_true("HR" %in% names(data.table::fread(result$files$adjusted)))
  expect_equal(nrow(result$skipped), 2L)
})
