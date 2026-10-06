immune_example <- function(n = 100L) {
  set.seed(812)
  genes <- unique(unlist(default_immune_signatures()))
  x <- matrix(rnorm(length(genes) * n, 4, 1), length(genes), n,
    dimnames = list(genes, paste0("S", seq_len(n))))
  clinical <- data.frame(sample_barcode = colnames(x), patient_id = paste0("P", seq_len(n)),
    cohort_label = "SIM", age = runif(n, 40, 80), sex = rep(c("F", "M"), length.out = n),
    OS.time = rexp(n, 0.001) + 1, OS.status = rep(c(0, 1, 1), length.out = n),
    PFI.time = rexp(n, 0.002) + 1, PFI.status = rep(c(1, 0), length.out = n))
  list(expression = x, clinical = clinical,
    marker = score_immune_signatures(x, clinical))
}

test_that("marker panels reproduce the documented mean-z formula and do not add a gene", {
  d <- immune_example()
  genes <- default_treg_signature()
  expected <- colMeans(t(scale(t(d$expression[genes, , drop = FALSE]))))
  got <- d$marker$values[d$marker$values$cell_type == "Treg", ]
  expect_equal(got$value, unname(expected), tolerance = 1e-12)
  expect_true(all(got$measurement == "marker_score" & got$denominator == "not_a_fraction"))
  expect_equal(nrow(d$marker$values), 500L)
  expect_true(all(d$marker$sample_qc$status == "ok"))
  expect_true(any(nzchar(d$marker$overlap$shared_genes)))
})

test_that("patient duplication and mismatched exact identifiers are blocked", {
  d <- immune_example(); m <- d$clinical; m$patient_id[2] <- m$patient_id[1]
  expect_error(score_immune_signatures(d$expression, m), "one specimen")
  m <- d$clinical; m$sample_barcode[1] <- "UNKNOWN"
  expect_error(score_immune_signatures(d$expression, m), "match expression")
  a <- d$marker$values; a$patient_id[1] <- "WRONG"
  expect_error(run_cell_abundance_survival(a, d$clinical, write_files = FALSE), "identities")
  expect_error(run_cell_abundance_survival(rbind(a, a[1, ]), d$clinical, write_files = FALSE), "Duplicate")
})

test_that("continuous Cox matches an independent patient-level survival fit", {
  d <- immune_example()
  r <- run_cell_abundance_survival(d$marker$values, d$clinical,
    covariates = c("age", "sex"), adjustment = "age_sex", write_files = FALSE)
  g <- r$groups[r$groups$cell_type == "Treg", ]
  m <- d$clinical; m$score <- g$abundance_z[match(m$sample_barcode, g$sample_barcode)]
  fit <- survival::coxph(survival::Surv(OS.time, OS.status) ~ score + age + sex, m)
  row <- r$results[r$results$cell_type == "Treg" & r$results$model == "continuous", ]
  expect_equal(row$status, "ok")
  expect_equal(row$HR, unname(exp(stats::coef(fit)["score"])), tolerance = 1e-10)
  expect_equal(row$wald_p, unname(summary(fit)$coefficients["score", "Pr(>|z|)"]), tolerance = 1e-10)
  expect_true(is.finite(row$PH_predictor_p))
  median_row <- r$results[r$results$cell_type == "Treg" & r$results$model == "median", ]
  m$group <- factor(ifelse(d$marker$values$value[d$marker$values$cell_type == "Treg"] > median_row$cutoff,
    "high", "low"), levels = c("low", "high"))
  median_fit <- survival::coxph(survival::Surv(OS.time, OS.status) ~ group + age + sex, m, x = TRUE)
  expect_equal(median_row$PH_predictor_p, unname(survival::cox.zph(median_fit)$table["group", "p"]))
})

test_that("cutoffs and SDs are independent of endpoint-specific missingness", {
  d <- immune_example(); d$clinical$OS.time[1:15] <- NA
  r <- run_cell_abundance_survival(d$marker$values, d$clinical,
    endpoints = c("OS", "PFI", "PFS"), write_files = FALSE)
  expect_true(all(r$results$n_patients[r$results$endpoint == "OS"] == 85L))
  expect_true(all(r$results$n_patients[r$results$endpoint == "PFI"] == 100L))
  expect_true(all(r$results$reason[r$results$endpoint == "PFS"] == "endpoint_not_measured"))
  i <- r$inputs[r$inputs$cell_type == "Treg", ]
  expect_equal(i$group[i$endpoint == "OS"], i$group[i$endpoint == "PFI"])
  expect_equal(unique(i$reference_sd), stats::sd(d$marker$values$value[d$marker$values$cell_type == "Treg"]))
})

test_that("tied fraction groups do not block an otherwise estimable continuous model", {
  d <- immune_example(); a <- d$marker$values[d$marker$values$cell_type == "Treg", ]
  a$measurement <- "estimated_fraction"; a$method <- "external_validated_method"
  a$unit <- "fraction"; a$denominator <- "all_cells"
  a$value <- c(rep(0, 80), seq(0.02, 0.2, length.out = 20))
  r <- run_cell_abundance_survival(a, d$clinical, min_group = 30, write_files = FALSE)
  expect_equal(r$results$status[r$results$model == "continuous"], "ok")
  expect_equal(r$results$reason[r$results$model == "median"], "insufficient_or_tied_groups")
  expect_equal(sum(r$groups$group == "high"), 20L)
  m <- d$clinical; m$fraction_10pp <- a$value / 0.1
  fit <- survival::coxph(survival::Surv(OS.time, OS.status) ~ fraction_10pp, m)
  expect_equal(r$results$HR_per_10pp[r$results$model == "continuous"], unname(exp(coef(fit))), tolerance = 1e-9)
  a$value <- 0
  r <- run_cell_abundance_survival(a, d$clinical, write_files = FALSE)
  expect_true(all(r$results$status == "not_estimable"))
  a$value[1] <- 1.1
  expect_error(run_cell_abundance_survival(a, d$clinical, write_files = FALSE), "in \\[0,1\\]")
})

test_that("both methods share planned FDR families while auxiliary outputs remain separate", {
  d <- immune_example(); a <- d$marker$values
  f <- a; f$measurement <- "estimated_fraction"; f$method <- "synthetic_fraction"
  f$unit <- "fraction"; f$denominator <- "all_cells"; f$value <- plogis(f$value) / 3
  aux <- f[f$cell_type == "CD4", ]; aux$cell_type <- "CD4_nonTreg"; aux$role <- "auxiliary"
  r <- run_cell_abundance_survival(rbind(a, f, aux), d$clinical, write_files = FALSE)
  for (model in c("continuous", "median")) {
    ii <- r$results$model == model & r$results$role == "primary"
    expect_equal(sum(ii), 10L)
    expect_equal(r$results$FDR[ii], p.adjust(r$results$wald_p[ii], "BH"))
  }
  expect_s3_class(plot_cell_abundance_forest(r, cohort = "SIM"), "ggplot")
  p <- plot_cell_abundance_km(r, "Treg", "marker_score", cohort = "SIM")
  expect_s3_class(p, "ggplot")
  expect_equal(attr(p, "risk_table")$n_risk[attr(p, "risk_table")$months == 0], c(50L, 50L))
  f$value[f$cell_type == "CD4"] <- 0
  failed <- run_cell_abundance_survival(rbind(a, f), d$clinical, write_files = FALSE)
  for (model in c("continuous", "median")) {
    ii <- failed$results$model == model
    expect_true(anyNA(failed$results$wald_p[ii]))
    expect_equal(failed$results$FDR[ii], p.adjust(failed$results$wald_p[ii], "BH", n = 10L))
  }
})

test_that("rank deficiency, missing covariates and insufficient event support remain visible", {
  d <- immune_example(); a <- d$marker$values[d$marker$values$cell_type == "Treg", ]
  d$clinical$copy_score <- as.numeric(scale(a$value))
  r <- run_cell_abundance_survival(a, d$clinical, covariates = "copy_score", write_files = FALSE)
  expect_equal(r$results$reason[r$results$model == "continuous"], "rank_deficient_design")
  d$clinical$age[1:9] <- NA
  r <- run_cell_abundance_survival(a, d$clinical, covariates = "age", write_files = FALSE,
    min_events_per_parameter = 100)
  expect_equal(sum(r$inputs$exclusion == "covariate_missing"), 9L)
  expect_true(all(r$results$reason == "insufficient_events_per_parameter"))
})

test_that("unified workflow inverts the explicit Xena scale and keeps score-only operation", {
  d <- immune_example(); out <- tempfile("immune_analysis_")
  r <- run_immune_abundance_analysis(d$expression, d$clinical,
    expression_scale = "xena_log2_tpm_0.001", run_deconvolution = FALSE, out_dir = out)
  expect_equal(r$module_status$status, c("completed", "not_requested"))
  expect_true(file.exists(file.path(out, "survival_results.tsv")))
  expect_equal(nrow(r$survival$results), 10L)
  bad <- d$expression; bad[1] <- -20
  expect_error(run_immune_abundance_analysis(bad, d$clinical,
    expression_scale = "xena_log2_tpm_0.001", run_deconvolution = FALSE, out_dir = out), "inconsistent")
})

test_that("quanTIseq wrapper matches the external implementation and retains the Other denominator", {
  skip_if_not_installed("quantiseqr")
  e <- new.env(); utils::data("dataset_racle", package = "quantiseqr", envir = e)
  x <- e$dataset_racle$expr_mat
  m <- data.frame(sample_barcode = colnames(x), patient_id = colnames(x), cohort_label = "SIM")
  r <- estimate_immune_fractions(x, m)
  direct <- suppressMessages(quantiseqr::run_quantiseq(x, is_tumordata = TRUE, return_se = FALSE))
  expect_equal(r$fractions, as.data.frame(direct))
  expect_equal(r$sample_qc$sum_fractions, rep(1, ncol(x)), tolerance = 1e-8)
  cd4 <- r$values$value[r$values$cell_type == "CD4"]
  expect_equal(cd4, direct$T.cells.CD4 + direct$Tregs)
  expect_equal(r$values$value[r$values$cell_type == "Macrophages"], direct$Macrophages.M1 + direct$Macrophages.M2)
  expect_true(all(r$values$denominator == "all_modeled_cells_including_Other"))
  expect_error(estimate_immune_fractions(log2(x + 1), m), "unsafe")
})

test_that("all-biotype loading is additive and keeps the historical protein-coding default", {
  gtf <- tempfile(fileext = ".gtf"); file <- tempfile(fileext = ".tsv")
  writeLines(c('1\tx\tgene\t1\t2\t.\t+\t.\tgene_id "E1"; gene_name "G1"; gene_type "protein_coding";',
    '1\tx\tgene\t3\t4\t.\t+\t.\tgene_id "E2"; gene_name "G2"; gene_type "lincRNA";'), gtf)
  writeLines(c("sample\tS1\tS2", "E1\t1\t2", "E2\t3\t4"), file)
  all <- read_tcga_expression_subset(file, gtf, c("S1", "S2"), gene_type = NULL)
  old <- read_tcga_expression_subset(file, gtf, c("S1", "S2"))
  expect_equal(rownames(all$expression), c("G1", "G2"))
  expect_equal(rownames(old$expression), "G1")
})
