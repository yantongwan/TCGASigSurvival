cell_gene_example <- function(n = 80L, measurement = "sorted") {
  set.seed(10)
  ids <- paste0("D", seq_len(n))
  x <- matrix(rnorm(100L * n), 100L, n, dimnames = list(c("ACLY", paste0("G", seq_len(99))), ids))
  meta <- data.frame(profile_id = ids, patient_id = ids, cancer_type = "SIM",
    cell_type = "Treg", tissue_status = "tumor", age = runif(n, 30, 70),
    cell_score = rnorm(n), OS.time = rexp(n, 0.001) + 1,
    OS.status = rep(c(0, 1), length.out = n), PFI.time = rexp(n, 0.002) + 1,
    PFI.status = rep(c(1, 1, 0), length.out = n))
  prepare_cell_expression_data(x, meta, measurement, provenance = "Synthetic donor data, not a discovery result")
}

test_that("cell data rejects mismatches, pseudoreplicates and unlabelled bulk inference", {
  d <- cell_gene_example()
  expect_error(prepare_cell_expression_data(d$expression, d$metadata[-1, ], provenance = "x"), "match exactly")
  m <- d$metadata; m$patient_id[2] <- m$patient_id[1]
  expect_error(prepare_cell_expression_data(d$expression, m, provenance = "x"), "Multiple profiles")
  proxy <- cell_gene_example(measurement = "bulk_proxy")
  expect_error(run_cell_gene_survival(proxy, write_files = FALSE), "not cell-specific")
  m <- d$metadata; m$tissue_status <- "normal"
  expect_error(prepare_cell_expression_data(d$expression, m, provenance = "x"), "Distinguish")
})

test_that("joint score matches the legacy target plus standardized marker formula", {
  d <- cell_gene_example(measurement = "bulk_proxy")
  r <- run_cell_gene_survival(d, endpoints = "OS", allow_proxy = TRUE, write_files = FALSE,
    score_mode = "joint", min_patients = 30, min_events = 5, min_group = 5)
  expected <- as.numeric(scale(d$expression["ACLY", ])) + as.numeric(scale(d$metadata$cell_score))
  expect_equal(r$groups$analysis_score, expected)
  expect_equal(as.character(r$groups$group), ifelse(expected >= median(expected), "high", "low"))
  independent <- survival::coxph(survival::Surv(OS.time, OS.status) ~ factor(ifelse(expected >= median(expected), "high", "low"), levels = c("low", "high")), d$metadata)
  expect_equal(r$results$HR[r$results$model == "median"], unname(exp(stats::coef(independent))), tolerance = 1e-8)
})

test_that("endpoint eligibility is independent and PFI is never relabelled PFS", {
  d <- cell_gene_example(); d$metadata$OS.time[1:10] <- NA
  r <- run_cell_gene_survival(d, endpoints = c("OS", "PFI", "PFS", "DSS"),
    write_files = FALSE, min_patients = 30, min_events = 5, min_group = 5)
  expect_equal(r$results$n_patients[r$results$endpoint == "PFI"], c(80L, 80L))
  expect_equal(r$results$n_patients[r$results$endpoint == "OS"], c(70L, 70L))
  expect_true(all(r$results$reason[r$results$endpoint %in% c("PFS", "DSS")] == "endpoint_not_measured"))
  g <- r$inputs$group[r$inputs$endpoint == "PFI"]
  expect_equal(g, r$inputs$group[r$inputs$endpoint == "OS"])
  med <- r$results$model == "median"
  expect_equal(r$results$FDR[med], p.adjust(r$results$wald_p[med], "BH"))
  expect_s3_class(plot_cell_gene_survival_grid(r), "ggplot")
})

test_that("adjacent controls are not healthy and pairs are donor-aware", {
  d <- cell_gene_example(40)
  normal <- d$metadata[1:10, ]; normal$profile_id <- paste0(normal$profile_id, "_N")
  normal$tissue_status <- "adjacent_normal"
  x <- cbind(d$expression, d$expression[, 1:10] - 0.3)
  colnames(x) <- c(d$metadata$profile_id, normal$profile_id)
  p <- prepare_cell_expression_data(x, rbind(d$metadata, normal), provenance = "Synthetic pairs")
  healthy <- compare_cell_gene_expression(p, write_files = FALSE)
  expect_equal(healthy$tests$status, "reference_not_available")
  expect_error(compare_cell_gene_expression(p, reference = "adjacent_normal", write_files = FALSE), "paired = TRUE")
  r <- compare_cell_gene_expression(p, reference = "adjacent_normal", paired = TRUE, write_files = FALSE)
  expect_equal(r$tests$tumor_n, 10L)
  expect_equal(r$tests$tumor_minus_reference, 0.3)
  expect_equal(sum(r$values$used_in_test), 20L)
  plot <- plot_cell_gene_expression(r)
  expect_identical(levels(plot$data$tissue_status), c("healthy", "adjacent_normal", "tumor"))
  built <- ggplot2::ggplot_build(plot)
  boxes <- built$data[[1]]
  expect_lt(as.numeric(boxes$x[boxes$fill == "#5F8BB4"]),
    as.numeric(boxes$x[boxes$fill == "#D95F59"]))
  points <- built$data[[2]]
  expect_lt(max(points$x[points$colour == "#5F8BB4"]),
    min(points$x[points$colour == "#D95F59"]))
})

test_that("pseudobulk counts aggregate cells by biological donor and audit cell QC", {
  counts <- matrix(seq_len(48), 4, 12, dimnames = list(paste0("G", 1:4), paste0("C", 1:12)))
  cells <- data.frame(cell_id = colnames(counts), patient_id = rep(c("P1", "P2", "P3"), c(5, 6, 1)),
    cancer_type = "SIM", cell_type = "Treg", tissue_status = "tumor")
  p <- aggregate_cell_pseudobulk(counts, cells, min_cells = 2, provenance = "Synthetic counts")
  expect_equal(ncol(p$expression), 2L)
  expect_equal(p$expression[, 1], rowSums(counts[, 1:5]))
  expect_equal(p$cell_qc$retained, c(TRUE, TRUE, FALSE))
  cells$age <- seq_len(12)
  expect_error(aggregate_cell_pseudobulk(counts, cells, min_cells = 2, provenance = "x"), "Conflicting")
})

test_that("DE and GSEA rank all tested genes and remove circular score components", {
  skip_if_not_installed("limma"); skip_if_not_installed("fgsea")
  d <- cell_gene_example(); d$cell_genes <- c("G1", "G2")
  sets <- list(SET_A = paste0("G", 1:30), SET_B = paste0("G", 31:60))
  old <- .Random.seed
  r <- run_cell_gene_enrichment(d, pathways = sets, score_mode = "joint",
    write_files = FALSE, min_group = 5)
  expect_identical(.Random.seed, old)
  expect_false(any(c("ACLY", "G1", "G2") %in% r$differential$gene))
  expect_equal(nrow(r$differential), 97L)
  expect_equal(length(r$ranks$SIM), 97L)
  expect_equal(r$qc$method, "limma_trend")
  expect_equal(nrow(r$enrichment), 2L)
})

test_that("raw donor counts use TMM voom and preserve input counts", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma"); skip_if_not_installed("fgsea")
  d <- cell_gene_example()
  x <- matrix(rpois(length(d$expression), lambda = 30), nrow(d$expression),
    dimnames = dimnames(d$expression))
  p <- prepare_cell_expression_data(x, d$metadata, "single_cell_pseudobulk", "counts", "Synthetic donor counts")
  r <- run_cell_gene_enrichment(p, pathways = list(A = paste0("G", 1:30)),
    min_group = 5, write_files = FALSE)
  expect_equal(r$qc$method, "TMM_limma_voom")
  expect_identical(p$expression, x)
  expect_true(r$qc$genes_tested > 50)
})

test_that("unified entry runs unaffected modules when healthy reference is absent", {
  d <- cell_gene_example()
  r <- run_cell_gene_analysis(d, out_dir = tempfile("cell_workflow_"), make_plots = FALSE,
    endpoints = c("OS", "PFS"))
  expect_equal(r$module_status$status[1], "reference_or_sample_blocked")
  expect_equal(r$module_status$status[2], "completed_with_qc")
  expect_equal(r$module_status$status[3], "not_requested")
})

test_that("GMT parsing preserves symbols and rejects duplicate pathways", {
  p <- tempfile(fileext = ".gmt")
  writeLines(c("A\tdescription\tACLY\tG1", "B\tdescription\tAcly\tG2"), p)
  expect_equal(read_gmt_pathways(p), list(A = c("ACLY", "G1"), B = c("Acly", "G2")))
  writeLines(c("A\tdescription\tG1", "A\tdescription\tG2"), p)
  expect_error(read_gmt_pathways(p), "Unique")
})
