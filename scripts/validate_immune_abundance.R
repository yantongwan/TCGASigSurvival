# Independent numerical checks against the saved patient-level analysis.
library(TCGASigSurvival)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Provide the completed output directory.")
out <- normalizePath(args[1])
res <- readRDS(file.path(out, "analysis/immune_abundance_analysis.rds"))
d <- readRDS(file.path(out, "PAAD_full_transcriptome_input.rds"))
r <- res$survival$results
checks <- list()
record <- function(name, ok, detail = "") {
  checks[[name]] <<- data.frame(check = name, pass = isTRUE(ok), detail = detail)
  if (!isTRUE(ok)) stop("Validation failed: ", name, " ", detail)
}
record("both_measurements_completed", all(res$module_status$status == "completed"))
record("one_sample_per_patient", !anyDuplicated(d$clinical$patient_id))
record("fraction_composition", all(abs(rowSums(res$fraction$fractions[, -1]) - 1) < 1e-6))
record("input_fingerprints_unchanged", identical(unname(tools::md5sum(d$manifest$path)), d$manifest$md5))
record("method_patient_pairing", all(data.table::fread(file.path(out, "method_pairing_qc.tsv"))$same_patients))
for (cohort in unique(res$marker$values$cohort_label)) for (cell in names(default_immune_signatures())) {
  m <- res$marker$values
  m <- m[m$cohort_label == cohort & m$cell_type == cell, ]
  qc <- res$marker$gene_qc
  genes <- qc$gene[qc$cohort_label == cohort & qc$cell_type == cell & qc$status == "usable"]
  x <- d$expression[genes, m$sample_barcode, drop = FALSE]
  z <- t(scale(t(x)))
  record(paste("direct_marker_formula", cohort, cell),
    isTRUE(all.equal(as.numeric(colMeans(z)), m$value, tolerance = 1e-10)))
}
fam <- interaction(r$cohort_label, r$endpoint, r$model, r$adjustment, r$role, drop = TRUE)
for (f in unique(fam)) {
  ii <- which(fam == f)
  record(paste("FDR", f), isTRUE(all.equal(r$FDR[ii], stats::p.adjust(r$wald_p[ii], "BH", n = length(ii)))))
}
for (i in which(r$status == "ok")) {
  row <- r[i, ]; x <- res$survival$inputs
  j <- x$cohort_label == row$cohort_label & x$measurement == row$measurement &
    x$cell_type == row$cell_type & x$endpoint == row$endpoint & x$adjustment == row$adjustment & x$exclusion == "included"
  x <- x[j, ]
  record(paste("patient_event_counts", i), nrow(x) == row$n_patients && sum(x$endpoint_status) == row$n_events)
  record(paste("fixed_groups", i), all(as.character(x$group) == ifelse(x$value > row$cutoff, "high", "low")))
  x$group <- factor(x$group, levels = c("low", "high"))
  predictor <- if (row$model == "continuous") "abundance_z" else "group"
  term <- if (row$model == "continuous") "abundance_z" else "grouphigh"
  cov <- if (nzchar(row$covariates_used)) strsplit(row$covariates_used, ";", fixed = TRUE)[[1]] else character()
  formula <- stats::as.formula(paste("survival::Surv(endpoint_time, endpoint_status) ~", paste(c(predictor, cov), collapse = " + ")))
  fit <- survival::coxph(formula, data = x)
  s <- summary(fit)
  expected <- unname(c(s$conf.int[term, c("exp(coef)", "lower .95", "upper .95")], s$coefficients[term, "Pr(>|z|)"]))
  got <- unname(unlist(row[, c("HR", "CI_lower", "CI_upper", "wald_p")]))
  record(paste("direct_cox", i), isTRUE(all.equal(expected, got, tolerance = 1e-8)))
}
figs <- data.table::fread(file.path(out, "figure_manifest.tsv"))
for (name in figs$figure[figs$status == "completed"]) {
  record(paste("all_exports", name), all(file.exists(file.path(out, "figures", paste0(name, c(".pdf", ".svg", ".png", ".tiff"))))))
  if (!grepl("_KM_", name, fixed = TRUE)) next
  risk <- data.table::fread(file.path(out, "figures", paste0(name, "_risk_table.tsv")))
  method <- if (grepl("marker_score$", name)) "marker_score" else "estimated_fraction"
  cell <- sub("^PAAD_OS_KM_", "", name); cell <- sub(paste0("_", method, "$"), "", cell)
  x <- res$survival$inputs
  x <- x[x$cohort_label == "PAAD" & x$cell_type == cell & x$measurement == method &
    x$endpoint == "OS" & x$adjustment == "unadjusted" & x$exclusion == "included", ]
  expected <- vapply(seq_len(nrow(risk)), function(i)
    sum(x$group == risk$group[i] & x$endpoint_time >= risk$months[i] * (365.25 / 12)), integer(1))
  record(paste("risk_counts", name), identical(expected, as.integer(risk$n_risk)))
}
data.table::fwrite(do.call(rbind, checks), file.path(out, "NUMERICAL_QA.tsv"), sep = "\t")
cat(length(checks), "numerical checks passed.\n")
