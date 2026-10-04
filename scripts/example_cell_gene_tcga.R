# Real local TCGA example. No patient data or gene sets are bundled in the package.
args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else getwd()
out <- if (length(args) > 1L) args[2] else file.path(project, "results", "ACLY_Treg_Figure_features_20261003")
stage <- if (length(args) > 2L) args[3] else "all"
stopifnot(stage %in% c("preflight", "survival", "enrichment", "all", "figures"))
library(TCGASigSurvival)
stopifnot(packageVersion("TCGASigSurvival") >= "0.4.0")
paths <- c(cache = file.path(project, "results", "panTCGA_ACLY_Treg_prepared_data.rds"),
  clinical = file.path(project, "data", "xena", "Survival_SupplementalTable_S1_20171025_xena_sp.tsv"),
  expression = file.path(project, "data", "xena", "tcga_RSEM_gene_tpm.gz"),
  annotation = file.path(project, "data", "xena", "gencode.v23.annotation.gtf.gz"))
stopifnot(all(file.exists(paths)))
cat("Stage: ", stage, "\n", paste(paths, collapse = "\n"), "\n", sep = "")
if (stage == "preflight") quit(status = 0)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
data <- prepare_tcga_cell_context(paths["cache"], paths["clinical"], cell_type = "Treg")
data.table::fwrite(data$specimen_qc, file.path(out, "specimen_qc.tsv"), sep = "\t")
data.table::fwrite(data$clinical_qc, file.path(out, "clinical_consensus_qc.tsv"), sep = "\t")
data.table::fwrite(data$marker_qc, file.path(out, "Treg_marker_qc.tsv"), sep = "\t")
metadata <- data$metadata
counts <- data.table::as.data.table(metadata)[, .(n_donors = data.table::uniqueN(patient_id)),
  by = .(cancer_type, tissue_status)]
data.table::fwrite(counts, file.path(out, "cohort_donor_counts.tsv"), sep = "\t")
manifest <- data.frame(role = names(paths), path = unname(paths), bytes = file.info(paths)$size,
  md5 = unname(tools::md5sum(paths)))
data.table::fwrite(manifest, file.path(out, "input_manifest.tsv"), sep = "\t")
saveRDS(data, file.path(out, "bulk_Treg_context_data.rds"))

if (stage %in% c("survival", "all")) {
  cat("Running legacy-compatible ACLY + Treg joint score\n")
  main <- run_cell_gene_analysis(data, reference = "adjacent_normal", paired = TRUE,
    endpoints = c("OS", "PFI", "DSS", "DFI", "PFS"), score_mode = "joint",
    min_patients = 80, min_events = 20, min_group = 10, allow_proxy = TRUE,
    out_dir = out, prefix = "ACLY_Treg_joint", make_plots = FALSE)
  healthy <- compare_cell_gene_expression(data, reference = "healthy", allow_proxy = TRUE,
    out_dir = out, prefix = "healthy_control_availability", score_mode = "joint")
  raw_paired <- compare_cell_gene_expression(data, reference = "adjacent_normal", paired = TRUE,
    allow_proxy = TRUE, out_dir = out, prefix = "bulk_ACLY_paired", score_mode = "gene")
  adjusted <- run_cell_gene_survival(data, score_mode = "joint", endpoints = c("OS", "PFI", "DSS", "DFI"),
    covariates = c("age", "sex"), min_patients = 80, min_events = 20, min_group = 10,
    allow_proxy = TRUE, out_dir = out, prefix = "ACLY_Treg_joint_age_sex")
  two <- run_cell_gene_survival(data, score_mode = "gene", context_subset = "high",
    endpoints = c("OS", "PFI", "DSS", "DFI"), covariates = c("age", "sex", "cell_score"),
    min_patients = 40, min_events = 10, min_group = 10, allow_proxy = TRUE,
    out_dir = out, prefix = "Treg_high_bulk_ACLY")
  saveRDS(list(main = main, adjusted = adjusted, two_layer = two,
    healthy = healthy, raw_paired = raw_paired), file.path(out, "survival_comparison.rds"))
  cat("Survival and comparison tables completed\n")
}

if (stage %in% c("enrichment", "all")) {
  # BRCA is preselected as the largest tumor cohort, not by significance or survival direction.
  meta <- metadata[metadata$cancer_type == "BRCA" & metadata$tissue_status == "tumor", ]
  cache <- file.path(out, "BRCA_protein_coding_expression.rds")
  if (file.exists(cache)) full <- readRDS(cache) else {
    cat("Streaming all protein-coding genes for BRCA donors\n")
    full <- read_tcga_expression_subset(paths["expression"], paths["annotation"], meta$profile_id)
    saveRDS(full, cache)
  }
  data.table::fwrite(full$gene_qc, file.path(out, "BRCA_gene_mapping_qc.tsv"), sep = "\t")
  brca <- prepare_cell_expression_data(full$expression, meta, "bulk_proxy", "log", data$provenance)
  brca$cell_genes <- default_treg_signature()
  pathway_cache <- file.path(out, "Hallmark_pathways_local.rds")
  if (file.exists(pathway_cache)) gs <- readRDS(pathway_cache) else {
    stopifnot(requireNamespace("msigdbr", quietly = TRUE))
    cat("Retrieving public human Hallmark gene sets via msigdbr\n")
    msig <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    gs <- list(pathways = split(msig$gene_symbol, msig$gs_name), version = unique(msig$db_version))
    saveRDS(gs, pathway_cache)
  }
  writeLines(c(paste("MSigDB version:", paste(gs$version, collapse = ";")),
    "Source: https://igordot.github.io/msigdbr/ ; https://www.gsea-msigdb.org/",
    "Gene sets used locally, not redistributed with package."), file.path(out, "Hallmark_provenance.txt"))
  cat("BRCA joint high/low DE and Hallmark GSEA\n")
  joint <- run_cell_gene_enrichment(brca, pathways = gs$pathways, score_mode = "joint",
    allow_proxy = TRUE, out_dir = out, prefix = "BRCA_joint")
  sensitivity <- run_cell_gene_enrichment(brca, pathways = gs$pathways, score_mode = "gene",
    context_subset = "high", covariates = c("age", "sex", "cell_score"),
    exclude_genes = c("ACLY", default_treg_signature()), allow_proxy = TRUE,
    out_dir = out, prefix = "BRCA_two_layer")
  saveRDS(list(joint = joint, two_layer = sensitivity), file.path(out, "enrichment_analysis.rds"))
  cat("DE and enrichment tables completed\n")
}

if (stage %in% c("figures", "all")) {
  stopifnot(requireNamespace("svglite", quietly = TRUE))
  cat("Rendering source-backed panels\n")
  s <- readRDS(file.path(out, "survival_comparison.rds"))
  e <- readRDS(file.path(out, "enrichment_analysis.rds"))
  figdir <- file.path(out, "figures")
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  save_plot <- function(p, name, width, height) {
    p <- p + ggplot2::theme(text = ggplot2::element_text(family = "Arial"),
      plot.caption = ggplot2::element_text(size = 8), plot.title = ggplot2::element_text(size = 12),
      plot.subtitle = ggplot2::element_text(size = 9))
    svglite::svglite(file.path(figdir, paste0(name, ".svg")), width = width, height = height)
    print(p)
    grDevices::dev.off()
    grDevices::cairo_pdf(file.path(figdir, paste0(name, ".pdf")), width = width, height = height,
      family = "Arial")
    print(p)
    grDevices::dev.off()
    ggplot2::ggsave(file.path(figdir, paste0(name, ".tiff")), p, width = width, height = height,
      dpi = 600, compression = "lzw")
    ggplot2::ggsave(file.path(figdir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  }
  save_plot(plot_cell_gene_expression(s$main$comparison), "A_joint_paired_adjacent", 11, 4.2)
  save_plot(plot_cell_gene_expression(s$raw_paired), "A_bulk_ACLY_paired_adjacent", 11, 4.2)
  save_plot(plot_cell_gene_enrichment(e$joint, "BRCA"), "B_BRCA_joint_Hallmark_ridges", 11, 8)
  save_plot(plot_cell_gene_enrichment(e$two_layer, "BRCA"), "B_BRCA_two_layer_Hallmark_ridges", 11, 8)
  # PFS stays in source/QC tables as not measured. The displayed fourth outcome is the actual PFI, not PFS.
  display <- s$main$survival
  display$results <- display$results[display$results$endpoint != "PFS", ]
  display$settings$endpoints <- c("OS", "PFI", "DSS", "DFI")
  save_plot(plot_cell_gene_survival_grid(display), "C_joint_survival_grid", 8, 9)
  save_plot(plot_cell_gene_survival_grid(s$adjusted), "C_joint_age_sex_grid", 8, 9)
  save_plot(plot_cell_gene_survival_grid(s$two_layer), "C_two_layer_survival_grid", 8, 9)
  # Fixed illustrative cancers, regardless of significance. No outcome-guided panel selection.
  for (cancer in c("LIHC", "LUAD", "PAAD", "SKCM")) {
    r <- s$main$survival$results
    if (any(r$cohort_label == cancer & r$endpoint == "OS" & r$model == "median" & r$status == "ok")) {
      save_plot(plot_cell_gene_km(s$main$survival, cancer, "OS"), paste0("D_", cancer, "_OS"), 7, 5)
    }
  }
  cat("Figures completed\n")
}
writeLines(c(paste("Completed stage:", stage), capture.output(sessionInfo())), file.path(out, paste0("session_", stage, ".txt")))
