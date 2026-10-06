if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".",
    "..keep",
    "cptac_cancer",
    "reason",
    "source_file",
    "spearman_FDR",
    "tcga_cancer"
    , "_PATIENT", "_primary_site", "_sample_type", "_study", "..probe_cols", "..samples"
    , "batch_FDR", "cancer type abbreviation", "cancer_plot", "cancer_type"
    , "cancer_type_patient", "cell_FDR", "cell_wald_p", "CI_lower", "CI_upper"
    , "cohort_label", "cutoff_method", "cutoff_value", "default_cohort"
    , "delta_mean_vs_overall", "delta_median_vs_overall", "detailed_category"
    , "direction", "FDR", "fdr_label", "gene_id", "gene_name", "gene_type", "group"
    , "group_events", "group_label", "group_mean", "group_median", "group_n", "has_os", "HR"
    , "is_tcga_barcode", "km_label", "log2_HR", "method_label", "order_hr", "order_value", "OS"
    , "OS.status", "OS.status_patient", "OS.time", "OS.time_patient", "overall_mean"
    , "overall_median", "patient_id", "phenotype_detailed_category", "phenotype_primary_site"
    , "phenotype_sample_type", "phenotype_study", "row_mean_expression", "sample_barcode"
    , "sample_mode", "sample_type", "sample_type_code", "signature_genes_missing"
    , "signature_genes_present", "signature_only_score", "signature_score_mean", "signature_z"
    , "significance", "significant", "survival_patient_id", "target_gene_expr", "target_gene_mean"
    , "target_protein", "target_rna", "target_signature_score", "target_z", "wald_p"
    , "analysis_score", "tissue_status", "endpoint", "log2_hr", "fdr_size", "raw_p"
    , "rank_statistic", "label", "padj", "endpoint_time", "endpoint_status", "months"
    , "cell_type", "measurement", "curve", "survival_estimate", "lower_survival", "upper_survival"
  ))
}
