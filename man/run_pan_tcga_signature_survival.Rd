\name{run_pan_tcga_signature_survival}
\alias{run_pan_tcga_signature_survival}
\title{Run pan-TCGA target-gene/signature survival analysis}
\usage{
run_pan_tcga_signature_survival(target_gene, signature_genes = default_treg_signature(), signature_name = "Treg", 
    project_dir = getwd(), data_dir = file.path(project_dir, "data", "xena"), out_dir = file.path(project_dir, 
        "results"), output_prefix = NULL, min_signature_genes = 6, min_patients = 80, min_events = 20, 
    include_skcm_primary_only = TRUE, cutoff_methods = c("optimized", "median"), optimized_minprop = 0.1, 
    make_km = TRUE, make_summary_plots = TRUE, make_expression_plot = TRUE, write_files = TRUE, 
    verbose = TRUE) 

}
\description{
Prepares selected target-gene and signature-gene expression from UCSC Xena pan-TCGA data, runs per-cancer survival analyses, and writes KM plots, forest plot, heatmap, result table, and QC report.
}
\value{
A list containing prepared data, survival results, and plot paths.
}
\arguments{
\item{target_gene}{Human target-gene symbol.}
\item{signature_genes}{Character vector of signature gene symbols.}
\item{signature_name}{Signature name used in labels and output filenames.}
\item{project_dir}{Project directory containing input and output folders.}
\item{data_dir}{Directory containing the public TCGA input files.}
\item{out_dir}{Directory for result tables and figures.}
\item{output_prefix}{Optional basename for generated output files.}
\item{min_signature_genes}{Minimum number of usable signature genes within a cohort.}
\item{min_patients}{Minimum number of patients for a cohort model.}
\item{min_events}{Minimum number of OS death events for a cohort model.}
\item{include_skcm_primary_only}{Include the overlapping SKCM primary-only sensitivity cohort.}
\item{cutoff_methods}{Median and/or optimized grouping methods. Optimized cutoffs are exploratory.}
\item{optimized_minprop}{Minimum group fraction in optimized cutoff selection.}
\item{make_km}{Whether to generate Kaplan-Meier PDFs.}
\item{make_summary_plots}{Whether to generate forest plots and heatmaps.}
\item{make_expression_plot}{Whether to generate the target-expression comparison for optimized groups.}
\item{write_files}{Whether to write tables and QC reports.}
\item{verbose}{Whether to print progress messages.}
}
