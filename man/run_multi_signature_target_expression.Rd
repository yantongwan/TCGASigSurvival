\name{run_multi_signature_target_expression}
\alias{run_multi_signature_target_expression}
\title{Compare target-gene expression across multiple cell signatures}
\usage{
run_multi_signature_target_expression(target_gene, signature_sets, project_dir = getwd(), data_dir = file.path(project_dir, 
    "data", "xena"), out_dir = file.path(project_dir, "results"), output_prefix = NULL, min_signature_genes = 6, 
    min_patients = 80, min_events = 20, include_skcm_primary_only = TRUE, cutoff_method = "optimized", 
    optimized_minprop = 0.1, make_km = FALSE, write_files = TRUE, verbose = TRUE) 

}
\description{
Runs optimized target-signature grouping for a named list of cell signatures and
plots target-gene bulk expression in signature-high and signature-low tumors
relative to each cancer type's overall expression baseline.
}
\value{
A list containing the prepared data, per-signature survival runs, expression
summaries, and combined output paths.
}
\arguments{
\item{target_gene}{Human target-gene symbol.}
\item{signature_sets}{Named list of gene-symbol vectors for multiple signatures.}
\item{project_dir}{Project directory containing input and output folders.}
\item{data_dir}{Directory containing the public TCGA input files.}
\item{out_dir}{Directory for result tables and figures.}
\item{output_prefix}{Optional basename for generated output files.}
\item{min_signature_genes}{Minimum number of usable signature genes within a cohort.}
\item{min_patients}{Minimum number of patients for a cohort model.}
\item{min_events}{Minimum number of OS death events for a cohort model.}
\item{include_skcm_primary_only}{Include the overlapping SKCM primary-only sensitivity cohort.}
\item{cutoff_method}{Grouping method to display, median or optimized.}
\item{optimized_minprop}{Minimum group fraction in optimized cutoff selection.}
\item{make_km}{Whether to generate Kaplan-Meier PDFs.}
\item{write_files}{Whether to write tables and QC reports.}
\item{verbose}{Whether to print progress messages.}
}
