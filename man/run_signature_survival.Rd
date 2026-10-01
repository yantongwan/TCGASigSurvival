\name{run_signature_survival}
\alias{run_signature_survival}
\title{Run survival analysis from prepared data}
\usage{
run_signature_survival(prepared, target_gene = NULL, signature_genes = NULL, signature_name = NULL, out_dir = file.path(getwd(), 
    "results"), output_prefix = NULL, min_signature_genes = 6, min_patients = 80, min_events = 20, 
    include_skcm_primary_only = TRUE, cutoff_methods = c("optimized", "median"), optimized_minprop = 0.1, 
    make_km = TRUE, write_files = TRUE, verbose = TRUE) 

}
\description{
Computes within-cancer target z-score, signature score, target-signature score, optimized and median cutoffs, KM plots, and Cox model summaries.
}
\value{
A list containing result, survival input, skipped cohort, and output path information.
}
\arguments{
\item{prepared}{Prepared-data list or path to its saved RDS file.}
\item{target_gene}{Human target-gene symbol.}
\item{signature_genes}{Character vector of signature gene symbols.}
\item{signature_name}{Signature name used in labels and output filenames.}
\item{out_dir}{Directory for result tables and figures.}
\item{output_prefix}{Optional basename for generated output files.}
\item{min_signature_genes}{Minimum number of usable signature genes within a cohort.}
\item{min_patients}{Minimum number of patients for a cohort model.}
\item{min_events}{Minimum number of OS death events for a cohort model.}
\item{include_skcm_primary_only}{Include the overlapping SKCM primary-only sensitivity cohort.}
\item{cutoff_methods}{Median and/or optimized grouping methods. Optimized cutoffs are exploratory.}
\item{optimized_minprop}{Minimum group fraction in optimized cutoff selection.}
\item{make_km}{Whether to generate Kaplan-Meier PDFs.}
\item{write_files}{Whether to write tables and QC reports.}
\item{verbose}{Whether to print progress messages.}
}
