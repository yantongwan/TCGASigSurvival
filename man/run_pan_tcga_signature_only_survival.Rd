\name{run_pan_tcga_signature_only_survival}
\alias{run_pan_tcga_signature_only_survival}
\alias{prepare_signature_only_data}
\alias{run_signature_only_survival}
\alias{plot_signature_only_forest}
\alias{plot_signature_only_heatmap}
\alias{plot_signature_only_results}
\title{Run pan-TCGA signature-only survival analysis}
\usage{
run_pan_tcga_signature_only_survival(signature_genes, signature_name = "Signature", project_dir = getwd(), data_dir = file.path(project_dir, 
    "data", "xena"), out_dir = file.path(project_dir, "results"), output_prefix = NULL, min_signature_genes = 6, 
    min_patients = 80, min_events = 20, include_skcm_primary_only = TRUE, cutoff_methods = c("optimized", 
        "median"), optimized_minprop = 0.1, make_km = TRUE, make_summary_plots = TRUE, write_files = TRUE, 
    verbose = TRUE) 

prepare_signature_only_data(signature_genes, signature_name = "Signature", project_dir = getwd(), data_dir = file.path(project_dir, 
    "data", "xena"), out_dir = file.path(project_dir, "results"), expression_file = file.path(data_dir, 
    "tcga_RSEM_gene_tpm.gz"), survival_file = file.path(data_dir, "Survival_SupplementalTable_S1_20171025_xena_sp.tsv"), 
    phenotype_file = file.path(data_dir, "TcgaTargetGTEX_phenotype.txt.gz"), gencode_file = file.path(data_dir, 
        "gencode.v23.annotation.gtf.gz"), output_prefix = NULL, write_files = TRUE, verbose = TRUE) 

run_signature_only_survival(prepared, signature_genes = NULL, signature_name = NULL, out_dir = file.path(getwd(), 
    "results"), output_prefix = NULL, min_signature_genes = 6, min_patients = 80, min_events = 20, 
    include_skcm_primary_only = TRUE, cutoff_methods = c("optimized", "median"), optimized_minprop = 0.1, 
    make_km = TRUE, write_files = TRUE, verbose = TRUE) 

plot_signature_only_forest(results = NULL, results_file = NULL, signature_name = NULL, out_dir = file.path(getwd(), 
    "results"), output_prefix = NULL, width = 10.5, height = NULL) 

plot_signature_only_heatmap(results = NULL, results_file = NULL, signature_name = NULL, out_dir = file.path(getwd(), 
    "results"), output_prefix = NULL, width = 7.5, height = NULL) 

plot_signature_only_results(results = NULL, results_file = NULL, signature_name = NULL, out_dir = file.path(getwd(), 
    "results"), output_prefix = NULL) 

}
\arguments{
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
\item{write_files}{Whether to write tables and QC reports.}
\item{verbose}{Whether to print progress messages.}
\item{expression_file}{Path to the Xena gene expression matrix.}
\item{survival_file}{Path to the Pan-Cancer Atlas clinical survival table.}
\item{phenotype_file}{Path to Xena sample phenotype metadata.}
\item{gencode_file}{Path to the GENCODE v23 GTF annotation.}
\item{prepared}{Prepared-data list or path to its saved RDS file.}
\item{results}{Survival result table, or a supported result object as described below.}
\item{results_file}{Path to an existing survival result TSV, used when results is NULL.}
\item{width}{Figure width in inches.}
\item{height}{Figure height in inches.}
}
\description{
Runs survival analysis using a cell-state or pathway signature score alone.
For each cancer cohort, each marker gene is z-scored, the signature score is
the mean of available marker-gene z-scores, and patients are stratified into
signature-high and signature-low groups. No target gene expression term is
added. This is appropriate for testing whether a trained-immunity, T cell,
or other signature score itself is associated with overall survival.
}
\details{
The optimized cutoff uses \code{survminer::surv_cutpoint} on the cohort-level
signature score with \code{minprop = optimized_minprop}. The median cutoff is
also available. Cox models are fit as \code{Surv(OS.time, OS.status) ~ group},
where the hazard ratio is high versus low.
}
\value{
A list with \code{prepared}, \code{survival}, and \code{plots}. Output files
are written to \code{out_dir} when \code{write_files = TRUE}.
}
\examples{
\dontrun{
state_genes <- read_signature_genes("state_signature_human.csv", gene_column = "Gene")

res <- run_pan_tcga_signature_only_survival(
  project_dir = "my_tcga_project",
  signature_genes = state_genes,
  signature_name = "StateSignature",
  min_signature_genes = 20
)
}
}
