\name{plot_target_expression_by_signature}
\alias{plot_target_expression_by_signature}
\title{Plot target-gene expression by optimized signature group}
\usage{
plot_target_expression_by_signature(survival = NULL, results = NULL, survival_input = NULL, results_file = NULL, survival_input_file = NULL, 
    target_gene = NULL, signature_name = NULL, cutoff_method = "optimized", default_only = TRUE, 
    out_dir = file.path(getwd(), "results"), output_prefix = NULL, write_files = TRUE, width = 9, 
    delta_height = NULL, boxplot_width = 11, boxplot_height = NULL) 

}
\description{
Uses optimized target-signature score cutoffs to split bulk TCGA tumors into
signature-high and signature-low groups, then compares target-gene bulk
expression in those groups against the overall expression baseline within each
cancer type. This does not estimate true cell-intrinsic gene expression.
}
\value{
A list containing the summary table, plot data, ggplot objects, and output paths.
}
\arguments{
\item{survival}{Returned survival-analysis list containing result and sample-score tables.}
\item{results}{Survival result table, or a supported result object as described below.}
\item{survival_input}{Patient-level survival input and score table.}
\item{results_file}{Path to an existing survival result TSV, used when results is NULL.}
\item{survival_input_file}{Path to a patient-level survival input TSV.}
\item{target_gene}{Human target-gene symbol.}
\item{signature_name}{Signature name used in labels and output filenames.}
\item{cutoff_method}{Grouping method to display, median or optimized.}
\item{default_only}{Whether to display only default cohorts, excluding sensitivity cohorts.}
\item{out_dir}{Directory for result tables and figures.}
\item{output_prefix}{Optional basename for generated output files.}
\item{write_files}{Whether to write tables and QC reports.}
\item{width}{Figure width in inches.}
\item{delta_height}{Height in inches for the mean-expression difference figure.}
\item{boxplot_width}{Width in inches for the target-expression boxplot.}
\item{boxplot_height}{Height in inches for the target-expression boxplot.}
}
