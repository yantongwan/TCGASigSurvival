\name{plot_signature_results}
\alias{plot_signature_results}
\alias{plot_signature_forest}
\alias{plot_signature_heatmap}
\title{Plot pan-TCGA signature survival results}
\usage{
plot_signature_results(results = NULL, results_file = NULL, target_gene = NULL, signature_name = NULL, out_dir = file.path(getwd(), 
    "results"), output_prefix = NULL) 

plot_signature_forest(results = NULL, results_file = NULL, target_gene = NULL, signature_name = NULL, out_dir = file.path(getwd(), 
    "results"), output_prefix = NULL, width = 10.5, height = NULL) 

plot_signature_heatmap(results = NULL, results_file = NULL, target_gene = NULL, signature_name = NULL, out_dir = file.path(getwd(), 
    "results"), output_prefix = NULL, width = 7.5, height = NULL) 

}
\description{
Writes forest plot PDF/PNG and heatmap PDF from an existing survival result table.
}
\value{
A list of plot objects and output paths.
}
\arguments{
\item{results}{Survival result table, or a supported result object as described below.}
\item{results_file}{Path to an existing survival result TSV, used when results is NULL.}
\item{target_gene}{Human target-gene symbol.}
\item{signature_name}{Signature name used in labels and output filenames.}
\item{out_dir}{Directory for result tables and figures.}
\item{output_prefix}{Optional basename for generated output files.}
\item{width}{Figure width in inches.}
\item{height}{Figure height in inches.}
}
