\name{summarize_cptac_proteogenomic_results}
\alias{summarize_cptac_proteogenomic_results}
\title{Summarize CPTAC Proteogenomic Validation Results}
\usage{
summarize_cptac_proteogenomic_results(results, fdr_cutoff = 0.05) 

}
\arguments{
\item{results}{Survival result table, or a supported result object as described below.}
\item{fdr_cutoff}{BH false-discovery-rate threshold used for summaries.}
}
\description{
Creates compact QC, RNA-protein concordance, survival-status, and TCGA-CPTAC
replication summary tables from loaded CPTAC proteogenomic validation results.
}
\value{
A list containing \code{qc_metrics}, \code{cohort_inventory},
\code{rna_protein_concordance}, \code{survival_status},
\code{tcga_cptac_replication}, and \code{interpretation_notes}.
}
\examples{
\dontrun{
cptac <- load_cptac_proteogenomic_results(
  project_dir = "my_tcga_project"
)
summ <- summarize_cptac_proteogenomic_results(cptac)
summ$qc_metrics
}
}
