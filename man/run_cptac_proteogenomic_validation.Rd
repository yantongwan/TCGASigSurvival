\name{run_cptac_proteogenomic_validation}
\alias{run_cptac_proteogenomic_validation}
\title{Run or Load CPTAC Proteogenomic Validation}
\usage{
run_cptac_proteogenomic_validation(project_dir = getwd(), proteogenomics_dir = tcgasig_cptac_default_dir(project_dir), 
    target_gene = "STRAP", signature_name = "Treg", signature_genes = NULL, rerun = FALSE, 
    runner = file.path(proteogenomics_dir, "scripts", "run_all.sh"), runner_args = "--all", 
    include_joint_tables = FALSE, verbose = TRUE) 

}
\arguments{
\item{project_dir}{Project directory containing input and output folders.}
\item{proteogenomics_dir}{Directory containing separately generated CPTAC workflow outputs.}
\item{target_gene}{Human target-gene symbol.}
\item{signature_name}{Signature name used in labels and output filenames.}
\item{signature_genes}{Character vector of signature gene symbols.}
\item{rerun}{Whether to execute the external CPTAC workflow before reading its outputs.}
\item{runner}{Path to an externally supplied CPTAC shell runner; not bundled with this package.}
\item{runner_args}{Arguments passed to the external shell runner.}
\item{include_joint_tables}{Whether to load patient-level CPTAC joint tables.}
\item{verbose}{Whether to print progress messages.}
}
\description{
Convenience wrapper for CPTAC proteogenomic validation. By default it reads
existing output tables. With \code{rerun = TRUE}, it first executes the external
\code{PanCancer_STRAP_Proteogenomics} workflow and then loads the updated
tables.
}
\details{
The CPTAC extension is a protein-level validation layer for the TCGA bulk RNA
signature workflow. It keeps the scientific boundary explicit: TCGA and CPTAC
are independent cohorts, and TCGA-CPTAC comparison is cancer/effect-size level
only. Patient-level RNA-protein matching happens only within CPTAC cohorts.
}
\value{
A \code{tcgasig_cptac_proteogenomics} list with an added \code{summary}
component.
}
\examples{
\dontrun{
cptac <- run_cptac_proteogenomic_validation(
  project_dir = "my_tcga_project",
  rerun = FALSE
)
cptac$summary$qc_metrics

cd8 <- run_cptac_proteogenomic_validation(
  project_dir = "my_tcga_project",
  target_gene = "STRAP",
  signature_name = "CD8",
  signature_genes = default_cd8_signature(),
  rerun = FALSE
)
}
}
