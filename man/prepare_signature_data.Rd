\name{prepare_signature_data}
\alias{prepare_signature_data}
\title{Prepare selected target and signature expression}
\usage{
prepare_signature_data(target_gene, signature_genes = default_treg_signature(), signature_name = "Treg", 
    project_dir = getwd(), data_dir = file.path(project_dir, "data", "xena"), out_dir = file.path(project_dir, 
        "results"), expression_file = file.path(data_dir, "tcga_RSEM_gene_tpm.gz"), survival_file = file.path(data_dir, 
        "Survival_SupplementalTable_S1_20171025_xena_sp.tsv"), phenotype_file = file.path(data_dir, 
        "TcgaTargetGTEX_phenotype.txt.gz"), gencode_file = file.path(data_dir, "gencode.v23.annotation.gtf.gz"), 
    output_prefix = NULL, write_files = TRUE, verbose = TRUE) 

}
\description{
Extracts selected target and signature gene rows from Xena expression data, joins survival and phenotype tables, and creates tumor cohort definitions.
}
\value{
A prepared-data list.
}
\arguments{
\item{target_gene}{Human target-gene symbol.}
\item{signature_genes}{Character vector of signature gene symbols.}
\item{signature_name}{Signature name used in labels and output filenames.}
\item{project_dir}{Project directory containing input and output folders.}
\item{data_dir}{Directory containing the public TCGA input files.}
\item{out_dir}{Directory for result tables and figures.}
\item{expression_file}{Path to the Xena gene expression matrix.}
\item{survival_file}{Path to the Pan-Cancer Atlas clinical survival table.}
\item{phenotype_file}{Path to Xena sample phenotype metadata.}
\item{gencode_file}{Path to the GENCODE v23 GTF annotation.}
\item{output_prefix}{Optional basename for generated output files.}
\item{write_files}{Whether to write tables and QC reports.}
\item{verbose}{Whether to print progress messages.}
}
