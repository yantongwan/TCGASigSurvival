library(TCGASigSurvival)

project_dir <- "my_tcga_project"
# Run download_tcga_data(project_dir) explicitly once before the real analysis.

run_pan_tcga_signature_survival(
  project_dir = project_dir,
  target_gene = "DYRK1A",
  signature_genes = default_treg_signature(),
  signature_name = "Treg"
)

cd8_signature <- c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GZMA", "IFNG")

run_pan_tcga_signature_survival(
  project_dir = project_dir,
  target_gene = "DYRK1A",
  signature_genes = cd8_signature,
  signature_name = "CD8T",
  min_signature_genes = 4
)
