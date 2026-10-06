\name{run_immune_abundance_analysis}
\alias{run_immune_abundance_analysis}
\alias{default_immune_signatures}
\alias{score_immune_signatures}
\alias{estimate_immune_fractions}
\alias{run_cell_abundance_survival}
\alias{prepare_tcga_abundance_data}
\alias{plot_cell_abundance_forest}
\alias{plot_cell_abundance_km}
\title{Immune Marker Scores and Estimated Fractions Against Patient Survival}
\description{Standalone abundance analyses without a target gene. Marker scores are not fractions; optional quanTIseq estimates require full-transcriptome linear TPM.}
\usage{
default_immune_signatures()
score_immune_signatures(expression, clinical,
  signatures = default_immune_signatures(), expression_scale = c("log", "tpm"),
  min_coverage = 0.7)
estimate_immune_fractions(tpm, clinical, min_signature_coverage = 0.8)
prepare_tcga_abundance_data(project_dir = getwd(), cancer_type = "PAAD",
  cache_file = NULL)
run_immune_abundance_analysis(expression, clinical,
  expression_scale = c("tpm", "xena_log2_tpm_0.001"),
  signatures = default_immune_signatures(), endpoints = "OS",
  covariate_sets = list(unadjusted = character()), run_deconvolution = TRUE,
  out_dir = file.path(getwd(), "results", "immune_abundance"),
  min_patients = 40L, min_events = 10L, min_group = 10L)
run_cell_abundance_survival(abundance, clinical, endpoints = "OS",
  covariates = character(), adjustment = "unadjusted", min_patients = 40L,
  min_events = 10L, min_group = 10L, min_events_per_parameter = 5,
  out_dir = file.path(getwd(), "results"), write_files = TRUE)
plot_cell_abundance_forest(result, cohort = "PAAD", endpoint = "OS",
  adjustment = "unadjusted", model = "continuous")
plot_cell_abundance_km(result, cell_type, measurement, cohort = "PAAD",
  endpoint = "OS")
}
\arguments{
\item{expression}{Numeric gene-by-sample matrix with unique human gene symbols and sample IDs.}
\item{tpm}{Full-transcriptome linear TPM, finite and nonnegative; not log expression or counts.}
\item{clinical}{Table with sample_barcode, patient_id, cohort_label and endpoint.time/endpoint.status columns. One independent specimen per patient within each cohort.}
\item{signatures}{Uniquely named marker lists; defaults are exploratory modules, not validated cell counters.}
\item{expression_scale}{Explicit scale. The scoring interface accepts log or tpm; the unified interface accepts tpm or xena_log2_tpm_0.001.}
\item{min_coverage}{Minimum usable marker fraction per sample.}
\item{min_signature_coverage}{Minimum exact-symbol effective TIL10 reference coverage.}
\item{project_dir}{Existing local TCGA project with data/xena inputs.}
\item{cancer_type}{TCGA abbreviation.}
\item{cache_file}{Optional new RDS cache; reuse requires matching input fingerprints/settings.}
\item{endpoints,endpoint}{Distinct clinical endpoints. PFI is never substituted for PFS.}
\item{covariate_sets}{Named list of clinical adjustment vectors.}
\item{run_deconvolution}{Run the external quantiseqr method if available. Failure is recorded without discarding marker results.}
\item{out_dir}{New output directory.}
\item{min_patients,min_events,min_group}{Patient/event eligibility gates; min_group applies only to median contrasts.}
\item{abundance}{Long table with sample_barcode, patient_id, cohort_label, cell_type, value, measurement, method, unit, denominator and role.}
\item{covariates}{Existing clinical adjustment columns.}
\item{adjustment}{Adjustment-set label.}
\item{min_events_per_parameter}{Exploratory computational reporting gate, not a guarantee of power.}
\item{write_files}{Write source tables.}
\item{result}{Survival result list, or results data frame for forest plots.}
\item{cohort}{Cohort label.}
\item{model}{Continuous or median contrast.}
\item{cell_type}{Cell label.}
\item{measurement}{marker_score or estimated_fraction.}
}
\value{Default signatures returns named gene lists. Preparation/scoring/estimation/analysis functions return lists of values, QC, settings and/or results. Plot functions return ggplot objects; KM at-risk counts are attached as risk_table.}
\details{Continuous Cox is per fixed reference-cohort SD; fraction effects additionally include HR per 10 percentage points. High means strictly above the reference-cohort median, ties low. Cutoffs are fixed before endpoint/covariate exclusions. FDR spans both measurements and all five primary cells within each cohort/endpoint/model/adjustment; auxiliary outputs are corrected separately. CD4 fraction is non-regulatory CD4 plus Treg; macrophages are model-defined M1 plus M2. Other cells remain in the denominator. No score or fraction is a direct cell count or causal effect.}
\seealso{\code{run_two_layer_survival}, \code{read_tcga_expression_subset}}
\examples{
default_immune_signatures()
}
