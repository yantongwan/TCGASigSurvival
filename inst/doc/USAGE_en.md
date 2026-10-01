# TCGASigSurvival 0.3.0 User Guide

## Installation

```r
install.packages("remotes")
remotes::install_github("yantongwan/TCGASigSurvival", upgrade = "never")
library(TCGASigSurvival)
```

For a local source release, install CRAN dependencies first and then use
`install.packages("TCGASigSurvival_0.3.0.tar.gz", repos = NULL, type = "source")`.
R >= 4.1 is required. Native analyses do not require Python. On Windows the
public input extractor falls back to R gzip streaming; Unix gzip/awk is used
when available. A GitHub release is not a CRAN submission.

## Offline Example

```r
demo <- tcgasig_demo_data()
prepared <- prepare_expression_survival_data(demo$expression, demo$clinical)
summarize_tcga_cohorts(prepared)
two <- run_two_layer_survival(
  prepared, state_signature_genes = demo$state_genes,
  cell_signature_genes = demo$cell_genes,
  state_name = "DemoState", cell_name = "NaiveCD4",
  cutoff_methods = "median",
  out_dir = file.path(tempdir(), "tcgasig_demo")
)
two$continuous
two$adjusted
two$cell_high
two$skipped
```

The demo is deterministic synthetic data; its state genes are illustrative,
not a validated trained-immunity signature. Patient data and unpublished gene
sets are not bundled. Random state is restored after example generation.

## Public TCGA Inputs

```r
project_dir <- "my_tcga_project"
tcga_data_sources()
download_tcga_data(project_dir = project_dir)
```

Four files are required in `data/xena`: `tcga_RSEM_gene_tpm.gz`,
`Survival_SupplementalTable_S1_20171025_xena_sp.tsv`,
`TcgaTargetGTEX_phenotype.txt.gz`, and `gencode.v23.annotation.gtf.gz`.
Installation/loading never downloads them. The downloader writes source URLs,
sizes and local MD5 fingerprints; these are not publisher checksum verification.
Existing complete files are reused; a failed `.part` is retained. Downloads are
not resumed from byte offsets. The original expression download was about
740 MB compressed.

Default cohorts use primary tumors, except SKCM uses primary/metastatic tumors
and provides a separate primary-only sensitivity cohort. Each patient contributes
the alphabetically first eligible sample within each cohort. Overlapping SKCM
cohorts must not be summed as independent patients.

## Target Plus Signature and Signature Alone

```r
joint <- run_pan_tcga_signature_survival(
  target_gene = "ZC3H12C", signature_genes = default_treg_signature(),
  signature_name = "Treg", project_dir = project_dir,
  cutoff_methods = "median"
)
genes <- read_signature_genes("human_signature.csv", gene_column = "Gene")
state <- run_pan_tcga_signature_only_survival(
  signature_genes = genes, signature_name = "MyState",
  project_dir = project_dir, cutoff_methods = "median"
)
joint$survival$results
state$survival$results
```

The historical joint model averages per-gene z scores to obtain a module score,
then sums the standardized module score and standardized target expression.
Its high/low groups are joint-score-defined tumors. The signature-only model
uses the module score alone. Neither isolates cells from bulk RNA.

`uppercase = TRUE` in `read_signature_genes()` normalizes letter case only.
Mouse-to-human orthologue mapping requires a verified external mapping.
Equal-weight mean scoring assumes genes encode a common direction; mixed
up/down signatures need an explicitly defined scoring strategy.

## Two Layers

```r
trained <- read_signature_genes("Gene_1_human.csv", "Gene")
two <- run_pan_tcga_two_layer_survival(
  state_signature_genes = trained,
  cell_signature_genes = default_naive_cd4_signature(),
  state_name = "TrainedImmunity", cell_name = "NaiveCD4",
  project_dir = project_dir, min_state_genes = 20,
  cutoff_methods = c("median", "optimized")
)
```

`continuous` estimates the HR per SD of state score adjusted for continuous
cell-context score. `adjusted` compares state high/low with cell-context
adjustment. `cell_high` repeats state grouping within tumors at or above the
cell-context median. These are independently scored modules, not a summed score.

Defaults are 80 patients/20 deaths in full cohorts and 40 patients/10 deaths
in cell-high subsets. Usable gene thresholds are applied per sample. Review
`scores`, `gene_qc`, `cohort_qc`, and `skipped`. Covariates are optional and must
be supplied explicitly; age, stage and purity are not automatically controlled.

The default Naive CD4 module contains T-cell and naive/central-memory markers
with overlapping expression in other lymphocyte states. It is a proxy context
score, not a uniquely Naive CD4 fraction. A bulk state association remains a
tumor-level association, even after context adjustment.

## Custom Cohorts and External Cell Scores

`expression` must be a named numeric gene-by-sample matrix. `clinical` requires
`sample_barcode`, `patient_id`, `cancer_type`, `OS.time` in days and `OS.status`
(0=censored, 1=dead). Additional columns are preserved.

```r
prepared <- prepare_expression_survival_data(expression, clinical,
  expression_scale = "log", duplicate_patients = "error")
two <- run_two_layer_survival(prepared, state_signature_genes = trained,
  covariates = c("age", "stage"), cutoff_methods = "median")
```

Use `expression_scale = "tpm"` for nonnegative TPM to apply log2(TPM+1).
Raw counts must be normalized appropriately first. Names must be unique;
`duplicate_patients = "first_sample"` is an explicit alternative. Invalid OS
records are excluded and counted in `prepared$input_qc`.

For previously computed deconvolution scores, pass a data frame containing
unique `sample_barcode` and numeric `cell_score`, and set `cell_score_method`:

```r
two <- run_two_layer_survival(prepared, state_signature_genes = trained,
  cell_scores = external_scores, cell_score_method = "documented_external_method",
  cutoff_methods = "median")
```

The package does not run CIBERSORT/xCell or infer cell-specific state expression.
Document the external algorithm/reference/version and exact sample mapping.

## Batch Analysis and Replotting

```r
batch <- run_pan_tcga_signature_batch(
  signature_sets = list(Treg = default_treg_signature(), CD8 = default_cd8_signature()),
  project_dir = project_dir, min_signature_genes = 4, cutoff_methods = "median"
)
gene_batch <- run_signature_batch(prepared,
  signature_sets = list(JUN = "JUN", ATF4 = "ATF4"),
  min_signature_genes = 1, cutoff_methods = "median", make_km = FALSE
)
plot_signature_results(results = joint$survival$results,
  target_gene = "ZC3H12C", signature_name = "Treg")
plot_target_expression_by_signature(survival = joint$survival,
  target_gene = "ZC3H12C", signature_name = "Treg", cutoff_method = "median")
```

Prepared inputs must include every requested gene. Batch public preparation
extracts the union once. `batch_FDR` adjusts across successful signature/cohort
tests within cutoff method and default/sensitivity cohort group. Historical
per-signature `FDR` is retained. Target-expression comparisons use tumor bulk
expression in historical joint-score-defined groups.

## Matched Proteogenomics

```r
matched <- match_proteogenomic_tables(rna_table, protein_table, clinical_table)
validation <- run_matched_proteogenomic_analysis(matched,
  target_gene = "ZC3H12C", signature_name = "Treg")
validation$correlations
validation$survival
validation$qc
validation$skipped
```

RNA and protein tables require `cancer_type`, `patient_id` and respectively
`target_rna`/`target_protein`; clinical requires the same keys and OS fields.
RNA may additionally contain `signature_score`. Inputs should be normalized
measurements from the same patients/study. Exact keys and duplicate rejection
prevent accidental joins but do not prove biological matching: source ID
semantics remain the user's responsibility.

Outputs include Spearman RNA-protein correlations, continuous/median RNA and
protein Cox models, and optional signature-adjusted models. TCGA and CPTAC are
independent validation cohorts; unrelated patients must not be joined.
The historical TCGA joint-score model differs from standalone RNA/protein Cox
effects, so effect definitions must be aligned before claiming replication.

Legacy CPTAC tables can be read with
`run_cptac_proteogenomic_validation(proteogenomics_dir = "pipeline", rerun = FALSE)`.
`rerun = TRUE` requires the separate shell/Python workflow and input data,
which are not bundled. Complete CPTAC downloading, drug screening and docking
are outside this R package's scope.

## Inference and Reproducibility

Prefer continuous models and prespecified median cutoffs. Optimized cutoffs
reuse the outcome for selection: Wald P/BH FDR do not adjust for this search.
They are exploratory. Convergence warnings and infinite estimates are logged
and excluded. Check PH_global_p, sparse events, covariate missingness,
cell-state correlation, marker overlap, purity and clinical confounding.

HR > 1 indicates higher modeled hazard; continuous HR is per 1 SD. Scores are
cohort-relative. New two-layer FDR families are model x cutoff x default_cohort;
legacy per-signature FDR includes all analyzed cohorts within cutoff method.
Keep parameter settings, gene lists, source fingerprints, QC/skipped tables,
package version and session information. See METHODS.md and
`help("run_two_layer_survival")`; cite via `citation("TCGASigSurvival")`.
