# Release Validation: 0.4.0

Validated locally on 2026-10-03 with R 4.4.3 on macOS arm64.

- Built and installed `TCGASigSurvival_0.4.0.tar.gz` in an isolated library.
- Final `R CMD check --no-manual --no-build-vignettes`: **Status: OK**,
  zero errors, warnings and notes. Network repository-index checks were unavailable;
  installed dependencies were checked locally.
- All legacy and new automated tests passed. Nine new test blocks cover donor
  identity, pseudobulk aggregation, legacy joint-score equivalence, endpoint-specific
  eligibility, PFI/PFS distinction, paired reference provenance, limma/voom,
  all-gene GSEA, score-component exclusion, RNG and partial-module execution.
- All five executable README tutorial blocks passed against 0.4.0.
- The real ACLY/Treg example completed preflight, survival, BRCA enrichment and
  figure stages: 32 selected cohorts, 9528 unique tumor patients, 20 eligible
  paired tumor/adjacent-normal cohorts and 676 donor pairs.
- The real example uses the legacy bulk joint score and a two-layer sensitivity;
  it does not measure purified-Treg ACLY. Healthy Treg controls and independent
  PFS endpoints were absent and are explicitly recorded, not imputed.
- Eleven single-panel figures have PDF/SVG/600-dpi TIFF/PNG exports. Source
  preflight and rendered-glyph geometry QA records accompany private local outputs;
  raw Cairo text-audit limitations and intentional KM censor marks are documented.
- Source archives include scripts and tutorials but not patient data, private
  gene sets or MSigDB definitions. The installable R tarball follows `.Rbuildignore`;
  command-line scripts are in the full source ZIP/checkout.
- The historical `v0.3.0` release remains unchanged. Local validation does not
  establish remote three-platform CI success; inspect GitHub Actions for the
  actual status of each published commit.

Real-data entry point (inputs must already be present):

```sh
Rscript scripts/example_cell_gene_tcga.R /path/to/tcga_project /path/to/new_output preflight
Rscript scripts/example_cell_gene_tcga.R /path/to/tcga_project /path/to/new_output all
```

See `inst/doc/CELL_GENE_zh.md` for the complete feature tutorial.

## Historical Validation: 0.3.0

Validated locally on 2026-10-01 with R 4.4.3 on macOS arm64.

- Built `TCGASigSurvival_0.3.0.tar.gz` and installed it into an isolated R library.
- `R CMD check ... --no-manual --no-build-vignettes`: **Status: OK**, zero
  package errors, warnings or notes. Online CRAN/Bioconductor index checks were
  unavailable in this local sandbox; all required dependencies were already installed.
- All automated tests and documented R examples passed, including direct
  independent Cox estimate comparisons, marker coverage, convergence failure,
  exact cancer/patient matching and readable empty output tables.
- All five executable README tutorial blocks passed in order, including
  state-only, two-layer, age-adjusted, batch and matched RNA/protein examples.
- Installed-package synthetic end-to-end validation exercised median/optimized
  two-layer models, KM plots, forests, target-expression comparisons, batch
  single-gene models and native proteogenomic analyses.
- A real public TCGA wrapper smoke test used canonical Naive-CD4-context/CD8
  modules: 32 default cancer cohorts, 21 successful continuous models,
  21 adjusted median models and 20 cell-high-subset median models.
- The real smoke test validates software execution; it is not a trained-immunity
  discovery result. The project's unpublished signature and patient-level
  outputs are not included in the release.

GitHub Actions is configured for Linux, macOS and Windows; local macOS success
does not by itself establish that remote cross-platform checks have passed.

Reproduce synthetic checks from a source checkout after installing the package:

```sh
Rscript scripts/validate_release.R /path/to/validation_output
Rscript scripts/validate_readme.R README.md /path/to/tutorial_output
```

Reproduce the real TCGA smoke test after preparing the four public input files:

```sh
Rscript scripts/validate_real_tcga.R /path/to/tcga_project
```
