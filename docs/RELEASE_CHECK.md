# Release Validation: 0.3.0

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
