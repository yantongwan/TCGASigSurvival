# TCGASigSurvival 0.3.0

- Promote cell-context/state two-layer analysis into general exported functions.
- Add custom expression/clinical preparation and exact external cell-score input.
- Add per-sample marker coverage, cohort/gene QC, separate sensitivity families,
  continuous adjusted models, PH diagnostics and model-failure logs.
- Exclude nonconvergent or infinite Cox estimates from successful inference.
- Add explicit public TCGA downloading, source/fingerprint manifests, signature
  file reading, patient inventory and deterministic synthetic examples.
- Add union-prepared signature/single-gene batches with an additional batch FDR.
- Add native R patient-matched RNA/protein correlation and survival validation.
- Retain historical interfaces and external CPTAC result loading/rerun wrapper;
  external Python workflows and patient data remain separate prerequisites.
- Add Chinese/English guides, method specification, full help pages, portable
  tests, MIT license text and macOS/Linux/Windows GitHub Actions configuration.
- New two-layer outputs may differ from the earlier project-specific script:
  failed models are excluded and sensitivity cohorts have separate BH families.
# TCGASigSurvival 0.4.0

- Add donor-level cell expression preparation and annotated single-cell pseudobulk aggregation.
- Add tumor-versus-healthy/adjacent-normal comparisons with explicit pairing and source distinctions.
- Add fixed-group OS/PFI/PFS/DSS/DFI Cox analyses without substituting PFI for PFS or filtering all patients on OS.
- Add limma-trend/voom differential expression, all-gene ranked fgsea, and leading-edge statistic ridgelines.
- Add opt-in TCGA bulk cell-context preparation, full-gene subset streaming, endpoint grids and KM plots.
- Keep measurement labels and bulk-proxy gates in every new analysis; no bundled patient data or silent downloads.
