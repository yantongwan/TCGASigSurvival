# Reproducible local PAAD workflow. No raw downloads or upstream sequencing runs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && (length(args) %% 2L || any(!grepl("^--", args[seq(1, length(args), 2)])))) {
  stop("Use --project-dir PATH --out-dir PATH --stage all|prepare|analyze|figures|report")
}
options <- if (length(args)) setNames(as.list(args[seq(2, length(args), 2)]),
  sub("^--", "", args[seq(1, length(args), 2)])) else list()
if (!is.null(options[["library-dir"]])) .libPaths(c(options[["library-dir"]], .libPaths()))
library(TCGASigSurvival)
if (utils::packageVersion("TCGASigSurvival") < "0.5.0") stop("Install TCGASigSurvival >= 0.5.0.")
project <- if (is.null(options[["project-dir"]])) getwd() else options[["project-dir"]]
out <- if (is.null(options[["out-dir"]])) file.path(project, "results", "PAAD_Immune_Abundance") else options[["out-dir"]]
stage <- if (is.null(options$stage)) "all" else options$stage
if (!stage %in% c("all", "prepare", "analyze", "figures", "report")) stop("Unknown stage.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
cache <- file.path(out, "PAAD_full_transcriptome_input.rds")
analysis_path <- file.path(out, "analysis", "immune_abundance_analysis.rds")
write_table <- function(d, name) data.table::fwrite(d, file.path(out, name), sep = "\t", na = "NA")
ledger <- function(step, state, detail = "") {
  data.table::fwrite(data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    stage = step, status = state, detail = detail), file.path(out, "run_ledger.tsv"),
    sep = "\t", append = file.exists(file.path(out, "run_ledger.tsv")))
}
ledger(stage, "started", paste("TCGASigSurvival", utils::packageVersion("TCGASigSurvival")))
if (stage %in% c("all", "prepare")) {
  d <- prepare_tcga_abundance_data(project, "PAAD", cache)
  for (name in c("manifest", "specimen_qc", "gene_qc", "transform_qc", "clinical"))
    write_table(d[[name]], paste0("input_", name, ".tsv"))
  metadata <- file.path(out, "tcga_RSEM_gene_tpm.metadata.json")
  if (!file.exists(metadata)) utils::download.file(
    "https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.json", metadata, mode = "wb", quiet = TRUE)
  metadata_text <- paste(readLines(metadata, warn = FALSE), collapse = "\n")
  if (!grepl('"unit"\\s*:\\s*"log2\\(tpm\\+0.001\\)"', metadata_text, perl = TRUE))
    stop("Official expression metadata does not confirm expected scale.")
  ledger("prepare", "completed", paste(nrow(d$clinical), "independent patients"))
}
if (stage %in% c("all", "analyze")) {
  if (file.exists(analysis_path)) stop("Analysis already exists. Use figures/report stages or a new versioned output directory.")
  d <- prepare_tcga_abundance_data(project, "PAAD", cache)
  strict <- d$clinical[d$clinical$histological_type == "Pancreas-Adenocarcinoma Ductal Type" &
    !is.na(d$clinical$histological_type), , drop = FALSE]
  strict$cohort_label <- "PAAD_Ductal"
  clinical <- rbind(d$clinical, strict)
  write_table(clinical, "analysis_clinical.tsv")
  hist <- as.data.frame(table(clinical$cohort_label, clinical$histological_type, useNA = "ifany"))
  names(hist) <- c("cohort_label", "histological_type", "patients")
  write_table(hist[hist$patients > 0, ], "histology_counts.tsv")
  stages <- as.data.frame(table(clinical$cohort_label, clinical$stage_group, useNA = "ifany"))
  names(stages) <- c("cohort_label", "stage_group", "patients")
  write_table(stages, "stage_counts.tsv")
  signatures <- default_immune_signatures()
  write_table(do.call(rbind, lapply(names(signatures), function(cell)
    data.frame(cell_type = cell, gene = signatures[[cell]]))), "marker_definitions.tsv")
  res <- run_immune_abundance_analysis(d$expression, clinical,
    expression_scale = d$expression_scale, signatures = signatures,
    endpoints = c("OS", "PFI", "DSS", "DFI", "PFS"),
    covariate_sets = list(unadjusted = character(), age_sex = c("age", "sex"),
      age_sex_stage = c("age", "sex", "stage_group")), out_dir = file.path(out, "analysis"))
  saveRDS(list(package_version = as.character(utils::packageVersion("TCGASigSurvival")),
    input_manifest = d$manifest, settings = res$settings, fraction_settings = res$fraction$settings),
    file.path(out, "run_settings.rds"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
  ledger("analyze", "completed", paste(nrow(res$survival$results), "explicit result rows; see module_status.tsv"))
}
if (stage %in% c("all", "figures", "report")) {
  if (!file.exists(analysis_path)) stop("Run analyze first.")
  res <- readRDS(analysis_path)
}
if (stage %in% c("all", "figures")) {
  figure_dir <- file.path(out, "figures"); dir.create(figure_dir, showWarnings = FALSE)
  figure_manifest <- list()
  export_plot <- function(p, name, width = 7.5, height = 5.2) {
    path <- file.path(figure_dir, name)
    ggplot2::ggsave(paste0(path, ".pdf"), p, width = width, height = height,
      device = grDevices::cairo_pdf, family = "Arial")
    ggplot2::ggsave(paste0(path, ".svg"), p, width = width, height = height, device = svglite::svglite)
    ggplot2::ggsave(paste0(path, ".png"), p, width = width, height = height, dpi = 300, device = ragg::agg_png)
    ggplot2::ggsave(paste0(path, ".tiff"), p, width = width, height = height, dpi = 600,
      device = ragg::agg_tiff, compression = "lzw")
    figure_manifest[[name]] <<- data.frame(figure = name, status = "completed", reason = "",
      width_in = width, height_in = height, source = "analysis/survival_results.tsv;analysis/survival_inputs.tsv",
      pdf_md5 = unname(tools::md5sum(paste0(path, ".pdf"))))
  }
  for (cohort in c("PAAD", "PAAD_Ductal")) for (adj in c("unadjusted", "age_sex", "age_sex_stage")) {
    export_plot(plot_cell_abundance_forest(res$survival, cohort = cohort, adjustment = adj),
      paste0(cohort, "_OS_continuous_", adj))
  }
  for (cell in names(default_immune_signatures())) for (method in c("marker_score", "estimated_fraction")) {
    name <- paste0("PAAD_OS_KM_", cell, "_", method)
    p <- tryCatch(plot_cell_abundance_km(res$survival, cell, method), error = function(e) e)
    if (inherits(p, "error")) {
      figure_manifest[[name]] <- data.frame(figure = name, status = "not_estimable", reason = conditionMessage(p),
        width_in = NA, height_in = NA, source = "analysis/survival_results.tsv", pdf_md5 = NA)
    } else {
      export_plot(p, name)
      data.table::fwrite(attr(p, "risk_table"), file.path(figure_dir, paste0(name, "_risk_table.tsv")), sep = "\t")
    }
  }
  write_table(do.call(rbind, figure_manifest), "figure_manifest.tsv")
  ledger("figures", "completed", paste(length(figure_manifest), "requested figures with explicit status"))
}
if (stage %in% c("all", "report")) {
  r <- res$survival$results
  counts <- unique(r[, c("cohort_label", "measurement", "endpoint", "adjustment", "n_patients", "n_events")])
  write_table(counts, "endpoint_patient_counts.tsv")
  write_table(as.data.frame(table(r$cohort_label, r$endpoint, r$status, r$reason)), "model_status_counts.tsv")
  agreement <- list(); pairing <- list()
  for (cohort in unique(r$cohort_label)) for (cell in names(default_immune_signatures())) {
    a <- res$abundance
    x <- a[a$cohort_label == cohort & a$cell_type == cell & a$measurement == "marker_score", ]
    y <- a[a$cohort_label == cohort & a$cell_type == cell & a$measurement == "estimated_fraction", ]
    m <- merge(x[, c("sample_barcode", "value")], y[, c("sample_barcode", "value")], by = "sample_barcode")
    test <- if (nrow(m) > 2 && stats::sd(m$value.x) > 0 && stats::sd(m$value.y) > 0)
      stats::cor.test(m$value.x, m$value.y, method = "spearman", exact = FALSE) else NULL
    agreement[[paste(cohort, cell)]] <- data.frame(cohort_label = cohort, cell_type = cell, n = nrow(m),
      rho = if (is.null(test)) NA else unname(test$estimate), p = if (is.null(test)) NA else test$p.value)
    for (ep in unique(r$endpoint)) for (adj in unique(r$adjustment)) {
      z <- res$survival$inputs
      z <- z[z$cohort_label == cohort & z$cell_type == cell & z$endpoint == ep &
        z$adjustment == adj & z$exclusion == "included", ]
      sx <- z$patient_id[z$measurement == "marker_score"]; sy <- z$patient_id[z$measurement == "estimated_fraction"]
      pairing[[paste(cohort, cell, ep, adj)]] <- data.frame(cohort_label = cohort, cell_type = cell,
        endpoint = ep, adjustment = adj, marker_n = length(sx), fraction_n = length(sy), same_patients = setequal(sx, sy))
    }
  }
  agreement <- do.call(rbind, agreement); agreement$FDR <- NA_real_
  for (cohort in unique(agreement$cohort_label)) {
    j <- agreement$cohort_label == cohort; agreement$FDR[j] <- stats::p.adjust(agreement$p[j], "BH")
  }
  write_table(agreement, "method_agreement.tsv"); write_table(do.call(rbind, pairing), "method_pairing_qc.tsv")
  main <- r[r$cohort_label == "PAAD" & r$endpoint == "OS" & r$role == "primary" & r$adjustment == "unadjusted", ]
  write_table(main, "main_OS_results.tsv")
  fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = 3, format = "g"))
  table_lines <- function(z) c("|细胞|测量|模型|患者/事件|HR (95%CI)|P|FDR|状态|", "|---|---|---|---|---|---|---|---|",
    vapply(seq_len(nrow(z)), function(i) paste0("|", z$cell_type[i], "|", z$measurement[i], "|", z$model[i], "|",
      z$n_patients[i], "/", z$n_events[i], "|", fmt(z$HR[i]), " (", fmt(z$CI_lower[i]), ", ", fmt(z$CI_upper[i]), ")|",
      fmt(z$wald_p[i]), "|", fmt(z$FDR[i]), "|", z$status[i], if (nzchar(z$reason[i])) paste0(":", z$reason[i]) else "", "|"), character(1)))
  ph <- r[r$status == "ok" & !is.na(r$PH_predictor_p) & r$PH_predictor_p < 0.05, ]
  write_table(ph, "PH_predictor_flags.tsv")
  files <- c("input_manifest.tsv", "analysis/survival_results.tsv", "analysis/survival_inputs.tsv",
    "analysis/marker_gene_qc.tsv", "analysis/full_fractions.tsv", "analysis/fixed_groups.tsv",
    "method_agreement.tsv", "method_pairing_qc.tsv", "endpoint_patient_counts.tsv", "figure_manifest.tsv")
  lines <- c("# PAAD 五类免疫细胞与生存：实际运行报告", "", paste0("运行版本：TCGASigSurvival ",
    utils::packageVersion("TCGASigSurvival"), "；输出根目录：`", out, "`。"), "",
    "## 样本和终点", "",
    paste0("主队列", max(r$initial_n[r$cohort_label == "PAAD"]), "位独立患者；严格导管型标签子集",
      max(r$initial_n[r$cohort_label == "PAAD_Ductal"]), "位。每患者仅一个原发肿瘤标本。两队列重叠，不是外部验证。"), "",
    paste(utils::capture.output(print(counts[counts$measurement == "marker_score" & counts$adjustment == "unadjusted", ], row.names = FALSE)), collapse = "\n"), "",
    "OS为主要终点，PFI/DSS/DFI为补充探索；本地CDR没有PFS，不以PFI替代。缺失或非正随访时间逐人记录。", "",
    "## 两种测量", "",
    "特征评分为癌种内逐基因z-score均值，不是比例；CD8含细胞毒性/NK共享基因，CD4含通用T细胞基因，髓系标记受炎症/极化影响。", "",
    "比例来自标准quanTIseq/TIL10：将官方log2(TPM+0.001)逆变换为线性全转录组TPM；微小负数只按舍入误差裁零。使用肿瘤设置、mRNA校正及标准去噪，不按结果调参。", "",
    "CD4主比例=非调节型CD4+Treg；巨噬细胞=M1+M2。Other保留在分母。五类比例不可直接相加（CD4包含Treg）。辅助结果另列非Treg CD4和M1/M2。", "",
    "## 主要OS结果：连续Cox", "",
    "HR为测量每升高一个参考队列SD的风险比。HR>1与较高风险相关，HR<1与较低风险相关。95%CI和FDR全部保留。", "",
    table_lines(main[main$model == "continuous", ]), "", "## 固定中位数高低组", "",
    "高组严格>中位数；相同值归低组。截点在终点与协变量剔除前固定。下表P为Cox Wald；KM图另标未校正log-rank P。高低组是描述性展示，不是最佳截点搜索。", "",
    table_lines(main[main$model == "median", ]), "", "## 校正、敏感性与反证", "",
    "另运行年龄+性别，以及年龄+性别+预设分期(I/II vs III/IV)校正。严格导管型子集另做全部模型；评分在该子集重新标准化，反卷积值使用同一批标准估计。", "",
    "不以单个未校正P<0.05宣布结论。BH主FDR覆盖每队列/终点/模型/校正方式的两种测量×五细胞10检验；辅助结果另组。方法间一致性是同一RNA的描述性比较，不是独立重复。", "",
    paste0("PH预测变量诊断P<0.05的模型共", nrow(ph), "个，源表见PH_predictor_flags.tsv；这些模型的固定HR需谨慎。"), "",
    "## 不能据此证明", "",
    "不能证明细胞功能或因果关系；分期晚期样本少、治疗信息和肿瘤纯度未控制。共享标记、组成依赖、低比例Treg/CD4区分及模型参考适配性会影响估计。", "",
    "低比例/零值导致分组不平衡不是错误，不随机拆分零值、不删去非显著细胞。DFI事件不足、PFS未测等以not_estimable记录，继续其他模块。", "",
    "## 文件", "", paste0("- `", file.path(out, files), "`"), "",
    "图位于figures：6张连续OS森林图及五细胞×两方法的10张KM图（不可估计时保留图状态，不伪造曲线）；PDF/SVG、300dpi PNG和600dpi TIFF。每张KM另附风险人数TSV。", "",
    "方法来源：[quanTIseq原始研究](https://doi.org/10.1186/s13073-019-0638-6)、[quantiseqr官方教程](https://github.com/federicomarini/quantiseqr/blob/devel/vignettes/using_quantiseqr.Rmd)、[Xena尺度元数据](https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.json)。")
  writeLines(enc2utf8(lines), file.path(out, "REPORT_zh.md"), useBytes = TRUE)
  ledger("report", "completed")
}
ledger(stage, "finished")
cat("Completed stage ", stage, ": ", out, "\n", sep = "")
