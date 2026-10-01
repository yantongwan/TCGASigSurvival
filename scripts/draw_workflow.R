library(grid)

base_dir <- getwd()
out_png <- file.path(base_dir, "results", "TCGASigSurvival_workflow.png")
out_pdf <- file.path(base_dir, "results", "TCGASigSurvival_workflow.pdf")
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

pal <- list(
  bg = "#fbfaf7",
  ink = "#1f2937",
  muted = "#475569",
  arrow = "#334155",
  blue = "#1d4ed8",
  orange = "#b45309",
  green = "#047857",
  rose = "#be123c",
  blue_fill = "#eff6ff",
  blue_border = "#93c5fd",
  orange_fill = "#fffbeb",
  orange_border = "#f59e0b",
  green_fill = "#ecfdf5",
  green_border = "#6ee7b7",
  rose_fill = "#fff1f2",
  rose_border = "#fda4af",
  note_fill = "#f8fafc",
  note_border = "#94a3b8",
  warn_fill = "#fef2f2",
  warn_border = "#f87171"
)

rect_round <- function(x, y, w, h, fill, border, lwd = 1.7, r = 0.008) {
  grid.roundrect(
    x = unit(x + w / 2, "npc"),
    y = unit(y - h / 2, "npc"),
    width = unit(w, "npc"),
    height = unit(h, "npc"),
    r = unit(r, "npc"),
    gp = gpar(fill = fill, col = border, lwd = lwd)
  )
}

header <- function(x, y, w, label, fill) {
  rect_round(x, y, w, 0.048, fill, fill, lwd = 0, r = 0.006)
  grid.text(
    label,
    x = unit(x + w / 2, "npc"),
    y = unit(y - 0.024, "npc"),
    gp = gpar(col = "white", fontsize = 27, fontface = "bold")
  )
}

box <- function(x, y, w, h, title, lines, fill, border,
                title_size = 20, body_size = 15, title_gap = 0.032) {
  rect_round(x, y, w, h, fill, border)
  grid.text(
    title,
    x = unit(x + 0.014, "npc"),
    y = unit(y - 0.021, "npc"),
    just = c("left", "top"),
    gp = gpar(col = "#111827", fontsize = title_size, fontface = "bold")
  )
  grid.text(
    paste(lines, collapse = "\n"),
    x = unit(x + 0.014, "npc"),
    y = unit(y - 0.021 - title_gap, "npc"),
    just = c("left", "top"),
    gp = gpar(col = pal$ink, fontsize = body_size, lineheight = 1.0)
  )
}

arrow_line <- function(x1, y1, x2, y2) {
  grid.lines(
    x = unit(c(x1, x2), "npc"),
    y = unit(c(y1, y2), "npc"),
    gp = gpar(col = pal$arrow, lwd = 2.2),
    arrow = arrow(type = "closed", length = unit(0.13, "inches"))
  )
}

draw_page <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = pal$bg, col = NA))

  grid.text(
    "TCGASigSurvival：数据、设计思路与 signature 推断流程",
    x = unit(0.5, "npc"),
    y = unit(0.955, "npc"),
    gp = gpar(col = "#202124", fontsize = 38, fontface = "bold")
  )
  grid.text(
    "bulk RNA-seq 目标基因表达 + 细胞 marker signature，高低分组后做 pan-TCGA 生存分析",
    x = unit(0.5, "npc"),
    y = unit(0.918, "npc"),
    gp = gpar(col = pal$muted, fontsize = 18)
  )

  x1 <- 0.04
  x2 <- 0.285
  x3 <- 0.53
  x4 <- 0.775
  w <- 0.205
  y_header <- 0.855

  header(x1, y_header, w, "1. 使用的数据", pal$blue)
  header(x2, y_header, w, "2. 预处理", pal$orange)
  header(x3, y_header, w, "3. Signature 推断", pal$green)
  header(x4, y_header, w, "4. 生存分析", pal$rose)

  box(x1, 0.77, w, 0.115, "表达矩阵", c(
    "UCSC Xena / Toil",
    "tcga_RSEM_gene_tpm.gz",
    "TCGA bulk RNA-seq"
  ), pal$blue_fill, pal$blue_border)
  box(x1, 0.625, w, 0.115, "生存信息", c(
    "PanCanAtlas Survival S1",
    "OS 与 OS.time",
    "癌种与患者 ID"
  ), pal$blue_fill, pal$blue_border)
  box(x1, 0.48, w, 0.115, "表型与注释", c(
    "TcgaTargetGTEX phenotype",
    "GENCODE v23 GTF",
    "symbol 转 gene_id"
  ), pal$blue_fill, pal$blue_border)

  box(x2, 0.77, w, 0.115, "选择基因", c(
    "目标基因：STRAP / DYRK1A",
    "signature：Treg markers",
    "只读取所需表达行"
  ), pal$orange_fill, pal$orange_border)
  box(x2, 0.625, w, 0.115, "整理表达", c(
    "判断是否已 log 转换",
    "必要时 log2(TPM + 1)",
    "重复 symbol 保留最高均值"
  ), pal$orange_fill, pal$orange_border)
  box(x2, 0.48, w, 0.13, "建立 cohort", c(
    "保留肿瘤样本 + OS",
    "常规癌种：primary tumor",
    "SKCM：01 + 06",
    "同一患者只保留 1 个样本"
  ), pal$orange_fill, pal$orange_border)

  box(x3, 0.77, w, 0.135, "默认 Treg markers", c(
    "FOXP3, IL2RA, CTLA4",
    "CCR8, TIGIT, IKZF2",
    "LAYN, ENTPD1, LRRC32",
    "BATF, IL2RB 等"
  ), pal$green_fill, pal$green_border)
  box(x3, 0.595, w, 0.14, "计算 Treg 信号", c(
    "每个癌种内标准化",
    "每个 marker 做 z-score",
    "signature_score =",
    "mean(marker z-scores)"
  ), pal$green_fill, pal$green_border)
  box(x3, 0.415, w, 0.145, "目标基因 + signature", c(
    "target_z：目标基因表达",
    "signature_z：Treg 信号",
    "综合分 = target_z",
    "+ signature_z"
  ), pal$green_fill, pal$green_border)
  box(x3, 0.255, w, 0.12, "方法边界", c(
    "这是 marker-based 推断",
    "不是细胞比例反卷积",
    "不是 Treg 内目标基因表达"
  ), pal$warn_fill, pal$warn_border)

  box(x4, 0.77, w, 0.115, "质控阈值", c(
    "每 cohort ≥ 80 例",
    "死亡事件 ≥ 20 个",
    "可用 signature 基因 ≥ 6"
  ), pal$rose_fill, pal$rose_border)
  box(x4, 0.625, w, 0.12, "high / low 分组", c(
    "Optimized cutoff",
    "surv_cutpoint minprop=0.10",
    "同时输出 median cutoff"
  ), pal$rose_fill, pal$rose_border)
  box(x4, 0.475, w, 0.125, "Cox 与 KM", c(
    "Surv(OS.time, OS.status) ~ group",
    "HR high vs low",
    "95% CI、Wald P、BH-FDR"
  ), pal$rose_fill, pal$rose_border)
  box(x4, 0.315, w, 0.125, "主要输出", c(
    "KM 曲线",
    "forest plot 与 heatmap",
    "survival_results.tsv",
    "QC_report.txt"
  ), pal$note_fill, pal$note_border)

  arrow_line(x1 + w + 0.012, 0.58, x2 - 0.012, 0.58)
  arrow_line(x2 + w + 0.012, 0.58, x3 - 0.012, 0.58)
  arrow_line(x3 + w + 0.012, 0.58, x4 - 0.012, 0.58)

  rect_round(0.10, 0.125, 0.80, 0.095, pal$note_fill, pal$note_border, r = 0.006)
  grid.text(
    "一句话解释",
    x = unit(0.5, "npc"),
    y = unit(0.106, "npc"),
    gp = gpar(col = "#111827", fontsize = 21, fontface = "bold")
  )
  grid.text(
    paste(
      "它问的是：bulk TCGA 肿瘤样本中，目标基因表达高 + Treg signature 强的样本，整体生存是否不同？",
      "所以结果是“目标基因与 Treg-rich 微环境共同定义的风险分组”，",
      "不能直接解释为真实 Treg 比例或 Treg 细胞内目标基因表达。",
      sep = "\n"
    ),
    x = unit(0.13, "npc"),
    y = unit(0.084, "npc"),
    just = c("left", "top"),
    gp = gpar(col = pal$ink, fontsize = 13.8, lineheight = 0.95)
  )

  grid.text(
    "TCGASigSurvival. Data: data/xena/*. Output: results/panTCGA_* and results/KM_plots/*.pdf.",
    x = unit(0.04, "npc"),
    y = unit(0.006, "npc"),
    just = c("left", "bottom"),
    gp = gpar(col = pal$muted, fontsize = 10)
  )
}

png(out_png, width = 3600, height = 2200, res = 150, type = "cairo")
draw_page()
dev.off()

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(out_pdf, width = 24, height = 14.67)
draw_page()
dev.off()

message("Wrote: ", out_png)
message("Wrote: ", out_pdf)
