# 细胞相关基因的表达、差异富集和多终点生存分析

适用版本：0.4.0。可以从 GitHub 当前源码或本地源码包安装；历史 `v0.3.0` 不含本章新增接口。

## 1. 先确认你实际比较什么

本包现在保留三种不同设计，而不是把它们混称为“细胞内表达”。

| 设计 | 分组依据 | 能回答的问题 | 不能直接证明 |
| --- | --- | --- | --- |
| 旧方法兼容联合评分 | `z(ACLY) + z(Treg_score)` | ACLY–Treg 相关联合高低与临床/转录程序的关联 | 纯 Treg 内 ACLY 高低、两个成分谁驱动关联 |
| 两层背景分析 | 先筛 Treg_score >= 本癌种中位数，再按 bulk ACLY 中位数分组 | Treg 相关背景较高患者中，bulk ACLY 高低的关联 | 测得 Treg 内 ACLY、真实 Treg 细胞数量 |
| 细胞特异数据 | sorted Treg 或按患者汇总的 Treg pseudobulk ACLY | 患者 Treg 表达高低的比较 | 转录关联等同于代谢功能或因果 |

`default_treg_signature()` 使用 FOXP3、IL2RA、CTLA4、CCR8、TIGIT、TNFRSF18、TNFRSF4、IKZF2、LAYN、ENTPD1、LRRC32、BATF、IL2RB。

Treg_score 是本癌种中逐标记基因标准化后的均值。联合模式再对这个均值做一次 z 标准化，然后与 ACLY z 相加，与旧包公式一致。它不是 GSVA/ssGSEA，也不是 CIBERSORT 等真正的细胞反卷积。除非确实使用了 GSVA，图例不标成 Higher/Lower GSVA。

所有新模块输出 `measurement` 和 `score_mode`。bulk 输入必须明确设置 `allow_proxy = TRUE`，防止源码用户误把代理分析解释为纯细胞表达。

## 2. 使用本地 TCGA 做完整示例

安装 0.4.0 包后：

```r
library(TCGASigSurvival)

project <- "/path/to/tcga_project"
out <- file.path(project, "results", "ACLY_Treg_figure_features")
cache <- file.path(project, "results", "panTCGA_ACLY_Treg_prepared_data.rds")
clinical <- file.path(project, "data", "xena",
  "Survival_SupplementalTable_S1_20171025_xena_sp.tsv")

# 如果还没有目标/标记缓存，使用已有接口准备；不会自动下载全部数据。
prepared <- prepare_signature_data(project_dir = project,
  target_gene = "ACLY", signature_genes = default_treg_signature(),
  signature_name = "Treg", out_dir = out, output_prefix = "new_ACLY_Treg")
# 或：prepared <- readRDS(cache)

data <- prepare_tcga_cell_context(prepared, clinical, cell_type = "Treg")

# 与旧包设计相同：联合评分高/低，不再混称为纯 Treg 内表达。
joint <- run_cell_gene_analysis(data,
  target_gene = "ACLY", cell_type = "Treg", score_mode = "joint",
  reference = "adjacent_normal", paired = TRUE,
  endpoints = c("OS", "PFI", "DSS", "DFI", "PFS"),
  min_patients = 80, min_events = 20, allow_proxy = TRUE,
  out_dir = out, prefix = "ACLY_Treg_joint", make_plots = FALSE)

joint$survival$results
joint$comparison$tests

# 两层对照：Treg 背景高，再比较单独 ACLY；不把细胞评分与 ACLY 相加。
two <- run_cell_gene_survival(data,
  target_gene = "ACLY", cell_type = "Treg", score_mode = "gene",
  context_subset = "high", covariates = c("age", "sex", "cell_score"),
  endpoints = c("OS", "PFI", "DSS", "DFI"),
  min_patients = 40, min_events = 10, allow_proxy = TRUE,
  out_dir = out, prefix = "Treg_high_ACLY")

plot_cell_gene_survival_grid(joint$survival)
plot_cell_gene_km(joint$survival, cohort = "LIHC", endpoint = "OS")
```

旧代码的 target_gene 改成 `"ACLY"` 后仍可直接调用，不影响已有 OS 接口。本章多终点分析不通过 `prepare_expression_survival_data()` 预先筛 OS，避免患者因 OS 缺失而失去 PFI/DSS/DFI 资格。

## 3. 肿瘤、邻癌正常、健康人必须分开

TCGA barcode 11 是实体组织正常样本，属于肿瘤患者的邻癌/癌旁正常组织；不是独立健康供者。当前本地 TCGA 文件没有健康供者 Treg 表达数据。

```r
paired <- compare_cell_gene_expression(data, "ACLY", "Treg",
  reference = "adjacent_normal", paired = TRUE, score_mode = "joint",
  allow_proxy = TRUE, out_dir = out)

# 没有健康供者时记录 reference_not_available，不伪造健康对照。
healthy <- compare_cell_gene_expression(data, "ACLY", "Treg",
  reference = "healthy", score_mode = "joint", allow_proxy = TRUE,
  out_dir = out, prefix = "healthy_availability")
```

配对分析只使用同一癌种同一患者同时存在的肿瘤和参考样本。非配对分析中若发现重叠供者，接口阻止将其当作独立样本。每个点是一名供者的一个组织状态；统计表记录配对数、效应、Wilcoxon P 和 BH FDR。

联合评分的肿瘤/正常比较均使用肿瘤样本的参考均值/SD，不能分别对两组做 z 后比较均值。评分是无量纲代理，不是 TPM。

## 4. 差异表达、GSEA 和“冲击图”

参考图 B 的形态是山脊图（ridgeline），不是 Sankey/冲积图。本包输出真正 GSEA 结果以及 leading-edge 基因 moderated t 的山脊分布，横轴明确标注统计量，不把基因密度伪称为 enrichment-score 分布。

需要全基因表达矩阵，不能只用 14 个 ACLY/Treg 标记基因进行全转录组差异和富集。建议先指定一个癌种示例，例如患者数最多的 BRCA，避免混合癌种的组织差异驱动结果。

```r
meta <- data$metadata[data$metadata$cancer_type == "BRCA" &
  data$metadata$tissue_status == "tumor", ]
full <- read_tcga_expression_subset(
  file.path(project, "data", "xena", "tcga_RSEM_gene_tpm.gz"),
  file.path(project, "data", "xena", "gencode.v23.annotation.gtf.gz"),
  sample_ids = meta$profile_id)

brca <- prepare_cell_expression_data(full$expression, meta,
  measurement = "bulk_proxy", expression_scale = "log",
  provenance = data$provenance)
brca$cell_genes <- default_treg_signature()

# 用户自行提供已合法获得的 GMT；不会在加载包时下载通路库。
pathways <- read_gmt_pathways("/path/to/hallmark.symbols.gmt")
enrich <- run_cell_gene_enrichment(brca, pathways = pathways,
  target_gene = "ACLY", cell_type = "Treg", score_mode = "joint",
  allow_proxy = TRUE, out_dir = out, prefix = "BRCA_joint")
plot_cell_gene_enrichment(enrich, cohort = "BRCA", top_n = 20)
```

- 高减低是统一方向，正 logFC / NES 偏向高评分组。
- 已归一化 log RNA 用 limma-trend；原始患者 pseudobulk counts 用 edgeR filterByExpr、TMM、limma-voom。
- GSEA 排名使用所有通过检测条件的基因的 moderated t，不只使用显著差异基因。
- 分组基因 ACLY 默认排除；联合模式还排除 `data$cell_genes`，避免把定义分组的基因再用作独立富集证据。
- 输出全部 DE 和通路表，保留未校正 P、BH FDR、NES、leading edge、基因宇宙和 QC；无隐含 fold-change 阈值。
- GSEA 固定随机种子、串行计算，不随机打散并列排名；记录并列数量。
- 山脊图展示前 20 个可估计通路，依据 FDR 排序，不等于这 20 个都显著。可设置 `significant_only = TRUE` 仅显示 raw P<0.05。
- TCGA bulk 分组后的 DE 是整块组织的程序差异，可能由细胞组成、肿瘤纯度和细胞内状态共同驱动，不是 Treg 内在 DE。

## 5. OS、PFI/PFS、DSS、DFI

TCGA-CDR 当前表含 OS、PFI、DSS、DFI。PFI 是 progression-free interval，不能直接改名为 PFS。接口支持用户队列中真实存在的 `PFS.time/PFS.status`；未测量时返回 `endpoint_not_measured`，不靠重命名制造一个终点。

默认中位数分组在所有符合组织/评分条件的肿瘤患者中一次性确定，随后各终点独立筛时间、事件、协变量。高组 `score >= median`；大量并列时不随机拆分，组过小则不可估计。该规则与旧 OS 实现在先筛 OS 后分组的样本集合上可能不同，属于多终点必要的显式差别，而不是修改评分公式。

每个模型记录输入人数、背景子集人数、纳入人数、事件数、高低组人数、截点、实际协变量、常量协变量、无法估计原因和逐患者排除记录。默认不搜索最显著截点。

Cox 输出 HR、95% CI、Wald P、PH GLOBAL P；KM/log-rank 是未调整描述，不能把校正 Cox P 标成 log-rank P。HR>1 表示高评分组风险更高，不等于某时刻生存率比。

`FDR` 主列对每一种模型的全部癌种×终点联合做 BH；`FDR_within_endpoint` 是分终点的辅助校正。不可估计行保留 NA，而不是 P=1 或“无差异”。网格图叉号代表不可估计，灰边圆点代表可估计但 raw P>=0.05。

由于各癌种随访质量与终点适用性不同，本包运行能力不表示所有癌种四终点都具有同等解释强度；请结合 TCGA-CDR 的癌种终点建议和源表事件数解读。

## 6. 接入真正患者级 Treg 数据

```r
# expr：gene × profile 的规范化 log 表达或原始 counts。
# meta：每个供者/组织状态一个 profile_id。
# 必须含 profile_id, patient_id, cancer_type, cell_type, tissue_status。
# tissue_status 只能明确填 tumor / healthy / adjacent_normal。
# 临床列例如 OS.time, OS.status, PFS.time, PFS.status；时间统一为天。
real <- prepare_cell_expression_data(expr, meta,
  measurement = "sorted", expression_scale = "log",
  provenance = "Study accession; sorted Treg gating; normalization; specimen origin")
res <- run_cell_gene_analysis(real, target_gene = "ACLY", cell_type = "Treg",
  pathways = pathways, reference = "healthy", score_mode = "gene",
  endpoints = c("OS", "PFS", "DSS", "DFI"), out_dir = out)
```

单细胞原始 counts 必须先按患者、癌种、组织状态汇总已注释 Treg，不能把每个细胞当患者。默认每位供者至少 20 个 Treg；排除记录在 `cell_qc`。

```r
# cells 只传入稳定的供者/临床列，不传逐细胞 UMAP/细胞周期等变化列。
pb <- aggregate_cell_pseudobulk(counts, cells, cell_type = "Treg", min_cells = 20,
  provenance = "Study accession; annotation/gating; raw count layer; doublet/QC rules")
res <- run_cell_gene_analysis(pb, target_gene = "ACLY", cell_type = "Treg",
  pathways = pathways, score_mode = "gene", out_dir = out)
```

已运行 CIBERSORTx 等方法的患者级细胞特异表达矩阵，可设置 `measurement = "deconvolved"` 并注明算法、参考图谱和高分辨率模式。仅有 Treg abundance 矩阵不等于具有 Treg 内 ACLY 表达矩阵。

## 7. 可复现运行入口与依赖

新增核心接口为 `run_cell_gene_analysis()`；真实本地 TCGA 批处理入口为 `scripts/example_cell_gene_tcga.R`，有 `preflight`、`survival`、`enrichment`、`figures`、`all` 阶段，原始输入只读，新增输出单独目录。

```sh
Rscript scripts/example_cell_gene_tcga.R /path/to/tcga_project /path/to/new_output preflight
Rscript scripts/example_cell_gene_tcga.R /path/to/tcga_project /path/to/new_output all
```

可选依赖：`BiocManager::install(c("limma", "edgeR", "fgsea", "BiocParallel"))`、`install.packages(c("ggridges", "msigdbr", "svglite"))`。示例脚本通过 `msigdbr` 明确获取公开 Hallmark 并记录库版本，通路定义只在本地缓存，不随包重分发。图形阶段导出可编辑 PDF/SVG、600 dpi TIFF 和 300 dpi PNG 预览；尺寸是独立教学图的展示尺寸，并非指定期刊的最终排版尺寸。

方法资料：

- TCGA-CDR：[Liu et al., Cell 2018](https://pubmed.ncbi.nlm.nih.gov/29625055/)
- 差异分析：[limma 官方文档](https://bioconductor.org/packages/release/bioc/html/limma.html)
- 预排序富集：[fgsea 官方文档](https://bioconductor.org/packages/release/bioc/html/fgsea.html)
- Hallmark 获取：[msigdbr 官方文档](https://igordot.github.io/msigdbr/reference/msigdbr.html)
- Hallmark 定义与许可：[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
