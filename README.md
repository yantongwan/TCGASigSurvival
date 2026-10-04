# TCGASigSurvival

**R package for gene/signature survival analysis, cell-context and state models, and matched RNA/protein validation.**

Version: **0.4.0** | R >= 4.1 | MIT license

新增：患者级细胞表达比较、全基因 DE/GSEA 山脊图和 OS/PFI/PFS/DSS/DFI 接口，保留旧 ACLY–Treg 联合评分并提供两层敏感性分析。完整教学见 [细胞相关基因功能指南](inst/doc/CELL_GENE_zh.md)。从当前 GitHub 源码安装 0.4.0 可使用新功能；历史 `v0.3.0` 不含这些新增接口。

[English usage guide](inst/doc/USAGE_en.md) | [中文使用指南](inst/doc/USAGE_zh.md) | [完整方法与统计定义](inst/doc/METHODS.md) | [本地发布检查](docs/RELEASE_CHECK.md)

这是面向研究者的可复现肿瘤生存关联分析工具。教程先用模拟数据跑通流程，再换成真实 TCGA 或自己的队列。

> **重要边界：**原有 TCGA 模块分析肿瘤整体 bulk RNA，不会先分离 Naive CD4、Treg 等细胞再测量其基因表达。细胞模块是背景代理评分，不是细胞计数或真正的反卷积；trained-immunity-derived score 不等于已经证明功能性训练免疫。0.4.0 新增接口支持用户自行提供的真实患者级细胞数据，并与 bulk 代理明确分开。

## 目录

- [完整模拟数据练习](#完整模拟数据练习)
- [安装与快速测试](#1-安装与快速测试)
- [TCGA 数据准备](#2-准备-tcga-数据)
- [目标基因与细胞标记](#3-目标基因细胞标记模块)
- [状态基因集](#4-仅检验状态或通路基因集)
- [两层分析及外部反卷积](#5-naive-cd4-背景trained-immunity-状态)
- [自定义输入与临床校正](#6-自定义队列与临床协变量)
- [批量筛查](#7-批量基因集与单基因)
- [表达比较](#8-表达比较与重绘)
- [RNA/蛋白匹配](#9-rna蛋白配对验证)
- [统计解释与排错](#10-输出统计解释与排错)
- [评分与模型的区别](#评分与模型的区别)
- [复现与引用](#复现与引用)

## 可以回答哪些问题

| 问题 | 入口 |
| --- | --- |
| 目标基因与细胞标记的联合信号是否与 OS 相关？ | `run_pan_tcga_signature_survival()` |
| 一个状态或通路模块本身是否与 OS 相关？ | `run_pan_tcga_signature_only_survival()` |
| 校正细胞背景后，状态是否仍与 OS 相关？ | `run_pan_tcga_two_layer_survival()` |
| 当前分析每癌种有多少患者和事件？ | `summarize_tcga_cohorts()` |
| 多个基因/模块如何批量分析？ | `run_pan_tcga_signature_batch()` |
| 严格匹配的 RNA 和蛋白是否一致？ | `match_proteogenomic_tables()` / `run_matched_proteogenomic_analysis()` |

分析单位是患者，按癌种独立拟合；OS 指总生存。人数取决于当前 prepared 的样本选择、OS 和模型完整病例条件，不是下载库的固定人数。

不包含：单细胞聚类注释、内置 CIBERSORT/xCell、自动同源基因转换、完整 CPTAC 下载、药物筛选、分子对接或临床预测模型。

## 完整模拟数据练习

先按第 1 节安装包，然后将本节代码块按顺序运行。两个模拟队列各 120 个虚拟患者；DemoState 的六个基因仅用于软件测试，不是推荐的训练免疫 signature。模拟 HR/P 不是生物学发现。

### A. 输入与人数检查

```r
# tutorial-demo
library(TCGASigSurvival)
demo <- tcgasig_demo_data(n = 120, seed = 42)
prepared <- prepare_expression_survival_data(
  expression = demo$expression, clinical = demo$clinical,
  expression_scale = "log"
)
tutorial_out <- file.path(tempdir(), "TCGASigSurvival_tutorial")
dim(demo$expression)
head(demo$clinical)
summarize_tcga_cohorts(prepared)
```

矩阵是基因 x 样本；临床每行对应样本。患者与样本在演示中一一对应，保留年龄 age。prepared 是后续分析共同使用的输入对象。

### B. 只看状态高低组

```r
# tutorial-demo
demo_state <- run_signature_only_survival(
  prepared = prepared,
  signature_genes = demo$state_genes,
  signature_name = "DemoState",
  cutoff_methods = "median",
  out_dir = file.path(tutorial_out, "state_only")
)
demo_state$results
```

状态高低组仅由状态模块决定，没有控制细胞背景。每基因在癌种内 z-score，再取等权均值。

### C. 分开计算细胞背景与状态

```r
# tutorial-demo
demo_two <- run_two_layer_survival(
  prepared = prepared,
  cell_signature_genes = demo$cell_genes,
  state_signature_genes = demo$state_genes,
  cell_name = "NaiveCD4", state_name = "DemoState",
  cutoff_methods = "median",
  out_dir = file.path(tutorial_out, "two_layer")
)
demo_two$continuous
demo_two$adjusted
demo_two$cell_high
demo_two$skipped
```

continuous 是状态每升高 1 SD 的校正效应；adjusted 是全队列状态高/低的校正效应；cell_high 是背景较高肿瘤患者中的状态高/低比较。不是提取出的 Naive CD4 细胞，不应只挑最小 P。

### D. 临床校正与单基因批量练习

```r
# tutorial-demo
demo_adjusted <- run_two_layer_survival(
  prepared = prepared,
  cell_signature_genes = demo$cell_genes,
  state_signature_genes = demo$state_genes,
  state_name = "DemoState", covariates = "age",
  cutoff_methods = "median",
  out_dir = file.path(tutorial_out, "age_adjusted")
)
demo_adjusted$continuous
demo_batch <- run_signature_batch(
  prepared = prepared,
  signature_sets = list(JUN = "JUN", ATF4 = "ATF4"),
  min_signature_genes = 1, cutoff_methods = "median",
  make_km = FALSE,
  out_dir = file.path(tutorial_out, "batch")
)
demo_batch$results
demo_batch$failures
```

年龄调整不会自动解决全部混杂。单基因是长度为 1 的基因集，需要 min_signature_genes = 1。

### E. 严格匹配 RNA 与蛋白

```r
# tutorial-demo
rna_demo <- demo$clinical[, c("cancer_type", "patient_id")]
rna_demo$target_rna <- as.numeric(demo$expression["ZC3H12C", ])
protein_demo <- rna_demo[, c("cancer_type", "patient_id")]
protein_demo$target_protein <- rna_demo$target_rna +
  0.2 * as.numeric(demo$expression["ATF4", ])
clinical_demo <- demo$clinical[, c("cancer_type", "patient_id", "OS.time", "OS.status")]
matched_demo <- match_proteogenomic_tables(rna_demo, protein_demo, clinical_demo)
matched_demo$qc
demo_protein <- run_matched_proteogenomic_analysis(
  data = matched_demo, target_gene = "ZC3H12C",
  out_dir = file.path(tutorial_out, "matched_demo")
)
demo_protein$correlations
demo_protein$survival
demo_protein$skipped
list.files(tutorial_out, recursive = TRUE)
```

这里蛋白值是人工构造的，不是 CPTAC。用途是学习表结构，不是验证任何生物学结论。临时目录适合练习，正式结果要保存到持久化项目目录。

---

本包将肿瘤 bulk RNA 表达、基因集评分和患者总生存信息连接起来，完成按癌种分析、结果汇总和配对 RNA/蛋白验证。路径均可替换为使用者自己的目录。

## 1. 安装与快速测试

```r
install.packages("remotes")
remotes::install_github("yantongwan/TCGASigSurvival", upgrade = "never")
library(TCGASigSurvival)
packageVersion("TCGASigSurvival")
```

本地安装源码包：

```r
install.packages(c("data.table", "ggplot2", "stringr", "survival", "survminer"))
install.packages("TCGASigSurvival_0.4.0.tar.gz", repos = NULL, type = "source")
```

要求 R ≥ 4.1。原生 R 分析不需要 Python。TCGA 数据在 macOS/Linux 上优先通过 gzip/awk 筛选，Windows 或缺少这些命令时使用 R gzip 流式读取。GitHub 安装不表示已经 CRAN 收录。

以下演示无需网络，也不包含真实患者或未公开 trained-immunity 基因集。演示状态基因仅供测试：

```r
demo <- tcgasig_demo_data(n = 120)
prepared <- prepare_expression_survival_data(demo$expression, demo$clinical)
summarize_tcga_cohorts(prepared)
demo_res <- run_two_layer_survival(
  prepared = prepared,
  cell_signature_genes = demo$cell_genes,
  state_signature_genes = demo$state_genes,
  cell_name = "NaiveCD4", state_name = "DemoState",
  cutoff_methods = "median",
  out_dir = file.path(tempdir(), "TCGASigSurvival_demo")
)
demo_res$continuous
demo_res$adjusted
demo_res$cell_high
demo_res$skipped
```

## 2. 准备 TCGA 数据

```r
project_dir <- "my_tcga_project"
tcga_data_sources()
download_tcga_data(project_dir = project_dir)
```

安装和载入包不会下载数据。下载函数默认复用已有文件、保留失败的 `.part`，写入本地 MD5、字节数和来源。MD5 是文件指纹，并非供应方校验和认证。函数没有断点续传；重跑时复用成功文件，未完成的大文件会重新下载。

| 输入文件 | 来源与用途 |
| --- | --- |
| `tcga_RSEM_gene_tpm.gz` | UCSC Xena Toil TCGA 表达，保留已对数化值 |
| `Survival_SupplementalTable_S1_20171025_xena_sp.tsv` | Pan-Cancer Atlas 临床和 OS |
| `TcgaTargetGTEX_phenotype.txt.gz` | 样本类型和项目等表型 |
| `gencode.v23.annotation.gtf.gz` | GENCODE v23：Ensembl ID 与人类符号映射 |

默认目录 `my_tcga_project/data/xena/`。原项目压缩表达矩阵约 740 MB；分析只提取请求的基因行。只下载小文件可设置 `roles = c("survival", "phenotype", "annotation")`，完整分析仍需表达矩阵。

默认取原发肿瘤 `01`；SKCM 默认合并 `01 + 06`，另建 `SKCM_primary_only` 敏感性队列。每队列每患者仅保留条形码字母顺序最前的样本。患者数与样本数分别报告，重叠的 SKCM 队列不能相加。

## 3. 目标基因＋细胞标记模块

```r
res <- run_pan_tcga_signature_survival(
  project_dir = project_dir,
  target_gene = "ZC3H12C",
  signature_genes = default_treg_signature(),
  signature_name = "Treg",
  cutoff_methods = c("median", "optimized"),
  min_patients = 80, min_events = 20
)
res$survival$results
summarize_tcga_cohorts(res$prepared)
```

历史接口：每癌种内各标记基因 z-score 的平均值为 signature score；再分别标准化该 score 与目标基因表达并相加。高低组由**联合评分**决定。它检验目标基因与 marker 模块的联合信号，不会从肿瘤中提取 Treg 再测量目标基因。

Treg 标记：`FOXP3, IL2RA, CTLA4, CCR8, TIGIT, TNFRSF18, TNFRSF4, IKZF2, LAYN, ENTPD1, LRRC32, BATF, IL2RB`。CD8 使用 `default_cd8_signature()`；自定义 DC 等基因集传入 `signature_genes`，按基因集大小设置 `min_signature_genes`。

## 4. 仅检验状态或通路基因集

输入 CSV 可使用下列结构；基因名是占位符，必须换成有来源依据的真实人类基因符号：

```csv
Gene
GENE_A
GENE_B
GENE_C
```

```r
genes <- read_signature_genes("my_signature.csv", gene_column = "Gene")
state_res <- run_pan_tcga_signature_only_survival(
  project_dir = project_dir, signature_genes = genes,
  signature_name = "MyState", min_signature_genes = 6,
  cutoff_methods = "median"
)
state_res$survival$results
```

评分是癌种内可用基因 z-score 的平均值，不叠加目标基因。

`read_signature_genes()` 去除空值与重复值，默认保留大小写。`uppercase = TRUE` 只做格式转换，**不能代替小鼠到人的同源基因映射**。例如 `Fasl` 对应 `FASLG`，部分 `Zfp` 对应 `ZNF`，不能只转成大写。

`Gene_1.csv` 未随公共包分发。使用者应提供自己有权使用的、已确认物种与方向的基因集。默认是等权、同向均值；若差异基因有相反方向，应先明确如何构成评分，不能全部视为“训练免疫增强”。

## 5. Naive CD4 背景＋trained immunity 状态

第一层默认使用 T-cell/naive-memory 背景标记：

```r
default_naive_cd4_signature()
# CD3D CD3E TRAC CD2 CD4 IL7R CCR7 TCF7 LEF1 SELL LTB MAL
```

这些基因也表达于其他淋巴细胞及 T 细胞状态，不能唯一界定 Naive CD4，也不是绝对细胞比例。第二层为独立定义的状态基因集。

```r
trained_genes <- read_signature_genes("Gene_1_human.csv", gene_column = "Gene")
two <- run_pan_tcga_two_layer_survival(
  project_dir = project_dir,
  cell_signature_genes = default_naive_cd4_signature(),
  state_signature_genes = trained_genes,
  cell_name = "NaiveCD4", state_name = "TrainedImmunity",
  min_state_genes = 20, min_cell_genes = 6,
  min_patients = 80, min_events = 20,
  min_subset_patients = 40, min_subset_events = 10,
  cutoff_methods = c("median", "optimized")
)
two$continuous
two$adjusted
two$cell_high
two$gene_qc
two$cohort_qc
two$skipped
```

| 表 | 回答的问题／内容 |
| --- | --- |
| `continuous` | 全队列：状态每升高 1 SD 的 HR，校正连续细胞背景分数 |
| `adjusted` | 全队列：状态高 vs 低 HR，校正连续细胞背景分数 |
| `cell_high` | 细胞背景不低于癌种中位数的患者中，状态高 vs 低 |
| `scores` | 患者、两个评分、可用基因数和背景分组 |
| `gene_qc` | 癌种、层和每个基因是否可评分 |
| `cohort_qc` | 排除数、两层相关性与重叠基因 |
| `skipped` | 样本／事件不足、常数、截点失败、不收敛等原因 |

两层独立评分、不相加。`cell_high` 中的状态截点在子集内重新计算。两层分数按癌种标准化，不能用其绝对值直接比较不同癌种。

连续模型建议作为主要证据，中位数分组便于绘图。优化截点在同一生存数据中搜索，Wald P 与 BH FDR 未校正该搜索，应列为探索性结果。默认未加入年龄、分期或肿瘤纯度；可以使用自定义 clinical 和 `covariates`。

建议表述：“Naive-CD4-derived 状态基因集的肿瘤 bulk 转录评分在校正独立细胞背景评分后仍与总生存相关。”细胞内训练免疫状态的结论仍需要单细胞、空间或细胞特异实验。

### 接入外部反卷积分数

本包没有内置 CIBERSORT/xCell。可以传入已计算好的结果：

```r
# prepared 包含所需状态基因和相同样本的生存信息
external <- read.csv("external_cell_scores.csv")
# 列：sample_barcode, cell_score；每个样本唯一
two <- run_two_layer_survival(
  prepared = prepared,
  state_signature_genes = trained_genes,
  cell_scores = external,
  cell_score_method = "my_documented_deconvolution_method",
  cell_name = "NaiveCD4", state_name = "TrainedImmunity",
  min_state_genes = 20, cutoff_methods = "median"
)
```

记录算法、版本、参考矩阵、输入尺度，并确认算法支持所需细胞类型。样本 ID 严格匹配；缺失评分会排除。控制细胞背景后仍不能直接得到该细胞内的状态表达。

## 6. 自定义队列与临床协变量

临床表最小结构如下，其他临床列可追加。事件必须显式编码，不能猜测标签含义或将其他生存终点冒充 OS：

| sample_barcode | patient_id | cancer_type | OS.time | OS.status | age |
| --- | --- | --- | --- | --- | --- |
| sample_01 | patient_01 | MyCancer | 365 | 1 | 60 |
| sample_02 | patient_02 | MyCancer | 730 | 0 | 55 |

```r
# expression：数值矩阵，行名为基因符号，列名为样本 ID
# clinical：sample_barcode patient_id cancer_type OS.time OS.status
# OS.time 为天，OS.status：0=删失、1=死亡
prepared <- prepare_expression_survival_data(
  expression, clinical, expression_scale = "log",
  duplicate_patients = "error"
)
```

`expression_scale = "tpm"` 明确使用 `log2(TPM+1)`。原始 counts 先进行恰当正规化。这个入口不猜测表达尺度。基因与样本名需唯一；多样本患者可显式选择 `duplicate_patients = "first_sample"`，检查 `prepared$input_qc`。

`age`、`stage` 等额外临床列会保留，可在两层模型指定 `covariates = c("age", "stage")`。分类变量正确编码为 factor；有缺失的患者排除并记录。`PH_global_p` 是比例风险诊断，假设不满足时需要进一步处理。

历史接口同样支持此 prepared：

```r
single <- run_signature_only_survival(prepared, signature_genes = genes,
  signature_name = "MyState", cutoff_methods = "median")
joint <- run_signature_survival(prepared, target_gene = "ZC3H12C",
  signature_genes = default_treg_signature(), signature_name = "Treg",
  cutoff_methods = "median")
```

## 7. 批量基因集与单基因

```r
batch <- run_pan_tcga_signature_batch(
  project_dir = project_dir,
  signature_sets = list(Treg = default_treg_signature(), CD8 = default_cd8_signature()),
  min_signature_genes = 4, cutoff_methods = "median"
)
batch$results
batch$failures

gene_batch <- run_pan_tcga_signature_batch(
  project_dir = project_dir,
  signature_sets = list(ZC3H12C = "ZC3H12C", DYRK1A = "DYRK1A", CTNNB1 = "CTNNB1"),
  min_signature_genes = 1, cutoff_methods = "median", make_km = FALSE
)
```

批量入口一次提取所有所需基因，避免重复读取全矩阵。保留历史单基因集 `FDR`，并添加跨基因集／癌种的 `batch_FDR`，按截点方法和 default_cohort 分族。应事先确定检验族，不能按显著性挑选。

## 8. 表达比较与重绘

```r
plot_signature_results(results = res$survival$results,
  target_gene = "ZC3H12C", signature_name = "Treg",
  out_dir = file.path(project_dir, "results"))
plot_target_expression_by_signature(
  survival = res$survival, target_gene = "ZC3H12C", signature_name = "Treg",
  cutoff_method = "median", out_dir = file.path(project_dir, "results")
)
multi <- run_multi_signature_target_expression(
  project_dir = project_dir, target_gene = "ZC3H12C",
  signature_sets = list(Treg = default_treg_signature(), CD8 = default_cd8_signature()),
  min_signature_genes = 4, make_km = FALSE
)
```

此处 `res` 应来自第 3 节目标＋signature 分析。历史表达图使用历史联合评分定义肿瘤高低组；signature-high 不意味着已分离该细胞。比较的是肿瘤整体目标表达，不能称为“Treg 内目标表达”。

## 9. RNA／蛋白配对验证

原生 R 接口接受已正规化、同癌种同研究同患者的表。RNA、蛋白都需 `cancer_type, patient_id`，另分别有 `target_rna`、`target_protein`；临床需 `OS.time, OS.status`。RNA 可含 `signature_score`。

```r
matched <- match_proteogenomic_tables(rna_table, protein_table, clinical_table)
matched$qc
cptac <- run_matched_proteogenomic_analysis(
  data = matched, target_gene = "ZC3H12C", signature_name = "Treg",
  out_dir = file.path(project_dir, "results", "cptac"),
  min_patients = 80, min_events = 20
)
cptac$correlations
cptac$survival
cptac$qc
cptac$skipped
```

严格使用癌种＋患者键，不修改 ID、不模糊匹配、拒绝重复键。使用者必须核实原始患者标识和样本对应关系；相同字符串不能证明是同一生物样本。

输出 RNA-蛋白 Spearman 相关、RNA/蛋白各自连续和中位数 Cox，以及可选 marker score 校正。CPTAC 作为独立队列验证，不将无关 TCGA RNA 与 CPTAC 蛋白拼成配对样本。旧 TCGA 联合评分与单独 RNA/蛋白模型效应定义不同，比较时必须核对。

也可读取既有外部 CPTAC 流程结果：

```r
legacy <- run_cptac_proteogenomic_validation(
  project_dir = project_dir, proteogenomics_dir = "my_cptac_pipeline",
  target_gene = "ZC3H12C", signature_name = "Treg", rerun = FALSE
)
legacy$summary$qc_metrics
legacy$summary$rna_protein_concordance
legacy$summary$tcga_cptac_replication
```

`rerun = TRUE` 要求已部署外部 shell/Python 流程、数据和依赖。该外部下载工作流未随 R 包分发，缺失文件不会自动补齐。公共包没有完整 CPTAC 自动抓取、药物筛选或分子对接模块。

## 10. 输出、统计解释与排错

默认输出至 `project_dir/results/`，不同分析名称可以避免覆盖。两层统一前缀 `panTCGA_<cell>_<state>_two_layer_*`：结果表、森林图 PDF/PNG、子集 KM PDF、QC 和 sessionInfo。`make_km = FALSE, make_summary_plots = FALSE` 可关闭图形输出。

HR > 1 表示高组死亡风险较高，HR < 1 表示较低；连续 HR 对应状态每升高 1 SD。CI 为 95% 区间。最低患者数只是筛选门槛，不能替代事件数量、模型稳定性和功效评估。

| 问题 | 处理 |
| --- | --- |
| Missing input file | 检查文件名、data_dir，或显式下载 |
| 基因找不到／覆盖不足 | 查看 gene_qc，核对物种、符号和 GENCODE 版本 |
| 某癌种没有结果 | 查看 skipped 的患者／事件／变异／模型失败原因 |
| 不收敛或无限 HR | 已记入跳过报告；不解释极端估计 |
| 与旧两层脚本略有差异 | 0.3.0 按默认／敏感性队列分族，并排除失败模型 |
| 需要临床调整 | 提供自定义 clinical 与显式 covariates |

这是研究用关联分析。细胞混合、肿瘤纯度、signature 来源、临床混杂、截点选择与比例风险假设都影响解释。高 trained-immunity-derived score 不等于功能性训练免疫已经证实。

帮助：`help("run_two_layer_survival")`。引用：`citation("TCGASigSurvival")`。模型公式和来源见 [METHODS](inst/doc/METHODS.md)。复现时保存版本、输入来源/指纹、基因集、参数、QC 与失败队列。

## 评分与模型的区别

对癌种 c、基因 g、患者 i，`z[g,i,c] = (expression[g,i,c] - mean[g,c]) / sd[g,c]`。使用 R 的样本 SD，忽略缺失、排除常数基因；模块是可用基因 z 的等权均值。新两层模型还检查每个患者的可用基因数。两层先按癌种评分，再在模型完整病例中标准化连续项。

| 分析 | 分组/公式 | 可以解释什么 |
| --- | --- | --- |
| 历史目标＋模块 | `J = z(target) + z(module)`；`Surv(...) ~ group(J)` | 联合信号的高低组差异 |
| 模块单独 | `Surv(...) ~ group(module)` | 模块整体转录评分的关联 |
| 两层连续 | `Surv(...) ~ state_score_z + cell_score_z + covariates` | 控制背景后的状态每 1 SD 效应 |
| 两层分组 | `Surv(...) ~ state_high_low + cell_score_z + covariates` | 控制背景后的高/低效应 |
| 背景较高子集 | `Surv(...) ~ state_high_low + covariates` | 背景较高肿瘤患者内的高/低差异 |

HR 是瞬时风险比，不是某时间点的存活率比或生存时间比。连续 HR 与高/低 HR 的单位不同。95% CI 横跨 1 表示仍与 HR = 1 相容。KM 是未调整的描述性曲线，即便对应 Cox 校正了协变量；Cox P 不应称为 log-rank P。

历史表达图的联合评分本身包含目标表达，因此高/低组的目标表达比较不独立于分组设计。不能据此称为“细胞内目标表达差异”。

不同癌种的评分绝对值不能直接比较。共线性、细胞组成、肿瘤纯度、年龄、分期、批次、signature 来源和方向都可能影响效应。连续/中位数模型可作为主要视角；优化截点需要独立验证或选择感知设计。研究用关联不等于因果、临床预测或细胞内机制。

## 复现与引用

```r
packageVersion("TCGASigSurvival")
sessionInfo()
citation("TCGASigSurvival")
system.file("doc", "USAGE_zh.md", package = "TCGASigSurvival")
```

从源码 checkout 开发与验证（先安装本包）：

```sh
R CMD build .
R CMD check TCGASigSurvival_0.4.0.tar.gz --no-manual
Rscript scripts/validate_readme.R README.md /path/to/tutorial_output
Rscript scripts/validate_release.R /path/to/validation_output
```

教程测试只提取带 `# tutorial-demo` 标记的 R 代码块并按顺序执行，不下载真实数据。真实 TCGA 流程测试见 `scripts/validate_real_tcga.R`。GitHub Actions 配置了 Linux/macOS/Windows 检查，远端通过情况以实际工作流为准；本地 macOS 已完成构建、安装、自动测试及真实 TCGA 流程验证。

公共仓库不包含患者数据、私人联系信息或未公开基因集。代码 MIT；外部数据及原始方法需要分别引用并遵守许可。

数据与方法来源：[UCSC Xena](https://xena.ucsc.edu/)、[GENCODE v23](https://www.gencodegenes.org/human/release_23.html)、[TCGA clinical resource](https://doi.org/10.1016/j.cell.2018.02.052)、[Toil expression resource](https://doi.org/10.1038/nbt.3772)、[survival](https://cran.r-project.org/package=survival)、[survminer cutoff](https://rpkgs.datanovia.com/survminer/reference/surv_cutpoint.html)。

提交问题：[GitHub Issues](https://github.com/yantongwan/TCGASigSurvival/issues)。请用去标识化的小型复现示例，不在公开 Issues 上传患者数据。
