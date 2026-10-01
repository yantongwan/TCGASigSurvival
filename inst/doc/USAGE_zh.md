# TCGASigSurvival 0.3.0 使用指南

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
install.packages("TCGASigSurvival_0.3.0.tar.gz", repos = NULL, type = "source")
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

帮助：`help("run_two_layer_survival")`。引用：`citation("TCGASigSurvival")`。模型公式和来源见同目录 `METHODS.md`。复现时保存版本、输入来源/指纹、基因集、参数、QC 与失败队列。
