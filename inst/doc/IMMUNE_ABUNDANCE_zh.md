# 五类免疫细胞与生存：双方法教学指南

适用于 TCGASigSurvival >= 0.5.0。本模块不需要指定 ACLY 等目标基因，也不先筛选细胞评分较高的患者。

## 1. 问题与边界

研究的是 CD8、Treg、中性粒细胞、巨噬细胞、CD4 的 bulk RNA 特征或估计丰度与患者生存的关系。

两条路线并列，不互相替代：

```text
独立患者的原发肿瘤 bulk RNA + 临床生存
  ├─ log表达 → 每基因z-score → 五类marker均值 → 连续Cox + 固定中位数KM
  └─ 全转录组线性TPM → quanTIseq/TIL10 → 五类估计比例 → 连续Cox + 固定中位数KM
```

marker score 不是百分比；quanTIseq 比例不是流式、病理计数或单细胞实测。高表达状态、细胞数量和组织纯度均可能改变 bulk 信号，关联不能直接证明细胞功能或因果关系。

## 2. 安装

本地源码：

```r
install.packages("/path/to/TCGASigSurvival", repos = NULL, type = "source")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("quantiseqr", ask = FALSE, update = FALSE)
library(TCGASigSurvival)
packageVersion("TCGASigSurvival")
```

quantiseqr 是外部可选依赖。比例分支缺失依赖或输入不合格时记录 blocked，特征评分分支继续。不要将这称为已完成两种方法。GitHub 安装能否取得 0.5.0 取决于实际远端发布，不能由本地版本推断。

## 3. 现有 TCGA 数据的一键运行

要求项目含三个原始文件：

```text
data/xena/tcga_RSEM_gene_tpm.gz
data/xena/gencode.v23.annotation.gtf.gz
data/xena/Survival_SupplementalTable_S1_20171025_xena_sp.tsv
```

入口不重新下载 FASTQ，不重跑上游测序流程，不修改这些输入。缓存按文件 MD5、文件大小、癌种和加载器版本核验，不一致时要求新建输出目录。

```sh
Rscript /path/to/TCGASigSurvival/scripts/analyze_paad_immune_abundance.R \
  --project-dir /path/to/TCGA_project \
  --out-dir /path/to/new_versioned_results \
  --stage all
```

可用阶段：prepare、analyze、figures、report、all。已存在分析结果时 all/analyze 拒绝覆盖；重绘或更新报告用对应阶段。隔离安装可增加 `--library-dir /path/to/R_library`。

该入口运行全部 PAAD 与严格 `Pancreas-Adenocarcinoma Ductal Type` 子集。后者与主队列重叠，不能称为外部验证。每患者选字典序最前的原发肿瘤 RNA 标本，其他标本保留排除原因。样本选择不依赖生存数据完整性。

## 4. 在 R 内调用两条方法

```r
d <- prepare_tcga_abundance_data(
  project_dir = "/path/to/TCGA_project", cancer_type = "PAAD",
  cache_file = "/path/to/new_results/PAAD_input.rds"
)
res <- run_immune_abundance_analysis(
  expression = d$expression, clinical = d$clinical,
  expression_scale = d$expression_scale,
  endpoints = c("OS", "PFI", "DSS", "DFI", "PFS"),
  covariate_sets = list(
    unadjusted = character(),
    age_sex = c("age", "sex"),
    age_sex_stage = c("age", "sex", "stage_group")
  ),
  out_dir = "/path/to/new_results/analysis"
)
res$module_status
res$marker$gene_qc
res$fraction$settings
res$survival$results
```

TCGA-CDR 的 PFI 不改名为 PFS；没有 PFS 列时返回 endpoint_not_measured。OS/DSS/PFI/DFI 可用人数和事件各不相同。`stage_group` 为预先固定的 I/II 对 III/IV，未知值作为缺失；晚期样本少时解释需谨慎。这里没有自动校正肿瘤纯度。

## 5. 评分的基因与计算

```r
default_immune_signatures()
```

Treg 与 CD8 沿用包原先的 `default_treg_signature()` 和 `default_cd8_signature()`。另外三个是可审计的探索性标记模块，而不是已经验证的细胞数量计：

|类别|默认额外模块|
|---|---|
|中性粒细胞|CEACAM8, FCGR3B, CSF3R, CXCR2, FPR1, S100A8, S100A9|
|巨噬细胞|C1QA, C1QB, C1QC, CSF1R, CD68, CD163, MSR1, MRC1|
|CD4|CD4, CD3D, CD3E, TRAC, CD2, IL7R, LTB|

每个基因在分析队列内计算 `(log表达 - 均值) / 标准差`，再对标记取等权均值。不使用目标基因乘积，不只分析 cell-high 子集。恒定基因不能标准化，标记缺失或不合格均记录。每样本至少达到 `ceiling(0.7 × 标记数)` 个可用基因，否则评分 NA，不以零代替缺失。

CD8 的细胞毒性基因可能来自 NK；CD4 含其他 T 细胞共享基因；Treg 活化标记不完全专一；髓系模块受炎症和极化影响。模块高分应解释为“相关 bulk 表达信号更高”。

自定义基因集应先给出独立理由，不能看到 P 值后选基因：

```r
panels <- default_immune_signatures()
panels$CD8 <- c("CD8A", "CD8B") # Example only; not a validated counter
scores <- score_immune_signatures(d$expression, d$clinical,
  signatures = panels, expression_scale = "log", min_coverage = 0.7)
```

对于普通线性 TPM，评分接口用 log2(TPM+1)。对于本地 Toil 文件，统一入口保留原有 log2(TPM+0.001) 作为评分输入，另做明确逆变换供反卷积；不重复 log，不猜测量纲。

## 6. quanTIseq 是如何估计比例的

使用外部 `quantiseqr::run_quantiseq()` 的 TIL10 参考表达矩阵，用约束最小二乘分解 bulk 表达。非负细胞贡献受总体约束；mRNA 含量校正后输出各参考细胞比例和 Other。不是把一个标记的表达直接转换成百分数。

本轮固定参数为 TIL10、lsei、RNA-seq、`is_tumordata=TRUE`、`scale_mRNA=TRUE`、标准 noisy/tumor-aberrant genes 去除。不寻找让生存 P 最小的设置。依赖版本、参考 MD5、有效参考覆盖率和每样本比例和写入结果。

输入要求完整人类基因符号矩阵和线性 TPM，至少 1000 基因，非负有限值；marker-only 矩阵、原始 counts 或 log TPM 不合格。有效参考基因精确符号覆盖至少 80%。外部实现另进行基因重注释及总量归一化。其自动尺度检测有局限，包装器拒绝最大值小于 50 的不安全输入，不静默自动修正。

五个主要输出：

|输出|定义|分母|
|---|---|---|
|CD8|T.cells.CD8|所有模型细胞，包括 Other|
|Treg|Tregs|同上|
|Neutrophils|Neutrophils|同上|
|Macrophages|Macrophages.M1 + Macrophages.M2|同上|
|CD4|T.cells.CD4 + Tregs，即总 CD4|同上|

另外输出非 Treg CD4、模型 M1 和 M2。M1/M2 是参考模型类别，不能作为功能极化的证明。CD4 包含 Treg，因此这五个主要比例不能简单相加。完整十一项模型输出（含 Other）才应合计为 1。

低 Treg/CD4 比例和相似表达背景会降低可分辨性；没有方法保证每个肿瘤的估计都准确。推荐后续用独立单细胞、病理或流式验证，而非据此宣称实测比例。

## 7. 生存模型和统计输出

每一种测量、每一种细胞分别拟合，不把五个相关变量同时塞入一个模型。临床校正集各单独运行。

- 连续 Cox 为主要分析，HR 表示增加一个参考队列 SD 的风险比。
- 比例另给增加 10 个百分点的 HR/95%CI；小比例 SD 很小时该尺度可能超出常见观察区间，不应过度外推。
- 高低组为补充：严格高于队列中位数为高组，等于中位数归低组。不优化截点，不随机拆分零值。
- 中位数与 SD 在终点/协变量排除前固定；终点或调整模型改变后，不重新分组。
- 默认至少 40 患者、10 事件；中位数组每组至少 10 人；每设计参数至少 5 事件是探索性输出门槛，不代表功效充足。
- HR>1 是较高测量与较高事件风险相关；HR<1 是较低风险相关，不是因果效应。
- `wald_p` 为对应 Cox 检验，`logrank_p` 为未校正高低组检验，不等同于调整后检验。
- BH `FDR` 覆盖两方法×五主要细胞；同一队列/终点/模型/调整集为一组。`FDR_within_method` 是该方法的五细胞范围。辅助细胞另组。
- PH_global_p、PH_predictor_p 是比例风险诊断。PH 不满足时固定 HR 的解释需谨慎，非显著诊断也不能证明假设成立。
- 缺失、未测、常量、分组不平衡、样本不足、事件不足、秩亏和拟合失败各保留原因。不能把 NA 当作不显著。

## 8. 图表

```r
p <- plot_cell_abundance_forest(res$survival, cohort = "PAAD",
  endpoint = "OS", adjustment = "unadjusted")
print(p)
km <- plot_cell_abundance_km(res$survival, cell_type = "Treg",
  measurement = "estimated_fraction", cohort = "PAAD", endpoint = "OS")
print(km)
attr(km, "risk_table")
```

KM 含 95% 生存区间、删失标记、组内患者/事件数、中位数截点和 log-rank P。风险人数附在对象属性，CLI 同时导出配对 TSV。CLI 导出 PDF/SVG、300dpi PNG、600dpi TIFF。森林图同时显示两方法连续 HR/95%CI；所有不显著结果保留。

## 9. 自己的队列 / 已有反卷积结果

临床要求 sample_barcode、patient_id、cohort_label；结局列如 OS.time（天）、OS.status（0/1）。矩阵列名与样本必须精确匹配，不能自动丢失不匹配患者。每队列每患者只能一个观测。

已有真实比例可直接整理为长表，避免再次反卷积：

```r
# abundance columns: sample_barcode, patient_id, cohort_label, cell_type,
# value, measurement, method, unit, denominator, role
# measurement="estimated_fraction", unit="fraction", value in [0,1]
# denominator must describe the original denominator; do not mix methods.
external <- run_cell_abundance_survival(abundance, clinical,
  endpoints = "OS", covariates = c("age", "sex"), adjustment = "age_sex")
```

不同方法的分母可能不同，不能混合为同一个比例变量。相同方法/细胞/队列中 unit、denominator、role 必须一致。不同细胞总量约束和共享标记会引起相关，不能将多模型 P 视作独立重复。

## 10. 输出与复现

`analysis/survival_results.tsv` 是全部主/辅助/敏感性/终点/调整模型源表。
`analysis/survival_inputs.tsv` 是逐患者的输入与排除账本。
`analysis/abundance.tsv` 是两种测量的全部数值。
`analysis/full_fractions.tsv` 是含 Other 的完整 quanTIseq 比例。
`marker_definitions.tsv`、marker_gene_qc、marker_sample_qc、marker_overlap 记录基因和覆盖率。
`fraction_gene_qc`、fraction_sample_qc 和 RDS settings 记录反卷积审计。
`input_manifest.tsv`、`sessionInfo.txt`、`run_settings.rds` 和 `run_ledger.tsv` 记录输入、版本及执行。
`method_agreement.tsv` 仅为同源 RNA 方法间描述性一致性；不是外部验证。

不要公开上传真实患者级表或原始表达，除非确认共享合规；公开源代码包不包含本轮患者数据。

方法依据：[quanTIseq 原始研究](https://doi.org/10.1186/s13073-019-0638-6)、[quantiseqr 官方教程](https://github.com/federicomarini/quantiseqr/blob/devel/vignettes/using_quantiseqr.Rmd)、[本地 Toil 表达尺度元数据](https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.json)。
