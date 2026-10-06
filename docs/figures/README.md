# README 例图说明

本目录仅提供 PDF 例图，保持原分析文件的完整页面、坐标、统计量与图形样式，不裁剪、不改写结果。GitHub README 使用相对链接打开 PDF 文件页面；不会使用失效的 `![...](...pdf)` 图片语法。

## 例图与来源

| 例图 | 原结果文件 | 源表 |
| --- | --- | --- |
| tcga_acly_tumor_adjacent.pdf | ACLY_Treg_Figure_features_20261003/figures/A_bulk_ACLY_paired_adjacent.pdf | bulk_ACLY_paired_tests.tsv、bulk_ACLY_paired_values.tsv |
| brca_acly_treg_hallmark.pdf | ACLY_Treg_Figure_features_20261003/figures/B_BRCA_joint_Hallmark_ridges.pdf | BRCA_joint_enrichment.tsv、BRCA_joint_leading_edge.tsv |
| pan_tcga_acly_treg_survival.pdf | ACLY_Treg_Figure_features_20261003/figures/C_joint_survival_grid.pdf | ACLY_Treg_joint_survival_results.tsv |
| paad_immune_abundance_forest.pdf | PAAD_Immune_Abundance_20261004/figures/PAAD_OS_continuous_unadjusted.pdf | analysis/survival_results.tsv |
| paad_cd4_marker_score_km.pdf | PAAD_Immune_Abundance_20261004/figures/PAAD_OS_KM_CD4_marker_score.pdf | analysis/survival_inputs.tsv、analysis/survival_results.tsv |
| paad_cd4_estimated_fraction_km.pdf | PAAD_Immune_Abundance_20261004/figures/PAAD_OS_KM_CD4_estimated_fraction.pdf | analysis/survival_inputs.tsv、analysis/survival_results.tsv |
| paad_ductal_immune_abundance_forest.pdf | PAAD_Immune_Abundance_20261004/figures/PAAD_Ductal_OS_continuous_unadjusted.pdf | analysis/survival_results.tsv |

原结果目录均位于运行者项目的 `results/` 下。源表在本地保留审计，不随公共例图目录分发。每张例图和原文件的字节数、MD5 一致，见 [manifest.tsv](manifest.tsv)。

## 方法与图注

**图1：配对表达比较。** 每个癌种至少3对同患者肿瘤/癌旁标本，20个癌种合计676对；正常在左、肿瘤在右。统计为供者级配对 Wilcoxon，原源表保留原始P与BH FDR。纵轴是正规化log表达，不是细胞比例。癌旁正常不是健康志愿者，bulk ACLY 不是 Treg 细胞内 ACLY。

**图2：BRCA 富集。** 1091位肿瘤患者，联合高组546、低组545。分组使用 `z(ACLY)+z(Treg_score)`；limma-trend 检验19609个蛋白编码基因，排除ACLY和13个分组标记。全部检验基因的 moderated t 排名进入 fgsea，不先筛显著差异基因；图中显示20个可估计 Hallmark。山脊是 leading-edge 基因 t 的描述性密度，不是 GSEA ES 分布，也不代表患者级不确定性。参考集合在本地使用，未在此分发。

**图3：多终点生存。** 32个癌种有显式状态，OS/PFI/DSS/DFI各自筛查患者和事件，不能将总数据库人数作为每个模型的n。颜色为log2(HR)，点大小为截顶后的-log10(FDR)，边框区分原始Wald P<0.05与否，叉号是不可估计。联合中位数模型的主FDR范围为全部可估计癌种×终点。PFS未测，不进入该图。

**图4：PAAD双方法。** 主队列178位独立原发肿瘤患者，OS完整177人/93事件；每细胞、每测量分别拟合未校正连续Cox，HR为每1参考队列SD、横线为95%CI。marker模块不是经验证的细胞数量计；quanTIseq使用完整线性TPM、TIL10/lsei、肿瘤设置及mRNA校正。比例分母包括Other，总CD4=非Treg CD4+Treg，巨噬细胞=模型M1+M2。

**图5和图6：固定中位数KM。** 截点在结局剔除前确定；严格高于中位数为高，等值归低。177人/93事件，两图均显示95%生存区间、删失标记和组内人数/事件。评分图log-rank P=0.238；比例图P=0.000707。P是未校正log-rank，不是临床校正后检验。比例图高组89人/35事件、低组88人/58事件；高低组Cox HR=0.491（0.323–0.748），BH FDR=0.00909。主FDR覆盖同队列/终点/模型/校正方式下两测量×五主要细胞。不能把它解释为CD4的因果抗肿瘤作用。

**图7：敏感性。** 仅保留明确导管型标签的147人，OS有效146人/85事件。该子集与主队列重叠；评分在子集重算标准化，比例复用同一标准估计。总CD4连续Cox HR=0.762（0.555–1.046）、P=0.0924、FDR=0.308；预测变量PH诊断P=0.0177，固定HR解释需谨慎。该敏感性没有通过主FDR，不能作为独立复现。

## 同步与发布

```sh
Rscript scripts/prepare_readme_examples.R /path/to/TCGA_project
Rscript scripts/validate_readme_assets.R .
```

同步脚本只复制已完成的PDF，不重跑分析。如果目标例图已存在但MD5不同，停止而不是覆盖。发布GitHub时需要同时提交README与此目录，不能只提交README。PDF是二进制文件，不应当作UTF-8文本上传；`scripts/github_manifest.R --paths` 可列出包括二进制例图的完整发布路径。
