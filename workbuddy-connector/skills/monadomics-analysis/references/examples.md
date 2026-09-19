# 调用示例

下面是参数写法，输入路径与分组名必须替换为当前用户的真实数据；这里不提供假分析结果。保持参数文件以便复跑。

## 原始 count 的疾病组 vs 对照组

先确认样本表中 `group` 列确实包含 `disease` 和 `control`，矩阵样本与表中样本一致。

```json
{
  "method": "deseq2",
  "matrix_type": "counts",
  "matrix_path": "counts.csv",
  "coldata_path": "samples.csv",
  "group_column": "group",
  "treat": "disease",
  "control": "control",
  "padj": 0.05,
  "log2fc": 1,
  "output_name": "disease_vs_control"
}
```

保存为 `deg.json`，在输入文件所在工作目录执行：

```bash
monadomics doctor --group deg
monadomics deg --params deg.json --output-dir analysis-results
```

## 样本 PCA

```json
{
  "type": "pca",
  "matrix_type": "counts",
  "matrix_path": "counts.csv",
  "coldata_path": "samples.csv",
  "colour_column": "group",
  "output_name": "sample_pca"
}
```

```bash
monadomics plot --params pca.json --output-dir analysis-results
```

火山图参数使用 `type: volcano`、`deg_path` 取前一步真实 `result_table`；热图使用 `type: heatmap`、原矩阵及从结果表选出的 `genes`。

## 富集分析

ORA 参数使用 `method: go`、`species`、`id_type` 和真实筛选出的 `genes` 数组。实际可入选的背景基因写入 `universe`，不要使用模型自行编写的基因清单。

GSEA 使用 `method: gsea`、`ranked_path` 指向完整排序表，`gene_column`/`metric_column` 填真实列名；完整声明 `species` 和 `id_type`。例如 `pvalue: 0.05, padj: 0.05, seed: 42`；除非分析计划需要额外筛 q-value，否则不传 `qvalue`。显著数量读取 `terms_significant` 并与 `significant_table` 核对，不数全表代替。

GSVA/ssGSEA 使用 `matrix_path`、`matrix_type`、物种/ID 类型和必填的 `gene_sets`。`gene_sets` 对象中的每一项为一个基因集名称及其成员数组，ID 应与矩阵一致。

```bash
monadomics enrich --params enrichment.json --output-dir analysis-results
```

## LASSO-Cox

```json
{
  "method": "lasso_cox",
  "data_path": "survival.csv",
  "time": "time",
  "event": "event",
  "predictors": ["gene1", "gene2", "gene3"],
  "seed": 42,
  "output_name": "prognostic_model"
}
```

```bash
monadomics survival --params survival.json --output-dir analysis-results
```

例中的候选变量名称仅表示列名写法，必须替换为数据中存在且有依据的候选变量。时间依赖 ROC 可使用返回的风险分数表，`risk_column` 填其中真实列名。
