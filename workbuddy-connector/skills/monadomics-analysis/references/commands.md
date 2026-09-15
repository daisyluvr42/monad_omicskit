# 命令与参数

使用 `monadomics <command> --params params.json`；参数必须为 JSON 对象。分析成功退出 0，stdout 为 `ok: true` 的 JSON；分析/参数错误退出 1，包含 `ok: false`、`error`、`message`。命令行用法错误退出 2。

公共选项：`--output-dir` 指定本地产物目录；`--timeout` 指定 R 运行秒数，默认分析 900、作图 600。文件路径字段支持本地 CSV/TSV，推荐绝对路径。未知参数字段会被拒绝；参数文件不包含公共命令行选项。

下列为 0.2.0 schema；运行 `monadomics schema <command>` 可核对当前安装版。

## deg

Differential expression with DESeq2, edgeR, or limma. DESeq2/edgeR require raw integer counts; use limma for microarray or already-normalised data. Returns the full result table plus the significant subset.

| 参数 | 类型 | 必填 | 默认值/允许值 |
|---|---|---|---|
| `method` | string | 按分析需要 | 默认 deseq2；deseq2, edger, limma |
| `matrix_type` | string | 是 | counts, normalized |
| `matrix_path` | string | 按分析需要 | — |
| `matrix` | array | 按分析需要 | — |
| `coldata_path` | string | 按分析需要 | — |
| `coldata` | array | 按分析需要 | — |
| `group_column` | string | 是 | — |
| `treat` | string | 是 | — |
| `control` | string | 是 | — |
| `covariates` | array | 按分析需要 | — |
| `log2fc` | number | 按分析需要 | 默认 1 |
| `padj` | number | 按分析需要 | 默认 0.05 |
| `padj_method` | string | 按分析需要 | 默认 bh |
| `voom` | boolean | 按分析需要 | 默认 false |
| `output_name` | string | 按分析需要 | — |

## enrich

Functional enrichment: GO/KEGG/Reactome over-representation, GSEA on a ranked list, or GSVA/ssGSEA per-sample pathway scores. Species is required because gene identifiers differ between organisms.

| 参数 | 类型 | 必填 | 默认值/允许值 |
|---|---|---|---|
| `method` | string | 按分析需要 | 默认 go；go, kegg, reactome, gsea, gsva, ssgsea |
| `species` | string | 是 | human, mouse |
| `id_type` | string | 是 | SYMBOL, ENSEMBL, ENTREZID |
| `genes` | array | 按分析需要 | — |
| `universe` | array | 按分析需要 | — |
| `ontology` | string | 按分析需要 | 默认 bp；BP, CC, MF, ALL |
| `ranked_path` | string | 按分析需要 | — |
| `ranked` | array | 按分析需要 | — |
| `gene_column` | string | 按分析需要 | — |
| `metric_column` | string | 按分析需要 | 默认 log2foldchange |
| `matrix_path` | string | 按分析需要 | — |
| `matrix` | array | 按分析需要 | — |
| `matrix_type` | string | 按分析需要 | counts, normalized |
| `gene_sets` | object | 按分析需要 | — |
| `pvalue` | number | 按分析需要 | 默认 0.05 |
| `qvalue` | number | 按分析需要 | 默认 0.2 |
| `top_n` | integer | 按分析需要 | 默认 10 |
| `output_name` | string | 按分析需要 | — |

## plot

Publication figures for expression analyses: volcano, heatmap, Venn, or PCA. Produces 300 dpi PNG plus editable SVG, and for Venn the per-region gene membership table.

| 参数 | 类型 | 必填 | 默认值/允许值 |
|---|---|---|---|
| `type` | string | 是 | volcano, heatmap, venn, pca |
| `deg_path` | string | 按分析需要 | — |
| `deg` | array | 按分析需要 | — |
| `matrix_path` | string | 按分析需要 | — |
| `matrix` | array | 按分析需要 | — |
| `matrix_type` | string | 按分析需要 | counts, normalized |
| `coldata_path` | string | 按分析需要 | — |
| `coldata` | array | 按分析需要 | — |
| `sets` | object | 按分析需要 | — |
| `genes` | array | 按分析需要 | — |
| `gene_column` | string | 按分析需要 | 默认 gene |
| `annotation_columns` | array | 按分析需要 | — |
| `colour_column` | string | 按分析需要 | — |
| `log2fc` | number | 按分析需要 | 默认 1 |
| `padj` | number | 按分析需要 | 默认 0.05 |
| `label_top` | integer | 按分析需要 | 默认 10 |
| `label_genes` | array | 按分析需要 | — |
| `label_samples` | boolean | 按分析需要 | 默认 false |
| `scale_rows` | boolean | 按分析需要 | 默认 true |
| `show_rownames` | boolean | 按分析需要 | — |
| `cluster_columns` | boolean | 按分析需要 | 默认 true |
| `top_variable` | integer | 按分析需要 | 默认 2000 |
| `title` | string | 按分析需要 | — |
| `width` | number | 按分析需要 | — |
| `height` | number | 按分析需要 | — |
| `output_name` | string | 按分析需要 | — |

## survival

Prognostic modelling: LASSO-Cox variable selection with risk score, time-dependent ROC, nomogram, bootstrap calibration, and decision curve analysis. Reports events-per-variable and flags optimistic training-set performance.

| 参数 | 类型 | 必填 | 默认值/允许值 |
|---|---|---|---|
| `method` | string | 是 | 默认 lasso_cox；lasso_cox, timeroc, nomogram, calibration, dca |
| `data_path` | string | 按分析需要 | — |
| `data` | array | 按分析需要 | — |
| `time` | string | 按分析需要 | 默认 time |
| `event` | string | 按分析需要 | 默认 event |
| `predictors` | array | 按分析需要 | — |
| `id_column` | string | 按分析需要 | — |
| `risk_column` | string | 按分析需要 | 默认 risk_score |
| `alpha` | number | 按分析需要 | 默认 1 |
| `nfolds` | integer | 按分析需要 | 默认 10 |
| `lambda` | string | 按分析需要 | 默认 1se；1se, min |
| `seed` | integer | 按分析需要 | 默认 42 |
| `times` | array | 按分析需要 | — |
| `thresholds` | array | 按分析需要 | — |
| `bootstrap` | integer | 按分析需要 | 默认 200 |
| `groups` | integer | 按分析需要 | — |
| `width` | number | 按分析需要 | — |
| `height` | number | 按分析需要 | — |
| `output_name` | string | 按分析需要 | — |

## 输入来源和返回值

- `matrix_path` / `matrix`、`coldata_path` / `coldata` 等分别选择文件或内联记录数组。通常优先文件，避免把大矩阵复制进对话。缺少具体方法需要的数据时 R 会返回可读错误。
- `deg`：`result_table`、`significant_table`、`significant`、`top_genes` 以及实际设计/比较信息。
- `enrich`：按方法返回富集表或通路分数表、图和映射信息；无显著条目不等于命令失败。
- `plot`：`figure.png`、`figure.svg`；PCA 还返回方差解释等，Venn 返回区域成员表。
- `survival`：按方法返回系数、风险分数、图及模型指标；检查所有 `warnings`。

依赖准备命令、输出目录和网络范围见 @references/installation.md。
