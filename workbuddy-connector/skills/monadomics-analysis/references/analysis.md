# 方法与边界

## 差异表达

原始整数 count 可使用 DESeq2、edgeR 或 limma-voom；已标准化表达使用 limma。`normalized` 不等于已 log2，也不等于已完成适当预处理，需结合数据来源、分布和处理记录判断。各组需有重复；不能只凭样本数机械选择方法。

`treat` 相对 `control` 决定 log2 fold change 正负。`covariates` 纳入设计公式。输出包括实际设计、样本信息、全部基因表、显著子集和预览。报告阈值、BH 等校正方法及上下调数量。

当前接口先取 treat/control 两组再拟合（`fit_scope: two_group_subset`），并非多组联合模型的 contrast。`sample_selection` 列出输入、纳入、其他组排除及仅在元数据中的样本；`group_levels` 按 control/treat 排列。缺失样本注释、组标签或协变量、混杂设计和无剩余自由度会报错；不要改标签绕过。矩阵空/重复基因 ID、重复样本列和非有限数值也会被拒绝，聚合必须在上游按有依据的规则完成并记录。

需要 gene symbol 时同时提供 `species` 和 `id_type`。返回 `annotation.mapping_table` 保存 INPUT/ENTREZID/SYMBOL 及缺失映射，记录数据库包版本；结果表的 `symbol` 只填唯一映射。未映射或一对多并不使统计值无效，保留原 ID，必要时查阅完整映射表。未请求注释时 `annotation` 为 null，不能猜测名称。

DESeq2 的 `results(alpha)` 与 `padj` 参数一致。显著定义为校正 P 严格小于阈值、绝对 log2FC 大于等于阈值且效应不为零。低计数预过滤与模型内独立过滤分开报告，使用返回的 `filtering`。CSV 空白统计值表示缺失，不能改填 1；`direction: unavailable` 不属于显著或可判定的不显著结果。选显著基因只取 `up`/`down`，不能简单筛 `direction != ns`。

## 富集

ORA（GO/KEGG/Reactome）使用候选基因及适当背景。未显式提供 `universe` 时，工具使用注释资源所支持的背景；对于靶向 panel 或过滤明显的实验，应提供实际可入选背景，避免选择偏倚。

工具核对映射后的候选集是否包含背景之外的基因。`annotation` 与 `background_annotation` 区分实际映射和未映射；`effective_counts` 来自检验的 GeneRatio/BgRatio 分母，GO ALL 可按 ontology 分别返回，无可用检验时为 null。分母不等于提交数，也不等于映射数。有效分析的空结果保存带表头的表并返回 0 条显著项，不要求阳性结果。上游零显著时跳过 ORA；全部 ID 无法映射时应修正物种/ID 类型，不能当作“无显著富集”。`directionality: unsigned_gene_set` 表示 ORA 本身不包含效应方向；即使分开分析上调和下调，富集也不能直接称为通路激活或机制证实。

GSEA 需要未经显著性筛选的完整、有方向的排序列表；说明实际使用的排序指标。不能用一个很短的 DEG 清单冒充完整排序。`pvalue` 控制原始 P，`padj` 控制 BH 校正 P，均默认 0.05，使用严格小于；显式传入 `qvalue` 时额外要求该值非缺失且小于阈值，不会忽略此参数。默认 `seed: 42`，返回随机种子与软件版本。相同环境和参数应可复跑；跨数据库/软件版本不能保证数值不变。

GSEA 的 `result_table` 为所有统计有效的检验条目，`significant_table` 为满足阈值的子集；`terms_tested`、`terms_invalid`、`terms_significant` 分别报告有效检验数、无效行数和显著数。点图和 `top_terms` 仅使用显著子集，横轴 NES 表示排序方向，颜色为校正 P。没有显著条目时交付空的显著表并明确说明，不强行生成阳性图。ORA 的 `pvalue` 同时用于原始 P 与 BH 校正 P 筛选，`qvalue` 默认 0.2，照实报告返回阈值。

GSVA/ssGSEA 使用逐样本表达矩阵，输出通路分数。自定义 `gene_sets` 为通路名到基因 ID 数组的对象，应与矩阵 ID 命名系统一致；工具返回 `kcdf` 等实际设置。通路分数并不直接等于通路激活的因果证据。

报告输入、映射与未映射基因数量；未映射 ID 不由模型猜测补齐。人类与小鼠之外的物种不在当前实现范围。

## 图

PCA 显示样本整体结构。热图默认逐行 z-score，用于展示相对表达差异；不等于差异表达检验。count 在 PCA/热图中先对完整矩阵做 edgeR TMM 归一化，再以 `prior.count=2` 计算 logCPM；返回每个样本的文库大小、归一化因子与有效文库大小。热图通过 `genes` 在变换之后选行，不能只给差异基因矩阵来估计归一化。PCA 在变换后选择高变基因，仅中心化，不把各基因缩放到相同方差。这些步骤不等同于批次校正，也不是 VST。

已有来源明确的 VST/rlog 或其他合适的变换矩阵，可用 `matrix_type: normalized` 原样作图，并记录外部预处理。原始 counts 仍用于 DESeq2/edgeR 模型，不用可视化矩阵替代。

Venn 支持 2–4 组，同时返回各区域的基因表。火山图使用 DEG 表的效应量和校正 P 值。图文件为 PNG 和 SVG，交付前实际打开检查。

## 预后

LASSO-Cox 返回所选变量、系数、风险分数和 KM。时间依赖 ROC、列线图、Bootstrap 校准、DCA 使用相应生存数据与参数。时间单位和事件编码由真实资料确定；`times` 必须与生存时间使用相同单位。

候选变量数相对事件数过多时，转述工具 EPV 警告并解释过拟合风险。交叉验证选择惩罚参数不等于独立验证；同一数据集拟合、分组和评估的表现有乐观偏倚。不得依据好看的训练曲线宣称具备临床应用价值。
