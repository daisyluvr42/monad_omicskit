---
name: monadomics-analysis
display_name: MonadOmics 生信分析
display_name_en: MonadOmics Bioinformatics
description: Prepare and verify tables from sequencing-company deliveries, then analyze bulk or pseudobulk expression data with local R tools for differential expression, enrichment, figures and prognostic models. Not for ordinary clinical statistics or full single-cell pipelines.
description_zh: 从测序公司交付材料提取规范表格，与原表及报告核验后，完成差异表达、功能富集、组学作图和预后建模。
description_en: Prepare and verify vendor data tables, then analyze expression, run enrichment, create figures and fit prognostic models locally.
version: 0.2.0
author: MonadOmics
---

# 生信分析

使用 `monadomics` CLI 调用 R。依据数据选择方法、解读结果；统计值、基因映射、通路条目、模型系数和图形必须来自实际计算。

## 第一步：数据整理与原材料校验

1. 接收公司交付的文件夹、压缩包、完整表格及原报告，结合用户目标读取 @references/data-preparation.md。先从已有材料确定物种、数据类型、样本与分组、比较方向；不要求用户重新手工制作材料中已有的表格。
2. 使用宿主文件工具或本地脚本提取数据，按规范保存本次分析需要的表格，保留原材料。此处的“归一化”仅指格式与字段标准化，不改变 counts、TPM、FPKM、对数尺度、缺失值或显著性数值。
3. **保存后重新读取规范表，再读取原始表格逐项核对 ID、样本与数值，并与原报告核对单位、分组、方法和完整性。** 不能只检查模型刚写出的预览，也不能用报告中的节选冒充全量数据。将来源、列映射及实际校验结果简记在 `input-check.md`。
4. 校验通过后继续分析，参数只引用已核验的表格或从这些表格读取的数据。发现不一致先修正并重验；仅在缺失信息或冲突无法由材料解决时询问用户，暂停受影响的分析。无需每次另请用户批准。用户已提供规范表时仍须核验，可直接复用而不重复转换。

## 第二步：调用分析

1. 会话中首次分析先运行 `monadomics doctor --group deg`（富集用 `enrich`，单独作图用 `plot`，预后用 `survival`）。检查 JSON 的 `ok` 和 `missing`。缺依赖时读取 @references/installation.md，按用户授权范围准备环境；未就绪时不生成分析结果。
2. 读取 @references/commands.md 中对应命令的参数；需要具体调用写法时读取 @references/examples.md。也可运行 `monadomics schema deg` 等查看当前安装版本的完整 schema。
3. 把参数写入 UTF-8 JSON 文件，再执行 `monadomics deg --params analysis.json --output-dir analysis-results`。其他分析使用 `enrich`、`plot`、`survival`。相对数据路径以命令的工作目录为准；实际任务优先写附件的绝对路径。
4. 等待命令完成，检查退出码、JSON 的 `ok`、产物是否存在。R 默认限时为分析 900 秒、作图 600 秒；按输入规模需要可用 `--timeout` 调整。宿主支持后台命令时用其等待机制继续跟踪同一个进程，避免因界面暂未返回重复启动分析。

用户询问能做什么时运行 `monadomics capabilities`，按需求介绍相关能力。

## 分析选择与关键约束

- **数据类型**：`matrix_type: counts` 仅用于原始非负整数 count；`normalized` 用于非 count 表达量。先确认是否已 log2，不能仅因数据是 TPM 就假定尺度已适合线性模型。DESeq2/edgeR 不能接收 TPM、FPKM 或 log2 数据。count 用 limma 时显式传 `voom: true`。
- **比较方向**：确认 `group_column`、`treat`、`control`。有批次等混杂时用 `covariates` 纳入设计；先检查 PCA 和样本注释，区分分组与批次结构。
- **富集**：物种仅支持 `human`/`mouse`，明确 `id_type: SYMBOL/ENSEMBL/ENTREZID`。ORA 的 `universe` 应反映实际进入检测/筛选范围的背景基因；GSEA 使用完整排序结果，不能只用显著 DEG。GSVA/ssGSEA 必须提供表达矩阵和来源明确的 `gene_sets`，当前工具不会自动下载默认基因集。
- **作图**：火山图使用真实 DEG 表。热图和 PCA 必须声明 `matrix_type`；count 会做 `log2(count + 1)`，normalized 不会再次转换。Venn 同时交付区域基因归属表。
- **预后**：生存时间单位明确，结局必须为 0/1，时间点与生存时间同单位。照实报告事件数、EPV 警告和训练集性能的乐观偏倚；训练集 C-index/AUC 不能表述为外部验证。

工具拒绝数据类型、分组或 ID 时，回到输入检查，不改标签绕过验证。具体方法判断和边界见 @references/analysis.md。

## 结果交付

打开并检查生成的图，核对比较方向、坐标、标签、样本注释和是否有遮挡。向用户提供规范输入表、`input-check.md`、结果表和图文件的可点击路径，并用实际输出说明方法、阈值、多重校正、主要结果及限制；保留提取脚本（如使用）和分析参数以便复跑。说明哪些结果沿用公司原分析、哪些由本次重新计算。

未映射 ID、低事件数和其他工具警告必须说明。富集支持关联解释，不能据此写成机制已证实；结果中出现的基因、通路和数值均须能追溯到输入或输出。

Seurat、CellChat、Monocle、SCENIC 完整流程不由本 kit 执行；如用户需要，可另行生成有参数说明的脚本，标明尚未运行。用户提供 pseudobulk count 后可使用本 kit 的差异分析。普通临床统计、文献检索和论文写作由其他适用工具处理，本 kit 不依赖其他连接器。
