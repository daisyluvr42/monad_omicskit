---
name: monadomics-analysis
display_name: MonadOmics 生信分析
display_name_en: MonadOmics Bioinformatics
description: Prepare and verify tables from sequencing-company deliveries, then analyze bulk or pseudobulk expression data with local R tools for differential expression, enrichment, figures and prognostic models. Not for ordinary clinical statistics or full single-cell pipelines.
description_zh: 从测序公司交付材料提取规范表格，与原表及报告核验后，完成差异表达、功能富集、组学作图和预后建模。
description_en: Prepare and verify vendor data tables, then analyze expression, run enrichment, create figures and fit prognostic models locally.
version: 0.2.3
author: MonadOmics
---

# 生信分析

使用 `monadomics` CLI 调用 R。依据数据选择方法、解读结果；统计值、基因映射、通路条目、模型系数和图形必须来自实际计算。

## 分析前：输入含义与比较目标是否明确

1. 接收公司交付的文件夹、压缩包、完整表格及原报告，先扫描文件、说明和已有对话，读取 @references/data-preparation.md。识别研究目标及本次任务需要的物种、数据类型、样本分组、比较方向和配对/批次等背景；已有明确依据的信息直接使用，不要求用户重新制作已有表格或重复回答。
2. **清洗和分析前，集中询问会影响本次方法或结论的缺失、冲突或仅能推测的信息。** 带上已识别的内容，用实验事实提问，不要求用户填写工具参数；研究目标不能由数据代替用户决定。必要信息随任务变化，不使用固定问卷。将来源和用户补充简记在 `input-check.md`，后续复用；用户也不清楚时保留未知，只暂停依赖它的步骤，继续文件盘点等不受影响的工作。
3. 使用宿主文件工具或本地脚本提取数据，按规范保存本次分析需要的表格，保留原材料。此处的“归一化”仅指格式与字段标准化，不改变 counts、TPM、FPKM、对数尺度、缺失值或显著性数值。
4. **保存后重新读取规范表及原材料，核对矩阵与样本表的字段含义、ID、分组和数值。** 原始表的行名、列名与说明共同决定字段含义；不能把“写出值与同一个内存对象一致”当成解析正确，也不能只检查计数而遗漏元数据。将来源、列映射及实际结果简记在 `input-check.md`。
5. 提取脚本修改后，在独立目录复跑并核对最终矩阵和样本表，再使用产物；临时命令修复也要回写脚本。已核验且未改动的材料直接复用，检查仅重做受改动影响的部分。后续发现新矛盾或信息缺口时，回到上述规则处理受影响的步骤；无需逐步审批。

## 调用分析

1. 会话中首次分析先运行 `monadomics doctor --group deg`（富集用 `enrich`，单独作图用 `plot`，预后用 `survival`）。检查 JSON 的 `ok` 和 `missing`。缺依赖时读取 @references/installation.md，按用户授权范围准备环境；未就绪时不生成分析结果。
2. 读取 @references/commands.md 中对应命令的参数；需要具体调用写法时读取 @references/examples.md。也可运行 `monadomics schema deg` 等查看当前安装版本的完整 schema。
3. 把参数写入 UTF-8 JSON 文件，再执行 `monadomics deg --params analysis.json --output-dir analysis-results`。其他分析使用 `enrich`、`plot`、`survival`。相对数据路径以命令的工作目录为准；实际任务优先写附件的绝对路径。
4. 等待命令完成，检查退出码、JSON 的 `ok`、产物是否存在。R 默认限时为分析 900 秒、作图 600 秒；按输入规模需要可用 `--timeout` 调整。宿主支持后台命令时用其等待机制继续跟踪同一个进程，避免因界面暂未返回重复启动分析。

## 工具返回后：实际执行与结果是否有效

工具自动检查结构和方法前提。模型对照任务核对返回的实际样本、设计、方向与预处理；参数可接受、命令成功不等于研究问题已回答。DEG 的 `fit_scope: two_group_subset` 表示排除其他组后拟合，`sample_selection` 给出实际样本；不能把它描述为所有组联合模型。需要联合模型时明确当前接口边界，不用多次两组调用冒充。

读取落盘结果按实际阈值核对显著子集与报告数量；零显著、PCA 分离不明显、部分未映射通常是可交付结果，不按预期基因、条目数或任意映射率门槛判失败。GSEA 使用 `significant_table` 和 `terms_significant`；ORA 的 `genes_mapped` 是已验证的 ID 映射数，实际注释分母见 `effective_counts`，不可混用。做对应分析前读取 @references/analysis.md 的相关小节。

用户询问能做什么时运行 `monadomics capabilities`，按需求介绍相关能力。

## 分析选择与关键约束

- **数据类型**：`matrix_type: counts` 仅用于原始非负整数 count；`normalized` 用于非 count 表达量。先确认是否已 log2，不能仅因数据是 TPM 就假定尺度已适合线性模型。DESeq2/edgeR 不能接收 TPM、FPKM 或 log2 数据。count 用 limma 时显式传 `voom: true`。
- **比较方向**：确认 `group_column`、`treat`、`control`。有批次等混杂时用 `covariates` 纳入设计；先检查 PCA 和样本注释，区分分组与批次结构。
- **基因名称**：需要使用 symbol 时，DEG 同时提供有依据的 `species`/`id_type`，读取返回的 `annotation` 与 `symbol`；只有唯一映射才填 symbol，歧义和未映射保留原 ID。也可使用已有且来源明确的注释表，不能凭印象补全名称。可选注释依赖用 `doctor --group enrich` 检查。
- **富集**：物种仅支持 `human`/`mouse`，明确 `id_type: SYMBOL/ENSEMBL/ENTREZID`。ORA 的 `universe` 应反映实际进入检测/筛选范围的背景基因；GSEA 使用完整排序结果，不能只用显著 DEG。GSVA/ssGSEA 必须提供表达矩阵和来源明确的 `gene_sets`，当前工具不会自动下载默认基因集。
- **作图**：火山图使用真实 DEG 表，缺失校正 P 值的基因不参与绘图并单独计数。热图和 PCA 必须声明 `matrix_type`；count 使用完整矩阵做 TMM 归一化和 logCPM 变换。热图用 `genes` 选行，不先截取少量基因再估计归一化因子。`normalized` 不会再次转换，须确认已经适合可视化。Venn 同时交付区域基因归属表。
- **预后**：生存时间单位明确，结局必须为 0/1，时间点与生存时间同单位。照实报告事件数、EPV 警告和训练集性能的乐观偏倚；训练集 C-index/AUC 不能表述为外部验证。

工具拒绝数据类型、分组或 ID 时，回到输入检查，不改标签绕过验证。具体方法判断和边界见 @references/analysis.md。

## 交付前：报告中的关键事实是否有依据

保存报告后读回，重点核对会改变结论的基因名称、比较方向、方法、数字和引用，分别对应注释、实际 CLI 输出、结果表或已读取的来源。缺少依据就保留原 ID、缩小表述或说明尚未验证；不把假设写成结果，不把未核验当成通过。发现矛盾时回到对应步骤修正，无需另建检查台账或机械重跑整套分析。

打开并检查生成的图，核对比较方向、坐标、标签、样本注释和是否有遮挡。向用户提供规范输入表、`input-check.md`、结果表和图文件的可点击路径，并用实际输出说明方法、阈值、多重校正、主要结果及限制；保留提取脚本（如使用）、分析参数和完整 CLI 返回 JSON 以便复跑。说明哪些结果沿用公司原分析、哪些由本次重新计算。

与原报告或论文对比时，先从对应实验的方法和表格确认单位、对数底数及比较方向，不由模型印象补全。统一方向后再计算一致性指标；回归必须使用图中实际的 x/y 方向，表、图、JSON 和正文由同一份核对结果生成。同源数据复算不能称作独立外部验证。

未映射 ID、低事件数和其他工具警告必须说明。混合上下调 ORA 只支持无方向的功能关联；方向解释需对应有方向的分析，上下调数量本身不能证明增殖或功能激活，富集也不证明机制。只有通过来源核对的事实可进入任务总结；写入记忆时仍保留适用范围和不确定性，不把临时排错办法当成通用规范。

Seurat、CellChat、Monocle、SCENIC 完整流程不由本 kit 执行；如用户需要，可另行生成有参数说明的脚本，标明尚未运行。用户提供 pseudobulk count 后可使用本 kit 的差异分析。普通临床统计、文献检索和论文写作由其他适用工具处理，本 kit 不依赖其他连接器。
