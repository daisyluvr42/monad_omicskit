---
name: omics-analysis
description: 当用户要做生信/组学分析——差异表达（DESeq2/edgeR/limma）、GO/KEGG/Reactome 富集、GSEA、GSVA/ssGSEA 通路打分、火山图、表达热图、基因集 Venn、样本 PCA、LASSO-Cox 预后模型、时间依赖 ROC、列线图、校准曲线、决策曲线，或 GEO/TCGA 表达矩阵挖掘时使用。Use for R-based omics analysis on expression matrices; not for ordinary clinical-table statistics, meta-analysis, or single-cell pipelines that must run outside the tool layer.
---

# 组学与生信分析

模型负责选方法、判断数据类型和解读结果；所有数值、富集条目、模型系数和图必须来自 Omics 的 R 工具。

## 模型偏好

优先 `deepseek-v4-pro`；其他模型也直接继续。

## 开始前必须确认

1. **物种**：human 还是 mouse。人鼠基因符号不通用，工具会强制要求这个参数，不要替用户假定。
2. **数据类型**：原始 count、TPM/FPKM、芯片信号值，还是已经 log2 转换过。**这一条决定方法选择，判断错了整条结果链都是错的。**
3. **矩阵方向**：基因在行、样本在列，第一列是基因 ID。
4. **分组信息**：样本表的哪一列是分组，比较方向是谁 vs 谁。
5. **基因 ID 类型**：SYMBOL、ENSEMBL 还是 ENTREZ。

会话中第一次分析前先调用 `omics_env`。缺包时把返回的 `install_command` 原样交给用户，不要自己猜安装命令。R 没装时如实说明，不要假装能算。

## 用户需要提供什么

| 分析 | 必须提供 | 可选 |
|---|---|---|
| 差异表达 | 表达矩阵（基因×样本）、样本分组表、比较方向 | 协变量、阈值 |
| 富集分析 | 基因列表（或 DEG 表）、物种 | 背景基因集、本体类别 |
| GSEA | **完整**排序基因列表、物种 | 排序指标列名 |
| 火山图 | DEG 表（含 log2FoldChange 和 padj） | 要标注的基因 |
| 热图 | 表达矩阵、基因列表 | 样本注释列 |
| 预后模型 | 每行一个样本：生存时间、结局(0/1)、候选变量 | 时间点、阈值范围 |

用户材料不全时，明确说缺什么，不要用默认值凑。

## 方法选择

### 差异表达（`omics_deg`）

- **原始 count**：`deseq2`（默认，样本量小时更稳）或 `edger`
- **芯片、TPM/FPKM、已标准化或已 log2**：`limma`
- **count 且样本量大**：`limma` 配 `voom: true`

工具会检查数值特征，发现把已转换数据喂给 DESeq2/edgeR 时直接报错。**报错时不要绕过，回去确认数据类型**——这个检查拦的是生信里最常见的致命错误。

每组至少 2 个重复。批次、性别、年龄等混杂通过 `covariates` 进模型，不要先分层再手工比较。

### 功能富集（`omics_enrich`）

- **有明确 DEG 列表**：`go` / `kegg` / `reactome`（过表达分析 ORA）
- **有完整排序基因列表**：`gsea`。GSEA 必须用全部基因的排序结果，**不能只喂筛选后的 DEG**，工具会拒绝过短的列表
- **要每个样本的通路活性**：`gsva` 或 `ssgsea`，输出后再做组间比较

ORA 默认背景是全基因组。如果实验只检测了部分基因（如靶向 panel），用 `universe` 传入实际检测到的基因，否则富集结果偏乐观。

工具会返回未能映射的基因符号。**如实报告这个数字**，不要悄悄丢掉再报总数。

### 图（`omics_plot`）

- `volcano`：DEG 表 → 火山图
- `heatmap`：表达矩阵 + 基因列表 → 热图，默认行 z-score
- `venn`：2–4 组基因列表 → Venn 图**加每个区域的基因归属表**
- `pca`：表达矩阵 → 样本 PCA

**出图前先看 PCA。** 如果样本按批次而不是按分组分开，先处理批次效应再解读差异。

### 预后模型（`omics_survival`）

- `lasso_cox`：候选变量筛选 + 风险评分 + 高低危 KM
- `timeroc`：时间依赖 ROC
- `nomogram`：列线图
- `calibration`：Bootstrap 校准曲线
- `dca`：决策曲线

事件数（event=1）必须足够。工具会计算 EPV（每变量事件数），低于 10 时返回警告，**必须把警告转达用户**，不能只报 C-index。

工具返回的训练集 C-index 和 AUC 带有乐观偏倚，输出里已附说明，**照实转达**。

## 单细胞与重流程

Seurat 全流程、CellChat、Monocle、SCENIC 运行时间小时级、内存 10GB 以上，**不在工具层运行**。正确做法：

1. 问清参数（组织、分组、marker 来源、QC 阈值、分辨率）
2. 生成参数化的、可直接 `Rscript` 执行的脚本交给用户
3. 用户在自己的环境跑完，把输出（metadata、DEG 表、聚类结果）拿回来
4. 用本 Skill 的工具做下游分析和解读

Pseudo-bulk 之后的差异表达可以用 `omics_deg` 真跑。**不要假装在工具里跑完了单细胞流程。**

## 解读顺序

1. 数据类型与预处理是否恰当
2. 样本量、分组和混杂控制
3. 差异基因数量是否合理（几万个"显著"通常意味着方法或阈值有问题）
4. 效应量与显著性分开讲
5. 富集结果的生物学合理性
6. 局限与需要的验证

## 学术诚信

- **基因符号、通路 ID 和条目名只能来自工具输出。** 绝不凭模型记忆写基因名——这是生信里最高频的事故。
- ID 转换走工具，不猜 ENSEMBL 与 SYMBOL 的对应。
- **富集是关联，不是机制。** 只能写"差异基因显著富集于该通路"，不能写"证实该通路驱动了表型"。
- 多重检验校正方法和阈值必须写进结果描述。
- 预后模型不能用同一批数据既建模又宣称已验证；必须提示需要外部验证。
- 相关不等于因果；共表达不等于调控。
- 公共数据（GEO/TCGA）分析必须写明数据集编号和纳排标准，不编样本量。
- 不伪造基因、通路、富集 P 值、模型系数或图。

## 与其他工具的分工

Omics Kit 只做组学分析。以下如果本机装了 Scholar（`monad_scholarkit`），转过去；没装就如实说明：

- 文献检索、精读、综述：Scholar 的 `sci-literature` / `sci-reading`
- 普通临床表格统计、基础 KM 与单变量 Cox：Scholar 的 `sci-stats`
- 非组学的通用数据图：Scholar 的 `sci-figure`
- Meta 分析：Scholar 的 `sci-meta`
- 把结果写成论文：Scholar 的 `sci-writing`
