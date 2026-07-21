# Omics Kit — 给 AI Agent 的生信分析工具箱

Omics Kit 是一套面向 AI agent 的组学分析 MCP 工具 + Skill。它让 agent 在生信流程中调用 R 完成差异表达、功能富集、出版级作图和预后建模，**所有数值和图都来自确定性计算，模型不许口算，也不许凭记忆写基因名**。

适配 WorkBuddy（腾讯 AI 工作台），也可作为独立 MCP server 接入任何支持 MCP 的 agent harness。

## 为什么单独做一个工具箱

生信分析的依赖（整个 Bioconductor 栈）、运行画像和适用人群都和通用科研工具差别很大。把它塞进通用工具箱会让后者变重。Omics Kit 独立安装、独立更新，和 [monad_scholarkit](https://github.com/daisyluvr42/monad_scholarkit) 互补：**Scholar 管文献、统计、写作、投稿；Omics 管组学分析。** 两者都装时会互相路由，只装一个也能独立工作。

## 设计原则

模型最容易在生信里犯的错不是算错，是**编**——编一个看起来合理的基因名，编一个像模像样的通路 ID。所以这里的分工是：

- **agent 负责**：判断数据类型、选方法、解读结果、把结论写清楚
- **R 工具负责**：所有数值、基因符号映射、富集条目、模型系数和图

工具还会主动拦截生信里最致命的输入错误：把已经 log2 转换或 TPM 标准化的数据喂给 DESeq2。这类错误不报错也能跑出结果，但整条结果链是错的。

## 能力概览

6 个 MCP 工具，21 项能力。

| 模块 | 覆盖 |
|---|---|
| 数据准备 | 表达矩阵核对、数据类型判断、样本 PCA 与批次结构检查 |
| 差异表达 | DESeq2、edgeR、limma/limma-voom、协变量校正、火山图、热图、Venn |
| 功能富集 | GO(BP/CC/MF)、KEGG、Reactome、GSEA、GSVA/ssGSEA |
| 预后模型 | LASSO-Cox 风险评分、高低危 KM、时间依赖 ROC、列线图、校准曲线、DCA |

## MCP 工具

| 工具 | 输入 | 输出 |
|---|---|---|
| `omics_env` | 无 | R 版本、各包状态、缺失包的安装命令 |
| `omics_deg` | 表达矩阵 + 样本表 + 比较方向 | 完整结果表、显著子集、上下调计数 |
| `omics_enrich` | 基因列表或排序列表 + 物种 | 富集表 + dotplot，含未映射基因清单 |
| `omics_plot` | DEG 表或表达矩阵 | 火山图/热图/Venn/PCA，PNG + SVG |
| `omics_survival` | 生存数据 + 候选变量 | 模型系数、风险评分、验证图，含 EPV 警告 |
| `omics_feature_menu` | 无 | 能力菜单与路由提示 |

所有生成文件默认写入 `~/.workbuddy/workspace/omics/`，按 `deg/`、`enrich/`、`figures/`、`survival/` 分类。

## 输入与输出示例

**差异表达**

> 用户：这是我从 GEO 下的 count 矩阵和样本表，帮我比较疾病组和对照组。
>
> 需准备：`counts.csv`（第一列基因符号，其余列为样本）、`samples.csv`（第一列样本名，含分组列）
>
> 得到：完整结果表 + 显著基因表（CSV），上调/下调基因数，前 25 个基因预览，实际使用的设计公式和校正方法

**功能富集**

> 用户：把上调的这批基因做 KEGG 富集。
>
> 需准备：基因符号列表、物种（human/mouse）
>
> 得到：富集通路表（含 P、校正 P、基因数、通路内基因）+ Top10 dotplot，以及未能映射的符号清单

**预后模型**

> 用户：用这 20 个基因建个预后模型。
>
> 需准备：每行一个样本，包含生存时间、结局（0/1）、各基因表达量
>
> 得到：LASSO 筛选后保留的变量及系数、风险评分表、交叉验证图、高低危生存曲线、C-index，以及 EPV 是否过低的警告

## 安装

### 前置条件

- Python 3.11+
- R 4.2+
- WorkBuddy 或任何支持 MCP 的 agent harness

### 步骤

```bash
git clone https://github.com/daisyluvr42/monad_omicskit.git ~/monad_omicskit
cd ~/monad_omicskit
python3 mcp/omics.py install workbuddy
```

安装器只写 MCP 配置和 Skill，**不会**自动装 R 包——首次 Bioconductor 编译可能超过 20 分钟，不应该让它有机会拖垮安装流程。安装完成后会报告 R 状态。

### 装 R 和依赖

```bash
brew install r                    # 没装 R 的话
cd ~/monad_omicskit
Rscript r/bootstrap.R             # 全部依赖，首次约 20 分钟
Rscript r/bootstrap.R --check     # 只看状态，不安装
Rscript r/bootstrap.R enrich      # 只装某一组
```

分组：`core` `deg` `enrich` `plot` `survival`。

R 装在非常规位置时，设 `OMICS_RSCRIPT` 指向 Rscript 路径。

## 更新与卸载

```bash
cd ~/monad_omicskit && python3 mcp/omics.py update workbuddy
```

`update` 自动执行 `git pull --ff-only`、刷新 MCP 配置和 Skill。也可以直接对 WorkBuddy 说"更新 Omics"。

```bash
python3 mcp/omics.py uninstall workbuddy
```

卸载只移除 WorkBuddy 配置和 Skill，不删项目代码、已装的 R 包或分析输出。

## 边界（Agent 必须遵守）

- **基因符号、通路 ID、模型系数只能来自工具输出**，不得由模型生成
- **物种必须显式确认**，人鼠基因符号不通用
- **数据类型判断错误会使整条结果链失效**：count 走 DESeq2/edgeR，芯片和已标准化数据走 limma
- **富集是关联不是机制**，不能写成"证实该通路驱动了表型"
- **训练集的 C-index 和 AUC 是乐观估计**，必须说明并提示外部验证
- **单细胞全流程不在工具层运行**：生成脚本交用户执行，再把下游结果拿回来分析
- 工具不可用或数据不足时明确说明，不口算、不编结果

## 与 Scholar 的分工

| 需求 | 去处 |
|---|---|
| 差异表达、富集、组学作图、预后建模 | Omics Kit |
| 文献检索、精读、综述 | Scholar `sci-literature` / `sci-reading` |
| 普通临床统计、基础 KM 与单变量 Cox | Scholar `sci-stats` |
| Meta 分析 | Scholar `sci-meta` |
| 论文写作、润色、投稿返修 | Scholar `sci-writing` / `sci-polish` / `sci-submit` |

## 测试

```bash
tests/smoke_mcp.sh                    # 不需要 R
python3 -m unittest discover -s tests # 需要 R 和依赖
```
