# MonadOmics 生信工具箱

MonadOmics 0.2.0 把原 Omics Skill + MCP 整理为 **一个本地 CLI + 一个主 Skill**，另附 WorkBuddy 市场连接器包。本机 Python CLI 调用 R，保留差异表达、功能富集、组学作图和预后建模等 21 项能力；数值和图来自实际计算。

**目前通过 GitHub 安装和更新。** 安装器直接部署仓库中的 CLI 运行文件和主 Skill，不需要 pip、PyPI 上架或腾讯市场审核。市场连接器包保留为后续可选发布方式。

## 能做什么

| 工作 | 实现 |
|---|---|
| 矩阵与样本检查 | 数据类型、样本注释、PCA |
| 差异表达 | DESeq2、edgeR、limma、limma-voom、协变量、多重校正 |
| 功能富集 | GO、KEGG、Reactome、GSEA、GSVA/ssGSEA |
| 组学图 | PCA、火山图、热图、2–4 组 Venn 与区域成员表 |
| 预后模型 | LASSO-Cox、风险评分、KM、时间依赖 ROC、列线图、校准、DCA |

用户可以提供检测公司的结果文件夹、压缩包或完整数据表，附已有报告和样本分组说明。主 Skill 先指导模型使用宿主工具提取规范表格，保存后与原表及报告校验，通过后调用分析。这里只统一文件和字段格式，不改变表达量单位或做统计归一化；详见 [输入整理规范](workbuddy-connector/skills/monadomics-analysis/references/data-preparation.md)。CLI 本身仍接收整理后的数据，不新增厂商导入命令。

这不是完整的 GEO/TCGA 下载器或 Seurat 单细胞流水线。单细胞 pseudobulk count 可以进入差异分析。

## GitHub 安装和第一次运行

需要 Python 3.11+；实际分析还需要 R 4.2+ 和对应 R 包。使用 GitHub CLI 克隆仓库并安装：

```bash
gh repo clone https://github.com/daisyluvr42/monad_omicskit.git
cd monad_omicskit
python3 install.py install
python3 install.py run --version
python3 install.py run doctor --group deg
```

没有 `gh` 时，第一步使用 `git clone https://github.com/daisyluvr42/monad_omicskit.git`。Windows 将 `python3` 换为 `py -3.12`。安装器只使用 Python 标准库，不下载 Python 依赖，也不自动安装 R/Bioconductor。

安装位置为 `~/.workbuddy/skills/monadomics-analysis`，包含主 Skill、输入校验规范和 Python/R 运行文件。安装后的 Skill 会注明本机解释器与 CLI 的完整路径，WorkBuddy 无需依赖终端 PATH 或激活虚拟环境。安装时使用的 Python 解释器需继续保留。安装后刷新技能或重启 WorkBuddy，再用自然语言提交分析任务。

### 从旧 MCP 版升级

在原仓库目录执行下面的命令，不再使用旧的 `mcp/omics.py update workbuddy`：

```bash
git pull --ff-only
python3 install.py install
```

安装器会备份原 `omics-analysis` Skill，并停用配置中指向旧 `mcp/omics_mcp.py` 的 `omics` MCP；其他连接器保持不变。旧配置和 Skill 保存在 `~/.workbuddy/monadomics-backups/`，分析数据与 R 包保留。重启 WorkBuddy 使旧 MCP 停用生效。

### 更新与卸载

```bash
python3 install.py update
python3 install.py uninstall
```

`update` 先检查仓库没有未保存的改动，再执行 `git pull --ff-only`，使用拉取后的新安装器刷新 CLI 与 Skill。下载源码压缩包的用户应下载新包后执行 `install`。卸载将本 kit 的 Skill 和运行文件移入备份目录，不删除 R 包、分析输出，也不会自动重新启用旧 MCP。

如需隔离安装位置，使用 `python3 install.py --workbuddy-dir <目录> install`，之后的 `run`/`update`/`uninstall` 使用同一选项；默认位置无需指定。

### 准备分析环境

R 本体从 [CRAN](https://cran.r-project.org/) 或对应系统渠道准备。R 不在 PATH 时，将 `OMICS_RSCRIPT` 设为实际 Rscript/Rscript.exe 的绝对路径。macOS 会额外识别官方 R 和常见 Homebrew 安装位置。

R 包与连接器初始化分开：

```bash
python3 install.py run setup-r deg
python3 install.py run doctor --group deg
```

模块：`core`、`deg`、`enrich`、`plot`、`survival`、`all`。setup-r 把缺失依赖安装到 R 的个人库 `R_LIBS_USER`，安装日志在 stderr；doctor 不安装任何东西。首次 Bioconductor 安装可能较久，不放进 WorkBuddy init。系统编译库仍由操作系统准备；未安装完不能进行相应分析。

## 调用分析

将参数保存为 `deg.json`，例如：

```json
{
  "method": "deseq2",
  "matrix_type": "counts",
  "matrix_path": "counts.csv",
  "coldata_path": "samples.csv",
  "group_column": "group",
  "treat": "disease",
  "control": "control",
  "output_name": "disease_vs_control"
}
```

```bash
python3 install.py run deg --params deg.json --output-dir analysis-results
python3 install.py run schema enrich
python3 install.py run capabilities
```

- `deg`、`enrich`、`plot`、`survival` 均读取 `--params` 指定的 JSON 对象；`--params -` 读取标准输入。
- 相对输入路径以命令工作目录为准；矩阵为基因×样本，首列基因 ID；样本表首列样本 ID。支持 CSV/TSV 或 schema 说明的内联记录，不直接读取 XLSX。
- 分析成功：退出 0，stdout 为 `ok: true` 的 JSON，含结果与文件路径。分析失败：退出 1，JSON 含 `ok: false` 与错误信息。命令用法错误：退出 2。
- `--output-dir` 控制产物目录；默认 `~/.workbuddy/workspace/omics`，也支持 `OMICS_OUTPUT_DIR`。
- 分析默认超时 900 秒，作图 600 秒；可按需要传 `--timeout`。长任务使用宿主命令进程的等待能力，不反复重跑同一分析。

参数和场景例子见 [主 Skill](workbuddy-connector/skills/monadomics-analysis/SKILL.md) 及其 references。

## 可选的 Python 包与 WorkBuddy 市场包

开发者也可在虚拟环境中执行 `python -m pip install .` 或安装本地 wheel，获得 PATH 中的 `monadomics` 命令；其参数与 `python3 install.py run` 一致。GitHub 安装不需要这一步。

```text
workbuddy-connector/
├── connector-meta.json
├── cli.json
├── icon.svg
└── skills/monadomics-analysis/
    ├── SKILL.md
    └── references/
```

`cli.json` 声明 WorkBuddy 托管 Python 3.12，最低 WorkBuddy 5.0.0；macOS/Linux/Windows 都使用固定版本的 pip 安装命令。无登录和 API Key，不配置虚构的 auth/status/unAuth。连接器只声明 CLI，不混入 MCP。

```bash
python -m pip install build
python -m build
python scripts/build_workbuddy_connector.py
```

生成 Python wheel、sdist，以及 `dist/monadomics-workbuddy-0.2.0.zip`。ZIP 根目录直接包含 connector-meta.json 等文件，符合 [WorkBuddy 连接器规范](https://open.workbuddy.cn/en/docs/connector) 和 [Skill 规范](https://open.workbuddy.cn/en/docs/skill)。

**仅市场连接器路线需要先发布 `monadomics==0.2.0` 到 PyPI，并通过 WorkBuddy 审核。** 本地 wheel 安装和 ZIP 校验不代表已经上架；当前验证范围见 [RELEASE.md](docs/RELEASE.md)。提交 WorkBuddy 的是连接器 ZIP，Python 包由 init 从 PyPI 安装。GitHub 路线直接部署源码，不依赖该市场初始化命令。

旧 0.1 MCP 代码在 Git 历史中保留。0.2 使用新 Skill `monadomics-analysis`，旧版迁移由上述 GitHub 安装器完成。

## 科学与数据边界

- 原始 count 与标准化/log2 表达必须区分；`normalized` 不代表已采用适当 log2 尺度。DESeq2/edgeR 拒绝非原始 count。
- 富集需确认物种、ID 类型和背景。GSEA 需要完整排序列表；GSVA/ssGSEA 需要提供来源明确的 `gene_sets`，不会自动下载默认基因集。
- 基因、通路、P 值、系数不能由模型编造。报告未映射 ID、校正方法和工具警告；富集是关联证据。
- 生存时间单位、事件编码、候选变量来源须明确。训练集 C-index/AUC 不代表外部验证或临床有效性。
- 本地计算不上传表达/生存表至另一分析服务器；R 安装、KEGG 等网络功能会访问公共资源。宿主读取的对话和工具输出按 WorkBuddy 自身规则处理。

## 验证

```bash
python -m pip install -e .
OMICS_OUTPUT_DIR="$PWD/test-output" python -m unittest discover -s tests -v
bash tests/smoke_cli.sh
python tests/run_acceptance.py --output-dir test-output/acceptance
```

测试数据均为固定种子的合成数据，不是公开患者数据或生物学发现。完整测试需要对应 R 包；缺包时后端测试会明确 skip，不能把 skip 当作分析通过。独立 wheel 验收应使用新环境安装 wheel，从源码目录之外运行 `tests/run_acceptance.py`，并检查生成的表和图。
