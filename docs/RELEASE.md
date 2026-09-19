# 0.2.1 发布说明

## 本次统计修正

- DESeq2 的独立过滤 alpha 与请求的校正 P 阈值一致；保留 NA，区分预过滤、统计值不可用与不显著结果。
- GSEA 全表与显著子集分开，显式执行校正 P 及可选 q-value 筛选；排除无效行，计数/预览/图使用同一子集，横轴显示 NES，记录随机种子。
- count 的 PCA/热图改用完整矩阵 TMM + logCPM；PCA 中心化且不按基因标准化，返回实际归一化信息。
- Skill 增加结果表与报告一致性、来源方法与比较方向检查。没有特定研究的基因名单、预期数量或专用分支。

这些修改会改变此前生成的结果和图，应按明确的参数重新分析，不沿用旧报告中的计数。

## 包含的文件

- `monadomics-0.2.1-py3-none-any.whl`：Python 运行包，包含 R 脚本，可本地安装。
- `monadomics-0.2.1.tar.gz`：Python 源码分发包。
- `monadomics-workbuddy-0.2.1.zip`：WorkBuddy CLI 连接器包，一个主 Skill。

## GitHub 发布与安装（当前路线）

2026-09-15：用户选择先通过 GitHub 分发，保留 CLI + 主 Skill 和数据整理校验步骤。仓库根目录的 `install.py` 提供 install、update、uninstall 和 run；安装只复制随仓库提供的 Python/R 文件，无 pip/PyPI 或腾讯审核依赖。使用方法见 README。

发布前验证：隔离 WorkBuddy 目录中的首次安装、旧 MCP/Skill 迁移、升级后的新安装器执行、卸载、错误返回和原文件保留；在安装后实际运行 R 分析。推送时使用正常 fast-forward，不覆盖远程变更；推送后从远程重新克隆并检查安装路径。

安装后的 Skill 增加本机绝对命令入口；其 references 和所带 Python/R 文件必须与源文件一致。仓库与市场包不得带开发者本机路径。旧 MCP 所在本机工作树不因发布动作而被强制更新，运行 `install.py install` 的用户才会触发本机迁移。

## 可选市场发布顺序

1. 在独立 Python 环境安装 wheel，检查 doctor，并运行 `tests/run_acceptance.py`。检查输出表和图；R/系统依赖未就绪不能宣称分析通过。
2. 发布固定版本 `monadomics==0.2.1` 到 PyPI。
3. 从新的环境按 `cli.json` 中的命令安装公开版本，检查版本和分析结果。
4. 向 WorkBuddy 开放平台提交连接器 ZIP，核对 `source: monadomics` 的唯一性，完成审核。
5. 从 WorkBuddy 市场安装，验证 Skill 触发、附件路径、CLI PATH、生成文件打开，以及更新/卸载行为。

GitHub 发布不执行上述市场步骤，不填写未获平台分配的 ID，也不声称经过市场审核。

## 当前验证

2026-09-19：0.2.1 在本机完成原有 47 项测试及新增统计回归测试。新增测试使用配对合成 counts，与官方 DESeq2 在两个阈值下逐值比较；检查缺失统计值、零效应、TMM/logCPM 与官方函数一致、先归一化后选热图基因、PCA 坐标、GSEA 原始 P/校正 P/q-value 筛选及无结果输出。另以小鼠 GO 注释和随机效应量验证固定种子重复运行一致；不指定任何研究数据集的目标数量或基因集合。

隔离环境安装 wheel 后，实际通过 doctor、PCA、DESeq2、火山图、热图、GO、GSEA 和 LASSO-Cox 八步验收。完整 CLI JSON、参数与产物保存于本地 `test-output/wheel-acceptance/`；统计测试日志位于 `test-output/statistical-final.log`。Python wheel、sdist 与 WorkBuddy ZIP 均构建成功。R 脚本解析和已变更图形的 PNG/SVG 检查通过。WorkBuddy 包按仓库专用校验器验证；Codex Skill 校验器不接受 WorkBuddy 的多语言与版本 frontmatter，不以它作为该格式的验收标准。

上述证据验证后端统计规则和本机安装流程，不等于完整生物学基准评测，也没有验证任意厂商 PDF/Excel 提取或模型生成报告的全部行为。现有分析输出不会随升级自动重算。

2026-09-15 GitHub 安装路线：macOS 上 47 项测试通过、0 skip，其中 6 项覆盖新安装器；在含中文及空格的隔离目录完成安装，通过安装器 `run` 实际执行 doctor、PCA 和 DESeq2。安装测试还检查移走源码后独立运行、最小 PATH、旧配置/输出保留，以及从本地 Git 远程拉取更新并执行新安装器。R 实算使用已有 R 与任务隔离包库，不代表新机器已具备 R 依赖。未通过真实 WorkBuddy 对话测试厂商表格提取，Linux/Windows 尚未实跑。

本轮的具体测试数、独立 wheel 验收和产物检查记录见交付目录中的 `验收说明.md`、`acceptance/acceptance.json` 和 `tests.log`。源码 README 中的复跑命令用于再次验证。

本机以外的 Linux/Windows 已提供相同的 Python init 命令，但尚需对应操作系统实测。R 本体、编译器和 Bioconductor 包不由 WorkBuddy 的 Python runtime 自动提供。

## 依据

2026-09-14 现场核对：[WorkBuddy Connector](https://open.workbuddy.cn/en/docs/connector)、[WorkBuddy Skill](https://open.workbuddy.cn/en/docs/skill)。采用文档已定义的 Python runtime 与 init；没有添加文档未给出完整结构的 versionCheck 配置。

无登录能力，因此不声明 auth/status/unAuth。R 缺依赖由 doctor/setup-r 处理，不伪装为“授权失败”。
