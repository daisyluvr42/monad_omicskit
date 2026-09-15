# MonadOmics WorkBuddy kit

2026-09-14：用户已确认采用一个 CLI 连接器 + 一个主 Skill，保留现有生信分析能力。

- 包名与命令：`monadomics`；显示名：MonadOmics 生信工具箱；版本：0.2.0。
- 把现有 Python 调用层和 R 脚本放进标准 `src/monadomics` 包。CLI 接替旧 MCP 接入，不另建任务队列、服务或数据库。
- CLI 提供 doctor、capabilities、schema、setup-r，以及 deg/enrich/plot/survival；分析参数从 JSON 文件读取，结果以 JSON 和本地文件返回。
- WorkBuddy 包只声明 CLI，使用托管 Python；R 本体和 R 包单独准备，连接器 init 不编译 Bioconductor。
- 主 Skill 负责输入判断、调用和结果解释；分析参数和依赖说明按需从 references 读取。
- 验收包括现有 R 回归、CLI 边界、wheel 独立安装、真实分析产物和 ZIP 结构。市场安装需公开可获取的固定版本，不能用本地构建冒充发布成功。

依据：[连接器规范](https://open.workbuddy.cn/en/docs/connector)、[Skill 规范](https://open.workbuddy.cn/en/docs/skill)。本轮现场核对日期为 2026-09-14。

2026-09-15：用户确认增加 Skill 前置的数据提取、规范表保存与原材料校验；随后决定先走 GitHub 安装，要求补齐安装/更新/卸载并推送。采用根目录 `install.py`，把标准库 Python CLI 与 R 文件部署到本机 Skill 的 scripts 目录，固定解释器入口；无需 pip。旧版迁移备份 Skill 和配置，只停用本 kit 的旧 MCP，保留其他连接器与分析数据。市场包仍独立保留。
