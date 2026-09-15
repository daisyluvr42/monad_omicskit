# 运行环境

GitHub 安装使用仓库根目录的 `python3 install.py install`，直接部署主 Skill 与 Python/R 运行文件，不使用 pip，也不需要市场审核。安装后的主 Skill 会给出本机命令入口，以下所有 `monadomics` 命令都使用该入口；从仓库终端也可用 `python3 install.py run` 执行相同子命令。Windows 使用 `py -3.12`。

GitHub 更新执行 `python3 install.py update`，卸载执行 `python3 install.py uninstall`。从旧版迁移时先 `git pull --ff-only`，再运行新安装器，随后重启 WorkBuddy；不要继续使用旧 `mcp/omics.py` 的更新命令。

可选的市场连接器声明 Python 3.12，初始化只安装 `monadomics` Python 包。两条路线都包含 R 脚本，不包含 R 可执行文件、R 包或系统编译库。

## 检查

```bash
monadomics --version
monadomics doctor --group deg
```

doctor 只检查，不安装、不写入授权状态；全部所需依赖就绪时退出 0 且 `ok: true`。未就绪时退出非 0，并返回缺失包或 R 路径问题。模块选择：`core`、`deg`、`enrich`、`plot`、`survival`、`all`。

## R 本体

最低 R 4.2，建议使用与当前 Bioconductor 配套的受支持 R 版本。

- macOS：可从 [CRAN](https://cran.r-project.org/bin/macosx/) 安装官方 R，以便使用可用的二进制包。Homebrew R 也能使用，但部分 Bioconductor 依赖需要源码编译。
- Windows：使用 [CRAN Windows](https://cran.r-project.org/bin/windows/base/)；将 `Rscript` 加入 PATH，或把环境变量 `OMICS_RSCRIPT` 设为实际 `Rscript.exe` 的绝对路径。
- Linux：按 [CRAN Linux](https://cran.r-project.org/bin/linux/) 和发行版说明准备 R 及编译库。

`OMICS_RSCRIPT` 指向已安装的可执行文件，不是目录。修改后在新进程中重试 doctor。不要在市场包内填写开发者本机绝对路径。

## R 包

用户已授权准备所需环境时，运行 doctor 返回的 `install_command`，例如：

```bash
monadomics setup-r deg
monadomics doctor --group deg
```

setup-r 检查当前 R 的库路径，把缺失依赖安装到 `R_LIBS_USER` 指定的个人库（含 R 默认的个人库路径），覆盖选定模块和公共依赖。它可能访问 CRAN/Bioconductor 并耗时较长，应等待同一次安装完成再检查；不要放进连接器 init 或并发启动多份 R 安装。安装日志写入 stderr，最终状态以 JSON 写入 stdout。不会自动删除 R 包的安装锁。

需要自定义 R 包库时，由用户的 R 配置或 `R_LIBS_USER` 设置；缺少编译器/系统库时，依据安装日志和平台文档处理。

## 数据与网络

表达矩阵、生存数据和分析产物在本机处理。安装依赖会联网，KEGG 等步骤也可能访问公共资源；离线时不能声称这些请求已成功。不会要求登录、API Key，或把数据送入另一个模型服务；WorkBuddy 本身读取的对话和工具输出遵循宿主的数据处理方式。

默认输出位置为 `~/.workbuddy/workspace/omics`；可用 `--output-dir` 或 `OMICS_OUTPUT_DIR` 指定。输出按 deg、enrich、figures、survival 分类。

## 安装范围

新版使用 `monadomics-analysis` 主 Skill，不再调用旧版 `omics` MCP。GitHub 安装器会备份旧 `omics-analysis` Skill，停用指向旧 `mcp/omics_mcp.py` 的 MCP 配置；其他连接器、既有分析文件和 R 包保留。安装后重启 WorkBuddy 使旧 MCP 停用生效。通过市场迁移时需自行停用旧版本。

安装器的 CLI 检查不代表 R 就绪；R 分析是否可用以 doctor 和实际分析结果为准。
