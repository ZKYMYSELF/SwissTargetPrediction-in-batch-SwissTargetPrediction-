# SwissTargetPrediction Windows Release

这是一个面向 Windows 的 SwissTargetPrediction 批量预测脚本。

## 文件说明

- `gui.py`: 图形化界面入口
- `swisstargetprediction_gui.py`: GUI 主程序
- `swisstargetprediction_batch.py`: 网站提交与结果解析核心逻辑
- `run_gui.bat`: 双击启动脚本
- `requirements.txt`: Python 依赖清单

## 使用前准备

- 安装 Python 3.9 或更高版本
- 确认安装 `requests`
- 如果想让双击直接启动，Windows 上建议保留 `py` 启动器，或者让 `pythonw` 在环境变量里可用

安装依赖：

```bash
pip install -r requirements.txt
```

## 启动方式

1. 双击 `run_gui.bat`
2. 或者运行：

```bash
python gui.py
```

## 输入方式

- 左边输入化合物名称
- 右边输入对应的 `SMILES`
- 两边数量必须一致，否则会报错
- 支持从 Excel 直接复制两列文本后，点击“从 Excel 剪贴板导入两列”

## 输出规则

- 每个化合物单独保存一个 `CSV`
- 文件名格式：

```
序号_化合物名称_时间戳.csv
```

- 表格 `A` 列为化合物名称
- 默认输出到程序目录下的 `输出结果` 文件夹

## 运行策略

- 程序会一个一个提交到 SwissTargetPrediction
- 如果某个化合物失败，只跳过该文件，继续处理后面的化合物
- 失败信息会显示在预览表和运行日志里

## 注意事项

- SwissTargetPrediction 可能会限制频繁请求，所以程序默认逐个提交，并保留间隔时间
- 如果网站返回错误页，程序会尽量抓取网页中的可读提示
- 如果担心IP被网站封锁，您可以使用代理而非自己的真实IP访问；如果很不幸被封了，可以尝试清除cookies或者直接联系SwissTarget官方，他们会在24小时之内回复你

## 开发者信息

KY Z  
E-mail：zkymyself@gmail.com