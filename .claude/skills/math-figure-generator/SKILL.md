---
name: math-figure-generator
description: Generate and render-verify publication-quality mathematical-modeling figures from saved evidence, using the approved figure plan, source data, claim, type, and consistent visual system.
---

# Preconditions

- Figure type, source artifacts, and target claim are known.
- Type 3 claim is human-confirmed.
- Submission figures use final/frozen evidence.

# References

Load only what the requested chart needs:

- `references/chart-patterns.md`
- `references/color-systems.md`
- `references/layout-guide.md`
- `references/render-check.md`

# Workflow

1. Verify source files and the exact variables/units to plot.
2. Choose the smallest chart form that communicates the claim.
3. Generate with deterministic code, preferably matplotlib.
4. Save editable source code and the requested output format.
5. Apply the shared color, typography, sizing, and labeling conventions.
6. Render the final output and inspect it visually.
7. Check clipping, overlap, illegible text, misleading axes, legends, empty panels, and source/claim mismatch.
8. Iterate until render checks pass.

# Output Locations

- Type 1/2 exploration: `results/Qx/experiments/roundN/figures/`
- Type 3/4 submission: `paper/figures/`

Use stable descriptive filenames. Do not copy Type 1 diagnostics into the paper directory.

# Figure Requirements

- Labels include units where applicable.
- Captions state what is shown and the evidence-backed takeaway without overstating causality.
- Baseline is visually distinct but not exaggerated.
- Uncertainty is shown when it is part of the claim.
- Type 3 raster output is at least 300 dpi; vector output is preferred when compatible.
- Accessibility and grayscale differentiation are considered.

# Rules

- Do not fabricate or manually alter plotted values.
- Do not use a chart type that hides concentration, uncertainty, or negative results.
- Do not truncate axes misleadingly.
- Do not create decorative 3D effects.
- Do not treat code execution as render verification.
- Keep diagnostic and paper roles separate.

## 图表与配色规范（强制）

### 图

1. **图内无标题**：图上方、图下方都不写标题，也不画面板编号，只保留数据；标题信息放图注。
2. **图题在图下方**（表题在表上方）。
3. **字号按印刷尺寸**：图内文字 8.5–11pt（最小 7.5pt），禁止默认字号整图缩放导致图大字小。
4. **布局不默认 2×2**：按面板数量与证据权重选 1×2、2×1、3×1、主图+辅助等。
5. **图型多样性**：同一类型图（折线/柱状/散点/热力图/箱线等）全篇 ≤3 张，同类图不连续 3 张。
6. **导出 300 DPI PNG**（矢量/灰度预览按需追加）。

### 配色（7 套期刊色板，默认 nature）

同一篇论文全文用同一套色板；严格色盲友好用 okabe-ito 或 wong。

| 色板 | 默认主色（按使用顺序） | 说明 |
|---|---|---|
| nature | `#0F6BBD` `#F26F21` `#9ECAE1` `#083C5F` `#7F8C8D` | 蓝主橙强调，2–3 组对比 |
| science | `#E41937` `#A31A30` `#2F6EB5` `#8FC1E1` `#8F8F8F` | 红为绝对主角 |
| cell | `#1E8449` `#186138` `#82E0AA` `#8E44AD` `#566573` | 绿系 + 低饱和 |
| ieee | `#0072BD` `#D95319` `#EDB120` `#7E2F8E` `#77AC30` `#4DBEEE` `#A2142F` | IEEE 常用 |
| okabe-ito | `#E69F00` `#56B4E9` `#009E73` `#F0E442` `#0072B2` `#D55E00` `#CC79A7` `#999999` | 色盲友好黄金标准 |
| wong | `#0077BB` `#EE7733` `#33BBEE` `#CC3311` `#009988` `#BBBBBB` | 色盲安全 |
| general | `#1f77b4` `#ff7f0e` `#2ca02c` `#d62728` `#9467bd` `#8c564b` `#e377c2` `#7f7f7f` | 通用 |

配色规则：
- 单张图类别色 ≤5，超过时拆图或合并次要类别；
- 不用颜色作唯一区分，配合线型/标记冗余编码（色盲与灰度打印均可读）；
- 连续数据用 viridis/magma，发散数据用 RdBu_r 且 0 值居中；
- 语义色全文一致（异常红、正常绿、警告橙黄）；
- 导出灰度预览，确认灰度下仍可区分。

### 表（三线表）

只用顶线、栏目线、底线；不画竖线、斜线；表题在上；单位标注在表头/列名，不重复写进单元格。

### 公式

公式独占一行居中、编号右对齐（`(1)(2)…` 连续）；正文用 `（n）` 引用；公式出现后解释每个符号。

## 知识库参考

绘图时按需读取 `knowledge/references/绘图参考/`、`knowledge/assets/期刊配色方案.md`。

# Verification

- Source, claim, type, and target section agree.
- Render inspection passed.
- Text is readable at final paper size.
- Legends, colors, markers, axes, units, and captions are consistent.
- Final output path exists and is recorded in the figure plan.
