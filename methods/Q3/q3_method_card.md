# Q3 Method Card

## Goal and success criteria

在 Q2 基础上，综合考虑身高、体重、年龄等多因素及检测误差，结合 Y 染色体浓度达标比例（浓度 >= 4% 的群体比例），按 BMI 分组给出每组最佳 NIPT 时点，使孕妇潜在风险最小。

成功标准：
- 多因素 GAMM 拟合优于仅含 BMI 的模型（LR 检验显著）
- 每个 BMI 组有明确的最佳检测时点和达标比例
- 蒙特卡洛误差分析给出时点推荐的置信区间
- 与纯 BMI 分组（Q2 思路）形成可比的改进论证

## Human constraints
- Output form: BMI 分组表（区间 + 最佳时点 + 达标比例 + 风险值）+ 误差影响分析
- Priority: 决策式优化（达标比例进目标函数，非事后报告）
- Unacceptable failure: 共线性导致模型不稳定；分组结果无法给出可操作的时点建议
- Experiment budget: 单次运行 <60s，总实验 <=5 轮

## 常规方案与短板分析

| # | 常规方案 | 短板 |
|---|---------|------|
| 1 | Q2 方案直接加入年龄/身高/体重作为额外协变量 | BMI/身高/体重 VIF>100，共线性导致系数不稳定、标准误膨胀 |
| 2 | 用 PCA 降维后建模 | 丢失物理解释性，无法按 BMI 分组 |
| 3 | 仅用 BMI + 年龄（丢弃身高体重） | 题目要求"综合考虑"，丢弃变量不符合题意 |

## 创新路线

### 路线 A：BMI + 残差分解 + 多因素 GAMM（推荐主方案）

- **创新点**：领域知识驱动的共线性处理——BMI = 体重/身高^2 已是身高体重的综合指标，额外提取"体重偏离 BMI 预期的残差"（weight_resid = 体重 - f(BMI, 身高)）作为独立信息；年龄直接纳入；达标比例作为群体级风险函数进入优化目标（决策式非报告式）；蒙特卡洛模拟传播测量误差
- **实现难度**：中
- **优缺点**：保留 BMI 的可解释性和分组功能，同时提取身高/体重的独立信息；残差分解有统计学依据（Frisch-Waugh 定理）。但残差变量的物理解释需要论证
- **论文落地**：强调共线性处理的创新（对比直接纳入 vs 残差分解的 VIF 改善）；达标比例概率化处理；蒙特卡洛误差传播

### 路线 B：BMI 分组 + 组内年龄交互

- **创新点**：在 Q2 基础上加入 BMI x 年龄交互项，不处理身高/体重共线性而是通过分组内分析规避
- **实现难度**：低
- **优缺点**：简单可解释，但未充分利用身高/体重信息
- **论文落地**：作为基线对比方案

## Shortlist

| ID | Role | Mathematical idea | Why eligible | Main risk | Implementation cost |
|----|------|-------------------|-------------|-----------|-------------------|
| M1 | main_candidate | 多因素 GAMM + 残差分解 + 达标比例风险优化 + MC 误差传播：y_conc ~ s(孕周) + s(BMI) + age + weight_resid + (1+孕周\|孕妇)，达标比例 P(y>=0.04\|t,group) 进入风险函数 | Q1 GAMM 已验证有效(AIC=-5240)；残差分解解决 VIF>100；达标比例是题目核心要求；MC 模拟覆盖误差分析 | 残差变量解释性需论证；MC 模拟增加计算时间 | 中 |
| M2 | usable_baseline | 仅 BMI 的 GAMM + 年龄交互：y_conc ~ s(孕周) + s(BMI) + age + BMI:age + (1+孕周\|孕妇)，同样的达标比例风险优化 | 与 M1 直接可比；不需共线性处理；验证多因素是否有实质改进 | 未利用身高/体重独立信息；交互项可能不显著 | 低 |

## Baseline validity
- Real task completed: 是，M2 输出 BMI 分组 + 最佳时点 + 达标比例 + 误差分析
- Comparable output/metric: 是，与 M1 直接比较 AIC/BIC、达标比例、推荐时点差异

## Risk-probe summary

| ID | Executability | Data/assumptions | Degeneracy | Sensitivity | Scale | Verdict |
|----|---------------|------------------|------------|-------------|-------|---------|
| M1 | PASS (收敛 0.69s, lbfgs) | PASS (残差分解后 VIF max=1.01) | PASS (1013 unique fitted / 1082 obs) | PASS (MC 500 iter, CI95 稳定) | PASS (全流程 <30s) | PASS |
| M2 | PASS (收敛 0.71s, lbfgs) | PASS (无共线性) | PASS (1015 unique fitted) | PASS | PASS | PASS |

## Fallback trigger
- Trigger: M1 残差分解后 weight_resid 的 VIF 仍 >10，或 weight_resid 系数 p>0.3 且 LR 检验 M1 vs M2 不显著（p>0.1）
- Evidence to evaluate: VIF 表 + LR 检验 M1 vs M2
- Round 1 状态: **已触发**。weight_resid p=0.685, age p=0.309，均不显著。但 VIF 分解成功（350->1.01），M1 AIC=-5214 优于 M2 AIC=-5209。结论：多因素（年龄、体重残差）对浓度预测无统计显著贡献，但残差分解本身作为方法论创新有效，且 BMI 分组 + 达标比例优化 + MC 误差分析框架完整可用。

## Compact history
- 2026-09-03: 初始方法筛选。BMI/身高/体重 VIF>100 需残差分解；年龄作为独立协变量；达标比例进入风险优化目标；MC 误差传播覆盖题目要求。
- 2026-09-03 round1: M1/M2 均收敛。VIF 分解成功（350->1.01）。weight_resid p=0.685, age p=0.309 均不显著——多因素对浓度无额外贡献，BMI 已充分捕捉体成分效应。BMI 分 4 组: [20,32) 最佳 12.8w(75.3%), [32,36) 13.3w(74.8%), [36,40) 13.3w(67.0%), [40,+inf) 24.8w(95.0%)。MC 误差分析显示前两组稳定（std<1w），高 BMI 组不确定性大（std=3-5w，样本量仅 4-16）。
