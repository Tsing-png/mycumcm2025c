# Q1 Method Card

## Goal and success criteria

建立男胎 Y 染色体浓度与孕周、BMI 等指标的关系模型，检验显著性。模型需正确处理纵向面板结构（267 人 × ~4 次测量），输出固定效应系数及 p 值，并为 Q2/Q3 提供可用的浓度预测基础。

成功标准：主要自变量（孕周、BMI）显著（p<0.05）；模型拟合优于忽略聚类的 pooled OLS；残差诊断无严重违反。

## Human constraints
- Output form: 固定效应系数表 + 显著性检验 + 模型拟合指标
- Priority: 解释性优先（描述性/推断为主，framing decision q1_framing_descriptive_vs_predictive = descriptive/inference）
- Unacceptable failure: 忽略重复测量结构导致标准误低估、p 值虚假显著
- Experiment budget: 单次运行 <30s，总实验 ≤5 轮

## 常规方案与短板分析

| # | 常规方案 | 短板 |
|---|---------|------|
| 1 | Pooled OLS 多元回归 | 忽略个体聚类（ICC=0.60），标准误严重低估，p 值不可信 |
| 2 | GEE（广义估计方程） | 处理聚类但仅估计群体平均效应，无法分离个体变异；小群组数（267）下工作相关矩阵估计不稳 |
| 3 | 固定效应面板模型 | 消除个体异质性但同时消除不随时间变化的协变量（BMI 近似不变），无法估计 BMI 效应 |

## 创新路线

### 路线 A：领域知识约束的混合效应模型（推荐主方案）

- **创新点**：将产科领域知识注入模型结构——孕周效应允许个体异质斜率（随机斜率），BMI 通过体脂稀释机制影响基线浓度（随机截距）；BMI/体重/身高的高度共线性（VIF>100）通过仅保留 BMI 作为体成分综合指标来处理，领域合理性约束要求孕周系数为正、BMI 系数为负
- **实现难度**：低（statsmodels/R lme4 直接支持）
- **优缺点**：正确处理纵向结构；可分离群体平均趋势与个体变异；直接输出显著性检验。局限是线性假设可能遗漏非线性效应
- **论文落地**：强调「面板数据不能用普通回归」的方法学论证 + ICC 证据 + 随机斜率的生物学解释（不同孕妇胎儿 DNA 释放速率不同）

### 路线 B：分段/样条混合效应模型

- **创新点**：在路线 A 基础上，用限制性三次样条（restricted cubic spline）捕捉孕周与浓度的非线性关系，将领域知识（孕早期浓度上升快、中后期趋缓）编码为样条节点位置选择的先验
- **实现难度**：中（需选择节点数和位置，模型自由度增加）
- **优缺点**：能捕捉非线性但增加过拟合风险；节点选择需要领域依据
- **论文落地**：与线性 LME 做嵌套似然比检验，量化非线性改进幅度

## Shortlist

| ID | Role | Mathematical idea | Why eligible | Main risk | Implementation cost |
|----|------|-------------------|-------------|-----------|-------------------|
| M1 | main_candidate | 广义加性混合模型（GAMM）：Y浓度 ~ s(孕周) + s(BMI) + (1+孕周\|孕妇), 样条拟合非线性+随机效应 | ICC=0.60 确认强聚类；孕周与浓度关系可能非线性（孕早期上升快、中后期趋缓）；GAMM 用数据驱动的样条自动捕捉非线性，避免强制线性假设；领域知识支持个体异质斜率 | 样条节点选择影响拟合；计算量略大于 LME；过拟合风险需交叉验证控制 | 中 |
| M2 | usable_baseline | 线性混合效应模型（LME）：Y浓度 ~ 孕周 + BMI + (1+孕周\|孕妇), 随机截距+随机斜率 | probe 已验证收敛（LL=2594），孕周(p<0.001)和BMI(p=0.025)均显著；与 GAMM 对比可量化非线性改进幅度 | 线性假设可能遗漏非线性效应；残差 skew=0.58 | 低 |
| M3 | conditional_fallback | Beta GLMM（logit 链接）：logit(Y浓度) ~ s(孕周) + BMI + (1+孕周\|孕妇) | Y浓度在(0,1)有界且右偏(skew=0.70)，Beta 族天然适配 | 收敛困难；解释性降低 | 中高 |

## Baseline validity
- Real task completed: 是，Pooled OLS 输出系数+显著性检验+R²，完成关系建模任务
- Comparable output/metric: 是，与 M1 直接比较系数估计、标准误、AIC/BIC

## Risk-probe summary

| ID | Executability | Data/assumptions | Degeneracy | Sensitivity | Scale | Verdict |
|----|---------------|------------------|------------|-------------|-------|---------|
| M1 | PASS (收敛, LL=2594) | CONDITIONAL (残差非正态 skew=0.58, BMI/体重 VIF>100 需排除体重) | PASS (系数非零, 随机效应方差>0) | PASS (BMI±1 扰动系数稳定) | PASS (<1s) | CONDITIONAL |
| M2 | PASS (OLS 无收敛问题) | CONDITIONAL (忽略 ICC=0.60 的聚类, 标准误低估) | PASS (R²=0.046 偏低但系数非零) | PASS | PASS | CONDITIONAL |

## Fallback trigger
- Trigger: M1 的 LME 残差 QQ 图严重偏离 + Shapiro p<0.001 + 残差方差随孕周系统性变化（异方差），且非线性诊断显示孕周二次项或样条项显著改善拟合（似然比检验 p<0.01）
- Evidence to evaluate: 残差诊断图 + 似然比检验 M1 vs M3

## Compact history
- 2026-09-03: 初始方法筛选完成。ICC=0.60 确认面板结构；LME probe 通过（孕周 p<0.001, BMI p=0.025）；Pooled OLS 作为 baseline（R²=0.046）；Beta GLMM 作为条件 fallback。
