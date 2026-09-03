# Q2 Code Plan

## Approved methods
- M1 (main): GAMM预测达标时间 + 连续风险函数优化 + 数据驱动BMI分组
- M2 (baseline): 临床经验BMI分组 + 组内统计描述

## Implementation strategy

### Step 1: 复用 Q1 GAMM 模型
- 重新拟合 Q1 的 GAMM: y_conc ~ bs(gest_week, df=4) + bs(bmi, df=3) + (1+gest_week|patient)
- 用固定效应预测: 给定 BMI，求 Y浓度 vs 孕周的群体平均曲线

### Step 2: 求达标时间函数 T*(BMI)
- 对每个 BMI 值，在 GAMM 群体平均曲线上找浓度首次达到 4% 的孕周
- 用二分法或 scipy.optimize.brentq 求解
- 输出: BMI → 达标孕周的映射曲线

### Step 3: 数据驱动 BMI 分组
- 方法: 在 T*(BMI) 曲线上用分段线性回归(piecewise regression)寻找自然断点
- 备选: 基于达标时间差异的层次聚类
- 与临床经验分组 [20,28),[28,32),[32,36),[36,40),40+ 对比

### Step 4: 风险函数构建与最优时点
- R_delay(t) = sigmoid 延迟风险: 越晚风险越高，12周后加速上升
- R_inaccuracy(t, BMI) = P(浓度 < 4% | t, BMI)，由 GAMM 预测 + 残差分布估计
- Total Risk = w1 * R_delay(t) + w2 * R_inaccuracy(t, BMI)
- 最优 NIPT 时点 = argmin_t Total_Risk(t, BMI)，在 [10, 25] 周范围内优化

### Step 5: 检测误差敏感性分析
- 浓度测量误差: 在 4% 阈值上加 ±0.5%, ±1%, ±2% 扰动
- 观察: 达标时间变化、最优时点变化、风险变化

### Step 6: M2 基线
- 用临床经验分组
- 组内统计达标时间的中位数和分位数
- 组内统计风险值

## Outputs
- figures/: 达标时间曲线、分组可视化、风险函数图、敏感性分析图
- tables/: 分组结果表、M1 vs M2 对比表
- metrics/: 核心数值指标
- run_summary.json
