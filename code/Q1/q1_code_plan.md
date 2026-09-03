# Q1 Code Plan

## Approved methods
- **M1 (main)**: GAMM via statsmodels MixedLM + natural spline basis for 孕周 (方案A)
- **M2 (baseline)**: LME via statsmodels MixedLM, linear fixed effects

## Implementation: 方案A — statsmodels MixedLM + spline basis

GAMM 近似策略：用 patsy 的 `cr()` (natural cubic spline) 对孕周做基展开，
将展开后的基列作为固定效应输入 statsmodels MixedLM，随机效应保留 (1+孕周|孕妇)。
这样固定效应部分捕捉非线性，随机效应部分处理面板聚类。

BMI 也用 `cr()` 做样条展开以捕捉可能的非线性（与 GAMM s(BMI) 对应）。

## Data
- Input: `workspace/data_clean/male_cleaned.csv`
- Key columns: 检测孕周_数值, BMI, Y染色体浓度, 孕妇代码

## Scripts
- `q1_main.py` — M1 GAMM (spline + MixedLM) + M2 LME baseline + comparison + diagnostics
- No `run_all.py` needed (single script)

## Outputs
- `results/Q1/experiments/round1/run_summary.json`
- `results/Q1/experiments/round1/figures/` — residual, QQ, partial effect plots
- `results/Q1/experiments/round1/tables/` — coefficient tables
- `results/Q1/experiments/round1/metrics/` — model comparison metrics

## Seeds and reproducibility
- Random seed: 42
- Python version: 3.x, deps: pandas, numpy, statsmodels, patsy, matplotlib, scipy
