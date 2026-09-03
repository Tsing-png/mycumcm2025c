---
name: method-selector
description: Build and risk-screen a compact role-based method shortlist for a mathematical-modeling subquestion. Use after problem framing and data profiling, before model code generation, to propose a main candidate, a usable baseline, and at most one conditional fallback without padding the pool.
---

# Purpose

Convert the framed problem and data profile into a small executable decision surface. Screen methods for load-bearing data, assumption, degeneracy, sensitivity, and scale risks before asking the human to choose.

This skill proposes and probes methods. The human chooses the method.

## 方法论创新规约（强制）

> 角色：方法论创新顾问。硬约束：底层求解器、基础算法固定不改，创新只发生在「建模层」，不碰「求解层」。生成创新路线供人工选择，仍遵循本 skill 的 shortlist 规则（主方案 + baseline + 至多一个 fallback），创新路线不全部实现。

### 一、创新方向（每次至少命中一个）

1. 模型结构改进：改变量定义、分层/解耦方式、分解结构。
2. 多模型融合创新：不同假设空间的模型集成/级联/迁移。
3. 目标函数 / 约束重构：隐含目标显式化、软约束写惩罚项、领域规律写正则或约束。
4. 评价指标体系创新：设计更贴问题本质的指标。
5. 领域知识适配改良：领域规律改造成先验、约束、特征或损失项。

### 二、创新方法论（先找点，再落笔）

1. 问题结构解剖：数据特性 / 约束特性 / 领域特性 / 规模特性。
2. 领域知识注入（创新第一来源）：追问「领域规律能不能变成方法的一部分，而非事后解释」。
3. 方法家族遍历：精确法 / 统计法 / 机器学习 / 仿真启发式，逐类写「为什么用/不用」。
4. 每个选择绑定「为什么」；记录试错（试过什么、为何放弃）。

### 二.5、领域知识注入的四种具体模式

把领域规律「变成方法的一部分」，有四种可落地模式（适用于物理、化学、农学、经济、工程等任何领域）：

1. **先验正则化**：把领域规律写成目标函数的正则/惩罚项，使模型在训练或求解时就受领域规律约束，而非事后解释。
2. **跨域/跨类迁移**：用可迁移的结构参数（如缩放系数）把在 A 类/域上训练的模型适配到 B 类/域，允许参数取负以处理方向反转。
3. **条件独立/图结构建模**：用精度矩阵、图模型等过滤数据固有结构引入的伪相关，提取直接关联。
4. **领域合理性约束**：把领域合理性写成对结果的检验或硬约束（如某物理量应落在某区间），不满足则标记为不合理。

选型时主动对照这四种，看领域规律能落到哪一种。

### 二.6、报告式 vs 决策式自查（强制）

凡问题要求「考虑 X / 权衡 X / 优化 X」（X 可为风险、成本、公平、稳健性、可持续性、多目标等），必须自查：

- ❌ **报告式**：先优化单一目标（如期望收益），求出方案后再报告 X 的数值。X 未影响方案选择。
- ✅ **决策式**：X 进入目标函数（如目标写成「收益 − 风险惩罚项」）或约束，X 的取值直接影响最终方案。

判定标准：删掉 X 的评估环节后，若方案完全不变，则 X 是「报告式」，尚未真正优化 X。

典型反例：「期望收益最大化 + 事后计算尾部风险」——风险被报告、但未参与决策；正确做法是把风险度量写入目标函数，让风险厌恶程度影响方案。

### 三、输出格式（写进方法卡/论文）

1. 列 2–3 种大众常规方案。
2. 逐条分析常规方案短板。
3. 给 2 套差异化创新路线。
4. 每条路线写：创新点 / 实现难度 / 优缺点 / 论文落地写作角度。

### 四、铁律

- 反羊群：不照搬网上范文方案；「大家都会做」的不算创新。
- 方法论层面：创新必须改模型/方法/指标本身，禁止只停留在可视化、后处理、调参。
- 领域知识进结构：领域规律只出现在「结果解读」而没进模型结构 = 没用上。
- 可验证：每条创新要能被对比实验或消融实验验证。

## 算法选型参考（五大题型速查）

### 五大题型识别

拿到题目先问：输出是一个数（优化）、一条序列（预测）、一个排名（评价）、一个类别（分类/聚类）还是一组微分方程（机理）？混合题型按主问题定主模型。

| 题型 | 核心问法 | 典型场景 |
|---|---|---|
| 优化 | 求最大/最小、资源分配 | 路径、调度、选址、装箱 |
| 预测 | 预测未来数值/趋势 | 销量、人口、价格 |
| 评价 | 多方案排序/打分 | 评级、比选 |
| 分类/聚类 | 判类别/自动分组 | 风险识别、分群 |
| 机理/动力学 | 系统随时间演化 | 传染病、生态、物理仿真 |

### 选型速查

- 预测：数据量 <15 → 灰色预测 GM(1,1)；纯时序 → ARIMA/指数平滑；多因子 → 回归/随机森林/XGBoost；非线性强 → LSTM；需不确定性 → Bootstrap。
- 评价：专家经验 → AHP；客观数据 → 熵权法；方案距离 → TOPSIS；模糊语言 → 模糊综合；效率 → DEA；指标高度相关 → PCA+TOPSIS。
- 优化：线性 → LP；含整数/0-1 → MIP；非线性 → 启发式；多目标 → 加权或 NSGA-II；大规模 → 遗传/模拟退火（固定种子、多次运行）。
- 分类/聚类：小样本可解释 → 逻辑回归/SVM；追求精度 → XGBoost/RF；不知 K → DBSCAN；指定 K → K-means。
- 机理：种群增长 → Logistic；传染病 → SIR/SEIR；捕食 → Lotka-Volterra；力学 → 牛顿方程。

### 防错速查

- 优化：不漏非负/整数/容量约束；连续松弛取整后重新验证可行性；非凸/启发式做多起点或多种子。
- 预测：时序按时间划分不随机打乱；scaler 只在训练集 fit；测试集不参与调参。
- 评价：正负向指标先统一方向；AHP 查一致性；TOPSIS 说明无量纲化并做权重敏感性。
- 图论：明确有向/无向；负权边不用 Dijkstra；TSP/VRP 防子回路。
- 统计/ML：类别不平衡用 F1/AUC 而非 Accuracy；特征工程有业务含义。

# Preconditions

- G1 problem framing passed.
- Required output and evaluation criteria are known.
- Relevant data inventory or audit exists.
- `planning/symbol_table.md` and `planning/model_assumptions.md` exist when the problem needs them.

If these are missing, return to the producer skill rather than guessing.

# Inputs

- Problem parse and classification.
- Data audit, including missingness, effective sample size, imbalance, cardinality, and distribution summaries.
- Literature analysis when available.
- Contest deadline, implementation language, interpretability needs, and compute limits.
- `planning/session_config.json`.
- Existing `methods/Qx/qx_method_card.md` and decision ledger when revising.

# Workflow

1. **Align the decision surface.**
   - Invoke `decision-prompt-builder` before generating an open-ended shortlist.
   - Ask about human-owned trade-offs, not algorithm names.
   - Reuse answers already present in the decision ledger.

2. **Derive method requirements.**
   - Start from required output, hard constraints, data characteristics, validation criteria, explanation burden, and experiment budget.
   - Identify the failure modes that would make a method unusable.

3. **Create a role-based shortlist.**
   - One `main_candidate`: best fit to the chosen trade-off.
   - One `usable_baseline`: completes the real task and yields directly comparable outputs.
   - At most one `conditional_fallback`: differs in a meaningful mathematical way and has an explicit activation trigger.
   - If a simple reference cannot complete the real task, label it `diagnostic_reference`; it does not satisfy the baseline requirement.
   - Do not add a method merely to reach a candidate count.

4. **Define method-specific risk checks.**
   - Use the contract in `references/risk-probe-contract.md`.
   - Select only relevant assumption checks.
   - Always check output degeneracy or concentration with metrics appropriate to the output.
   - Bound probe runtime rather than source-line count.

5. **Run the risk probe on the main candidate and usable baseline.**
   - Use a representative slice or full-data diagnostic as appropriate; never rely only on the first rows.
   - The probe may use reusable scripts and may save detailed metrics, but its canonical output is one compact summary.
   - Probe the fallback only enough to establish that its trigger and risk profile are credible. Do not fully implement it.

6. **Write canonical artifacts.**
   - `methods/Qx/qx_method_card.md`
   - `methods/Qx/probes/risk_probe_summary.json`
   - Update `planning/manifests/Qx.json` if present.

7. **Ask for the method choice.**
   - Present the probe evidence through a choice card.
   - After the user answers, hand the exact answer to `modeler-decision-logger` for append-only capture in `methods/Qx/qx_decisions.jsonl`.
   - If no answer is available, stop. Do not create a placeholder decision file.

# Method Card Contract

`qx_method_card.md` stays compact and contains:

```markdown
# Qx Method Card

## Goal and success criteria

## Human constraints
- Output form:
- Priority:
- Unacceptable failure:
- Experiment budget:

## Shortlist
| ID | Role | Mathematical idea | Why eligible | Main risk | Implementation cost |

## Baseline validity
- Real task completed:
- Comparable output/metric:
- If no, classification: diagnostic_reference

## Risk-probe summary
| ID | Executability | Data/assumptions | Degeneracy | Sensitivity | Scale | Verdict |

## Fallback trigger
- Trigger:
- Evidence to evaluate:

## Compact history
- One line per material change, with decision_id when human-owned.
```

Do not maintain a separate iteration log for new work.

# Probe Verdicts

- `PASS`: eligible for the human choice.
- `CONDITIONAL`: eligible only with a stated mitigation or fallback trigger.
- `FAIL`: not offered as a selectable main or baseline.

A method fails screening when a load-bearing assumption fails, the output degenerates, it cannot produce a legal result, or its cost violates the user's budget. A method does not fail merely because an irrelevant generic diagnostic is unavailable.

# Output and Handoff

After G2 screening:

- If the human choice is absent: return the evidence-backed choice card.
- If G2.5 is decided: hand the method card, probe summary, chosen IDs, and experiment budget to `model-code-analyzer`.
- Instruct code generation to implement only the approved main method and usable baseline.
- Keep the fallback dormant until its recorded trigger fires.

# Rules

- Do not use a fixed candidate count.
- Do not use source-line count as validation quality.
- Do not invent missing data fields, constraints, labels, or evaluation metrics.
- Do not call a nonfunctional toy method a baseline.
- Do not fully implement all shortlisted methods.
- Do not select the method or write the human rationale.
- Keep AI suggestions visibly separate from the human decision.

# Compatibility

When revising an older workspace, read:

- `methods/Qx/qx_method_candidates.md`
- `methods/Qx/qx_method_iteration_log.md`
- `methods/Qx/poc/`

Migrate material evidence into the method card and probe summary. Do not require new legacy PoCs or iteration logs.

# References

- Risk checks and summary schema: `references/risk-probe-contract.md`
- Method-family routing cues: `references/method-family-guide.md`

## 知识库参考

选型时按需读取 `knowledge/references/模型选型树.md`、`算法索引.md`，以及 `knowledge/assets/01-07 算法说明`。

# Verification

- Shortlist contains a main candidate and a genuinely usable baseline.
- Optional fallback has a concrete trigger.
- Main and baseline have evidence-backed probe verdicts.
- Output-degeneracy checks are present.
- Method card and probe summary exist.
- No per-skill pending decision file was created.
- No code-generation handoff occurs before a human method choice is recorded.
