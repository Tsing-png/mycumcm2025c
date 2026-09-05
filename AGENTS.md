# Core Philosophy

- **The AI owns mechanical correctness; the human owns modeling judgment.**
- Start from goals, objects, constraints, data, outputs, variables, relationships, and checkable conclusions.
- Do not start from model names or favorite techniques.
- Separate assumptions, observations, derivations, and validated conclusions.
- Preserve evidence that changes a decision; do not create files merely to prove that a skill ran.

# Configuration

`planning/session_config.json` has two independent controls:

```json
{
  "interaction_mode": "learning",
  "rigor_profile": "lean"
}
```

- `interaction_mode`: `learning` or `speed`. It changes question density and when AI suggestions are shown.
- `rigor_profile`: `lean` or `submission`. It changes artifact and audit density, never the human-judgment boundary.
- Default to `learning + lean` in a fresh workspace.
- Use `lean` while exploring and iterating. Switch to `submission` only when preparing writer handoff or final assembly.
- For compatibility, read legacy `{ "mode": "learning" | "speed" }` as `interaction_mode`.

# Repository Skill Copies

- `.codex/skills/` and `.claude/skills/` are two complete, independently usable skill trees.
- Every skill and referenced local resource required at runtime must exist in both trees; neither tree may depend on a wrapper, symlink, or path into the other.
- When a shared skill contract changes, update and validate both copies in the same change.
- Runtime-specific wording may differ only when necessary, but each copy must remain standalone and behaviorally consistent with this policy.
- `plugins/mathmodeling-skills/skills/` is the generated distribution copy used by both native plugin manifests. After the two standalone trees agree, refresh it with `scripts/sync-plugin.sh` and verify it with `scripts/sync-plugin.sh --check`.
- Keep both plugin manifests and both marketplace catalogs aligned for every release. Bump the version in both plugin manifests and the Claude marketplace entry together.

# Workflow Discipline

- Parse before classifying; classify before screening methods.
- Ask the modeler about output form, priority, unacceptable failure, and experiment budget before creating a method shortlist.
- Build a role-based shortlist rather than filling a quota:
  - one `main_candidate`;
  - one `usable_baseline`;
  - at most one `conditional_fallback`.
- Allow only a main candidate plus baseline when no genuine fallback exists.
- A trivial reference that cannot complete the real task is `diagnostic_reference`, not a baseline.
- Fully implement the human-approved main method and usable baseline only. Activate a fallback only when its recorded trigger fires.
- Keep changes minimal, traceable, and reviewable.

# Human Decision Convention

Human decisions are captured in one append-only ledger per subquestion:

`methods/Qx/qx_decisions.jsonl`

Use `planning/framing_decisions.jsonl` for global or pre-subquestion framing decisions made before a Qx method directory exists.

Each line is a JSON object with at least:

```json
{
  "decision_id": "q2_method_choice",
  "decision_type": "method_choice",
  "status": "DECIDED",
  "decided_by": "human",
  "captured_in_mode": "learning",
  "choice": "M2",
  "rationale": "M2 is selected because ...",
  "evidence_refs": ["methods/Q2/probes/risk_probe_summary.json"],
  "decided_at": "ISO-8601 timestamp"
}
```

- The AI may present evidence and options but must not originate the human's choice, rationale, confidence, physical interpretation, or submission authorization.
- The AI may append the user's answer verbatim or faithfully structure it; it must not strengthen or invent the rationale.
- Do not create per-skill `*_modeler_decision.md` files for new work.
- Existing decision Markdown files remain readable during migration but are not required for new work.
- A decision passes only when it is human-authored, evidence-linked, non-empty, and contains no placeholder.

# Choice Cards

Use choice cards only at modeling-judgment points, normally twice per subquestion:

1. Before method screening: output form, interpretability/performance priority, unacceptable failure, experiment budget.
2. After the first meaningful experiment: proceed, adjust, or activate the fallback.

An optional third card may be used before final freeze for claim scope and confidence. Do not ask users to decide mechanically checkable matters.

# Workflow Gates

## G1 — PROBLEM_FRAMED

- Parse, classification, data inventory, success criteria, and human framing exist.

## G2 — METHOD_SCREENED

- `methods/Qx/qx_method_card.md` defines the main candidate, usable baseline, and optional conditional fallback.
- `methods/Qx/probes/risk_probe_summary.json` exists.
- The main candidate and usable baseline pass the applicable risk checks.
- Any fallback has an explicit activation trigger.
- No fixed candidate count or source-line limit is used.

## G2.5 — METHOD_CHOSEN_BY_HUMAN

- `qx_decisions.jsonl` contains a `DECIDED` human `method_choice` record citing probe evidence.
- Code generation is allowed only when G2 and G2.5 both pass.

## G3 — CODE_AND_EXPERIMENT_REVIEWED

- The approved main method and usable baseline ran.
- `results/Qx/experiments/roundN/run_summary.json` records configuration, seed, metrics, outputs, and failures.
- A language review artifact contains the required named checks:
  - `syntax`
  - `input_contract`
  - `method_alignment`
  - `reproducibility`
  - `output_contract`
- New review artifacts use `code/Qx/reviews/qx_<lang>_review.json`. Legacy Markdown reviews may be read during migration.

## G4 — RESULTS_JUDGED_AND_FROZEN

- The human decision ledger contains result, stability, and claim-scope verdicts tied to computed evidence.
- Final result analysis and robustness report exist.
- In `submission` profile, the solution package and immutable `frozen_numbers.json` exist and are current.

## G5 — PAPER_SECTION_READY

- The writer uses the solution package as the primary source.
- Numerical claims come from `frozen_numbers.json`.
- Physical/domain interpretation and contribution claims are human-confirmed.
- Every paper figure passes render verification.

## G6 — FINAL_AUDIT_PASSED

Run only in `submission` profile. All three must pass:

- cross-media consistency;
- semantic completeness;
- final quality assurance.

# Risk Probe Contract

The risk probe replaces universal ≤30-line PoCs. It is time-bounded, method-specific, and may use reusable scripts.

`methods/Qx/probes/risk_probe_summary.json` must contain:

- `executability`: can the method produce a legal result?
- `data_coverage`: missingness, effective sample size, imbalance, cardinality, and distribution coverage.
- `assumption_checks`: only checks relevant to the method, such as stationarity, multicollinearity, identifiability, clusterability, or constraint feasibility.
- `output_degeneracy`: variance, unique-output count, top-k mass, entropy/Gini, score or rank concentration, and constraint slack where applicable.
- `perturbation_sensitivity`: response to a small justified perturbation.
- `scale_check`: runtime and memory at representative sizes.
- `verdict`: `PASS`, `CONDITIONAL`, or `FAIL`, with evidence and fallback trigger when conditional.

Do not reject a method merely because an irrelevant generic test is unavailable. Do reject or condition it when a load-bearing assumption fails or its output degenerates.

# Lean Artifact Contract

During exploration, keep only:

```text
planning/session_config.json
planning/framing_decisions.jsonl       # only when global framing decisions exist
planning/manifests/Qx.json
methods/Qx/qx_method_card.md
methods/Qx/qx_decisions.jsonl
methods/Qx/probes/risk_probe_summary.json
results/Qx/experiments/roundN/run_summary.json
```

- `planning/manifests/Qx.json` is the machine-readable state source.
- Derive dashboards from manifests; do not rewrite a large dashboard after every state transition.
- `qx_method_card.md` contains roles, assumptions, risks, fallback triggers, and a compact decision history. Do not maintain a separate iteration log for new work.
- Successful runs store summaries and artifact paths. Persist full console logs only for failures or when needed to reproduce an anomaly.
- Ordinary rounds do not require a Markdown experiment report. Generate one only at a human decision point or for the final round.

# Submission Artifact Contract

Before writer handoff, add:

```text
methods/Qx/qx_final_method_explanation.md
code/Qx/reviews/qx_<lang>_review.json
results/Qx/reports/qx_final_result_analysis.md
robustness/Qx/qx_robustness_report.md
results/Qx/reports/qx_solution_package_for_writer.md
results/Qx/reports/frozen_numbers.json
```

The three critical writer rules remain:

1. No final method explanation, no paper section.
2. No final result analysis, no writer handoff.
3. The writer reads the solution package rather than guessing from scattered results.

# Change Impact and Auditing

Classify a change before auditing:

- `NONE`: scratch files, formatting, comments, non-semantic documentation. No consistency audit.
- `LOCAL`: exploratory code or method-card changes before freeze. Run local tests/review only.
- `CANONICAL`: data schema/units, symbols, equations, parameters, official result values, or figure paths. Run a scoped consistency check for affected Qx.
- `FROZEN`: anything that can change a frozen number or paper claim. Log the thaw, update the canonical source, rerun affected experiments, re-freeze, then run scoped consistency.

Do not run a full-workspace audit merely because multiple files changed. Always run the full three-auditor layer once in `submission` profile before final assembly.

# Frozen Numbers

- Numbers flow code → results → freeze → paper.
- Never edit `frozen_numbers.json` by hand.
- To change a frozen value: **解冻 → 修改 canonical source → 重跑 affected work → 重冻结**.
- Record the reason in `results/Qx/reports/freeze_change_log.md`.
- A freeze is stale when a referenced canonical source is newer than `frozen_at`.

# Experiment Output

Every executed round writes:

```text
results/Qx/experiments/roundN/
├── figures/
├── tables/
├── metrics/
└── run_summary.json
```

Create `logs/` only when a failure, warning, or reproducibility need justifies it.

`run_summary.json` records question, round, approved methods, role, status, inputs, outputs, metric summary, seed, environment, warnings, and fallback-trigger state.

# Modeling and Coding Rules

- Match methods to output, data, interpretability, time, and contest constraints.
- Do not choose complexity for appearance.
- Do not invent data, assumptions, evidence, results, or references.
- Keep assumptions explicit and distinguish necessary from simplifying assumptions.
- Maintain `planning/symbol_table.md`; define every symbol and unit before use.
- Use fixed random seeds.
- Save formal outputs to files; console output alone is not a deliverable.
- Keep raw data untouched under `workspace/data_raw/`; write cleaned copies under `workspace/data_clean/`.

# Figures and Paper

- Type 1 diagnostic: internal only.
- Type 2 comparison: may appear in paper.
- Type 3 paper: must support a main claim and pass publication-quality render checks.
- Type 4 appendix: supplementary and referenced from the main text.
- Paper claims must remain proportional to tested evidence.
- Mention eliminated methods only when the record helps explain a real trade-off; do not manufacture breadth.

# Verification

- In `lean`, verify the current gate and only the affected artifacts.
- In `submission`, verify all required artifacts, frozen-number lineage, figure rendering, references, and the three independent audits.
- A review or audit passes by completing its named semantic checks, not by reaching an arbitrary bullet count.
- Flag uncertainty and blocking issues explicitly.
- Do not approve final assembly while any G6 auditor fails.

# 论文写作避免出现

- 一、语言层面（最直观）

  1. 句式过于规整、排比化、书面腔过重
     AI 习惯用长难句、对称句式，比如 “综上所述，本模型具备较强的可行性与鲁棒性，能够为实际决策提供理论支撑与数据参考”，人工写会更口语化、有侧重，不会通篇都是套话。
  2. 空话套话密集，无具体细节
     高频模板句：“具有重要的现实意义”“效果显著提升”“为后续研究奠定基础”，没有绑定本题的具体数据、参数、场景。
  3. 用词过于精致学术，脱离工程实际
     人工写数模会直白说 “迭代 100 次收敛”，AI 会写成 “经由多轮迭代优化，模型实现了稳定收敛，误差满足预设精度要求”。
  4. 段落结构高度同质化
     每个小节开头结尾都是固定模板，引言、模型假设、模型求解的行文节奏完全一致，没有人工写作的自然起伏。
- 二、逻辑与内容层面（评委最敏感）

  1. 泛泛而谈，缺少本题专属细节
     AI 生成的内容可以套用到任何数模题，没有针对本题的变量、约束、数据特点做针对性分析，比如 pareto 优化只说原理，不说本题目标函数怎么构建、约束怎么设定。
  2. 因果逻辑薄弱，跳跃感强
     直接给出结论，缺少推导过程、试错过程，比如 “采用多目标规划求解”，不解释为什么选这个方法、对比了哪些其他方法、舍弃的原因。
  3. 数据和分析脱节
     罗列图表数据，但解读部分高度模板化，不会结合数据异常点、波动情况做针对性解释，人工写作会盯着数据细节分析。
  4. 创新点空泛
     AI 编的创新点大多是 “改进 XX 算法、引入智能优化策略”，没有具体的改进点、量化的提升指标，显得很虚。
- 三、格式与排版层面

  1. 段落长度高度均匀
     AI 生成的段落长短基本一致，人工写作会有长短错落，重点段落会写得更详细。
  2. 参考文献、公式引用生硬
     公式堆砌，符号说明敷衍，参考文献格式标准但和正文内容匹配度低，很多是 AI 随便匹配的文献。
  3. 摘要、关键词模板化严重
     摘要结构固定：问题分析→模型构建→算法求解→结果验证→结论，话术高度统一，缺少本题独有的研究亮点。
- 四、数模场景专属 AI 痕迹

  1. 算法部分只讲通用原理，不写本题的参数调参过程、代码实现细节；
  2. 灵敏度分析、鲁棒性分析走个过场，没有针对性的扰动实验设计；
  3. 问题改进建议部分千篇一律，没有结合题目实际场景提出落地的优化方向。

# 快速淡化 AI 味的修改技巧

- 删掉冗余套话，把长句拆成短句，加入口语化的工程表述；
- 每一段都绑定本题的具体数据、参数、模型变量，去掉可套用的通用内容；
- 补充试错过程，比如 “最初尝试 XX 算法，因 XX 问题收敛慢，最终改用 pareto 优化”，增加真实感；
- 调整段落长短，打乱过于规整的行文节奏。

# 去 AI 味核心规则（整合自 humanizer-zh）

## 五条核心原则

1. **删除填充短语** — 去除开场白和强调性拐杖词（”值得注意的是””综上所述””具有重要意义”）
2. **打破公式结构** — 避免二元对比、戏剧性分段、修辞性设置、否定式排比（”不仅……而且……”）
3. **变化节奏** — 混合句子长度，两项优于三项，段落结尾要多样化
4. **信任读者** — 直接陈述事实，跳过软化、辩解和过度解释
5. **删除金句** — 如果听起来像可引用的语句，重写它

## 必须检测并修复的 AI 模式

- **夸大的象征意义**：”作为……的证明””标志着……的里程碑” → 直接说做了什么
- **宣传性语言**：”无缝、直观和强大””革命性的” → 用具体功能和数据替代
- **三段式法则**：三个形容词/名词并列 → 减到两个或用具体描述替换
- **模糊归因**：”业内专家认为””研究表明” → 给出具体来源或删除
- **AI 高频词汇**：此外、综上所述、关键作用、不断演变的格局、奠定基础 → 直接删除或用口语替代
- **否定式排比**：”这不仅仅是X，而是Y” → 直接说是什么
- **过多连接性短语**：此外、与此同时、值得注意的是、需要指出的是 → 大部分可以直接删除

## 注入真实感的方法

- **有观点**：不要只报告事实，对结果做出反应，写出选择背后的真实考量
- **承认复杂性**：真实的分析会有犹豫和权衡，”效果不错但计算开销偏大”比”效果显著”更真实
- **允许一些混乱**：完美的结构感觉像算法，适当的试错叙述和半成型想法是人性的体现
- **对数据要具体**：不是”误差满足预设精度”，而是”MSE 降到 0.003，比基线低 40%”

## 改写示例

改写前（AI 味）：
> 新的软件更新作为公司致力于创新的证明。此外，它提供了无缝、直观和强大的用户体验——确保用户能够高效地完成目标。

改写后（人性化）：
> 软件更新添加了批处理、键盘快捷键和离线模式。来自测试用户的早期反馈是积极的，大多数报告任务完成速度更快。

# 整体硬性排版要求

- 摘要：限制一页篇幅。每个问题摘要统一以针对问题 X：开头；每个问题段落内部使用首先、其次、然后、最后串联行文（首先说明题目特点和数据特征，其次说明方法选择依据，然后描述具体做法，最后给出关键结果数值），注意这些过渡词用在同一个问题段落内部，不是用来串联不同问题之间。摘要末尾罗列模型创新点，创新点带有序号，序号不必强制换行。在说用了什么方法前，可以适当说明题目特点，给方法选择提供依据。摘要中不放数学公式，数值用行内文字表述。
- 每个问题段落内需要写明加粗展示模型名称，同时展示实验结果，关键数值加粗标注。摘要总字数控制 700–800 字，正文常规字体。
- 正文章节命名：第一大点命名为问题重述，第二大点命名为问题分析；各个小题独立分段分析。
- 图文严格对应：正文提及、描述完成的图像紧跟文字后方，图片环境添加 [H] 固定位置。所有正文段落首行缩进两字符。
- 全文禁止使用破折号，不准使用冒号；流程图支持树状结构等多种样式排版；全文总页数控制在 30 页以内。
- 统一文中人称。
