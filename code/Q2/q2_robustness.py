"""
Q2 鲁棒性分析:
1. BMI分组断点敏感性 (±1, ±2)
2. 风险函数权重敏感性 (w_delay / w_inaccuracy 组合)
3. 4%阈值敏感性 (3.5%, 4.0%, 4.5%)
Seed: 42
"""

import os, sys, json, warnings, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, optimize
from patsy import dmatrix, build_design_matrices

SEED = 42
np.random.seed(SEED)
warnings.filterwarnings('ignore')

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_ROOT, 'workspace', 'data_clean', 'male_cleaned.csv')
OUT_DIR = os.path.join(PROJECT_ROOT, 'robustness', 'Q2')
FIG_DIR = os.path.join(OUT_DIR, 'figures')
for d in [OUT_DIR, FIG_DIR]:
    os.makedirs(d, exist_ok=True)

for fn in ['SimHei', 'Microsoft YaHei', 'STHeiti']:
    try:
        rcParams['font.sans-serif'] = [fn] + rcParams['font.sans-serif']
        break
    except Exception:
        pass
rcParams['axes.unicode_minus'] = False

# ── Load & refit GAMM (same as Q2 main) ─────────────────────────────
print("[Q2 Robustness] Loading data & refitting GAMM...")
df = pd.read_csv(DATA_PATH, encoding='utf-8-sig')
cols = list(df.columns)
df = df.rename(columns={
    cols[1]: 'patient', cols[9]: 'gest_week_raw', cols[10]: 'bmi',
    cols[21]: 'y_conc', cols[31]: 'gest_week', cols[2]: 'age', cols[6]: 'ivf',
})
df = df[['patient','gest_week','bmi','y_conc','age','ivf']].dropna(
    subset=['gest_week','bmi','y_conc','patient'])
df['patient'] = df['patient'].astype(str)
print(f"   Rows={len(df)}, Patients={df['patient'].nunique()}")

# Refit GAMM
spline_rhs = "bs(gest_week, df=4, include_intercept=False) + bs(bmi, df=3, include_intercept=False)"
X_design = dmatrix(spline_rhs, data=df, return_type='dataframe')
design_info = X_design.design_info
safe_cols = []
for c in X_design.columns:
    if c == 'Intercept': safe_cols.append('sp_intercept')
    elif 'gest_week' in c: safe_cols.append(f'gw_s{len([x for x in safe_cols if x.startswith("gw_s")])}')
    else: safe_cols.append(f'bmi_s{len([x for x in safe_cols if x.startswith("bmi_s")])}')
X_design.columns = safe_cols
if 'sp_intercept' in X_design.columns:
    X_design = X_design.drop(columns=['sp_intercept'])
    safe_cols = [c for c in safe_cols if c != 'sp_intercept']
df_gamm = pd.concat([df.reset_index(drop=True), X_design.reset_index(drop=True)], axis=1)
fe_formula = "y_conc ~ " + " + ".join(safe_cols)
mdf_gamm = None
for opt in ['lbfgs', 'powell', 'cg']:
    try:
        md_g = smf.mixedlm(fe_formula, data=df_gamm, groups=df_gamm["patient"],
                           re_formula="~gest_week")
        fit_g = md_g.fit(reml=True, method=opt)
        if fit_g.converged:
            mdf_gamm = fit_g; break
    except Exception:
        continue
if mdf_gamm is None:
    md_g = smf.mixedlm(fe_formula, data=df_gamm, groups=df_gamm["patient"])
    mdf_gamm = md_g.fit(reml=True, method='lbfgs')
fe = mdf_gamm.fe_params
resid_std = float(mdf_gamm.resid.std())
GW_MIN = float(df['gest_week'].min())
GW_MAX = float(df['gest_week'].max())
print(f"   GAMM converged={mdf_gamm.converged}, resid_std={resid_std:.6f}")

def rename_pred(X):
    orig = list(X.columns); new = []; gi = bi = 0
    for c in orig:
        if c == 'Intercept': new.append('_drop')
        elif 'gest_week' in c: new.append(f'gw_s{gi}'); gi += 1
        else: new.append(f'bmi_s{bi}'); bi += 1
    X.columns = new
    if '_drop' in X.columns: X = X.drop(columns=['_drop'])
    X.insert(0, 'Intercept', 1.0)
    return X

def predict_yconc(gw_arr, bmi_arr):
    gw_c = np.clip(gw_arr, GW_MIN, GW_MAX)
    bmi_c = np.clip(bmi_arr, df['bmi'].min(), df['bmi'].max())
    pdf = pd.DataFrame({'gest_week': gw_c, 'bmi': bmi_c})
    Xp = build_design_matrices([design_info], pdf, return_type='dataframe')[0]
    Xp = rename_pred(Xp.copy())
    return (Xp[fe.index] @ fe).values

# Per-patient empirical threshold time
patient_rows = []
for pid, grp in df.groupby('patient'):
    grp_s = grp.sort_values('gest_week')
    bmi_val = grp_s['bmi'].iloc[0]
    reached = grp_s[grp_s['y_conc'] >= 0.04]
    t_emp = reached['gest_week'].min() if len(reached) > 0 else np.nan
    patient_rows.append({'patient': pid, 'bmi': bmi_val, 't_empirical': t_emp})
df_pts = pd.DataFrame(patient_rows)
df_valid = df_pts.dropna(subset=['t_empirical']).copy()

# Risk functions
def risk_delay(t):
    return 1.0 / (1.0 + np.exp(-0.5 * (t - 18.0)))

def risk_inaccuracy(t, bmi_val, threshold=0.04, noise_std=None):
    if noise_std is None: noise_std = resid_std
    pred = predict_yconc(np.array([t]), np.array([bmi_val]))[0]
    z = (threshold - pred) / noise_std
    return float(stats.norm.cdf(z))

def total_risk(t, bmi_val, w_delay=0.4, w_inacc=0.6, threshold=0.04):
    return w_delay * risk_delay(t) + w_inacc * risk_inaccuracy(t, bmi_val, threshold)

def optimal_nipt_time(bmi_val, gw_range=(10.0, 25.0), **kwargs):
    lo = max(gw_range[0], GW_MIN)
    hi = min(gw_range[1], GW_MAX)
    res = optimize.minimize_scalar(lambda t: total_risk(t, bmi_val, **kwargs),
                                   bounds=(lo, hi), method='bounded')
    return float(res.x), float(res.fun)

# Reference breakpoints and groups
REF_BREAKS = [30.0, 36.0]
bmi_lo_data = df['bmi'].min()
bmi_hi_data = df['bmi'].max()

def compute_group_results(breaks, threshold=0.04, w_delay=0.4, w_inacc=0.6):
    edges = [bmi_lo_data] + list(breaks) + [bmi_hi_data + 0.01]
    results = []
    for i in range(len(edges)-1):
        lo, hi = edges[i], edges[i+1]
        mask = (df_valid['bmi'] >= lo) & (df_valid['bmi'] < hi)
        sub = df_valid[mask]
        if len(sub) == 0: continue
        bmi_mid = float(sub['bmi'].mean())
        t_opt, r_opt = optimal_nipt_time(bmi_mid, threshold=threshold,
                                          w_delay=w_delay, w_inacc=w_inacc)
        results.append({
            'bmi_lo': round(lo, 1), 'bmi_hi': round(hi, 1),
            'n': int(mask.sum()), 'bmi_mean': round(bmi_mid, 1),
            'optimal_week': round(t_opt, 2), 'min_risk': round(r_opt, 4),
        })
    return results

ref_results = compute_group_results(REF_BREAKS)
print("   Reference results (breaks=[30,36]):")
for r in ref_results:
    print(f"     [{r['bmi_lo']},{r['bmi_hi']}): n={r['n']}, opt={r['optimal_week']}w")

# ══════════════════════════════════════════════════════════════════════
# 1. BMI分组断点敏感性 (±1, ±2)
# ══════════════════════════════════════════════════════════════════════
print("\n[1/3] BMI breakpoint sensitivity...")
bp_sensitivity = []
for delta1 in [-2, -1, 0, 1, 2]:
    for delta2 in [-2, -1, 0, 1, 2]:
        bp1 = 30.0 + delta1
        bp2 = 36.0 + delta2
        if bp2 <= bp1 + 2:
            continue  # ensure minimum gap
        breaks = [bp1, bp2]
        res = compute_group_results(breaks)
        opt_weeks = [r['optimal_week'] for r in res]
        bp_sensitivity.append({
            'bp1': bp1, 'bp2': bp2, 'delta1': delta1, 'delta2': delta2,
            'n_groups': len(res),
            'group_sizes': [r['n'] for r in res],
            'optimal_weeks': opt_weeks,
            'min_risks': [r['min_risk'] for r in res],
        })

# Compare to reference
ref_weeks = [r['optimal_week'] for r in ref_results]
for bp in bp_sensitivity:
    if len(bp['optimal_weeks']) == len(ref_weeks):
        max_week_diff = max(abs(a - b) for a, b in zip(bp['optimal_weeks'], ref_weeks))
    else:
        max_week_diff = float('nan')
    bp['max_week_diff_vs_ref'] = round(max_week_diff, 2)

print(f"   Tested {len(bp_sensitivity)} breakpoint combinations")
diffs = [b['max_week_diff_vs_ref'] for b in bp_sensitivity if not np.isnan(b['max_week_diff_vs_ref'])]
print(f"   Max week diff range: [{min(diffs):.2f}, {max(diffs):.2f}]")
bp_stable = all(d < 1.5 for d in diffs)

# ══════════════════════════════════════════════════════════════════════
# 2. 风险函数权重敏感性
# ══════════════════════════════════════════════════════════════════════
print("\n[2/3] Risk function weight sensitivity...")
weight_combos = [(0.2, 0.8), (0.3, 0.7), (0.4, 0.6), (0.5, 0.5), (0.6, 0.4), (0.7, 0.3)]
weight_results = []
for wd, wi in weight_combos:
    res = compute_group_results(REF_BREAKS, w_delay=wd, w_inacc=wi)
    opt_weeks = [r['optimal_week'] for r in res]
    weight_results.append({
        'w_delay': wd, 'w_inaccuracy': wi,
        'optimal_weeks': opt_weeks,
        'min_risks': [r['min_risk'] for r in res],
    })
    print(f"   w_d={wd}, w_i={wi}: opt_weeks={opt_weeks}")

# Check stability
all_opt_g1 = [w['optimal_weeks'][0] for w in weight_results if len(w['optimal_weeks']) >= 1]
all_opt_g2 = [w['optimal_weeks'][1] for w in weight_results if len(w['optimal_weeks']) >= 2]
weight_range_g1 = round(max(all_opt_g1) - min(all_opt_g1), 2) if all_opt_g1 else float('nan')
weight_range_g2 = round(max(all_opt_g2) - min(all_opt_g2), 2) if all_opt_g2 else float('nan')
print(f"   G1 optimal week range across weights: {weight_range_g1}")
print(f"   G2 optimal week range across weights: {weight_range_g2}")

# ══════════════════════════════════════════════════════════════════════
# 3. 阈值敏感性 (3.5%, 4.0%, 4.5%)
# ══════════════════════════════════════════════════════════════════════
print("\n[3/3] Threshold sensitivity (3.5%, 4.0%, 4.5%)...")
threshold_results = []
for thr_pct in [3.0, 3.5, 4.0, 4.5, 5.0]:
    thr = thr_pct / 100.0
    res = compute_group_results(REF_BREAKS, threshold=thr)
    opt_weeks = [r['optimal_week'] for r in res]
    threshold_results.append({
        'threshold_pct': thr_pct,
        'optimal_weeks': opt_weeks,
        'min_risks': [r['min_risk'] for r in res],
    })
    print(f"   threshold={thr_pct}%: opt_weeks={opt_weeks}")

thr_weeks_g1 = [t['optimal_weeks'][0] for t in threshold_results if len(t['optimal_weeks']) >= 1]
thr_range_g1 = round(max(thr_weeks_g1) - min(thr_weeks_g1), 2) if thr_weeks_g1 else float('nan')
print(f"   G1 optimal week range across thresholds: {thr_range_g1}")

# ══════════════════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════════════════
print("\nGenerating figures...")

# Fig 1: Breakpoint sensitivity heatmap (G2 optimal week)
bp1_vals = sorted(set(b['bp1'] for b in bp_sensitivity))
bp2_vals = sorted(set(b['bp2'] for b in bp_sensitivity))
heat = np.full((len(bp1_vals), len(bp2_vals)), np.nan)
for b in bp_sensitivity:
    if len(b['optimal_weeks']) >= 2:
        i = bp1_vals.index(b['bp1'])
        j = bp2_vals.index(b['bp2'])
        heat[i, j] = b['optimal_weeks'][1]  # G2 optimal week

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(heat, cmap='YlOrRd', aspect='auto')
ax.set_xticks(range(len(bp2_vals)))
ax.set_xticklabels([f'{v:.0f}' for v in bp2_vals])
ax.set_yticks(range(len(bp1_vals)))
ax.set_yticklabels([f'{v:.0f}' for v in bp1_vals])
ax.set_xlabel('Breakpoint 2 (upper)')
ax.set_ylabel('Breakpoint 1 (lower)')
for i in range(len(bp1_vals)):
    for j in range(len(bp2_vals)):
        if not np.isnan(heat[i, j]):
            ax.text(j, i, f'{heat[i,j]:.1f}', ha='center', va='center', fontsize=8)
plt.colorbar(im, ax=ax, label='G2 Optimal NIPT Week')
ax.set_title('Q2: Breakpoint Sensitivity (G2 Optimal Week)')
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q2_breakpoint_sensitivity.png'), dpi=150)
plt.close(fig)

# Fig 2: Weight sensitivity
fig, ax = plt.subplots(figsize=(8, 5))
wd_vals = [w['w_delay'] for w in weight_results]
for gi in range(min(len(w['optimal_weeks']) for w in weight_results)):
    weeks = [w['optimal_weeks'][gi] for w in weight_results]
    ax.plot(wd_vals, weeks, 'o-', label=f'Group {gi+1}', lw=2)
ax.set_xlabel('w_delay')
ax.set_ylabel('Optimal NIPT Week')
ax.set_title('Q2: Risk Weight Sensitivity')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q2_weight_sensitivity.png'), dpi=150)
plt.close(fig)

# Fig 3: Threshold sensitivity
fig, ax = plt.subplots(figsize=(8, 5))
thr_vals = [t['threshold_pct'] for t in threshold_results]
for gi in range(min(len(t['optimal_weeks']) for t in threshold_results)):
    weeks = [t['optimal_weeks'][gi] for t in threshold_results]
    ax.plot(thr_vals, weeks, 's-', label=f'Group {gi+1}', lw=2)
ax.set_xlabel('Y-conc Threshold (%)')
ax.set_ylabel('Optimal NIPT Week')
ax.set_title('Q2: Detection Threshold Sensitivity')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q2_threshold_sensitivity.png'), dpi=150)
plt.close(fig)

# ══════════════════════════════════════════════════════════════════════
# Save summary
# ══════════════════════════════════════════════════════════════════════
summary = {
    "question": "Q2",
    "robustness_checks": [
        {
            "check": "bmi_breakpoint_sensitivity",
            "claim": "Optimal NIPT timing is stable under ±2 BMI breakpoint shifts",
            "perturbation": "bp1 in [28,32], bp2 in [34,38], step=1",
            "n_tested": len(bp_sensitivity),
            "max_week_diff_range": [min(diffs), max(diffs)],
            "results": bp_sensitivity,
            "status": "PASS" if bp_stable else "CONDITIONAL",
            "limitation": "Only tested integer shifts around data-driven breakpoints"
        },
        {
            "check": "risk_weight_sensitivity",
            "claim": "Optimal timing is moderately sensitive to risk weight allocation",
            "perturbation": "w_delay in [0.2, 0.7], w_inaccuracy = 1 - w_delay",
            "results": weight_results,
            "g1_week_range": weight_range_g1,
            "g2_week_range": weight_range_g2,
            "status": "PASS" if max(weight_range_g1, weight_range_g2) < 3.0 else "CONDITIONAL",
            "limitation": "Linear weight combinations only; nonlinear risk aggregation not tested"
        },
        {
            "check": "threshold_sensitivity",
            "claim": "Optimal timing shifts predictably with Y-conc threshold",
            "perturbation": "threshold in {3.0%, 3.5%, 4.0%, 4.5%, 5.0%}",
            "results": threshold_results,
            "g1_week_range": thr_range_g1,
            "status": "PASS" if thr_range_g1 < 3.0 else "CONDITIONAL",
            "limitation": "Threshold perturbation does not account for measurement error in threshold itself"
        },
    ],
    "seed": SEED,
    "figures": ["q2_breakpoint_sensitivity.png", "q2_weight_sensitivity.png", "q2_threshold_sensitivity.png"],
}

with open(os.path.join(OUT_DIR, 'q2_robustness_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print(f"\nQ2 robustness done. Outputs in {OUT_DIR}")
for chk in summary['robustness_checks']:
    print(f"   {chk['check']}: {chk['status']}")
