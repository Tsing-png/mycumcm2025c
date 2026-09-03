"""
Q2: BMI分组 + 最佳NIPT时点优化 + 检测误差敏感性分析
M1 (main): 混合策略——GAMM连续预测 + 经验达标时间回归 + 风险函数优化 + 数据驱动BMI分组
M2 (baseline): 临床经验BMI分组 + 组内统计
依赖: Q1 GAMM 模型结构
Seed: 42
"""

import os, sys, json, time, warnings
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
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ── Paths ────────────────────────────────────────────────────────────
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_ROOT, 'workspace', 'data_clean', 'male_cleaned.csv')
OUT_DIR = os.path.join(PROJECT_ROOT, 'results', 'Q2', 'experiments', 'round1')
FIG_DIR = os.path.join(OUT_DIR, 'figures')
TBL_DIR = os.path.join(OUT_DIR, 'tables')
MET_DIR = os.path.join(OUT_DIR, 'metrics')
for d in [FIG_DIR, TBL_DIR, MET_DIR]:
    os.makedirs(d, exist_ok=True)

for fn in ['SimHei', 'Microsoft YaHei', 'STHeiti']:
    try:
        rcParams['font.sans-serif'] = [fn] + rcParams['font.sans-serif']
        break
    except Exception:
        pass
rcParams['axes.unicode_minus'] = False

# ── 1. Load & prep ──────────────────────────────────────────────────
print("[1/8] Loading data...")
t_start = time.time()
df = pd.read_csv(DATA_PATH, encoding='utf-8-sig')
cols = list(df.columns)
df = df.rename(columns={
    cols[1]: 'patient', cols[9]: 'gest_week_raw', cols[10]: 'bmi',
    cols[21]: 'y_conc', cols[31]: 'gest_week', cols[2]: 'age', cols[6]: 'ivf',
})
df = df[['patient','gest_week','bmi','y_conc','age','ivf']].dropna(
    subset=['gest_week','bmi','y_conc','patient'])
df['patient'] = df['patient'].astype(str)
N_OBS = len(df)
N_PAT = df['patient'].nunique()
print(f"   Rows={N_OBS}, Patients={N_PAT}")
print(f"   BMI: [{df['bmi'].min():.1f}, {df['bmi'].max():.1f}]")
print(f"   GW:  [{df['gest_week'].min():.1f}, {df['gest_week'].max():.1f}]")

# ── 2. Refit Q1 GAMM ────────────────────────────────────────────────
print("\n[2/8] Refitting Q1 GAMM for continuous prediction...")
t0 = time.time()

spline_rhs = "bs(gest_week, df=4, include_intercept=False) + bs(bmi, df=3, include_intercept=False)"
X_design = dmatrix(spline_rhs, data=df, return_type='dataframe')
design_info = X_design.design_info

safe_cols = []
for i, c in enumerate(X_design.columns):
    if c == 'Intercept':
        safe_cols.append('sp_intercept')
    elif 'gest_week' in c:
        safe_cols.append(f'gw_s{len([x for x in safe_cols if x.startswith("gw_s")])}')
    else:
        safe_cols.append(f'bmi_s{len([x for x in safe_cols if x.startswith("bmi_s")])}')
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
            mdf_gamm = fit_g
            print(f"   GAMM converged with '{opt}' in {time.time()-t0:.2f}s")
            break
    except Exception:
        continue
if mdf_gamm is None:
    md_g = smf.mixedlm(fe_formula, data=df_gamm, groups=df_gamm["patient"])
    mdf_gamm = md_g.fit(reml=True, method='lbfgs')
    print(f"   GAMM (intercept-only RE) in {time.time()-t0:.2f}s")

fe = mdf_gamm.fe_params
resid_std = float(mdf_gamm.resid.std())
print(f"   Residual std = {resid_std:.6f}")

GW_MIN = float(df['gest_week'].min())
GW_MAX = float(df['gest_week'].max())

def rename_pred(X):
    orig = list(X.columns); new = []; gi = bi = 0
    for c in orig:
        if c == 'Intercept':      new.append('_drop')
        elif 'gest_week' in c:    new.append(f'gw_s{gi}'); gi += 1
        else:                     new.append(f'bmi_s{bi}'); bi += 1
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

# ── 3. Per-patient empirical threshold time ──────────────────────────
print("\n[3/8] Computing per-patient empirical threshold times...")
patient_rows = []
for pid, grp in df.groupby('patient'):
    grp_s = grp.sort_values('gest_week')
    bmi_val = grp_s['bmi'].iloc[0]
    reached = grp_s[grp_s['y_conc'] >= 0.04]
    t_emp = reached['gest_week'].min() if len(reached) > 0 else np.nan
    patient_rows.append({'patient': pid, 'bmi': bmi_val, 't_empirical': t_emp})
df_pts = pd.DataFrame(patient_rows)
n_reached = df_pts['t_empirical'].notna().sum()
print(f"   {n_reached}/{len(df_pts)} patients reached 4% threshold")

# ── 4. Empirical regression: T_emp ~ f(BMI) ─────────────────────────
print("\n[4/8] Fitting empirical threshold-time model T*(BMI)...")
df_valid = df_pts.dropna(subset=['t_empirical']).copy()

# Fit a robust regression: T_emp ~ BMI + BMI^2 (allows nonlinearity)
df_valid['bmi2'] = df_valid['bmi'] ** 2
X_reg = sm.add_constant(df_valid[['bmi', 'bmi2']])
rlm_fit = sm.RLM(df_valid['t_empirical'], X_reg, M=sm.robust.norms.HuberT()).fit()
print(f"   RLM coefficients: const={rlm_fit.params['const']:.4f}, bmi={rlm_fit.params['bmi']:.4f}, bmi2={rlm_fit.params['bmi2']:.6f}")

bmi_grid = np.linspace(df['bmi'].min(), df['bmi'].max(), 300)
t_star_empirical = rlm_fit.params['const'] + rlm_fit.params['bmi'] * bmi_grid + rlm_fit.params['bmi2'] * bmi_grid**2
print(f"   T* range (empirical model): [{t_star_empirical.min():.2f}, {t_star_empirical.max():.2f}] weeks")

# ── 5. Data-driven BMI grouping ─────────────────────────────────────
print("\n[5/8] Data-driven BMI grouping via optimal risk-based segmentation...")

# Strategy: use the empirical T*(BMI) curve to find natural breakpoints.
# Test piecewise-constant fits (each segment gets one representative T*)
# to find groupings that minimize within-group T* variance.
from itertools import combinations

def segmented_variance(bmi_vals, t_vals, breaks):
    """Total within-segment variance of t_vals given BMI breakpoints."""
    segs = np.digitize(bmi_vals, breaks)
    total_var = 0.0
    for s in range(len(breaks) + 1):
        mask = segs == s
        if mask.sum() > 1:
            total_var += np.var(t_vals[mask]) * mask.sum()
        # single-element segments contribute 0
    return total_var

# Use per-patient empirical data for grouping
bmi_pts = df_valid['bmi'].values
t_pts = df_valid['t_empirical'].values

# Candidate breakpoints: every 1.0 BMI unit
candidates = np.arange(24, 45, 1.0)
bmi_p5, bmi_p95 = np.percentile(bmi_pts, [5, 95])
candidates = candidates[(candidates > bmi_p5) & (candidates < bmi_p95)]

best_bic = np.inf
best_breaks = None
n_pts = len(bmi_pts)
MIN_WIDTH = 3.0   # minimum BMI range per group
MIN_COUNT = 5     # minimum patients per group

for n_bp in [2, 3, 4]:
    for combo in combinations(candidates, n_bp):
        breaks = sorted(combo)
        # Enforce minimum width between consecutive breakpoints and boundaries
        all_edges = [bmi_pts.min()] + list(breaks) + [bmi_pts.max()]
        widths_ok = all(all_edges[j+1] - all_edges[j] >= MIN_WIDTH for j in range(len(all_edges)-1))
        if not widths_ok:
            continue
        segs = np.digitize(bmi_pts, breaks)
        counts = [np.sum(segs == s) for s in range(n_bp + 1)]
        if min(counts) < MIN_COUNT:
            continue
        sv = segmented_variance(bmi_pts, t_pts, breaks)
        k = 2 * (n_bp + 1)
        bic = n_pts * np.log(sv / n_pts + 1e-12) + k * np.log(n_pts)
        if bic < best_bic:
            best_bic = bic
            best_breaks = list(breaks)

print(f"   Best breakpoints (BIC={best_bic:.2f}): {best_breaks}")

# Also try the clinical reference for comparison
clinical_edges_inner = [28, 32, 36, 40]
sv_clinical = segmented_variance(bmi_pts, t_pts, clinical_edges_inner)
bic_clinical = n_pts * np.log(sv_clinical / n_pts + 1e-12) + 2*5 * np.log(n_pts)
print(f"   Clinical grouping BIC: {bic_clinical:.2f}")

# Build final groups
bmi_lo_data = df['bmi'].min()
bmi_hi_data = df['bmi'].max()
dd_edges = [bmi_lo_data] + best_breaks + [bmi_hi_data + 0.01]
n_groups = len(dd_edges) - 1
dd_labels = []
for i in range(n_groups):
    lo_r = round(dd_edges[i], 1)
    hi_r = round(dd_edges[i+1], 1) if i < n_groups - 1 else round(dd_edges[i+1] - 0.01, 1)
    dd_labels.append(f"G{i+1}: [{lo_r}, {hi_r})")
print(f"   Data-driven groups ({n_groups}): {dd_labels}")

# Clinical groups
clinical_edges = [20, 28, 32, 36, 40, 50]
clinical_labels = ['[20,28)', '[28,32)', '[32,36)', '[36,40)', '>=40']

# ── 6. Risk function & optimal NIPT timing ──────────────────────────
print("\n[6/8] Building risk function & optimizing NIPT timing...")

def risk_delay(t):
    """Delay risk: sigmoid increasing after 12 weeks."""
    return 1.0 / (1.0 + np.exp(-0.5 * (t - 18.0)))

def risk_inaccuracy(t, bmi_val, threshold=0.04, noise_std=None):
    """P(Y-conc < threshold) at week t for given BMI, from GAMM + residual."""
    if noise_std is None:
        noise_std = resid_std
    pred = predict_yconc(np.array([t]), np.array([bmi_val]))[0]
    z = (threshold - pred) / noise_std
    return float(stats.norm.cdf(z))

def total_risk(t, bmi_val, w_delay=0.4, w_inacc=0.6, threshold=0.04, noise_std=None):
    return w_delay * risk_delay(t) + w_inacc * risk_inaccuracy(t, bmi_val, threshold, noise_std)

def optimal_nipt_time(bmi_val, gw_range=(10.0, 25.0), **kwargs):
    lo = max(gw_range[0], GW_MIN)
    hi = min(gw_range[1], GW_MAX)
    res = optimize.minimize_scalar(lambda t: total_risk(t, bmi_val, **kwargs),
                                   bounds=(lo, hi), method='bounded')
    return float(res.x), float(res.fun)

# Compute optimal time on BMI grid
opt_times = np.array([optimal_nipt_time(b)[0] for b in bmi_grid])
opt_risks = np.array([optimal_nipt_time(b)[1] for b in bmi_grid])

# ── M1 group-level results ──────────────────────────────────────────
m1_results = []
for i in range(n_groups):
    lo, hi = dd_edges[i], dd_edges[i+1]
    pt_mask = (df_valid['bmi'] >= lo) & (df_valid['bmi'] < hi)
    sub = df_valid[pt_mask]
    bmi_mid = sub['bmi'].mean() if len(sub) > 0 else (lo + hi) / 2
    t_opt, r_opt = optimal_nipt_time(bmi_mid)
    t_emp_med = sub['t_empirical'].median() if len(sub) > 0 else np.nan
    t_emp_mean = sub['t_empirical'].mean() if len(sub) > 0 else np.nan
    # All patients (incl. those not reaching threshold)
    all_mask = (df_pts['bmi'] >= lo) & (df_pts['bmi'] < hi)
    n_all = int(all_mask.sum())
    n_reached = int(pt_mask.sum())
    m1_results.append({
        'group': dd_labels[i],
        'bmi_lo': round(lo, 1), 'bmi_hi': round(hi, 1),
        'n_patients': n_all, 'n_reached_threshold': n_reached,
        'bmi_mean': round(float(bmi_mid), 1),
        'empirical_threshold_median_week': round(float(t_emp_med), 2) if pd.notna(t_emp_med) else None,
        'empirical_threshold_mean_week': round(float(t_emp_mean), 2) if pd.notna(t_emp_mean) else None,
        'optimal_nipt_week': round(t_opt, 2),
        'min_total_risk': round(r_opt, 4),
        'risk_delay_at_opt': round(float(risk_delay(t_opt)), 4),
        'risk_inaccuracy_at_opt': round(float(risk_inaccuracy(t_opt, bmi_mid)), 4),
    })

df_m1 = pd.DataFrame(m1_results)
print("\n   M1 Results (data-driven groups):")
print(df_m1[['group','n_patients','bmi_mean','empirical_threshold_median_week',
             'optimal_nipt_week','min_total_risk']].to_string(index=False))

# ── M2 baseline: clinical groups ────────────────────────────────────
print("\n   M2 Baseline (clinical groups):")
m2_results = []
for i in range(len(clinical_edges)-1):
    lo, hi = clinical_edges[i], clinical_edges[i+1]
    pt_mask_all = (df_pts['bmi'] >= lo) & (df_pts['bmi'] < hi)
    pt_mask = (df_valid['bmi'] >= lo) & (df_valid['bmi'] < hi)
    sub = df_valid[pt_mask]
    if pt_mask_all.sum() == 0:
        continue
    bmi_mid = sub['bmi'].mean() if len(sub) > 0 else (lo + hi) / 2
    t_opt, r_opt = optimal_nipt_time(bmi_mid)
    t_emp_med = sub['t_empirical'].median() if len(sub) > 0 else np.nan
    m2_results.append({
        'group': clinical_labels[i],
        'bmi_lo': lo, 'bmi_hi': hi,
        'n_patients': int(pt_mask_all.sum()),
        'n_reached_threshold': int(pt_mask.sum()),
        'bmi_mean': round(float(bmi_mid), 1),
        'empirical_threshold_median_week': round(float(t_emp_med), 2) if pd.notna(t_emp_med) else None,
        'optimal_nipt_week': round(t_opt, 2),
        'min_total_risk': round(r_opt, 4),
    })
df_m2 = pd.DataFrame(m2_results)
print(df_m2[['group','n_patients','bmi_mean','empirical_threshold_median_week',
             'optimal_nipt_week','min_total_risk']].to_string(index=False))

# ── 7. Error sensitivity analysis ────────────────────────────────────
print("\n[7/8] Detection error sensitivity analysis...")

error_levels = [0.0, 0.005, 0.01, 0.02]
bmi_test_points = []
for row in m1_results:
    bmi_test_points.append(row['bmi_mean'])

sensitivity_results = []
for bmi_val in bmi_test_points:
    for err in error_levels:
        eff_thresh = 0.04 + err
        t_opt, r_opt = optimal_nipt_time(bmi_val, threshold=eff_thresh)
        sensitivity_results.append({
            'bmi': bmi_val, 'error_type': 'threshold_shift',
            'error_pct': round(err * 100, 1),
            'effective_threshold_pct': round(eff_thresh * 100, 2),
            'optimal_nipt_week': round(t_opt, 2),
            'min_risk': round(r_opt, 4),
        })
    for noise_mult in [1.0, 1.5, 2.0, 3.0]:
        t_opt_n, r_opt_n = optimal_nipt_time(bmi_val, noise_std=resid_std * noise_mult)
        sensitivity_results.append({
            'bmi': bmi_val, 'error_type': 'noise_inflation',
            'noise_multiplier': noise_mult,
            'effective_noise_std': round(resid_std * noise_mult, 6),
            'optimal_nipt_week': round(t_opt_n, 2),
            'min_risk': round(r_opt_n, 4),
        })

df_sens = pd.DataFrame(sensitivity_results)
print("   Threshold shift effect:")
sub_ts = df_sens[df_sens['error_type'] == 'threshold_shift']
print(sub_ts[['bmi','error_pct','optimal_nipt_week','min_risk']].to_string(index=False))
print("\n   Noise inflation effect:")
sub_ni = df_sens[df_sens['error_type'] == 'noise_inflation']
print(sub_ni[['bmi','noise_multiplier','optimal_nipt_week','min_risk']].to_string(index=False))

# ── 8. Figures & outputs ─────────────────────────────────────────────
print("\n[8/8] Generating figures and saving outputs...")

# --- Fig 1: T*(BMI) empirical + fitted curve + group boundaries ---
fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(df_valid['bmi'], df_valid['t_empirical'], alpha=0.4, s=20, c='gray',
           label='Per-patient empirical T*', zorder=2)
ax.plot(bmi_grid, t_star_empirical, 'b-', lw=2.5, label='Robust regression T*(BMI)', zorder=3)
for bp in best_breaks:
    ax.axvline(bp, color='red', ls='--', alpha=0.8, lw=1.2)
# annotate groups
for i in range(n_groups):
    mid_x = (dd_edges[i] + dd_edges[i+1]) / 2
    ax.text(mid_x, ax.get_ylim()[1] * 0.95, dd_labels[i], ha='center', fontsize=8,
            color='red', fontweight='bold')
ax.axhline(12, color='green', ls=':', alpha=0.5, label='12w boundary')
ax.set_xlabel('BMI', fontsize=12)
ax.set_ylabel('Gestational Week at Y-conc >= 4%', fontsize=12)
ax.set_title('Threshold-Reaching Time vs BMI', fontsize=13)
ax.legend(fontsize=9, loc='upper left')
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'M1_threshold_time_vs_bmi.png'), dpi=150)
plt.close(fig)

# --- Fig 2: Risk function components ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
gw_range = np.linspace(GW_MIN, GW_MAX, 200)

# 2a: Risk components for a mid-group BMI
bmi_show = m1_results[len(m1_results)//2]['bmi_mean']
r_del = np.array([risk_delay(t) for t in gw_range])
r_ina = np.array([risk_inaccuracy(t, bmi_show) for t in gw_range])
r_tot = np.array([total_risk(t, bmi_show) for t in gw_range])
axes[0].plot(gw_range, r_del, 'r--', label='R_delay')
axes[0].plot(gw_range, r_ina, 'b--', label='R_inaccuracy')
axes[0].plot(gw_range, r_tot, 'k-', lw=2, label='Total Risk')
t_opt_show, _ = optimal_nipt_time(bmi_show)
axes[0].axvline(t_opt_show, color='green', ls=':', label=f'Optimal={t_opt_show:.1f}w')
axes[0].set_xlabel('Gestational Week'); axes[0].set_ylabel('Risk')
axes[0].set_title(f'Risk Components (BMI={bmi_show:.1f})'); axes[0].legend(fontsize=8)

# 2b: Total risk for each group midpoint
colors_grp = plt.cm.tab10(np.linspace(0, 0.8, n_groups))
for idx, row in enumerate(m1_results):
    bv = row['bmi_mean']
    r = np.array([total_risk(t, bv) for t in gw_range])
    t_o = row['optimal_nipt_week']
    axes[1].plot(gw_range, r, color=colors_grp[idx],
                 label=f'{row["group"]} (opt={t_o:.1f}w)')
axes[1].set_xlabel('Gestational Week'); axes[1].set_ylabel('Total Risk')
axes[1].set_title('Total Risk by BMI Group'); axes[1].legend(fontsize=7)

# 2c: Optimal NIPT time vs BMI (continuous)
axes[2].plot(bmi_grid, opt_times, 'b-', lw=2)
axes[2].fill_between(bmi_grid, opt_times - 0.5, opt_times + 0.5, alpha=0.15, color='blue')
for bp in best_breaks:
    axes[2].axvline(bp, color='red', ls='--', alpha=0.7)
axes[2].set_xlabel('BMI'); axes[2].set_ylabel('Optimal NIPT Week')
axes[2].set_title('Optimal NIPT Timing vs BMI (M1)')
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'M1_risk_function_analysis.png'), dpi=150)
plt.close(fig)

# --- Fig 3: Error sensitivity ---
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for bv in bmi_test_points:
    sub = sub_ts[sub_ts['bmi'] == bv]
    axes[0].plot(sub['error_pct'], sub['optimal_nipt_week'], 'o-', label=f'BMI={bv:.0f}')
axes[0].set_xlabel('Concentration Error (%)'); axes[0].set_ylabel('Optimal NIPT Week')
axes[0].set_title('Threshold Shift Effect'); axes[0].legend(fontsize=8)

for bv in bmi_test_points:
    sub = sub_ni[sub_ni['bmi'] == bv]
    axes[1].plot(sub['noise_multiplier'], sub['min_risk'], 'o-', label=f'BMI={bv:.0f}')
axes[1].set_xlabel('Noise Multiplier'); axes[1].set_ylabel('Min Achievable Risk')
axes[1].set_title('Measurement Noise Effect'); axes[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'M1_error_sensitivity.png'), dpi=150)
plt.close(fig)

# --- Fig 4: GAMM predicted Y-conc curves by group ---
fig, ax = plt.subplots(figsize=(10, 6))
gw_plot = np.linspace(GW_MIN, GW_MAX, 200)
for idx, row in enumerate(m1_results):
    bv = row['bmi_mean']
    yp = predict_yconc(gw_plot, np.full_like(gw_plot, bv))
    ax.plot(gw_plot, yp * 100, color=colors_grp[idx], lw=2,
            label=f'{row["group"]} (opt={row["optimal_nipt_week"]:.1f}w)')
    ax.axvline(row['optimal_nipt_week'], color=colors_grp[idx], ls=':', alpha=0.5)
ax.axhline(4, color='black', ls='--', lw=1, label='4% threshold')
ax.set_xlabel('Gestational Week', fontsize=12)
ax.set_ylabel('Y Chromosome Concentration (%)', fontsize=12)
ax.set_title('GAMM Predicted Curves by BMI Group (M1)', fontsize=13)
ax.legend(fontsize=8); ax.set_xlim(GW_MIN, GW_MAX)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'M1_conc_curves_by_group.png'), dpi=150)
plt.close(fig)

# --- Fig 5: M1 vs M2 comparison bar chart ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
# M1
x1 = range(len(m1_results))
axes[0].bar(x1, [r['optimal_nipt_week'] for r in m1_results], color='steelblue', alpha=0.8)
axes[0].set_xticks(x1)
axes[0].set_xticklabels([r['group'] for r in m1_results], rotation=30, ha='right', fontsize=8)
axes[0].set_ylabel('Optimal NIPT Week')
axes[0].set_title('M1: Data-Driven Groups')
for xi, r in zip(x1, m1_results):
    axes[0].text(xi, r['optimal_nipt_week'] + 0.2, f"{r['optimal_nipt_week']:.1f}w", ha='center', fontsize=9)
# M2
x2 = range(len(m2_results))
axes[1].bar(x2, [r['optimal_nipt_week'] for r in m2_results], color='darkorange', alpha=0.8)
axes[1].set_xticks(x2)
axes[1].set_xticklabels([r['group'] for r in m2_results], rotation=30, ha='right', fontsize=8)
axes[1].set_ylabel('Optimal NIPT Week')
axes[1].set_title('M2: Clinical Reference Groups')
for xi, r in zip(x2, m2_results):
    axes[1].text(xi, r['optimal_nipt_week'] + 0.2, f"{r['optimal_nipt_week']:.1f}w", ha='center', fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'M1_vs_M2_comparison.png'), dpi=150)
plt.close(fig)

# ── Save tables ──────────────────────────────────────────────────────
df_m1.to_csv(os.path.join(TBL_DIR, 'M1_group_results.csv'), index=False, encoding='utf-8-sig')
df_m2.to_csv(os.path.join(TBL_DIR, 'M2_baseline_results.csv'), index=False, encoding='utf-8-sig')
df_sens.to_csv(os.path.join(TBL_DIR, 'error_sensitivity.csv'), index=False, encoding='utf-8-sig')

# ── Save metrics ─────────────────────────────────────────────────────
metrics = {
    'M1_data_driven_groups': m1_results,
    'M2_clinical_groups': m2_results,
    'data_driven_breakpoints': [round(b, 1) for b in best_breaks],
    'data_driven_bic': round(best_bic, 2),
    'clinical_bic': round(bic_clinical, 2),
    'empirical_regression': {
        'const': round(float(rlm_fit.params['const']), 4),
        'bmi': round(float(rlm_fit.params['bmi']), 4),
        'bmi2': round(float(rlm_fit.params['bmi2']), 6),
    },
    'risk_weights': {'w_delay': 0.4, 'w_inaccuracy': 0.6},
    'residual_std': round(resid_std, 6),
    'gamm_converged': bool(mdf_gamm.converged),
}
with open(os.path.join(MET_DIR, 'q2_metrics.json'), 'w', encoding='utf-8') as f:
    json.dump(metrics, f, indent=2, ensure_ascii=False)

# ── Run summary ──────────────────────────────────────────────────────
elapsed = round(time.time() - t_start, 2)
run_summary = {
    "question": "Q2", "round": 1,
    "decision_id": "q2_method_choice",
    "methods": [
        {
            "id": "M1", "role": "main_candidate",
            "name": "GAMM连续预测 + 经验达标时间回归 + 风险函数优化 + 数据驱动BMI分组",
            "status": "completed", "time_s": elapsed,
            "key_results": {
                "n_groups": n_groups,
                "breakpoints": [round(b, 1) for b in best_breaks],
                "groups": m1_results,
            },
        },
        {
            "id": "M2", "role": "usable_baseline",
            "name": "临床经验BMI分组 + 组内统计",
            "status": "completed",
            "key_results": {"groups": m2_results},
        },
    ],
    "risk_function": {
        "type": "weighted_sum",
        "components": ["sigmoid_delay_risk", "normal_cdf_inaccuracy_risk"],
        "weights": {"w_delay": 0.4, "w_inaccuracy": 0.6},
        "design": "decisional (risk enters objective, not post-hoc report)",
    },
    "error_sensitivity": {
        "threshold_shift_tested_pct": [0, 0.5, 1.0, 2.0],
        "noise_multipliers_tested": [1.0, 1.5, 2.0, 3.0],
    },
    "seed": SEED, "n_obs": N_OBS, "n_patients": N_PAT,
    "environment": {
        "python": sys.version,
        "statsmodels": sm.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    },
    "inputs": [DATA_PATH],
    "outputs": {
        "figures": [
            'M1_threshold_time_vs_bmi.png', 'M1_risk_function_analysis.png',
            'M1_error_sensitivity.png', 'M1_conc_curves_by_group.png',
            'M1_vs_M2_comparison.png',
        ],
        "tables": ['M1_group_results.csv', 'M2_baseline_results.csv', 'error_sensitivity.csv'],
        "metrics": ['q2_metrics.json'],
    },
    "degeneracy_check": {
        "unique_optimal_times": int(len(np.unique(np.round(opt_times, 2)))),
        "risk_range": [round(float(opt_risks.min()), 4), round(float(opt_risks.max()), 4)],
    },
    "fallback_trigger": {
        "triggered": False,
        "criteria": "T* < 8w or T* > 30w, or no valid breakpoints found",
        "t_star_range": [round(float(t_star_empirical.min()), 2),
                         round(float(t_star_empirical.max()), 2)],
    },
    "warnings": [], "errors": [],
}

with open(os.path.join(OUT_DIR, 'run_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(run_summary, f, indent=2, ensure_ascii=False)

print(f"\nDone in {elapsed:.1f}s. Outputs in {OUT_DIR}")
print(f"   M1: {n_groups} groups, breaks at {[round(b,1) for b in best_breaks]}")
for r in m1_results:
    print(f"     {r['group']}: n={r['n_patients']}, opt={r['optimal_nipt_week']}w, risk={r['min_total_risk']}")
print(f"   M2: {len(m2_results)} clinical groups")
for r in m2_results:
    print(f"     {r['group']}: n={r['n_patients']}, opt={r['optimal_nipt_week']}w, risk={r['min_total_risk']}")
