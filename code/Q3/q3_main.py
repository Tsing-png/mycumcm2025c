"""
Q3: 多因素GAMM + 残差分解 + 达标比例风险优化 + MC误差传播
M1 (main): 扩展GAMM — y_conc ~ s(孕周) + s(BMI) + age + weight_resid + (1+孕周|孕妇)
           达标比例 P(y>=0.04|t,group) 进入风险函数; MC模拟误差传播
M2 (baseline): y_conc ~ s(孕周) + s(BMI) + age + BMI:age + (1+孕周|孕妇)
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
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats
from scipy.optimize import minimize_scalar
from patsy import dmatrix, build_design_matrices

SEED = 42
np.random.seed(SEED)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ── Paths ────────────────────────────────────────────────────────────
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_ROOT, 'workspace', 'data_clean', 'male_cleaned.csv')
OUT_DIR = os.path.join(PROJECT_ROOT, 'results', 'Q3', 'experiments', 'round1')
FIG_DIR = os.path.join(OUT_DIR, 'figures')
TBL_DIR = os.path.join(OUT_DIR, 'tables')
MET_DIR = os.path.join(OUT_DIR, 'metrics')
for d in [FIG_DIR, TBL_DIR, MET_DIR]:
    os.makedirs(d, exist_ok=True)

# CJK font
for fn in ['SimHei', 'Microsoft YaHei', 'STHeiti']:
    try:
        rcParams['font.sans-serif'] = [fn] + rcParams['font.sans-serif']
        break
    except Exception:
        pass
rcParams['axes.unicode_minus'] = False

# ── 1. Load & prepare ────────────────────────────────────────────────
print("[1/10] Loading data...")
df = pd.read_csv(DATA_PATH, encoding='utf-8-sig')
df = df.rename(columns={
    '检测孕周_数值': 'gest_week', 'BMI': 'bmi',
    'Y染色体浓度': 'y_conc', '孕妇代码': 'patient',
    '年龄': 'age', '身高': 'height', '体重': 'weight',
})
cols_need = ['patient','gest_week','bmi','y_conc','age','height','weight']
df = df[cols_need].dropna(subset=cols_need)
df['patient'] = df['patient'].astype(str)
print(f"   Rows={len(df)}, Patients={df['patient'].nunique()}")
print(f"   Age range: {df['age'].min()}-{df['age'].max()}")
print(f"   Height range: {df['height'].min()}-{df['height'].max()}")
print(f"   Weight range: {df['weight'].min()}-{df['weight'].max()}")
print(f"   BMI range: {df['bmi'].min():.1f}-{df['bmi'].max():.1f}")

# ── 2. Residual decomposition for collinearity ──────────────────────
print("\n[2/10] Residual decomposition: weight_resid = weight - f(BMI, height)...")
# weight ~ BMI + height => residual captures weight info independent of BMI
X_wr = sm.add_constant(df[['bmi', 'height']])
wr_model = sm.OLS(df['weight'], X_wr).fit()
df['weight_resid'] = wr_model.resid
print(f"   weight ~ bmi + height: R2={wr_model.rsquared:.4f}")
print(f"   weight_resid std={df['weight_resid'].std():.3f}, mean={df['weight_resid'].mean():.6f}")

# VIF check: before and after decomposition
print("\n[3/10] VIF check...")
vif_before = pd.DataFrame({
    'variable': ['bmi', 'height', 'weight', 'age'],
})
X_vif_before = sm.add_constant(df[['bmi', 'height', 'weight', 'age']])
vif_before['VIF'] = [variance_inflation_factor(X_vif_before.values, i+1) for i in range(4)]

X_vif_after = sm.add_constant(df[['bmi', 'age', 'weight_resid']])
vif_after = pd.DataFrame({
    'variable': ['bmi', 'age', 'weight_resid'],
    'VIF': [variance_inflation_factor(X_vif_after.values, i+1) for i in range(3)]
})
print("   VIF BEFORE decomposition:")
print(vif_before.to_string(index=False))
print("   VIF AFTER decomposition:")
print(vif_after.to_string(index=False))

vif_before.to_csv(os.path.join(TBL_DIR, 'vif_before_decomposition.csv'), index=False)
vif_after.to_csv(os.path.join(TBL_DIR, 'vif_after_decomposition.csv'), index=False)

# ── 4. M1 Main: Multi-factor GAMM with residual decomposition ────────
print("\n[4/10] Fitting M1 (multi-factor GAMM + residual decomposition)...")

# Build spline basis (same as Q1: bs(gest_week,df=4) + bs(bmi,df=3))
spline_formula_rhs = "bs(gest_week, df=4, include_intercept=False) + bs(bmi, df=3, include_intercept=False)"
X_design = dmatrix(spline_formula_rhs, data=df, return_type='dataframe')
design_info = X_design.design_info

# Rename spline columns to safe names
col_map = {}
safe_cols = []
for i, c in enumerate(X_design.columns):
    if c == 'Intercept':
        safe = 'sp_intercept'
    elif 'gest_week' in c:
        safe = f'gw_s{len([x for x in safe_cols if x.startswith("gw_s")])}'
    else:
        safe = f'bmi_s{len([x for x in safe_cols if x.startswith("bmi_s")])}'
    col_map[c] = safe
    safe_cols.append(safe)
X_design.columns = safe_cols
if 'sp_intercept' in X_design.columns:
    X_design = X_design.drop(columns=['sp_intercept'])
    safe_cols = [c for c in safe_cols if c != 'sp_intercept']

gw_spline_cols = [c for c in safe_cols if c.startswith('gw_s')]
bmi_spline_cols = [c for c in safe_cols if c.startswith('bmi_s')]

df_m1 = pd.concat([df.reset_index(drop=True), X_design.reset_index(drop=True)], axis=1)

# M1 formula: splines + age + weight_resid
fe_cols_m1 = safe_cols + ['age', 'weight_resid']
fe_formula_m1 = "y_conc ~ " + " + ".join(fe_cols_m1)

t0 = time.time()
mdf_m1 = None
for opt in ['lbfgs', 'powell', 'cg', 'nm']:
    try:
        md = smf.mixedlm(fe_formula_m1, data=df_m1,
                         groups=df_m1["patient"], re_formula="~gest_week")
        fit = md.fit(reml=True, method=opt)
        if fit.converged:
            mdf_m1 = fit
            print(f"   M1 converged with optimizer '{opt}'")
            break
        elif mdf_m1 is None:
            mdf_m1 = fit
    except Exception as e:
        print(f"   M1 optimizer '{opt}' failed: {e}")
        continue

if mdf_m1 is None:
    print("   M1 all optimizers failed with random slope. Trying random intercept only...")
    md = smf.mixedlm(fe_formula_m1, data=df_m1, groups=df_m1["patient"])
    mdf_m1 = md.fit(reml=True, method='lbfgs')

m1_time = time.time() - t0
print(f"   M1 fit in {m1_time:.2f}s, converged={mdf_m1.converged}")
print(mdf_m1.summary())

# ── 5. M2 Baseline: GAMM + age + BMI:age interaction ────────────────
print("\n[5/10] Fitting M2 baseline (GAMM + age + BMI:age)...")

# M2 uses same spline basis but adds age and bmi_raw*age interaction
df_m2 = df_m1.copy()
df_m2['bmi_age'] = df_m2['bmi'] * df_m2['age']
fe_cols_m2 = safe_cols + ['age', 'bmi_age']
fe_formula_m2 = "y_conc ~ " + " + ".join(fe_cols_m2)

t0 = time.time()
mdf_m2 = None
for opt in ['lbfgs', 'powell', 'cg', 'nm']:
    try:
        md = smf.mixedlm(fe_formula_m2, data=df_m2,
                         groups=df_m2["patient"], re_formula="~gest_week")
        fit = md.fit(reml=True, method=opt)
        if fit.converged:
            mdf_m2 = fit
            print(f"   M2 converged with optimizer '{opt}'")
            break
        elif mdf_m2 is None:
            mdf_m2 = fit
    except Exception as e:
        print(f"   M2 optimizer '{opt}' failed: {e}")
        continue

if mdf_m2 is None:
    print("   M2 all optimizers failed with random slope. Trying random intercept only...")
    md = smf.mixedlm(fe_formula_m2, data=df_m2, groups=df_m2["patient"])
    mdf_m2 = md.fit(reml=True, method='lbfgs')

m2_time = time.time() - t0
print(f"   M2 fit in {m2_time:.2f}s, converged={mdf_m2.converged}")
print(mdf_m2.summary())

# ── 6. Model comparison ──────────────────────────────────────────────
print("\n[6/10] Model comparison...")

def model_metrics(fit, name):
    ll = fit.llf
    n_fe = len(fit.fe_params)
    cov_re_flat = fit.cov_re.values[np.triu_indices_from(fit.cov_re.values)]
    k = n_fe + len(cov_re_flat) + 1
    n = fit.nobs
    aic = -2*ll + 2*k
    bic = -2*ll + np.log(n)*k
    resid = fit.resid
    rmse = np.sqrt(np.mean(resid**2))
    mae = np.mean(np.abs(resid))
    return {
        'model': name, 'log_likelihood': round(float(ll),2),
        'n_fe': n_fe, 'n_params': int(k),
        'AIC': round(float(aic),2), 'BIC': round(float(bic),2),
        'RMSE': round(float(rmse),6), 'MAE': round(float(mae),6),
    }

met_m1 = model_metrics(mdf_m1, 'M1_MultiGAMM')
met_m2 = model_metrics(mdf_m2, 'M2_BaselineGAMM')
comparison = [met_m1, met_m2]
print(pd.DataFrame(comparison).to_string(index=False))

# LR test M1 vs M2
lr_stat = 2*(met_m1['log_likelihood'] - met_m2['log_likelihood'])
lr_df = abs(met_m1['n_params'] - met_m2['n_params'])
lr_p = float('nan')
if lr_df > 0 and lr_stat > 0:
    lr_p = 1 - stats.chi2.cdf(lr_stat, lr_df)
elif lr_df > 0 and lr_stat <= 0:
    lr_p = 1.0  # M1 not better
print(f"   LR test M1 vs M2: chi2={lr_stat:.2f}, df={lr_df}, p={lr_p:.6f}")

pd.DataFrame(comparison).to_csv(os.path.join(TBL_DIR, 'model_comparison.csv'), index=False)

def save_coef(fit, path):
    d = pd.DataFrame({'coef': fit.fe_params, 'std_err': fit.bse_fe,
                       'z': fit.tvalues, 'p_value': fit.pvalues})
    d.index.name = 'variable'
    d.to_csv(path)
    return d

coef_m1 = save_coef(mdf_m1, os.path.join(TBL_DIR, 'M1_coefficients.csv'))
coef_m2 = save_coef(mdf_m2, os.path.join(TBL_DIR, 'M2_coefficients.csv'))
print("\nM1 coefficients:")
print(coef_m1.to_string())
print("\nM2 coefficients:")
print(coef_m2.to_string())

# Residual diagnostics
def residual_diag(fit):
    resid = fit.resid
    sw_stat, sw_p = stats.shapiro(resid) if len(resid) <= 5000 else (float('nan'), float('nan'))
    return {
        'shapiro_stat': round(float(sw_stat),4), 'shapiro_p': round(float(sw_p),6),
        'skewness': round(float(stats.skew(resid)),4),
        'kurtosis': round(float(stats.kurtosis(resid)),4),
        'resid_mean': round(float(resid.mean()),6), 'resid_std': round(float(resid.std()),6),
    }

diag_m1 = residual_diag(mdf_m1)
diag_m2 = residual_diag(mdf_m2)

with open(os.path.join(MET_DIR, 'model_comparison.json'), 'w', encoding='utf-8') as f:
    json.dump({'models': comparison,
               'lr_test': {'chi2': round(lr_stat,4), 'df': lr_df, 'p_value': round(lr_p,6) if not np.isnan(lr_p) else None},
               'diagnostics': {'M1': diag_m1, 'M2': diag_m2},
               'vif_after': vif_after.to_dict('records')}, f, indent=2, ensure_ascii=False)

# ── 7. Prediction function & BMI grouping + optimal timepoint ────────
print("\n[7/10] Building prediction & optimization pipeline...")

# Training data range for clipping predictions
gw_train_min = float(df['gest_week'].min())
gw_train_max = float(df['gest_week'].max())
bmi_train_min = float(df['bmi'].min())
bmi_train_max = float(df['bmi'].max())

# Vectorized prediction: build spline basis for batch (gw, bmi) pairs
def predict_m1_batch(gw_arr, bmi_arr, age_arr, wr_arr):
    """Predict y_conc for arrays of inputs using M1 fixed effects (population-average)."""
    gw_c = np.clip(gw_arr, gw_train_min, gw_train_max)
    bmi_c = np.clip(bmi_arr, bmi_train_min, bmi_train_max)
    pred_df = pd.DataFrame({'gest_week': gw_c, 'bmi': bmi_c})
    X_sp = build_design_matrices([design_info], pred_df, return_type='dataframe')[0]
    orig = list(X_sp.columns)
    new_names = []
    gi, bi = 0, 0
    for c in orig:
        if c == 'Intercept':
            new_names.append('_patsy_intercept')
        elif 'gest_week' in c:
            new_names.append(f'gw_s{gi}'); gi += 1
        else:
            new_names.append(f'bmi_s{bi}'); bi += 1
    X_sp.columns = new_names
    if '_patsy_intercept' in X_sp.columns:
        X_sp = X_sp.drop(columns=['_patsy_intercept'])
    X_sp.insert(0, 'Intercept', 1.0)
    X_sp['age'] = age_arr
    X_sp['weight_resid'] = wr_arr
    fe = mdf_m1.fe_params
    return (X_sp[fe.index] @ fe).values

# Residual standard deviation for MC simulation
resid_std_m1 = float(mdf_m1.resid.std())
re_intercept_std = float(np.sqrt(mdf_m1.cov_re.iloc[0, 0]))

# BMI grouping: patient-level
df_patient = df.groupby('patient').first().reset_index()
bmi_vals = df_patient['bmi'].values

# Clinical-informed grouping with minimum group size
clinical_bounds = [20, 24, 28, 32, 36, 40]

def make_groups(bounds, bmi_arr, min_size=15):
    valid_bounds = [bounds[0]]
    for b in bounds[1:]:
        cnt = np.sum((bmi_arr >= valid_bounds[-1]) & (bmi_arr < b))
        if cnt >= min_size:
            valid_bounds.append(b)
    valid_bounds.append(np.inf)
    groups = []
    for i in range(len(valid_bounds)-1):
        lo, hi = valid_bounds[i], valid_bounds[i+1]
        mask = (bmi_arr >= lo) & (bmi_arr < hi)
        if mask.sum() > 0:
            groups.append({'lo': lo, 'hi': hi, 'n': int(mask.sum()),
                          'bmi_mean': float(bmi_arr[mask].mean())})
    return groups

groups = make_groups(clinical_bounds, bmi_vals, min_size=10)
print(f"   BMI groups: {len(groups)}")
for g in groups:
    hi_str = f"{g['hi']:.0f}" if g['hi'] < 100 else "inf"
    print(f"     [{g['lo']:.0f}, {hi_str}): n={g['n']}, mean_bmi={g['bmi_mean']:.1f}")

# ── 8. Optimal timepoint via compliance-rate risk function ───────────
print("\n[8/10] Computing optimal NIPT timepoints per group (vectorized)...")

THRESHOLD = 0.04
GW_MIN, GW_MAX = 11.0, 25.0
N_GW = 150  # grid points for compliance curve
gw_fine = np.linspace(GW_MIN, GW_MAX, N_GW)
N_SIM_PER_PATIENT = 50  # MC draws per patient per timepoint

def late_risk_arr(t_arr):
    """Vectorized piecewise late risk."""
    r = np.zeros_like(t_arr, dtype=float)
    m1 = (t_arr > 12) & (t_arr <= 20)
    m2 = t_arr > 20
    r[m1] = (t_arr[m1] - 12) / 8.0 * 0.5
    r[m2] = 0.5 + (t_arr[m2] - 20) / 5.0 * 0.5
    return r

alpha_risk = 0.3
beta_comply = 0.7

# Precompute: for each patient, predict y_conc at every gw in gw_fine (batch call)
# Then add MC noise to get compliance rates
print("   Precomputing patient-level predictions across timepoints...")
n_patients_total = len(df_patient)

# Build large batch: n_patients * N_GW rows
gw_tile = np.tile(gw_fine, n_patients_total)           # repeat gw for each patient
bmi_rep = np.repeat(df_patient['bmi'].values, N_GW)    # repeat each patient's bmi N_GW times
age_rep = np.repeat(df_patient['age'].values, N_GW)
wr_rep  = np.repeat(df_patient['weight_resid'].values, N_GW)

pred_all = predict_m1_batch(gw_tile, bmi_rep, age_rep, wr_rep)
# Reshape to (n_patients, N_GW)
pred_matrix = pred_all.reshape(n_patients_total, N_GW)
print(f"   Prediction matrix shape: {pred_matrix.shape}")

# For each group, compute compliance via vectorized MC
rng = np.random.RandomState(SEED)
optimal_results = []

for g in groups:
    label = f"[{g['lo']:.0f},{g['hi']:.0f})" if g['hi'] < 100 else f"[{g['lo']:.0f},+inf)"
    mask = (df_patient['bmi'].values >= g['lo']) & (df_patient['bmi'].values < g['hi'])
    grp_pred = pred_matrix[mask]  # (n_grp, N_GW)
    n_grp = grp_pred.shape[0]

    # MC: for each patient at each timepoint, draw N_SIM_PER_PATIENT samples
    # y_sim = pred + RE_intercept + residual_noise
    # Shape: (n_grp, N_GW, N_SIM_PER_PATIENT)
    re_draws = rng.normal(0, re_intercept_std, (n_grp, 1, N_SIM_PER_PATIENT))
    noise_draws = rng.normal(0, resid_std_m1, (n_grp, N_GW, N_SIM_PER_PATIENT))
    y_sim = grp_pred[:, :, np.newaxis] + re_draws + noise_draws  # broadcast

    # Compliance rate at each timepoint: fraction of (patient, sim) pairs >= threshold
    comply_matrix = (y_sim >= THRESHOLD).mean(axis=(0, 2))  # shape (N_GW,)

    # Risk curve
    late_r = late_risk_arr(gw_fine)
    risk_curve = alpha_risk * late_r + beta_comply * (1 - comply_matrix)

    idx_opt = np.argmin(risk_curve)
    gw_opt = gw_fine[idx_opt]
    comply_opt = comply_matrix[idx_opt]
    risk_opt = risk_curve[idx_opt]

    optimal_results.append({
        'group': label, 'n': g['n'], 'bmi_mean': round(g['bmi_mean'], 1),
        'optimal_week': round(float(gw_opt), 1),
        'compliance_rate': round(float(comply_opt), 3),
        'min_risk': round(float(risk_opt), 4),
        'risk_curve': risk_curve.tolist(),
        'compliance_curve': comply_matrix.tolist(),
    })
    print(f"   {label}: optimal_week={gw_opt:.1f}, compliance={comply_opt:.1%}, risk={risk_opt:.4f}")

opt_df = pd.DataFrame([{
    'BMI_group': r['group'], 'n': r['n'], 'bmi_mean': r['bmi_mean'],
    'optimal_week': r['optimal_week'],
    'compliance_rate': f"{r['compliance_rate']:.1%}",
    'min_risk': r['min_risk']
} for r in optimal_results])
opt_df.to_csv(os.path.join(TBL_DIR, 'q3_optimal_timepoints.csv'), index=False, encoding='utf-8-sig')
print("\n   Optimal timepoints table:")
print(opt_df.to_string(index=False))

# ── 9. Monte Carlo error propagation ─────────────────────────────────
print("\n[9/10] Monte Carlo error propagation analysis (vectorized)...")

MC_ITERS = 500
sigma_meas = resid_std_m1

mc_results = []
for g_idx, g in enumerate(groups):
    label = optimal_results[g_idx]['group']
    gw_opt_base = optimal_results[g_idx]['optimal_week']

    mask = (df_patient['bmi'].values >= g['lo']) & (df_patient['bmi'].values < g['hi'])
    grp_pred = pred_matrix[mask]  # (n_grp, N_GW)
    n_grp = grp_pred.shape[0]

    # For each MC iteration, vary noise scale and recompute optimal week
    opt_weeks_mc = np.zeros(MC_ITERS)
    comply_mc = np.zeros(MC_ITERS)
    rng_mc = np.random.RandomState(SEED + g_idx)

    for mc_i in range(MC_ITERS):
        noise_scale = rng_mc.uniform(0.5, 1.5) * sigma_meas
        re_draws = rng_mc.normal(0, re_intercept_std, (n_grp, 1))
        noise_draws = rng_mc.normal(0, noise_scale, (n_grp, N_GW))
        y_sim = grp_pred + re_draws + noise_draws  # (n_grp, N_GW)
        comply_curve = (y_sim >= THRESHOLD).mean(axis=0)  # (N_GW,)
        risk_curve = alpha_risk * late_risk_arr(gw_fine) + beta_comply * (1 - comply_curve)
        idx_opt = np.argmin(risk_curve)
        opt_weeks_mc[mc_i] = gw_fine[idx_opt]
        comply_mc[mc_i] = comply_curve[idx_opt]

    mc_result = {
        'group': label,
        'optimal_week_mean': round(float(np.mean(opt_weeks_mc)), 2),
        'optimal_week_std': round(float(np.std(opt_weeks_mc)), 2),
        'optimal_week_ci95': [round(float(np.percentile(opt_weeks_mc, 2.5)), 2),
                              round(float(np.percentile(opt_weeks_mc, 97.5)), 2)],
        'compliance_mean': round(float(np.mean(comply_mc)), 3),
        'compliance_std': round(float(np.std(comply_mc)), 3),
        'compliance_ci95': [round(float(np.percentile(comply_mc, 2.5)), 3),
                           round(float(np.percentile(comply_mc, 97.5)), 3)],
        'sigma_meas': round(sigma_meas, 6),
        'mc_opt_weeks_samples': opt_weeks_mc.tolist(),  # store for histogram
    }
    mc_results.append(mc_result)
    print(f"   {label}: week={mc_result['optimal_week_mean']}+/-{mc_result['optimal_week_std']}, "
          f"CI95={mc_result['optimal_week_ci95']}, "
          f"compliance={mc_result['compliance_mean']:.1%}+/-{mc_result['compliance_std']:.1%}")

mc_df = pd.DataFrame([{
    'BMI_group': r['group'],
    'optimal_week_mean': r['optimal_week_mean'],
    'optimal_week_std': r['optimal_week_std'],
    'CI95_lo': r['optimal_week_ci95'][0],
    'CI95_hi': r['optimal_week_ci95'][1],
    'compliance_mean': f"{r['compliance_mean']:.1%}",
    'compliance_std': f"{r['compliance_std']:.1%}",
} for r in mc_results])
mc_df.to_csv(os.path.join(TBL_DIR, 'q3_mc_error_analysis.csv'), index=False, encoding='utf-8-sig')

# Strip large sample arrays before saving JSON
mc_results_json = [{k: v for k, v in r.items() if k != 'mc_opt_weeks_samples'} for r in mc_results]
with open(os.path.join(MET_DIR, 'mc_error_propagation.json'), 'w', encoding='utf-8') as f:
    json.dump({'n_iterations': MC_ITERS, 'sigma_meas': round(sigma_meas, 6),
               'results': mc_results_json}, f, indent=2, ensure_ascii=False)

# ── 10. Figures ──────────────────────────────────────────────────────
print("\n[10/10] Generating figures...")

# Fig 1: Compliance rate curves by BMI group
fig, ax = plt.subplots(figsize=(10, 6))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
for i, r in enumerate(optimal_results):
    c = colors[i % len(colors)]
    ax.plot(gw_fine, r['compliance_curve'], '-', color=c, lw=2, label=r['group'])
    ax.axvline(r['optimal_week'], color=c, ls='--', alpha=0.5, lw=1)
    ax.plot(r['optimal_week'], r['compliance_rate'], 'o', color=c, ms=8)
ax.axhline(0.9, color='gray', ls=':', lw=1, alpha=0.5)
ax.set_xlabel('Gestational Week')
ax.set_ylabel('Compliance Rate (Y conc >= 4%)')
ax.set_title('Q3: Compliance Rate by BMI Group (Multi-factor GAMM)')
ax.legend(title='BMI Group', loc='lower right')
ax.set_xlim(GW_MIN, GW_MAX)
ax.set_ylim(0, 1.05)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q3_compliance_curves.png'), dpi=150)
plt.close(fig)

# Fig 2: Risk curves by BMI group
fig, ax = plt.subplots(figsize=(10, 6))
for i, r in enumerate(optimal_results):
    c = colors[i % len(colors)]
    ax.plot(gw_fine, r['risk_curve'], '-', color=c, lw=2, label=r['group'])
    ax.plot(r['optimal_week'], r['min_risk'], '*', color=c, ms=12)
ax.set_xlabel('Gestational Week')
ax.set_ylabel('Total Risk')
ax.set_title('Q3: Risk Function by BMI Group (alpha=0.3 late + beta=0.7 non-compliance)')
ax.legend(title='BMI Group', loc='upper left')
ax.set_xlim(GW_MIN, GW_MAX)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q3_risk_curves.png'), dpi=150)
plt.close(fig)

# Fig 3: MC error analysis — optimal week distributions
fig, axes = plt.subplots(1, len(mc_results), figsize=(4*len(mc_results), 4), squeeze=False)
axes = axes[0]
for i, r in enumerate(mc_results):
    ax = axes[i]
    samples = np.array(r['mc_opt_weeks_samples'])
    ci = r['optimal_week_ci95']
    mean_w = r['optimal_week_mean']
    ax.hist(samples, bins=30, density=True, alpha=0.7, color=colors[i % len(colors)], edgecolor='black')
    ax.axvline(mean_w, color='red', lw=2, label=f'Mean={mean_w:.1f}')
    ax.axvline(ci[0], color='red', ls='--', lw=1)
    ax.axvline(ci[1], color='red', ls='--', lw=1)
    ax.set_title(f'BMI {r["group"]}')
    ax.set_xlabel('Optimal Week')
    ax.legend(fontsize=8)
fig.suptitle('Q3: MC Error Propagation — Optimal Timepoint Distributions', y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q3_mc_timepoint_distributions.png'), dpi=150, bbox_inches='tight')
plt.close(fig)

# Fig 4: Residual diagnostics for M1
resid_m1 = mdf_m1.resid
fitted_m1 = mdf_m1.fittedvalues
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].scatter(fitted_m1, resid_m1, alpha=0.3, s=10, edgecolors='none')
axes[0].axhline(0, color='red', lw=1, ls='--')
axes[0].set_xlabel('Fitted'); axes[0].set_ylabel('Residuals')
axes[0].set_title('M1: Residuals vs Fitted')

stats.probplot(resid_m1, dist='norm', plot=axes[1])
axes[1].set_title('M1: Q-Q Plot')

axes[2].hist(resid_m1, bins=40, density=True, alpha=0.7, edgecolor='black')
xr = np.linspace(resid_m1.min(), resid_m1.max(), 200)
axes[2].plot(xr, stats.norm.pdf(xr, resid_m1.mean(), resid_m1.std()), 'r-', lw=2)
axes[2].set_xlabel('Residuals'); axes[2].set_ylabel('Density')
axes[2].set_title('M1: Residual Distribution')
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'M1_diagnostics.png'), dpi=150)
plt.close(fig)

# Fig 5: VIF comparison
fig, ax = plt.subplots(figsize=(8, 5))
x_pos = np.arange(len(vif_before))
x_pos2 = np.arange(len(vif_after))
ax.bar(x_pos - 0.2, vif_before['VIF'].values, 0.35, label='Before decomposition', color='#d62728', alpha=0.8)
# For after, align with matching variables
after_map = {r['variable']: r['VIF'] for _, r in vif_after.iterrows()}
after_vals = [after_map.get(v, 0) for v in vif_before['variable']]
ax.bar(x_pos + 0.2, after_vals, 0.35, label='After decomposition', color='#2ca02c', alpha=0.8)
ax.set_xticks(x_pos)
ax.set_xticklabels(vif_before['variable'])
ax.set_ylabel('VIF')
ax.set_title('Q3: VIF Before vs After Residual Decomposition')
ax.axhline(10, color='red', ls='--', lw=1, label='VIF=10 threshold')
ax.legend()
ax.set_yscale('log')
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q3_vif_comparison.png'), dpi=150)
plt.close(fig)

# ── Run summary ──────────────────────────────────────────────────────
print("\nWriting run summary...")

def safe_re_var(fit):
    cre = fit.cov_re
    d = {"intercept": round(float(cre.iloc[0,0]),6)}
    if cre.shape[0] > 1:
        d["gest_week"] = round(float(cre.iloc[1,1]),6)
    return d

# Fallback evaluation
weight_resid_p = float(mdf_m1.pvalues.get('weight_resid', 1.0))
max_vif_after = float(vif_after['VIF'].max())
fallback_triggered = (max_vif_after > 10) or (weight_resid_p > 0.3 and (lr_p > 0.1 if not np.isnan(lr_p) else True))

warn_list = []
if not mdf_m1.converged:
    warn_list.append("M1 optimizer did not fully converge")
if not mdf_m2.converged:
    warn_list.append("M2 optimizer did not fully converge")
if weight_resid_p > 0.1:
    warn_list.append(f"weight_resid marginally significant (p={weight_resid_p:.4f})")

# Strip large sample arrays before saving run summary
for mc_r in mc_results:
    mc_r.pop('mc_opt_weeks_samples', None)

run_summary = {
    "question": "Q3", "round": 1,
    "decision_id": "q3_method_choice",
    "methods": [
        {
            "id": "M1", "role": "main_candidate",
            "name": "Multi-factor GAMM + residual decomposition",
            "formula": "y_conc ~ bs(gest_week,df=4) + bs(bmi,df=3) + age + weight_resid + (1+gest_week|patient)",
            "status": "converged" if mdf_m1.converged else "converged_with_warnings",
            "time_s": round(m1_time, 2),
            "metrics": met_m1, "residual_diagnostics": diag_m1,
            "fixed_effects": {k: round(float(v),6) for k,v in mdf_m1.fe_params.items()},
            "random_effects_variance": safe_re_var(mdf_m1),
        },
        {
            "id": "M2", "role": "usable_baseline",
            "name": "GAMM + age + BMI:age interaction",
            "formula": "y_conc ~ bs(gest_week,df=4) + bs(bmi,df=3) + age + bmi*age + (1+gest_week|patient)",
            "status": "converged" if mdf_m2.converged else "converged_with_warnings",
            "time_s": round(m2_time, 2),
            "metrics": met_m2, "residual_diagnostics": diag_m2,
            "fixed_effects": {k: round(float(v),6) for k,v in mdf_m2.fe_params.items()},
            "random_effects_variance": safe_re_var(mdf_m2),
        },
    ],
    "model_comparison": {
        "lr_test": {"chi2": round(lr_stat,4), "df": lr_df,
                    "p_value": round(lr_p,6) if not np.isnan(lr_p) else None},
        "m1_better_aic": met_m1['AIC'] < met_m2['AIC'],
        "m1_better_bic": met_m1['BIC'] < met_m2['BIC'],
    },
    "vif_analysis": {
        "before_decomposition": vif_before.to_dict('records'),
        "after_decomposition": vif_after.to_dict('records'),
        "max_vif_after": round(max_vif_after, 2),
    },
    "bmi_groups": [{
        "group": r['group'], "n": r['n'], "bmi_mean": r['bmi_mean'],
        "optimal_week": r['optimal_week'],
        "compliance_rate": r['compliance_rate'],
        "min_risk": r['min_risk'],
    } for r in optimal_results],
    "mc_error_analysis": {
        "n_iterations": MC_ITERS,
        "sigma_meas": round(sigma_meas, 6),
        "results": [{k: v for k, v in r.items()} for r in mc_results],
    },
    "risk_function": {
        "alpha_late_risk": alpha_risk,
        "beta_non_compliance": beta_comply,
        "threshold": THRESHOLD,
        "description": "R(t,g) = 0.3*late_risk(t) + 0.7*(1-compliance_rate(t,g))"
    },
    "seed": SEED, "n_obs": int(len(df)), "n_patients": int(df['patient'].nunique()),
    "environment": {
        "python": sys.version,
        "statsmodels": sm.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    },
    "inputs": [DATA_PATH],
    "outputs": {
        "figures": [
            'q3_compliance_curves.png', 'q3_risk_curves.png',
            'q3_mc_timepoint_distributions.png', 'M1_diagnostics.png',
            'q3_vif_comparison.png',
        ],
        "tables": [
            'model_comparison.csv', 'M1_coefficients.csv', 'M2_coefficients.csv',
            'vif_before_decomposition.csv', 'vif_after_decomposition.csv',
            'q3_optimal_timepoints.csv', 'q3_mc_error_analysis.csv',
        ],
        "metrics": ['model_comparison.json', 'mc_error_propagation.json'],
    },
    "degeneracy_check": {
        "m1_resid_var": round(float(mdf_m1.resid.var()), 6),
        "m2_resid_var": round(float(mdf_m2.resid.var()), 6),
        "m1_unique_fitted": int(len(np.unique(np.round(mdf_m1.fittedvalues, 6)))),
        "m2_unique_fitted": int(len(np.unique(np.round(mdf_m2.fittedvalues, 6)))),
    },
    "fallback_trigger": {
        "triggered": fallback_triggered,
        "criteria": "max_vif_after>10 OR (weight_resid_p>0.3 AND lr_p>0.1)",
        "current_max_vif": round(max_vif_after, 2),
        "current_weight_resid_p": round(weight_resid_p, 4),
    },
    "warnings": warn_list, "errors": [],
}

with open(os.path.join(OUT_DIR, 'run_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(run_summary, f, indent=2, ensure_ascii=False)

print(f"\nDone. Outputs in {OUT_DIR}")
print(f"   M1 AIC={met_m1['AIC']:.2f}  M2 AIC={met_m2['AIC']:.2f}")
print(f"   M1 better AIC: {met_m1['AIC'] < met_m2['AIC']}")
print(f"   M1 better BIC: {met_m1['BIC'] < met_m2['BIC']}")
print(f"   LR test p={lr_p}")
print(f"   Fallback triggered: {fallback_triggered}")
print(f"   weight_resid p={weight_resid_p:.4f}")
print(f"   Max VIF after decomposition: {max_vif_after:.2f}")
