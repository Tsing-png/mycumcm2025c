"""
Q1: Y染色体浓度与孕周/BMI关系建模
M1 (main): GAMM — statsmodels MixedLM + natural spline basis
M2 (baseline): LME — statsmodels MixedLM, linear fixed effects
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
from scipy import stats
from patsy import dmatrix, build_design_matrices

SEED = 42
np.random.seed(SEED)
warnings.filterwarnings('ignore', category=FutureWarning)

# ── Paths ────────────────────────────────────────────────────────────
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_ROOT, 'workspace', 'data_clean', 'male_cleaned.csv')
OUT_DIR = os.path.join(PROJECT_ROOT, 'results', 'Q1', 'experiments', 'round1')
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

# ── 1. Load ──────────────────────────────────────────────────────────
print("[1/7] Loading data...")
df = pd.read_csv(DATA_PATH, encoding='utf-8-sig')
df = df.rename(columns={
    '检测孕周_数值': 'gest_week', 'BMI': 'bmi',
    'Y染色体浓度': 'y_conc', '孕妇代码': 'patient',
    '年龄': 'age', 'IVF妊娠': 'ivf',
})
df = df[['patient','gest_week','bmi','y_conc','age','ivf']].dropna(
    subset=['gest_week','bmi','y_conc','patient'])
df['patient'] = df['patient'].astype(str)
print(f"   Rows={len(df)}, Patients={df['patient'].nunique()}")

# ── 2. ICC ───────────────────────────────────────────────────────────
print("[2/7] Computing ICC...")
grand_mean = df['y_conc'].mean()
grps = df.groupby('patient')['y_conc']
k0 = len(df) / df['patient'].nunique()
ms_b = sum(grps.count() * (grps.mean() - grand_mean)**2) / (df['patient'].nunique() - 1)
ms_w = sum(grps.apply(lambda x: sum((x - x.mean())**2))) / (len(df) - df['patient'].nunique())
icc = (ms_b - ms_w) / (ms_b + (k0 - 1) * ms_w)
print(f"   ICC = {icc:.4f}")

# ── 3. M2 Baseline: LME ─────────────────────────────────────────────
print("[3/7] Fitting M2 baseline (LME)...")
t0 = time.time()
md_lme = smf.mixedlm("y_conc ~ gest_week + bmi", data=df,
                      groups=df["patient"], re_formula="~gest_week")
mdf_lme = md_lme.fit(reml=True, method='lbfgs')
lme_time = time.time() - t0
print(f"   LME fit in {lme_time:.2f}s, converged={mdf_lme.converged}")
print(mdf_lme.summary())

# ── 4. M1 Main: GAMM (spline + MixedLM) ─────────────────────────────
print("\n[4/7] Fitting M1 main (GAMM via spline basis + MixedLM)...")

# Build spline basis using patsy bs() — much better conditioning than cr()
# bs(df=4) for gest_week, bs(df=3) for bmi
spline_formula_rhs = "bs(gest_week, df=4, include_intercept=False) + bs(bmi, df=3, include_intercept=False)"
X_design = dmatrix(spline_formula_rhs, data=df, return_type='dataframe')
# Save design_info for prediction later
design_info = X_design.design_info

# Rename columns to safe names for statsmodels formula
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
# Drop the patsy intercept col — statsmodels will add its own
if 'sp_intercept' in X_design.columns:
    X_design = X_design.drop(columns=['sp_intercept'])
    safe_cols = [c for c in safe_cols if c != 'sp_intercept']

gw_spline_cols = [c for c in safe_cols if c.startswith('gw_s')]
bmi_spline_cols = [c for c in safe_cols if c.startswith('bmi_s')]

df_gamm = pd.concat([df.reset_index(drop=True),
                      X_design.reset_index(drop=True)], axis=1)

fe_formula = "y_conc ~ " + " + ".join(safe_cols)

# Try multiple optimizers for convergence
t0 = time.time()
mdf_gamm = None
for opt in ['lbfgs', 'powell', 'cg', 'nm']:
    try:
        md_g = smf.mixedlm(fe_formula, data=df_gamm,
                            groups=df_gamm["patient"],
                            re_formula="~gest_week")
        fit_g = md_g.fit(reml=True, method=opt)
        if fit_g.converged:
            mdf_gamm = fit_g
            print(f"   GAMM converged with optimizer '{opt}'")
            break
        elif mdf_gamm is None:
            mdf_gamm = fit_g  # keep best non-converged
    except Exception as e:
        print(f"   optimizer '{opt}' failed: {e}")
        continue

if mdf_gamm is None:
    # Fallback: random intercept only (simpler RE structure)
    print("   All optimizers failed with random slope. Trying random intercept only...")
    md_g = smf.mixedlm(fe_formula, data=df_gamm,
                        groups=df_gamm["patient"])
    mdf_gamm = md_g.fit(reml=True, method='lbfgs')

gamm_time = time.time() - t0
gamm_converged = mdf_gamm.converged
print(f"   GAMM fit in {gamm_time:.2f}s, converged={gamm_converged}")
print(mdf_gamm.summary())

# ── 5. Model comparison ──────────────────────────────────────────────
print("\n[5/7] Model comparison...")

def model_metrics(fit, name):
    ll = fit.llf
    # count fixed-effect params + unique RE variance params
    n_fe = len(fit.fe_params)
    cov_re_flat = fit.cov_re.values[np.triu_indices_from(fit.cov_re.values)]
    k = n_fe + len(cov_re_flat) + 1  # +1 for residual variance
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

met_lme = model_metrics(mdf_lme, 'M2_LME')
met_gamm = model_metrics(mdf_gamm, 'M1_GAMM')
comparison = [met_gamm, met_lme]
print(pd.DataFrame(comparison).to_string(index=False))

lr_stat = 2*(met_gamm['log_likelihood'] - met_lme['log_likelihood'])
lr_df = met_gamm['n_params'] - met_lme['n_params']
lr_p = float('nan')
if lr_df > 0 and lr_stat > 0:
    lr_p = 1 - stats.chi2.cdf(lr_stat, lr_df)
print(f"   LR test: chi2={lr_stat:.2f}, df={lr_df}, p={lr_p:.6f}")

# Save tables
pd.DataFrame(comparison).to_csv(os.path.join(TBL_DIR, 'model_comparison.csv'), index=False)

def save_coef(fit, path):
    d = pd.DataFrame({'coef': fit.fe_params, 'std_err': fit.bse_fe,
                       'z': fit.tvalues, 'p_value': fit.pvalues})
    d.index.name = 'variable'
    d.to_csv(path)
    return d

save_coef(mdf_lme, os.path.join(TBL_DIR, 'M2_LME_coefficients.csv'))
save_coef(mdf_gamm, os.path.join(TBL_DIR, 'M1_GAMM_coefficients.csv'))

with open(os.path.join(MET_DIR, 'model_comparison.json'), 'w', encoding='utf-8') as f:
    json.dump({'models': comparison,
               'lr_test': {'chi2': round(lr_stat,4), 'df': lr_df, 'p_value': round(lr_p,6)},
               'icc': round(icc,4)}, f, indent=2, ensure_ascii=False)

# ── 6. Diagnostic plots ──────────────────────────────────────────────
print("\n[6/7] Generating diagnostic plots...")

def residual_diagnostics(fit, name, tag):
    resid = fit.resid
    fitted = fit.fittedvalues
    # Residuals vs Fitted
    fig, ax = plt.subplots(figsize=(7,5))
    ax.scatter(fitted, resid, alpha=0.3, s=10, edgecolors='none')
    ax.axhline(0, color='red', lw=1, ls='--')
    ax.set_xlabel('Fitted values'); ax.set_ylabel('Residuals')
    ax.set_title(f'{name}: Residuals vs Fitted')
    fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, f'{tag}_resid_vs_fitted.png'), dpi=150); plt.close(fig)
    # QQ
    fig, ax = plt.subplots(figsize=(6,6))
    stats.probplot(resid, dist='norm', plot=ax)
    ax.set_title(f'{name}: Q-Q Plot')
    fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, f'{tag}_qq.png'), dpi=150); plt.close(fig)
    # Histogram
    fig, ax = plt.subplots(figsize=(7,5))
    ax.hist(resid, bins=40, density=True, alpha=0.7, edgecolor='black')
    xr = np.linspace(resid.min(), resid.max(), 200)
    ax.plot(xr, stats.norm.pdf(xr, resid.mean(), resid.std()), 'r-', lw=2)
    ax.set_xlabel('Residuals'); ax.set_ylabel('Density')
    ax.set_title(f'{name}: Residual Distribution')
    fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, f'{tag}_resid_hist.png'), dpi=150); plt.close(fig)
    # Stats
    sw_stat, sw_p = stats.shapiro(resid) if len(resid)<=5000 else (float('nan'), float('nan'))
    return {
        'shapiro_stat': round(float(sw_stat),4), 'shapiro_p': round(float(sw_p),6),
        'skewness': round(float(stats.skew(resid)),4),
        'kurtosis': round(float(stats.kurtosis(resid)),4),
        'resid_mean': round(float(resid.mean()),6), 'resid_std': round(float(resid.std()),6),
    }

diag_lme = residual_diagnostics(mdf_lme, 'M2 LME', 'M2_LME')
diag_gamm = residual_diagnostics(mdf_gamm, 'M1 GAMM', 'M1_GAMM')
print(f"   LME  diag: {diag_lme}")
print(f"   GAMM diag: {diag_gamm}")

# ── Partial effect plots ─────────────────────────────────────────────
print("   Generating partial effect plots...")

gw_grid = np.linspace(df['gest_week'].min(), df['gest_week'].max(), 200)
bmi_grid = np.linspace(df['bmi'].min(), df['bmi'].max(), 200)
bmi_med = df['bmi'].median()
gw_med = df['gest_week'].median()

# Use saved design_info to build prediction matrices (avoids refit-knot error)
pred_gw_df = pd.DataFrame({'gest_week': gw_grid, 'bmi': bmi_med})
pred_bmi_df = pd.DataFrame({'gest_week': gw_med, 'bmi': bmi_grid})

X_pred_gw = build_design_matrices([design_info], pred_gw_df, return_type='dataframe')[0]
X_pred_bmi = build_design_matrices([design_info], pred_bmi_df, return_type='dataframe')[0]

# Rename to safe cols matching model, drop patsy intercept
def rename_pred(X):
    # Map patsy columns back to safe names used in the model
    orig_cols = list(X.columns)
    new = []
    gi, bi = 0, 0
    for c in orig_cols:
        if c == 'Intercept':
            new.append('_patsy_intercept')
        elif 'gest_week' in c:
            new.append(f'gw_s{gi}'); gi += 1
        else:
            new.append(f'bmi_s{bi}'); bi += 1
    X.columns = new
    if '_patsy_intercept' in X.columns:
        X = X.drop(columns=['_patsy_intercept'])
    # Add intercept with the name statsmodels expects
    X.insert(0, 'Intercept', 1.0)
    return X

X_pred_gw = rename_pred(X_pred_gw.copy())
X_pred_bmi = rename_pred(X_pred_bmi.copy())

fe = mdf_gamm.fe_params
pred_y_gw = X_pred_gw[fe.index] @ fe
pred_y_bmi = X_pred_bmi[fe.index] @ fe

# GAMM partial effect: gest_week
fig, ax = plt.subplots(figsize=(8,5))
ax.scatter(df['gest_week'], df['y_conc'], alpha=0.15, s=8, color='gray', label='Observed')
ax.plot(gw_grid, pred_y_gw.values, 'b-', lw=2, label=f'GAMM | BMI={bmi_med:.1f}')
ax.set_xlabel('Gestational Week'); ax.set_ylabel('Y Chromosome Concentration')
ax.set_title('Partial Effect: Gestational Week (GAMM)'); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, 'M1_GAMM_partial_gest_week.png'), dpi=150); plt.close(fig)

# GAMM partial effect: BMI
fig, ax = plt.subplots(figsize=(8,5))
ax.scatter(df['bmi'], df['y_conc'], alpha=0.15, s=8, color='gray', label='Observed')
ax.plot(bmi_grid, pred_y_bmi.values, 'r-', lw=2, label=f'GAMM | GW={gw_med:.1f}')
ax.set_xlabel('BMI'); ax.set_ylabel('Y Chromosome Concentration')
ax.set_title('Partial Effect: BMI (GAMM)'); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, 'M1_GAMM_partial_bmi.png'), dpi=150); plt.close(fig)

# LME vs GAMM comparison
lp = mdf_lme.fe_params
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
y_lme_gw = lp['Intercept'] + lp['gest_week']*gw_grid + lp['bmi']*bmi_med
axes[0].scatter(df['gest_week'], df['y_conc'], alpha=0.15, s=8, color='gray')
axes[0].plot(gw_grid, y_lme_gw, 'g-', lw=2, label='LME linear')
axes[0].plot(gw_grid, pred_y_gw.values, 'b--', lw=2, label='GAMM spline')
axes[0].set_xlabel('Gestational Week'); axes[0].set_ylabel('Y Conc.')
axes[0].set_title('Gest. Week: LME vs GAMM'); axes[0].legend()

y_lme_bmi = lp['Intercept'] + lp['gest_week']*gw_med + lp['bmi']*bmi_grid
axes[1].scatter(df['bmi'], df['y_conc'], alpha=0.15, s=8, color='gray')
axes[1].plot(bmi_grid, y_lme_bmi, 'g-', lw=2, label='LME linear')
axes[1].plot(bmi_grid, pred_y_bmi.values, 'r--', lw=2, label='GAMM spline')
axes[1].set_xlabel('BMI'); axes[1].set_ylabel('Y Conc.')
axes[1].set_title('BMI: LME vs GAMM'); axes[1].legend()
fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, 'M1_vs_M2_partial_effects.png'), dpi=150); plt.close(fig)

# ── 7. Run summary ───────────────────────────────────────────────────
print("\n[7/7] Writing run summary...")

fallback_triggered = (
    abs(diag_gamm['skewness']) > 1.0
    and diag_gamm['shapiro_p'] < 0.001
    and (lr_p < 0.01 if not np.isnan(lr_p) else False)
)

gamm_status = "converged" if gamm_converged else "converged_with_warnings"
warn_list = []
if not gamm_converged:
    warn_list.append("GAMM optimizer did not fully converge; results may be on boundary")

def safe_re_var(fit):
    cre = fit.cov_re
    d = {"intercept": round(float(cre.iloc[0,0]),6)}
    if cre.shape[0] > 1:
        d["gest_week"] = round(float(cre.iloc[1,1]),6)
    return d

run_summary = {
    "question": "Q1", "round": 1,
    "decision_id": "q1_method_choice",
    "methods": [
        {
            "id": "M1", "role": "main_candidate",
            "name": "GAMM (spline basis + MixedLM)",
            "formula": "y_conc ~ bs(gest_week,df=4) + bs(bmi,df=3) + (1+gest_week|patient)",
            "status": gamm_status, "time_s": round(gamm_time,2),
            "metrics": met_gamm, "residual_diagnostics": diag_gamm,
            "fixed_effects": {k: round(float(v),6) for k,v in mdf_gamm.fe_params.items()},
            "random_effects_variance": safe_re_var(mdf_gamm),
        },
        {
            "id": "M2", "role": "usable_baseline",
            "name": "LME (linear MixedLM)",
            "formula": "y_conc ~ gest_week + bmi + (1+gest_week|patient)",
            "status": "converged", "time_s": round(lme_time,2),
            "metrics": met_lme, "residual_diagnostics": diag_lme,
            "fixed_effects": {k: round(float(v),6) for k,v in mdf_lme.fe_params.items()},
            "random_effects_variance": safe_re_var(mdf_lme),
        },
    ],
    "model_comparison": {
        "lr_test": {"chi2": round(lr_stat,4), "df": lr_df, "p_value": round(lr_p,6)},
        "gamm_better_aic": met_gamm['AIC'] < met_lme['AIC'],
        "gamm_better_bic": met_gamm['BIC'] < met_lme['BIC'],
    },
    "icc": round(icc,4),
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
            'M1_GAMM_resid_vs_fitted.png','M1_GAMM_qq.png','M1_GAMM_resid_hist.png',
            'M2_LME_resid_vs_fitted.png','M2_LME_qq.png','M2_LME_resid_hist.png',
            'M1_GAMM_partial_gest_week.png','M1_GAMM_partial_bmi.png',
            'M1_vs_M2_partial_effects.png',
        ],
        "tables": ['model_comparison.csv','M1_GAMM_coefficients.csv','M2_LME_coefficients.csv'],
        "metrics": ['model_comparison.json'],
    },
    "degeneracy_check": {
        "gamm_resid_var": round(float(mdf_gamm.resid.var()),6),
        "lme_resid_var": round(float(mdf_lme.resid.var()),6),
        "gamm_unique_fitted": int(len(np.unique(np.round(mdf_gamm.fittedvalues,6)))),
        "lme_unique_fitted": int(len(np.unique(np.round(mdf_lme.fittedvalues,6)))),
    },
    "fallback_trigger": {
        "triggered": fallback_triggered,
        "criteria": "abs(skew)>1.0 AND shapiro_p<0.001 AND LR_p<0.01",
        "current_skew": diag_gamm['skewness'],
        "current_shapiro_p": diag_gamm['shapiro_p'],
    },
    "warnings": warn_list, "errors": [],
}

with open(os.path.join(OUT_DIR, 'run_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(run_summary, f, indent=2, ensure_ascii=False)

print(f"\nDone. Outputs in {OUT_DIR}")
print(f"   GAMM AIC={met_gamm['AIC']:.2f}  LME AIC={met_lme['AIC']:.2f}")
print(f"   GAMM better AIC: {met_gamm['AIC'] < met_lme['AIC']}")
print(f"   GAMM better BIC: {met_gamm['BIC'] < met_lme['BIC']}")
print(f"   LR test p={lr_p:.6f}")
print(f"   Fallback triggered: {fallback_triggered}")
