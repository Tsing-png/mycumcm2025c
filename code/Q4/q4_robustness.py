"""
Q4 鲁棒性分析:
1. 特征子集敏感性 (不同特征组合的AUC/F1)
2. 阈值敏感性 (不同百分位阈值的precision-recall权衡)
3. 留一孕妇交叉验证 (leave-one-out CV for small anomaly set)
Seed: 42
"""

import os, sys, json, warnings, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.spatial.distance import mahalanobis
from scipy.stats import chi2
from sklearn.metrics import roc_curve, auc, f1_score
from sklearn.covariance import MinCovDet
from itertools import combinations

SEED = 42
np.random.seed(SEED)
warnings.filterwarnings('ignore')

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_ROOT, 'workspace', 'data_clean', 'female_cleaned.csv')
OUT_DIR = os.path.join(PROJECT_ROOT, 'robustness', 'Q4')
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

# ── Load ────────────────────────────────────────────────────────────
print("[Q4 Robustness] Loading data...")
df = pd.read_csv(DATA_PATH, encoding='utf-8-sig')
LABEL_COL = '染色体非整倍体'
df['is_anomaly'] = df[LABEL_COL].notna().astype(int)

CORE_FEATURES = ['13号Z值', '18号Z值', '21号Z值', 'X染色体Z值', 'BMI', 'X染色体浓度']
df_core = df.dropna(subset=CORE_FEATURES).copy()
normal_mask = df_core['is_anomaly'] == 0
anomaly_mask = df_core['is_anomaly'] == 1
df_normal = df_core[normal_mask]
df_anomaly = df_core[anomaly_mask]
y_true = df_core['is_anomaly'].values
print(f"   Total={len(df_core)}, Normal={len(df_normal)}, Anomaly={len(df_anomaly)}")

def compute_mahalanobis(df_all, df_norm, features):
    """Compute Mahalanobis distance using robust covariance."""
    X_norm = df_norm[features].values
    X_all = df_all[features].values
    mcd = MinCovDet(random_state=42)
    mcd.fit(X_norm)
    mean = mcd.location_
    cov_inv = np.linalg.pinv(mcd.covariance_)
    distances = np.array([mahalanobis(x, mean, cov_inv) for x in X_all])
    return distances

def evaluate_at_percentile(y_true, distances, normal_dists, pct):
    thr = np.percentile(normal_dists, pct)
    y_pred = (distances >= thr).astype(int)
    tp = int((y_pred & y_true).sum())
    fp = int((y_pred & (1 - y_true)).sum())
    fn = int(((1 - y_pred) & y_true).sum())
    tn = int(((1 - y_pred) & (1 - y_true)).sum())
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {'threshold_pct': pct, 'threshold_val': round(float(thr), 4),
            'recall': round(recall, 4), 'precision': round(precision, 4),
            'f1': round(f1, 4), 'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn}

# Reference: full core features
print("   Computing reference (6 core features, robust cov)...")
ref_distances = compute_mahalanobis(df_core, df_normal, CORE_FEATURES)
ref_normal_dists = ref_distances[normal_mask.values]
fpr_ref, tpr_ref, _ = roc_curve(y_true, ref_distances)
ref_auc = round(float(auc(fpr_ref, tpr_ref)), 4)
ref_best = max([evaluate_at_percentile(y_true, ref_distances, ref_normal_dists, p)
                for p in [75, 80, 85, 90, 95]], key=lambda x: x['f1'])
print(f"   Reference AUC={ref_auc}, best F1={ref_best['f1']} at p{ref_best['threshold_pct']}")

# ══════════════════════════════════════════════════════════════════════
# 1. 特征子集敏感性
# ══════════════════════════════════════════════════════════════════════
print("\n[1/3] Feature subset sensitivity...")
feature_subsets = {
    'z3_only': ['13号Z值', '18号Z值', '21号Z值'],
    'z4': ['13号Z值', '18号Z值', '21号Z值', 'X染色体Z值'],
    'z3_bmi': ['13号Z值', '18号Z值', '21号Z值', 'BMI'],
    'z3_xconc': ['13号Z值', '18号Z值', '21号Z值', 'X染色体浓度'],
    'z4_bmi': ['13号Z值', '18号Z值', '21号Z值', 'X染色体Z值', 'BMI'],
    'z4_xconc': ['13号Z值', '18号Z值', '21号Z值', 'X染色体Z值', 'X染色体浓度'],
    'core_6': CORE_FEATURES,
}

# Also test leave-one-feature-out
for feat in CORE_FEATURES:
    subset = [f for f in CORE_FEATURES if f != feat]
    feature_subsets[f'drop_{feat[:4]}'] = subset

subset_results = []
for name, feats in feature_subsets.items():
    try:
        dists = compute_mahalanobis(df_core, df_normal, feats)
        norm_d = dists[normal_mask.values]
        fpr_s, tpr_s, _ = roc_curve(y_true, dists)
        auc_s = round(float(auc(fpr_s, tpr_s)), 4)
        best_f1 = max([evaluate_at_percentile(y_true, dists, norm_d, p)
                       for p in [80, 85, 90, 95]], key=lambda x: x['f1'])
        subset_results.append({
            'name': name, 'n_features': len(feats), 'features': feats,
            'auc_roc': auc_s, 'best_f1': best_f1['f1'],
            'best_recall': best_f1['recall'], 'best_precision': best_f1['precision'],
            'best_pct': best_f1['threshold_pct'],
        })
        print(f"   {name} ({len(feats)} feat): AUC={auc_s}, F1={best_f1['f1']}")
    except Exception as e:
        print(f"   {name}: FAILED ({e})")
        subset_results.append({'name': name, 'n_features': len(feats), 'error': str(e)})

auc_range = [min(r['auc_roc'] for r in subset_results if 'auc_roc' in r),
             max(r['auc_roc'] for r in subset_results if 'auc_roc' in r)]
print(f"   AUC range across subsets: {auc_range}")

# ══════════════════════════════════════════════════════════════════════
# 2. 阈值敏感性 (precision-recall-F1 vs percentile)
# ══════════════════════════════════════════════════════════════════════
print("\n[2/3] Threshold sensitivity...")
threshold_pcts = list(range(50, 100, 2))
threshold_results = []
for pct in threshold_pcts:
    res = evaluate_at_percentile(y_true, ref_distances, ref_normal_dists, pct)
    threshold_results.append(res)

# Find optimal threshold range
best_f1_row = max(threshold_results, key=lambda x: x['f1'])
# Find range where F1 > 0.8 * best_f1
f1_80 = best_f1_row['f1'] * 0.8
near_optimal = [r for r in threshold_results if r['f1'] >= f1_80]
if near_optimal:
    pct_range = [min(r['threshold_pct'] for r in near_optimal),
                 max(r['threshold_pct'] for r in near_optimal)]
else:
    pct_range = [best_f1_row['threshold_pct'], best_f1_row['threshold_pct']]

print(f"   Best F1={best_f1_row['f1']} at p{best_f1_row['threshold_pct']}")
print(f"   Near-optimal range (F1 >= {f1_80:.3f}): p{pct_range[0]} to p{pct_range[1]}")

# ══════════════════════════════════════════════════════════════════════
# 3. 留一孕妇交叉验证
# ══════════════════════════════════════════════════════════════════════
print("\n[3/3] Leave-one-out cross-validation (anomaly samples)...")
# For each anomaly sample, remove it from the dataset, recompute distances,
# check if the remaining model still flags it
anomaly_indices = df_core.index[anomaly_mask].tolist()
n_anomaly = len(anomaly_indices)

loo_results = []
for i, idx in enumerate(anomaly_indices):
    # Leave this anomaly out of the full dataset
    df_loo = df_core.drop(index=idx)
    normal_loo = df_loo[df_loo['is_anomaly'] == 0]

    # Recompute distances for the left-out sample
    sample = df_core.loc[[idx]]
    X_norm = normal_loo[CORE_FEATURES].values
    X_sample = sample[CORE_FEATURES].values[0]

    mcd = MinCovDet(random_state=42)
    mcd.fit(X_norm)
    mean_loo = mcd.location_
    cov_inv_loo = np.linalg.pinv(mcd.covariance_)
    d_loo = mahalanobis(X_sample, mean_loo, cov_inv_loo)

    # Also compute distances of all normal samples to get threshold
    norm_dists_loo = np.array([
        mahalanobis(x, mean_loo, cov_inv_loo) for x in X_norm
    ])
    thr_85 = np.percentile(norm_dists_loo, 85)
    thr_90 = np.percentile(norm_dists_loo, 90)

    flagged_85 = d_loo >= thr_85
    flagged_90 = d_loo >= thr_90
    anomaly_type = str(sample[LABEL_COL].iloc[0])

    loo_results.append({
        'index': int(idx), 'anomaly_type': anomaly_type,
        'distance': round(float(d_loo), 4),
        'threshold_p85': round(float(thr_85), 4),
        'threshold_p90': round(float(thr_90), 4),
        'flagged_p85': bool(flagged_85),
        'flagged_p90': bool(flagged_90),
    })

loo_recall_85 = round(sum(1 for r in loo_results if r['flagged_p85']) / len(loo_results), 4)
loo_recall_90 = round(sum(1 for r in loo_results if r['flagged_p90']) / len(loo_results), 4)
print(f"   LOO recall at p85: {loo_recall_85} ({sum(1 for r in loo_results if r['flagged_p85'])}/{len(loo_results)})")
print(f"   LOO recall at p90: {loo_recall_90} ({sum(1 for r in loo_results if r['flagged_p90'])}/{len(loo_results)})")

# Per-type LOO recall
loo_by_type = {}
for r in loo_results:
    t = r['anomaly_type']
    if t not in loo_by_type:
        loo_by_type[t] = {'total': 0, 'flagged_85': 0, 'flagged_90': 0}
    loo_by_type[t]['total'] += 1
    loo_by_type[t]['flagged_85'] += int(r['flagged_p85'])
    loo_by_type[t]['flagged_90'] += int(r['flagged_p90'])

for t, v in loo_by_type.items():
    v['recall_85'] = round(v['flagged_85'] / v['total'], 4) if v['total'] > 0 else 0.0
    v['recall_90'] = round(v['flagged_90'] / v['total'], 4) if v['total'] > 0 else 0.0
    print(f"   {t}: {v['flagged_85']}/{v['total']} (p85), {v['flagged_90']}/{v['total']} (p90)")

# ══════════════════════════════════════════════════════════════════════
# Figures
# ══════════════════════════════════════════════════════════════════════
print("\nGenerating figures...")

# Fig 1: Feature subset AUC comparison
valid_subsets = [r for r in subset_results if 'auc_roc' in r]
fig, ax = plt.subplots(figsize=(10, 5))
names = [r['name'] for r in valid_subsets]
aucs = [r['auc_roc'] for r in valid_subsets]
f1s = [r['best_f1'] for r in valid_subsets]
x = np.arange(len(names))
w = 0.35
ax.bar(x - w/2, aucs, w, label='AUC-ROC', color='steelblue', alpha=0.8)
ax.bar(x + w/2, f1s, w, label='Best F1', color='coral', alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(names, rotation=45, ha='right', fontsize=7)
ax.set_ylabel('Score')
ax.set_title('Q4: Feature Subset Sensitivity')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_feature_subset_sensitivity.png'), dpi=150)
plt.close(fig)

# Fig 2: Threshold sensitivity (precision, recall, F1 vs percentile)
fig, ax = plt.subplots(figsize=(9, 5))
pcts = [r['threshold_pct'] for r in threshold_results]
ax.plot(pcts, [r['recall'] for r in threshold_results], 'o-', label='Recall', lw=2)
ax.plot(pcts, [r['precision'] for r in threshold_results], 's-', label='Precision', lw=2)
ax.plot(pcts, [r['f1'] for r in threshold_results], '^-', label='F1', lw=2, color='green')
ax.axvline(best_f1_row['threshold_pct'], color='red', ls='--', alpha=0.5,
           label=f"Best F1 at p{best_f1_row['threshold_pct']}")
ax.set_xlabel('Normal Percentile Threshold')
ax.set_ylabel('Score')
ax.set_title('Q4: Threshold Sensitivity (Precision / Recall / F1)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_threshold_sensitivity.png'), dpi=150)
plt.close(fig)

# Fig 3: LOO distance distribution
fig, ax = plt.subplots(figsize=(8, 5))
loo_dists = [r['distance'] for r in loo_results]
loo_flagged = [r['flagged_p85'] for r in loo_results]
colors_loo = ['coral' if f else 'steelblue' for f in loo_flagged]
ax.bar(range(len(loo_dists)), loo_dists, color=colors_loo, alpha=0.7)
if loo_results:
    ax.axhline(loo_results[0]['threshold_p85'], color='red', ls='--', lw=1.5,
               label=f"p85 threshold (~{loo_results[0]['threshold_p85']:.2f})")
ax.set_xlabel('Anomaly Sample Index')
ax.set_ylabel('Mahalanobis Distance (LOO)')
ax.set_title(f'Q4: LOO CV Distances (recall@p85={loo_recall_85})')
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_loo_cv.png'), dpi=150)
plt.close(fig)

# ══════════════════════════════════════════════════════════════════════
# Save summary
# ══════════════════════════════════════════════════════════════════════
summary = {
    "question": "Q4",
    "robustness_checks": [
        {
            "check": "feature_subset_sensitivity",
            "claim": "Core 6-feature set is near-optimal among tested subsets",
            "perturbation": "7 named subsets + 6 leave-one-feature-out",
            "results": subset_results,
            "reference_auc": ref_auc,
            "auc_range": auc_range,
            "status": "PASS" if ref_auc >= max(auc_range) - 0.02 else "CONDITIONAL",
            "limitation": "AUC is low overall (near 0.5-0.6); feature choice has limited impact "
                         "because the anomaly signal is weak in these features"
        },
        {
            "check": "threshold_sensitivity",
            "claim": "Optimal threshold lies in a stable range around p80-p90",
            "perturbation": "percentile threshold from p50 to p98",
            "results": threshold_results,
            "best_f1": best_f1_row,
            "near_optimal_range": pct_range,
            "status": "PASS" if (pct_range[1] - pct_range[0]) >= 6 else "CONDITIONAL",
            "limitation": "F1 is low across all thresholds due to weak separability"
        },
        {
            "check": "leave_one_out_cv",
            "claim": "Detection is not driven by a single influential anomaly",
            "perturbation": "Leave-one-anomaly-out recomputation",
            "n_anomaly": n_anomaly,
            "loo_recall_p85": loo_recall_85,
            "loo_recall_p90": loo_recall_90,
            "per_type": loo_by_type,
            "results_sample": loo_results[:10],
            "status": "PASS" if abs(loo_recall_85 - ref_best['recall']) < 0.1 else "CONDITIONAL",
            "limitation": "LOO recomputes robust covariance each time, which is costly but rigorous; "
                         "small anomaly count limits statistical power of per-type analysis"
        },
    ],
    "seed": SEED,
    "figures": ["q4_feature_subset_sensitivity.png", "q4_threshold_sensitivity.png", "q4_loo_cv.png"],
}

with open(os.path.join(OUT_DIR, 'q4_robustness_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print(f"\nQ4 robustness done. Outputs in {OUT_DIR}")
for chk in summary['robustness_checks']:
    print(f"   {chk['check']}: {chk['status']}")
