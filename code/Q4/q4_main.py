"""
Q4 Main Method: Multivariate Mahalanobis Distance Anomaly Detection
Decision: q4_method_choice (M1)

Detects chromosomal aneuploidy in female fetal NIPT data by computing
Mahalanobis distance from the normal-sample covariance baseline.
"""

import os
import sys
import json
import time
import warnings
import numpy as np
import pandas as pd
from scipy.spatial.distance import mahalanobis
from scipy.stats import chi2
from sklearn.metrics import (
    precision_score, recall_score, f1_score,
    roc_curve, auc, precision_recall_curve, confusion_matrix
)
from sklearn.covariance import MinCovDet
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
np.random.seed(42)

# === Paths ===
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_PATH = os.path.join(PROJECT_ROOT, 'workspace', 'data_clean', 'female_cleaned.csv')
OUT_DIR = os.path.join(PROJECT_ROOT, 'results', 'Q4', 'experiments', 'round1')
FIG_DIR = os.path.join(OUT_DIR, 'figures')
TABLE_DIR = os.path.join(OUT_DIR, 'tables')
METRIC_DIR = os.path.join(OUT_DIR, 'metrics')

for d in [FIG_DIR, TABLE_DIR, METRIC_DIR]:
    os.makedirs(d, exist_ok=True)

# === Load Data ===
print("Loading data...")
df = pd.read_csv(DATA_PATH, encoding='utf-8-sig')
print(f"  Shape: {df.shape}")

# Label: AB column (染色体非整倍体) — non-null = anomaly, null = normal
LABEL_COL = '染色体非整倍体'
df['is_anomaly'] = df[LABEL_COL].notna().astype(int)
print(f"  Anomaly distribution: {df['is_anomaly'].value_counts().to_dict()}")
print(f"  Anomaly types: {df[LABEL_COL].value_counts().to_dict()}")

# === Feature sets ===
CORE_FEATURES = ['13号Z值', '18号Z值', '21号Z值', 'X染色体Z值', 'BMI', 'X染色体浓度']
EXTENDED_FEATURES = CORE_FEATURES + [
    'GC含量', '13号GC含量', '18号GC含量', '21号GC含量',
    '原始读段数', '唯一比对读段数', '比对比例', '重复读段比例', '被过滤读段比例'
]

# Drop rows with missing values in core features
df_core = df.dropna(subset=CORE_FEATURES).copy()
print(f"  After dropping NaN in core features: {len(df_core)} rows")

# Split normal / anomaly
normal_mask = df_core['is_anomaly'] == 0
anomaly_mask = df_core['is_anomaly'] == 1
df_normal = df_core[normal_mask]
df_anomaly = df_core[anomaly_mask]
print(f"  Normal: {len(df_normal)}, Anomaly: {len(df_anomaly)}")

# === Mahalanobis Distance Computation ===
def compute_mahalanobis(df_all, df_norm, features, method_name="empirical"):
    """Compute Mahalanobis distance for all samples using normal-sample covariance."""
    X_norm = df_norm[features].values
    X_all = df_all[features].values

    if method_name == "robust":
        # Robust covariance estimation (MinCovDet)
        mcd = MinCovDet(random_state=42)
        mcd.fit(X_norm)
        mean = mcd.location_
        cov_inv = np.linalg.pinv(mcd.covariance_)
    else:
        mean = X_norm.mean(axis=0)
        cov = np.cov(X_norm, rowvar=False)
        cov_inv = np.linalg.pinv(cov)

    distances = np.array([
        mahalanobis(x, mean, cov_inv) for x in X_all
    ])
    return distances


def evaluate_threshold(y_true, distances, threshold):
    """Evaluate binary classification at a given distance threshold."""
    y_pred = (distances >= threshold).astype(int)
    if y_pred.sum() == 0:
        return {'threshold': threshold, 'recall': 0.0, 'precision': 0.0, 'f1': 0.0,
                'tp': 0, 'fp': 0, 'fn': int(y_true.sum()), 'tn': int((1-y_true).sum()),
                'fpr': 0.0, 'n_flagged': 0}
    tp = int((y_pred & y_true).sum())
    fp = int((y_pred & (1 - y_true)).sum())
    fn = int(((1 - y_pred) & y_true).sum())
    tn = int(((1 - y_pred) & (1 - y_true)).sum())
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    return {'threshold': round(threshold, 4), 'recall': round(recall, 4),
            'precision': round(precision, 4), 'f1': round(f1, 4),
            'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
            'fpr': round(fpr, 4), 'n_flagged': int(y_pred.sum())}


# --- Experiment 1: Core features, empirical covariance ---
print("\n=== M1 Experiment 1: Core features, empirical covariance ===")
distances_core = compute_mahalanobis(df_core, df_normal, CORE_FEATURES, "empirical")
df_core['mahal_core'] = distances_core

y_true = df_core['is_anomaly'].values
# Percentile thresholds from NORMAL samples only
normal_dists = distances_core[normal_mask.values]
percentiles = [50, 75, 80, 85, 90, 95, 97, 99]
thresholds_pct = {p: np.percentile(normal_dists, p) for p in percentiles}

# Also try chi2 thresholds (Mahalanobis^2 ~ chi2(p) under normality)
p = len(CORE_FEATURES)
chi2_alphas = [0.10, 0.05, 0.01, 0.005]
thresholds_chi2 = {f"chi2_{a}": np.sqrt(chi2.ppf(1 - a, p)) for a in chi2_alphas}

print("\n  Percentile thresholds (normal sample):")
results_core = []
for pct, thr in thresholds_pct.items():
    res = evaluate_threshold(y_true, distances_core, thr)
    res['threshold_type'] = f'p{pct}'
    results_core.append(res)
    print(f"    p{pct}: thr={thr:.3f}, recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

print("\n  Chi2 thresholds:")
for name, thr in thresholds_chi2.items():
    res = evaluate_threshold(y_true, distances_core, thr)
    res['threshold_type'] = name
    results_core.append(res)
    print(f"    {name}: thr={thr:.3f}, recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

# --- Experiment 2: Core features, robust covariance ---
print("\n=== M1 Experiment 2: Core features, robust covariance ===")
distances_robust = compute_mahalanobis(df_core, df_normal, CORE_FEATURES, "robust")
df_core['mahal_robust'] = distances_robust

normal_dists_r = distances_robust[normal_mask.values]
results_robust = []
for pct in percentiles:
    thr = np.percentile(normal_dists_r, pct)
    res = evaluate_threshold(y_true, distances_robust, thr)
    res['threshold_type'] = f'p{pct}'
    results_robust.append(res)
    print(f"    p{pct}: thr={thr:.3f}, recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

# --- Experiment 3: Extended features ---
print("\n=== M1 Experiment 3: Extended features ===")
df_ext = df_core.dropna(subset=EXTENDED_FEATURES).copy()
normal_ext = df_ext[df_ext['is_anomaly'] == 0]
anomaly_ext = df_ext[df_ext['is_anomaly'] == 1]
print(f"  Extended feature set: {len(df_ext)} rows (normal={len(normal_ext)}, anomaly={len(anomaly_ext)})")

distances_ext = compute_mahalanobis(df_ext, normal_ext, EXTENDED_FEATURES, "empirical")
df_ext['mahal_ext'] = distances_ext
y_true_ext = df_ext['is_anomaly'].values
normal_dists_ext = distances_ext[df_ext['is_anomaly'].values == 0]

results_ext = []
for pct in percentiles:
    thr = np.percentile(normal_dists_ext, pct)
    res = evaluate_threshold(y_true_ext, distances_ext, thr)
    res['threshold_type'] = f'p{pct}'
    results_ext.append(res)
    print(f"    p{pct}: thr={thr:.3f}, recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

# --- Experiment 4: Z-values only (Z13, Z18, Z21) ---
print("\n=== M1 Experiment 4: Z-values only (Z13, Z18, Z21) ===")
Z_ONLY = ['13号Z值', '18号Z值', '21号Z值']
distances_z3 = compute_mahalanobis(df_core, df_normal, Z_ONLY, "empirical")
df_core['mahal_z3'] = distances_z3
normal_dists_z3 = distances_z3[normal_mask.values]

results_z3 = []
for pct in percentiles:
    thr = np.percentile(normal_dists_z3, pct)
    res = evaluate_threshold(y_true, distances_z3, thr)
    res['threshold_type'] = f'p{pct}'
    results_z3.append(res)
    print(f"    p{pct}: thr={thr:.3f}, recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

# --- Per-anomaly-type analysis (best experiment) ---
print("\n=== Per-anomaly-type detection (core features, p85 threshold) ===")
best_thr = np.percentile(normal_dists, 85)
df_core['flagged_core'] = (distances_core >= best_thr).astype(int)
anomaly_types = df_core[df_core['is_anomaly'] == 1].groupby(LABEL_COL).agg(
    total=('is_anomaly', 'count'),
    detected=('flagged_core', 'sum')
).reset_index()
anomaly_types['recall'] = (anomaly_types['detected'] / anomaly_types['total']).round(4)
print(anomaly_types.to_string(index=False))
anomaly_types.to_csv(os.path.join(TABLE_DIR, 'q4_m1_per_type_detection.csv'), index=False)

# === ROC Curves and Figures ===
print("\n=== Generating ROC curves ===")

plt.rcParams.update({
    'font.sans-serif': ['SimHei', 'Microsoft YaHei', 'DejaVu Sans'],
    'axes.unicode_minus': False,
    'figure.dpi': 150
})

fpr_core, tpr_core, _ = roc_curve(y_true, distances_core)
auc_core = auc(fpr_core, tpr_core)
fpr_rob, tpr_rob, _ = roc_curve(y_true, distances_robust)
auc_rob = auc(fpr_rob, tpr_rob)
fpr_z3, tpr_z3, _ = roc_curve(y_true, distances_z3)
auc_z3 = auc(fpr_z3, tpr_z3)
fpr_ext, tpr_ext, _ = roc_curve(y_true_ext, distances_ext)
auc_ext = auc(fpr_ext, tpr_ext)

# Figure 1: ROC comparison
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(fpr_core, tpr_core, label=f'M1 Core 6-feat (AUC={auc_core:.3f})', linewidth=2)
ax.plot(fpr_rob, tpr_rob, label=f'M1 Core Robust (AUC={auc_rob:.3f})', lw=2, ls='--')
ax.plot(fpr_z3, tpr_z3, label=f'M1 Z3-only (AUC={auc_z3:.3f})', lw=2, ls='-.')
ax.plot(fpr_ext, tpr_ext, label=f'M1 Extended 15-feat (AUC={auc_ext:.3f})', lw=2, ls=':')
ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate (Recall)')
ax.set_title('Q4 M1 Mahalanobis Distance: ROC Comparison')
ax.legend(loc='lower right')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_m1_roc_comparison.png'))
plt.close(fig)
print("  Saved q4_m1_roc_comparison.png")

# Figure 2: Distance distribution
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(distances_core[normal_mask.values], bins=50, alpha=0.6, label='Normal', density=True, color='steelblue')
ax.hist(distances_core[anomaly_mask.values], bins=30, alpha=0.6, label='Anomaly', density=True, color='coral')
for pct in [85, 90, 95]:
    thr = np.percentile(normal_dists, pct)
    ax.axvline(thr, color='gray', ls='--', alpha=0.7, label=f'Normal p{pct}={thr:.2f}')
ax.set_xlabel('Mahalanobis Distance')
ax.set_ylabel('Density')
ax.set_title('Q4 M1: Mahalanobis Distance Distribution (Core Features)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_m1_distance_distribution.png'))
plt.close(fig)
print("  Saved q4_m1_distance_distribution.png")

# Figure 3: Precision-Recall curve
prec_curve, rec_curve, _ = precision_recall_curve(y_true, distances_core)
pr_auc = auc(rec_curve, prec_curve)
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(rec_curve, prec_curve, linewidth=2, label=f'M1 Core (PR-AUC={pr_auc:.3f})')
ax.set_xlabel('Recall')
ax.set_ylabel('Precision')
ax.set_title('Q4 M1: Precision-Recall Curve (Core Features)')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_m1_precision_recall.png'))
plt.close(fig)
print("  Saved q4_m1_precision_recall.png")

# Figure 4: Feature contribution scatter
X_norm_vals = df_normal[CORE_FEATURES].values
mean_normal = X_norm_vals.mean(axis=0)
cov_normal = np.cov(X_norm_vals, rowvar=False)
cov_inv_normal = np.linalg.pinv(cov_normal)
feature_contributions = np.zeros((len(df_core), len(CORE_FEATURES)))
X_all = df_core[CORE_FEATURES].values
diff = X_all - mean_normal
for i in range(len(CORE_FEATURES)):
    for j in range(len(CORE_FEATURES)):
        feature_contributions[:, i] += diff[:, i] * diff[:, j] * cov_inv_normal[i, j]
top_contrib_idx = np.argsort(np.abs(feature_contributions).mean(axis=0))[::-1][:2]
feat1, feat2 = CORE_FEATURES[top_contrib_idx[0]], CORE_FEATURES[top_contrib_idx[1]]

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(df_core.loc[normal_mask, feat1], df_core.loc[normal_mask, feat2],
           alpha=0.3, s=20, label='Normal', color='steelblue')
ax.scatter(df_core.loc[anomaly_mask, feat1], df_core.loc[anomaly_mask, feat2],
           alpha=0.7, s=40, label='Anomaly', color='coral', edgecolors='darkred', linewidth=0.5)
ax.set_xlabel(feat1)
ax.set_ylabel(feat2)
ax.set_title('Q4 M1: Top Contributing Features')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_m1_feature_scatter.png'))
plt.close(fig)
print("  Saved q4_m1_feature_scatter.png")

# === Save M1 metrics ===
best_core = max(results_core, key=lambda x: x['f1'])
best_robust = max(results_robust, key=lambda x: x['f1'])
best_ext = max(results_ext, key=lambda x: x['f1'])
best_z3 = max(results_z3, key=lambda x: x['f1'])

m1_metrics = {
    'method': 'M1_mahalanobis_distance',
    'experiments': {
        'core_empirical': {
            'features': CORE_FEATURES, 'n_features': len(CORE_FEATURES),
            'best_threshold': best_core, 'all_thresholds': results_core,
            'auc_roc': round(auc_core, 4), 'auc_pr': round(pr_auc, 4),
        },
        'core_robust': {
            'features': CORE_FEATURES, 'n_features': len(CORE_FEATURES),
            'best_threshold': best_robust, 'all_thresholds': results_robust,
            'auc_roc': round(auc_rob, 4),
        },
        'extended': {
            'features': EXTENDED_FEATURES, 'n_features': len(EXTENDED_FEATURES),
            'best_threshold': best_ext, 'all_thresholds': results_ext,
            'auc_roc': round(auc_ext, 4),
        },
        'z3_only': {
            'features': Z_ONLY, 'n_features': len(Z_ONLY),
            'best_threshold': best_z3, 'all_thresholds': results_z3,
            'auc_roc': round(auc_z3, 4),
        },
    }
}
with open(os.path.join(METRIC_DIR, 'q4_m1_metrics.json'), 'w', encoding='utf-8') as f:
    json.dump(m1_metrics, f, ensure_ascii=False, indent=2)
print("\nSaved q4_m1_metrics.json")

df_core[['is_anomaly', LABEL_COL, 'mahal_core', 'mahal_robust', 'mahal_z3']].to_csv(
    os.path.join(TABLE_DIR, 'q4_m1_distances.csv'), index=False, encoding='utf-8-sig')

print("\n=== M1 Summary ===")
print(f"  Best core empirical: F1={best_core['f1']}, recall={best_core['recall']}, "
      f"precision={best_core['precision']} at {best_core['threshold_type']}")
print(f"  Best robust: F1={best_robust['f1']}, recall={best_robust['recall']} at {best_robust['threshold_type']}")
print(f"  Best extended: F1={best_ext['f1']}, recall={best_ext['recall']} at {best_ext['threshold_type']}")
print(f"  Best Z3-only: F1={best_z3['f1']}, recall={best_z3['recall']} at {best_z3['threshold_type']}")
print(f"  AUC-ROC: core={auc_core:.4f}, robust={auc_rob:.4f}, ext={auc_ext:.4f}, z3={auc_z3:.4f}")

if best_core['recall'] < 0.15 and best_robust['recall'] < 0.15:
    print("\n  WARNING: Fallback trigger CHECK - recall < 15% for both empirical and robust.")
    fallback_triggered = True
else:
    fallback_triggered = False

print("\nM1 done.")
