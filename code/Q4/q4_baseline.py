"""
Q4 Baseline Method: Multivariate Z-value Joint Threshold Rules
Decision: q4_method_choice (M2)

Detects chromosomal aneuploidy using per-chromosome Z-value thresholds.
Any chromosome exceeding threshold => flagged as anomaly.
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc, precision_recall_curve
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

plt.rcParams.update({
    'font.sans-serif': ['SimHei', 'Microsoft YaHei', 'DejaVu Sans'],
    'axes.unicode_minus': False,
    'figure.dpi': 150
})

# === Load Data ===
print("Loading data...")
df = pd.read_csv(DATA_PATH, encoding='utf-8-sig')
LABEL_COL = '染色体非整倍体'
df['is_anomaly'] = df[LABEL_COL].notna().astype(int)
y_true = df['is_anomaly'].values
print(f"  Shape: {df.shape}, Anomaly: {y_true.sum()}, Normal: {(1-y_true).sum()}")

# Z-value columns
Z_COLS = ['13号Z值', '18号Z值', '21号Z值']
ZX_COL = 'X染色体Z值'

# === M2: Z-value Threshold Rules ===

def evaluate_threshold_binary(y_true, y_pred):
    """Evaluate binary predictions."""
    tp = int((y_pred & y_true).sum())
    fp = int((y_pred & (1 - y_true)).sum())
    fn = int(((1 - y_pred) & y_true).sum())
    tn = int(((1 - y_pred) & (1 - y_true)).sum())
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    return {'recall': round(recall, 4), 'precision': round(precision, 4),
            'f1': round(f1, 4), 'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
            'fpr': round(fpr, 4), 'n_flagged': int(y_pred.sum())}


# --- Strategy 1: Absolute Z-value threshold (any chromosome) ---
print("\n=== M2 Strategy 1: |Z| > threshold (any of Z13, Z18, Z21) ===")
thresholds_abs = [1.5, 2.0, 2.5, 3.0, 3.5]
results_abs = []
for thr in thresholds_abs:
    flags = np.zeros(len(df), dtype=int)
    for col in Z_COLS:
        flags |= (np.abs(df[col].values) > thr).astype(int)
    res = evaluate_threshold_binary(y_true, flags)
    res['threshold'] = thr
    res['strategy'] = 'abs_any'
    results_abs.append(res)
    print(f"  |Z|>{thr}: recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

# --- Strategy 2: Directional threshold (Z > thr, since trisomy => positive Z) ---
print("\n=== M2 Strategy 2: Z > threshold (directional, any of Z13, Z18, Z21) ===")
thresholds_dir = [1.0, 1.5, 2.0, 2.5, 3.0]
results_dir = []
for thr in thresholds_dir:
    flags = np.zeros(len(df), dtype=int)
    for col in Z_COLS:
        flags |= (df[col].values > thr).astype(int)
    res = evaluate_threshold_binary(y_true, flags)
    res['threshold'] = thr
    res['strategy'] = 'directional_any'
    results_dir.append(res)
    print(f"  Z>{thr}: recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

# --- Strategy 3: Include X chromosome Z-value ---
print("\n=== M2 Strategy 3: |Z| > threshold (Z13, Z18, Z21, ZX) ===")
Z_COLS_X = Z_COLS + [ZX_COL]
results_with_x = []
for thr in thresholds_abs:
    flags = np.zeros(len(df), dtype=int)
    for col in Z_COLS_X:
        flags |= (np.abs(df[col].values) > thr).astype(int)
    res = evaluate_threshold_binary(y_true, flags)
    res['threshold'] = thr
    res['strategy'] = 'abs_any_with_X'
    results_with_x.append(res)
    print(f"  |Z|>{thr} (+X): recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

# --- Strategy 4: Weighted combination score ---
print("\n=== M2 Strategy 4: Weighted Z-score combination ===")
# Sum of squared Z-values as a continuous score
df['z_score_sum'] = df['13号Z值']**2 + df['18号Z值']**2 + df['21号Z值']**2
df['z_score_sum_x'] = df['z_score_sum'] + df['X染色体Z值']**2

# Normal-sample percentile thresholds
normal_mask = df['is_anomaly'] == 0
percentiles = [50, 75, 80, 85, 90, 95, 97, 99]
results_weighted = []
for pct in percentiles:
    thr = np.percentile(df.loc[normal_mask, 'z_score_sum'].values, pct)
    flags = (df['z_score_sum'].values >= thr).astype(int)
    res = evaluate_threshold_binary(y_true, flags)
    res['threshold'] = round(float(thr), 4)
    res['threshold_type'] = f'p{pct}'
    res['strategy'] = 'z_squared_sum'
    results_weighted.append(res)
    print(f"  Z2sum p{pct} (thr={thr:.3f}): recall={res['recall']:.3f}, precision={res['precision']:.3f}, "
          f"f1={res['f1']:.3f}, fpr={res['fpr']:.3f}, flagged={res['n_flagged']}")

# --- Per-chromosome analysis ---
print("\n=== Per-chromosome detection at |Z|>2 ===")
per_chrom = {}
for col in Z_COLS:
    flags = (np.abs(df[col].values) > 2.0).astype(int)
    res = evaluate_threshold_binary(y_true, flags)
    res['chromosome'] = col
    per_chrom[col] = res
    print(f"  {col} |Z|>2: recall={res['recall']:.3f}, precision={res['precision']:.3f}, flagged={res['n_flagged']}")

# === ROC curve for M2 (using z_score_sum as continuous score) ===
print("\n=== Generating M2 ROC ===")
fpr_m2, tpr_m2, _ = roc_curve(y_true, df['z_score_sum'].values)
auc_m2 = auc(fpr_m2, tpr_m2)
fpr_m2x, tpr_m2x, _ = roc_curve(y_true, df['z_score_sum_x'].values)
auc_m2x = auc(fpr_m2x, tpr_m2x)
print(f"  M2 Z2sum AUC={auc_m2:.4f}, M2 Z2sum+X AUC={auc_m2x:.4f}")

# Figure: M2 ROC
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(fpr_m2, tpr_m2, label=f'M2 Z2-sum (AUC={auc_m2:.3f})', linewidth=2)
ax.plot(fpr_m2x, tpr_m2x, label=f'M2 Z2-sum+X (AUC={auc_m2x:.3f})', lw=2, ls='--')
ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate (Recall)')
ax.set_title('Q4 M2 Z-value Threshold: ROC Curves')
ax.legend(loc='lower right')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_m2_roc.png'))
plt.close(fig)
print("  Saved q4_m2_roc.png")

# Figure: Per-chromosome Z-value distributions
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for idx, col in enumerate(Z_COLS):
    ax = axes[idx]
    ax.hist(df.loc[normal_mask, col], bins=40, alpha=0.6, label='Normal', density=True, color='steelblue')
    ax.hist(df.loc[~normal_mask, col], bins=20, alpha=0.6, label='Anomaly', density=True, color='coral')
    ax.axvline(2.0, color='gray', ls='--', alpha=0.7, label='|Z|=2')
    ax.axvline(-2.0, color='gray', ls='--', alpha=0.7)
    ax.axvline(3.0, color='red', ls='--', alpha=0.5, label='|Z|=3')
    ax.axvline(-3.0, color='red', ls='--', alpha=0.5)
    ax.set_title(col)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
fig.suptitle('Q4 M2: Z-value Distributions by Chromosome', y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'q4_m2_z_distributions.png'), bbox_inches='tight')
plt.close(fig)
print("  Saved q4_m2_z_distributions.png")

# === Save M2 metrics ===
best_abs = max(results_abs, key=lambda x: x['f1'])
best_dir = max(results_dir, key=lambda x: x['f1'])
best_wx = max(results_with_x, key=lambda x: x['f1'])
best_weighted = max(results_weighted, key=lambda x: x['f1'])

m2_metrics = {
    'method': 'M2_z_value_threshold_rules',
    'strategies': {
        'absolute_any': {'best': best_abs, 'all': results_abs},
        'directional_any': {'best': best_dir, 'all': results_dir},
        'absolute_any_with_X': {'best': best_wx, 'all': results_with_x},
        'z_squared_sum': {'best': best_weighted, 'all': results_weighted},
    },
    'per_chromosome': per_chrom,
    'auc_roc_z2sum': round(auc_m2, 4),
    'auc_roc_z2sum_x': round(auc_m2x, 4),
}
with open(os.path.join(METRIC_DIR, 'q4_m2_metrics.json'), 'w', encoding='utf-8') as f:
    json.dump(m2_metrics, f, ensure_ascii=False, indent=2)
print("\nSaved q4_m2_metrics.json")

print("\n=== M2 Summary ===")
print(f"  Best |Z|>thr (any): F1={best_abs['f1']}, recall={best_abs['recall']} at thr={best_abs['threshold']}")
print(f"  Best directional: F1={best_dir['f1']}, recall={best_dir['recall']} at thr={best_dir['threshold']}")
print(f"  Best Z2-sum: F1={best_weighted['f1']}, recall={best_weighted['recall']} at {best_weighted['threshold_type']}")
print(f"  AUC-ROC: Z2sum={auc_m2:.4f}, Z2sum+X={auc_m2x:.4f}")

print("\nM2 done.")
