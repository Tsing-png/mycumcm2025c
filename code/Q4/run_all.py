"""
Q4 Run All: Execute M1 (main) and M2 (baseline), then generate
comparative ROC figure and canonical run_summary.json.
"""

import os
import sys
import json
import time
import subprocess

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(PROJECT_ROOT, 'results', 'Q4', 'experiments', 'round1')
FIG_DIR = os.path.join(OUT_DIR, 'figures')
METRIC_DIR = os.path.join(OUT_DIR, 'metrics')

for d in [FIG_DIR, METRIC_DIR]:
    os.makedirs(d, exist_ok=True)

def run_script(name, path):
    print(f"\n{'='*60}")
    print(f"Running {name}: {path}")
    print('='*60)
    t0 = time.time()
    result = subprocess.run(
        [sys.executable, path],
        capture_output=False, text=True, cwd=PROJECT_ROOT
    )
    elapsed = time.time() - t0
    status = 'success' if result.returncode == 0 else 'failed'
    print(f"\n{name} finished: {status} ({elapsed:.1f}s)")
    return {'status': status, 'elapsed_s': round(elapsed, 2), 'returncode': result.returncode}

# === Run M1 and M2 ===
t_start = time.time()
r1 = run_script('M1_main', os.path.join(CODE_DIR, 'q4_main.py'))
r2 = run_script('M2_baseline', os.path.join(CODE_DIR, 'q4_baseline.py'))
total_time = time.time() - t_start

# === Comparative ROC figure ===
print("\n=== Generating comparative ROC (M1 vs M2) ===")
try:
    import numpy as np
    import pandas as pd
    from sklearn.metrics import roc_curve, auc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        'font.sans-serif': ['SimHei', 'Microsoft YaHei', 'DejaVu Sans'],
        'axes.unicode_minus': False, 'figure.dpi': 150
    })

    DATA_PATH = os.path.join(PROJECT_ROOT, 'workspace', 'data_clean', 'female_cleaned.csv')
    df = pd.read_csv(DATA_PATH, encoding='utf-8-sig')
    LABEL_COL = '染色体非整倍体'
    df['is_anomaly'] = df[LABEL_COL].notna().astype(int)
    y = df['is_anomaly'].values

    # M1: recompute core Mahalanobis distance
    from scipy.spatial.distance import mahalanobis
    CORE = ['13号Z值', '18号Z值', '21号Z值', 'X染色体Z值', 'BMI', 'X染色体浓度']
    df_c = df.dropna(subset=CORE).copy()
    y_c = df_c['is_anomaly'].values
    X_n = df_c.loc[df_c['is_anomaly']==0, CORE].values
    mean_n = X_n.mean(axis=0)
    cov_inv = np.linalg.pinv(np.cov(X_n, rowvar=False))
    d_m1 = np.array([mahalanobis(x, mean_n, cov_inv) for x in df_c[CORE].values])

    # M2: z-squared sum
    df_c['z2sum'] = df_c['13号Z值']**2 + df_c['18号Z值']**2 + df_c['21号Z值']**2

    fpr1, tpr1, _ = roc_curve(y_c, d_m1)
    auc1 = auc(fpr1, tpr1)
    fpr2, tpr2, _ = roc_curve(y_c, df_c['z2sum'].values)
    auc2 = auc(fpr2, tpr2)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(fpr1, tpr1, label=f'M1 Mahalanobis (AUC={auc1:.3f})', linewidth=2)
    ax.plot(fpr2, tpr2, label=f'M2 Z-squared sum (AUC={auc2:.3f})', linewidth=2, ls='--')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate (Recall)')
    ax.set_title('Q4: M1 vs M2 ROC Comparison')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'q4_m1_vs_m2_roc.png'))
    plt.close(fig)
    print("  Saved q4_m1_vs_m2_roc.png")
    comparison_auc = {'m1_mahalanobis': round(auc1, 4), 'm2_z2sum': round(auc2, 4)}
except Exception as e:
    print(f"  Comparison figure failed: {e}")
    comparison_auc = {}

# === Build run_summary.json ===
m1_path = os.path.join(METRIC_DIR, 'q4_m1_metrics.json')
m2_path = os.path.join(METRIC_DIR, 'q4_m2_metrics.json')
m1_met = json.load(open(m1_path, encoding='utf-8')) if os.path.exists(m1_path) else {}
m2_met = json.load(open(m2_path, encoding='utf-8')) if os.path.exists(m2_path) else {}

# Extract best metrics
m1_best = m1_met.get('experiments', {}).get('core_empirical', {}).get('best_threshold', {})
m2_best = m2_met.get('strategies', {}).get('z_squared_sum', {}).get('best', {})

run_summary = {
    'decision_id': 'q4_method_choice',
    'round': 1,
    'methods': {
        'M1_main': {
            'id': 'M1',
            'role': 'main_candidate',
            'name': 'Multivariate Mahalanobis Distance Anomaly Detection',
            'status': r1['status'],
            'elapsed_s': r1['elapsed_s'],
            'best_result': m1_best,
            'auc_roc': m1_met.get('experiments', {}).get('core_empirical', {}).get('auc_roc'),
            'auc_pr': m1_met.get('experiments', {}).get('core_empirical', {}).get('auc_pr'),
        },
        'M2_baseline': {
            'id': 'M2',
            'role': 'usable_baseline',
            'name': 'Multivariate Z-value Joint Threshold Rules',
            'status': r2['status'],
            'elapsed_s': r2['elapsed_s'],
            'best_result': m2_best,
            'auc_roc': m2_met.get('auc_roc_z2sum'),
        }
    },
    'comparison_auc': comparison_auc,
    'inputs': ['workspace/data_clean/female_cleaned.csv'],
    'outputs': {
        'figures': [
            'q4_m1_roc_comparison.png', 'q4_m1_distance_distribution.png',
            'q4_m1_precision_recall.png', 'q4_m1_feature_scatter.png',
            'q4_m2_roc.png', 'q4_m2_z_distributions.png', 'q4_m1_vs_m2_roc.png'
        ],
        'tables': ['q4_m1_distances.csv', 'q4_m1_per_type_detection.csv'],
        'metrics': ['q4_m1_metrics.json', 'q4_m2_metrics.json']
    },
    'seed': 42,
    'environment': {'python': sys.version.split()[0]},
    'total_elapsed_s': round(total_time, 2),
    'degeneracy_check': {
        'all_normal': bool(m1_best.get('recall', 0) == 0),
        'all_anomaly': bool(m1_best.get('fpr', 0) > 0.5),
        'verdict': 'PASS' if m1_best.get('recall', 0) > 0 and m1_best.get('fpr', 1) < 0.5 else 'WARN'
    },
    'fallback_trigger': {
        'condition': 'M1 recall < 15% and FPR > 30% at best threshold',
        'm1_recall': m1_best.get('recall'),
        'triggered': m1_best.get('recall', 0) < 0.15
    },
    'warnings': [],
    'errors': []
}

if r1['status'] == 'failed':
    run_summary['errors'].append('M1 execution failed')
if r2['status'] == 'failed':
    run_summary['errors'].append('M2 execution failed')
if m1_best.get('recall', 0) < 0.15:
    run_summary['warnings'].append(f"M1 best recall={m1_best.get('recall')} < 15% threshold")

with open(os.path.join(OUT_DIR, 'run_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(run_summary, f, ensure_ascii=False, indent=2)

print(f"\n{'='*60}")
print("Q4 Round 1 Complete")
print(f"  Total time: {total_time:.1f}s")
print(f"  M1 status: {r1['status']}, M2 status: {r2['status']}")
print(f"  M1 best: recall={m1_best.get('recall')}, F1={m1_best.get('f1')}")
print(f"  M2 best: recall={m2_best.get('recall')}, F1={m2_best.get('f1')}")
print(f"  Fallback triggered: {run_summary['fallback_trigger']['triggered']}")
print(f"  run_summary.json saved to {OUT_DIR}")
print('='*60)
