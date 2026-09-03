"""Data audit and cleaning script for 附件.xlsx"""
import pandas as pd
import numpy as np
import re
import json
from pathlib import Path
import sys
import io

# Force UTF-8 output
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

BASE = Path(r"C:\Users\HUAWEI\mycumcm2025c")
RAW = BASE / "workspace" / "data_raw" / "附件.xlsx"
CLEAN_DIR = BASE / "workspace" / "data_clean"
DATA_DIR = BASE / "workspace" / "data"
CLEAN_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Actual column names from the Excel file
# We'll rename them to shorter names for convenience
RENAME_MAP = {
    "序号": "序号",
    "孕妇代码": "孕妇代码",
    "年龄": "年龄",
    "身高": "身高",
    "体重": "体重",
    "末次月经": "末次月经",
    "IVF妊娠": "IVF妊娠",
    "检测日期": "检测日期",
    "检测抽血次数": "检测抽血次数",
    "检测孕周": "检测孕周",
    "孕妇BMI": "BMI",
    "原始读段数": "原始读段数",
    "在参考基因组上比对的比例": "比对比例",
    "重复读段的比例": "重复读段比例",
    "唯一比对的读段数  ": "唯一比对读段数",
    "GC含量": "GC含量",
    "13号染色体的Z值": "13号Z值",
    "18号染色体的Z值": "18号Z值",
    "21号染色体的Z值": "21号Z值",
    "X染色体的Z值": "X染色体Z值",
    "Y染色体的Z值": "Y染色体Z值",
    "Y染色体浓度": "Y染色体浓度",
    "X染色体浓度": "X染色体浓度",
    "13号染色体的GC含量": "13号GC含量",
    "18号染色体的GC含量": "18号GC含量",
    "21号染色体的GC含量": "21号GC含量",
    "被过滤掉读段数的比例": "被过滤读段比例",
    "染色体的非整倍体": "染色体非整倍体",
    "怀孕次数": "怀孕次数",
    "生产次数": "生产次数",
    "胎儿是否健康": "胎儿是否健康",
}

# PLACEHOLDER_FUNCTIONS

def parse_gestational_week(val):
    if pd.isna(val):
        return np.nan
    s = str(val).strip()
    m = re.match(r"(\d+)\s*[wW周]\s*\+?\s*(\d+)", s)
    if m:
        return int(m.group(1)) + int(m.group(2)) / 7.0
    m2 = re.match(r"(\d+)\s*[wW周]$", s)
    if m2:
        return float(m2.group(1))
    try:
        return float(s)
    except Exception:
        return np.nan


def audit_sheet(df, sheet_name, report_lines):
    r = report_lines
    r.append(f"\n{'='*60}")
    r.append(f"Sheet: {sheet_name} ({len(df)} rows, {len(df.columns)} cols)")
    r.append("=" * 60)

    # --- Missing values ---
    r.append("\n## 缺失值统计")
    r.append(f"{'列名':<20} {'缺失数':>8} {'缺失比例':>10}")
    r.append("-" * 42)
    for c in df.columns:
        n_miss = int(df[c].isna().sum())
        if n_miss > 0:
            r.append(f"{c:<20} {n_miss:>8} {n_miss/len(df)*100:>9.1f}%")
    total_miss = int(df.isna().sum().sum())
    r.append(f"总缺失单元格数: {total_miss}")

    # --- Numeric stats ---
    num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    r.append(f"\n## 数值列基本统计 ({len(num_cols)} 列)")
    for c in num_cols:
        s = df[c].dropna()
        if len(s) == 0:
            continue
        r.append(f"\n  {c}:")
        r.append(f"    有效值={len(s)}  均值={s.mean():.4f}  标准差={s.std():.4f}")
        r.append(f"    最小={s.min():.4f}  中位数={s.median():.4f}  最大={s.max():.4f}")

    # --- 检测孕周解析 ---
    r.append("\n## 检测孕周解析")
    df["检测孕周_数值"] = df["检测孕周"].apply(parse_gestational_week)
    parsed = df["检测孕周_数值"].dropna()
    fail_count = len(df) - len(parsed)
    r.append(f"  成功解析: {len(parsed)}/{len(df)}")
    if fail_count > 0:
        failed = df[df["检测孕周_数值"].isna()]["检测孕周"]
        r.append(f"  解析失败: {fail_count}条, 样例: {list(failed.head(5))}")
    if len(parsed) > 0:
        r.append(f"  范围: {parsed.min():.2f} ~ {parsed.max():.2f} 周")
        r.append(f"  均值: {parsed.mean():.2f}, 中位数: {parsed.median():.2f}")
        r.append(f"  <10周: {(parsed < 10).sum()}, 10-14周: {((parsed >= 10) & (parsed <= 14)).sum()}, >14周: {(parsed > 14).sum()}")

    # --- BMI 分布 ---
    r.append("\n## BMI 分布")
    bmi = df["BMI"].dropna()
    if len(bmi) > 0:
        r.append(f"  有效值: {len(bmi)}, 范围: {bmi.min():.2f} ~ {bmi.max():.2f}")
        bins = [0, 18.5, 24, 28, 100]
        labels = ["偏瘦(<18.5)", "正常(18.5-24)", "超重(24-28)", "肥胖(>=28)"]
        cats = pd.cut(bmi, bins=bins, labels=labels, right=False)
        for lab in labels:
            cnt = int((cats == lab).sum())
            r.append(f"  {lab}: {cnt} ({cnt/len(bmi)*100:.1f}%)")

    # --- 同一孕妇多次检测 ---
    r.append("\n## 同一孕妇多次检测")
    vc = df["孕妇代码"].value_counts()
    multi = vc[vc > 1]
    r.append(f"  唯一孕妇数: {len(vc)}")
    r.append(f"  多次检测孕妇数: {len(multi)}")
    if len(multi) > 0:
        r.append(f"  最大检测次数: {multi.max()}")
        freq_dist = vc.value_counts().sort_index()
        r.append(f"  检测次数分布: {dict(zip(freq_dist.index.astype(int), freq_dist.values.astype(int)))}")

    # --- Y染色体浓度 (男胎) ---
    if "男" in sheet_name:
        r.append("\n## Y染色体浓度分布 (男胎)")
        y_conc = df["Y染色体浓度"].dropna()
        if len(y_conc) > 0:
            r.append(f"  有效值: {len(y_conc)}")
            r.append(f"  范围: {y_conc.min():.6f} ~ {y_conc.max():.6f}")
            r.append(f"  均值: {y_conc.mean():.6f}, 中位数: {y_conc.median():.6f}")
            r.append(f"  <=0: {int((y_conc <= 0).sum())}条")
            # percentiles
            pcts = y_conc.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
            r.append(f"  分位数: P5={pcts[0.05]:.6f} P25={pcts[0.25]:.6f} P50={pcts[0.5]:.6f} P75={pcts[0.75]:.6f} P95={pcts[0.95]:.6f}")

    # --- 异常值检测 ---
    r.append("\n## 异常值检测")
    anomalies = {}

    # GC含量 outside 40%-60%
    for gc_col in ["GC含量", "13号GC含量", "18号GC含量", "21号GC含量"]:
        if gc_col in df.columns:
            vals = df[gc_col].dropna()
            out = vals[(vals < 0.40) | (vals > 0.60)]
            if len(out) > 0:
                anomalies[gc_col + "_异常"] = int(len(out))
                r.append(f"  {gc_col} 不在[0.40,0.60]: {len(out)}条")
                r.append(f"    范围: {out.min():.4f} ~ {out.max():.4f}")

    # X染色体浓度负值
    x_conc = df["X染色体浓度"].dropna()
    x_neg_count = int((x_conc < 0).sum())
    if x_neg_count > 0:
        anomalies["X染色体浓度_负值"] = x_neg_count
        neg_vals = x_conc[x_conc < 0]
        r.append(f"  X染色体浓度 < 0: {x_neg_count}条")
        r.append(f"    范围: {neg_vals.min():.6f} ~ {neg_vals.max():.6f}")

    # 比例列 outside [0,1]
    for pct_col in ["比对比例", "重复读段比例", "被过滤读段比例"]:
        if pct_col in df.columns:
            vals = df[pct_col].dropna()
            out = vals[(vals < 0) | (vals > 1)]
            if len(out) > 0:
                anomalies[pct_col + "_越界"] = int(len(out))
                r.append(f"  {pct_col} 不在[0,1]: {len(out)}条")

    # BMI extreme
    if len(bmi) > 0:
        bmi_extreme = bmi[(bmi < 14) | (bmi > 50)]
        if len(bmi_extreme) > 0:
            anomalies["BMI_极端值"] = int(len(bmi_extreme))
            r.append(f"  BMI 极端值(<14或>50): {len(bmi_extreme)}条, 值: {sorted(bmi_extreme.values)}")

    # 年龄
    age = df["年龄"].dropna()
    age_out = age[(age < 15) | (age > 55)]
    if len(age_out) > 0:
        anomalies["年龄_极端值"] = int(len(age_out))
        r.append(f"  年龄极端值(<15或>55): {len(age_out)}条")

    # Z值极端 (|Z| > 5 for chr 13, 18, 21)
    for zc in ["13号Z值", "18号Z值", "21号Z值"]:
        if zc in df.columns:
            zv = df[zc].dropna()
            extreme = zv[zv.abs() > 5]
            if len(extreme) > 0:
                anomalies[zc + "_极端"] = int(len(extreme))
                r.append(f"  {zc} |Z|>5: {len(extreme)}条, 值: {sorted(extreme.values)}")

    if not anomalies:
        r.append("  未发现明显异常值")

    # --- 染色体非整倍体分布 ---
    r.append("\n## 染色体非整倍体分布")
    ab = df["染色体非整倍体"]
    r.append(f"  缺失: {int(ab.isna().sum())}")
    vc_ab = ab.value_counts(dropna=False)
    for val, cnt in vc_ab.items():
        label = str(val) if pd.notna(val) else "NaN"
        r.append(f"  [{label}]: {cnt} ({cnt/len(df)*100:.1f}%)")

    # --- 胎儿是否健康 ---
    r.append("\n## 胎儿是否健康分布")
    ae = df["胎儿是否健康"]
    vc_ae = ae.value_counts(dropna=False)
    for val, cnt in vc_ae.items():
        label = str(val) if pd.notna(val) else "NaN"
        r.append(f"  [{label}]: {cnt} ({cnt/len(df)*100:.1f}%)")

    # --- IVF妊娠分布 ---
    r.append("\n## IVF妊娠分布")
    ivf = df["IVF妊娠"]
    vc_ivf = ivf.value_counts(dropna=False)
    for val, cnt in vc_ivf.items():
        label = str(val) if pd.notna(val) else "NaN"
        r.append(f"  [{label}]: {cnt} ({cnt/len(df)*100:.1f}%)")

    return df, anomalies


# ============================================================
# Main
# ============================================================
print("Reading raw data...")
df_male_raw = pd.read_excel(RAW, sheet_name="男胎检测数据")
df_female_raw = pd.read_excel(RAW, sheet_name="女胎检测数据")

# Strip whitespace from column names
df_male_raw.columns = [c.strip() for c in df_male_raw.columns]
df_female_raw.columns = [c.strip() for c in df_female_raw.columns]

# Female sheet has Unnamed:20, Unnamed:21 for Y染色体Z值, Y染色体浓度 (empty for females)
female_rename_extra = {}
for c in df_female_raw.columns:
    if "Unnamed: 20" in str(c):
        female_rename_extra[c] = "Y染色体的Z值"
    elif "Unnamed: 21" in str(c):
        female_rename_extra[c] = "Y染色体浓度"
if female_rename_extra:
    df_female_raw = df_female_raw.rename(columns=female_rename_extra)

# Build rename map (strip whitespace from keys too)
rename_clean = {}
for k, v in RENAME_MAP.items():
    rename_clean[k.strip()] = v

df_male = df_male_raw.rename(columns=rename_clean)
df_female = df_female_raw.rename(columns=rename_clean)

print(f"男胎: {df_male.shape}, 女胎: {df_female.shape}")
print(f"男胎列: {list(df_male.columns)}")

report = []
report.append("# 数据审计报告")
report.append(f"生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}")
report.append(f"原始文件: workspace/data_raw/附件.xlsx")
report.append(f"男胎记录数: {len(df_male)}, 女胎记录数: {len(df_female)}")

df_male_clean, anom_male = audit_sheet(df_male, "男胎检测数据", report)
df_female_clean, anom_female = audit_sheet(df_female, "女胎检测数据", report)

# --- Add anomaly flag columns ---
for gc_col in ["GC含量", "13号GC含量", "18号GC含量", "21号GC含量"]:
    for df_ in [df_male_clean, df_female_clean]:
        if gc_col in df_.columns:
            df_[f"{gc_col}_异常"] = df_[gc_col].apply(
                lambda x: 1 if pd.notna(x) and (x < 0.40 or x > 0.60) else 0)

for df_ in [df_male_clean, df_female_clean]:
    df_["X染色体浓度_负值标记"] = (df_["X染色体浓度"] < 0).astype(int)

# --- Save cleaned data ---
df_male_clean.to_csv(CLEAN_DIR / "male_cleaned.csv", index=False, encoding="utf-8-sig")
df_female_clean.to_csv(CLEAN_DIR / "female_cleaned.csv", index=False, encoding="utf-8-sig")
print(f"\nCleaned data saved to {CLEAN_DIR}")

# --- Save report ---
report_text = "\n".join(report)
report_path = DATA_DIR / "data_report.md"
with open(report_path, "w", encoding="utf-8") as f:
    f.write(report_text)
print(f"Report saved to {report_path}")

# --- Build data_profile.json ---
def safe_int(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    return x

def safe_dict(d):
    """Recursively convert numpy types for JSON serialization."""
    if isinstance(d, dict):
        return {str(k): safe_dict(v) for k, v in d.items()}
    if isinstance(d, list):
        return [safe_dict(i) for i in d]
    if isinstance(d, (np.integer,)):
        return int(d)
    if isinstance(d, (np.floating,)):
        return float(d)
    if isinstance(d, np.bool_):
        return bool(d)
    if isinstance(d, float) and np.isnan(d):
        return None
    return d

male_miss = {c: int(df_male[c].isna().sum()) for c in df_male.columns if df_male[c].isna().sum() > 0}
female_miss = {c: int(df_female[c].isna().sum()) for c in df_female.columns if df_female[c].isna().sum() > 0}

profile = {
    "schema_version": 1,
    "raw_files": ["workspace/data_raw/附件.xlsx"],
    "attachment_mapping": [
        {"file": "附件.xlsx", "sheet": "男胎检测数据", "mapped_to": "shared", "rows": len(df_male), "cols": 31},
        {"file": "附件.xlsx", "sheet": "女胎检测数据", "mapped_to": "shared", "rows": len(df_female), "cols": 31}
    ],
    "fields": list(df_male.columns[:31]),
    "quality": {
        "missingness": {"male": male_miss, "female": female_miss},
        "duplicates": {
            "male_duplicate_rows": int(df_male.duplicated().sum()),
            "female_duplicate_rows": int(df_female.duplicated().sum())
        },
        "impossible_values": {
            "male": anom_male,
            "female": anom_female
        },
        "outliers": {}
    },
    "coverage": {
        "rows": len(df_male) + len(df_female),
        "male_rows": len(df_male),
        "female_rows": len(df_female),
        "effective_sample_size": {
            "male_unique_patients": int(df_male["孕妇代码"].nunique()),
            "female_unique_patients": int(df_female["孕妇代码"].nunique())
        },
        "time_range": None,
        "time_gaps": None
    },
    "distribution_risks": {
        "class_imbalance": {
            "male_vs_female_ratio": f"{len(df_male)}:{len(df_female)} ({len(df_male)/len(df_female):.2f}:1)"
        },
        "rare_categories": [],
        "high_cardinality": [],
        "redundancy_warnings": [],
        "concentration_metrics": {}
    },
    "per_question_readiness": {
        "shared": "ready_with_warnings",
        "notes": "GC content anomalies and negative X chromosome concentrations flagged; gestational week parsed; column names standardized"
    },
    "cleaned_files": [
        "workspace/data_clean/male_cleaned.csv",
        "workspace/data_clean/female_cleaned.csv"
    ],
    "unresolved_risks": [
        "Some GC content values outside expected 40-60% range - flagged, not removed",
        "Negative X chromosome concentrations present - flagged, not removed",
        "Male-to-female sample ratio ~1.8:1 imbalance"
    ]
}

profile = safe_dict(profile)

with open(DATA_DIR / "data_profile.json", "w", encoding="utf-8") as f:
    json.dump(profile, f, ensure_ascii=False, indent=2, default=str)
print(f"Profile saved to {DATA_DIR / 'data_profile.json'}")

# Print report
print("\n" + "=" * 60)
print(report_text)
