import io
import os
import re
import subprocess
from pathlib import Path
from glob import glob

import pandas as pd
import numpy as np

try:
    import cupy as cp
except ModuleNotFoundError:
    cp = None

from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score, f1_score
from sklearn.calibration import calibration_curve
from sklearn.metrics import confusion_matrix

import matplotlib.pyplot as plt
from matplotlib import font_manager
from xgboost import XGBClassifier


def cuda_available() -> bool:
    if cp is None:
        return False
    try:
        return cp.cuda.runtime.getDeviceCount() > 0
    except Exception:
        return False


USE_CUDA = cuda_available()


def backend_asarray(a, dtype=None):
    # Keep host-side numpy arrays and let XGBoost move data to cuda:0.
    # CuPy in this env is cuda12x and requires libnvrtc.so.12, while ctl206 has
    # CUDA 13.1 system libraries. XGBoost's numpy input path works on the GPU.
    return np.asarray(a, dtype=dtype)


def backend_to_numpy(a):
    return np.asarray(a)


def configure_plot_font():
    font_candidates = [
        Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
        Path("/usr/share/fonts/truetype/msttcorefonts/Arial.TTF"),
        Path("/usr/share/fonts/truetype/msttcorefonts/arial.ttf"),
        Path.home() / ".local/share/fonts/arial-corefonts/Arial.TTF",
    ]

    font_family = "DejaVu Sans"
    for font_path in font_candidates:
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
            font_family = font_manager.FontProperties(fname=str(font_path)).get_name()
            break

    plt.rcParams["font.family"] = font_family
    plt.rcParams["font.sans-serif"] = [font_family]


configure_plot_font()


# =========================
# 0) Paths
# =========================
DATASET_NAME = "WirbelJ_2019"

ROOT_DATA5 = "/path/to/data2/CRC/CCDC2"
ROOT_DATA = "/path/to/data2/CRC/CCDC2_val"
ROOT_DATA2 = "/path/to/data3/CRC_DATA/CCDC2_data"

DATA_ROOT = os.path.join(ROOT_DATA2, DATASET_NAME)
DATA_DIR = os.path.join(ROOT_DATA, DATASET_NAME, "DiTASiC", "abundance")
SLURM_LOG_DIR = os.path.join(ROOT_DATA, DATASET_NAME, "slurm_out_optimized")
RUN_GROUP_FILE = os.path.join(DATA_ROOT, "PRJEB27928_PE_runs_grouped.txt")
GROUP_FILE = f"{DATASET_NAME}_CRC_Group.txt"

BESTK_TXT = os.path.join(ROOT_DATA5, "ML_results", "bestK_features_fullrange_xgb_gpu.txt")
USE_PVAL_FILTER = False

OUT_DIR = os.path.join(ROOT_DATA, "ML_results")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PREFIX = os.path.join(OUT_DIR, f"{DATASET_NAME}_bestK")
READS_CACHE_FILE = os.path.join(OUT_DIR, f"{DATASET_NAME}_sample_total_reads.tsv")

SEQKIT_BIN = "/path/to/conda/envs/CRC/bin/seqkit"
SEQKIT_THREADS = 8


def split_runs(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(";") if item.strip()]


def fastq_pair_exists(run: str) -> bool:
    return os.path.isfile(os.path.join(DATA_ROOT, run, f"{run}_1.fastq.gz")) and os.path.isfile(
        os.path.join(DATA_ROOT, run, f"{run}_2.fastq.gz")
    )


def normalize_group(group: str) -> str:
    group_norm = str(group).strip()
    return "CRC" if group_norm.lower() == "crc" else "control"


def load_group_by_run(path: str) -> dict[str, str]:
    group_df_raw = pd.read_csv(path, sep="\t", header=0)
    group_df_raw.columns = ["Sample", "Group"]

    group_by_run = {}
    for _, row in group_df_raw.iterrows():
        group = normalize_group(row["Group"])
        for run in split_runs(row["Sample"]):
            previous = group_by_run.get(run)
            if previous is not None and previous != group:
                raise ValueError(f"Conflicting group labels for {run}: {previous} vs {group}")
            group_by_run[run] = group
    return group_by_run


def build_sample_table() -> pd.DataFrame:
    group_by_run = load_group_by_run(GROUP_FILE)
    rows = []

    with open(RUN_GROUP_FILE, "r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            pe_runs = split_runs(line)
            if not pe_runs:
                continue

            downloaded_runs = [run for run in pe_runs if fastq_pair_exists(run)]
            if not downloaded_runs:
                continue

            missing_labels = [run for run in downloaded_runs if run not in group_by_run]
            if missing_labels:
                raise ValueError(
                    f"Missing group labels for downloaded PE runs on line {line_no}: "
                    + ";".join(missing_labels)
                )

            labels = {group_by_run[run] for run in downloaded_runs}
            if len(labels) != 1:
                raise ValueError(
                    f"Downloaded PE runs on line {line_no} have inconsistent labels: "
                    + ";".join(downloaded_runs)
                )

            rows.append(
                {
                    "Sample": downloaded_runs[0],
                    "Group": labels.pop(),
                    "Runs": ";".join(downloaded_runs),
                    "NumRuns": len(downloaded_runs),
                }
            )

    if not rows:
        raise ValueError("No downloaded paired-end WirbelJ_2019 sample groups were found.")

    sample_table = pd.DataFrame(rows)
    if sample_table["Sample"].duplicated().any():
        dupes = sample_table.loc[sample_table["Sample"].duplicated(), "Sample"].tolist()
        raise ValueError(f"Duplicate sample IDs after PE filtering: {dupes}")

    sample_table["label"] = (sample_table["Group"] == "CRC").astype(int)
    sample_table = sample_table.set_index("Sample")

    print(f"[INFO] WirbelJ sample groups: {sample_table.shape[0]}")
    print(f"[INFO] Label counts: {sample_table['Group'].value_counts().to_dict()}")
    print(f"[INFO] Multi-run sample groups: {(sample_table['NumRuns'] > 1).sum()}")

    return sample_table


def normalize_feat(x: str) -> str:
    x = os.path.basename(str(x).strip())
    for suf in (".fa.gz", ".fasta.gz", ".fa", ".fasta"):
        if x.endswith(suf):
            x = x[: -len(suf)]
            break
    return x


def fastq_paths_for_table(sample_table: pd.DataFrame) -> tuple[list[str], dict[str, str]]:
    paths = []
    run_to_sample = {}
    for sample_id, row in sample_table.iterrows():
        for run in split_runs(row["Runs"]):
            run_to_sample[run] = sample_id
            paths.append(os.path.join(DATA_ROOT, run, f"{run}_1.fastq.gz"))
            paths.append(os.path.join(DATA_ROOT, run, f"{run}_2.fastq.gz"))
    return paths, run_to_sample


def parse_int_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(",", "", regex=False).astype(int)


def read_cached_total_reads(sample_table: pd.DataFrame) -> pd.Series | None:
    if not os.path.exists(READS_CACHE_FILE):
        return None

    cache = pd.read_csv(READS_CACHE_FILE, sep="\t")
    if not {"Sample", "num_reads"}.issubset(cache.columns):
        return None

    reads = cache.set_index("Sample")["num_reads"].astype(float)
    if set(sample_table.index).issubset(set(reads.index)):
        print(f"[INFO] Loaded cached sample read counts: {READS_CACHE_FILE}")
        return reads.loc[sample_table.index]
    return None


def read_total_reads_from_kneaddata_logs(sample_table: pd.DataFrame) -> pd.Series | None:
    if not os.path.isdir(SLURM_LOG_DIR):
        return None

    reads_by_sample = {}
    missing = []

    for sample_id in sample_table.index:
        log_paths = sorted(
            Path(SLURM_LOG_DIR).glob(f"{sample_id}.*.kneaddata.log"),
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
        if not log_paths:
            missing.append(sample_id)
            continue

        raw_pair_counts = {}
        for line in log_paths[0].read_text(errors="replace").splitlines():
            match = re.search(r"READ COUNT: raw pair([12]) .*: ([0-9.]+)$", line)
            if match:
                raw_pair_counts[match.group(1)] = float(match.group(2))

        if {"1", "2"}.issubset(raw_pair_counts):
            reads_by_sample[sample_id] = raw_pair_counts["1"] + raw_pair_counts["2"]
        else:
            missing.append(sample_id)

    if missing:
        print(
            "[WARN] Could not read raw read counts from kneaddata logs for "
            f"{len(missing)} samples; falling back to seqkit. Examples: {missing[:10]}"
        )
        return None

    reads = pd.Series(reads_by_sample, dtype=float).loc[sample_table.index]
    reads.name = "num_reads"
    reads.reset_index().rename(columns={"index": "Sample"}).to_csv(
        READS_CACHE_FILE,
        sep="\t",
        index=False,
    )
    print(f"[INFO] Loaded sample read counts from kneaddata logs: {SLURM_LOG_DIR}")
    print(f"[INFO] Saved sample read count cache: {READS_CACHE_FILE}")
    return reads


def count_reads_with_seqkit(sample_table: pd.DataFrame) -> pd.Series:
    fastq_paths, run_to_sample = fastq_paths_for_table(sample_table)
    missing = [path for path in fastq_paths if not os.path.isfile(path)]
    if missing:
        raise FileNotFoundError(f"Missing FASTQ files, first examples: {missing[:5]}")

    if not os.path.isfile(SEQKIT_BIN):
        raise FileNotFoundError(f"seqkit not found: {SEQKIT_BIN}")

    print(f"[INFO] Computing sample read counts with seqkit for {len(fastq_paths)} FASTQ files")
    completed = subprocess.run(
        [SEQKIT_BIN, "stats", "-T", "-j", str(SEQKIT_THREADS), *fastq_paths],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    seqkit_df = pd.read_csv(io.StringIO(completed.stdout), sep="\t")
    seqkit_df["run"] = seqkit_df["file"].map(lambda p: Path(str(p)).parent.name)
    seqkit_df["Sample"] = seqkit_df["run"].map(run_to_sample)
    seqkit_df["num_seqs"] = parse_int_series(seqkit_df["num_seqs"])

    reads = seqkit_df.groupby("Sample")["num_seqs"].sum().astype(float)
    reads = reads.loc[sample_table.index]
    reads.rename("num_reads").reset_index().to_csv(READS_CACHE_FILE, sep="\t", index=False)
    print(f"[INFO] Saved sample read count cache: {READS_CACHE_FILE}")
    return reads


sample_table = build_sample_table()


# =========================
# 0.5) bestK features
# =========================
if not os.path.isfile(BESTK_TXT):
    raise FileNotFoundError(f"[ERROR] BESTK_TXT not found: {BESTK_TXT}")

best_raw = []
with open(BESTK_TXT, "r", encoding="utf-8") as f:
    for line in f:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        best_raw.append(s.split()[0])

seen = set()
best_raw_uniq = []
for x in best_raw:
    if x not in seen:
        seen.add(x)
        best_raw_uniq.append(x)

best_norm2name = {}
best_norm_list = []
for x in best_raw_uniq:
    nx = normalize_feat(x)
    if nx not in best_norm2name:
        best_norm2name[nx] = x
        best_norm_list.append(nx)

print(f"[INFO] bestK features loaded: {len(best_raw_uniq)} (unique)")


# =========================
# 1) Read counts
# =========================
total_reads_per_sample = read_cached_total_reads(sample_table)
if total_reads_per_sample is None:
    total_reads_per_sample = read_total_reads_from_kneaddata_logs(sample_table)
if total_reads_per_sample is None:
    total_reads_per_sample = count_reads_with_seqkit(sample_table)


# =========================
# 2) Abundance matrix
# =========================
file_list = glob(os.path.join(DATA_DIR, "abundance_*.txt"))

file_samples = set(
    os.path.basename(p).replace("abundance_", "").replace(".txt", "")
    for p in file_list
)
group_samples = set(sample_table.index)
missing_abundance_files = sorted(group_samples - file_samples)
print(f"[INFO] group samples without abundance files: {len(missing_abundance_files)}")
print(f"[INFO] examples: {missing_abundance_files[:20]}")

sample_data = []
kept_samples = []

skipped_not_in_group = []
skipped_no_total_reads = []
skipped_empty_after_filter = []

feature_present_counter = {k: 0 for k in best_raw_uniq}

for file_path in file_list:
    sample_id = os.path.basename(file_path).replace("abundance_", "").replace(".txt", "")

    if sample_id not in sample_table.index:
        skipped_not_in_group.append(sample_id)
        continue

    if sample_id not in total_reads_per_sample.index:
        skipped_no_total_reads.append(sample_id)
        continue

    df = pd.read_csv(file_path, sep="\t")
    df["_norm"] = df["taxa.name"].astype(str).map(normalize_feat)
    df = df[df["_norm"].isin(best_norm_list)]

    if USE_PVAL_FILTER and "raw.pval" in df.columns:
        df = df[df["raw.pval"] <= 0.05]

    if df.empty:
        skipped_empty_after_filter.append(sample_id)
        continue

    df["Feature"] = df["_norm"].map(best_norm2name)

    for feat in df["Feature"].unique().tolist():
        if feat in feature_present_counter:
            feature_present_counter[feat] += 1

    s = df.set_index("Feature")["count.estimate"]
    s.name = sample_id

    sample_data.append(s)
    kept_samples.append(sample_id)

print(f"[INFO] kept (used) samples: {len(kept_samples)}")
print(f"[INFO] skipped_not_in_group: {len(skipped_not_in_group)}")
print(f"[INFO] skipped_no_total_reads: {len(skipped_no_total_reads)}")
print(f"[INFO] skipped_empty_after_filter: {len(skipped_empty_after_filter)}")

if not sample_data:
    raise ValueError("No usable abundance files were found for WirbelJ_2019.")

abundance_matrix = pd.DataFrame(sample_data).fillna(0)
abundance_matrix = abundance_matrix.reindex(columns=best_raw_uniq, fill_value=0)

never_seen = [k for k, c in feature_present_counter.items() if c == 0]
never_seen_path = f"{OUT_PREFIX}_features_never_seen.txt"
with open(never_seen_path, "w", encoding="utf-8") as f:
    for k in never_seen:
        f.write(k + "\n")
print(f"[INFO] bestK features never seen in any sample: {len(never_seen)}")
print(f"[INFO] Saved never-seen feature list to: {never_seen_path}")

reads_aligned = total_reads_per_sample.loc[abundance_matrix.index].astype(float)

valid_reads_mask = reads_aligned > 0
if (~valid_reads_mask).any():
    bad_samples = reads_aligned.index[~valid_reads_mask].tolist()
    print(f"[WARN] reads<=0 samples will be removed: {len(bad_samples)} (first 20: {bad_samples[:20]})")
    abundance_matrix = abundance_matrix.loc[valid_reads_mask]
    reads_aligned = reads_aligned.loc[valid_reads_mask]

abundance_rel = abundance_matrix.div(reads_aligned, axis=0)
abundance_rel = abundance_rel.dropna(axis=0, how="all")
abundance_rel = abundance_rel.loc[(abundance_rel.fillna(0) != 0).any(axis=1)]

ordered_samples = [sample for sample in sample_table.index if sample in abundance_rel.index]
abundance_rel = abundance_rel.loc[ordered_samples]
labels = sample_table.loc[abundance_rel.index, "label"]

X = abundance_rel.to_numpy(dtype=np.float32)
y = labels.to_numpy(dtype=int)

feature_names = abundance_rel.columns.tolist()
n_features = len(feature_names)

print(f"[INFO] Dataset: {DATASET_NAME}")
print(f"[INFO] Samples: {X.shape[0]}, Features: {X.shape[1]}")
print(f"[INFO] XGBoost backend: {'cuda:0' if USE_CUDA else 'cpu'}")


# =========================
# 3) Leave-One-Out CV
# =========================
loo = LeaveOneOut()

y_true, y_scores = [], []
train_accuracies, val_accuracies = [], []

importance_sum = np.zeros(n_features, dtype=float)
importance_sq_sum = np.zeros(n_features, dtype=float)
n_folds = 0

for train_idx, test_idx in loo.split(X, y):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    X_train_backend = backend_asarray(X_train, dtype=np.float32)
    y_train_backend = backend_asarray(y_train)
    X_test_backend = backend_asarray(X_test, dtype=np.float32)

    clf = XGBClassifier(
        eval_metric="logloss",
        random_state=42,
        tree_method="hist",
        device="cuda:0" if USE_CUDA else "cpu",
        importance_type="gain",
    )

    clf.fit(X_train_backend, y_train_backend)

    fold_imp = np.asarray(clf.feature_importances_, dtype=float)
    if fold_imp.shape[0] != n_features:
        raise ValueError(f"Feature importance length mismatch: got {fold_imp.shape[0]}, expected {n_features}")

    importance_sum += fold_imp
    importance_sq_sum += fold_imp ** 2
    n_folds += 1

    train_pred = backend_to_numpy(clf.predict(X_train_backend))
    test_pred = backend_to_numpy(clf.predict(X_test_backend))

    train_accuracy = accuracy_score(y_train, train_pred)
    val_accuracy = accuracy_score(y_test, test_pred)

    train_accuracies.append(train_accuracy)
    val_accuracies.append(val_accuracy)

    prob = backend_to_numpy(clf.predict_proba(X_test_backend)[:, 1])[0]
    y_scores.append(float(prob))
    y_true.append(int(y_test[0]))


# =========================
# 4) Mean feature importance
# =========================
importance_mean = importance_sum / max(n_folds, 1)
importance_var = (importance_sq_sum / max(n_folds, 1)) - (importance_mean ** 2)
importance_var = np.maximum(importance_var, 0)
importance_std = np.sqrt(importance_var)

importance_df = pd.DataFrame(
    {
        "Feature": feature_names,
        "MeanImportance": importance_mean,
        "StdImportance": importance_std,
    }
)

importance_df["ImportanceRank"] = importance_df["MeanImportance"].rank(ascending=False, method="min").astype(int)
importance_df = importance_df.sort_values("MeanImportance", ascending=False)

importance_outfile = f"{OUT_PREFIX}_feature_mean_importance.csv"
importance_df.to_csv(importance_outfile, index=False)

print(f"[INFO] Saved mean feature importance to: {importance_outfile}")
print("[INFO] Top 20 features by mean importance:")
print(importance_df.head(20))


# =========================
# 5) ROC / AUC + 95% CI
# =========================
fpr, tpr, _ = roc_curve(y_true, y_scores)
roc_auc = auc(fpr, tpr)

n_bootstraps = 1000
rng = np.random.RandomState(42)
bootstrapped_scores = []
bootstrapped_tprs = []
fpr_grid = np.linspace(0, 1, 201)

y_true_arr = np.array(y_true)
y_scores_arr = np.array(y_scores)
y_pred_arr = (y_scores_arr >= 0.5).astype(int)

for _ in range(n_bootstraps):
    indices = rng.randint(0, len(y_scores_arr), len(y_scores_arr))
    if len(np.unique(y_true_arr[indices])) < 2:
        continue
    fpr_b, tpr_b, _ = roc_curve(y_true_arr[indices], y_scores_arr[indices])
    score = auc(fpr_b, tpr_b)
    bootstrapped_scores.append(score)

    interp_tpr = np.interp(fpr_grid, fpr_b, tpr_b)
    interp_tpr[0] = 0.0
    interp_tpr[-1] = 1.0
    bootstrapped_tprs.append(interp_tpr)

ci_lower = np.percentile(bootstrapped_scores, 2.5) if bootstrapped_scores else float("nan")
ci_upper = np.percentile(bootstrapped_scores, 97.5) if bootstrapped_scores else float("nan")

tpr_ci_lower = np.percentile(bootstrapped_tprs, 2.5, axis=0) if bootstrapped_tprs else np.zeros_like(fpr_grid)
tpr_ci_upper = np.percentile(bootstrapped_tprs, 97.5, axis=0) if bootstrapped_tprs else np.zeros_like(fpr_grid)

accuracy = accuracy_score(y_true_arr, y_pred_arr)
precision = precision_score(y_true_arr, y_pred_arr, zero_division=0)
recall = recall_score(y_true_arr, y_pred_arr, zero_division=0)
f1 = f1_score(y_true_arr, y_pred_arr, zero_division=0)
tn, fp, fn, tp = confusion_matrix(y_true_arr, y_pred_arr).ravel()
specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

roc_color = "#ff5a36"
roc_fill_color = "#ffd9cf"

fig, ax = plt.subplots(figsize=(7.4, 7.0))
ax.set_facecolor("#fffdfb")

ax.fill_between(
    fpr_grid * 100,
    tpr_ci_lower * 100,
    tpr_ci_upper * 100,
    step="post",
    color=roc_fill_color,
    alpha=0.45,
    linewidth=0,
)
ax.step(fpr * 100, tpr * 100, where="post", color=roc_color, linewidth=2.8)
ax.plot([0, 100], [0, 100], linestyle="--", color="#cfcfcf", linewidth=1.2)

ax.set_xlim(-5, 105)
ax.set_ylim(-5, 105)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xlabel("False Positive Percentage", fontsize=21)
ax.set_ylabel("True Positive Percentage", fontsize=21)
ax.set_title(f"CRC {DATASET_NAME} (bestK)", fontsize=20, pad=18)

for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
    tick_label.set_fontsize(16)

ax.grid(False)
for spine in ax.spines.values():
    spine.set_linewidth(1.2)
    spine.set_color("#bdbdbd")

metric_text = "\n".join(
    [
        f"AUROC: {roc_auc * 100:.2f}%",
        f"95% CI: {ci_lower * 100:.2f}-{ci_upper * 100:.2f}%",
        f"Accuracy: {accuracy * 100:.2f}%",
        f"Precision: {precision * 100:.2f}%",
        f"Recall: {recall * 100:.2f}%",
        f"Specificity: {specificity * 100:.2f}%",
        f"F1 score: {f1 * 100:.2f}%",
    ]
)

ax.text(
    0.98,
    0.07,
    metric_text,
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=13.5,
    color=roc_color,
    linespacing=1.35,
    bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="none", alpha=0.8),
)

plt.tight_layout()

roc_pdf = f"{OUT_PREFIX}_roc_curve.pdf"
plt.savefig(roc_pdf, bbox_inches="tight")
plt.close()


# =========================
# 6) Calibration curve
# =========================
prob_true, prob_pred = calibration_curve(y_true, y_scores, n_bins=10)
plt.figure(figsize=(6, 6))
plt.plot(prob_pred, prob_true, marker="o", label="Calibration curve")
plt.plot([0, 1], [0, 1], linestyle="--", label="Perfectly calibrated")
plt.xlabel("Mean predicted probability")
plt.ylabel("Fraction of positives")
plt.title("Calibration Curve (bestK)")
plt.legend(loc="lower right")
plt.tight_layout()

cal_pdf = f"{OUT_PREFIX}_calibration_curve.pdf"
plt.savefig(cal_pdf)
plt.close()


# =========================
# 7) DCA net benefit
# =========================
def calculate_net_benefit_model(thresh_group, y_pred_score, y_label):
    net_benefit_model = np.array([])
    for thresh in thresh_group:
        y_pred_label = y_pred_score > thresh
        tn, fp, fn, tp = confusion_matrix(y_label, y_pred_label).ravel()
        n = len(y_label)
        net_benefit = (tp / n) - (fp / n) * (thresh / (1 - thresh))
        net_benefit_model = np.append(net_benefit_model, net_benefit)
    return net_benefit_model


def calculate_net_benefit_all(thresh_group, y_label):
    net_benefit_all = np.array([])
    tn, fp, fn, tp = confusion_matrix(y_label, y_label).ravel()
    total = tp + tn
    for thresh in thresh_group:
        net_benefit = (tp / total) - (tn / total) * (thresh / (1 - thresh))
        net_benefit_all = np.append(net_benefit_all, net_benefit)
    return net_benefit_all


thresh_group = np.arange(0, 1, 0.01)
net_benefit_model = calculate_net_benefit_model(thresh_group, np.array(y_scores), np.array(y_true))
net_benefit_all = calculate_net_benefit_all(thresh_group, np.array(y_true))

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(thresh_group, net_benefit_model, color="crimson", label="Model")
ax.plot(thresh_group, net_benefit_all, color="black", label="Treat all")
ax.plot((0, 1), (0, 0), color="black", linestyle=":", label="Treat none")

y2 = np.maximum(net_benefit_all, 0)
y1 = np.maximum(net_benefit_model, y2)
ax.fill_between(thresh_group, y1, y2, color="crimson", alpha=0.2)

ax.set_xlim(0, 1)
ax.set_ylim(min(net_benefit_model.max() - 0.15, -0.18), net_benefit_model.max() + 0.15)
ax.set_xlabel("Threshold Probability", fontsize=15)
ax.set_ylabel("Net Benefit", fontsize=15)
ax.grid("major")
ax.spines["right"].set_color((0.8, 0.8, 0.8))
ax.spines["top"].set_color((0.8, 0.8, 0.8))
ax.legend(loc="upper right")
plt.tight_layout()

dca_pdf = f"{OUT_PREFIX}_dca_curve.pdf"
plt.savefig(dca_pdf)
plt.close()


# =========================
# 8) Per-sample predictions + metrics
# =========================
sample_names = abundance_rel.index.tolist()

results_df = pd.DataFrame(
    {
        "Sample": sample_names,
        "True Label": y_true,
        "Predicted Probability": y_scores,
    }
)

results_df["Predicted Label"] = (results_df["Predicted Probability"] >= 0.5).astype(int)
results_df["Correct"] = results_df["True Label"] == results_df["Predicted Label"]

metrics_df = pd.DataFrame(
    {
        "Metric": [
            "Accuracy",
            "Precision",
            "Recall",
            "Specificity",
            "F1 Score",
            "AUROC",
            "AUROC 95% CI Lower",
            "AUROC 95% CI Upper",
        ],
        "Score": [accuracy, precision, recall, specificity, f1, roc_auc, ci_lower, ci_upper],
    }
)

metrics_csv = f"{OUT_PREFIX}_metrics.csv"
pred_csv = f"{OUT_PREFIX}_prediction_results.csv"
sample_table_csv = f"{OUT_PREFIX}_sample_table.tsv"

metrics_df.to_csv(metrics_csv, index=False)
results_df.to_csv(pred_csv, index=False)
sample_table.reset_index().to_csv(sample_table_csv, sep="\t", index=False)

print(metrics_df)
print("[INFO] Misclassified samples:")
print(results_df[~results_df["Correct"]])

print(f"[INFO] Saved ROC: {roc_pdf}")
print(f"[INFO] Saved Calibration: {cal_pdf}")
print(f"[INFO] Saved DCA: {dca_pdf}")
print(f"[INFO] Saved metrics: {metrics_csv}")
print(f"[INFO] Saved predictions: {pred_csv}")
print(f"[INFO] Saved sample table: {sample_table_csv}")
