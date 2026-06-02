import os

sample_list_file = '/path/to/data1/FengQ_2015/SRR_Acc_List_CRC.txt'
output_job_dir = '34_C_TCG_DiTASiC_jobs/FengQ_2015'

os.makedirs(output_job_dir, exist_ok=True)

with open(sample_list_file, 'r') as f:
    samples = [line.strip() for line in f if line.strip()]

template = r'''#!/usr/bin/env bash

#SBATCH --job-name=__SAMPLE__
#SBATCH --nodelist=cn204
#SBATCH --chdir=/path/to/scratch/tmp_slurm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --output=/path/to/scratch/tmp_slurm/FengQ_2015/slurm_logs/%x.%j.out
#SBATCH --error=/path/to/scratch/tmp_slurm/FengQ_2015/slurm_logs/%x.%j.err

set -Eeuo pipefail
trap 'echo "[ERROR] 命令失败：$BASH_COMMAND (行号 $LINENO)" >&2' ERR
umask 002
ulimit -n 4096

THREADS="${SLURM_CPUS_PER_TASK:-16}"
JOB_TAG="${SLURM_JOB_ID:-manual_$(date +%Y%m%d%H%M%S)}"
JOB_NAME="${SLURM_JOB_NAME:-__SAMPLE__}"
BASE_NAME="__SAMPLE__"

REMOTE_HOST="192.168.2.206"
REMOTE_USER_HOST="user@${REMOTE_HOST}"
REMOTE_READ1="/path/to/data1/FengQ_2015/__SAMPLE__/__SAMPLE___1.fastq.gz"
REMOTE_READ2="/path/to/data1/FengQ_2015/__SAMPLE__/__SAMPLE___2.fastq.gz"
REMOTE_DITASIC_ROOT="/path/to/data2/CRC/CCDC2/DiTASiC"

BASE_TMP="/path/to/scratch/tmp_slurm/FengQ_2015/DiTASiC"
LOG_DIR="/path/to/scratch/tmp_slurm/FengQ_2015/slurm_logs"
LOCAL_SHARED_DITASIC_ROOT="${BASE_TMP}/shared_from_206"
LOCAL_INDEX="${LOCAL_SHARED_DITASIC_ROOT}/kallisto_index"
LOCAL_REMOTE_REF_PATHS="${LOCAL_SHARED_DITASIC_ROOT}/ref_paths.txt"
LOCAL_REF_PATHS="${LOCAL_SHARED_DITASIC_ROOT}/ref_paths.local.txt"
LOCAL_SIMILARITY="${LOCAL_SHARED_DITASIC_ROOT}/similarity_matrix.npy"
LOCAL_REF_CACHE_ROOT="${LOCAL_SHARED_DITASIC_ROOT}/reference_cache"
LOCAL_REF_FILES_FROM="${LOCAL_SHARED_DITASIC_ROOT}/ref_paths.files_from.txt"
LOCAL_MAPPING_DIR="${LOCAL_SHARED_DITASIC_ROOT}/ditasic_mapping"
LOCAL_ABUNDANCE_DIR="${LOCAL_SHARED_DITASIC_ROOT}/abundance"
LOCAL_SHARED_RUN_DIR="${LOCAL_SHARED_DITASIC_ROOT}/job_runs/${BASE_NAME}/${JOB_TAG}"
SYNC_LOCK_FILE="${BASE_TMP}/.shared_from_206.lock"

LOCAL_ROOT="${BASE_TMP}/${BASE_NAME}/${JOB_TAG}"
LOCAL_INPUT_DIR="${LOCAL_ROOT}/input"
LOCAL_WORK_DIR="${LOCAL_ROOT}/work"
LOCAL_KNEAD_DIR="${LOCAL_WORK_DIR}/non_human_reads"
LOCAL_TMP_DIR="${LOCAL_WORK_DIR}/ditasic_mapping_tmp"
LOCAL_RUNLOG="${LOCAL_ROOT}/run.log"
LOCAL_META="${LOCAL_ROOT}/meta.txt"
SLURM_STDOUT="${LOG_DIR}/${JOB_NAME}.${JOB_TAG}.out"
SLURM_STDERR="${LOG_DIR}/${JOB_NAME}.${JOB_TAG}.err"

CONDA_BIN="/path/to/conda/bin/conda"
CRC_ENV_PREFIX="/path/to/conda/envs/CRC"
CRC_PYTHON="${CRC_ENV_PREFIX}/bin/python"
KNEADDATA_BIN="${CRC_ENV_PREFIX}/bin/kneaddata"
KNEADDATA_DB="/path/to/scratch/tools/bowtie2-2.5.4-linux-x86_64/hg38_index"
TRIMMOMATIC_DIR="${CRC_ENV_PREFIX}/share/trimmomatic"
DITASIC_BIN_DIR="/path/to/ditasic"
KALLISTO_BIN_DIR="${DITASIC_BIN_DIR}/bin"
KALLISTO_BIN="${KALLISTO_BIN_DIR}/kallisto"
DITASIC_BIN="${DITASIC_BIN_DIR}/ditasic"
DITASIC_MAPPING_PY="${DITASIC_BIN_DIR}/ditasic_mapping.py"
RSYNC_RSH="ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new"

mkdir -p \
  "${LOG_DIR}" \
  "${LOCAL_REF_CACHE_ROOT}" \
  "${LOCAL_MAPPING_DIR}" \
  "${LOCAL_ABUNDANCE_DIR}" \
  "${LOCAL_SHARED_RUN_DIR}" \
  "${LOCAL_INPUT_DIR}" \
  "${LOCAL_KNEAD_DIR}" \
  "${LOCAL_TMP_DIR}"
exec > >(tee -a "${LOCAL_RUNLOG}") 2>&1

write_meta() {
  local rc="$1"
  cat > "${LOCAL_META}" <<META
job_id=${SLURM_JOB_ID:-manual}
job_name=${JOB_NAME}
node=$(hostname)
rc=${rc}
remote_host=${REMOTE_HOST}
remote_read1=${REMOTE_READ1}
remote_read2=${REMOTE_READ2}
remote_ditasic_root=${REMOTE_DITASIC_ROOT}
local_ditasic_root=${LOCAL_SHARED_DITASIC_ROOT}
local_ref_cache_root=${LOCAL_REF_CACHE_ROOT}
local_ref_paths=${LOCAL_REF_PATHS}
local_abundance_dir=${LOCAL_ABUNDANCE_DIR}
local_root=${LOCAL_ROOT}
timestamp=$(date -Is)
META
}

cleanup() {
  local rc=$?
  set +e
  write_meta "$rc"
  cp -f "${LOCAL_RUNLOG}" "${LOCAL_SHARED_RUN_DIR}/run.log" 2>/dev/null || true
  cp -f "${LOCAL_META}" "${LOCAL_SHARED_RUN_DIR}/meta.txt" 2>/dev/null || true
  if [ "$rc" -eq 0 ]; then
    echo "[INFO] 任务完成"
  else
    echo "[INFO] 任务失败，保留 204 本地现场以便排查"
  fi
  echo "[INFO] 本地工作目录: ${LOCAL_ROOT}"
  echo "[INFO] 本地 DiTASiC 目录: ${LOCAL_SHARED_DITASIC_ROOT}"
  echo "[INFO] Slurm 日志: ${SLURM_STDOUT} ${SLURM_STDERR}"
  exit "$rc"
}
trap cleanup EXIT

command -v ssh >/dev/null || { echo "[ERROR] ssh 不存在" >&2; exit 10; }
command -v rsync >/dev/null || { echo "[ERROR] rsync 不存在" >&2; exit 11; }
command -v flock >/dev/null || { echo "[ERROR] flock 不存在" >&2; exit 12; }
[ -x "${CONDA_BIN}" ] || { echo "[ERROR] conda 不存在或不可执行: ${CONDA_BIN}" >&2; exit 13; }
[ -x "${CRC_PYTHON}" ] || { echo "[ERROR] CRC 环境 python 不存在: ${CRC_PYTHON}" >&2; exit 14; }
[ -x "${KNEADDATA_BIN}" ] || { echo "[ERROR] kneaddata 不存在: ${KNEADDATA_BIN}" >&2; exit 15; }
[ -d "${KNEADDATA_DB}" ] || { echo "[ERROR] KneadData 参考库目录不存在: ${KNEADDATA_DB}" >&2; exit 16; }
[ -d "${TRIMMOMATIC_DIR}" ] || { echo "[ERROR] Trimmomatic 目录不存在: ${TRIMMOMATIC_DIR}" >&2; exit 17; }
[ -x "${KALLISTO_BIN}" ] || { echo "[ERROR] kallisto 不存在: ${KALLISTO_BIN}" >&2; exit 18; }
[ -x "${DITASIC_BIN}" ] || { echo "[ERROR] ditasic 不存在: ${DITASIC_BIN}" >&2; exit 19; }
[ -r "${DITASIC_MAPPING_PY}" ] || { echo "[ERROR] ditasic_mapping.py 不存在: ${DITASIC_MAPPING_PY}" >&2; exit 20; }

printf '[INFO] 节点: %s\n' "$(hostname)"
printf '[INFO] 线程数: %s\n' "${THREADS}"
printf '[INFO] 本地工作目录: %s\n' "${LOCAL_ROOT}"
printf '[INFO] 远端输入: %s\n' "${REMOTE_READ1}"
printf '[INFO] 远端输入: %s\n' "${REMOTE_READ2}"
printf '[INFO] 远端 DiTASiC 根目录: %s\n' "${REMOTE_DITASIC_ROOT}"
printf '[INFO] 本地 DiTASiC 根目录: %s\n' "${LOCAL_SHARED_DITASIC_ROOT}"
printf '[INFO] 本地参考基因组缓存: %s\n' "${LOCAL_REF_CACHE_ROOT}"

ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new \
  "${REMOTE_USER_HOST}" \
  "test -r '${REMOTE_READ1}' && test -r '${REMOTE_READ2}' && test -r '${REMOTE_DITASIC_ROOT}/kallisto_index' && test -r '${REMOTE_DITASIC_ROOT}/ref_paths.txt' && test -r '${REMOTE_DITASIC_ROOT}/similarity_matrix.npy'"

echo "[INFO] 从 206 拉取原始输入到 204 本地临时目录"
rsync -av --partial -e "${RSYNC_RSH}" \
  "${REMOTE_USER_HOST}:${REMOTE_READ1}" \
  "${REMOTE_USER_HOST}:${REMOTE_READ2}" \
  "${LOCAL_INPUT_DIR}/"

echo "[INFO] 从 206 同步 DiTASiC 资源到 204 本地目录"
exec 9>"${SYNC_LOCK_FILE}"
flock 9
rsync -av --partial -e "${RSYNC_RSH}" \
  "${REMOTE_USER_HOST}:${REMOTE_DITASIC_ROOT}/kallisto_index" \
  "${REMOTE_USER_HOST}:${REMOTE_DITASIC_ROOT}/ref_paths.txt" \
  "${REMOTE_USER_HOST}:${REMOTE_DITASIC_ROOT}/similarity_matrix.npy" \
  "${LOCAL_SHARED_DITASIC_ROOT}/"

sed -e 's/\r$//' -e '/^[[:space:]]*$/d' -e 's#^/##' "${LOCAL_REMOTE_REF_PATHS}" > "${LOCAL_REF_FILES_FROM}"
[ -s "${LOCAL_REF_FILES_FROM}" ] || { echo "[ERROR] ref_paths 为空: ${LOCAL_REMOTE_REF_PATHS}" >&2; exit 21; }

echo "[INFO] 从 206 同步参考基因组到 204 本地缓存"
rsync -av --partial --files-from="${LOCAL_REF_FILES_FROM}" -e "${RSYNC_RSH}" \
  "${REMOTE_USER_HOST}:/" \
  "${LOCAL_REF_CACHE_ROOT}/"

awk -v cache_root="${LOCAL_REF_CACHE_ROOT}" '
  {
    sub(/\r$/, "", $0)
    if ($0 == "") next
    ref = $0
    if (substr(ref, 1, 1) == "/") {
      print cache_root ref
    } else {
      print cache_root "/" ref
    }
  }
' "${LOCAL_REMOTE_REF_PATHS}" > "${LOCAL_REF_PATHS}"
flock -u 9

LOCAL_READ1="${LOCAL_INPUT_DIR}/$(basename "${REMOTE_READ1}")"
LOCAL_READ2="${LOCAL_INPUT_DIR}/$(basename "${REMOTE_READ2}")"

[ -r "${LOCAL_INDEX}" ] || { echo "[ERROR] 本地 DiTASiC index 不存在: ${LOCAL_INDEX}" >&2; exit 22; }
[ -r "${LOCAL_REMOTE_REF_PATHS}" ] || { echo "[ERROR] 本地 DiTASiC 原始 ref_paths 不存在: ${LOCAL_REMOTE_REF_PATHS}" >&2; exit 23; }
[ -r "${LOCAL_REF_PATHS}" ] || { echo "[ERROR] 本地改写 ref_paths 不存在: ${LOCAL_REF_PATHS}" >&2; exit 24; }
[ -r "${LOCAL_SIMILARITY}" ] || { echo "[ERROR] 本地 DiTASiC similarity_matrix 不存在: ${LOCAL_SIMILARITY}" >&2; exit 25; }

while IFS= read -r ref; do
  [ -r "${ref}" ] || { echo "[ERROR] 本地参考基因组不存在: ${ref}" >&2; exit 26; }
done < "${LOCAL_REF_PATHS}"

echo "[INFO] 激活 CRC 环境"
eval "$("$CONDA_BIN" shell.bash hook)"
set +u
conda activate "${CRC_ENV_PREFIX}"
set -u
export PATH="${CRC_ENV_PREFIX}/bin:${KALLISTO_BIN_DIR}:${DITASIC_BIN_DIR}:$PATH"

echo "[INFO] 运行 kneaddata"
extra_opts=()
if [ -n "${TRIMMOMATIC_DIR}" ]; then
  extra_opts+=(--trimmomatic "${TRIMMOMATIC_DIR}")
fi

"${KNEADDATA_BIN}" \
  --input1 "${LOCAL_READ1}" \
  --input2 "${LOCAL_READ2}" \
  --output "${LOCAL_KNEAD_DIR}" \
  --output-prefix "${BASE_NAME}" \
  --reference-db "${KNEADDATA_DB}" \
  --threads "${THREADS}" \
  --max-memory 30000M \
  --remove-intermediate-output \
  "${extra_opts[@]}"

paired1=$(find "${LOCAL_KNEAD_DIR}" -maxdepth 1 -type f \( -name "${BASE_NAME}_paired_1.fastq" -o -name "${BASE_NAME}_paired_1.fastq.gz" \) | head -n 1 || true)
paired2=$(find "${LOCAL_KNEAD_DIR}" -maxdepth 1 -type f \( -name "${BASE_NAME}_paired_2.fastq" -o -name "${BASE_NAME}_paired_2.fastq.gz" \) | head -n 1 || true)

[ -s "${paired1}" ] || { echo "[ERROR] KneadData 输出不存在: ${paired1}" >&2; exit 27; }
[ -s "${paired2}" ] || { echo "[ERROR] KneadData 输出不存在: ${paired2}" >&2; exit 28; }
printf '[INFO] KneadData 输出: %s\n' "${paired1}"
printf '[INFO] KneadData 输出: %s\n' "${paired2}"

LOCAL_COMBINED_FQ="${LOCAL_WORK_DIR}/${BASE_NAME}_combined.fq"
if [[ "${paired1}" == *.gz ]]; then
  zcat "${paired1}" > "${LOCAL_COMBINED_FQ}"
else
  cat "${paired1}" > "${LOCAL_COMBINED_FQ}"
fi
if [[ "${paired2}" == *.gz ]]; then
  zcat "${paired2}" >> "${LOCAL_COMBINED_FQ}"
else
  cat "${paired2}" >> "${LOCAL_COMBINED_FQ}"
fi

MAPPED_COUNTS_FILE="${BASE_NAME}_combined_mapped_counts.npy"
TOTAL_COUNTS_FILE="${BASE_NAME}_combined_total.npy"
FINAL_MAPPED_COUNTS="${LOCAL_MAPPING_DIR}/${MAPPED_COUNTS_FILE}"
FINAL_TOTAL_COUNTS="${LOCAL_MAPPING_DIR}/${TOTAL_COUNTS_FILE}"
FINAL_ABUNDANCE="${LOCAL_ABUNDANCE_DIR}/abundance_${BASE_NAME}.txt"

echo "[INFO] 运行 ditasic_mapping.py"
pushd "${LOCAL_WORK_DIR}" >/dev/null
"${CRC_PYTHON}" "${DITASIC_MAPPING_PY}" \
  -i "${LOCAL_INDEX}" \
  -l 100 \
  -t "${LOCAL_TMP_DIR}/tmp_${BASE_NAME}" \
  "${LOCAL_REF_PATHS}" \
  "${LOCAL_COMBINED_FQ}"

[ -f "${MAPPED_COUNTS_FILE}" ] || { echo "[ERROR] 缺少映射结果: ${MAPPED_COUNTS_FILE}" >&2; exit 29; }
[ -f "${TOTAL_COUNTS_FILE}" ] || { echo "[ERROR] 缺少总计结果: ${TOTAL_COUNTS_FILE}" >&2; exit 30; }

mv -f "${MAPPED_COUNTS_FILE}" "${FINAL_MAPPED_COUNTS}"
mv -f "${TOTAL_COUNTS_FILE}" "${FINAL_TOTAL_COUNTS}"

echo "[INFO] 运行 ditasic"
"${DITASIC_BIN}" \
  -r "${LOCAL_REF_PATHS}" \
  -a "${LOCAL_SIMILARITY}" \
  -x "${FINAL_MAPPED_COUNTS}" \
  -n "${FINAL_TOTAL_COUNTS}" \
  -o "${FINAL_ABUNDANCE}"
popd >/dev/null

cat > "${LOCAL_SHARED_RUN_DIR}/stage_meta.txt" <<META
job_id=${SLURM_JOB_ID:-manual}
job_name=${JOB_NAME}
node=$(hostname)
source_host=${REMOTE_HOST}
source_read1=${REMOTE_READ1}
source_read2=${REMOTE_READ2}
source_ditasic_root=${REMOTE_DITASIC_ROOT}
local_ditasic_root=${LOCAL_SHARED_DITASIC_ROOT}
local_ref_cache_root=${LOCAL_REF_CACHE_ROOT}
local_ref_paths=${LOCAL_REF_PATHS}
local_root=${LOCAL_ROOT}
mapped_counts=${FINAL_MAPPED_COUNTS}
total_counts=${FINAL_TOTAL_COUNTS}
abundance=${FINAL_ABUNDANCE}
timestamp=$(date -Is)
META

echo "[DONE] ${BASE_NAME} 完成"
echo "[DONE] abundance: ${FINAL_ABUNDANCE}"
'''

for sample in samples:
    script_path = os.path.join(output_job_dir, f"{sample}.sh")
    script_text = template.replace('__SAMPLE__', sample)
    with open(script_path, 'w') as script:
        script.write(script_text)
    os.chmod(script_path, 0o755)
