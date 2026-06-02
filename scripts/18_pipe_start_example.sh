#!/usr/bin/env bash
#SBATCH --job-name=18.丰度计算FengQ_2015
#SBATCH --nodelist=cn204
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --chdir=/path/to/scratch/tmp_slurm
#SBATCH --output=/path/to/scratch/tmp_slurm/FengQ_2015/slurm_logs/%x.%j.out
#SBATCH --error=/path/to/scratch/tmp_slurm/FengQ_2015/slurm_logs/%x.%j.err

set -Eeuo pipefail
trap 'echo "[ERROR] 命令失败：$BASH_COMMAND (行号 $LINENO)" >&2' ERR

REMOTE_HOST="192.168.2.206"
REMOTE_USER_HOST="user@${REMOTE_HOST}"
REMOTE_SAMPLE_LIST="/path/to/data1/FengQ_2015/SRR_Acc_List_CRC.txt"
PIPELINE_SCRIPT="/path/to/scratch/pipeline/2.pipe.sh"
DITASIC_HOME="/path/to/ditasic"
DITASIC_BIN="${DITASIC_HOME}/bin"
PYTHON_BIN="/path/to/conda/bin/python"
JOB_TAG="${SLURM_JOB_ID:-manual_$(date +%Y%m%d%H%M%S)}"
LOCAL_ROOT="/path/to/scratch/tmp_slurm/FengQ_2015/pipe18/${JOB_TAG}"
LOCAL_SAMPLE_LIST="${LOCAL_ROOT}/SRR_Acc_List.txt"
LOCAL_INPUT_VIEW="${LOCAL_ROOT}/input_view"
LOCAL_BIN="${LOCAL_ROOT}/bin"
RAW_INPUT_ROOT="/path/to/scratch/tmp_slurm/FengQ_2015"
OUTPUT_DIR="/path/to/scratch/tmp_slurm/FengQ_2015/Results"
RSYNC_RSH="ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new"

mkdir -p /path/to/scratch/tmp_slurm/FengQ_2015/slurm_logs "${LOCAL_ROOT}" "${LOCAL_BIN}" "${LOCAL_INPUT_VIEW}" "${OUTPUT_DIR}"
ln -sf "${PYTHON_BIN}" "${LOCAL_BIN}/python"
export PATH="${LOCAL_BIN}:${DITASIC_BIN}:${DITASIC_HOME}:${PATH}"

[ -x "${PIPELINE_SCRIPT}" ] || { echo "[ERROR] 204 管道脚本不存在或不可执行: ${PIPELINE_SCRIPT}" >&2; exit 11; }
[ -x "${DITASIC_BIN}/kallisto" ] || { echo "[ERROR] kallisto 不存在: ${DITASIC_BIN}/kallisto" >&2; exit 12; }
[ -f "${DITASIC_HOME}/ditasic_matrix.py" ] || { echo "[ERROR] ditasic_matrix.py 不存在: ${DITASIC_HOME}/ditasic_matrix.py" >&2; exit 13; }
[ -x "${PYTHON_BIN}" ] || { echo "[ERROR] python 不存在: ${PYTHON_BIN}" >&2; exit 14; }

ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new "${REMOTE_USER_HOST}" "test -r '${REMOTE_SAMPLE_LIST}'"
rsync -a -e "${RSYNC_RSH}" "${REMOTE_USER_HOST}:${REMOTE_SAMPLE_LIST}" "${LOCAL_SAMPLE_LIST}"

while read -r sample; do
    [ -n "$sample" ] || continue
    sample_dir="${RAW_INPUT_ROOT}/${sample}"
    [ -d "$sample_dir" ] || continue

    latest_job=""
    for meta in "$sample_dir"/*/meta.txt; do
        [ -f "$meta" ] || continue
        rc=$(awk -F= '/^rc=/{print $2}' "$meta")
        [ "$rc" = "0" ] || continue
        job=$(awk -F= '/^job_id=/{print $2}' "$meta")
        if [ -z "$latest_job" ] || [ "$job" -gt "$latest_job" ]; then
            latest_job="$job"
        fi
    done

    [ -n "$latest_job" ] || continue
    actual_output="${sample_dir}/${latest_job}/output"
    [ -d "$actual_output" ] || continue
    mkdir -p "${LOCAL_INPUT_VIEW}/${sample}"
    ln -sfn "$actual_output" "${LOCAL_INPUT_VIEW}/${sample}/output"
done < "${LOCAL_SAMPLE_LIST}"

bash "${PIPELINE_SCRIPT}" "${LOCAL_SAMPLE_LIST}" "${LOCAL_INPUT_VIEW}" "${OUTPUT_DIR}"
