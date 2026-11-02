#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=./PBS/DAP_%x.out
#SBATCH --error=./PBS/DAP_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate diffraction
    echo "[INFO] Using conda environment diffraction."
else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

bed_a=$1     # ../test_data/ChIPseq/GM12878_peaks.bed
bed_b=$2     # ../test_data/ChIPseq/K562_peaks.bed
name_a=$3    # GM12878
name_b=$4    # K562
mark=$5     # H3K27Ac
overlap=${6:-0.5}

output_dir="./tissue_comparison/DAP/"
mkdir -p ${output_dir}

bash $(dirname $0)/dap.sh \
    --bed_a   ${bed_a} --bed_b   ${bed_b} \
    --name_a  ${name_a} --name_b  ${name_b} --mark ${mark} \
    --output  ${output_dir} --overlap ${overlap}