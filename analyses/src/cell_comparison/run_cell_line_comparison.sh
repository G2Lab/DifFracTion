#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=95G
#SBATCH --output=./PBS/cell_line_comparison_%x.out
#SBATCH --error=./PBS/cell_line_comparison_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction
    echo "[INFO] Using conda environment diffraction."

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

hic_path1=$1
hic_path2=$2
chromosome=$3
analysis_name=$4
resolution=$5
pvalue_cutoff=$6
log2_fold_change_cutoff=$7

cd ..

python -m analyses.src.cell_comparison.cell_line_comparison --hic_path_1 ${hic_path1} --hic_path_2 ${hic_path2} \
        --chromosome ${chromosome} --resolution ${resolution} \
        --analysis_name ${analysis_name} --pvalue_cutoff ${pvalue_cutoff} \
        --log2_fold_change_cutoff ${log2_fold_change_cutoff}

#chr1 α = -1.16, intercept = 21.31
# α = -1.22, intercept = 20.70