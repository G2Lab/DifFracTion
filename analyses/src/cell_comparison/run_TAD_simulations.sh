#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=./PBS/TAD_simulations_%x.out
#SBATCH --error=./PBS/TAD_simulations_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction
    echo "[INFO] Using conda environment diffraction."

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

while getopts "d:r:l:t:s:" opt; do
  case $opt in
    d) diffraction_type="$OPTARG" ;;
    r) resolution="$OPTARG" ;;
    l) log_fold="$OPTARG" ;;
    t) type_analysis="$OPTARG" ;; #cell_line or tissue
    s) sample="$OPTARG" ;; # which cell line or tissue to analyze
    *) echo "Usage: $0 -d <diffraction_type> -r <resolution> -l <log_fold> -t <type_analysis> -s <cell_line|tissue>"
       exit 1 ;;
  esac
done

if [ "$type_analysis" == "cell_line" ]; then
    extra_args="--cell_line $sample"
elif [ "$type_analysis" == "tissue" ]; then
    extra_args="--tissue $sample"
fi

#submiting from analyses, so cd sends back to the parent directory
cd ..
python -m analyses.src.cell_comparison.TAD_simulations \
  --diffraction_type $diffraction_type \
  --resolution $resolution --log_fold $log_fold \
  --type $type_analysis $extra_args