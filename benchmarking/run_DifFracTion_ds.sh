#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=90G
#SBATCH --output=./PBS/ds_diffraction_%x.out
#SBATCH --error=./PBS/ds_diffraction_%x.err


if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction_benchmarking
    echo "[INFO] Using conda environment diffraction_benchmarking. R 4.5.3"

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

usage() {
    echo ""
    echo "Usage: $0 --matrixA <matrix_A_file> --matrixB <matrix_B_file> --chrom <chromosome> --resolution <bp> --pval <p-value threshold> [options]"
    echo "Options:"
    echo "  -matrixA,  --matrixA         Input matrix A file"
    echo "  -matrixB,  --matrixB         Input matrix B file"
    echo "  -c,  --chrom       Chromosome to process (e.g., 1)"
    echo "  -r,  --resolution  Resolution (bin size) in bp"
    echo "  -p,  --pval        P-value threshold for significance."
    echo "  -n,  --type_norm     Type of normalization to apply alpha or iterative (default: alpha)"
    echo "  -a,  --adjusted_pvalues_method     Method for adjusting p-values (default: distance)"
    echo "  -o,  --output_dir    Directory to save the results (/performance)"
    echo ""
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -matrixA|--matrixA) matrix_A="$2"; shift 2 ;;
        -matrixB|--matrixB) matrix_B="$2"; shift 2 ;;
        -c|--chrom) chrom="$2"; shift 2 ;;
        -r|--resolution) resolution="$2"; shift 2 ;;
        -p|--pval) p_val_threshold="$2"; shift 2 ;;
        -n|--type_norm) type_norm="$2"; shift 2 ;;
        -a|--adjusted_pvalues_method) adjusted_pvalues_method="$2"; shift 2 ;;
        -o|--output_dir) output_dir="$2"; shift 2 ;;

        *) echo "[DifFracTion] [ERROR] Unknown parameter: $1"; usage ;;
    esac
done


cd ..
python -m benchmarking.run_DifFracTion_ds --matrixA $matrix_A \
                                       --matrixB $matrix_B \
                                       --chrom $chrom \
                                       --resolution $resolution \
                                       --pval $p_val_threshold \
                                       --type_norm $type_norm \
                                       --adjusted_pvalues_method $adjusted_pvalues_method \
                                       --output_dir $output_dir

##
# Example usage:
#sbatch run_DifFracTion.sh --matrixA /users/jcc2340/g2lab/projects/DifFracTion/results_spikeins/DifFracTion/22_100000_4/input_files/GM12878_chr22_100kb_spike_4.npz \
#                           --matrixB /users/jcc2340/g2lab/projects/DifFracTion/results_spikeins/DifFracTion/22_100000_4/input_files/GM12878_chr22_100kb.npz \
#                            --chrom 22 \
#                            --resolution 100000 \
#                            --pval 0.01 \
#                            --type_norm alpha \
#                            --output_dir /users/jcc2340/g2lab/projects/DifFracTion/results_spikeins/DifFracTion/22_100000_4/performance/

