#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --output=./PBS/ds_benchmarking_%x.out
#SBATCH --error=./PBS/ds_benchmarking_%x.err

## This script runs the benchmarking of DifFracTion against other tools (diffHiC, HiCcompare, HiCDCPlus, multiHiCcompare) using the generated input files. It evaluates the performance of each tool based on the generated spike-ins and their neighbors.
## It differs from run becnhmarking because here we are focusing on the performance (FP) of tools when the data is downsampled not when we introduce spike ins.

#mamba env create -f ../environment_benchmarking.yml

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction_benchmarking
    echo "[INFO] Using conda environment diffraction_benchmarking. R 4.5.3"

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

if ! command -v Rscript &> /dev/null; then
    echo "[ERROR] Rscript is not available. Please ensure R is installed and Rscript is in your PATH."
    exit 1
fi

check_packs=false

if $check_packs; then
    missing_pkg=false
    set +H  # disable bash history expansion so ! in R expressions works
    # Check required R packages are installed in their lib dirs
    declare -A R_PACKAGES=(
        ["diffHic"]="diffHic_lib"
        ["HiCcompare"]="HiCcompare_lib"
        ["HiCDCPlus"]="HiCDCPlus_lib"
        ["multiHiCcompare"]="multiHiCcompare_lib"
    )
    for pkg in "${!R_PACKAGES[@]}"; do
        if ! Rscript -e "if (!requireNamespace('$pkg', quietly=TRUE)) quit(status=1)" &> /dev/null; then
            echo "[ERROR] R package '$pkg' not found. Make sure to activate the correct conda environment ."
            missing_pkg=true
        else
            echo "[INFO] R package '$pkg' is available."
        fi
    done
    if $missing_pkg; then exit 1; fi
fi

usage() {
    echo ""
    echo "Usage: $0 -f <HiC_file> -c <chromosome> -r <resolution> -d <downsample_factor> [options]"
    echo ""
    echo "Options:"
    echo "  -f,  --hic              Input Hi-C file (.hic format)"
    echo "  -c,  --chrom            Chromosome to process (e.g., 1)"
    echo "  -r,  --resolution       Resolution (bin size) in bp"
    echo "  -p,  --pval             P-value threshold for significance."
    echo "  -d,  --downsample_factor  Downsample factor. Default: 1.0"
    echo "  -s,  --tool             Run only a specific tool: diffraction, diffHic, HiCcompare, HiCDCPlus, multiHiCcompare. Default: all."
    echo "  -t,  --type_norm        DifFracTion norm: alpha, iterative, or both. Only used with --tool diffraction."
    echo ""
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--hic) hic="$2"; shift 2 ;;
        -c|--chrom) chrom="$2"; shift 2 ;;
        -r|--resolution) resolution="$2"; shift 2 ;;
        -p|--pval) p_val_threshold="$2"; shift 2 ;;
        -d|--downsample_factor) downsample_factor="$2"; shift 2 ;;
        -s|--tool) tool="$2"; shift 2 ;;
        -t|--type_norm) type_norm="$2"; shift 2 ;;
        *) echo "[DifFracTion] [ERROR] Unknown parameter: $1"; usage ;;
    esac
done

if [[ -z "$hic" || -z "$chrom" || -z "$resolution" || -z "$p_val_threshold" ]]; then
    echo "[DifFracTion] [ERROR] Missing required parameters."
    usage
fi

if [[ -z "$downsample_factor" ]]; then
    downsample_factor=1.0
fi

echo "[DifFracTion] [INFO] Starting benchmarking with the following parameters:"
echo "  Hi-C file: $hic"
echo "  Chromosome: $chrom"
echo "  Resolution: $resolution"
echo "  Downsample factor: $downsample_factor"
echo "  P-value threshold: $p_val_threshold"
echo "  Tool: ${tool:-all}"
echo "  Directory: ../results_downsample/{tool_name}/${chrom}_${resolution}_${downsample_factor}/"
echo "----------------------------------------"

# Returns true if this tool should run (no --tool flag = run all) (aka we did not specify a tool)
run_tool() { [[ -z "$tool" || "$tool" == "$1" ]]; }

echo "[DifFracTion] [INFO] Generating input files for benchmarking..."
cd ..
python -m benchmarking.generate_inputs_downsample --hic $hic --chrom $chrom --resolution $resolution --downsample_factor $downsample_factor ${tool:+--tool $tool}
cd benchmarking
echo "[DifFracTion] [INFO] Input files generation completed."
echo "----------------------------------------"

# DifFracTion
if run_tool diffraction; then
    echo "[DifFracTion] [INFO] Evaluating DifFracTion (norm: ${type_norm:-both})..."
    echo "----------------------------------------"
    diffraction_dir=$(realpath "../results_downsample/DifFracTion/${chrom}_${resolution}_${downsample_factor}/input_files/")
    matrix_A=$(find ${diffraction_dir} -name "*downsampled*.npz" | head -1)
    matrix_B=$(find ${diffraction_dir} -name "*kb.npz" | head -1)
    echo "  matrixA (downsampled): ${matrix_A}"
    echo "  matrixB (original):    ${matrix_B}"
    echo "  output_dir: $(dirname ${diffraction_dir})/performance/"
    if [[ -z "$type_norm" || "$type_norm" == "alpha" || "$type_norm" == "both" ]]; then
        echo "  [norm: alpha]"
        sbatch run_DifFracTion_ds.sh --matrixA ${matrix_A} --matrixB ${matrix_B} \
            --chrom $chrom --resolution $resolution --pval $p_val_threshold \
            --type_norm alpha --adjusted_pvalues_method distance \
            --output_dir $(dirname ${diffraction_dir})/performance/
    fi
    if [[ -z "$type_norm" || "$type_norm" == "iterative" || "$type_norm" == "both" ]]; then
        echo "  [norm: iterative]"
        sbatch run_DifFracTion_ds.sh --matrixA ${matrix_A} --matrixB ${matrix_B} \
            --chrom $chrom --resolution $resolution --pval $p_val_threshold \
            --type_norm iterative --adjusted_pvalues_method distance \
            --output_dir $(dirname ${diffraction_dir})/performance/
    fi
    echo "----------------------------------------"
fi

# diffHic
if run_tool diffHic; then
    echo "[DifFracTion] [INFO] Evaluating diffHiC..."
    echo "----------------------------------------"
    diffHiC_dir=$(realpath "../results_downsample/diffHic/${chrom}_${resolution}_${downsample_factor}/input_files/")
    echo "  table: ${diffHiC_dir}/diffHic_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table"
    Rscript run_diffHiC_ds.R ${diffHiC_dir}/diffHic_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table ${p_val_threshold}
    echo "----------------------------------------"
fi

# HiCcompare
if run_tool HiCcompare; then
    echo "[DifFracTion] [INFO] Evaluating HiCcompare..."
    echo "----------------------------------------"
    HiCcompare_dir=$(realpath "../results_downsample/HiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/")
    echo "  table: ${HiCcompare_dir}/HiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table"
    Rscript run_HiCcompare_ds.R ${HiCcompare_dir}/HiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table ${p_val_threshold}
    echo "----------------------------------------"
fi

# HiCDCPlus
if run_tool HiCDCPlus; then
    echo "[DifFracTion] [INFO] Evaluating HiCDCPlus..."
    echo "----------------------------------------"
    HiCDCPlus_dir=$(realpath "../results_downsample/HiCDCPlus/${chrom}_${resolution}_${downsample_factor}/input_files/")
    echo "  A1: ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A1.table"
    Rscript run_HiCDCPlus_ds.R \
        ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A1.table \
        ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A2.table \
        ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B1.table \
        ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B2.table \
        ${p_val_threshold}
    echo "----------------------------------------"
fi

# multiHiCcompare
if run_tool multiHiCcompare; then
    echo "[DifFracTion] [INFO] Evaluating multiHiCcompare..."
    echo "----------------------------------------"
    multiHiCcompare_dir=$(realpath "../results_downsample/multiHiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/")
    echo "  A1: ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_A1.table"
    Rscript runmultiHiCcompare_ds.R \
        ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_A1.table \
        ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_A2.table \
        ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_B1.table \
        ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_B2.table \
        ${p_val_threshold}
    echo "----------------------------------------"
fi