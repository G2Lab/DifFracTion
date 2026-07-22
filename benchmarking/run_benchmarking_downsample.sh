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
    echo "Usage: $0 --hic <HiC_file> --chrom <chromosome> --resolution <bp> --downsample_factor <factor> [options]"
    echo ""
    echo "Options:"
    echo "  -h,  --hic         Input Hi-C file (.hic format)"
    echo "  -c,  --chrom       Chromosome to process (e.g., 1)"
    echo "  -r,  --resolution  Resolution (bin size) in bp"
    echo "  -p,  --pval        P-value threshold for significance."
    echo "  -d,  --downsample_factor  Downsample factor. Default: 1.0"

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

        *) echo "[DifFracTion] [ERROR] Unknown parameter: $1"; usage ;;
    esac
done
# Check required parameters
if [[ -z "$hic" || -z "$chrom" || -z "$resolution" || -z "$p_val_threshold" ]]; then
    echo "[DifFracTion] [ERROR] Missing required parameters."
    usage
fi  

# Set default value for downsample_factor if not provided   
if [[ -z "$downsample_factor" ]]; then
    downsample_factor=1.0
fi
echo "[DifFracTion] [INFO] Starting benchmarking with the following parameters:"
echo "  Hi-C file: $hic"
echo "  Chromosome: $chrom"
echo "  Resolution: $resolution"
echo "  Downsample factor: $downsample_factor"
echo "  P-value threshold: $p_val_threshold"
echo "  Input files generation for tools will be stored in:"
echo "  Directory: ../results_downsample/{tool_name}/${chrom}_${resolution}_${downsample_factor}/input_files/"

echo "[DifFracTion] [INFO] Running benchmarking..."
echo "----------------------------------------"
echo "[DifFracTion] [INFO] Generating input files for benchmarking..."
echo "----------------------------------------"


run_all_tools=false
generate_inputs=false

if $generate_inputs; then
    echo "[DifFracTion] [INFO] Generating input files for benchmarking..."
    cd ..
    python -m benchmarking.generate_inputs_downsample --hic $hic --chrom $chrom --resolution $resolution --downsample_factor $downsample_factor
    cd benchmarking
fi

echo "[DifFracTion] [INFO] Input files generation completed."
echo "----------------------------------------"


if $run_all_tools; then
     echo "[DifFracTion] [INFO] Starting evaluation of all tools (DifFraction, diffHiC, HiCcompare, HiCDCPlus, multiHiCcompare) for spike in detection..."
     echo "----------------------------------------"

     # --- DifFracTion ---
     echo "[DifFracTion] [INFO] Evaluating DifFracTion..."
     echo "[DifFracTion] [INFO] ../results_downsample/DifFracTion/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"
     diffraction_dir=$(realpath "../results_downsample/DifFracTion/${chrom}_${resolution}_${downsample_factor}/input_files/")
     matrix_A=$(find ${diffraction_dir} -name "*downsampled*.npz" | head -1)
     matrix_B=$(find ${diffraction_dir} -name "*kb.npz" | head -1)
     echo "[DifFracTion] [INFO] DifFracTion inputs:"
     echo "  matrixA (downsampled): ${matrix_A}"
     echo "  matrixB (original):    ${matrix_B}"
     echo "  chrom: ${chrom}  resolution: ${resolution}  pval: ${p_val_threshold}"
     echo "  output_dir: $(dirname ${diffraction_dir})/performance/"
     echo "  [norm: alpha]"
     echo "  [norm: iterative]"
     sbatch run_DifFracTion_ds.sh --matrixA ${matrix_A}  \
                              --matrixB ${matrix_B} \
                              --chrom $chrom \
                              --resolution $resolution \
                              --pval $p_val_threshold \
                              --type_norm alpha \
                              --adjusted_pvalues_method distance \
                              --output_dir $(dirname ${diffraction_dir})/performance/
     sbatch run_DifFracTion_ds.sh --matrixA ${matrix_A}  \
                              --matrixB ${matrix_B} \
                              --chrom $chrom \
                              --resolution $resolution \
                              --pval $p_val_threshold \
                              --type_norm iterative \
                              --adjusted_pvalues_method distance \
                              --output_dir $(dirname ${diffraction_dir})/performance/
     echo "----------------------------------------"

     # --- diffHiC ---
     echo "[DifFracTion] [INFO] Evaluating diffHiC..."
     echo "[DifFracTion] [INFO] ../results_downsample/diffHiC/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"
     diffHiC_dir=$(realpath "../results_downsample/diffHic/${chrom}_${resolution}_${downsample_factor}/input_files/")
     echo "[DifFracTion] [INFO] diffHiC inputs:"
     echo "  table: ${diffHiC_dir}/diffHic_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table"
     echo "  pval:  ${p_val_threshold}"
     Rscript run_diffHiC_ds.R ${diffHiC_dir}/diffHic_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table ${p_val_threshold}
     echo "----------------------------------------"

     # --- HiCcompare ---
     echo "[DifFracTion] [INFO] Evaluating HiCcompare."
     echo "[DifFracTion] [INFO] ../results_downsample/HiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"
     HiCcompare_dir=$(realpath "../results_downsample/HiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/")
     echo "[DifFracTion] [INFO] HiCcompare inputs:"
     echo "  table: ${HiCcompare_dir}/HiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table"
     echo "  pval:  ${p_val_threshold}"
     Rscript run_HiCcompare_ds.R ${HiCcompare_dir}/HiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table ${p_val_threshold}
     echo "----------------------------------------"

     # --- HiCDCPlus ---
     echo "[DifFracTion] [INFO] Evaluating HiCDCPlus."
     echo "[DifFracTion] [INFO] ../results_downsample/HiCDCPlus/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"
     HiCDCPlus_dir=$(realpath "../results_downsample/HiCDCPlus/${chrom}_${resolution}_${downsample_factor}/input_files/")
     echo "[DifFracTion] [INFO] HiCDCPlus inputs:"
     echo "  A1: ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A1.table"
     echo "  A2: ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A2.table"
     echo "  B1: ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B1.table"
     echo "  B2: ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B2.table"
     echo "  pval: ${p_val_threshold}"
     Rscript run_HiCDCPlus_ds.R ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A1.table \
      ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A2.table \
      ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B1.table \
      ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B2.table \
      ${p_val_threshold}
     echo "----------------------------------------"

     # --- multiHiCcompare ---
     echo "[DifFracTion] [INFO] Evaluating multiHiCcompare."
     echo "[DifFracTion] [INFO] ../results_downsample/multiHiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"
     multiHiCcompare_dir=$(realpath "../results_downsample/multiHiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/")
     echo "[DifFracTion] [INFO] multiHiCcompare inputs:"
     echo "  A1: ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_A1.table"
     echo "  A2: ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_A2.table"
     echo "  B1: ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_B1.table"
     echo "  B2: ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_B2.table"
     echo "  pval: ${p_val_threshold}"
     Rscript runmultiHiCcompare_ds.R ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_A1.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_A2.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_B1.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_IF_B2.table \
          ${p_val_threshold}
     echo "----------------------------------------"
else 
    echo "[DifFracTion] [INFO] Evaluating DifFracTion. Set run_all_tools=true to evaluate all tools."
         # --- DifFracTion ---
     echo "[DifFracTion] [INFO] Evaluating DifFracTion..."
     echo "[DifFracTion] [INFO] ../results_downsample/DifFracTion/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"
     diffraction_dir=$(realpath "../results_downsample/DifFracTion/${chrom}_${resolution}_${downsample_factor}/input_files/")
     matrix_A=$(find ${diffraction_dir} -name "*downsampled*.npz" | head -1)
     matrix_B=$(find ${diffraction_dir} -name "*kb.npz" | head -1)
     echo "[DifFracTion] [INFO] DifFracTion inputs:"
     echo "  matrixA (downsampled): ${matrix_A}"
     echo "  matrixB (original):    ${matrix_B}"
     echo "  chrom: ${chrom}  resolution: ${resolution}  pval: ${p_val_threshold}"
     echo "  output_dir: $(dirname ${diffraction_dir})/performance/"
     echo "  [norm: alpha]"
     echo "  [norm: iterative]"
     sbatch run_DifFracTion_ds.sh --matrixA ${matrix_A}  \
                              --matrixB ${matrix_B} \
                              --chrom $chrom \
                              --resolution $resolution \
                              --pval $p_val_threshold \
                              --type_norm alpha \
                              --adjusted_pvalues_method distance \
                              --output_dir $(dirname ${diffraction_dir})/performance/
     sbatch run_DifFracTion_ds.sh --matrixA ${matrix_A}  \
                              --matrixB ${matrix_B} \
                              --chrom $chrom \
                              --resolution $resolution \
                              --pval $p_val_threshold \
                              --type_norm iterative \
                              --adjusted_pvalues_method distance \
                              --output_dir $(dirname ${diffraction_dir})/performance/
fi
