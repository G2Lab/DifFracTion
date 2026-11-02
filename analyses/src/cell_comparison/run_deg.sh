#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --output=./PBS/DEG_%x.out
#SBATCH --error=./PBS/DEG_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate diffraction_deg
    echo "[INFO] Using conda environment diffraction_deg."
else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

counts_a=$1   # e.g. ../test_data/RNAseq/GM12878.rawCounts
counts_b=$2   # e.g. ../test_data/RNAseq/K562.rawCounts
name_a=$3     # e.g. GM12878
name_b=$4     # e.g. K562
lfc=${5:-1.0}
padj=${6:-0.05}
type=${7:-"cell_line"}  # default to cell_line if not provided

if [ "$type" == "cell_line" ]; then
    output_dir="./cell_comparison/DEG/"
elif [ "$type" == "tissue" ]; then
    output_dir="./tissue_comparison/DEG/"
else
    echo "[ERROR] Unknown type: ${type}. Must be 'cell_line' or 'tissue'."
    exit 1
fi


mkdir -p ${output_dir}
output="${output_dir}/${name_a}_vs_${name_b}_DEGs.tsv"

# Inspect the format of the input files
# If the file contains a formated header like GeneID	GM12878.rep1	GM12878.rep2 then it is assumed to be in the correct format for DESeq2. If not, it will be reformatted.
if ! head -n 1 ${counts_a} | grep -q "GeneID"; then
    echo "[INFO] Generating pseudo-replicates for ${counts_a} and setting DESeq2 format."
    echo "[INFO] Assuming this format comes from ENCODE..."
    cut -f1,5 ${counts_a} | awk -v name=${name_a} -v seed=42 '
BEGIN {
  OFS = "\t"
  srand(seed)
}
NR==1 { print "GeneID", name".rep1", name".rep2"; next }
{
  n = $2 + 0
  if (n < 1) {
    r1 = 0
  } else {
    u1 = rand(); u2 = rand()
    if (u1 < 1e-12) u1 = 1e-12          # avoid log(0)
    z = sqrt(-2*log(u1)) * cos(6.283185307179586 * u2)   # Box-Muller standard normal
    r1 = int(n/2 + sqrt(n)/2 * z + 0.5)  # normal approx to Binomial(n, 0.5)
    if (r1 < 0) r1 = 0
    if (r1 > n) r1 = n
  }
  r2 = n - r1
  if ($1 ~ "^ENSG") {
    print $1, r1, r2
  }
}' > ${counts_a}.tmp
counts_a=${counts_a}.tmp
fi

if ! head -n 1 ${counts_b} | grep -q "GeneID"; then
    echo "[INFO] Generating pseudo-replicates for ${counts_b} and setting DESeq2 format."
    echo "[INFO] Assuming this format comes from ENCODE..."
    cut -f1,5 ${counts_b} | awk -v name=${name_b} -v seed=42 '
BEGIN {
  OFS = "\t"
  srand(seed)
}
NR==1 { print "GeneID", name".rep1", name".rep2"; next }
{
  n = $2 + 0
  if (n < 1) {
    r1 = 0
  } else {
    u1 = rand(); u2 = rand()
    if (u1 < 1e-12) u1 = 1e-12          # avoid log(0)
    z = sqrt(-2*log(u1)) * cos(6.283185307179586 * u2)   # Box-Muller standard normal
    r1 = int(n/2 + sqrt(n)/2 * z + 0.5)  # normal approx to Binomial(n, 0.5)
    if (r1 < 0) r1 = 0
    if (r1 > n) r1 = n
  }
  r2 = n - r1
  if ($1 ~ "^ENSG") {
    print $1, r1, r2
  }
}' > ${counts_b}.tmp
counts_b=${counts_b}.tmp

fi

Rscript src/cell_comparison/deg.R \
    --counts_a ${counts_a} \
    --counts_b ${counts_b} \
    --name_a ${name_a} \
    --name_b ${name_b} \
    --output ${output} \
    --lfc ${lfc} \
    --padj ${padj}