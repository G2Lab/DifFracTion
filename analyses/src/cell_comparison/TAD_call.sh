#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=90G
#SBATCH --output=./PBS/TAD_%x.out
#SBATCH --error=./PBS/TAD_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction
    echo "[INFO] Using conda environment diffraction."

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi
hic_file=$1
name=$(basename ${hic_file} .hic)
resolution=$2
window_size=$3

type_analysis=${4:-"cell"}


#./ refering where is being run from
if [ "$type_analysis" == "cell" ]; then
    results_path="./cell_comparison/TADs"
    output_path=${results_path}
elif [ "$type_analysis" == "tissue" ]; then
    results_path="./tissue_comparison/TADs"
    output_path=${results_path}
else
    echo "[ERROR] Invalid type_analysis value. Use 'cell' or 'tissue'."
    exit 1
fi

mkdir -p ${results_path}

HIC_FILE="${hic_file}" COOL_OUT="${output_path}/${name}_${resolution}.cool" RESOLUTION="${resolution}" python3 - <<'PYEOF'
import hicstraw, os

hic_file   = os.environ['HIC_FILE']
cool_file  = os.environ['COOL_OUT']
resolution = int(os.environ['RESOLUTION'])

hic = hicstraw.HiCFile(hic_file)

assert resolution in hic.getResolutions(), \
    f"{resolution} not in available resolutions: {hic.getResolutions()}"

chromosomes = [c for c in hic.getChromosomes() if c.name != "All"]

size_file = cool_file.replace('.cool', '.chrom.sizes')
txt_file  = cool_file.replace('.cool', '.txt')

with open(size_file, 'w') as fsize:
    for chrom in chromosomes:
        fsize.write(f"{chrom.name}\t{chrom.length}\n")

with open(txt_file, 'w') as fo:
    for i in range(len(chromosomes)):
        for j in range(i, len(chromosomes)):
            chrom1 = chromosomes[i].name
            chrom2 = chromosomes[j].name
            result = hicstraw.straw('observed', 'NONE', hic_file, chrom1, chrom2, 'BP', resolution)
            for r in result:
                fo.write(f"{chrom1}\t{r.binX}\t{r.binX + resolution}\t"
                         f"{chrom2}\t{r.binY}\t{r.binY + resolution}\t{r.counts}\n")

os.system(f"cooler load -f bg2 {size_file}:{resolution} {txt_file} {cool_file}")
os.remove(txt_file)
os.remove(size_file)
PYEOF

cooler balance ${output_path}/${name}_${resolution}.cool
cooltools insulation ${output_path}/${name}_${resolution}.cool ${window_size} -o ${results_path}/${name}_TADs.tsv
     