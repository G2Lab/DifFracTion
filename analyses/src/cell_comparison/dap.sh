#!/bin/bash

usage() {
    echo "Usage: dap.sh --bed_a <file> --bed_b <file> --name_a <name> --name_b <name> --output <dir> [--overlap <fraction>]"
    echo ""
    echo "  --bed_a    BED file for condition A (e.g. GM12878_peaks.bed)"
    echo "  --bed_b    BED file for condition B (e.g. K562_peaks.bed)"
    echo "  --name_a   Name for condition A (e.g. GM12878)"
    echo "  --name_b   Name for condition B (e.g. K562)"
    echo "  --mark     ChIP-seq mark (e.g. H3K27Ac)"
    echo "  --output   Output directory"
    echo "  --overlap  Reciprocal overlap fraction cutoff (default: 0.5)"
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bed_a)   bed_a="$2";   shift ;;
        --bed_b)   bed_b="$2";   shift ;;
        --name_a)  name_a="$2";  shift ;;
        --name_b)  name_b="$2";  shift ;;
        --mark)    mark="$2";    shift ;;
        --output)  output="$2";  shift ;;
        --overlap) overlap="$2"; shift ;;
        *) echo "[ERROR] Unknown parameter: $1"; usage ;;
    esac
    shift
done

if [[ -z "$bed_a" || -z "$bed_b" || -z "$name_a" || -z "$name_b" || -z "$output" || -z "$mark" ]]; then
    echo "[ERROR] Missing required arguments."
    usage
fi

if [[ ! -f "$bed_a" ]]; then echo "[ERROR] BED file not found: $bed_a"; exit 1; fi
if [[ ! -f "$bed_b" ]]; then echo "[ERROR] BED file not found: $bed_b"; exit 1; fi

mkdir -p "$output"

total_a=$(grep -vc "^#\|^track\|^browser" "$bed_a" || true)
total_b=$(grep -vc "^#\|^track\|^browser" "$bed_b" || true)

echo "[INFO] Total peaks in ${name_a}: ${total_a}"
echo "[INFO] Total peaks in ${name_b}: ${total_b}"
echo "[INFO] Reciprocal overlap fraction cutoff: ${overlap}"

# Peaks unique to A (not overlapping B by the given reciprocal fraction)
unique_a="${output}/${name_a}_${mark}_unique.bed"
bedtools intersect -a "$bed_a" -b "$bed_b" -v -f "$overlap" -r | grep -v "chrX" | grep -v "chrY" > "$unique_a"

# Peaks unique to B (not overlapping A)
unique_b="${output}/${name_b}_${mark}_unique.bed"
bedtools intersect -a "$bed_b" -b "$bed_a" -v -f "$overlap" -r | grep -v "chrX" | grep -v "chrY" > "$unique_b"

# Shared peaks (A peaks that overlap B)
shared="${output}/${name_a}_${name_b}_${mark}_shared.bed"
bedtools intersect -a "$bed_a" -b "$bed_b" -u -f "$overlap" -r | grep -v "chrX" | grep -v "chrY" > "$shared"

n_unique_a=$(wc -l < "$unique_a")
n_unique_b=$(wc -l < "$unique_b")
n_shared=$(wc -l < "$shared")

echo "[INFO] Peaks unique to ${name_a}: ${n_unique_a}"
echo "[INFO] Peaks unique to ${name_b}: ${n_unique_b}"
echo "[INFO] Shared peaks: ${n_shared}"
echo "[INFO] Results written to ${output}"