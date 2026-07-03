import hicstraw
import sys,os,shutil
import importlib
import numpy as np
import pandas as pd
import argparse
import scipy.sparse

import diffraction_edits
from diffraction_edits import src as DifFracTion
from diffraction_edits import utils as DifFracTion_utils
from diffraction_edits import spikein as DifFracTion_spikein
from diffraction_edits import benchmarking as DifFracTion_bench

'''REAL DATA BENCHMARKING - SPIKE-INS
no downsampling
STARTED : 04 May 2026
STATUS: FULLY FUNCTIONAL'''

'''This script generates input files for benchmarking differential Hi-C tools by introducing spike-ins into Hi-C data.
It processes a specified chromosome at a given resolution and spike-in fold change factor (k).
The script creates altered Hi-C matrices with spike-ins and saves the necessary input files for each tool'''


def parse_args():

     parser = argparse.ArgumentParser(description="Plot the distance distribution of Hi-C contacts")
     parser.add_argument("-f", "--hic", required=True, help="Input Hi-C file (hic format)")
     parser.add_argument("-c", "--chrom", required=True, help="Chromosome to process (e.g., 1)")
     parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution (bin size) in bp")
     parser.add_argument("-k","--k_spike", type=int, default=5, required=True, help="Spike-in fold change factor k. Default: 5")
     return parser.parse_args()

args = parse_args()
main_path=DifFracTion_utils.main_path
print("[INFO] Parameters:")
print(f"  Hi-C file: {args.hic}")
print(f"  Chromosome: {args.chrom}")
print(f"  Resolution: {args.resolution}")
print(f"  K-spike: {args.k_spike}")
print("Main path of the project: ", main_path)

# EXAMPLE USAGE:
# Run in main path of the project, i.e. DifFracTion/
# python -m benchmarking.generate_inputs   -f test_data/GM12878-HRC.hic   -c 10   -r 50000   -k 5

def main():
     # Fixed parameters
     neighbor_degrees = 2
     tool_names = ['DifFracTion', 'HiCcompare', 'multiHiCcompare', 'diffHic', 'HiCDCPlus']
     # we are on benchmarking/
     dictionary_files = {}
     matrices_dictionary = {}
     # Create output directories
     for tool_name in tool_names:
          output_path = f'{main_path}/results_spikeins/{tool_name}/'
          os.makedirs(output_path, exist_ok=True)

          # Directory Inputs for R script generation
          os.makedirs(f"{output_path}/{args.chrom}_{args.resolution}_{args.k_spike}/", exist_ok=True)
          input_files_output_dir = f"{output_path}/{args.chrom}_{args.resolution}_{args.k_spike}/input_files/"
          os.makedirs(input_files_output_dir, exist_ok=True)

          coordinate_file =  f"{input_files_output_dir}spikein_and_neighbors_coordinates_k{args.k_spike}.txt"
          
          dictionary_files[tool_name] = coordinate_file

     # Generate spike-ins once — all tools use the same altered matrix
     # We use DifFractTion because is the only tool that requires both matrices so generate_and_save_spike_ins 
     # saves the original matrix on DifFraction directory
     first_tool = tool_names[0]

     coordinate_file_ref = dictionary_files[first_tool]
     matrix_file_output_ref = os.path.dirname(coordinate_file_ref)
     altered_matrix, raw_matrix, _, name_matrix = DifFracTion_bench.generate_and_save_spike_ins(
          args.hic, args.chrom, args.resolution, args.k_spike,
          neighbor_degrees, matrix_file_output_ref, coordinate_file_ref
     )


     for tool_name in tool_names:

          coordinate_file = dictionary_files[tool_name]
          matrix_file_output = os.path.dirname(coordinate_file)
          matrices_dictionary[tool_name] = name_matrix # assigns DifFractions matrix file, all tools use the same

          if tool_name != first_tool: #to avoid the shutil copy for the first tool since we already generated the spike-ins and saved the coordinate file
               shutil.copy(coordinate_file_ref, coordinate_file)
   

          if tool_name == 'DifFracTion':
              name_altered=os.path.basename(name_matrix).split(".npz")[0]
              altered_name_matrix = f"{matrix_file_output}/{name_altered}_spike_{args.k_spike}.npz"
              scipy.sparse.save_npz(altered_name_matrix, scipy.sparse.csr_matrix(altered_matrix))
              print("=" * 50)
              print(f"[DifFracTion] Input Altered: {altered_name_matrix}")
              print(f"[DifFracTion] Input Original: {name_matrix}")
              print(f"[DifFracTion] Spike-in coordinates: {coordinate_file_ref}")
              print("=" * 50)
          if tool_name == 'HiCcompare':
               save_path = os.path.dirname(coordinate_file)
               save_path = f"{save_path}/HiCcompare_input_chr{args.chrom}_res{args.resolution}_k{args.k_spike}.table"
               
               hiccompare,hiccompare_out= DifFracTion_bench.generate_HiCcompare_input(altered_matrix, raw_matrix, 
                                   args.chrom, args.resolution, args.k_spike,
                                   os.path.dirname(coordinate_file))
               
               print(f"[HiCcompare] Input Table (altered and original): {save_path}")
               print(f"[HiCcompare] Spike-in coordinates: {coordinate_file}")
               print("=" * 50)
          if tool_name == 'diffHic':
               diffHic_out,diffHic_path = DifFracTion_bench.generate_diffHic_input(altered_matrix, raw_matrix,
                                        args.chrom, args.resolution, args.k_spike,
                                        os.path.dirname(coordinate_file))
               print(f"[diffHic] Input Table (altered and original): {diffHic_path}")
               print(f"[diffHic] Spike-in coordinates: {coordinate_file}")
               print("=" * 50)
          if tool_name == 'multiHiCcompare':
               multiHiCcompare_outA1, multiHiCcompare_outA2, multiHiCcompare_outB1, multiHiCcompare_outB2 = DifFracTion_bench.generate_multiHiCcompare_input(altered_matrix, raw_matrix,
                                        args.chrom, args.resolution, args.k_spike,
                                        os.path.dirname(coordinate_file))
               print(f"[multiHiCcompare] Output A1: {multiHiCcompare_outA1}")
               print(f"[multiHiCcompare] Output A2: {multiHiCcompare_outA2}")
               print(f"[multiHiCcompare] Output B1: {multiHiCcompare_outB1}")
               print(f"[multiHiCcompare] Output B2: {multiHiCcompare_outB2}")
               print(f"[multiHiCcompare] Spike-in coordinates: {coordinate_file}")
               print("=" * 50)
          if tool_name == 'HiCDCPlus':

               HiCDCPlus_outA,HiCDCPlus_outA2,HiCDCPlus_outB,HiCDCPlus_outB2 = DifFracTion_bench.generate_HiCDCPlus_input(altered_matrix, raw_matrix,
                                        args.chrom, args.resolution, args.k_spike,
                                        os.path.dirname(coordinate_file))
               print(f"[HiCDCPlus] Output A: {HiCDCPlus_outA}")
               print(f"[HiCDCPlus] Output A2: {HiCDCPlus_outA2}")
               print(f"[HiCDCPlus] Output B: {HiCDCPlus_outB}")
               print(f"[HiCDCPlus] Output B2: {HiCDCPlus_outB2}")
               print(f"[HiCDCPlus] Spike-in coordinates: {coordinate_file}")
               print("=" * 50)
          #                              
if __name__ == "__main__":
    main()