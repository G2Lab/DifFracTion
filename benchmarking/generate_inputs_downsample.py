import hicstraw
import sys,os

import numpy as np
import pandas as pd
import argparse
import scipy.sparse


from diffraction_edits import src as DifFracTion
from diffraction_edits import utils as DifFracTion_utils
from diffraction_edits import benchmarking as DifFracTion_bench

'''REAL DATA BENCHMARKING
STARTED : 04 May 2026
STATUS: FULLY FUNCTIONAL AS OF May 2026'''

'''This script generates input files for benchmarking differential Hi-C tools 
by downsampling Hi-C matrices. It is not related to spike ins'''

#It will generate downsampled matrices for a given chromosome and resolution and generate results for the following tools: HiCcompare, multiHiCcompare, diffHic, HiCDCPlus.
#it now incorporates run_synthetic_validation_resolution.py that basically does the same but for DifFraction
def parse_args():

     parser = argparse.ArgumentParser(description="Plot the distance distribution of Hi-C contacts")
     parser.add_argument("-f", "--hic", required=True, help="Input Hi-C file (hic format)")
     parser.add_argument("-c", "--chrom", required=True, help="Chromosome to process (e.g., 1)")
     parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution (bin size) in bp")
     parser.add_argument("-d","--downsample_factor", type=float, default=1.0, required=True, help="Downsample factor. Default: 1.0")
     parser.add_argument("-t", "--tool", choices=['diffraction', 'HiCcompare', 'multiHiCcompare', 'diffHic', 'HiCDCPlus'], default=None, help="Tool to generate input files for. Default: all tools")
     return parser.parse_args()

args = parse_args()
# Project directory
main_path=DifFracTion_utils.main_path

# EXAMPLE USAGE:
# Run in main path of the project, i.e. DifFracTion/
# python -m benchmarking.generate_inputs_downsample   -f test_data/GM12878-HRC.hic   -c 1   -r 50000   -d 0.5
# FULLY FUNCTIONAL

def load_matrix(
     hic_path: str,
     chromosome: str,
     resolution: int,
     matrix_save_path: str
):   
     #mainpath
     resolution_kb = resolution // 1000
     #has to be on mainpath/generated_data
     cell_type=os.path.basename(hic_path).split("-")[0]
     name_matrix=f"{matrix_save_path}/{cell_type}_chr{chromosome}_{resolution_kb}kb.npz"
     if not os.path.exists(f"{name_matrix}"):
          DifFracTion_utils.save_sparse_matrix(hic_path,chromosome,resolution,f"{name_matrix}")

     matrix = DifFracTion_utils.load_dense_matrix(f"{name_matrix}")
     
     return matrix,name_matrix



def main():
     raw_matrix,name_matrix = load_matrix(args.hic, args.chrom, args.resolution, f"{main_path}/generated_data")
     tool_names = ['DifFracTion', 'HiCcompare', 'multiHiCcompare', 'diffHic', 'HiCDCPlus']
     # we are on benchmarking/
     dictionary_paths = {}
     # Create output directories
     for tool_name in tool_names:
          output_path = f'{main_path}/results_downsample/{tool_name}/'
          os.makedirs(output_path, exist_ok=True)

          # Directory Inputs for R script generation
          os.makedirs(f"{output_path}/{args.chrom}_{args.resolution}_{args.downsample_factor}/", exist_ok=True)
          input_files_output_dir = f"{output_path}/{args.chrom}_{args.resolution}_{args.downsample_factor}/input_files/"
          os.makedirs(input_files_output_dir, exist_ok=True)
          dictionary_paths[tool_name] = input_files_output_dir


     
     #Now we have the matrices saved so we just load them
     if args.tool == 'DifFracTion' or args.tool is None:
          #change npz for {ds} and .npz for {ds}_downsampled
          downsampled_diffraction_B_name = dictionary_paths['DifFracTion'] + os.path.splitext(os.path.basename(name_matrix))[0] + f"_downsampled_{args.downsample_factor}.npz"
          downsampled_diffraction_A_name = dictionary_paths['DifFracTion'] + os.path.splitext(os.path.basename(name_matrix))[0] + ".npz"
         
          print(downsampled_diffraction_B_name)
          if not os.path.exists(downsampled_diffraction_B_name):
               downsampled_diffraction_B = DifFracTion_utils.synthetic_datasets(raw_matrix, args.downsample_factor)
               downsampled_diffraction_B = scipy.sparse.csr_matrix(downsampled_diffraction_B)
               scipy.sparse.save_npz(downsampled_diffraction_B_name, downsampled_diffraction_B)
          if not os.path.exists(downsampled_diffraction_A_name):
               downsampled_diffraction_A = scipy.sparse.csr_matrix(raw_matrix)
               scipy.sparse.save_npz(downsampled_diffraction_A_name, downsampled_diffraction_A)
               #When loading do the same matrix = DifFracTion_utils.load_dense_matrix("path/to/output.npz")
     
     if args.tool == 'HiCcompare' or args.tool is None:
          save_path = (
               f"{dictionary_paths['HiCcompare']}/HiCcompare_input_chr{args.chrom}_res{args.resolution}_ds{args.downsample_factor}.table"
          )
          print(save_path)
          if not os.path.exists(save_path):
               hiccompare,hiccompare_out= DifFracTion_bench.generate_HiCcompare_input_ds(raw_matrix, 
                                        args.chrom, args.resolution, args.downsample_factor,
                                        dictionary_paths['HiCcompare'])
          
     if args.tool == 'diffHic' or args.tool is None:
          save_path = (
               f"{dictionary_paths['diffHic']}/diffHic_input_chr{args.chrom}_res{args.resolution}_ds{args.downsample_factor}.table"
          )
          if not os.path.exists(save_path):
               diffHic_out = DifFracTion_bench.generate_diffHic_input_ds(raw_matrix,
                                        args.chrom, args.resolution, args.downsample_factor,
                                        dictionary_paths['diffHic'])
     if args.tool == 'multiHiCcompare' or args.tool is None:
          # All files will have the downsample tag but IF_A1,IF_A2 are full depth
          save_path = (
               f"{dictionary_paths['multiHiCcompare']}/multiHiCcompare_input_chr{args.chrom}_res{args.resolution}_ds{args.downsample_factor}_IF_A1.table"
          )
          if not os.path.exists(save_path):
               multiHiCcompare_outA1, multiHiCcompare_outA2, multiHiCcompare_outB1, multiHiCcompare_outB2 = DifFracTion_bench.generate_multiHiCcompare_input_ds(raw_matrix,
                                   args.chrom, args.resolution, args.downsample_factor,
                                   dictionary_paths['multiHiCcompare'])

     if args.tool == 'HiCDCPlus' or args.tool is None:
          # All files will have the downsample tag but A1,A2 are full depth
          save_path = (
               f"{output_path}/HiCDCPlus_input_chr{args.chrom}_res{args.resolution}_ds{args.downsample_factor}_B2.table"
          )
          if not os.path.exists(save_path):
               HiCDCPlus_outA,HiCDCPlus_outA2,HiCDCPlus_outB,HiCDCPlus_outB2 = DifFracTion_bench.generate_HiCDCPlus_input_ds(raw_matrix,
                                   args.chrom, args.resolution, args.downsample_factor,
                                   dictionary_paths['HiCDCPlus'])

     
                                        
if __name__ == "__main__":
    main()