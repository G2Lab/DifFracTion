import hicstraw
import sys,os
import importlib
import numpy as np
import pandas as pd
import argparse

import diffraction
from diffraction import src as DifFracTion
from diffraction import utils as DifFracTion_utils
from diffraction import spikein as DifFracTion_spikein
from diffraction import benchmarking as DifFracTion_bench

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

def main():
     # Fixed parameters
     neighbor_degrees = 2
     tool_names = ['HiCcompare', 'multiHiCcompare', 'diffHic', 'HiCDCPlus']
     # we are on benchmarking/
     dictionary_files = {}
     # Create output directories
     for tool_name in tool_names:
          output_path = f'./results/{tool_name}/'
          os.makedirs(output_path, exist_ok=True)

          # Directory Inputs for R script generation
          os.makedirs(f"{output_path}/{args.chrom}_{args.resolution}_{args.k_spike}/", exist_ok=True)
          input_files_output_dir = f"{output_path}/{args.chrom}_{args.resolution}_{args.k_spike}/input_files/"
          os.makedirs(input_files_output_dir, exist_ok=True)
          coordinate_file =  f"{input_files_output_dir}spikein_and_neighbors_coordinates_k{args.k_spike}.txt"
          dictionary_files[tool_name] = coordinate_file

     for tool_name in tool_names:
	
          coordinate_file = dictionary_files[tool_name]
          altered_matrix, raw_matrix, spikein_all_lists = DifFracTion_bench.generate_and_save_spike_ins(args.hic,args.chrom,args.resolution,args.k_spike,
                                                                 neighbor_degrees,coordinate_file)
               
          if tool_name == 'HiCcompare':
               hiccompare,hiccompare_out= DifFracTion_bench.generate_HiCcompare_input(altered_matrix, raw_matrix, 
                              args.chrom, args.resolution, args.k_spike,
                              os.path.dirname(coordinate_file))
               
          if tool_name == 'diffHic':
               diffHic_out = DifFracTion_bench.generate_diffHic_input(altered_matrix, raw_matrix,
                                        args.chrom, args.resolution, args.k_spike,
                                        os.path.dirname(coordinate_file))
               
          if tool_name == 'multiHiCcompare':
               multiHiCcompare_outA1, multiHiCcompare_outA2, multiHiCcompare_outB1, multiHiCcompare_outB2 = DifFracTion_bench.generate_multiHiCcompare_input(altered_matrix, raw_matrix,
                                        args.chrom, args.resolution, args.k_spike,
                                        os.path.dirname(coordinate_file))
          if tool_name == 'HiCDCPlus':
               HiCDCPlus_outA,HiCDCPlus_outA2,HiCDCPlus_outB,HiCDCPlus_outB2 = DifFracTion_bench.generate_HiCDCPlus_input(altered_matrix, raw_matrix,
                                        args.chrom, args.resolution, args.k_spike,
                                        os.path.dirname(coordinate_file))

                                        
if __name__ == "__main__":
    main()