import hicstraw
import sys,os

import numpy as np
import pandas as pd
import argparse

import diffraction
from diffraction import src as DifFracTion
from diffraction import utils as DifFracTion_utils
from diffraction import spikein as DifFracTion_spikein
from diffraction import benchmarking as DifFracTion_bench


'''This script generates input files for benchmarking differential Hi-C tools by downsampling Hi-C matrices.'''


def parse_args():

     parser = argparse.ArgumentParser(description="Plot the distance distribution of Hi-C contacts")
     parser.add_argument("-f", "--hic", required=True, help="Input Hi-C file (hic format)")
     parser.add_argument("-c", "--chrom", required=True, help="Chromosome to process (e.g., 1)")
     parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution (bin size) in bp")
     parser.add_argument("-d","--downsample_factor", type=float, default=1.0, required=True, help="Downsample factor. Default: 1.0")
     return parser.parse_args()

args = parse_args()

def load_matrix(
    hic_path: str,
    chromosome: str,
    resolution: int,
):
    chromosome_length = DifFracTion_utils.get_chromosome_length(hic_path, chromosome)
    matrixZoom = DifFracTion_utils.prepare_MatrixZoomData(hic_path, chromosome, resolution)
    raw_matrix = DifFracTion_utils.rawMatrix(matrixZoom, chromosome_length)
    return raw_matrix


     
def main():
     raw_matrix = load_matrix(args.hic, args.chrom, args.resolution)
     tool_names = ['HiCcompare', 'multiHiCcompare', 'diffHic', 'HiCDCPlus']
     # we are on benchmarking/
     dictionary_paths = {}
     # Create output directories
     for tool_name in tool_names:
          output_path = f'./results_downsample/{tool_name}/'
          os.makedirs(output_path, exist_ok=True)

          # Directory Inputs for R script generation
          os.makedirs(f"{output_path}/{args.chrom}_{args.resolution}_{args.downsample_factor}/", exist_ok=True)
          input_files_output_dir = f"{output_path}/{args.chrom}_{args.resolution}_{args.downsample_factor}/input_files/"
          os.makedirs(input_files_output_dir, exist_ok=True)
          dictionary_paths[tool_name] = input_files_output_dir


     for tool_name in tool_names:
          if tool_name == 'HiCcompare':
               hiccompare,hiccompare_out= DifFracTion_bench.generate_HiCcompare_input_ds(raw_matrix, 
                                        args.chrom, args.resolution, args.downsample_factor,
                                        dictionary_paths['HiCcompare'])
               
          if tool_name == 'diffHic':
               diffHic_out = DifFracTion_bench.generate_diffHic_input_ds(raw_matrix,
                                        args.chrom, args.resolution, args.downsample_factor,
                                        dictionary_paths['diffHic'])
               
          if tool_name == 'multiHiCcompare':
               multiHiCcompare_outA1, multiHiCcompare_outA2, multiHiCcompare_outB1, multiHiCcompare_outB2 = DifFracTion_bench.generate_multiHiCcompare_input_ds(raw_matrix,
                                        args.chrom, args.resolution, args.downsample_factor,
                                        dictionary_paths['multiHiCcompare'])

          if tool_name == 'HiCDCPlus':
               HiCDCPlus_outA,HiCDCPlus_outA2,HiCDCPlus_outB,HiCDCPlus_outB2 = DifFracTion_bench.generate_HiCDCPlus_input_ds(raw_matrix,
                                        args.chrom, args.resolution, args.downsample_factor,
                                        dictionary_paths['HiCDCPlus'])
	
          
                                        
if __name__ == "__main__":
    main()