import shutil
import hicstraw
import os

import numpy as np
import pandas as pd
import argparse
import scipy.sparse


from diffraction_edits import src as DifFracTion
from diffraction_edits import utils as DifFracTion_utils
from diffraction_edits import benchmarking as DifFracTion_bench

'''REAL DATA BENCHMARKING
STARTED : 25 May 2026
STATUS:  On development as of May 2026'''

'''This script generates input files for benchmarking differential Hi-C tools 
by downsampling Hi-C matrices. While also introducing spike ins'''

#It will generate downsampled matrices for a given chromosome and resolution and generate results for the following tools: HiCcompare, multiHiCcompare, diffHic, HiCDCPlus.
#it now incorporates run_synthetic_validation_resolution.py that basically does the same but for DifFraction
def parse_args():

     parser = argparse.ArgumentParser(description="Plot the distance distribution of Hi-C contacts")
     parser.add_argument("-f", "--hic", required=True, help="Input Hi-C file (hic format)")
     parser.add_argument("-c", "--chrom", required=True, help="Chromosome to process (e.g., 1)")
     parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution (bin size) in bp")
     parser.add_argument("-d","--downsample_factor", type=float, default=1.0, required=True, help="Downsample factor. Default: 1.0")
     parser.add_argument("-k","--k_spikes", type=int, default=0, required=False, help="Number of spike ins to introduce. Default: 0")
     return parser.parse_args()

args = parse_args()
# Project directory
main_path=DifFracTion_utils.main_path



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
     neighbor_degrees = 2
     tool_names = ['DifFracTion', 'HiCcompare', 'multiHiCcompare', 'diffHic', 'HiCDCPlus']
     # we are on benchmarking/
     dictionary_paths = {}
     # Create output directories
     for tool_name in tool_names:
          output_path = f'{main_path}/results_downsample_n_spikes/{tool_name}/'
          os.makedirs(output_path, exist_ok=True)

          # Directory Inputs for R script generation
          os.makedirs(f"{output_path}/{args.chrom}_{args.resolution}_{args.downsample_factor}_{args.k_spikes}/", exist_ok=True)
          input_files_output_dir = f"{output_path}/{args.chrom}_{args.resolution}_{args.downsample_factor}_{args.k_spikes}/input_files/"
          os.makedirs(input_files_output_dir, exist_ok=True)
          dictionary_paths[tool_name] = input_files_output_dir

     first_tool = tool_names[0]
     coordinate_file_ref = dictionary_paths[first_tool] + f"spikein_and_neighbors_coordinates_k{args.k_spikes}.txt"
     matrix_file_output_ref = os.path.dirname(coordinate_file_ref)

     altered_matrix, raw_matrix, _, name_matrix = DifFracTion_bench.generate_and_save_spike_ins(
          args.hic, args.chrom, args.resolution, args.k_spikes,
          neighbor_degrees, matrix_file_output_ref, coordinate_file_ref
     )

     downsampled = DifFracTion_utils.synthetic_datasets(raw_matrix, args.downsample_factor)
     diffraction_altered = DifFracTion_utils.synthetic_datasets(altered_matrix, 1)

     for tool_name in tool_names:
          if tool_name != first_tool:
               coordinate_file = dictionary_paths[tool_name] + f"spikein_and_neighbors_coordinates_k{args.k_spikes}.txt"
               shutil.copy(coordinate_file_ref, coordinate_file)

          #Now we have the matrices saved so we just load them
          if tool_name == 'DifFracTion':
               
               #change npz for {ds} and .npz for {ds}_downsampled
               downsampled_diffraction_B_name = dictionary_paths[tool_name] + os.path.splitext(os.path.basename(name_matrix))[0] + f"_downsampled_{args.downsample_factor}.npz"
               downsampled_diffraction_A_name = dictionary_paths[tool_name] + os.path.splitext(os.path.basename(name_matrix))[0] + ".npz"
              
               print(downsampled_diffraction_B_name)
               #B matrix will always be the one that is downsampled
               downsampled_diffraction_B = downsampled
               downsampled_diffraction_B = scipy.sparse.csr_matrix(downsampled_diffraction_B)
               scipy.sparse.save_npz(downsampled_diffraction_B_name, downsampled_diffraction_B)

               #A will always be the original matrix with spike ins if k_spikes > 0, or the original matrix if k_spikes = 0
               downsampled_diffraction_A = scipy.sparse.csr_matrix(diffraction_altered)
               scipy.sparse.save_npz(downsampled_diffraction_A_name, downsampled_diffraction_A)

                    #When loading do the same matrix = DifFracTion_utils.load_dense_matrix("path/to/output.npz")
          
          if tool_name == 'HiCcompare':
               save_path = (
                    f"{dictionary_paths['HiCcompare']}/HiCcompare_input_chr{args.chrom}_res{args.resolution}_k{args.k_spikes}.table"
               )
               print(save_path)
               hiccompare,hiccompare_out= DifFracTion_bench.generate_HiCcompare_input(altered_matrix, downsampled,
                                   args.chrom, args.resolution, args.k_spikes,
                                   os.path.dirname(save_path))
               print(f"[HiCcompare] Input Table (altered and original): {save_path}")
               print(f'[HiCcompare] Spike-in coordinates: {dictionary_paths["HiCcompare"]}spikein_and_neighbors_coordinates_k{args.k_spikes}.txt')               
          
          if tool_name == 'diffHic':
               save_path = (
                    f"{dictionary_paths['diffHic']}/diffHic_input_chr{args.chrom}_res{args.resolution}_k{args.k_spikes}.table"
               )
               # Altered matrix will generate two matrices (replicates ) so for downsampled
               diffHic_out,diffHic_path = DifFracTion_bench.generate_diffHic_input(altered_matrix, downsampled,
                                    args.chrom, args.resolution, args.k_spikes,
                                    os.path.dirname(save_path))
          if tool_name == 'multiHiCcompare':
               # All files will have the downsample tag but IF_A1,IF_A2 are full depth
               save_path = (
                    f"{dictionary_paths['multiHiCcompare']}/multiHiCcompare_input_chr{args.chrom}_res{args.resolution}_k{args.k_spikes}_IF_A1.table"
               )
               multiHiCcompare_outA1, multiHiCcompare_outA2, multiHiCcompare_outB1, multiHiCcompare_outB2 = DifFracTion_bench.generate_multiHiCcompare_input(altered_matrix, downsampled,
                                    args.chrom, args.resolution, args.k_spikes,
                                    os.path.dirname(save_path))
#
          if tool_name == 'HiCDCPlus':
               # All files will have the downsample tag but A1,A2 are full depth
               save_path = (
                    f"{dictionary_paths['HiCDCPlus']}/HiCDCPlus_input_chr{args.chrom}_res{args.resolution}_k{args.k_spikes}_B2.table"
               )
               HiCDCPlus_outA,HiCDCPlus_outA2,HiCDCPlus_outB,HiCDCPlus_outB2 = DifFracTion_bench.generate_HiCDCPlus_input(altered_matrix, downsampled,
                                    args.chrom, args.resolution, args.k_spikes,
                                    os.path.dirname(save_path))
          
                                        
if __name__ == "__main__":
    main()