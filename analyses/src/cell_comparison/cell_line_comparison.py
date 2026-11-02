import hicstraw
import scipy
from scipy.sparse import coo_matrix
from scipy import ndimage
from scipy.stats import linregress
from statsmodels.stats.multitest import multipletests

import numpy as np
import pandas as pd
import os
from typing import Any,Literal
import argparse

import matplotlib.pyplot as plt
from matplotlib import colors

from diffraction import utils as DifFracTion_utils
from diffraction import src as DifFracTion

def parse_args():
    
    parser = argparse.ArgumentParser(description='Iterative correction of Hi-C matrices based on distance decay.')
    parser.add_argument('--hic_path_1', type=str, required=True, help='Path to the input Hi-C file (hic format).')
    parser.add_argument('--hic_path_2', type=str, required=True, help='Path to the second input Hi-C file (hic format).')
    parser.add_argument('--chromosome', type=str, required=True, help='Chromosome to process (e.g., 1).')
    parser.add_argument('--resolution', type=int, default=10000, help='Resolution of the Hi-C matrices (default: 10000).')
    parser.add_argument('--pvalue_cutoff', type=float, default=0.05, help='P-value cutoff for identifying significant interactions (default: 0.05).')
    parser.add_argument('--log2_fold_change_cutoff', type=float, default=1.0, help='Log2 fold change cutoff for identifying significant interactions (default: 1.0).')
    parser.add_argument('--analysis_name', type=str, default='default', help='Directory to save the corrected matrices (default: corrected_matrices).')
    
    return parser.parse_args()

args= parse_args()
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
     matrix_A, name_matrix_A = load_matrix(args.hic_path_1, args.chromosome, args.resolution, f"{main_path}/generated_data")
     matrix_B, name_matrix_B = load_matrix(args.hic_path_2, args.chromosome, args.resolution, f"{main_path}/generated_data")

     output_path = f'{main_path}/{args.analysis_name}/'
     os.makedirs(output_path, exist_ok=True)

     normA,normB =DifFracTion.alpha_normalization(matrix_A, matrix_B, args.resolution, metric='median', desired_alpha=1.08)
     #save the normalized matrices
     np.savez_compressed(f"{output_path}/chr{args.chromosome}_{args.resolution}_alpha_normalized_matrix_A.npz", normA)
     np.savez_compressed(f"{output_path}/chr{args.chromosome}_{args.resolution}_alpha_normalized_matrix_B.npz", normB)

     i_normA,i_normB = DifFracTion.iterative_normalization(matrix_A, matrix_B, args.resolution, metric='median', weight_factor=-1.08, eta=0.15)
     #save the normalized matrices
     np.savez_compressed(f"{output_path}/chr{args.chromosome}_{args.resolution}_iterative_normalized_matrix_A.npz", i_normA)
     np.savez_compressed(f"{output_path}/chr{args.chromosome}_{args.resolution}_iterative_normalized_matrix_B.npz", i_normB)

     alpha_int,_=DifFracTion.identify_differential_interactions(normA, normB, args.resolution, 
                                                    filter_by='padjusted', pvalue_cutoff=args.pvalue_cutoff
                                                    , log2fc_cutoff=args.log2_fold_change_cutoff, neighbor_support=True, bayesian=True, permutations=1000, upper_limit=7e6, lower_limit=5e4, adjusted_pvalues_method='distance')
     

     iterative_int,_=DifFracTion.identify_differential_interactions(i_normA, i_normB, args.resolution,
                                                    filter_by='padjusted', pvalue_cutoff=args.pvalue_cutoff
                                                    , log2fc_cutoff=args.log2_fold_change_cutoff, neighbor_support=True, bayesian=True, permutations=1000, upper_limit=7e6, lower_limit=5e4, adjusted_pvalues_method='distance')
     
     alpha_int.to_csv(f"{output_path}/chr{args.chromosome}_{args.resolution}_alpha_normalization_differential_interactions.csv", index=False,sep='\t')
     iterative_int.to_csv(f"{output_path}/chr{args.chromosome}_{args.resolution}_iterative_normalization_differential_interactions.csv", index=False,sep='\t')

if __name__ == "__main__":
     main()