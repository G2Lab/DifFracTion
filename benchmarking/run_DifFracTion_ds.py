import hicstraw
import scipy
from scipy.sparse import coo_matrix
from scipy import ndimage
from scipy.stats import linregress
from statsmodels.stats.multitest import multipletests

import numpy as np
import pandas as pd
import os
import glob

import time
import tracemalloc
import argparse

import matplotlib.pyplot as plt
from matplotlib import colors

from diffraction_edits import utils as DifFracTion_utils
from diffraction_edits import src as DifFracTion

def parse_args():
    parser = argparse.ArgumentParser(description='Run DifFracTion on a given dataset')
    parser.add_argument('--matrixA', type=str, required=True, help='Path to the dataset (matrix or cooler file)')
    parser.add_argument('--matrixB', type=str, required=True, help='Path to the second dataset (matrix or cooler file)')
    parser.add_argument('--chrom', type=str, required=True, help='Chromosome to analyze') 
    parser.add_argument('--resolution', type=int, default=10000, help='Resolution for the analysis (default: 10000)')
    parser.add_argument('--pval', type=float, default=0.05, help='FDR threshold for significant interactions (default: 0.05)')
    parser.add_argument('--type_norm', type=str, default='alpha', help='Type of normalization to apply alpha or iterative (default: alpha)')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results')
    parser.add_argument('--adjusted_pvalues_method', type=str, default='distance', help='Method for adjusting p-values (default: distance)')
    
    return parser.parse_args()


if __name__ == "__main__":
     args = parse_args()
     print(f"Running DifFracTion with the following parameters:\nMatrix A: {args.matrixA}\nMatrix B: {args.matrixB}\nChromosome: {args.chrom}\nResolution: {args.resolution}\nP-value threshold: {args.pval}\nNormalization type: {args.type_norm}\nAdjusted p-values method: {args.adjusted_pvalues_method}\nOutput directory: {args.output_dir}")
     matrix_A=DifFracTion_utils.load_dense_matrix(args.matrixA)
     matrix_B=DifFracTion_utils.load_dense_matrix(args.matrixB)

     

     if args.type_norm == 'alpha':
        tracemalloc.start()
        start_time = time.time()

        norm_matrix_A, norm_matrix_B = DifFracTion.alpha_normalization(matrix_A,matrix_B,
                                                                         args.resolution,metric='median',
                                                                         desired_alpha=1.08)
        normalization_time = time.time() - start_time
        _, norm_mem_peak = tracemalloc.get_traced_memory()
        mem_mb_norm = norm_mem_peak / 1024 / 1024

        alpha_padjusted_interactions,n_bins = DifFracTion.identify_differential_interactions(norm_matrix_A, norm_matrix_B,
                                                                      args.resolution, filter_by='padjusted', pvalue_cutoff=args.pval,
                                                                      log2fc_cutoff=0, neighbor_support=True, adjusted_pvalues_method=args.adjusted_pvalues_method)
        fully_elapsed_sec = time.time() - start_time
        _, mem_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        mem_mb = mem_peak / 1024 / 1024

        perf_log = pd.DataFrame([{'tool': 'DifFracTion_alpha', 'fully_elapsed_sec': fully_elapsed_sec, 'fully_mem_mb': mem_mb,'normalization_time': normalization_time, 'norm_mem_mb': mem_mb_norm}])
        os.makedirs(args.output_dir, exist_ok=True) #/performance
        perf_log.to_csv(f"{args.output_dir}/runtime_DifFracTion_alpha.txt", sep='\t', index=False)
        
        FP=len(alpha_padjusted_interactions)
        FP_support=alpha_padjusted_interactions[alpha_padjusted_interactions['neighbor_support'] == True]
        FPR=FP/n_bins if FP > 0 else 0
        FPR_support=len(FP_support)/n_bins if len(FP_support) > 0 else 0
        final=pd.DataFrame([{'FP': FP, 'FP_with_neighbor_support': len(FP_support), 'FPR': FPR, 'FPR_with_neighbor_support': FPR_support,'n_tests': n_bins}])
        final.to_csv(f"{args.output_dir}/FPR_alpha_padjusted_interactions_summary.txt", sep='\t', index=False)
        #alpha_padjusted_interactions.to_csv(f"{args.output_dir}/alpha_padjusted_interactions.txt", sep='\t', index=False)
     elif args.type_norm == 'iterative':
        tracemalloc.start()
        start_time = time.time()

        norm_matrix_A_iterative, norm_matrix_B_iterative = DifFracTion.iterative_normalization(matrix_A, matrix_B,
                                                                                            args.resolution, metric='median',
                                                                                            weight_factor=-1.08)
        
        normalization_time = time.time() - start_time
        _, norm_mem_peak = tracemalloc.get_traced_memory()
        mem_mb_norm = norm_mem_peak / 1024 / 1024

        iterative_padjusted_interactions,n_bins= DifFracTion.identify_differential_interactions(norm_matrix_A_iterative, norm_matrix_B_iterative,
                                                                      args.resolution, filter_by='padjusted', pvalue_cutoff=args.pval,
                                                                      log2fc_cutoff=0, neighbor_support=True)
        fully_elapsed_sec = time.time() - start_time
        _, mem_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        mem_mb = mem_peak / 1024 / 1024

        perf_log = pd.DataFrame([{'tool': 'DifFracTion_iterative', 'fully_elapsed_sec': fully_elapsed_sec, 'fully_mem_mb': mem_mb, 'normalization_time': normalization_time, 'norm_mem_mb': mem_mb_norm}])
        os.makedirs(args.output_dir, exist_ok=True)
        perf_log.to_csv(f"{args.output_dir}/runtime_DifFracTion_iterative.txt", sep='\t', index=False)
        FP=len(iterative_padjusted_interactions)
        FP_support=iterative_padjusted_interactions[iterative_padjusted_interactions['neighbor_support'] == True]
        FPR=FP/n_bins if FP > 0 else 0
        FPR_support=len(FP_support)/n_bins if len(FP_support) > 0 else 0
        final=pd.DataFrame([{'FP': FP, 'FP_with_neighbor_support': len(FP_support), 'FPR': FPR, 'FPR_with_neighbor_support': FPR_support,'n_tests': n_bins}])
        final.to_csv(f"{args.output_dir}/FPR_iterative_padjusted_interactions_summary.txt", sep='\t', index=False)
     else:
        print(f"Normalization type {args.type_norm} not recognized. Please choose 'alpha' or 'iterative'.")
