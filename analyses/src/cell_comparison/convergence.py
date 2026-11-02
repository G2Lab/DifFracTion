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
    parser.add_argument('--hic_path', type=str, required=True, help='Path to the input Hi-C file (hic format).')
    parser.add_argument('--chromosome', type=str, required=True, help='Chromosome to process (e.g., 1).')
    parser.add_argument('--resolution', type=int, default=10000, help='Resolution of the Hi-C matrices (default: 10000).')
    parser.add_argument('--downsample_factor', type=float, default=1.0, help='Downsample factor for synthetic dataset generation (default: 1.0).')
    parser.add_argument('--metric', type=str, choices=['median', 'mean'], default='median', help='Metric to use for error calculation (default: median).')
    parser.add_argument('--weight_factor', type=float, default=-1.08, help='Weight factor for correction (default: -1.08).')
    parser.add_argument('--analysis_name', type=str, default='default', help='Directory to save the corrected matrices (default: corrected_matrices).')
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

def iterative_correction(matrix_A, matrix_B, resolution,metric='median',weight_factor=-1.08, eta=0.15, tolerance=1e-6, plot=False):
	best_error  = float("inf")
	updated_factor = 1.0
	min_window = 5 # minimum number of distances to check per cycle

	all_distances, counts_A1, counts_B1 = DifFracTion_utils.get_distances_and_counts(matrix_A, matrix_B, resolution)
	ref_counts, to_scale_counts = DifFracTion_utils.set_reference(counts_A1, counts_B1)

	# Remember DifFracTion iterative sets up the alpha of the reference, we do not change it to 1.08
	ref_values, to_scale_values = DifFracTion_utils.calculate_distance_decay(
		all_distances, ref_counts, to_scale_counts, verbose=False)

	n_distances = len(ref_values[0])
	max_i = n_distances - 1  # start with full range

	improved = True
	cycle = 1
	factor_history = [1.0]  # to track the cumulative factor applied across cycles
	while improved:
		no_improvement = 0
		improved = False
		last_improvement_idx = -1

		#i is the current distance index for the correction step
		# So basically it iterates and corrects by distance
		# But the WMSE takes into consideration the effect of the refitting in all the distances, not just the one being corrected
		# Convergence will always be reached because we only update the best error and the counts if there is a significant improvement, so eventually we will reach a point where no further improvements are possible and the loop will stop. The tolerance threshold ensures that we don't get stuck in an infinite loop due to very small improvements that may not be meaningful.
		for i in range(max_i + 1):
			trial_counts, trial_values, trial_error, factor, trial_decay, ref_decay = DifFracTion_utils.correction(
				ref_values, to_scale_values, to_scale_counts, all_distances, i,
				metric=metric, weight_factor=weight_factor, eta=eta)
			# After applying the correction for distance i we check if the error has improved compared to the best error so far. We use a tolerance threshold to determine if the improvement is significant enough to update our best error and continue iterating.
			if trial_error + tolerance < best_error:
				print(f"Improvement at distance {ref_values[0][i]:.2f}: error improved from {best_error:.6f} to {trial_error:.6f} (factor: {factor:.4f})")
				to_scale_counts = trial_counts
				to_scale_values = trial_values
				best_error = trial_error
				updated_factor *= factor
				factor_history.append(updated_factor)
				improved = True
				last_improvement_idx = i
				if plot:
					plt.figure(figsize=(4, 3))
					plt.plot(np.log(ref_values[0]), np.log(ref_decay), label='Reference Decay')
					plt.plot(np.log(ref_values[0]), np.log(trial_decay), label='Scaled Decay')
					plt.xlabel('Log Distance')
					plt.ylabel('Log Value')
					plt.title(f'Cycle {cycle}, Distance {ref_values[0][i]:.2f}, Error: {trial_error:.6f}')
					plt.legend()
					plt.show()
			else:
				no_improvement += 1

		# Window for next cycle to last index that improved, with a minimum of min_window
		if improved:
			max_i = max(last_improvement_idx, min_window - 1)
			print(f"Next cycle will search up to index {max_i} (distance {ref_values[0][max_i]:.2f})")

		cycle += 1
		print(f"Out of {max_i + 1} distances checked, {no_improvement} showed no improvement in this cycle.")
		if not improved:
			print("No improvement in this cycle, stopping correction.")
			break
	
	matrix_B_norm_iterative = matrix_B.copy() * updated_factor
	matrix_A_norm_iterative = matrix_A.copy() * 1.0

	return matrix_A_norm_iterative, matrix_B_norm_iterative, factor_history

def main():
	raw_matrix,name_matrix=load_matrix(args.hic_path, args.chromosome, args.resolution, f"{main_path}/generated_data/")
	output_path = f"{main_path}/{args.analysis_name}/"
	os.makedirs(output_path, exist_ok=True)
	print(f"Results will be saved in: {output_path}")
	matrix_A = raw_matrix.copy()
	matrix_B = DifFracTion_utils.synthetic_datasets(raw_matrix, args.downsample_factor)

	norm_matrix_A_iterative, norm_matrix_B_iterative,factor_history = iterative_correction(matrix_A, matrix_B,
                                                                                            args.resolution, metric='median',
                                                                                            weight_factor=-1.08)
	# save the cycle history of the factors applied
	factor_history_df = pd.DataFrame({'Cycle': range(len(factor_history)), 'Cumulative_Factor': factor_history})
	factor_history_df.to_csv(f"{output_path}/{args.chromosome}_{args.resolution}_{args.downsample_factor}_factor_history", index=False,sep="\t")

if __name__ == "__main__":
	main()
