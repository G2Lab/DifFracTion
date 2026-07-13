import os
import hicstraw 
import cooler

import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from statistics import mean, median
from scipy.stats import linregress
import scipy.sparse
from scipy.sparse import coo_matrix
import math
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from scipy.stats import bootstrap

import hicstraw
from scipy.sparse import coo_matrix
from scipy import ndimage
from scipy.stats import linregress

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import colors


#### DifFracTion functions ####

main_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#  ── Preprocessing - MA related ───────────────────────────────────────────────────────────────

def compute_log2fc(matrix_A,matrix_B,resolution,upper_limit=7e6,lower_limit=5e4):
	n = matrix_A.shape[0]
	#Minumum and maximum distances 
	min_d = int(max(1, lower_limit // resolution)) # minimum distance in bins, avoid main diagonal d=0
	max_d = int(min(n, upper_limit // resolution))  # don't go beyond the upper limit in bins

	distances_list, log2fc_list, rows_in_list, cols_in_list,IF_list = [], [], [], [], []
	distances_out_list, log2fc_out_list, rows_out_list, cols_out_list,IF_out_list = [], [], [], [], []
	
	for d in range(1, n): # avoid main diagonal
		diag_A = np.asarray(matrix_A.diagonal(d), dtype=np.float32).ravel()
		diag_B = np.asarray(matrix_B.diagonal(d), dtype=np.float32).ravel()	
		valid = (diag_A > 0) & (diag_B > 0)
		valid = valid.ravel()
		rows_valid = np.where(valid)[0]
		
		cols_valid = rows_valid + d  # Assuming square matrix
		if not np.any(valid):
			continue

		log2fc = np.log2((diag_A[valid] + 0.01) / (diag_B[valid] + 0.01))
		mean_IF = (diag_A[valid] + diag_B[valid]) / 2
		d_bp = d * resolution

		if min_d <= d < max_d:   # within range
			log2fc_list.append(log2fc)
			distances_list.append(np.full(log2fc.size, d_bp, dtype=np.float32))
			IF_list.append(np.full(log2fc.size, mean_IF, dtype=np.float32))
			rows_in_list.append(rows_valid)
			cols_in_list.append(cols_valid)
		else:                    # outside range
			log2fc_out_list.append(log2fc)
			distances_out_list.append(np.full(log2fc.size, d_bp, dtype=np.float32))
			IF_out_list.append(np.full(log2fc.size, mean_IF, dtype=np.float32))
			rows_out_list.append(rows_valid)
			cols_out_list.append(cols_valid)
	obs_rows = np.concatenate(rows_in_list + rows_out_list) if (rows_in_list or rows_out_list) else np.array([])
	obs_cols = np.concatenate(cols_in_list + cols_out_list) if (cols_in_list or cols_out_list) else np.array([])
	distances     = np.concatenate(distances_list)
	log2fc_values = np.concatenate(log2fc_list)
	distances_out     = np.concatenate(distances_out_list)     if distances_out_list     else np.array([])
	log2fc_values_out = np.concatenate(log2fc_out_list)        if log2fc_out_list        else np.array([])
	IF_values         = np.concatenate(IF_list)               if IF_list               else np.array([])
	IF_values_out     = np.concatenate(IF_out_list)           if IF_out_list           else np.array([])

	return distances, log2fc_values, distances_out, log2fc_values_out, obs_rows, obs_cols, IF_values, IF_values_out


def plot_MA(matrix_A: np.ndarray, 
		  matrix_B: np.ndarray, 
		  resolution: int, 
		  log2_fc_cutoff: float, 
		  upper_limit: float = 7e6, 
		  lower_limit: float = 5e4, 
		  xscale: str = 'log', 
		  bayesian: bool = True):
	if bayesian:
		distances, log2_fc, distances_out, log2_fc_out, obs_rows, obs_cols, IF_values, IF_values_out = compute_log2fc_bayesian(matrix_A, matrix_B, resolution, lower_limit, upper_limit)
	else:
		distances, log2_fc, distances_out, log2_fc_out, obs_rows, obs_cols, IF_values, IF_values_out = compute_log2fc(matrix_A, matrix_B, resolution, upper_limit, lower_limit)

	# If there are no points outside the range, we can just plot the points within the range
	all_distances = np.concatenate([distances, distances_out]) if distances_out.size else distances
	all_log2fc    = np.concatenate([log2_fc,   log2_fc_out])   if log2_fc_out.size   else log2_fc
	all_IF        = np.concatenate([IF_values, IF_values_out]) if IF_values_out.size else IF_values

	IF_90th_percentile = np.percentile(all_IF, 90)

	fig = plt.figure(figsize=(8, 6), dpi=200)

	plt.scatter(all_distances, all_log2fc, alpha=0.5, s=1, 
		   		c=all_IF, cmap='RdYlBu_r',vmax=IF_90th_percentile)
	plt.colorbar(label='Mean Interaction Frequency (IF)')
	plt.hlines(0,  xmin=0, xmax=all_distances.max(), colors='red',  linestyles='dashed')
	plt.hlines( log2_fc_cutoff, xmin=0, xmax=all_distances.max(), colors='gray', linestyles='dashed')
	plt.hlines(-log2_fc_cutoff, xmin=0, xmax=all_distances.max(), colors='gray', linestyles='dashed')
	plt.vlines(upper_limit, ymin=all_log2fc.min(), ymax=all_log2fc.max(), colors='red', linestyles='dotted', linewidth=0.5)
	if lower_limit > 0:
		plt.vlines(lower_limit, ymin=all_log2fc.min(), ymax=all_log2fc.max(), colors='red', linestyles='dotted', linewidth=0.5)
	plt.xlabel('Genomic Distance (bp)')
	if xscale == 'log':
		plt.xscale('log')
	
	plt.ylabel('Log2 Fold Change (Sample A / Sample B)')

	return fig







# ── Preprocessing data functions ───────────────────────────────────────────────────────────────

def save_sparse_matrix(hic_path,chromosome,resolution,output_path):
	'Load from .hic file'
	mdz = prepare_MatrixZoomData(hic_path,chromosome,resolution)
	chr_len = get_chromosome_length(hic_path,chromosome)

	n_bins = chr_len //	resolution + 1

	records = mdz.getRecords(0, chr_len,0,chr_len)

	rows,cols,values = [],[],[]
	for record in records:
		# Might have to implement if it is on coordinates or bins
		i,j = record.binX//resolution, record.binY//resolution
		rows.append(i); cols.append(j); values.append(record.counts)
	sparse_matrix = scipy.sparse.csr_matrix(
		(values, (rows, cols)), shape=(n_bins, n_bins)
	)

	scipy.sparse.save_npz(output_path, sparse_matrix)

def load_sparse_matrix(npz_path):
	return scipy.sparse.load_npz(npz_path)

def load_dense_matrix(npz_path):
	sparse_matrix = load_sparse_matrix(npz_path)
	dense_matrix = np.asarray(sparse_matrix.todense(), dtype=np.float32)
	return dense_matrix

def get_chromosome_length(hic_path,chromosome):
	"""

	helper function to get the length of a chromosome from a .hic file	

	
	:param hic_path: Local path to the .hic file
	:param chromosome: Chromosome name (e.g., 'chr1')
	"""
	hic = hicstraw.HiCFile(hic_path)
	chromosomes = hic.getChromosomes()
	for chrom in chromosomes:
		if chrom.name == chromosome:
			chromosome_length = chrom.length
			return int(chromosome_length)
		else:
			continue
	raise ValueError(f"Chromosome {chromosome} not found in the Hi-C file.")

def prepare_MatrixZoomData(hic_path,chromosome,resolution,normalization="KR"):
    """
    Prepare a hicstraw MatrixZoomData object for a given chromosome, resolution, and normalization.
    
    Parameters:
    hic_path (str): Path to the .hic file.
    chromosome (str): Chromosome name (e.g., 'chr1').
    resolution (int): Resolution in base pairs (e.g., 10000 for 10kb).
    normalization (str): Normalization method (default is "KR").
    
    Returns:
    MatrixZoomData: The prepared MatrixZoomData object.
    """
    straw = hicstraw.HiCFile(hic_path)
    matrix_data = straw.getMatrixZoomData(chromosome, chromosome, "observed", normalization,"BP", resolution)
    return matrix_data


# ────── Synthetic datasets functions ────────────────────────────────────────────────────────────────	
# Synthetic data generation
def synthetic_datasets(matrix, ds_factor):
	"""downsampling method
	ds_factor = downsampling factor, e.g. 0.5 means half the reads"""
	if ds_factor == 1:
		return matrix.copy()
	else:
		n_reads = get_reads_count(matrix) * ds_factor
		tag_mat= dense2tag(matrix) #tag_len is the number of non-zero elements in the matrix
		sample_idx = np.random.choice(len(tag_mat), int(n_reads)) #from tag len randomly select n_reads number of items
		sample_tag = tag_mat[sample_idx]# getting coordinates from the coo matrix of the selected items
		down_mat = tag2dense(sample_tag, matrix.shape[0])
		down_mat = np.matrix(down_mat) # convert to numpy matrix
		return down_mat
	
def get_reads_count(matrix):
	'''Get the total number of reads (upper triangle) in a contact matrix.'''
	total_reads = int(np.sum(np.triu(matrix, k=0)))
	return total_reads

def dense2tag(matrix):
    """Convert dense square matrix to COO-based tag matrix (row, col per count)."""
    matrix = np.round(matrix).astype(int)
    matrix = np.triu(matrix)  # keep upper triangle
    coo_mat = coo_matrix(matrix)
    row, col, data = coo_mat.row, coo_mat.col, coo_mat.data

    # repeat rows and cols by data counts
    tag_mat = np.column_stack((np.repeat(row, data), np.repeat(col, data)))
    return tag_mat

def tag2dense(tag, nsize):
    """Convert tag matrix (row, col repeated) back to dense symmetric square matrix."""
    row, col = tag[:, 0], tag[:, 1]
    data = np.ones(len(tag), dtype=int)
    # build upper-triangular COO matrix directly
    dense_mat = coo_matrix((data, (row, col)), shape=(nsize, nsize)).toarray()
    # mirror to make symmetric
    dense_mat = dense_mat 
    return dense_mat

def matrix_downsampling(matrix, n_reads):
	"""From an initial matrix, downsample n_reads and return the downsampled matrix with the same distribution of reads."""
	tag_mat= dense2tag(matrix) #tag_len is the number of non-zero elements in the matrix
	sample_idx = np.random.choice(len(tag_mat), int(n_reads)) #from tag len randomly select n_reads number of items
	sample_tag = tag_mat[sample_idx]# getting coordinates from the coo matrix of the selected items
	down_mat = tag2dense(sample_tag, matrix.shape[0])
	down_mat = np.matrix(down_mat)
	return down_mat


# ── Alpha Normalization functions ──────────────────────────────────────────────────────────────
def _CI_bootstrap(data, n_resamples=100, confidence_level=0.95) -> bootstrap.BootstrapResult:
	'''Calculate the confidence interval for the median of the data using bootstrap resampling.
	:param data: 1D array of data points
	:param n_resamples: number of bootstrap resamples. Default is 500
	:param confidence_level: confidence level for the interval. Default is 0.95 
	'''
	data = data[data > 0]  # Filter out non-positive or zero values
	R = bootstrap((data,), np.median, confidence_level=confidence_level, n_resamples=n_resamples, method='percentile')

	return R

def _calculate_cutoffs(matrix_A: np.ndarray, matrix_B: np.ndarray, resolution: int,upper_limit: float = 7e6):
	'''Calculate cutoffs for matrices A and B based on the upper limit. It uses a confidence interval bootstraping technique
	which allow us to get a robust estimate of the median of the counts outside the distance decay, therefore allowing us
	to set up a "noise" cutoff.
	Args:
	:param matrix_A: dense matrix A
	:param matrix_B: dense matrix B
	:param resolution: resolution of the matrices
	:param upper_limit: upper limit for the distance
	'''
	n=matrix_A.shape[0]

	min_d = int(min(n, upper_limit // resolution)) # because we want to get the bins outside the distance decay range
	max_d = int(max(n, upper_limit // resolution)) # because we want to get the bins outside the distance decay range
	out_A, out_B = [], []
	
	# Get all counts outside the distance decay range (overdispersion)
	for d in range(min_d, max_d):
		diag_A = np.asarray(matrix_A.diagonal(d), dtype=np.float32).ravel()
		diag_B = np.asarray(matrix_B.diagonal(d), dtype=np.float32).ravel()
		out_A.append(diag_A)
		out_B.append(diag_B)


	out_A = np.concatenate(out_A)
	out_B = np.concatenate(out_B)

	# We have each matrix out of distance points, we calculate a CI
	# Calculate CI
	R_A = _CI_bootstrap(out_A[out_A > 0])
	R_B = _CI_bootstrap(out_B[out_B > 0])

	cutoff_A = R_A.confidence_interval.high if R_A.confidence_interval.high > 0 else 0
	cutoff_B = R_B.confidence_interval.high if R_B.confidence_interval.high > 0 else 0

	return cutoff_A, cutoff_B

def get_distances_and_counts_paired(matrix_A: np.ndarray, matrix_B: np.ndarray, resolution: int,
                                     upper_limit: float = 7e6, lower_limit: float = 5e4, filter: bool = True
                                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    
    if filter:
        cutoff_A, cutoff_B = _calculate_cutoffs(matrix_A, matrix_B, resolution, upper_limit)
    else:
        cutoff_A, cutoff_B = 0, 0

    n = matrix_A.shape[0]

    min_d = int(max(1, lower_limit // resolution))
    max_d = int(min(n, upper_limit // resolution))

    # Everything inside the range but that also has counts above the cutoff for both matrices
    # This is just for modeling the distance decay, we havent removed any bins yet
    distances_list, counts_A, counts_B = [], [], []
    for d in range(min_d, max_d):
        diag_A = np.asarray(matrix_A.diagonal(d), dtype=np.float32).ravel()
        diag_B = np.asarray(matrix_B.diagonal(d), dtype=np.float32).ravel()

        valid = (diag_A > cutoff_A) & (diag_B > cutoff_B)
        if not np.any(valid):
            continue

        counts_A.append(diag_A[valid])
        counts_B.append(diag_B[valid])
        distances_list.append(np.full(valid.sum(), d * resolution, dtype=np.float32))

    return (np.concatenate(distances_list), np.concatenate(counts_A), np.concatenate(counts_B))

def get_distances_and_counts_individual(matrix: np.ndarray, resolution: int, 
								upper_limit: float = 7e6, 
								lower_limit: float = 5e4) -> tuple[np.ndarray, np.ndarray]:
	n = matrix.shape[0]
	min_d = int(max(1, lower_limit // resolution))
	max_d = int(min(n, upper_limit // resolution))

	distances_list, counts_list = [], []
	for d in range(min_d, max_d):
		diag = np.asarray(matrix.diagonal(d), dtype=np.float32).ravel()
		valid = diag > 0
		if not np.any(valid):
			continue

		counts_list.append(diag[valid])
		distances_list.append(np.full(valid.sum(), d * resolution, dtype=np.float32))

	return np.concatenate(distances_list), np.concatenate(counts_list)


def calculate_distance_decay_individual(distances: np.ndarray, 
							     counts_array: np.ndarray, 
								metric: str='median', verbose=True):
	'''Calculate distance decay for a single matrix
	Args:
		:param distances: array of distances (one value per bin), within the genomic range
		:param counts_array: array of counts (must be the same length as distances), within the genomic range
		:param metric: metric to use ('median' or 'mean') (what to use to calculate the decay per distance)
		:param verbose: whether to print verbose output
	Returns: 
		# Distances, - linear 0
		# Observed Metric Values, - linear 1
		# Expected Fitted values from power law, -linear 2
		# alpha - scalar (log) 3
		# intercept - log 4
		# ------------
		# Log - log plots
		# -------------
		# log distances 5
		# log values 6
'''
	unique_distances=np.unique(distances) 
	labels = np.searchsorted(unique_distances,distances) #All distances are assigned a label corresponding to their unique distance index

	if metric == 'median':
		idx = np.arange(len(unique_distances))
		#Map each count to its corresponding unique distance index
		metric_decay = ndimage.median(counts_array,labels=labels,index=idx)

	elif metric == 'mean':
		# sum of counts (weights let us sum the counts, otherwise is occurrences of labels) for each unique distance/ label 
		A = np.bincount(labels,weights=counts_array)
		n = np.bincount(labels)
		valid = n > 0
		metric_decay = A[valid] / n[valid]

	logx   = np.log(unique_distances)   # x   = log(d)
	logy   = np.log(metric_decay)      # y   = log(P(d))

	valid_values = np.isfinite(logy) & np.isfinite(logx)
	logx_v, logy_v = logx[valid_values], logy[valid_values]

	slope, intercept, r, *_ = linregress(logx_v, logy_v)

	alpha = slope
	fit = np.exp(intercept) * unique_distances[valid_values] ** alpha

	if verbose:
		print(f"α = {alpha:.2f}, intercept = {intercept:.2f}")

	return [unique_distances[valid_values], metric_decay[valid_values], fit, alpha, intercept, logx_v, logy_v]

def refit_intercept(metric_values: np.ndarray, distances_array: np.ndarray, desired_alpha: float, form='log') -> float:
    '''Refit the decay curve to a desired alpha and return the expected counts.
    :param form='log'    : inputs already log-transformed (from calculate_distance_decay [5,6])
    :param form='linear' : inputs in linear space         (from calculate_distance_decay [0,1])
    '''
    m = -abs(desired_alpha)   # normalize sign — slope is always negative

    if form == 'log':
        logy, logd = metric_values, distances_array
    elif form == 'linear':
        logy = np.log(metric_values)
        logd = np.log(distances_array)

    b = np.mean(logy - m * logd)
    return b

def shared_intercept(A_values: list, B_values: list, desired_alpha: float) -> float:
	'''Find the single intercept b that minimizes error for both A and B jointly.
	Pools both log-count arrays and solves for b at the desired slope simultaneously.
	This is the least-squares solution for b given fixed slope m, weighted by number of points.
	When datasets are the same size, this equals the mean of the two individual intercepts.
	When datasets differ in size, it naturally weights toward the larger/denser dataset.
	Args:
		:param A_values: list from calculate_distance_decay for dataset A
		:param B_values: list from calculate_distance_decay for dataset B
		:param desired_alpha: desired power-law exponent (positive or negative convention both accepted)
	Returns:
		b_shared: shared intercept in log space
	'''
	m = -abs(desired_alpha)   # normalize sign — slope is always negative

	# Pool residuals from both datasets: b_i = logy_i - m·logx_i
	# Taking the mean across all points jointly is the least-squares solution for b
	b_shared = np.mean(np.concatenate([
		A_values[6] - m * A_values[5],   # logy_A - m·logx  (residuals for A)
		B_values[6] - m * B_values[5]    # logy_B - m·logx  (residuals for B)
	]))
	return b_shared

def refit_decay(b: float, desired_alpha: float, distances_array: np.ndarray, form='log') -> np.ndarray:
	'''Given the new intercept b and the desired alpha, return the expected counts for each distance in the array.
	Args:
		:param b: new intercept in log space
		:param desired_alpha: desired power-law exponent
		:param distances_array: array of distances
		:param form: 'log' or 'linear'
	Returns:
		fit: array of expected counts
	'''
	m = -abs(desired_alpha)   # normalize sign
	if form == 'log':
		logd = distances_array
	elif form == 'linear':
		logd = np.log(distances_array)
	return np.exp(b + m * logd)   # P(d) = e^b · d^m  back in linear space

def get_factors(original:list, fit:list) -> np.ndarray:
	'''Calculate the factors by which we need to multiply the original counts to get the fitted counts.
	Args:
		original: list of original counts
		fit: list of fitted counts
	Returns:
		factors: array of factors
	'''
	factors = fit / original
	return factors

def apply_decay_factors(matrix: np.ndarray, factors: np.ndarray, distances: np.ndarray, resolution: int, outside_range='interp'):
	'''Apply decay factors to a matrix based on the distance of the bins and the factors per distance.
	Args:
		matrix: dense matrix to which we want to apply the factors
		factors: array of factors per distance (same size as distances)
		distances: array of unique distances in bp corresponding to the factors inside the fitted range
		resolution: bin resolution in base pairs
		outside_range: how to handle distances outside the fitted range
			'mean'  : apply the mean of all computed factors (neutral depth correction,
			          no shape assumption — recommended for distances where power law doesn't hold)
			'interp': clamp to the nearest boundary factor (assumes decay trend continues)
	Returns:
		normalized: dense matrix with the factors applied'''

	n = matrix.shape[0]
	normalized = matrix.astype(np.float32).copy()

	# All diagonal distances in bp: d=1  resolution, d=2  2*resolution, ...
	# 0 indexed but starting on d=1
	d_bp_all = np.arange(1, n, dtype=np.float64) * resolution   # shape (n-1,)

	if outside_range == 'mean':
		# Within the fitted range → interpolate factor from the decay model
		# Outside the fitted range → apply mean factor (neutral depth correction)
		mean_factor = np.mean(factors)
		d_min, d_max = distances[0], distances[-1]

		within = (d_bp_all >= d_min) & (d_bp_all <= d_max)
		#Initialize the array of factors with the mean factor, then override the values within the fitted range with the interpolated factors
		all_factors = np.full(n - 1, mean_factor, dtype=np.float32)         # default: mean
		all_factors[within] = np.interp(d_bp_all[within], distances, factors)  # override within range
	elif outside_range == 'median':
		# Within the fitted range → interpolate factor from the decay model
		# Outside the fitted range → apply median factor (neutral depth correction)
		median_factor = np.median(factors)
		d_min, d_max = distances[0], distances[-1]

		within = (d_bp_all >= d_min) & (d_bp_all <= d_max)
		#Initialize the array of factors with the median factor, then override the values within the fitted range with the interpolated factors
		all_factors = np.full(n - 1, median_factor, dtype=np.float32)         # default: median
		all_factors[within] = np.interp(d_bp_all[within], distances, factors)  # override within range
	elif outside_range == 'interp':
		# the data to extrapolate, the data we know (x), and the factors we know (y)
		# So all distances will have a factor
		last_factor = factors[-1]
		first_factor = factors[0]
		d_min, d_max = distances[0], distances[-1]

		within = (d_bp_all >= d_min) & (d_bp_all <= d_max)
		#Initialize the array of factors with the last factor, then override the values within the fitted range with the interpolated factors
		all_factors = np.full(n - 1, last_factor, dtype=np.float32)         # default: last factor

		all_factors[within] = np.interp(d_bp_all[within], distances, factors).astype(np.float32)
		all_factors[d_bp_all < d_min] = first_factor  # clamp to first factor for distances below the fitted range
		all_factors[d_bp_all > d_max] = last_factor   # clamp to last factor for distances above the fitted range

	# Apply factors — one factor per diagonal, vectorized interp already done above
	for d in range(1, n): # never touch the main diagonal d=0
		i = np.arange(n - d) # We get a vector of indices for the diagonal d
		# adding i+d allow us to get the corresponding column
		# d = 2, n =10
		# i = [0,1,2,3,4,5,6,7]  # rows
		# i+d = [2,3,4,5,6,7,8,9]  # columns, which allow us to get the diagonal d=2
		normalized[i, i + d] *= all_factors[d - 1]
		
	normalized = np.asarray(normalized, dtype=np.float32).copy()

	return normalized

# ── Iterative Normalization functions ──────────────────────────────────────────────────────────
def set_reference(counts_A,counts_B):
	if counts_A.sum() >= counts_B.sum():
		ref_counts = counts_A
		to_scale_counts = counts_B
	else:
		ref_counts = counts_B
		to_scale_counts = counts_A
	return ref_counts, to_scale_counts

def get_distances_and_counts(matrix_A, matrix_B, resolution, upper_limit=7e6, lower_limit=5e4):
    n = matrix_A.shape[0]
    #min_d max between 1 and the lower limit in bins, to avoid main diagonal d=0
    min_d = int(max(1, lower_limit // resolution))
    #max_d min between n and the upper limit in bins, to avoid going beyond the upper limit
    max_d = int(min(n, upper_limit // resolution))

    distances_list, counts_A, counts_B = [], [], []
    for d in range(min_d, max_d):
        diag_A = np.asarray(matrix_A.diagonal(d), dtype=np.float32).ravel()
        diag_B = np.asarray(matrix_B.diagonal(d), dtype=np.float32).ravel()

        # Each matrix contributes only bins where IT has real counts
        # Using intersection: both must have real data to be included in the decay fit
        valid = (diag_A > 0) & (diag_B > 0)
        if not np.any(valid):
            continue

        counts_A.append(diag_A[valid])
        counts_B.append(diag_B[valid])
        distances_list.append(np.full(valid.sum(), d * resolution, dtype=np.float32))

    return np.concatenate(distances_list), np.concatenate(counts_A), np.concatenate(counts_B)

def calculate_distance_decay(distances,counts_A_array,counts_B_array,metric='median',verbose=True):
	'''Calculates distance decay based on counts that lay between the distance decay limits.
		By default it uses median to calculate it'''
	
	unique_distances=np.unique(distances)
	#Map distances to dense integer labels i.e. the first distance
	#in unique distances get 0, and it is assigned to all elements in distances that are that distance
	labels = np.searchsorted(unique_distances,distances)

	if metric == 'median':
		#Index for distance
		idx = np.arange(len(unique_distances))
		# Calculates median of the values of an array over labeled regions
		# Since labels already link the information of counts to distances 
		# i.e. same distance have the same label, we avoid using a loop
		metric_decay_A = ndimage.median(counts_A_array,labels=labels,index=idx)
		metric_decay_B = ndimage.median(counts_B_array,labels=labels,index=idx)

	elif metric == 'mean':
	
		#bincount Counts the number of occurrences on labels, how many bins we have per distance
		#Which corresponds to the size of unique distances
		#Since each element in labels corresponds to one count we can use weights 
		#Basically is gonna sum the counts that have the same label 
		# A and B are arrays with a lenght equal to the amount of unique distances
		# The values inside each array are the sum of all bins that lay a given distance
		A = np.bincount(labels,weights=counts_A_array)
		B = np.bincount(labels,weights=counts_B_array)

		n = np.bincount(labels)
		valid = n > 0 # just in case, we don't want to encounter problems by dividing by 0

		# For mean it will return an array which contains a mean per distance
		metric_decay_A = A[valid] / n[valid]
		metric_decay_B = B[valid] / n[valid]

	# ── Log-log regression to extract chromatin alpha ──────────────────────────
	# Power law:        P(d) = c · d^α
	# Linearized:   log P(d) = α · log(d) + log(c)
	# Linear form:          y = m · x     + b
	# ───────────────────────────────────────────────────────────────────────────

	logx   = np.log(unique_distances)   # x   = log(d)
	logy_A = np.log(metric_decay_A)     # y_A = log(P_A(d))
	logy_B = np.log(metric_decay_B)     # y_B = log(P_B(d))

	# Remove points where log is undefined (counts were 0 → log(0) = -inf)
	valid_values = np.isfinite(logy_A) & np.isfinite(logy_B) & np.isfinite(logx)
	logx_v, logy_A_v, logy_B_v = logx[valid_values], logy_A[valid_values], logy_B[valid_values]

	# Fit y = m·x + b  in log-log space
	slope_A, intercept_A, r_A, *_ = linregress(logx_v, logy_A_v)  # slope_A = α_A,  intercept_A = log(c_A)
	slope_B, intercept_B, r_B, *_ = linregress(logx_v, logy_B_v)  # slope_B = α_B,  intercept_B = log(c_B)

	alpha_A, alpha_B = slope_A, slope_B   # α = m  (slope IS the power-law exponent)

	# Back-transform fitted line to original space: P(d) = c · d^α
	# intercept = log(c)  →  c = e^(intercept)
	# fit = e^(intercept) · d^α  =  c · d^α
	fit_A = np.exp(intercept_A) * unique_distances[valid_values] ** alpha_A
	fit_B = np.exp(intercept_B) * unique_distances[valid_values] ** alpha_B

	if verbose:
		print(f"α_A = {alpha_A:.2f}, intercept_A = {intercept_A:.2f}")
		print(f"α_B = {alpha_B:.2f}, intercept_B = {intercept_B:.2f}")
	
	#Returns: Distances, - linear 0
			# Observed Metric Values, - linear 1
			# Expected Fitted values from power law, -linear 2
			# alpha - scalar (log) 3
			# intercept - log 4
			# ------------
			# Log - log plots
			# -------------
			# log distances 5
			# log values 6

	A_values = [unique_distances[valid_values], 
		   		metric_decay_A[valid_values], 
				fit_A, 
				alpha_A, 
				intercept_A,
				logx_v,
				logy_A_v]
	B_values = [unique_distances[valid_values], 
		   		metric_decay_B[valid_values], 
				fit_B, 
				alpha_B, 
				intercept_B,
				logx_v,
				logy_B_v]
	return A_values,B_values

def iterative_correction(matrix_A, matrix_B, resolution,metric='median',
					weight_factor=-1.08, eta=0.15, tolerance=1e-6, filter=True,plot=False):
	'''Iteratively corrects the distance decay of matrix_B to match that of matrix_A.
	Args:
		:param matrix_A: reference matrix (dense)
		:param matrix_B: matrix to be corrected (dense)
		:param resolution: bin resolution in base pairs
		:param metric: metric to use for distance decay ('median' or 'mean')
		:param weight_factor: factor for weighting the correction. Default is -1.08, which prioritizes short-range interactions.
		:param eta: learning rate for the correction. Default is 0.15, which controls the step size of the correction.
		:param tolerance: tolerance for convergence. Default is 1e-6, which determines when to stop iterating based on the improvement in error.
		:param filter: whether to filter the matrices based on CI cutoffs. Default is True, which removes low-count bins that may introduce noise in the distance decay estimation.
		:param plot: whether to plot the correction
	'''
	best_error  = float("inf")
	updated_factor = 1.0
	min_window = 10 # minimum number of distances to check per cycle

	all_distances, counts_A1, counts_B1 = get_distances_and_counts_paired(matrix_A, matrix_B, resolution,filter=filter)
	ref_counts, to_scale_counts = set_reference(counts_A1, counts_B1)
	

	ref_values, to_scale_values = calculate_distance_decay(
		all_distances, ref_counts, to_scale_counts, verbose=False)

	n_distances = len(ref_values[0])
	max_i = n_distances - 1  # start with full range

	improved = True
	cycle = 1
	while improved:
		no_improvement = 0
		improved = False
		last_improvement_idx = -1

		#i is the current distance index for the correction step
		# So basically it iterates and corrects by distance
		# But the WMSE takes into consideration the effect of the refitting in all the distances, not just the one being corrected
		# Convergence will always be reached because we only update the best error and the counts if there is a significant improvement, so eventually we will reach a point where no further improvements are possible and the loop will stop. The tolerance threshold ensures that we don't get stuck in an infinite loop due to very small improvements that may not be meaningful.
		for i in range(max_i + 1):
			trial_counts, trial_values, trial_error, factor, trial_decay, ref_decay = correction(
				ref_values, to_scale_values, to_scale_counts, all_distances, i,
				metric=metric, weight_factor=weight_factor, eta=eta)
			# After applying the correction for distance i we check if the error has improved compared to the best error so far. We use a tolerance threshold to determine if the improvement is significant enough to update our best error and continue iterating.
			if trial_error + tolerance < best_error:
				print(f"Improvement at distance {ref_values[0][i]:.2f}: error improved from {best_error:.6f} to {trial_error:.6f} (factor: {factor:.4f})")
				to_scale_counts = trial_counts
				to_scale_values = trial_values
				best_error = trial_error
				updated_factor *= factor
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
	
	return matrix_A_norm_iterative, matrix_B_norm_iterative

def correction(ref_values, to_scale_values, to_scale_counts, distances, i,
			   metric='median', weight_factor=-1.08, eta=0.15):
	unique_distances = ref_values[0]
	ratios = ref_values[1] / to_scale_values[1]
	weights = unique_distances ** weight_factor
	weights = weights / weights.max()

	raw_factor = ratios[i] * weights[i]
	factor = 1 + eta * (raw_factor - 1)

	trial_counts = to_scale_counts * factor
	trial_values = calculate_distance_decay_individual(
		distances, trial_counts, metric=metric, verbose=False)

	ref_decay = ref_values[2]
	trial_decay = trial_values[2]

	# this is the error of all distances, weighted by distance to prioritize fitting the short-range decay
	def error(scaled_decay, ref_decay, error_type='mse'):
		if error_type == 'mse':
			se = ((np.log(scaled_decay) - np.log(ref_decay)) ** 2)
			w = unique_distances ** weight_factor
			w = w / w.sum()
			return np.average(se, weights=w)
		else:
			raise ValueError(f"Unsupported error type: {error_type}")

	trial_error = error(trial_decay, ref_decay)
	return trial_counts, trial_values, trial_error, factor, trial_decay, ref_decay

# ── Differential Interactions ──────────────────────────────────────────────────────────────────
# MAIN
def differential_interactions(norm_matrix_A: np.ndarray, norm_matrix_B: np.ndarray, resolution: float, permutations_n: int = 1000, pvalue_cutoff: float = 0.05, log2fc_cutoff: float = 0,
						filter_by: str = 'pvalue', neighbor_support: bool = True, bayesian: bool = True,
						upper_limit: float = 7e6, lower_limit: float = 5e4, adjusted_pvalues_method: str = 'distance', do_IF_cutoff: bool = True):
	'''
	Docstring for differential_interactions
	
	:param norm_matrix_A: Description
	:type norm_matrix_A: np.ndarray
	:param norm_matrix_B: Description
	:type norm_matrix_B: np.ndarray
	:param matrix_A: Description
	:type matrix_A: np.ndarray
	:param matrix_B: Description
	:type matrix_B: np.ndarray
	:param resolution: Description
	:type resolution: float
	:param permutations_n: Description
	:type permutations_n: int
	:param pvalue_cutoff: Description
	:type pvalue_cutoff: float
	:param log2fc_cutoff: Description
	:type log2fc_cutoff: float
	:param filter_by: Description
	:type filter_by: str
	:param neighbor_support: Description
	:type neighbor_support: bool
	:param bayesian: Description
	:type bayesian: bool
	:param upper_limit: Description
	:type upper_limit: float
	:param lower_limit: Description
	:type lower_limit: float
	:param adjusted_pvalues_method: Description
	:type adjusted_pvalues_method: str
	:param do_IF_cutoff: Description
	:type do_IF_cutoff: bool
	'''
	if do_IF_cutoff:
		rows, cols, log2fc, pvalues, distances = get_differential_interactions(norm_matrix_A, norm_matrix_B, resolution, 
															  permutations_n=permutations_n, bayesian=bayesian,
																 upper_limit=upper_limit, lower_limit=lower_limit,do_IF_cutoff=True)
	else:
		rows, cols, log2fc, pvalues, distances = get_differential_interactions(norm_matrix_A, norm_matrix_B, resolution, 
															  permutations_n=permutations_n, bayesian=bayesian,
																 upper_limit=upper_limit, lower_limit=lower_limit,do_IF_cutoff=False)
	print(f"[1] Total bin pairs tested                          : {len(rows):,}")
	total_tested = len(rows)

	# If the size of all the elements is not the same, quit
	if not (len(rows) == len(cols) == len(log2fc) == len(pvalues) == len(distances)):
		raise ValueError("The lengths of rows, cols, log2fc, pvalues, and distances must be the same.")
		sys.exit(1)
		
	if adjusted_pvalues_method == 'distance':
		print(f"Adjusting p-values by distance using method: {adjusted_pvalues_method}")
		# Size is the same as pvalues, i.e, the number of bins that passed the cutoff
		adjusted_pvalues = adjust_pvalues_by_distance(pvalues, distances, method='fdr_bh')
		print(adjusted_pvalues.shape)
	elif adjusted_pvalues_method == 'global':
		print(f"Adjusting p-values globally using method: {adjusted_pvalues_method}")
		adjusted_pvalues = adjust_pvalues_global(pvalues, method='fdr_bh')
		print(adjusted_pvalues.shape)
	else:
		raise ValueError(f"Unsupported adjusted p-values method: {adjusted_pvalues_method}")

	# Filter by raw pvalue (used when filter_by='pvalue')
	significant_rows, significant_cols, significant_log2fc, significant_pvalues, significant_adjusted_pvalues, significant_distances = filter_significant_interactions(
	rows, cols, log2fc, pvalues, adjusted_pvalues, distances, pvalue_cutoff=pvalue_cutoff, log2fc_cutoff=log2fc_cutoff
	)
	print(f"[2] Significant by pvalue<={pvalue_cutoff} & |log2fc|>={log2fc_cutoff} : {len(significant_rows):,}")

	# Filter by adjusted pvalue independently on the full set (not pre-filtered by raw pvalue)
	mask_padj_full = (adjusted_pvalues <= pvalue_cutoff) & (np.abs(log2fc) >= log2fc_cutoff)
	padj_rows      = rows[mask_padj_full]
	padj_cols      = cols[mask_padj_full]
	padj_log2fc    = log2fc[mask_padj_full]
	padj_pvalues   = pvalues[mask_padj_full]
	padj_adjusted  = adjusted_pvalues[mask_padj_full]
	padj_distances = distances[mask_padj_full]
	print(f"[3] Significant by padjusted<={pvalue_cutoff} & |log2fc|>={log2fc_cutoff}: {mask_padj_full.sum():,}")

	#All of these are just significant by pvalue and logfold
	sig_by_pval = pd.DataFrame({
	'row': significant_rows,
	'col': significant_cols,
	'log2fc': significant_log2fc,
	'pvalue': significant_pvalues,
	'adjusted_pvalue': significant_adjusted_pvalues,
	'distance': significant_distances
	})

	sig_by_pval = neighbor_support_f(sig_by_pval)
	n_support_pval = sig_by_pval['neighbor_support'].sum()
	print(f"[4] With neighbor support (of pval-significant)     : {n_support_pval:,}")

	# Build padjusted DataFrame independently from the full set
	sig_by_padj = pd.DataFrame({
		'row': padj_rows, 'col': padj_cols, 'log2fc': padj_log2fc,
		'pvalue': padj_pvalues, 'adjusted_pvalue': padj_adjusted, 'distance': padj_distances
	})
	sig_by_padj = neighbor_support_f(sig_by_padj)
	n_support_padj = sig_by_padj['neighbor_support'].sum()
	print(f"[5] With neighbor support (of padjusted-significant): {n_support_padj:,}")
	print(f"    --> Filtering by: {filter_by} | neighbor_support in plot={neighbor_support}")

	if filter_by == 'pvalue':
		if len(sig_by_pval) > 0:
			if neighbor_support:
				mask = sig_by_pval['neighbor_support'].values.astype(bool)
				plot_MA_significant(norm_matrix_A, norm_matrix_B, resolution, log2fc_cutoff, upper_limit=upper_limit, lower_limit=lower_limit, xscale='log', sig_rows=significant_rows[mask], sig_cols=significant_cols[mask])
			else:
				plot_MA_significant(norm_matrix_A, norm_matrix_B, resolution, log2fc_cutoff, upper_limit=upper_limit, lower_limit=lower_limit, xscale='log', sig_rows=significant_rows, sig_cols=significant_cols)
		return sig_by_pval,total_tested

	elif filter_by == 'padjusted':
		if len(sig_by_padj) > 0: #because if we don't have significant interactions then we cannot get the neighbor support, and we cannot plot, so we need to check if we have significant interactions before trying to get the neighbor support and plot
			if neighbor_support:
				mask_n = sig_by_padj['neighbor_support'].values.astype(bool)
				plot_MA_significant(norm_matrix_A, norm_matrix_B, resolution, log2fc_cutoff, upper_limit=upper_limit, lower_limit=lower_limit, xscale='log', sig_rows=padj_rows[mask_n], sig_cols=padj_cols[mask_n])
			else:
				plot_MA_significant(norm_matrix_A, norm_matrix_B, resolution, log2fc_cutoff, upper_limit=upper_limit, lower_limit=lower_limit, xscale='log', sig_rows=padj_rows, sig_cols=padj_cols)
		return sig_by_padj,total_tested

def get_differential_interactions(norm_matrix_A: np.ndarray, norm_matrix_B: np.ndarray, resolution: int, permutations_n: int,
						    upper_limit=7e6, lower_limit=5e4, bayesian=True, do_IF_cutoff=True):
	''''Identify differential interactions between two normalized matrices using a permutation-based approach.
	Args:
	norm_matrix_A: normalized contact matrix for sample A
	norm_matrix_B: normalized contact matrix for sample B
	resolution: bin resolution in base pairs
	permutations_n: number of permutations to perform for significance testing
	cut_off: p-value cutoff for significance
	upper_limit: upper distance limit for testing interactions (default 7Mb)
	lower_limit: lower distance limit for testing interactions (default 50kb)
	bayesian: whether to use Bayesian estimation for log2 fold changes (default True)
	Returns:
	significant_rows: array of row indices for significant interactions
	significant_cols: array of column indices for significant interactions
	significant_log2fc: array of log2 fold changes for significant interactions
	pvalues_significant: array of p-values for significant interactions
	distances_significant: array of genomic distances for significant interactions'''
	#### _______ Initial normalized matrices _______

	
	if bayesian:
		#Obs rows and cols are the indices of the bins that have valid counts in both matrices, so we can compute log2fc only for those bins. The distances and log2fc arrays will be filtered to only include those bins.
		distances, log2_fc, distances_out, log2_fc_out, obs_rows, obs_cols, IF_values, IF_values_out = compute_log2fc_bayesian(norm_matrix_A,norm_matrix_B,resolution,lower_limit, upper_limit)
		# For permutation test, use raw log2fc so the comparison is not biased by Bayesian shrinkage
		_, log2_fc_raw, _, log2_fc_raw_out, *_ = compute_log2fc(norm_matrix_A, norm_matrix_B, resolution, upper_limit, lower_limit)
		raw_log2fc_for_perm = np.concatenate([log2_fc_raw, log2_fc_raw_out]) if log2_fc_raw_out.size else log2_fc_raw
	else:
		distances, log2_fc, distances_out, log2_fc_out , obs_rows, obs_cols, IF_values, IF_values_out = compute_log2fc(norm_matrix_A, norm_matrix_B, resolution,lower_limit, upper_limit)

	# Get cutoffs
	#Must be roughly the same because we have normalized
	cutA,cutB = _calculate_cutoffs(norm_matrix_A, norm_matrix_B,
						      resolution, lower_limit, upper_limit)

	IF_cutoff = np.mean([cutA, cutB])
	print("Cutoff: ", IF_cutoff)

	IF_matrix = (norm_matrix_A + norm_matrix_B) / 2
	mask_matrix = (IF_matrix >= IF_cutoff)

	# This is to get the amount of reads to sample, we cannot make the cutoff here because we need to model dispersion in the permutation test
	# These have to use the original counts, because they define the variance of the null distribution (how spread the fc are)
	# Normalization does not change the total number of reads, it only rescales the values.
	# If we use an inflated n, we will get a narrow null distribution, which leaves small room, and the multinomial becomes too precise, everything will be significant.
	# Changed to use original counts
	n_A,n_B=get_reads_count(norm_matrix_A), get_reads_count(norm_matrix_B)

	distances = np.concatenate([distances, distances_out]) if distances_out.size else distances
	observed_log2fc    = np.concatenate([log2_fc,   log2_fc_out])   if log2_fc_out.size   else log2_fc
	
	# Use raw log2fc for permutation test when bayesian=True to keep the comparison symmetric
	obs_log2fc = raw_log2fc_for_perm if bayesian else observed_log2fc

	# Filter all arrays to bins where average cutoff is met, if requested
	if do_IF_cutoff: 
		# whether the specific bin passed the cutoff
		if_mask = mask_matrix[obs_rows, obs_cols]
		# Keep only those bins
		obs_rows_f       = obs_rows[if_mask]
		obs_cols_f       = obs_cols[if_mask]
		distances      = distances[if_mask]
		observed_log2fc = observed_log2fc[if_mask]
		obs_log2fc    = obs_log2fc[if_mask]

	# _______ Create background matrix for permutation testing _______
	# ______ Keep all info because we are modeling dispersion in the background too
	
	obs_A = np.asarray(norm_matrix_A,dtype=np.float32)
	obs_B = np.asarray(norm_matrix_B,dtype=np.float32)
	background = obs_A + obs_B
	# We keep all bins because we need to model the dispersion
	background = np.asarray(background, dtype=np.float64)
	background = background[obs_rows, obs_cols] #


	# The trick is here, probabilities are always fixed, because we divided them, so it does not matter the initial number of reads 
	# Bin_i will always have the same probability of being sampled, and the multinomial will take care of the variance, because it is a random sampling process
	probs = background / background.sum()

	# Same size as the initial log fold change array
	# The size changes because we will filter inside the loop to only keep bins with mask
	exceed_counts = np.zeros(len(observed_log2fc), dtype=np.int32)

	for p in range(permutations_n):
		print("Subsample: ", p)
		# Downsample A and B from the background matrix using multinomial distribution
		counts_A = np.random.multinomial(n_A, probs)
		counts_B = np.random.multinomial(n_B, probs)

		# Here is where we apply the IF cutoff to the counts, so that we only keep bins that pass the cutoff for the permutation test. This ensures that we are comparing apples to apples, and not including bins that are too sparse to be reliable.
		# Because we need to compare the same bins that are in the observed log2fc, we need to apply the same mask to the counts. This is important because we want to make sure that we are comparing the same set of bins in both the observed and permuted data.
		counts_A = counts_A[if_mask] if do_IF_cutoff else counts_A
		counts_B = counts_B[if_mask] if do_IF_cutoff else counts_B
		# I don't know if we should also do the bayesian estimation here, or just compute the log2fc directly from the counts. For now, I will just compute the log2fc directly from the counts, but we can change it later if we want to be more consistent with the observed log2fc.
		#if a pvalue is small after the bayesian estimation, then comparing against the normal values will 
		# always be small vs 'bigger'
		# Lets say i got 0.5 after bayesian
		# then i do all logs on the subsampling tests, most likely these are gonna be bigger than the bayesian one
		# so our numerator will be close to the number of permutations
		# and the pvalue of those small logfc will be close to 1
		# Which makes more sense because if the log2fc was pushed to zero by the bayesian stimation
		# then it should come from sparse counts
		# therefore it is not confident calling, then it makes sense that the log2fc values from the subsampling will be bigger than the observed one, 
		# and therefore the pvalue will be close to 1, which is what we expect for low confidence calls
		log2_fc_perm = np.log2((counts_A ) / (counts_B ))

		exceed_counts += (np.abs(log2_fc_perm) > np.abs(obs_log2fc)).astype(np.int32)

	# Go back to coordinates of the bins
	p_values = exceed_counts / permutations_n

	if do_IF_cutoff:
		obs_rows = obs_rows_f
		obs_cols = obs_cols_f
	
	return obs_rows, obs_cols, observed_log2fc, p_values, distances

def adjust_pvalues_by_distance(pvalues, distances, method='fdr_bh'):
	'''Adjust p-values for multiple testing while accounting for genomic distance.
	Args:
	pvalues: array of p-values to adjust
	distances: array of genomic distances corresponding to each p-value
	method: method for multiple testing correction (default 'fdr_bh' for Benjamini-Hochberg)
	Returns:
	adjusted_pvalues: array of adjusted p-values'''

	unique_distances = np.unique(distances)
	adjusted_pvalues = np.ones_like(pvalues)
	for d in unique_distances:
		# The key is generating this mask that will help us retrieve each entry 
		# corresponding to a given distance
		mask = distances == d
		# So we can get pvalues that correspond to a given distance with the mask
		pvalues_d = pvalues[mask]
		if len(pvalues_d) > 0:
			# And apply the multiple testing correction method to the p-values of that distance
			_, adj_pvalues_d, _, _ = multipletests(pvalues_d, method=method)
			# We can just populate the already defined array of adjusted pvalues with the mask
			adjusted_pvalues[mask] = adj_pvalues_d

	return adjusted_pvalues

def adjust_pvalues_global(pvalues, method='fdr_bh'):
	'''Adjust p-values for multiple testing across all interactions without accounting for distance.
	Args:
	pvalues: array of p-values to adjust
	method: method for multiple testing correction (default 'fdr_bh' for Benjamini-Hochberg)
	Returns:
	adjusted_pvalues: array of adjusted p-values'''
	_, adjusted_pvalues, _, _ = multipletests(pvalues, method=method)
	return adjusted_pvalues

def filter_significant_interactions(rows, cols, log2fc, pvalues, adjusted_pvalues, distances, pvalue_cutoff=0.05, log2fc_cutoff=1):
	'''Filter interactions based on  p-value and log2 fold change cutoffs.
	Args:
	rows: array of row indices for interactions
	cols: array of column indices for interactions
	log2fc: array of log2 fold changes for interactions
	pvalues: array of raw p-values for interactions
	adjusted_pvalues: array of adjusted p-values for interactions
	distances: array of genomic distances for interactions
	pvalue_cutoff: significance cutoff for adjusted p-values (default 0.05)
	log2fc_cutoff: minimum absolute log2 fold change for significance (default 1)
	Returns:
	significant_rows: array of row indices for significant interactions based on pvalue cutoff and log2fc cutoff
	significant_cols: array of column indices for significant interactions based on pvalue cutoff and log2fc cutoff
	significant_log2fc: array of log2 fold changes for significant interactions based on pvalue cutoff and log2fc cutoff
	significant_pvalues: array of adjusted p-values for significant interactions
	significant_distances: array of genomic distances for significant interactions'''
	significant_mask = (pvalues <= pvalue_cutoff) & (np.abs(log2fc) >= log2fc_cutoff)
	return rows[significant_mask], cols[significant_mask], log2fc[significant_mask], pvalues[significant_mask], adjusted_pvalues[significant_mask],	 distances[significant_mask]

def plot_MA_significant(matrix_A, matrix_B, resolution, log2_fc_cutoff, upper_limit=7e6, lower_limit=5e4, xscale='log',sig_rows=None, sig_cols=None):
    distances, log2_fc, distances_out, log2_fc_out, obs_rows, obs_cols, IF_values, IF_values_out = compute_log2fc(matrix_A, matrix_B, resolution, upper_limit, lower_limit)

    all_distances = np.concatenate([distances, distances_out]) if distances_out.size else distances
    all_log2fc    = np.concatenate([log2_fc,   log2_fc_out])   if log2_fc_out.size   else log2_fc
    all_rows      = np.concatenate([obs_rows[:len(distances)],  obs_rows[len(distances):]]) if distances_out.size else obs_rows
    all_cols      = np.concatenate([obs_cols[:len(distances)],  obs_cols[len(distances):]]) if distances_out.size else obs_cols

    fig = plt.figure(figsize=(8, 6), dpi=200)
    if distances_out.size:
        plt.scatter(distances_out, log2_fc_out, alpha=0.3, s=1, color='lightgray', label='Outside range')
    plt.scatter(distances, log2_fc, alpha=0.5, s=1, color='black', label='Within range')

    if sig_rows is not None and sig_cols is not None:
        # Build a set of significant (row, col) pairs for fast lookup
        sig_set = set(zip(sig_rows, sig_cols))
        sig_mask = np.array([( r, c) in sig_set for r, c in zip(all_rows, all_cols)])
        plt.scatter(all_distances[sig_mask], all_log2fc[sig_mask],
                    alpha=0.8, s=3, color='red', label='Significant', zorder=5)

    plt.hlines(0,  xmin=all_distances.min(), xmax=all_distances.max(), colors='red',  linestyles='dashed')
    plt.hlines( log2_fc_cutoff, xmin=all_distances.min(), xmax=all_distances.max(), colors='gray', linestyles='dashed')
    plt.hlines(-log2_fc_cutoff, xmin=all_distances.min(), xmax=all_distances.max(), colors='gray', linestyles='dashed')
    plt.vlines(upper_limit, ymin=all_log2fc.min(), ymax=all_log2fc.max(), colors='red', linestyles='dotted', linewidth=0.5)
    if lower_limit > 0:
        plt.vlines(lower_limit, ymin=all_log2fc.min(), ymax=all_log2fc.max(), colors='red', linestyles='dotted', linewidth=0.5)
    plt.xlabel('Genomic Distance (bp)')
    if xscale == 'log':
        plt.xscale('log')
    plt.ylabel('Log2 Fold Change (Sample A / Sample B)')
    plt.legend(markerscale=5)
    return fig



# ── Neighbor Support ───────────────────────────────────────────────────────
def get_neighbors(i, j):
    
    """Return direct (4-connected) neighbors of a matrix element.
	d is diagonal of the spike in
	row, column,value"""
    neighbors = []
    n_rows, n_cols = 1000000000, 1000000000

    if i + 1 < n_rows: #ensures we don't go out of bounds on the rows (lower bound)
        neighbors.append((i + 1, j))
    if i - 1 >= 0: #ensures we don't go out of bounds on the rows (upper bound)
        neighbors.append((i - 1, j))
    if j + 1 < n_cols: #ensures we don't go out of bounds on the columns (right bound)
        neighbors.append((i, j + 1))
    if j - 1 >= 0: #ensures we don't go out of bounds on the columns (left bound)
        neighbors.append((i, j - 1))

    return neighbors

def neighbor_support_f(final_results):
	sig_interactions = set(
		zip(final_results['row'], final_results['col'])
	)

	k = 1

	neighbor_support = []

	if k > 0:
		for a, b in zip(final_results['row'], final_results['col']):
			neighbors = get_neighbors_coordinates([(a, b)], K=1) #The logic is that even the outer neighbors will have direct neighbors that can support the interaction, if not then is most likely noise
			has_support = False
			for degree in range(k): #degree 0 is first degree neighbors
				for neighbor in neighbors[degree]:
					if (neighbor[0], neighbor[1]) in sig_interactions:
						has_support = True
						break
				if has_support:
					break

			neighbor_support.append(has_support)

	else:
		neighbor_support = [False] * len(final_results)
		
	final_results['neighbor_support'] = neighbor_support
	return final_results

def get_neighbors_coordinates(center, K=2, verbose=False):
    """
    center: list of (i, j) coordinates to start from
    K: number of degrees of neighbors to expand
    General K-degree neighbor expansion.
    Returns:
        A list of lists:
            neighbors_per_degree[d] = list of (i, j) for degree d+1
    """
    seen = set()
    neighbors_per_degree = [[] for _ in range(K)]


    frontier = []  # list of (i, j)

    for (i0, j0) in center:
        seen.add((i0, j0))
        first_neighbors = get_neighbors(i0, j0)

        for i, j in first_neighbors:
            if (i, j) not in seen:
                frontier.append((i, j))
                neighbors_per_degree[0].append((i, j))
                seen.add((i, j))

    if verbose:
        print(f"Degree 1 neighbors: {len(neighbors_per_degree[0])}")

    for degree in range(1, K):
        next_frontier = []
        for (i, j) in frontier:
            new_neighbors = get_neighbors(i, j)
            for ni, nj in new_neighbors:
                if (ni, nj) not in seen:
                    neighbors_per_degree[degree].append((ni, nj))
                    next_frontier.append((ni, nj))
                    seen.add((ni, nj))

        if verbose:
            print(f"Degree {degree+1} neighbors: {len(neighbors_per_degree[degree])}")

        frontier = next_frontier

        # No more expansion possible
        if len(frontier) == 0:
            break

    return (*neighbors_per_degree,)


# ───── Bayesian approach ──────────────────────────────────────────────────────
def bayesian_stimation(counts_A, counts_B):
		theta_grid = np.linspace(-6, 6, counts_A.shape[0]) # grid of possible true log2 fold changes
		#For each theta in the grid compute a prior probability based on a standard normal distribution
		prior = np.exp(-theta_grid**2 / 2) # standard normal prior
		prior /= np.trapezoid(prior, theta_grid) # Normalize the prior over the grid

		#transform all thetas into probabilities
		p_theta = 2**theta_grid / (1 + 2**theta_grid)
		#[:,None] adds a new axis to make it a rows vector, 
		# so that we can broadcast it with p_theta which is a column vector
		logL = (counts_A[:,None] * np.log(p_theta[None, :]) + counts_B[:,None] * np.log(1 - p_theta[None, :]))
		#We normalize the log likelihoods by subtracting the maximum log likelihood for each entry, to prevent numerical underflow when we exponentiate later. This does not change the relative probabilities, just rescales them.
		logL -= logL.max(axis=1, keepdims=True)      
		L = np.exp(logL)
		L /= np.trapezoid(L, theta_grid, axis=1)[:, None] # Normalize the likelihoods over the grid for each entry

		posterior = L * prior[None, :]
		posterior /= np.trapezoid(posterior, theta_grid, axis=1)[:, None] # Normalize the posterior over the grid for each entry
		post_mean = np.trapezoid(theta_grid[None, :] * posterior, theta_grid, axis=1)
		return post_mean


def compute_IF_cutoff(if_in: np.ndarray, if_out: np.ndarray, option: str = "out"):
	"""Compute the cutoff for interaction frequencies (IF) based on the specified option."""
	if option == "all":
		all_if = np.concatenate([if_in, if_out]) if if_out.size else if_in
		cutoff = np.percentile(all_if, 15) if all_if.size else None
	elif option == "out":
		cutoff = np.percentile(if_out, 15) if if_out.size else None
	else:
		raise ValueError("Invalid option. Please choose 'all' or 'out'.")

	return cutoff

def _generate_IF_from_matrices(matrix_A, matrix_B, resolution,upper_limit=7e6, lower_limit=5e4):
	n = matrix_A.shape[0]
	min_d= int(max(1, lower_limit // resolution))
	max_d= int(min(n, upper_limit // resolution))

	IF_list, IF_out_list = [], []

	for d in range(1,n):
		diag_A=np.asarray(matrix_A.diagonal(d),dtype=np.float32).ravel()
		diag_B=np.asarray(matrix_B.diagonal(d),dtype=np.float32).ravel()
		
		valid = (diag_A > 0) & (diag_B > 0)
		valid = valid.ravel()
		if not np.any(valid):
			continue
		# Compute mean interaction frequency for valid entries
		mean_IF = (diag_A[valid] + diag_B[valid]) / 2

		if min_d <= d < max_d:
			IF_list.append(mean_IF)
		else:
			IF_out_list.append(mean_IF)

	IF_values = np.concatenate(IF_list) if IF_list else np.array([])
	out_IF_values = np.concatenate(IF_out_list) if IF_out_list else np.array([])

	return IF_values, out_IF_values

def compute_log2fc_bayesian(matrix_A,matrix_B,resolution,lower_limit=5e4, upper_limit=7e6):
	#The normalized matrices still include all bins, we just removed them to avoid normalizing on dispersed values 
	n = matrix_A.shape[0]
	min_d= int(max(1, lower_limit // resolution))
	max_d= int(min(n, upper_limit // resolution))

    
	distances_list, log2fc_list, rows_in_list, cols_in_list, IF_list= [], [], [], [], []
	distances_out_list, log2fc_out_list, rows_out_list, cols_out_list, out_IF_list = [], [], [], [], []

	for d in range(1,n):
		diag_A=np.asarray(matrix_A.diagonal(d),dtype=np.float32).ravel()
		diag_B=np.asarray(matrix_B.diagonal(d),dtype=np.float32).ravel()
		
		valid = (diag_A > 0) & (diag_B > 0)
		valid = valid.ravel()
		rows_valid = np.where(valid)[0] #because it is being treated as a 1D array after raveling
		cols_valid = rows_valid + d
		if not np.any(valid):
			continue
		# Compute log2 fold change for valid entries
		#For each valid entry we have the observed counts in A and B, each element correspond to the respective bin
		#shape (N,) where N is the number of valid entries in the diagonal
		diag_A_valid = (diag_A[valid]).astype(np.float32) # add a small pseudocount to avoid log2(0)
		diag_B_valid = (diag_B[valid]).astype(np.float32) # add a small pseudocount to avoid log2(0)
		
		#Now this is the log2fc
		post_mean = bayesian_stimation(diag_A_valid, diag_B_valid)
		log2fc = post_mean
		mean_IF = (diag_A_valid + diag_B_valid) / 2
		d_bp = d * resolution

		if min_d <= d < max_d:
			log2fc_list.append(log2fc)
			distances_list.append(np.full(log2fc.size, d_bp, dtype=np.float32))
			rows_in_list.append(rows_valid) #keep track of the row indices of the valid entries in the diagonal
			cols_in_list.append(cols_valid) #keep track of the column indices of the valid entries in the diagonal
			IF_list.append(np.full(log2fc.size, mean_IF, dtype=np.float32))
		else:
			log2fc_out_list.append(log2fc)
			distances_out_list.append(np.full(log2fc.size, d_bp, dtype=np.float32))
			rows_out_list.append(rows_valid) #keep track of the row indices of the valid entries in the diagonal
			cols_out_list.append(cols_valid) #keep track of the column indices of the valid entries in the diagonal
			out_IF_list.append(np.full(log2fc.size, mean_IF, dtype=np.float32))
		
	obs_rows = np.concatenate(rows_in_list + rows_out_list) if (rows_in_list or rows_out_list) else np.array([], dtype=int)
	obs_cols = np.concatenate(cols_in_list + cols_out_list) if (cols_in_list or cols_out_list) else np.array([], dtype=int)

	distances = np.concatenate(distances_list) if distances_list else np.array([])
	log2fc_values = np.concatenate(log2fc_list) if log2fc_list else np.array([])
	distances_out     = np.concatenate(distances_out_list)     if distances_out_list     else np.array([])
	log2fc_values_out = np.concatenate(log2fc_out_list)        if log2fc_out_list        else np.array([])
	
	IF_values = np.concatenate(IF_list) if IF_list else np.array([])
	out_IF_values = np.concatenate(out_IF_list) if out_IF_list else np.array([])

	return distances, log2fc_values, distances_out, log2fc_values_out, obs_rows, obs_cols, IF_values, out_IF_values



# ------- Plot decays
def _plot_distance_decay(matrix_A: np.ndarray, matrix_B: np.ndarray, resolution: int, type_plot: str = "log", filter: bool = True):
	# It should then just plot the decay based on the filtered values because we already filtered inside get_distances_and_counts_paired
	distances_A, counts_A, counts_B = get_distances_and_counts_paired(matrix_A, matrix_B, resolution, filter=filter)
	distances_B = distances_A  # Since we're dealing with paired matrices, the distances should be the same

	A_values=calculate_distance_decay_individual(distances_A, counts_A, metric='median')
	B_values=calculate_distance_decay_individual(distances_B, counts_B, metric='median')
	#Arrays here have elements per distance, not per bin
	
	if type_plot == "linear":
		plt.figure(figsize=(6, 4), dpi=200)
		plt.plot(A_values[0], A_values[2], label='Matrix A', color='blue')
		plt.plot(B_values[0], B_values[2], label='Matrix B', color='orange')
		plt.xlabel('Distance (bp)')
		plt.ylabel('Distance Decay (median)')
		plt.title('Distance Decay Comparison')
		plt.legend()
		plt.show()
	elif type_plot == "log":
		plt.figure(figsize=(6, 4), dpi=200)
		plt.plot(A_values[5], A_values[6], label='Matrix A', color='blue')
		plt.plot(B_values[5], B_values[6], label='Matrix B', color='orange')
		plt.xlabel('Distance (bp)')
		plt.ylabel('Distance Decay (median)')
		plt.title('Distance Decay Comparison')
		plt.legend()
		plt.show()
	else:
		raise ValueError("Invalid type_plot. Choose 'linear' or 'log'.")
	