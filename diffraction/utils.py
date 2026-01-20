import os
import hicstraw 
import cooler

import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from statistics import mean, median
from scipy.stats import linregress
from scipy.sparse import coo_matrix
import math
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

#### DifFracTion functions ####

permutations_n = 1000

## Main function to normalize
def alpha_normalization(matrixA,matrixB,chromosome_length,resolution,desired_alpha=-1.08,upper_limit=7e6,lower_limit=0,metric='median'):
	'''
	Docstring for alpha_normalization
	
	:param matrixA: nxn numpy array representing the first contact matrix
	:param matrixB: nxn numpy array representing the second contact matrix
	:param chromosome_length: Length of the chromosome in base pairs
	:param resolution: Bin resolution in base pairs
	:param desired_alpha: Desired alpha value for normalization
	:param upper_limit: Upper distance limit for normalization, default is 7e6
	:param lower_limit: Lower distance limit for normalization, default is 0
	:param metric: Metric to use for distance decay calculation (e.g., 'median' or 'mean')
	'''

	matrixA_counts,matrixB_counts=get_counts_by_distance_both_densematrices(matrixA,matrixB,chromosome_length,resolution)
	d_1, y_1, m_1, b_1, alpha_1 = calculate_distance_decay(matrixA_counts, metric=metric,
                              lower_limit=lower_limit,upper_limit=upper_limit)
	d_2, y_2, m_2, b_2, alpha_2 = calculate_distance_decay(matrixB_counts, metric=metric,
                              lower_limit=lower_limit,upper_limit=upper_limit)
     																

	scaled_counts_1,scaled_counts_2,full_factors_r,full_factors_s=pre_counts(matrixA_counts,matrixB_counts,d_1,d_2,y_1,y_2,desired_alpha,resolution)
     
	#### Apply factors to dense matrices
	matrixA_adj=apply_dict_factors_to_dense(matrixA,full_factors_r,resolution)
	matrixB_adj=apply_dict_factors_to_dense(matrixB,full_factors_s,resolution)
	return matrixA_adj,matrixB_adj

def iterative_normalization(matrixA,matrixB,chromosome_length,resolution,upper_limit=7e6,lower_limit=0,metric='median', verbose=True):
	matrixA_counts,matrixB_counts=get_counts_by_distance_both_densematrices(matrixA,matrixB,chromosome_length,resolution)
	matrixB_counts_scaled, best_error, updated_factor = cycle_correction(matrixA_counts, matrixB_counts,metric=metric,
											error_type='mse',upper_limit=upper_limit,lower_limit=lower_limit)
	matrixA_counts_scaled=matrixA_counts
	#### Apply factors to dense matrices
	matrixA_adj=apply_factors_to_dense(matrixA,1,resolution)
	matrixB_adj=apply_factors_to_dense(matrixB,updated_factor,resolution)
	if verbose:
		print(f"Final scaling factor applied to matrix B: {updated_factor:.4f}")
		print(f"Final error after iterative normalization: {best_error:.4f}")
	return matrixA_adj,matrixB_adj








#
### Individual fitting
def get_individual_factors(y, expected_fit):
    """
    Iteratively bring y closer to expected_fit without collapsing into a straight line.
    Returns: the new y values after the correction with their respective factors
    """
    y = np.asarray(y, dtype=float)
    ef = np.asarray(expected_fit, dtype=float)
    factors = ef / y
    average_factor = np.mean(factors)
    factors = np.append(factors,average_factor)
    return factors

def apply_factors_to_dict_counts(counts,factors,resolution):
	upper_limit = 7000000
	"""
	Apply factors to counts
	Args:
		counts: dict {distance_in_bp: list_of_counts}
		factors: list of factors per distance and overall average as the last element
		resolution: bin resolution in base pairs
	Returns:
		new_counts: dict {distance_in_bp: list_of_rescaled_counts}
	"""
	
	n_distances = len(counts)

	factors_dict={}

	keys_for_fits = [resolution*(k) if resolution*(k) > 0 else None for k in range(0, (upper_limit//resolution)+1)]
	keys_for_fits.remove(None)
	
	for k,f in zip(keys_for_fits, factors):
		factors_dict[k]=f
	for k in counts.keys():
		if k not in factors_dict:
			factors_dict[k]=factors[-1]

	new_counts = defaultdict(list)

	for k,d in counts.items():
		counts_rescaled=(np.array(d, dtype=float) * factors_dict[k]).tolist()
		new_counts[k] = counts_rescaled
	
	return new_counts,factors_dict

def apply_factors_to_dict_counts_v1(counts,factors,resolution):
	"""
	Apply factors to counts
	Args:
		counts: dict {distance_in_bp: list_of_counts}
		factors: list of factors per distance and overall average as the last element
		resolution: bin resolution in base pairs
	Returns:
		new_counts: dict {distance_in_bp: list_of_rescaled_counts}
	"""
	keys = list(counts.keys())
	values = list(counts.values())
	n_distances = len(keys)
	per_distance = factors[:-1]
	overall_average = factors[-1]

	full_factors = np.empty(n_distances)
	full_factors[:len(per_distance)] = per_distance
	full_factors[len(per_distance):] = overall_average
 
	# New version of inferring the unit
	support_bp= 0
	support_index = 0
	for i in range(n_distances):
		basis = keys[i] % resolution
		if basis == 0:
			support_bp += 1
		elif basis != 0:
			support_index += 1
	if support_bp > support_index:
		unit = 'bp'
	elif support_index > support_bp:
		unit = 'index' # it will never be this case because we force the dict to have distances in bp
	else:
		unit = 'unknown'
		
	new_counts = defaultdict(list)
	if unit == 'bp':
		for k,c in counts.items():
			index = (k // resolution)-1
			factor = full_factors[index]
			counts_rescaled = (np.array(c, dtype=float) * factor).tolist()
			new_counts[k] = counts_rescaled

	elif unit == 'index':
		for k,c in counts.items():
			index = k
			factor = full_factors[index]
			counts_rescaled = (np.array(c, dtype=float) * factor).tolist()
			new_counts[k] = counts_rescaled
	else:
		raise ValueError("Could not infer bin unit (bp or index).")

	return new_counts,full_factors

### Data Preparation 
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

def infer_bins_unit(records, resolution, chromosome_length, sample_n=1000):
    import math
    """Infer whether the binX values in records are in base pairs (bp) or bin indices.
    Parameters:
    records (iterable): An iterable of records with a 'binX' attribute.
    resolution (int): The resolution in base pairs.
    chromosome_length (int): The length of the chromosome in base pairs.
    sample_n (int): Number of records to sample for inference (default is 1000)."""
    
    support_bp= 0
    support_index = 0
    xs = [r.binX for r in records[:sample_n]] if hasattr(records, '__getitem__') else []
    if not xs:
        for k, r in zip(range(sample_n), records):
            xs.append(r.binX)
    for x in xs:
        basis = x % resolution

        if basis == 0:
            support_bp += 1
        elif basis != 0:
            support_index += 1
    if support_bp > support_index:
        return 'bp'
    elif support_index > support_bp:
        return 'index'	
    else:
    	return 'unknown'
 
def get_counts_by_distance_from_matrixZoom(matrix,chromosome_length,resolution):
	'''Get contact counts grouped by genomic distance from a MatrixZoomData object.
	Args: 
		matrix: MatrixZoomData object from hicstraw
		resolution: resolution in base pairs
		chromosome_length: length of the chromosome in base pairs
	Returns:
		counts: dict with distances as keys and list of counts as values
	'''
	recs = matrix.getRecords(0, chromosome_length, 0, chromosome_length) 
	unit = infer_bins_unit(recs, resolution, chromosome_length)
	counts = defaultdict(list)
	for r in recs:
		d_bins = abs(r.binX - r.binY)
		if d_bins == 0:
			continue  # Skip diagonal
		if unit == 'bp':
			d_bp = d_bins
		elif unit == 'index':
			d_bp = d_bins * resolution
		else:
			raise ValueError("Could not infer bin unit (bp or index).")
		if r.counts > 0:
			counts[d_bp].append(r.counts)
	return counts  

def get_counts_by_distance_from_dense(matrix, resolution):
    """
    Get contact counts grouped by genomic distance from a dense Hi-C matrix.

    Args:
        matrix: 2D numpy array (square contact matrix)
        resolution: bin resolution in base pairs

    Returns:
        counts: dict {distance_in_bp: list_of_counts}
    """
    n_bins = matrix.shape[0]
    counts = defaultdict(list)
	# It will always be index because we are passing a dense matrix and we are iterating over n_bins range
    for i in range(n_bins):
        for j in range(i+1, n_bins):   # skip diagonal, only upper triangle WE NEVER INCLUDE 0
            d_bins = j - i

            d_bp = d_bins * resolution
            val = matrix[i, j]
            if val > 0:
                counts[d_bp].append(val)

    return counts

def get_shared_distances(counts_1,counts_2):
	'''Get shared distances between two counts dictionaries.
	Args:
		counts_1: dict with distances as keys and list of counts as values
		counts_2: dict with distances as keys and list of counts as values'''
	shared_keys = set(counts_1.keys()).intersection(set(counts_2.keys()))
	counts_1 = {k: counts_1[k] for k in shared_keys}
	counts_2 = {k: counts_2[k] for k in shared_keys}
	counts_1 = dict(sorted(counts_1.items()))
	counts_2 = dict(sorted(counts_2.items()))
	return counts_1, counts_2

def get_counts_by_distance_both_matricesZoom(matrix_1,matrix_2,chromosome_length,resolution):
	counts_1 = get_counts_by_distance_from_matrixZoom(matrix_1,chromosome_length,resolution)
	counts_2 = get_counts_by_distance_from_matrixZoom(matrix_2,chromosome_length,resolution)
	counts_1, counts_2 = get_shared_distances(counts_1,counts_2)
	return counts_1, counts_2

def get_counts_by_distance_both_densematrices(matrix_1,matrix_2,chromosome_length,resolution):
	counts_1 = get_counts_by_distance_from_dense(matrix_1,resolution)
	counts_2 = get_counts_by_distance_from_dense(matrix_2,resolution)
	counts_1, counts_2 = get_shared_distances(counts_1,counts_2)
	return counts_1, counts_2

def rawMatrix(matrix,chromosome_length):
	'''Get the raw contact matrix from a MatrixZoomData object.
	Args:
		matrix: MatrixZoomData object from hicstraw
		chromosome_length: length of the chromosome in base pairs
	Returns:
		raw_matrix: contact matrix as a 2D numpy array getRecordsAsMatrix
	'''
	raw_matrix = matrix.getRecordsAsMatrix(0,chromosome_length,0,chromosome_length)
	return raw_matrix


### Prepare matrix to remove distances outside limits
def filter_matrix_by_distance(matrix,resolution,upper_limit,lower_limit):
    '''
    Set to zero the contacts in the matrix that are outside the specified distance limits.
    Args:
	   matrix: 2D numpy array (square contact matrix)
	   resolution: bin resolution in base pairs
	   upper_limit: upper distance limit in base pairs
	   lower_limit: lower distance limit in base pairs
	Returns:
		filtered_matrix: 2D numpy array with contacts outside the distance limits set to zero
    '''
    n_bins = matrix.shape[0]
    filtered_matrix = matrix.copy()
    for i in range(n_bins):
        for j in range(n_bins):
            d_bins = abs(j - i)
            d_bp = d_bins * resolution
            if d_bp > upper_limit or d_bp < lower_limit:
                filtered_matrix[i, j] = 0
    return filtered_matrix

### Decay Fitting
def extract_metric_decay(counts,metric='mean'):
	'''Get decay data (mean per distance) from counts dictionary.
	Args:
 		counts: dict with distances as keys and list of counts as values
		metric: 'mean' or 'median' (default is 'mean')	
	Returns:
		decay_dict: dict with distances as keys and mean/median counts as values (one value per key)
    '''
	if metric == 'median':
		decay_dict = {d: median(c) for d, c in counts.items()}
	elif metric == 'mean':
		decay_dict = {d: mean(c) for d, c in counts.items()}
	return decay_dict

def calculate_distance_decay(counts,metric='mean',upper_limit=7e6,lower_limit=5e5,verbose=False):
	'''Calculate distance decay and fit a power-law model.
	Args:
		counts: dict with distances as keys and list of counts as values
		metric: 'mean' or 'median' (default is 'mean')
		upper_limit: upper distance limit for fitting (default is 7e6)
		lower_limit: lower distance limit for fitting (default is 5e5)
		verbose: if True, print alpha and intercept (default is False) 
	Returns:
		d: list of distances used in fitting
		c: list of mean/median counts used in fitting
		m: slope of the fitted line in log-log space
		b: intercept of the fitted line in log-log space
		alpha: exponent of the power-law decay (alpha = -m)
    '''
	
	decay_dict = extract_metric_decay(counts,metric)

	#Filter by distance, only the distances that are involved in the fitting
	decay_filtered={d: c for d,c in decay_dict.items() if d <= upper_limit and d >= lower_limit}
	d = list(decay_filtered.keys())
	c = list(decay_filtered.values())
	#Flag
	logx, logy = np.log(d), np.log(c)
	regression = linregress(logx, logy) #log-log regression
	m,b = regression.slope, regression.intercept
	alpha = -m
	if verbose:
		print(f"α = {alpha:.2f}, intercept = {b:.2f}")
	return d,c,m,b,alpha

def refit_intercept(y, d, m):
	
	"""
	We need to refit the b because the previous one was based on a given slope
	Fit intercept b given a fixed slope m.
	y: observed counts (log10 scale)
	d: distances (log10 scale)
	"""
	#Flag
	logy = np.log(y)
	logd = np.log(d)
	b = np.mean(logy - m*logd)
	return b

def linear_distance_decay(b, m, distances):
	'''Based on the fitted parameters, get the expected counts at given distances.
	Args:
		b: intercept of the fitted line in log-log space
		m: slope of the fitted line in log-log space
		distances: list of distances to calculate expected counts for
	Returns:
		expected_counts: list of expected counts at the given distances
	'''
	#Flag 
	d = np.asarray(distances, dtype=float)
	linear_fit = np.exp(b + m*np.log(d))
	return linear_fit

### Correction
def cycle_correction(ref_counts, scaled_counts,metric='median',error_type='mse',upper_limit=7e6,lower_limit=0,verbose=False):
	best_error    = float("inf")
	max_cycles=20
	tol=1e-6
	w_factor=1.08
	updated_factor=1.0
	for cycle in range(max_cycles):
		if verbose:
			print(f"\nCycle {cycle+1}")
		
		improved = False
		set_usable_distances = {d: c for d,c in ref_counts.items() if d <= upper_limit } # Just the set of distances that are actually used for decay
		distances=list(range(1,len(set_usable_distances.keys())))
		for i in distances:
			trial_scaled, error,factor = correction(ref_counts, scaled_counts, i-1,upper_limit=upper_limit,
											lower_limit=lower_limit,metric=metric,error_type=error_type,
            									weight_factor=-1*w_factor,plot=False)
			if error + tol < best_error:
				scaled_counts = trial_scaled
				best_error = error
				updated_factor *= factor  # accumulate factor
				if verbose:
					print(f"    Error: {error:.4f}")
					print("     New best error, updating counts")
				improved = True
			else:
				if verbose:
					print("")			
		#If at least one of the distances make an improvement, then improved=True
		if not improved: #if none of the distances made an improvement then we stop
			break
	return scaled_counts, best_error, updated_factor

def correction(ref_counts, counts_to_scale, i,
                         lower_limit=5e5,upper_limit=7e6,
                         metric='mean',error_type='mse',weight_factor=-1.08,
                         plot=False,verbose=False):
	do_weights = True
 
	# Decays. Metric has to be consistent with the one used on previous steps
	decay_ref = extract_metric_decay(ref_counts, metric=metric)
	decay_2   = extract_metric_decay(counts_to_scale, metric=metric)

	# Values at distances
	y_1 = list(decay_ref.values())
	y_2 = list(decay_2.values())
    
	distances = np.array(list(decay_ref.keys()))
	ratios = np.array(y_1) / np.array(y_2)
    

    
    # The logic here is to scale the counts at distance i by a factor 
    # that is a weighted average of the ratios at all distances,
    # with weights decreasing with distance

	if do_weights:
		weights = distances ** weight_factor # weights decrease with distance
 
	# normalize weights
	# normalizing by the maximum just rescales the whole array 
  	# to be between 0 and 1, ratios are intact
     
		weights = weights / weights.max() # normalizing one weight at the time,
		raw_factor = ratios[i] * weights[i]
  
	else:
      
		raw_factor = ratios[i]

	# Smooth the factor a bit, 0.2 means 20% of the way from 1 to raw_factor
	# so if raw_factor is 2, the final factor will be 1.2
	#Think of η as how aggressively you steer a car back into the lane:	
	# (raw_factor - 1)is the correction away from 1, if raw factor is 2 then 2 - 1 = 1
	# eta * (raw_factor - 1) is the fraction of that correction we apply
	
	eta = 0.15
	factor = 1 + eta * (raw_factor - 1)

	if verbose:
		print("Scale factor at distance index", i, ":", factor)
	scaled_counts = defaultdict(list)
	
	#Populate the scaled counts, all distances values are scaled by the same factor
	#We apply the same factor because the HiC Matrix is a whole polymer and changing one distance only is not possible
	for d, c in counts_to_scale.items(): 
		arr = np.array(c, dtype=float)
		scaled_counts[d] = (arr * factor).tolist()
	
	#0 to upper limit fitting
	d_r, y_r, m_r, b_r, alpha_r = calculate_distance_decay(ref_counts,lower_limit=lower_limit,upper_limit=upper_limit)
	d_s, y_s, m_s, b_s, alpha_s = calculate_distance_decay(scaled_counts,lower_limit=lower_limit,upper_limit=upper_limit)

	#Each iteration
	if plot:
		decay_boxPlot(d_r, ref_counts, m_r, b_r,
					d_s, scaled_counts, m_s, b_s)
     
	decay_ref_fit    = linear_distance_decay(b_r, m_r, d_r)
	decay_scaled_fit = linear_distance_decay(b_s, m_s, d_s)

	error = measure_error_by_distance(decay_scaled_fit, decay_ref_fit, method=error_type)
	
	return scaled_counts, error,factor

def measure_error_by_distance(y_obs,y_exp,method='mse'):
	#Error is measured after fitting the distance decay
	# Now we want to decide which correction wins
	# if we weight shorter distances more then:
	     # when short distances get closer the weighter error will drop even if the large distances diverge
	     # When only large distances improve, the weighted error barely moves 
    
	y_obs=np.asarray(y_obs, dtype=float)
	y_exp=np.asarray(y_exp, dtype=float)
	mask = (y_obs > 0) & (y_exp > 0)
    
	# They are obviously in order thats why we can use the relation
	distances = np.arange(1,len(y_obs)+1)    
	distances = np.asarray(distances)[mask]
	if method == 'mse':
		se=((np.log(y_obs[mask])-np.log(y_exp[mask]))**2)
		weights = distances ** -1.08
		# weighted average across all distances, all weights sum to 1
		# if we were to normalize by max, we would be scaling the whole error
		#That’s why you normalize by the sum, so the weights behave like probabilities.
		weights = weights / weights.sum()
		mse = np.average(se, weights=weights)
		return mse
	elif method == 'mae':
		ae = np.abs(np.log(y_obs[mask]) - np.log(y_exp[mask]))
		weights = distances ** -1.08
		weights = weights / weights.sum()
		mae = np.average(ae, weights=weights)
		return mae

### Apply factors
def apply_factors_to_records(matrix,resolution, chromosome_length,factor):
     
    '''
    Apply a scaling factor to all counts in a MatrixZoomData object and return adjusted records.
    Args: 
    matrix: a matrix in the form getMatrixZoomData from hic-straw
    resolution: resolution of the matrix
    chromosome_length: length of the chromosome
    factor: final factor to be applied to all counts
    Returns:
    adjusted_records: list of tuples (i,j,adjusted_count)
    '''
    # Get all records
    records = matrix.getRecords(0, chromosome_length, 0, chromosome_length)
    basis = infer_bins_unit(records, resolution, chromosome_length)
        
    adjusted_records = []
    for r in records:
        if basis == 'bp':
            format_bins = abs(r.binX - r.binY)
        elif basis == 'index':
            format_bins = abs(r.binX - r.binY) * resolution
        r.counts = r.counts * factor
        adjusted_records.append((r.binX, r.binY, r.counts))
        adjusted_records.append((r.binY, r.binX, r.counts))
    return adjusted_records, basis

def records_to_coomatrix(adj_records,basis,chromosome_length,resolution):
    '''
    	Arguments: 
     adj_records: list of tuples (i,j,adjusted_count)'
     chromosome_length: length of the chromosome
     resolution: resolution of the matrix
     Returns:
     sparse matrix in coo format
    '''
    num_bins = math.ceil(chromosome_length / resolution)
    rows, cols, data = [], [], []
    for i,j,val in adj_records:
        if basis == 'bp':
            i = i // resolution
            j = j // resolution
        elif basis == 'index':
            i = i
            j = j
        rows.append(i)
        cols.append(j)
        data.append(val)
    sparse_matrix = coo_matrix((data, (rows, cols)), shape=(num_bins, num_bins))
    return sparse_matrix
  
def coo2dense(coo_matrix):
    dense_matrix = coo_matrix.toarray()  
    return dense_matrix      

def get_corrected_matrix(matrix,resolution,chromosome_length,factor):
    '''
    Apply a scaling factor to all counts in a MatrixZoomData object and return the corrected dense matrix.
    Args:
    matrix: a matrix in the form getMatrixZoomData from hic-straw
    resolution: resolution of the matrix
    chromosome_length: length of the chromosome
    factor: final factor to be applied to all counts
    Returns:
    dense_matrix: corrected contact matrix as a 2D numpy array
    '''
    adjusted_records, basis = apply_factors_to_records(matrix,resolution, chromosome_length,factor)
    coo_matrix = records_to_coomatrix(adjusted_records,basis,chromosome_length,resolution)	
    dense_matrix = coo2dense(coo_matrix)
    dense_matrix = np.triu(dense_matrix,k=1)  # Keep only upper triangle
    return dense_matrix


def apply_factors_to_dense(matrix, factor,resolution ):
    """
	Apply factors to a dense Hi-C matrix.
	Args:
		matrix: 2D numpy array (triangular contact matrix)
		resolution: bin resolution in base pairs
    """
    new_matrix = matrix.copy()  # <- make a copy so original is untouched
    n_bins = matrix.shape[0]
	# It will always be index because we are passing a dense matrix and we are iterating over n_bins range
    for i in range(n_bins):
        for j in range(i+1, n_bins):   # skip diagonal, only upper triangle
            new_matrix[i, j] = matrix[i, j] * factor
    return np.triu(new_matrix)

def apply_list_factors_to_dense(matrix, list_factor,resolution ):
	"""
	Apply factors to a dense Hi-C matrix.
	Args:
		matrix: 2D numpy array (triangular contact matrix)
		resolution: bin resolution in base pairs
	"""
	new_matrix = matrix.copy()  # <- make a copy so original is untouched

	n_bins = matrix.shape[0]
	max_index = len(list_factor) - 1 
	# It will always be index because we are passing a dense matrix and we are iterating over n_bins range
	for i in range(n_bins):
		for j in range(i+1, n_bins):   # skip diagonal, only upper triangle
			index = np.abs(j - i) - 1
			if index > max_index:
				factor = list_factor[-1]
			else:
				factor = list_factor[index]
			new_matrix[i, j] = matrix[i, j] * factor

	return np.triu(new_matrix)

def apply_dict_factors_to_dense(matrix, dict_factor,resolution ):
	"""
	Apply factors to a dense Hi-C matrix.
	Args:
		matrix: 2D numpy array (triangular contact matrix)
		resolution: bin resolution in base pairs
	"""
	new_matrix = matrix.copy()  # <- make a copy so original is untouched

	n_bins = matrix.shape[0]
	max_key = max(dict_factor.keys())
	# It will always be index because we are passing a dense matrix and we are iterating over n_bins range
	for i in range(n_bins):
		for j in range(i+1, n_bins):   # skip diagonal, only upper triangle
			index = np.abs(j - i) * resolution
			if index > max_key or (index not in dict_factor):
				factor = dict_factor[max_key]
				
			else:
				factor = dict_factor[index]
	
			new_matrix[i, j] = matrix[i, j] * factor

	return np.triu(new_matrix)
### Plots

def decay_boxPlot(distances1, counts_1, m1, b1,
               distances2, counts_2, m2, b2,upper_limit=7e6):
    
	'''Plot distance decay boxplots for two datasets with fitted lines.
	Args:
 	distances1/distances2: list of distances used in fitting for dataset 1/2
  	counts_1/counts_2: dict with distances as keys and list of counts as values for dataset 1/2 (all distances)
  	m1/m2: slope of the fitted line in log-log space for dataset 1/2
  	b1/b2: intercept of the fitted line in log-log space for dataset 1/2
  	upper_limit: upper distance limit for fitting (default is 7e6)
	Returns:
		None (displays the plot)

	'''
	fit_y = linear_distance_decay(b1, m1, distances1) #Line
	fit_y2 = linear_distance_decay(b2, m2, distances2) 

	keys_1 = list(counts_1.keys())
 
	values_1 = list(counts_1.values())
	values_2 = [counts_2[k] for k in keys_1]
 
	width = 0.04 * np.array(keys_1)
	shift = width / 2
 
	fig = plt.figure(figsize=(10,5), dpi=240)

	#Yellow
	plt.boxplot(
		values_1,
		positions=np.array(keys_1) - shift,   
		widths=width,
		showfliers=False,
		patch_artist=True,
  		boxprops=dict(facecolor="#F2A3A3", color="#C94C4C", linewidth=0.5),
   		whiskerprops=dict(color="#C94C4C", linewidth=0.5),
    		capprops=dict(color="#C94C4C", linewidth=0.5),
    		medianprops=dict(color="#7A0A0A", linewidth=0.8)
		)
 
	plt.plot(distances1, fit_y, color='#492E00', label=f'Dataset A | α={-m1:.2f}', linewidth=0.5)
 
	#Teal
	plt.boxplot(
		values_2,
		positions=np.array(keys_1) + shift,   # shift right
		widths=width,
		showfliers=False,
		patch_artist=True,
		boxprops=dict(facecolor="#D6D8DB", color="#8A8D91", linewidth=0.5),
		whiskerprops=dict(color="#8A8D91", linewidth=0.5),
		capprops=dict(color="#8A8D91", linewidth=0.5),
		medianprops=dict(color="#4A4C4E", linewidth=0.8)
		)
		
	plt.plot(distances2, fit_y2, color='#536562', label=f'Dataset B | α={-m2:.2f}',linewidth=0.5)
	
 
	# Dinamically add y-axis vline
	# two mins and two maxs, because we are getting the min of each list but then we need the general min
	ymin = min(min(min(values_1)),min(min(values_2)))
	ymax = max(max(max(values_1)),max(max(values_2)))

 
	plt.vlines(upper_limit, ymin=ymin, ymax=ymax, 
            colors='grey', linestyles='dotted', linewidth=0.5)
 
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
 
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Genomic distance (bp)')
	plt.ylabel('Interaction Frequency')
	plt.legend().remove()
	plt.tight_layout()
	plt.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.3)
	plt.show()
	return fig

def plot_both_matrices(raw_matrix1,raw_matrix2):
	fig, axes = plt.subplots(1, 3, figsize=(12, 5),dpi=250)
	#Generate full matrix
	raw_matrix1 = raw_matrix1 + raw_matrix1.T
	raw_matrix2 = raw_matrix2 + raw_matrix2.T
	im1 = axes[0].imshow(raw_matrix1, cmap="OrRd", 
                    vmin=0, vmax=np.percentile(raw_matrix1, 97))
	axes[0].set_title("Sample A")
	plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)

	im2 = axes[1].imshow(raw_matrix2, cmap="OrRd", 
                    vmin=0, vmax=np.percentile(raw_matrix2, 97))
 
	axes[1].set_title("Sample B")
	plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)

	diff = np.array(raw_matrix1) - np.array(raw_matrix2)
	vmax = np.percentile(diff, 99)  

	im3 = axes[2].imshow(diff, cmap="binary", vmin=0, vmax=vmax)
	axes[2].set_title("Absolute difference (Sample A - Sample B)")
	plt.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04)

	plt.tight_layout()
	plt.show()

def plot_individual_matrix(matrix):
	plt.figure(figsize=(3,3), dpi=250)
	plt.imshow(matrix, cmap='RdYlBu_r',vmin=0,vmax=np.percentile(matrix, 99))
	plt.colorbar(fraction=0.046, pad=0.04)
	plt.show()


### Visualization 
def matrix2df(matrix):
	'''Args:
	matrix: 2D numpy array resultant from getRecordsAsMatrix 
	Returns:
	DataFrame with three columns: bin1, bin2, count'''
	matrix = np.triu(matrix,k=1)  # Keep only upper triangle
	df = pd.DataFrame([(i, j, matrix[i, j]) for i in range(matrix.shape[0]) for j in range(i+1,matrix.shape[1])])
	return df

def matrix2longdf(matrix,resolution=None):
	'''Args:
	matrix: 2D numpy array resultant from getRecordsAsMatrix 
	Returns:
	DataFrame with five columns: start1, end1, start2, end2, count [0,1,2,3,4]'''
	if resolution is None:
		resolution = 1

	matrix = np.triu(matrix,k=1)
	df = pd.DataFrame([(i * resolution, (i + 1) * resolution, j * resolution, (j + 1) * resolution, matrix[i, j]) for i in range(matrix.shape[0]) for j in range(i+1,matrix.shape[0]) ])
	return df

def merge_matrices_df_log2fc(df1,df2,resolution):
	'''Merge two contact count DataFrames and calculate distance and log2 fold-change.
		Args:
	 		df1: DataFrame with columns ['0', '1', '2'] for dataset 1
			df2: DataFrame with columns ['0', '1', '2'] for dataset 2
   			resolution: resolution in base pairs	
	Returns:
		both_df: DataFrame with columns ['bin1', 'bin2', 'count_1', 'count_2', 'distance', 'log2_fc']
	'''
	both_df = pd.merge(df1, df2, on=[0, 1], suffixes=('_1', '_2'))
	both_df.columns = ['bin1', 'bin2', 'count_1', 'count_2']
 
	both_df['distance'] = resolution * (np.abs(both_df['bin1'] - both_df['bin2']))
	both_df['log2_fc'] = np.log2((both_df['count_1'] + 0.01) / (both_df['count_2'] + 0.01))
	return both_df

def merge_matrices_df_difference(df1,df2,resolution):
	'''Merge two contact count DataFrames and calculate distance and log2 fold-change.
		Args:
	 		df1: DataFrame with columns ['0', '1', '2'] for dataset 1
			df2: DataFrame with columns ['0', '1', '2'] for dataset 2
   			resolution: resolution in base pairs	
	Returns:
		both_df: DataFrame with columns ['bin1', 'bin2', 'count_1', 'count_2', 'distance', 'log2_fc']
	'''
	both_df = pd.merge(df1, df2, on=[0, 1], suffixes=('_1', '_2'))
	both_df.columns = ['bin1', 'bin2', 'count_1', 'count_2']
 
	both_df['distance'] = resolution * (np.abs(both_df['bin1'] - both_df['bin2']))
	both_df['diff'] = (both_df['count_1'] + 1) - (both_df['count_2'] +1)
	return both_df

def define_plot(merged_df,log2_fc_cutoff,upper_limit):
	lg2fc_cutoff_pos = log2_fc_cutoff
	lg2fc_cutoff_neg = -log2_fc_cutoff
	both_df = merged_df.copy()
	'''Plot a MA plot comparing a processed merged DataFrame.
	'''
	fig = plt.figure(figsize=(8, 6), dpi=200)
	plt.scatter(both_df['distance'], both_df['log2_fc'], alpha=0.5, s=1, color='black')
	plt.hlines(0, xmin=0, xmax=both_df['distance'].max(), colors='red', linestyles='dashed')
	plt.hlines(lg2fc_cutoff_pos, xmin=0, xmax=both_df['distance'].max(), colors='gray', linestyles='dashed')
	plt.hlines(lg2fc_cutoff_neg, xmin=0, xmax=both_df['distance'].max(), colors='gray', linestyles='dashed')
	
	plt.vlines(upper_limit, ymin=both_df['log2_fc'].min(), ymax=both_df['log2_fc'].max(),
		  colors='red', linestyles='dotted', linewidth=0.5)
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
 
	plt.xlabel('Genomic Distance (bp)')
	plt.xscale('log') #the values stay the same but the axis is log scaled
	#The tick marks are ploted on an equal visual distance
	plt.ylabel('Log2 Fold Change (Sample A / Sample B)')
	return fig

def plot_MA(matrix1,matrix2,resolution,log2_fc_cutoff,upper_limit=7e6):
	'''Plot a MA plot comparing two contact matrices.
	Args:
	matrix1: 2D numpy array resultant from getRecordsAsMatrix for dataset 1
	matrix2: 2D numpy array resultant from getRecordsAsMatrix for dataset 2
	resolution: resolution in base pairs
	log2_fc_cutoff: log2 fold-change cutoff for highlighting significant changes'''
	df1 = matrix2df(matrix1)
	df2 = matrix2df(matrix2)
	merged_df = merge_matrices_df_log2fc(df1,df2,resolution)
	fig = define_plot(merged_df,log2_fc_cutoff,upper_limit)
	plt.tight_layout()
	return fig


### Statistical Significance
#This won't work becayse HiC maps arent independent, nearby bins share signal due domain continuity.
# If we only randomize each diagonal, we are removing this correlation and increase variance per bin, making the background uniform and not correlated with noise 



def create_background_matrix(matrix1,matrix2):
	'''Generate a background matrix by adding two contact matrices.'''
	background_matrix = (matrix1 + matrix2)
	return background_matrix


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
	"""downsampling method"""
	tag_mat= dense2tag(matrix) #tag_len is the number of non-zero elements in the matrix
	sample_idx = np.random.choice(len(tag_mat), int(n_reads)) #from tag len randomly select n_reads number of items
	sample_tag = tag_mat[sample_idx]# getting coordinates from the coo matrix of the selected items
	down_mat = tag2dense(sample_tag, matrix.shape[0])
	return down_mat

def get_reads_count(matrix):
	'''Get the total number of reads (upper triangle) in a contact matrix.'''
	total_reads = int(np.sum(np.triu(matrix, k=1)))
	return total_reads

def perform_permutation_test(matrix1, matrix2, resolution,method='log2fc', permutations_n=1000):
	df_obs1 = matrix2df(matrix1)
	df_obs2 = matrix2df(matrix2)
	if method == 'log2fc':
		merged_obs_df = merge_matrices_df_log2fc(df_obs1, df_obs2, resolution=resolution)
	elif method == 'difference':
		merged_obs_df = merge_matrices_df_difference(df_obs1, df_obs2, resolution=resolution)
	background_matrix = create_background_matrix(matrix1, matrix2)
	n_reads1 = get_reads_count(matrix1)
	n_reads2 = get_reads_count(matrix2)
	all_perms = merged_obs_df.copy()

	for p in range(permutations_n):
		print('Starting subsampling', p+1)
		down_mat1 = matrix_downsampling(background_matrix, n_reads1)
		down_mat2 = matrix_downsampling(background_matrix, n_reads2)

		df1 = matrix2df(down_mat1)
		df2 = matrix2df(down_mat2)
		if method == 'difference':
			permutation_dataframe = merge_matrices_df_difference(df1, df2, resolution=resolution)  
			permutation_dataframe = permutation_dataframe[['bin1', 'bin2', 'diff']].rename(
							columns={'diff': f'diff_{p+1}'}
			)	
		elif method == 'log2fc':
			permutation_dataframe = merge_matrices_df_log2fc(df1, df2, resolution=resolution)
			permutation_dataframe = permutation_dataframe[['bin1', 'bin2', 'log2_fc']].rename(
							columns={'log2_fc': f'log2_fc_{p+1}'}
			)
		all_perms = pd.merge(all_perms, permutation_dataframe, on=['bin1', 'bin2'], suffixes=('_obs', '_perm'))
	
	#Here is the filter before performing p-value calculations
	#all_perms = all_perms[all_perms['distance'] <= upper_limit].reset_index(drop=True)
	all_perms = all_perms[(all_perms['count_1'] > 0.0) & (all_perms['count_2'] > 0.0)].reset_index(drop=True)
	
	return all_perms

def get_p_values(merged_df,method='log2fc'):
	'''Calculate empirical p-values
	Args:
		merged_df: DataFrame with observed and permutation differences
	Returns:
		merged_df: DataFrame with an additional 'p_value' column
	'''
	if method == "difference":
		num_permutations = len([col for col in merged_df.columns if col.startswith('diff_')])
	elif method == "log2fc":
		num_permutations = len([col for col in merged_df.columns if col.startswith('log2_fc_')])
	p_values = []

	for _, row in merged_df.iterrows():
		if method == "log2fc":
			obs_diff = abs(row['log2_fc'])
			# count how many permuted absolute log2 fcs are >= observed absolute log2 fc
			count = sum(abs(row[f'log2_fc_{i+1}']) > obs_diff for i in range(num_permutations))
		elif method == "difference":
			obs_diff = abs(row['diff'])  # take absolute value of observed
			# count how many permuted absolute diffs are >= observed absolute diff
			count = sum(abs(row[f'diff_{i+1}']) > obs_diff for i in range(num_permutations))
			
		# empirical p-value 
		p_value = (count) / (num_permutations)
		p_values.append(p_value)

	merged_df['p_value'] = p_values
	return merged_df
 
def reformat_permutation_df(df, resolution, upper_limit=7e6, lower_limit=0):

    df = df.copy()  # FORCE real copy, avoids hidden duplication
    
    df.columns = ['bin1', 'bin2', 'count_1', 'count_2', 'difference', 'p_value']
    
    # keep upper triangle
    df = df.loc[df['bin2'] > df['bin1']].copy()
    
    # distance
    df['distance'] = resolution * (df['bin2'] - df['bin1']).abs()
    
    # distance filter
    df = df.loc[
        (df['distance'] > lower_limit)
    ].copy()
    
    # count filter
    df = df.loc[
        (df['count_1'] > 0.0) &
        (df['count_2'] > 0.0)
    ].copy()
    
    # log2 fold change
    df['log2_fc'] = np.log2((df['count_1'] + 0.01) /
                            (df['count_2'] + 0.01))
    
    return df[['bin1','bin2','count_1','count_2','log2_fc','distance','p_value']]

def adjust_p_values(merged_df):
	'''Apply Benjamini-Hochberg correction to p-values in the
	merged_df: DataFrame with a 'p_value' column
	Returns:
		merged_df: DataFrame with an additional 'p_value_adj' and FDR column
	'''
	#Filtering should be done before to avoid adjusting p-values for tests we are not interested in.
	#If we include them then we will be saying they are part on the 5% FDR when in reality we don't care about them.
	#	An adjusted p-value is the minimum FDR threshold at which that test would still be significant.
	pvals = merged_df['p_value'].values
	_, pvals_adj, _, FDR = multipletests(pvals, method='fdr_bh')
	merged_df['p_value_adj'] = pvals_adj
	merged_df['FDR'] = FDR
	return merged_df

def prior_probability(counts_1,counts_2,true_fc,tau2=1.0):
	theta_vals = np.linspace(true_fc.min(), true_fc.max(),len(counts_1))
	# Assume a normal prior with a mean of 0 and variance tau2
	prior = np.exp(-theta_vals**2 / (2 * tau2)) / np.sqrt(2 * np.pi * tau2)
	prior /= np.trapz(prior, theta_vals)  # Normalize the prior
	return prior, theta_vals

def likelihood(counts_1,counts_2,theta_vals):
	p_theta = 2**theta_vals / (1 + 2**theta_vals)
	logL = (counts_1 * np.log(p_theta) + counts_2 * np.log(1 - p_theta))
	logL -= logL.max(axis=1, keepdims=True)  # For numerical stability
	L = np.exp(logL)
	L /= np.trapz(L, theta_vals, axis=1)  # Normalize the likelihoods
	return L

def posterior_distribution(prior, L, theta_vals):
	posterior = L * prior[None, :]
	posterior /= np.trapz(posterior, theta_vals, axis=1)  # Normalize the posterior
	return posterior

def posterior_summary(posterior, theta_vals):
	post_mean = np.trapz(theta_vals * posterior, theta_vals, axis=1)
	i_map = np.argmax(posterior, axis=1)
	theta_map = theta_vals[i_map]
	return post_mean, theta_map

def bayesian_fc_estimation(permutations_df,tau2=1.0):
	"""
	Estimate the posterior mean and variance of the log2 fold change using a Bayesian approach.
	
	Parameters:
	permutations_df (pd.DataFrame): DataFrame containing observed log2 fold changes and their variances.
	tau2 (float): Prior variance for the log2 fold change.
	
	Returns:
	pd.DataFrame: Updated DataFrame with posterior mean and variance of log2 fold change.
	"""
	counts_1 = permutations_df['count_1'].to_numpy()[:,None]
	counts_2 = permutations_df['count_2'].to_numpy()[:,None]
	true_log2fc = permutations_df['log2_fc'].to_numpy()[:,None]
	prior, theta_vals = prior_probability(counts_1, counts_2, true_log2fc, tau2)
	L = likelihood(counts_1, counts_2, theta_vals)
	posterior = posterior_distribution(prior, L, theta_vals)
	post_mean, theta_map = posterior_summary(posterior, theta_vals)

	final_df = permutations_df.copy()
	final_df['post_mean_log2fc'] = post_mean
	final_df['post_map_log2fc'] = theta_map
	final_df = final_df[['bin1','bin2','count_1','count_2','log2_fc','post_mean_log2fc','post_map_log2fc','distance','p_value','p_value_adj']]


	return final_df

 
def significant_interactions(permutations_df,pval,lg2_fc_cutoff=0,column='p_value_adj'):
	'''Filter significant interactions based on a p-value threshold.
	Args:
		permutations_df
		pval: significance threshold for calling a difference significant
		column: user can choose between p_value or p_value_adj, default is 'p_value_adj'
	Returns:
		predicted_positive: DataFrame with significant interactions
	'''
	# Total number of tests
	permutations_df['tests'] = len(permutations_df)
	predicted_positive = permutations_df[(permutations_df[column]< pval) & (permutations_df['distance'] > 0)]
	predicted_positive = predicted_positive[(predicted_positive['post_map_log2fc'] >= lg2_fc_cutoff) | (predicted_positive['post_map_log2fc'] <= -lg2_fc_cutoff)]
	predicted_positive = predicted_positive[['bin1','bin2','count_1','count_2','log2_fc','post_mean_log2fc','post_map_log2fc','distance','p_value','p_value_adj']]

	return predicted_positive

def significant_interactions_no_bayesian(permutations_df,pval,log2_fc_cutoff=0,column='p_value_adj'):
	'''Filter significant interactions based on a p-value threshold.
	Args:
		permutations_df
		pval: significance threshold for calling a difference significant
		column: user can choose between p_value or p_value_adj, default is 'p_value_adj'
	Returns:
		predicted_positive: DataFrame with significant interactions
	'''
	# Total number of tests
	permutations_df['tests'] = len(permutations_df)
	predicted_positive = permutations_df[(permutations_df[column]< pval) & (permutations_df['distance'] > 0)]
	predicted_positive = predicted_positive[(predicted_positive['log2_fc'] >= log2_fc_cutoff) | (predicted_positive['log2_fc'] <= -log2_fc_cutoff)]
	predicted_positive = predicted_positive[['bin1','bin2','count_1','count_2','log2_fc','distance','p_value','p_value_adj']]

	return predicted_positive


# Neigborhood functions

def get_neighbors_coordinates(center, K=5, verbose=False):
    """
    General K-degree neighbor expansion.
    Returns:
        A list of lists:
            neighbors_per_degree[d] = list of (i, j) for degree d+1
    """
    seen = set()
    neighbors_per_degree = [[] for _ in range(K)]

    # --- First degree neighbors ---
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

    # --- Higher degree neighbors ---
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

def get_neighbors_M(M,i, j):
    
    """Return direct (4-connected) neighbors of a matrix element.
	d is diagonal of the spike in
	row, column,value"""
    neighbors = []
    n_rows, n_cols = M.shape

    if i + 1 < n_rows: #ensures we don't go out of bounds on the rows (lower bound)
        neighbors.append((i + 1, j))
    if i - 1 >= 0: #ensures we don't go out of bounds on the rows (upper bound)
        neighbors.append((i - 1, j))
    if j + 1 < n_cols: #ensures we don't go out of bounds on the columns (right bound)
        neighbors.append((i, j + 1))
    if j - 1 >= 0: #ensures we don't go out of bounds on the columns (left bound)
        neighbors.append((i, j - 1))

    return neighbors

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

def neighbor_support(final_results):
	sig_interactions = set(
		zip(final_results['bin1'], final_results['bin2'])
	)

	k = 1

	neighbor_support = []

	if k > 0:
		for a, b in zip(final_results['bin1'], final_results['bin2']):
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

### Synthetic datasets 

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
		return down_mat



#Others

def get_chromosome_length(hic_path,chromosome):
	hic = hicstraw.HiCFile(hic_path)
	chromosomes = hic.getChromosomes()
	for chrom in chromosomes:
		if chrom.name == chromosome:
			chromosome_length = chrom.length
			return int(chromosome_length)
		else:
			continue
	raise ValueError(f"Chromosome {chromosome} not found in the Hi-C file.")

def pre_counts(synthetic_counts1,synthetic_counts2,d_1,d_2,y_1,y_2,desired_alpha,resolution):
    
	counts_1_v2 = synthetic_counts1.copy()
	counts_2_v2 = synthetic_counts2.copy()

	b_refit_1 = refit_intercept(y_1, d_1, desired_alpha)
	b_refit_2 = refit_intercept(y_2, d_2, desired_alpha)

	b_mean = (b_refit_1 + b_refit_2) / 2

	desired_fit_1 = linear_distance_decay(b_mean,desired_alpha,d_1)
	desired_fit_2 = linear_distance_decay(b_mean,desired_alpha,d_2)

	factors_r = get_individual_factors(y_1,desired_fit_1)
	factors_s = get_individual_factors(y_2,desired_fit_2)

	#Apply factors to the dictionary counts
	scaled_counts_1,full_factors_r=apply_factors_to_dict_counts(counts_1_v2,factors_r,resolution)
	scaled_counts_2,full_factors_s=apply_factors_to_dict_counts(counts_2_v2,factors_s,resolution)
	
	return scaled_counts_1,scaled_counts_2,full_factors_r,full_factors_s

## Higher resolutions

def generate_cool_matrix(matrix,chromosome,resolution,output_path):
	n = matrix.shape[0]

	#bin table is always sqared
	bins = pd.DataFrame({
		'chrom': [chromosome]*n,
		'start': np.arange(0, n*resolution, resolution),
		'end': np.arange(resolution, (n+1)*resolution, resolution),
		"KR": np.ones(n, dtype=np.float32) # Dummy normalization vector

	})

	i,j = np.triu_indices(n)

	# Pixels table is a 2-D array of every pixel in the upper triangle of the matrix
	pixels = pd.DataFrame({
		'bin1_id': i,
		'bin2_id': j,
		'count': matrix[i,j].astype(int)
	})
	pixels = pixels[pixels["count"] > 0]

	cooler.create_cooler(
		output_path,
		bins,
		pixels,
		dtypes={"count": "int32"},
		ordered=True, symmetric_upper=True
	)

	#mcool 

	mcool_name = output_path.replace('.cool', '.mcool')
	
	if os.path.exists(mcool_name):
		os.remove(mcool_name)
	
	cooler.zoomify_cooler(output_path,chunksize=256,
		outfile=mcool_name,
		resolutions=[resolution] 
	)

	return mcool_name
