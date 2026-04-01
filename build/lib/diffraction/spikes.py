import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from . import utils as DifFracTion_utils

def plot_spike_boxplots(data, positions, labels):
    plt.figure(figsize=(5, 3))
    plt.boxplot(data, 
                positions=positions, 
                widths=0.22, patch_artist=True, showfliers=False,
                labels=labels,
                boxprops=dict(facecolor="#565FAF", color="#4C71BF"), 
                medianprops=dict(color="#D3DDEE"))

    m, b,_ = fit_line(data, positions)
    fit_y = m * positions + b
    plt.plot(positions, fit_y, color="#000000", linestyle='-', linewidth=1, label='Fitted line')
    
    plt.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.6)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title(f"y = {m:.6f}x + {b:.6f}")
    plt.tight_layout()
    plt.ylim(0, None)
    plt.show()     
    plt.close()   
    
def fit_line(data,positions):
	# median or mean?
	medians = np.array([np.median(d) for d in data])
	slope, intercept = np.polyfit(positions, medians, 1) #1 is linear
	return slope, intercept,medians

def highlight_spikes(matrix, spikes):
	fig, ax = plt.subplots(figsize=(4,4), dpi=200)
	percentile_99 = np.percentile(matrix, 99)
	ax.imshow(matrix, cmap="RdYlBu_r", vmin=0, vmax=percentile_99)

	for (i, j) in spikes:
		circle = Circle((j, i), radius=1, edgecolor='white', facecolor='none', lw=0.2)
		ax.add_patch(circle)
	plt.show()
	return fig
	
def get_spikes_signal(M, spike,verbose=False):
	'''Return a list of lists containing the spike values from the matrix M at the given spike coordinates.'''
	
	
	individual_spike_values = [] #Add a list for each spike in so we can keep track of the respective values and neighbors
	individual_spike_values.append(M[spike[0], spike[1]])
	
	if verbose:
		print("Spike values:", individual_spike_values)
	return individual_spike_values

def get_spikes_coordinates(M, spike,verbose=False):
	'''Return a list of lists containing the spike values from the matrix M at the given spike coordinates.'''
	
	
	individual_spike_values = [] #Add a list for each spike in so we can keep track of the respective values and neighbors
	individual_spike_values.append((spike[0], spike[1], M[spike[0], spike[1]]))
	
	if verbose:
		print("Spike values:", individual_spike_values)
	return individual_spike_values

def get_neighbors_coordinates(M, spikes, K=5, verbose=False):
    """
    General K-degree neighbor expansion.
    Returns:
        A list of lists:
            neighbors_per_degree[d] = list of (i, j, val) for degree d+1
    """
    seen = set()
    neighbors_per_degree = [[] for _ in range(K)]

    # --- First degree neighbors ---
    frontier = []  # list of (i, j, val)

    for (i0, j0) in spikes:
        seen.add((i0, j0))
        first_neighbors = get_neighbors(M, i0, j0)

        for i, j, val in first_neighbors:
            if (i, j) not in seen:
                frontier.append((i, j, val))
                neighbors_per_degree[0].append((i, j, val))
                seen.add((i, j))

    if verbose:
        print(f"Degree 1 neighbors: {len(neighbors_per_degree[0])}")

    # --- Higher degree neighbors ---
    for degree in range(1, K):
        next_frontier = []
        for (i, j, val) in frontier:
            new_neighbors = get_neighbors(M, i, j)
            for ni, nj, nval in new_neighbors:
                if (ni, nj) not in seen:
                    neighbors_per_degree[degree].append((ni, nj, nval))
                    next_frontier.append((ni, nj, nval))
                    seen.add((ni, nj))

        if verbose:
            print(f"Degree {degree+1} neighbors: {len(neighbors_per_degree[degree])}")

        frontier = next_frontier

        # No more expansion possible
        if len(frontier) == 0:
            break

    return (*neighbors_per_degree,)


def analyze_neighbors(M, spikes, K=5, verbose=False):
    """
    General k-degree neighbor extraction.
    Returns: list of lists (values_per_degree)
        values_per_degree[d] = list of values for degree d+1
    """
    seen = set()
    values_per_degree = [[] for _ in range(K)]

    # Initialize frontier with all first-degree neighbors for each spike
    frontier = []  # list of (i,j,val)

    # Collect all first-degree neighbors from all spikes
    for (i0, j0) in spikes:
        seen.add((i0, j0))
        first_neighbors = get_neighbors(M, i0, j0)
        for i, j, val in first_neighbors:
            if (i, j) not in seen:
                frontier.append((i, j, val))
                seen.add((i, j))
        # Save first-degree values
        values_per_degree[0].extend([val for _, _, val in first_neighbors])

    if verbose:
        print(f"Degree 1 neighbors: {len(values_per_degree[0])}")

    # Iterate through higher degrees
    for degree in range(1, K):
        next_frontier = []
        for (i, j, val) in frontier:
            new_neighbors = get_neighbors(M, i, j)
            for ni, nj, nval in new_neighbors:
                if (ni, nj) not in seen:
                    next_frontier.append((ni, nj, nval))
                    seen.add((ni, nj))
                    values_per_degree[degree].append(nval)

        if verbose:
            print(f"Degree {degree+1} neighbors: {len(values_per_degree[degree])}")

        frontier = next_frontier  # Move frontier outward

        # If no new neighbors at this degree, stop early
        if len(frontier) == 0:
            break

    return (*values_per_degree,) # unpack the list of lists into separate lists

def analyze_neighbors_distance(M, spikes, dc,K=5, verbose=False):
    """
    General k-degree neighbor extraction.
    Returns: list of lists (values_per_degree)
        values_per_degree[d] = list of values for degree d+1
    """
    seen = set()
    values_per_degree = [[] for _ in range(K)]

    # Initialize frontier with all first-degree neighbors for each spike
    frontier = []  # list of (i,j,val, diag)

    # Collect all first-degree neighbors from all spikes
    for (i0, j0) in spikes:
        seen.add((i0, j0))
        first_neighbors = get_neighbors_with_diagonal(M, i0, j0)
        for i, j, val, diag in first_neighbors:
            if (i, j) not in seen:
                frontier.append((i, j, val, diag)) # First neighbors with diagonal
                seen.add((i, j))
        # Save first-degree values
        values_per_degree[0].extend([(i,j,val,dc,diag) for i, j, val, diag in first_neighbors])

    if verbose:
        print(f"Degree 1 neighbors: {len(values_per_degree[0])}")

    # Iterate through higher degrees
    for degree in range(1, K):
        next_frontier = []
        for (i, j, val, diag) in frontier:
            new_neighbors = get_neighbors_with_diagonal(M, i, j)
            for ni, nj, nval, ndiag in new_neighbors:
                if (ni, nj) not in seen:
                    next_frontier.append((ni, nj, nval, ndiag)) # Keep diagonal info
                    seen.add((ni, nj))
                    values_per_degree[degree].append((ni,nj,nval,dc, ndiag))

        if verbose:
            print(f"Degree {degree+1} neighbors: {len(values_per_degree[degree])}")

        frontier = next_frontier  # Move frontier outward

        # If no new neighbors at this degree, stop early
        if len(frontier) == 0:
            break

    return (*values_per_degree,) # unpack the list of lists into separate lists


def get_neighbors(M, i, j):
    
    """Return direct (4-connected) neighbors of a matrix element.
	d is diagonal of the spike in
	row, column,value"""
    neighbors = []
    n_rows, n_cols = M.shape

    if i + 1 < n_rows:
        neighbors.append((i + 1, j, M[i + 1][j]))
    if i - 1 >= 0:
        neighbors.append((i - 1, j, M[i - 1][j]))
    if j + 1 < n_cols:
        neighbors.append((i, j + 1, M[i][j + 1]))
    if j - 1 >= 0:
        neighbors.append((i, j - 1, M[i][j - 1]))

    return neighbors

def get_neighbors_with_diagonal(M, i, j):
    
    """Return direct (4-connected) neighbors of a matrix element.
	d is diagonal of the spike in
	row, column,value, diagonal of spikein """
    neighbors = []
    n_rows, n_cols = M.shape
    def diag(r, c):
        return np.abs(r - c)
    
    if i + 1 < n_rows:
        neighbors.append((i + 1, j, M[i + 1][j],diag(i + 1,j)))
    if i - 1 >= 0:
        neighbors.append((i - 1, j, M[i - 1][j],diag(i - 1,j)))
    #if j + 1 < n_cols:
    #    neighbors.append((i, j + 1, M[i][j + 1],diag(i,j + 1)))
    if j - 1 >= 0:
        neighbors.append((i, j - 1, M[i][j - 1],diag(i,j - 1)))

    return neighbors






# Functions to model spike in signal


def spike_alteration(k,original_spike_signal,slope):
	''' Modify spike signal by a factor k and adjust intercept to maintain slope
	Args:
		k: multiplication factor for the spike signal
		original_spike_signal: list of original spike signal values
		slope: slope of the original line'''
	y_1_modified = [f * k for f in original_spike_signal] #for each spike value f, multiply by k
	y_1_median_modified = np.median(y_1_modified) #get median of modified spike values
	# Adjust intercept to maintain slope
	b_prime = y_1_median_modified - slope * 1 #assuming x=1 for spike position
	return slope, b_prime

def get_factors_from_medians(medians,m_prime,b_prime,positions):
	target_y = m_prime * positions + b_prime

	return target_y / medians


def rescale_line(slope, b_prime, positions):
	''' Rescale line using new slope and intercept
	Args:
		slope: slope of the line which is the same as original
		b_prime: new intercept of the line
		positions: list of x positions [1,2,3...] to calculate y values for'''
	return [slope * x + b_prime for x in positions]

def get_factors(slope,b,b_prime,positions):
	''' Get rescaling factors for each position based on original and modified line'''
	return [(slope * x + b_prime)/(slope * x + b) for x in positions]


def apply_perturbations(data, k_factors):
	'''Apply perturbations to data based on k_factors
	Args:
		data: list of lists, each sublist contains signal values (spike, neighbors1, neighbors2, etc)
		k_factors: list of multiplication factors corresponding to each sublist in data'''
	new_data = []
	for i in range(len(data)):
		new_data.append([val * k_factors[i] for val in data[i]])
	return new_data

def apply_perturbations_to_matrix(M, spikes, k_factors, K=None):
    """
    Apply multiplicative perturbations to Hi-C matrix M based on k_factors.

    Args:
        M : 2D numpy array (Hi-C matrix)
        spikes : list of (i, j) coordinates
        k_factors : list of multiplicative factors
                        k_factors[0] → spike
                        k_factors[1] → degree-1 neighbors
                        k_factors[2] → degree-2 neighbors
                        ...
        K : optional number of degrees (default = len(k_factors) - 1)

    Returns:
        Modified copy of M with perturbations applied
    """
    M_perturbed = M.copy()

    if K is None:
        K = len(k_factors) - 1  # spikes = degree 0

    # Set of visited coordinates to avoid double perturbation
    seen = set()

    for (i0, j0) in spikes:
        # Apply perturbation to spike itself
        M_perturbed[i0, j0] *= k_factors[0]
        seen.add((i0, j0))

        # ----- FIRST DEGREE NEIGHBORS -----
        frontier = []
        first_neighbors = get_neighbors(M, i0, j0)
        for i, j, val in first_neighbors:
            if (i, j) not in seen:
                M_perturbed[i, j] *= k_factors[1]
                frontier.append((i, j))
                seen.add((i, j))

        # ----- HIGHER DEGREES -----
        for degree in range(2, K+1):
            next_frontier = []

            for (i, j) in frontier:
                neighbors = get_neighbors(M, i, j)

                for ni, nj, val in neighbors:
                    if (ni, nj) not in seen:
                        M_perturbed[ni, nj] *= k_factors[degree]
                        next_frontier.append((ni, nj))
                        seen.add((ni, nj))

            frontier = next_frontier
            if len(frontier) == 0:
                break  # no more neighbors to expand into

    return M_perturbed

def apply_perturbations_to_matrix_distance_dependent(M, spikes, k_f,k_c, K=None):
    """
    Apply multiplicative perturbations to Hi-C matrix M based on k_factors.

    Args:
        M : 2D numpy array (Hi-C matrix)
        spikes : list of (i, j) coordinates
        k_f : list of multiplicative factors for farther distances
                        k_f[0] → spike
                        k_f[1] → degree-1 neighbors
                        k_f[2] → degree-2 neighbors
                        ...
        k_c : list of multiplicative factors for closer distances

        K : optional number of degrees (default = len(k_f) - 1)

    Returns:
        Modified copy of M with perturbations applied
    """
    M_perturbed = M.copy()

    if K is None:
        K = len(k_f) - 1  # spikes = degree 0

    # Set of visited coordinates to avoid double perturbation
    seen = set()

    for (i0, j0) in spikes:
        # Apply perturbation to spike itself
        M_perturbed[i0, j0] *= k_f[0]
        seen.add((i0, j0))

        dc=abs(i0 - j0)  # diagonal of the spikein

        # ----- FIRST DEGREE NEIGHBORS -----
        frontier = []
        first_neighbors = get_neighbors_with_diagonal(M, i0, j0)
        for i, j, val, diag in first_neighbors:
            if (i, j) not in seen:
                if diag > dc:
                    M_perturbed[i, j] *= k_f[1]  # farther
                elif diag < dc:
                    M_perturbed[i, j] *= k_c[1]  # closer
                else:
                    pass
                frontier.append((i, j))
                seen.add((i, j))

        # ----- HIGHER DEGREES -----
        for degree in range(2, K+1):
            next_frontier = []

            for (i, j) in frontier:
                neighbors = get_neighbors_with_diagonal(M, i, j)

                for ni, nj, val, diag in neighbors:
                    if (ni, nj) not in seen:
                        if diag > dc:
                            M_perturbed[ni, nj] *= k_f[degree]
                        elif diag < dc:
                            M_perturbed[ni, nj] *= k_c[degree]
                        else:
                            pass
                        next_frontier.append((ni, nj))
                        seen.add((ni, nj))

            frontier = next_frontier
            if len(frontier) == 0:
                break  # no more neighbors to expand into

    return M_perturbed

# Other utility functions 

def save_spike_and_neighbors_coordinates(spikes, all_neighbors, output_file):
	'''Save spike coordinates to a text file
	Args:
		spikes: list of tuples containing spike coordinates (i,j)
		all_neighbors: dictionary of neighbor lists per degree, with degree as key and list of (i,j,value) tuples as values
		output_file: path to the output text file'''
	
	with open(output_file, 'w') as f:
		for (i, j,v) in spikes:
			f.write(f"{i}\t{j}\t{v}\tspike\n")
		for deg, neighbors in all_neighbors.items():
			for (i, j,v) in neighbors:
				f.write(f"{i}\t{j}\t{v}\tneighbor{deg+1}\n")
                        

def get_diagonal_bin_signal(M, diagonal_position):
    '''Return a list of values from the matrix M at the given diagonal position.'''
    
    diagonal_values = []
    n_rows, n_cols = M.shape

    for i in range(n_rows):
        j = i + diagonal_position
        if 0 <= j < n_cols:
            if M[i, j] != 0:
                diagonal_values.append(M[i, j])
    
    return diagonal_values

def get_diagonal_position(coordinates):
	row, col = coordinates
	return abs(row - col)

def get_genomic_distance(diagonal_position, resolution):
	return diagonal_position * resolution

def run_spikes(spikein_coordinates, k, raw_matrix_sub,resolution,neighbor_degrees):
	non_modified_raw_matrix=raw_matrix_sub.copy()
	counts_dict= DifFracTion_utils.get_counts_by_distance_from_dense(non_modified_raw_matrix,resolution)
	_,_,m_fit,b_fit,a_fit= DifFracTion_utils.calculate_distance_decay(counts_dict,'median')

	all_neighbors = {deg: [] for deg in range(neighbor_degrees)}   # K = 8
	spikein_all_lists = []


	for spikein in spikein_coordinates:
		print(f'Analyzing spike-in at position: {spikein}')
		spike_diagonal_position = get_diagonal_position(spikein)
		spike_genomic_distance = get_genomic_distance(spike_diagonal_position, resolution)
		#m = -alpha
		expected_signal = np.exp(b_fit) * (spike_genomic_distance ** (-1.08))
		
		spike_signal=get_spikes_signal(non_modified_raw_matrix, spikein, verbose=False)
		
		n_f=spike_signal[0]/expected_signal
		
		print(f'Expected signal: {expected_signal}, Observed signal: {spike_signal[0]}, Normalization factor: {n_f}')
		
		k_prime = k/n_f

		print(k, k_prime)

		spikes_list=get_spikes_coordinates(non_modified_raw_matrix,spikein, verbose=False)
		spike_dup = []

		# It takes the current bin signal
		for i in spike_signal:
			spike_dup.extend([i, i]) # Duplicating spike signal for boxplot visualization
		
		#Coords are (row,col,value)
		signal_neigh=analyze_neighbors(raw_matrix_sub,spikes=[spikein],K=neighbor_degrees)
		coords_neigh=get_neighbors_coordinates(raw_matrix_sub,spikes=[spikein],K=neighbor_degrees)
		
        #access list by the index of the degree of neighborhood

		for deg in range(neighbor_degrees):
			all_neighbors[deg].extend(coords_neigh[deg])
		spikein_all_lists.extend(spikes_list)

		data = [spike_dup] + [signal_neigh[d] for d in range(neighbor_degrees)]
		labels = [f'Spike-in'] + [f'{i+1}' for i in range(neighbor_degrees)]
		positions = np.arange(1,len(data)+1)
		
		# Fit line to medians
		m, b , medians = fit_line(data, positions)

		if m > 0:
			print(f'Warning: Positive slope detected (m={m}). Check spike-in definition and neighborhood extraction.')
			for _ in spikes_list:
				spikein_all_lists.pop()
			for deg in range(neighbor_degrees):
				all_neighbors[deg].pop()
			continue

		plot_spike_boxplots(data,positions,labels)

		#Rescaling by k, slope should be the same, but intercept changes
		#### We have to modify this to take into consideration the diagonal position of the spike
		m_prime, b_prime = spike_alteration(k_prime, spike_signal, m)
		print(f'Original line: y = {m}x + {b}')

		# Get factors from actual medians
		k_factors = get_factors_from_medians(medians,m_prime,b_prime,positions)

		new_data = apply_perturbations(data, k_factors)
		
		m_check, b_check, _ = fit_line(new_data, positions)

		raw_matrix_sub = apply_perturbations_to_matrix(raw_matrix_sub, spikes=[spikein], k_factors=k_factors)

		plot_spike_boxplots(new_data,positions,labels)
		print(f'k_factors per position: {k_factors}')
	#raw_matrix_sub now contains all spike ins applied
	return raw_matrix_sub, spikein_all_lists, all_neighbors
