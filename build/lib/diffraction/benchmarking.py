from . import utils as DifFracTion_utils
from . import spikein as DifFracTion_spikes
import numpy as np
import pandas as pd
from functools import reduce
import os


def generate_and_save_spike_ins(
    hic_path: str,
    chromosome: str,
    resolution: int,
    spikes_k: int,
    neighbor_degrees: list,
    coordinate_file: str
):
    chromosome_length = DifFracTion_utils.get_chromosome_length(hic_path, chromosome)
    matrixZoom = DifFracTion_utils.prepare_MatrixZoomData(hic_path, chromosome, resolution)
    raw_matrix = DifFracTion_utils.rawMatrix(matrixZoom, chromosome_length)

    _, spikein_coordinates = DifFracTion_spikes.select_spikeins(
        raw_matrix,
        resolution,
        n_spikes=50,
        fold_change=spikes_k
    )

    altered_matrix, spikein_all_lists, all_neighbors = DifFracTion_spikes.apply_spikeins(
        raw_matrix,
        spikein_coordinates,
        neighbor_degrees,
        spikes_k
    )

    DifFracTion_spikes.save_spike_and_neighbors_coordinates(
        spikein_all_lists,
        all_neighbors,
        coordinate_file
    )

    return altered_matrix, raw_matrix, spikein_all_lists

def generate_HiCDCPlus_input(altered_matrix,
    raw_matrix,
    chromosome: str,
    resolution: int,
    spikes_k: float,
    output_path: str
):  
    #Made based on requirements from HiCDCPlus
    def generate_individual_frame(matrix, resolution):
        m = DifFracTion_utils.matrix2df(matrix)  # upper triangle already
        m.columns = ['startI','startJ','counts']
        m['startI'] = m['startI'] * resolution
        m['startJ'] = m['startJ'] * resolution
        m['chr'] = f'chr{chromosome}'
        m = m[['chr','startI','startJ','counts']]
        return m
    
    def save_individual_frame(m, resolution, spikes_k, condition):
        save_path = (
            f"{output_path}/HiCDCPlus_input_chr{chromosome}_res{resolution}_k{spikes_k}_{condition}.table"
        )
        m.to_csv(save_path, sep='\t', index=False, header=True)
        return save_path
    
    # convert to long format
    M1_1 = generate_individual_frame(altered_matrix, resolution)
    # A condition
    save_path_A = save_individual_frame(M1_1, resolution, spikes_k, "A1")
    
    altered_matrix_rep = DifFracTion_utils.synthetic_datasets(altered_matrix, 0.9999999) # this is because we already have spike-ins
    M1_2 = generate_individual_frame(altered_matrix_rep, resolution)
    save_path_A2 = save_individual_frame(M1_2, resolution, spikes_k, "A2")

     # B condition
    synthetic = DifFracTion_utils.synthetic_datasets(raw_matrix, 1)
    M2_1 = generate_individual_frame(synthetic, resolution)
    save_path_B = save_individual_frame(M2_1, resolution, spikes_k, "B1")
    
    synthetic_rep = DifFracTion_utils.synthetic_datasets(raw_matrix, 0.9999999)
    M2_2 = generate_individual_frame(synthetic_rep, resolution)
    save_path_B2 = save_individual_frame(M2_2, resolution, spikes_k, "B2")
    return save_path_A, save_path_A2, save_path_B, save_path_B2
    

def generate_HiCDCPlus_input_ds(
    raw_matrix,
    chromosome: str,
    resolution: int,
    downsample_factor: float,
    output_path: str
):  
    #Made based on requirements from HiCDCPlus
    def generate_individual_frame(matrix, resolution):
        m = DifFracTion_utils.matrix2df(matrix)  # upper triangle already
        m.columns = ['startI','startJ','counts']
        m['startI'] = m['startI'] * resolution
        m['startJ'] = m['startJ'] * resolution
        m['chr'] = f'chr{chromosome}'
        m = m[['chr','startI','startJ','counts']]
        return m
    
    def save_individual_frame(m, resolution, downsample_factor, condition):
        save_path = (
            f"{output_path}/HiCDCPlus_input_chr{chromosome}_res{resolution}_ds{downsample_factor}_{condition}.table"
        )
        m.to_csv(save_path, sep='\t', index=False, header=True)
        return save_path
    
    # convert to long format
    M1 = DifFracTion_utils.synthetic_datasets(raw_matrix,0.9999999)
    M1_1 = generate_individual_frame(M1, resolution)
    # A condition
    save_path_A = save_individual_frame(M1_1, resolution, downsample_factor, "A1")
    
    M1_2 = DifFracTion_utils.synthetic_datasets(raw_matrix, 0.9999999) # this is because we already have spike-ins
    M1_2 = generate_individual_frame(M1_2, resolution)
    save_path_A2 = save_individual_frame(M1_2, resolution, downsample_factor, "A2")


     # B condition
    synthetic = DifFracTion_utils.synthetic_datasets(raw_matrix, downsample_factor)
    M2_1 = generate_individual_frame(synthetic, resolution)
    save_path_B = save_individual_frame(M2_1, resolution, downsample_factor, "B1")
    
    synthetic_rep = DifFracTion_utils.synthetic_datasets(raw_matrix,downsample_factor)
    M2_2 = generate_individual_frame(synthetic_rep, resolution)
    save_path_B2 = save_individual_frame(M2_2, resolution, downsample_factor, "B2")
    return save_path_A, save_path_A2, save_path_B, save_path_B2


def generate_diffHic_input(altered_matrix,
    raw_matrix,
    chromosome: str,
    resolution: int,
    spikes_k: float,
    output_path: str):
    #A condition
    # Rep1
    m1 = DifFracTion_utils.matrix2longdf(altered_matrix, resolution)  # upper triangle already
    # Rep2
    altered_matrix_rep = DifFracTion_utils.synthetic_datasets(altered_matrix, 0.999999999999) # this is because we already have spike-ins
    m1_2 = DifFracTion_utils.matrix2longdf(altered_matrix_rep, resolution)

    #B condition
    #Rep1 
    synthetic = DifFracTion_utils.synthetic_datasets(raw_matrix, 1)
    m2 = DifFracTion_utils.matrix2longdf(synthetic, resolution)
    # Rep2
    synthetic_rep = DifFracTion_utils.synthetic_datasets(raw_matrix, 0.999999999999)
    m2_2 = DifFracTion_utils.matrix2longdf(synthetic_rep, resolution)


    # Merge 
    m1    = m1.rename(columns={4: "sampleA_rep1"})
    m1_2  = m1_2.rename(columns={4: "sampleA_rep2"})
    m2    = m2.rename(columns={4: "sampleB_rep1"})
    m2_2  = m2_2.rename(columns={4: "sampleB_rep2"})

    matrices_dfs = [m1, m1_2, m2, m2_2]
    matrices_merged = reduce( lambda left, right:
                             pd.merge(left,right,
                                      on=[0,1,2,3],
                                      how='inner').fillna(0),
                             matrices_dfs)

    
    matrices_merged.columns = ['start1','end1','start2','end2','sampleA_rep1','sampleA_rep2','sampleB_rep1','sampleB_rep2']
    
    # add chromosome fields
    matrices_merged['chr1'] = chromosome
    matrices_merged['chr2'] = chromosome

    matrices_merged = matrices_merged[
        ['chr1','start1','end1','chr2','start2','end2','sampleA_rep1','sampleA_rep2','sampleB_rep1','sampleB_rep2']
    ]


    save_path = (
        f"{output_path}/diffHic_input_chr{chromosome}_res{resolution}_k{spikes_k}.table"
    )

    matrices_merged.to_csv(save_path, sep='\t', index=False, header=True)
    return matrices_merged, save_path

def generate_diffHic_input_ds(
	raw_matrix,
	chromosome: str,
	resolution: int,
	downsample_factor: float,
    save_path: str ):
	# A condition 2 replicates full depth
	synthetic_matrix1=DifFracTion_utils.synthetic_datasets(raw_matrix,0.999999999999)
	synthetic_matrix1_2=DifFracTion_utils.synthetic_datasets(raw_matrix,0.999999999999)
	m1=DifFracTion_utils.matrix2longdf(synthetic_matrix1)
	m1_2=DifFracTion_utils.matrix2longdf(synthetic_matrix1_2)

	# B condition 2 replicates
	synthetic_matrix2=DifFracTion_utils.synthetic_datasets(raw_matrix,downsample_factor)
	synthetic_matrix2_2=DifFracTion_utils.synthetic_datasets(raw_matrix,downsample_factor)
	m2=DifFracTion_utils.matrix2longdf(synthetic_matrix2)
	m2_2=DifFracTion_utils.matrix2longdf(synthetic_matrix2_2)

	# Merge 
	m1    = m1.rename(columns={4: "sampleA_rep1"})
	m1_2  = m1_2.rename(columns={4: "sampleA_rep2"})
	m2    = m2.rename(columns={4: "sampleB_rep1"})
	m2_2  = m2_2.rename(columns={4: "sampleB_rep2"})
	matrices_dfs = [m1, m1_2, m2, m2_2]
	matrices_merged = reduce( lambda left, right:
							pd.merge(left,right,
									on=[0,1,2,3],
									how='inner').fillna(0),
								matrices_dfs)
		
	matrices_merged.columns = ['start1','end1','start2','end2','sampleA_rep1','sampleA_rep2','sampleB_rep1','sampleB_rep2']

	# add chromosome fields
	matrices_merged['chr1'] = chromosome
	matrices_merged['chr2'] = chromosome
	matrices_merged = matrices_merged[
		['chr1','start1','end1','chr2','start2','end2','sampleA_rep1','sampleA_rep2','sampleB_rep1','sampleB_rep2']
	]
	# Save
	save_path = f'{save_path}/diffHic_input_chr{chromosome}_res{resolution}_ds{downsample_factor}.table'
	matrices_merged.to_csv(save_path, sep='\t', index=False, header=True)
	return matrices_merged, save_path

def generate_multiHiCcompare_input(altered_matrix,
    raw_matrix,
    chromosome: str,
    resolution: int,
    spikes_k: float,
    output_path: str
):
    
    ###. Re do, now we have to save each file individually  
    # and just keep chromosome, start, end, and IF
    # make_hicexp(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4, groups = c(1, 1, 2, 2))
    # Takes care of buiolding the object internally (R)
    """
    altered_matrix : matrix with spike-ins applied
    raw_matrix     : original real matrix
    DifFracTion_utils  : module namespace (DifFracTion_utils) so function works in scripts too
    """

    def generate_individual_frame(matrix, resolution, col_name):
        m = DifFracTion_utils.matrix2longdf(matrix, resolution)  # upper triangle already
        m = m.rename(columns={4: col_name})
        m.columns = ['start1','end1','start2','end2',col_name]
        m['chr1'] = chromosome
        m['chr2'] = chromosome
        m = m[['chr1','start1','end1','chr2','start2','end2',col_name]]
        return m

    def save_individual_frame(m, resolution, spikes_k, condition_rep):
        save_path = (
            f"{output_path}/multiHiCcompare_input_chr{chromosome}_res{resolution}_k{spikes_k}_{condition_rep}.table"
        )
        m.to_csv(save_path, sep='\t', index=False, header=True)
        return save_path
    
    # A condition
    # Rep1
    m1 = generate_individual_frame(altered_matrix, resolution, "IF_A1")
    # Rep2
    altered_matrix_rep = DifFracTion_utils.synthetic_datasets(altered_matrix, 0.9999999) # this is because we already have spike-ins
    m1_2 = generate_individual_frame(altered_matrix_rep, resolution, "IF_A2")
    
    # B condition
    # Rep1
    synthetic = DifFracTion_utils.synthetic_datasets(raw_matrix, 1)
    m2 = generate_individual_frame(synthetic, resolution, "IF_B1")
    # Rep2
    synthetic_rep = DifFracTion_utils.synthetic_datasets(raw_matrix, 0.9999999)
    m2_2 = generate_individual_frame(synthetic_rep, resolution, "IF_B2")

    save_path_A1 = save_individual_frame(m1, resolution, spikes_k, "A1")
    save_path_A2 = save_individual_frame(m1_2, resolution, spikes_k, "A2")
    save_path_B1 = save_individual_frame(m2, resolution, spikes_k, "B1")
    save_path_B2 = save_individual_frame(m2_2, resolution, spikes_k, "B2")

    return save_path_A1, save_path_A2, save_path_B1, save_path_B2

def generate_multiHiCcompare_input_ds(
    raw_matrix,
    chromosome: str,
    resolution: int,
    downsample_factor: float,
    output_path: str
):
    
    ###. Re do, now we have to save each file individually  
    # and just keep chromosome, start, end, and IF
    # make_hicexp(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4, groups = c(1, 1, 2, 2))
    # Takes care of buiolding the object internally (R)
    """
    altered_matrix : matrix with spike-ins applied
    raw_matrix     : original real matrix
    DifFracTion_utils  : module namespace (DifFracTion_utils) so function works in scripts too
    """

    def generate_individual_frame(matrix, resolution, col_name):
        m = DifFracTion_utils.matrix2longdf(matrix, resolution)  # upper triangle already
        m = m.rename(columns={4: col_name})
        m.columns = ['start1','end1','start2','end2',col_name]
        m['chr1'] = chromosome
        m['chr2'] = chromosome
        m = m[['chr1','start1','end1','chr2','start2','end2',col_name]]
        return m

    def save_individual_frame(m, resolution, downsample_factor, condition_rep):
        save_path = (
            f"{output_path}/multiHiCcompare_input_chr{chromosome}_res{resolution}_ds{downsample_factor}_{condition_rep}.table"
        )
        m.to_csv(save_path, sep='\t', index=False, header=True)
        return save_path
    
    # A condition
    # Rep1
    m1 = DifFracTion_utils.synthetic_datasets(raw_matrix,0.9999999)
    m1 = generate_individual_frame(m1, resolution, "IF_A1")
    # Rep2
    m1_2 = DifFracTion_utils.synthetic_datasets(raw_matrix,0.9999999)
    m1_2 = generate_individual_frame(m1_2, resolution, "IF_A2")
    # B condition
    # Rep1
    m2 = DifFracTion_utils.synthetic_datasets(raw_matrix, downsample_factor) # this is because we already have spike-ins
    m2 = generate_individual_frame(m2, resolution, "IF_B1")
    m2_2 = DifFracTion_utils.synthetic_datasets(raw_matrix, downsample_factor)
    m2_2 = generate_individual_frame(m2_2, resolution, "IF_B2")

    save_path_A1 = save_individual_frame(m1, resolution, downsample_factor, "A1")
    save_path_A2 = save_individual_frame(m1_2, resolution, downsample_factor, "A2")
    save_path_B1 = save_individual_frame(m2, resolution, downsample_factor, "B1")
    save_path_B2 = save_individual_frame(m2_2, resolution, downsample_factor, "B2")

    return save_path_A1, save_path_A2, save_path_B1, save_path_B2

def generate_HiCcompare_input(altered_matrix,
    raw_matrix,
    chromosome: str,
    resolution: int,
    spikes_k: float,
    output_path: str
):
    """
    altered_matrix : matrix with spike-ins applied
    raw_matrix     : original real matrix
    DifFracTion_utils  : module namespace (DifFracTion_utils) so function works in scripts too
    """

    # convert to long format
    m1 = DifFracTion_utils.matrix2longdf(altered_matrix, resolution)  # upper triangle already
    synthetic = DifFracTion_utils.synthetic_datasets(raw_matrix, 1)
    m2 = DifFracTion_utils.matrix2longdf(synthetic, resolution)

    # merge IF1 / IF2
    matrices_merged = m1.merge(m2, on=[0,1,2,3], how='inner').fillna(0)
    matrices_merged.columns = ['start1','end1','start2','end2','IF1','IF2']

    # add chromosome fields
    matrices_merged['chr1'] = chromosome
    matrices_merged['chr2'] = chromosome

    # distance in bins
    matrices_merged['D'] = abs(matrices_merged['start1'] - matrices_merged['start2']) // resolution

    # M = log2 fold change
    matrices_merged['M'] = np.log2((matrices_merged['IF2'] + 1) / (matrices_merged['IF1'] + 1))

    matrices_merged = matrices_merged[
        ['chr1','start1','end1','chr2','start2','end2','IF1','IF2','D','M']
    ]


    save_path = (
        f"{output_path}/HiCcompare_input_chr{chromosome}_res{resolution}_k{spikes_k}.table"
    )

    matrices_merged.to_csv(save_path, sep='\t', index=False, header=True)

    return matrices_merged, save_path

def generate_HiCcompare_input_ds(
    raw_matrix,
    chromosome: str,
    resolution: int,
    downsample_factor: float,
    output_path: str
):
    """
    altered_matrix : matrix with spike-ins applied
    raw_matrix     : original real matrix
    DifFracTion_utils  : module namespace (DifFracTion_utils) so function works in scripts too
    """

    m1 = DifFracTion_utils.synthetic_datasets(raw_matrix,0.9999999)
    m1 = DifFracTion_utils.matrix2longdf(m1, resolution)  # upper triangle already

    m2 = DifFracTion_utils.synthetic_datasets(raw_matrix, downsample_factor)
    m2 = DifFracTion_utils.matrix2longdf(m2, resolution)

    # merge IF1 / IF2
    matrices_merged = m1.merge(m2, on=[0,1,2,3], how='inner').fillna(0)
    matrices_merged.columns = ['start1','end1','start2','end2','IF1','IF2']

    # add chromosome fields
    matrices_merged['chr1'] = chromosome
    matrices_merged['chr2'] = chromosome

    # distance in bins
    matrices_merged['D'] = abs(matrices_merged['start1'] - matrices_merged['start2']) // resolution

    # M = log2 fold change
    matrices_merged['M'] = np.log2((matrices_merged['IF2'] + 1) / (matrices_merged['IF1'] + 1))

    matrices_merged = matrices_merged[
        ['chr1','start1','end1','chr2','start2','end2','IF1','IF2','D','M']
    ]


    save_path = (
        f"{output_path}/HiCcompare_input_chr{chromosome}_res{resolution}_ds{downsample_factor}.table"
    )

    matrices_merged.to_csv(save_path, sep='\t', index=False, header=True)

    return matrices_merged, save_path

