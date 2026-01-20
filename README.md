# DiFracTion Installation

DifFraction can be installed in the following ways

     git clone https://github.com/g2lab/DifFracTion.git
     cd DifFracTion
     pip install .


Data used in this study corresponds to cell line GM12878 and was retrieved from GSE63525 accession number:

     mkdir test_data/
     cd test_data/
     wget [https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63525&format=file&file=GSE63525%5FGM12878%5Fcombined%5Finterchromosomal%5Fcontact%5Fmatrices%2Etar%2Egz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Finterchromosomal%5Fcontact%5Fmatrices.tar.gz)
     

The DiFfracTion_Tutorial notebook offers a walkthough to the main DifFracTion functions:
     
     DifFracTion.alpha_normalization()
     DifFracTion.iterative_normalization()
     DifFracTion.identify_differential_interactions()

## Normalization

The inputs to alpha_normalization() and iterative_normalization() are two $n x n $ chromosome-level balanced Hi-C matrices (e.g., KR or any other balancing method), referred to as matrixA and matrixB, along with the length of the corresponding chromosome and the matrix resolution in base pairs.
     
     matrix_A_norm,matrix_B_norm,_,_ = DifFracTion.alpha_normalization(matrix_A,matrix_B,chromosome_length,resolution)

     matrix_A_norm_iter,matrix_B_norm_iter,_,_ = DifFracTion.iterative_normalization(matrix_A,matrix_B,chromosome_length,resolution)

Both functions output the fully normalized upper triangles of both matrices in plain-text format.

## DCIs

For identify_differential_interactions(), we strongly recommend using a high-performance computing (HPC) system, as the subsampling procedure generates 
$2N$ additional Hi-C matrices to construct the null distribution, which is later used to test for significance. Ideal results are observed with $N=1000$ subsampling iterations.

     dcis = DifFracTion.identify_differential_interactions(matrix_A_norm,matrix_B_norm,resolution,method='log2fc', subsample_n=20,log2_fc_cutoff=0,p_value_threshold=0.01,bayesian=False)


This function takes as input two fully normalized Hi-C contact matrices (matrix\_A\_norm and matrix\_B\_norm) at the same resolution, along with the matrix resolution in base pairs. Differential interactions can be quantified using a log2 fold-change or difference based operation (method='log2fc'). 

Interactions are reported if they pass a user-defined log2 fold-change threshold (log2_fc_cutoff) and a significance cutoff (p_value_threshold). 

Output format:

     bin1	bin2	count_1	count_2	log2_fc   distance	p_value	p_value_adj	neighbor_support
     0	34	249.8892  178.0     0.4893    3400000   0.0  0.0  True
     0	35	218.7481  146.0     0.5832    3500000	0.0  0.0  True

The bayesian flag controls whether a Bayesian framework is used for significance estimation.

When the Bayesian framework is enabled, interactions are filtered based on the maximum a posteriori (MAP) probability rather than the log2 fold-change statistic.

Output format:

     bin1	bin2	count_1	count_2	log2_fc	post_mean_log2fc	post_map_log2fc	distance	p_value	p_value_adj	neighbor_support
     0	34	249.8892	178.0	0.4893	25312.7513	0.4594	3400000	0.0	0.0	True
     0	35	218.74814	146.0	0.5832	135823.0750	1.8127	3500000	0.0	0.0	True

The *neighbor_support* column is a boolean flag indicating whether a given interaction has at least one statistically significant neighboring interaction. Although both `True` and `False` values are reported, only interactions with `neighbor_support = True` are retained in the final results.
