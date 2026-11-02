# DiFracTion Installation
```bash
     git clone https://github.com/g2lab/DifFracTion.git

     cd DifFracTion

     conda env create -f environment.yml
     conda activate diffraction
     pip install .
```

## Validation Data

Pre-processed chr19 Hi-C matrices at 10kb resolution are provided for quick validation of the DifFracTion pipeline. These files are included in the repository under `generated_data/` and can be used directly as inputs without downloading the full `.hic` files.

| File | Cell type / Tissue |
|------|--------------------|
| `generated_data/K562_chr19_10kb.npz` | K562 (chronic myelogenous leukemia) |
| `generated_data/GM12878_chr19_10kb.npz` | GM12878 (lymphoblastoid) |
| `generated_data/TransverseColon_chr19_10kb.npz` | Transverse colon |
| `generated_data/GastrocnemiusMedialis_chr19_10kb.npz` | Gastrocnemius medialis (skeletal muscle) |

Data used in this study corresponds to cell line GM12878 and was retrieved from GSE63525 accession number:

```bash
     mkdir test_data/
     cd test_data/

     Read counts 1,844,107,778
     wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_insitu_primary%2Breplicate_combined_30.hic
     mv GSE63525_GM12878_insitu_primary+replicate_combined_30.hic GM12878-HRC.hic

     Read counts 735,043,093
     wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_IMR90_combined_30.hic 
     mv GSE63525_IMR90_combined_30.hic IMR90-HRC.hic
     
     Read counts 641,402,880
     wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_K562_combined_30.hic
     mv GSE63525_K562_combined_30.hic K562-HRC.hic

```

The DiFfracTion_Tutorial notebook offers a walkthough to the main DifFracTion functions:

```python
     DifFracTion.alpha_normalization()
     DifFracTion.iterative_normalization()
     DifFracTion.identify_differential_interactions()
```
## Normalization

The inputs to alpha_normalization() and iterative_normalization() are two $n x n $ chromosome-level balanced Hi-C matrices (e.g., KR or any other balancing method), referred to as matrixA and matrixB, along with the length of the corresponding chromosome and the matrix resolution in base pairs.
     
     norm_matrix_A, norm_matrix_B = DifFracTion.alpha_normalization(matrix_A,matrix_B,resolution,metric='median',desired_alpha=1.08,filter=True)

        
     norm_matrix_A_iterative, norm_matrix_B_iterative = DifFracTion.iterative_normalization(matrix_A, matrix_B, 
                                                       resolution, metric='median', weight_factor=-1.08)
        

Both functions output the fully normalized upper triangles of both matrices in plain-text format.

Please refer to the function vignettes for parameter details.

## DCIs

For identify_differential_interactions(), we strongly recommend using a high-performance computing (HPC) system on high resolutions (10kb), as the subsampling procedure generates 
$2N$ additional Hi-C matrices to construct the null distribution, which is later used to test for significance. Ideal results are observed with $N=1000$ subsampling iterations. Refer to our paper supplementary materials for execution time and memory details.

     alpha_padjusted_interactions = DifFracTion.identify_differential_interactions(norm_matrix_A, norm_matrix_B,
                                   resolution, filter_by='padjusted', pvalue_cutoff=0.01,
                                   log2fc_cutoff=0, permutations=100, bayesian=True)

       
This function takes as input two fully normalized Hi-C contact matrices (matrix\_A\_norm and matrix\_B\_norm) at the same resolution, along with the matrix resolution in base pairs. 


Interactions are reported if they pass a user-defined log2 fold-change threshold (log2_fc_cutoff) and a significance cutoff (p_value_threshold). 

Output format:

     row	col	condition_A_counts	condition_B_counts	log2fc	pvalue	adjusted_pvalue	distance	neighbor_support         
     12819	12831	281.6880798339844	53.76190001719214	2.2957815112963003	0.0	0.0	120000.0	True
     12820	12832	140.6009521484375	14.022085522013567	2.931904096404026	0.0	0.0	120000.0	True


The *neighbor_support* column is a boolean flag indicating whether a given interaction has at least one statistically significant neighboring interaction. Although both `True` and `False` values are reported, only interactions with `neighbor_support = True` are retained in the final results.

# Paper results

To replicate any figure included on the paper refer to the corresponding folder and follow the corresponding README instructions.
