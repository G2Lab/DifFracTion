import hicstraw 
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from . import utils as DifFracTion_utils
import seaborn as sns

upper_limit = 7_000_000  # 7Mb

def load_significant_contacts(search_path, method, p_val, resolution, option='all'):
    '''Returns a list of dictionaries with significant contacts dataframes and parameters used.'''
    if option == 'all':
        string_search = "significant_contacts*tsv"
    elif option == 'supported':
        string_search = "supported*_contacts*tsv"
    else:
        raise ValueError("option must be either 'all' or 'supported'")

    significant_files = [
        f for f in search_path.rglob(string_search)
        if method in str(f.parent)
        and str(p_val) in str(f.name)
        and str(resolution) in str(f.name)
        and f.is_file()
    ]

    significant_files.sort()
    records = []

    for file in significant_files:
        name = file.parent.name
        parts = name.split('_')

        if method == 'alpha_based':
            chrstr, pvalstr, restr, dsstr, _, _ = parts
        elif method == 'cycle':
            chrstr, pvalstr, restr, dsstr, _ = parts
        else:
            raise ValueError("Unknown method")

        df = pd.read_csv(file, sep='\t')

        records.append({
            'chromosome': chrstr,
            'p_value_threshold': float(pvalstr),
            'resolution': int(restr),
            'proportion': float(dsstr),
            'df': df
        })

    return records

def get_significant_contacts_fpr(
    search_path, method, p_val, resolution, hic_file,option='all'):
    '''Returns a data frame generated from a list of dictionaries with FPR calculations.'''
    records = load_significant_contacts(
        search_path, method, p_val, resolution, option
    )
    rows = []
    for rec in records: #iterate over the dictionaries of the list
        
        #rec is the individual dictionary == records[i]
        total_tests = get_number_tests(
            hic_file,
            rec['chromosome'],
            resolution)
        print(rec['chromosome'])
        print(rec['resolution'])
        print(total_tests)

        FP = len(rec['df'])  
        FPR = FP / total_tests if total_tests > 0 else 0.0

        #This is a list of dictionaries that will be converted to a dataframe
        rows.append({
            'chromosome': rec['chromosome'],
            'p_value_threshold': rec['p_value_threshold'],
            'resolution': rec['resolution'],
            'proportion': rec['proportion'],
            'total_tests': total_tests,
            'FP': FP,
            'FPR': FPR
        })

    return pd.DataFrame(rows)

def get_number_tests(hic_file,chromosome,resolution):
     chromosome_length = DifFracTion_utils.get_chromosome_length(hic_file, chromosome)
     matrixZoom=DifFracTion_utils.prepare_MatrixZoomData(hic_file, chromosome, resolution)
     raw_matrix = DifFracTion_utils.rawMatrix(matrixZoom,chromosome_length)
     raw_matrix = np.triu(raw_matrix)
     df_m=DifFracTion_utils.matrix2df(raw_matrix)
     df_m['d'] = (df_m[1] - df_m[0])*resolution
     df_m= df_m[df_m[2] > 0] # counts > 0
     df_m = df_m[df_m['d'] > 0]
     n_tests = len(df_m)
     return n_tests

def plot_results(results, method,p_val,column='FPR'):
    fig =plt.figure(figsize=(5,3), dpi=400)
    results = results[results['method']==method]
    results = results[results['p_value_threshold']==p_val]
    colors = ["#E9C46A", "#2a9d8f", "#e76f51", "#669bbc", "#f4a261"]
    sns.set_theme()
    sns.set_style("ticks")
    sns.barplot(
          data=results,
          x='resolution',
          y=column,
          hue='proportion',
          palette=colors,
          errorbar="se",
          linewidth=0.3,
          edgecolor='black',
          capsize=0.15,err_kws={ "linewidth": 0.5},
     )

    plt.ylabel(column)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    sns.despine()
    plt.tight_layout()
    return fig

def generate_latex_table(metrics_df, method, p_val):
     
     metrics_df = metrics_df[metrics_df['method']==method]
     metrics_df = metrics_df[metrics_df['p_value_threshold']==p_val]
     metrics_df = metrics_df[['chromosome', 'p_value_threshold', 'resolution', 'proportion', 'total_tests', 'FP', 'FPR']]
     latex_str = metrics_df.to_latex(
    	index=False,
    	float_format="%.5f",
    	column_format="lrrrrrr" #left align for chromosome, right align for numbers
     )
     return latex_str