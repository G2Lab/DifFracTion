import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import seaborn as sns
from . import utils as DifFracTion_utils


def get_TP_list(search_path,method,p_val,resolution):
	TP_files = [f for f in search_path.rglob(f"spikein_and_neighbors_coordinates*txt")
		   if method in str(f.parent) and str(p_val) in str(f.name) and f"_{resolution}" in str(f.name)]
	TP_files = [f for f in TP_files if f.is_file()]
	TP_files.sort()
	df_list = []
	for file in TP_files:
		name = file.parent.name
		parts = name.split('_')
		if method == 'alpha_based':
			chrstr,pvalstr,restr,dsstr,k,_,_ = parts
		elif method == 'cycle':
			chrstr,pvalstr,restr,dsstr,k,_ = parts
		df = pd.read_csv(file, sep='\t',
			    header=None,index_col=False,names=['bin1','bin2','count'])
		df['chromosome'] = chrstr
		df['p_value_threshold'] = float(pvalstr) #adjusted pval
		df['resolution'] = int(restr)
		df['proportion'] = float(dsstr)
		df['k'] = int(k)
		df_list.append(df)
		
	return df_list

def get_significant_contacts(search_path,method,p_val,resolution,option='all'):
     if option=='all':
          string_search = "significant_contacts*tsv"
     elif option=='supported':
          string_search = "supported*_contacts*tsv"
     else:
          raise ValueError("option must be either 'all' or 'supported'")
	
     significant_files = [f for f in search_path.rglob(string_search)
		   if method in str(f.parent) and str(p_val) in str(f.name) and f"_{resolution}" in str(f.name)]
     significant_files = [f for f in significant_files if f.is_file()]
	
     significant_files.sort()
     df_list = []
     for file in significant_files:
          name = file.parent.name
          parts = name.split('_')
          if method == 'alpha_based':
               chrstr,pvalstr,restr,dsstr,k,_,_ = parts
          elif method == 'cycle':
               chrstr,pvalstr,restr,dsstr,k,_ = parts
          df = pd.read_csv(file, sep='\t')
          df['chromosome'] = chrstr
          df['p_value_threshold'] = float(pvalstr) #adjusted pval
          df['resolution'] = int(restr)
          df['proportion'] = float(dsstr)
          df['k'] = int(k)
          df_list.append(df)
		
     return df_list

def get_number_tests(hic_file,chromosome,resolution):
     chromosome_length = DifFracTion_utils.get_chromosome_length(hic_file, chromosome)
     matrixZoom=DifFracTion_utils.prepare_MatrixZoomData(hic_file, chromosome, resolution)
     raw_matrix = DifFracTion_utils.rawMatrix(matrixZoom,chromosome_length)
     raw_matrix = np.triu(raw_matrix) # because only upper triangle is used
     df_m=DifFracTion_utils.matrix2df(raw_matrix)
     df_m['d'] = (df_m[1] - df_m[0])*resolution
     df_m= df_m[df_m[2] > 0] # counts > 0
     df_m = df_m[df_m['d'] > 0]
     n_tests = len(df_m)
     return n_tests

def get_confusion_matrix(tp_df_list, sig_df_list,hic_file):
     sig_lookup = {}
     full_df = pd.DataFrame()
	
	# What we found to be significant
     for df in sig_df_list:
          key = (
            df['chromosome'].iloc[0],
            df['p_value_threshold'].iloc[0],
            df['resolution'].iloc[0],
            df['proportion'].iloc[0],
            df['k'].iloc[0],
        )
          sig_lookup[key] = df

     # What we should have found to be significant
     for tp_df_c in tp_df_list:
          chromosome = tp_df_c['chromosome'].iloc[0]
          pvalue    = tp_df_c['p_value_threshold'].iloc[0]
          resolution = tp_df_c['resolution'].iloc[0]
          proportion = tp_df_c['proportion'].iloc[0]
          k         = tp_df_c['k'].iloc[0]   
          
          n_tests= get_number_tests(hic_file,chromosome,resolution)

          key = (chromosome, pvalue, resolution, proportion, k)  
		
          # Get matching sig dataframe if it exists
          sig_df_c = sig_lookup.get(key, None) #what the tool found
    
          K_TP = set(zip(tp_df_c['bin1'], tp_df_c['bin2'])) # Total true positives inserted (spike-ins)
          N = set(zip(sig_df_c['bin1'], sig_df_c['bin2'])) if sig_df_c is not None else set() # Total significant contacts found by the tool
    
          if sig_df_c is not None:
               TP = len(K_TP & N) # True positives found by the tool
               FN = len(K_TP - N) # False negatives, what the tool missed
               FP = len(N - K_TP) # False positives, what the tool incorrectly found
			# True negatives, what the tool correctly identified as not significant
               TN = n_tests -  TP - FP - FN
          else:
               TP = 0
               FN = len(K_TP)
               FP = 0
               TN = 0
          
          full_df = pd.concat([full_df, pd.DataFrame({
               'chromosome': [chromosome],
               'p_value_threshold': [pvalue],
               'resolution': [resolution],
               'proportion': [proportion],
               'k': [k],
               'TP': [TP],
               'FP': [FP],
               'TN': [TN],
               'FN': [FN],
               'Total_TP': [len(K_TP)]
          })], ignore_index=True)
          
     return full_df
    
def get_metrics(confusion_matrix):
	confusion_matrix['sensitivity'] = confusion_matrix['TP'] / (confusion_matrix['TP'] + confusion_matrix['FN'])
	confusion_matrix['specificity'] = confusion_matrix['TN'] / (confusion_matrix['TN'] + confusion_matrix['FP'])
	confusion_matrix['precision'] = confusion_matrix['TP'] / (confusion_matrix['TP'] + confusion_matrix['FP'])
	confusion_matrix['accuracy'] = (confusion_matrix['TP'] + confusion_matrix['TN']) / (confusion_matrix['TP'] + confusion_matrix['TN'] + confusion_matrix['FP'] + confusion_matrix['FN'])
	confusion_matrix = confusion_matrix.replace([np.inf, -np.inf], np.nan).fillna(0)
	return confusion_matrix

def generate_latex_table(metrics_df, chromosome):
     metrics_df = metrics_df[metrics_df['chromosome']==chromosome]
     metrics_df = metrics_df[metrics_df['k'].isin([2,3,4])]
     latex_str = metrics_df.to_latex(
    index=False,
    float_format="%.2f",
    column_format="lrrrrrrrrrrrrr"
     )
     return latex_str

def plot_results(results, resolution, column):
    fig =plt.figure(figsize=(5,3), dpi=400)

    colors = ["#E9C46A", "#2a9d8f", "#e76f51", "#669bbc", "#f4a261"]
    sns.set_theme()
    sns.set_style("ticks")
    
    data = results[results['resolution'] == resolution]
    p_values = data['p_value_threshold'].unique()
    pval=p_values[0]
    data = data[data['k'].isin([2,3,4])]
    sns.barplot(
          data=data,
          x='k',
          y=column,
          hue='proportion',
          palette=colors,
          errorbar="se",
          linewidth=0.3,
          edgecolor='black',
          capsize=0.15,err_kws={ "linewidth": 0.5},
     )

    plt.ylabel(column)
    plt.xlabel("k")
    plt.title(f"Resolution: {resolution} - p-value: {pval}")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    sns.despine()
    plt.tight_layout()
    return fig