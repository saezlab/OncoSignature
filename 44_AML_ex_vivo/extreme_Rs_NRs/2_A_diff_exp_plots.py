import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from data_tools.plots import volcano

# Defining current working directory
os.chdir(#'/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/44_AML_ex_vivo/
	'results/2_diff_exp/')

# Loading the results from the differential expression analysis
files = [f for f in os.listdir(os.getcwd()) if f.endswith('_ttop.csv')]

for f in files:
    name = ' vs. '.join(f.split('_')[0:2])

    df = pd.read_csv(f, index_col=0)

    # Volcano plot
    volcano(df['logFC'], -np.log10(df['P.Value']), title=name,
            filename='%s_%s.pdf' %tuple(f.split('_')[0:2]))

    # p-value histogram
    fig, ax = plt.subplots()

    ax.hist(df['P.Value'].dropna(), bins=100)
    ax.set_xlim(0, 1)
    ax.set_xlabel('P.Value')
    ax.set_ylabel('Frequency')
    ax.set_title('%s : p-value histogram' %name)

    fig.tight_layout()

    fig.savefig('%s_%s_pval_hist.pdf' %tuple(f.split('_')[0:2]))
