import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from data_tools.plots import venn


os.chdir(#'/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/2+2_AML_patients/
	'results/3_gsea')

comp = ('R2_R1', 'NR2_NR1') # NOTE: Order matters for Venn diagram

df = pd.read_csv('results_gsea_viper.csv')
df['dif'] = abs(df[comp[0]] - df[comp[1]])
df.sort_values(by='dif', ascending=False, inplace=True)
df.head()

# SUMMARY HEATMAP
n = 50

subdf = df.loc[:, comp].iloc[:n, :]
subdf.index = df.iloc[:n, 0]
subdf.index.name = None

fig, ax = plt.subplots(figsize=(12, 15))
im = ax.imshow(subdf.values, interpolation='none', cmap='bwr',
               aspect='auto', origin='lower')

for i in range(n):
    ax.text(0, i, np.round(subdf.iloc[i, :][comp[0]], 3), fontdict={'ha':'center'})

for i in range(n):
    ax.text(1, i, np.round(subdf.iloc[i, :][comp[1]], 3), fontdict={'ha':'center'})

ax.set_xticks(range(2))
ax.set_xticklabels(subdf.columns)

ax.set_yticks(range(n))
ax.set_yticklabels(subdf.index)

fig.colorbar(im)

fig.tight_layout()
fig.savefig('top_dif_viper_scores.pdf')

# VENN DIAGRAM
r_up = set([df.iloc[i, 0] for i in range(len(df)) if df.iloc[i, :][comp[0]] > 0])
r_dw = set([df.iloc[i, 0] for i in range(len(df)) if df.iloc[i, :][comp[0]] < 0])
nr_up = set([df.iloc[i, 0] for i in range(len(df)) if df.iloc[i, :][comp[1]] > 0])
nr_dw = set([df.iloc[i, 0] for i in range(len(df)) if df.iloc[i, :][comp[1]] < 0])

venn([r_up, r_dw, nr_up, nr_dw], labels=['R up', 'R down', 'NR up', 'NR down'],
     filename='venn_viper.pdf')
