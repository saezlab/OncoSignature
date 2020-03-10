# Copyright (C) 2020 Nicolàs Palacio
#
# Contact: nicolas.palacio@bioquant.uni-heidelberg.de
#
# GNU-GLPv3:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# A full copy of the GNU General Public License can be found on
# http://www.gnu.org/licenses/.
#
# Differential expression plots
# =============================
#
# This script generates the volcano plots to visualize the differential
# expression results

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from data_tools.plots import volcano
from data_tools.iterables import chunk_this
from data_tools.databases import up_map

#----------------------------------- INPUT -----------------------------------#
parent_dir = 'results'
#-----------------------------------------------------------------------------#
usedirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir)
           if d.startswith('2_diff_exp')]

# Parameters to plot arrows between points and labels
adj_txt_kwargs = dict(arrowprops=dict(zorder=4.5, color='k', arrowstyle='-'))

# Downloading mapping table UniProt AC to GeneSymbol
somefile = [f for f  in os.listdir(usedirs[0])
            if (f.endswith('_ttop.csv') and not f.startswith('sig'))][0]
somedf = pd.read_csv(os.path.join(usedirs[0], somefile), index_col=0)
ids = list(set([i.split('_')[0] for i in somedf.index]))

mapping = pd.Series(index=ids, name='GENENAME')
# We chunk the requests in 1000-sized ids because server rejects otherwise
for p in chunk_this(ids, 1000):
    aux = up_map(p)
    mapping[aux['ACC']] = aux['GENENAME']

pd.isna(mapping).sum()
# If no mapping available, go back to UniProt AC
no_map = pd.isna(mapping)
mapping[no_map] = mapping.index[no_map].tolist()

for dir_ in usedirs:
    # Loading the results from the differential expression analysis
    files = [f for f in os.listdir(dir_) if (f.endswith('_ttop.csv')
                                             and not f.startswith('sig'))]

    for f in files:
        fname = f.replace('_ttop.csv', '')
        name = ' vs. '.join(fname.split('vs'))

        df = pd.read_csv(os.path.join(dir_, f), index_col=0)

        # Converting ids to genename
        upids, residues = zip(*[i.split('_')[:2] for i in df.index])
        mapped = [mapping[i] for i in upids]
        df.index = ['_'.join([str(a), b]) for a, b in zip(mapped, residues)]

        # Volcano plot
        volcano(df['logFC'], -np.log10(df['P.Value']), title=name,
                labels=df.index.values,
                filename=os.path.join(dir_, '%s.pdf' % fname),
                adj_txt_kwargs=adj_txt_kwargs)

        # p-value histogram
        fig, ax = plt.subplots()

        ax.hist(df['P.Value'].dropna(), bins=100)
        ax.set_xlim(0, 1)
        ax.set_xlabel('P.Value')
        ax.set_ylabel('Frequency')
        ax.set_title('%s : p-value histogram' % name)

        fig.tight_layout()

        fig.savefig(os.path.join(dir_, '%s_pval_hist.pdf' % fname))
        plt.close('all')
