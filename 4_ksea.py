# Copyright (C) 2019 Nicolàs Palacio
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
# Kinase-substrate enrichment analysis
# ====================================
#
# This script computes the kinase enrichment based on the differential
# phosphoproteomics data

import os

import pandas as pd
import numpy as np

import kinact
from data_tools.plots import volcano
from data_tools.databases import op_kinase_substrate
from data_tools.databases import up_map

#----------------------------------- INPUT -----------------------------------#
data_dir = 'data'
dex_dir = 'results/2_diff_exp'
out_dir = 'results/4_ksea'

# If exists, reads it, otherwise writes it
adj_file = 'ks_adj_281119.csv'
#-----------------------------------------------------------------------------#

# Loading kinase-substrate network
adj_path = os.path.join(data_dir, adj_file)

if os.path.exists(adj_path):
    ks_adj = pd.read_csv(adj_path)

# Download the kinase-substrate interactions table from OmniPath
else:
    ks_tab = op_kinase_substrate(gsymbols=False, incl_phosphatases=True)

    # Join phosphosite identifiers
    ks_tab['psite'] = ['{}_{}{}'.format(*r[1].values[1:4])
                       for r in ks_tab.iterrows()]

    # Building adjacency matrix
    unique_kinases = ks_tab['enzyme'].sort_values().unique()
    unique_psites = ks_tab['psite'].sort_values().unique()

    ks_adj = pd.DataFrame(index=unique_psites, columns=unique_kinases)

    for k in unique_kinases:
        targets = ks_tab.loc[ks_tab['enzyme'] == k, 'psite']
        ks_adj.loc[targets, k] = 1

    ks_adj.to_csv(adj_path)

ttop_files = [f for f in os.listdir(dex_dir)
              if (f.endswith('_ttop.csv') and not f.startswith('sig_'))]
ttops = dict([f.replace('_ttop.csv', ''), pd.read_csv(os.path.join(dex_dir, f),
                                                      index_col=0,
                                                      usecols=[0, 1])]
              for f in ttop_files)

results = {}

for k, v in ttops.items():
    aux = v.dropna().copy()
    aux.index = [i.split('__')[0] for i in aux.index]

    score, pval = kinact.ksea.ksea_mean(aux, ks_adj, minimum_set_size=5)
    results[k] = pd.DataFrame({'score':score, 'pval':pval})

pval_thr = 0.05

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

for k, v in results.items():
    title = 'KSEA {} vs. {}'.format(*k.split('vs'))
    print('Significant kinases of', title)

    idmap = up_map([i.split('_')[0] for i in v.index])
    v.index = [i.replace(i.split('_')[0],
                         idmap.loc[idmap['ACC'] == i.split('_')[0],
                                   'GENENAME'].values[0]) for i in v.index]

    print(v[v['pval'] <= pval_thr])

    plot = volcano(v['score'], -np.log10(v['pval']), thr_fc=1, thr_pval=pval_thr,
                   title=title)
    ax = plot.gca()
    ax.set_xlabel('KSEA score')
    plot.savefig(os.path.join(out_dir, 'ksea_%s.pdf' %k))

    v.to_csv(os.path.join(out_dir, '%s_KSEA_results.csv' %k))
