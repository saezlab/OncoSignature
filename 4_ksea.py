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
from data_tools.databases import up_map

#----------------------------------- INPUT -----------------------------------#
data_dir = 'data'
p_dex_dir = 'results'
p_out_dir = 'results/4_ksea'

ks_file = 'opath_ptms_filtered_100320.csv'

pval_thr = 0.05
# Arguments for labels
adj_txt_kwargs = dict(arrowprops=dict(zorder=4.5, color='k', arrowstyle='-'))
#-----------------------------------------------------------------------------#

# Loading kinase-substrate network
ks_path = os.path.join(data_dir, ks_file)

ks_net = pd.read_csv(ks_path)

# Join phosphosite identifiers
ks_net['psite'] = ['{}_{}{}'.format(*r[1].values[[1, 4, 5]])
                   for r in ks_net.iterrows()]

# Building adjacency matrix
unique_kinases = ks_net['enzyme'].sort_values().unique()
unique_psites = ks_net['psite'].sort_values().unique()

ks_adj = pd.DataFrame(index=unique_psites, columns=unique_kinases)

for k in unique_kinases:
    targets = ks_net.loc[ks_net['enzyme'] == k, 'psite']
    ks_adj.loc[targets, k] = 1

ks_adj.to_csv(ks_path.replace('.csv', '_adj.csv'))

usedirs = [os.path.join(p_dex_dir, d) for d in os.listdir(p_dex_dir)
           if d.startswith('2_diff_exp')]

for dex_dir in usedirs:
    ttop_files = [f for f in os.listdir(dex_dir)
                  if (f.endswith('_ttop.csv') and not f.startswith('sig_'))]
    ttops = dict([f.replace('_ttop.csv', ''),
                  pd.read_csv(os.path.join(dex_dir, f), index_col=0,
                              usecols=[0, 1])] for f in ttop_files)

    results = {}

    for k, v in ttops.items():
        aux = v.dropna().copy()
        aux.index = [i.split('__')[0] for i in aux.index]

        score, pval = kinact.ksea.ksea_mean(aux, ks_adj, minimum_set_size=3)
        results[k] = pd.DataFrame({'score':score, 'pval':pval})

    out_dir = os.path.join(p_out_dir, dex_dir.split('_')[-1])

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

        plot = volcano(v['score'], -np.log10(v['pval']), thr_fc=1,
                       thr_pval=pval_thr, title=title, labels=v.index.tolist(),
                       adj_txt_kwargs=adj_txt_kwargs, maxlabels=50)
        ax = plot.gca()
        ax.set_xlabel('KSEA score')
        plot.savefig(os.path.join(out_dir, 'ksea_%s.pdf' %k))

        v.to_csv(os.path.join(out_dir, '%s_KSEA_results.csv' %k))
