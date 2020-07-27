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
# Dose-response combination therapy
# =================================
#
# This script generates the plots regarding the data normalization.

import os
from itertools import product

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from data_tools.plots import cluster_hmap

data_dir = 'data'
out_dir = 'results'
# First is raw (before), second normalized (after) (not relevant tho)
prefixes = ['raw', 'norm_data']
suffixes = ['cl', 'ex']
# Here be the batches annotated
annotf = 'response.csv'

batches = pd.read_csv(os.path.join(data_dir, annotf),
                      index_col=0).batch.to_dict()
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for p, s in product(prefixes, suffixes):
    path = os.path.join(data_dir, '%s_%s.csv' % (p, s))

    df = pd.read_csv(path, index_col=0)
    num_batch = len(set(batches[c] for c in df.columns))

    corr = df.corr()
    #print(corr.min().min())
    dtype = 'Raw' if p == 'raw' else 'Normalized'
    extype = 'cell line' if s == 'cl' else 'ex-vivo'
    # Min correlation overall is 0.434...
    fig = cluster_hmap(corr.values, xlabels=corr.columns, ylabels=corr.index,
                       hmap_kwargs={'vmin': 0.4, 'vmax': 1},
                       title='%s data of %s samples' % (dtype, extype), figsize=(9, 9))
    ax = fig.axes[1]

    for x in ax.get_xticklabels():
        idx = batches[x.get_text()]
        x.set_color(colors[idx if idx <= 4 else idx - 4])

    for y in ax.get_yticklabels():
        idx = batches[y.get_text()]
        y.set_color(colors[idx if idx <= 4 else idx - 4])

    fig.savefig(os.path.join(out_dir, '%s_%s_cor_hmap.pdf') % (p, s))
