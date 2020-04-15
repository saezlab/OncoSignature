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

#----------------------------------- INPUT -----------------------------------#
parent_dir = 'results/perseus'
params = dict([('Cell_Lines', dict([('xlim', (-2.25, 2.25)),
                                    ('ylim', (-0.25, 8.25))])),
               ('ExVivo', dict([('xlim', (-2.25, 2.25)),
                                ('ylim', (-0.25, 20))]))])
#-----------------------------------------------------------------------------#
usedirs = os.listdir(parent_dir)

for dir_ in usedirs:
    fulldir = os.path.join(parent_dir, dir_)
    # Loading the results from the differential expression analysis
    dex_files = [f for f in os.listdir(fulldir) if f.startswith('DEX_')]

    for dexf in dex_files:
        scvf = dexf.replace('DEX', 'Scurve')
        name = dexf.lstrip('DEX_').rstrip('.txt')

        df = pd.read_csv(os.path.join(fulldir, dexf), sep='\t')
        df.set_index('Fusion', drop=True, inplace=True)
        df.index.name = None
        df = df.loc[:, ['Significant', '(-)log10(P)', 'Difference']]
        is_sig = np.array([i == '+' for i in df['Significant'].values])

        curve = pd.read_csv(os.path.join(fulldir, scvf), sep='\t')

        #print(name)
        #print(df['Difference'].max(), df['Difference'].min())
        #print(df['(-)log10(P)'].max())

        fig, ax = plt.subplots()
        # Significant points
        ax.scatter(df.loc[is_sig, 'Difference'],
                   df.loc[is_sig, '(-)log10(P)'],
                   c='C1', s=5, alpha=0.3, label='Significant')
        # Non-significant points
        ax.scatter(df.loc[~is_sig, 'Difference'],
                   df.loc[~is_sig, '(-)log10(P)'],
                   c='C0', s=5, alpha=0.3, label='Non-significant')
        # S-curve
        ax.plot(curve.x, curve.y, 'k--', alpha=0.5)

        ax.set_xlim(*params[dir_]['xlim'])
        ax.set_ylim(*params[dir_]['ylim'])

        ax.set_title(name)
        ax.set_xlabel(r'$\log_2$(FC)')
        ax.set_ylabel(r'$-\log_{10}$(p-val)')

        fig.savefig(os.path.join(fulldir, 'volcano_%s.svg' % name))
