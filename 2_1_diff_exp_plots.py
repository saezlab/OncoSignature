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

#----------------------------------- INPUT -----------------------------------#
parent_dir = 'results'
#-----------------------------------------------------------------------------#
usedirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir)
           if d.startswith('2_diff_exp')]

for dir_ in usedirs:
    # Loading the results from the differential expression analysis
    files = [f for f in os.listdir(dir_) if (f.endswith('_ttop.csv')
                                             and not f.startswith('sig'))]

    for f in files:
        fname = f.replace('_ttop.csv', '')
        name = ' vs. '.join(fname.split('vs'))

        df = pd.read_csv(os.path.join(dir_, f), index_col=0)

        # Volcano plot
        volcano(df['logFC'], -np.log10(df['P.Value']), title=name,
                filename=os.path.join(dir_, '%s.pdf' % fname))

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
