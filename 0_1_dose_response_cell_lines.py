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
# This script generates the plots for the survival experiments with
# combination therapy.

import os

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

#----------------------------------- INPUT -----------------------------------#
data_dir = 'data'
res_dir = 'results'
#-----------------------------------------------------------------------------#

usefiles = [f for f in os.listdir(data_dir) if f.startswith('survival_combi')]
for f in usefiles:
    data = pd.read_csv(os.path.join(data_dir, f))
    name = f.split('_')[-1].rstrip('.csv')

    com_doses = sorted(data['MK2206 dose (nM)'].unique())[::-1]

    cmap = matplotlib.cm.get_cmap('cool')
    colors = list(map(cmap,
                      np.linspace(1, 0, len(com_doses))))

    rep_cols = [i for i in data.columns if i.startswith('Raw')]

    control_row = [(r[1]['Selinexor dose (nM)'] == 0
                    and r[1]['MK2206 dose (nM)'] == 0) for r in data.iterrows()]

    ref = data.loc[control_row, rep_cols].values

    norm_data = data.copy()
    norm_data[rep_cols] = norm_data[rep_cols] / ref

    fig, ax = plt.subplots()

    for i, d in enumerate(com_doses):
        subdf = norm_data.loc[norm_data['MK2206 dose (nM)'] == d, :]
        subdf.set_index('Selinexor dose (nM)', drop=True, inplace=True)
        subdf = subdf[rep_cols]

        ax.errorbar(subdf.index, subdf.mean(axis=1), yerr=subdf.std(axis=1),
                    c=colors[i], label='%.1f nM MK2206' %d)
        ax.scatter(np.repeat(subdf.index.values, subdf.shape[1]),
                   subdf.values.reshape(-1), c=colors[i], alpha=0.5,
                   marker='o')

    ax.set_xscale('log')
    ax.legend(loc=0)
    ax.set_title(name)

    fig.savefig(os.path.join(res_dir, '%s_combi_dose_response.png' %name))
