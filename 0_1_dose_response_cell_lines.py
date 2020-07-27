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
from itertools import product

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#----------------------------------- INPUT -----------------------------------#
data_dir = 'data'
res_dir = 'results'
#-----------------------------------------------------------------------------#

usefiles = [f for f in os.listdir(data_dir) if f.startswith('survival_combi')]

for f in usefiles:
    data = pd.read_csv(os.path.join(data_dir, f))
    name = f.split('_')[-1].rstrip('.csv')

    # Unique dose list
    treat_doses = sorted(data['Selinexor dose (nM)'].unique())
    com_doses = sorted(data['MK2206 dose (nM)'].unique())

    # Sequential colors corresponding to each dose
    cmap = matplotlib.cm.get_cmap('cividis')
    colors = list(map(cmap,
                      np.linspace(1, 0, len(com_doses))))

    rep_cols = [i for i in data.columns if i.startswith('Raw')]

    # Reference values for normalization (0nM/0nM)
    control_row = [(r[1]['Selinexor dose (nM)'] == 0
                    and r[1]['MK2206 dose (nM)'] == 0) for r in data.iterrows()]
    ref = data.loc[control_row, rep_cols].values

    norm_data = data.copy()
    norm_data[rep_cols] = norm_data[rep_cols] / ref
    norm_data['Average'] = norm_data[rep_cols].mean(axis=1)
    norm_data['stdev'] = norm_data[rep_cols].std(axis=1)

    # Dose-response curves
    fig, ax = plt.subplots()

    for i, d in enumerate(com_doses):
        subdf = norm_data.loc[norm_data['MK2206 dose (nM)'] == d, :]
        subdf.set_index('Selinexor dose (nM)', drop=True, inplace=True)
        #subdf = subdf[rep_cols]

        ax.errorbar(subdf.index, subdf['Average'], yerr=subdf['stdev'],
                    c=colors[i], label='%.1f nM MK2206' %d)
        ax.scatter(np.repeat(subdf.index.values, len(rep_cols)),
                   subdf[rep_cols].values.reshape(-1), c=colors[i], alpha=0.5,
                   marker='o')

    ax.plot([treat_doses[0], treat_doses[-1]], [1, 1], 'k--', alpha=0.5)
    ax.plot([treat_doses[0], treat_doses[-1]], [0.5, 0.5], 'r--', alpha=0.5)

    ax.set_xscale('log')
    ax.legend(loc=0)
    ax.set_xlabel('Selinexor dose (nM)')
    ax.set_ylabel('Survival')
    ax.set_title(name)
    fig.tight_layout()

    fig.savefig(os.path.join(res_dir, '%s_combi_dose_response.png' %name))

    # Heatmaps
    doses = product(treat_doses, com_doses)
    z = pd.DataFrame(index=com_doses, columns=treat_doses)

    for i, j in doses:
        loc = [(r[1]['Selinexor dose (nM)'] == i
                and r[1]['MK2206 dose (nM)'] == j)
               for r in norm_data.iterrows()]
        z.loc[j, i] = float(norm_data.loc[loc, 'Average'].values)

    x, y = np.meshgrid(treat_doses, com_doses)
    x = np.log10(x)
    y = np.log10(y)
    z = z.values.astype(float)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x, y, z, cmap='viridis')

    ax.view_init(30, 45)
    ax.set_xlabel('Selinexor dose log10(nM)')
    ax.set_ylabel('MK2206 dose log10(nM)')
    ax.set_zlabel('Survival')
    ax.set_title(name)
    fig.tight_layout()

    fig.savefig(os.path.join(res_dir, '%s_combi_3d.png' %name))
