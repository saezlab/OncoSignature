# -*- coding: utf-8 -*-
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

from data_tools.models import DoseResponse

#----------------------------------- INPUT -----------------------------------#
data_dir = 'data'
res_dir = 'results'
#-----------------------------------------------------------------------------#

usefiles = [f for f in os.listdir(data_dir) if f.startswith('survival_combi')]
dataframes = dict()

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

    dataframes[name] = norm_data

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

    fig.savefig(os.path.join(res_dir, '%s_combi_dose_response.pdf' %name))

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

    fig.savefig(os.path.join(res_dir, '%s_combi_3d.pdf' % name))

# Checking 'synergy'
x0 = [1e3, 1, -1] # Initial guess
bounds = ([0, 0, -np.inf], [np.inf, np.inf, 0]) # Parameter boundaries
models = dict()
asmatrix = np.zeros((len(dataframes), len(com_doses)))

fig, ax = plt.subplots()

for i, (k, v) in enumerate(dataframes.items()): # For each cell line
    models[k] = dict()

    for n, df in v.groupby('MK2206 dose (nM)'): # For each MK2206 dose
        models[k][n] = DoseResponse(df['Selinexor dose (nM)'].values[:-1],
                                    df['Average'].values[:-1],
                                    x0=x0, bounds=bounds)

    ecs = [m.ec() for m in models[k].values()]
    doses = [d for d in models[k].keys()]

    asmatrix[i, :] = ecs

    rng = range(len(ecs))
    ax.plot(rng, ecs, label=k)

ax.set_xticks(rng)
ax.set_xticklabels(doses, rotation=90)
#ax.set_xlim(0, len(ecs))
ax.set_xlabel('MK2206 dose (nM)')
ax.set_ylabel('EC50 of Selinexor (nM)')
ax.legend()

fig.tight_layout()
fig.savefig(os.path.join(res_dir, 'combi_ec50s.pdf'))

fig, ax = plt.subplots()
im = ax.imshow(asmatrix, aspect='auto')
ax.set_yticks(range(asmatrix.shape[0]))
ax.set_yticklabels([k for k in dataframes.keys()])
ax.set_xticks(range(asmatrix.shape[1]))
ax.set_xticklabels(doses, rotation=90)
ax.set_xlabel('MK2206 dose (nM)')
ax.set_title('EC50 of Selinexor (nM)')
fig.colorbar(im)

fig.tight_layout()
fig.savefig(os.path.join(res_dir, 'ec50s_all.pdf'))

def isobole(x, k_x, k_y, m_x, m_y, p, q):
    return ((k_y**q/(2*m_y - 1)) ** (1/q) - (k_y**q*m_x*x**p/(k_x**p*m_y
            - m_x*x**p + m_y*x**p)) ** (1/q))

def resp_y(y, k_y, m_y, q):
    return m_y * y ** q / (k_y ** q + y ** q)

def y_of_x(x, k_x, k_y, m_x, m_y, p, q):
    return (k_y**q*m_x*x**p/(k_x**p*m_y - m_x*x**p + m_y*x**p)) ** (1/q)

# Generating the dose combination matrices
measurements = dict()
predicted = dict()

fig, ax = plt.subplots()
# Computing synergy
for k, v in dataframes.items():
    # Rows are Slx and columns are MK doses
    measurements[k] = pd.DataFrame(index=treat_doses, columns=com_doses)
    predicted[k] = pd.DataFrame(index=treat_doses, columns=com_doses)

    slx_dose, slx_resp = v.loc[v['MK2206 dose (nM)'] == 0,
                                 ['Selinexor dose (nM)', 'Average']].values.T
    mk_dose, mk_resp = v.loc[v['Selinexor dose (nM)'] == 0,
                                 ['MK2206 dose (nM)', 'Average']].values.T

    slx_mod = DoseResponse(slx_dose[:-1], slx_resp[:-1], x0=x0, bounds=bounds)
    k_x, m_x, p = slx_mod.params
    mk_mod = DoseResponse(mk_dose[:-1], mk_resp[:-1], x0=x0, bounds=bounds)
    k_y, m_y, q = mk_mod.params

    x = np.linspace(0, slx_mod.ec(), 1000)
    ax.plot(x, isobole(x, k_x, k_y, m_x, m_y, p, q), label=k)

    # Computing predicted response
    for i, r in v.iterrows():
        # Retrieve the doses of each drug
        x, y = r[['Selinexor dose (nM)', 'MK2206 dose (nM)']].values
        # Store the real response value for that combination
        measurements[k].loc[x, y] = r['Average']
        # Theoretical y
        pred_y = y_of_x(x, k_x, k_y, m_x, m_y, p, q) + y
        # Predicted response of theoretical y
        predicted[k].loc[x, y] = resp_y(pred_y, k_y, m_y, q)

ax.set_xlabel('Selinexor dose (nM)')
ax.set_ylabel('MK2206 dose (nM)')
ax.set_title(r'Isobole of EC$_{50}$')
ax.legend()
ax.set_xlim(0, 180)
ax.set_ylim(0, 4500)
fig.savefig(os.path.join(res_dir, 'isoboles.pdf'))

# Optmized Loewe (Loewe adapted to our DR model)
for k in measurements.keys():
    print(k)
    dif = predicted[k] - measurements[k]
    print(dif.mean().mean())

    fig, ax = plt.subplots()
    im = ax.imshow(dif.values.astype(float), cmap='coolwarm',
                   vmin=-0.5, vmax=0.5)
    ax.set_title('%s (%.3f)' % (k, dif.mean().mean()))
    fig.colorbar(im)
    ax.set_xticks(range(len(dif.columns)))
    ax.set_xticklabels(dif.columns, rotation=90)
    ax.set_yticks(range(len(dif.index)))
    ax.set_yticklabels(dif.index)

    ax.set_xlabel('MK2206 dose (nM)')
    ax.set_ylabel('Selinexor dose (nM)')
    fig.tight_layout()
    fig.savefig(os.path.join(res_dir, 'resp_dif_%s.pdf') % k)

# Bliss
for k, v in measurements.items():
    predicted = np.outer(v.iloc[:, 0], v.iloc[0, :])
    dif = predicted - v
    print(k)
    print(dif.mean().mean())
