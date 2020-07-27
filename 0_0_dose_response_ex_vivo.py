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
# Dose-response model
# ===================
#
# This script fits a dose-response curve to the response data from our
# ex-vivo samples to classify them as responders or non-responders.

import os

import numpy as np
import pandas as pd
from scipy.stats import t
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import colors

from data_tools.models import DoseResponse
from data_tools.iterables import chunk_this

#----------------------------------- INPUT -----------------------------------#
data_dir = 'data'
res_dir = 'results'
resp_data = 'raw_survival_data.tab'

n = 4 # Number of columns on the multiplots
# Selected samples
use_samples = [6, 8, 9, 13, 14, 15, 18, 19, 24, 25,
               27, 28, 29, 34, 36, 39, 40, 41, 42, 44]
#-----------------------------------------------------------------------------#

def plot(ax, df):
    '''
    Generates a dose-response plot on the given the axes (ax) and data (df)
    '''

    ax.errorbar(df.index, df.mean(axis=1), yerr=df.std(axis=1), c='k')
    ax.scatter(np.repeat(df.index.values, df.shape[1]),
               df.values.reshape(-1),
               c=list(map('C{}'.format, list(range(df.shape[1])) * len(df))))


def hill(x, k, m, n):
    '''
    Hill function, where x is the dependant variable (dose) and k, m and n
    are the parameters.
    '''

    return m * x ** n / (k ** n + x ** n)


df_raw = pd.read_csv(os.path.join(data_dir, resp_data), sep='\t')

rep_cols = [i for i in df_raw.columns if i.startswith('Raw')]

norm_samples = {}

for s in df_raw.groupby('Sample #'):
    sample, data = s

    ref = data.loc[data['slx dose (nM)'] == 0, rep_cols]
    norm = data[rep_cols].values / ref.values

    df = pd.DataFrame(norm[2:],
                      index=data['slx dose (nM)'].values[2:],
                      columns=map('R{}'.format, range(1, 5)))

    norm_samples[sample] = df

for k, v in norm_samples.items():
    if (v.dropna().values > 2).any():
        print('Sample #%d' % k)
        print(v.dropna()[v.dropna().values > 1.5])

# We remove outlier measurements (e.g. after normalizing survival on 0nM
# there can't be 200% increase, either measurement error or contamination
# of the samples)
norm_samples[1][norm_samples[1] > 20] = np.nan
norm_samples[33][norm_samples[33] > 2] = np.nan

xmin = min(df_raw['slx dose (nM)'].unique()[2:])
xmax = max(df_raw['slx dose (nM)'].unique()[2:])

groups = chunk_this(use_samples, n)
m = len(groups)

fig, ax = plt.subplots(m, n, figsize=(3 * n, 3 * m),
                       sharex=True, sharey=True)

# j is vertical index, g is labels of group samples
for j, g in enumerate(groups):
    # i is horizontal index, s is sample id
    for i, s in enumerate(g):
        df = norm_samples[s]

        plot(ax[j, i], df)

        # Plot dashed lines at S=1 and S=0.5 for reference
        ax[j, i].plot([xmin, xmax], [1, 1],
                      '--k', alpha=.5)
        ax[j, i].plot([xmin, xmax], [.5, .5],
                      '--r', alpha=.5)

        ax[j, i].set_title('Sample #%d' %s)

        if j + 1 == m:
            ax[j, i].set_xlabel(r'$D$')

        if i == 0:
            ax[j, i].set_ylabel(r'$S$')

# Remove axis of empty plots in last line
idxs = [i for i in range(n) if len(ax[m - 1, i].get_lines()) == 0]

for i in idxs:
    ax[m - 1, i].axis('off')

# Dose-response figures
fig.legend([Patch(color='C%d' %i) for i in range(df.shape[1])],
           df.columns)
fig.tight_layout()
fig.savefig(os.path.join(res_dir, 'dose_response.pdf'))

# MODEL FITTING
#  We fit a Hill function to the response data (average across replicates)
# Model figures
fig, ax = plt.subplots(m, n, figsize=(3 * n, 3 * m),
                       sharex=True, sharey=True)
rng = np.linspace(xmin, xmax, 1000)

# Parameter order: k, m, n
x0 = [1e3, 1, -1] # Initial guess
bounds = ([0, 0, -np.inf], [np.inf, np.inf, 0]) # Parameter boundaries

fit_params = {}
ec50s = {}

# j is vertical index, g is labels of group samples
for j, g in enumerate(groups):
    # i is horizontal index, s is sample id
    for i, s in enumerate(g):
        xdata = norm_samples[s].index.values[:-1]
        ydata = norm_samples[s].mean(axis=1).values[:-1]

        model = DoseResponse(xdata, ydata, x0=x0, bounds=bounds)

        fit_params[s] = model.params
        ec50s[s] = model.ec()

        ax[j, i].scatter(xdata, ydata, label='Data')
        ax[j, i].plot(rng, hill(rng, *model.params), 'k', label='Fit')

        # Plot dashed lines at S=1 and S=0.5 for reference
        ax[j, i].plot([xmin, xmax], [1, 1],
                      '--k', alpha=.5)
        ax[j, i].plot([xmin, xmax], [.5, .5],
                      '--r', alpha=.5)

        ax[j, i].set_title('Sample #%d' %s)
        ax[j, i].text(0, 1.2, r'$k=%.1f$ | $m=%.2f$ | $n=%.2f$'
                      %tuple(model.params), fontsize=9)

        if j + 1 == m:
            ax[j, i].set_xlabel(r'$D$')

        if i == 0:
            ax[j, i].set_ylabel(r'$S$')

# Remove axis of empty plots in last line
idxs = [i for i in range(n) if len(ax[m - 1, i].get_lines()) == 0]

for i in idxs:
    ax[m - 1, i].axis('off')

# Barplot of EC50s
fig.tight_layout()
fig.savefig(os.path.join(res_dir, 'dose_response_fit.pdf'))
fig, ax = plt.subplots(figsize=(8, 5))

xmin = min(fit_params.keys())
xmax = max(fit_params.keys())

rng = range(len(fit_params))

for i, (k, v) in enumerate(fit_params.items()):

    val = ec50s[k]

    if val > 1e3:
        c = 'r'

    else:
        c = 'g'

    ax.bar(i, val if val != np.inf else 2e96, color=c)

ax.plot([rng[0] - 1, rng[-1] + 1], [1e3, 1e3], 'k--')
ax.set_xlim(rng[0] - 1, rng[-1] + 1)

ax.set_yscale('log')
ax.set_ylim(0, 1e5)
ax.set_xticks(rng)
ax.set_xticklabels(fit_params.keys(),
                   rotation=90, ha='center',
                   va='top')

ax.set_xlabel('Sample #')
ax.set_ylabel(r'$EC_{50}$')

ax.legend([Patch(color=i) for i in ['r', 'g']],
          [r'$EC_{50}>10^{3}$', r'$EC_{50}<10^{3}$'], loc=0)

fig.tight_layout()
fig.savefig(os.path.join(res_dir, 'ec50s.pdf'))

# Generating annotation file for further analyses
sample_class = pd.DataFrame(index=range(1, len(ec50s)),
                            columns=['sample', 'condition'])

for k, v in ec50s.items():
    sample_class.loc[k, 'sample'] = str(k)
    sample_class.loc[k, 'condition'] = 'R' if v <= 1e3 else 'NR'

missing = [1, 4, 7, 12, 26]
samples = pd.DataFrame(ec50s.values(), index=ec50s.keys(), columns=['EC50'])

for s in missing:

    try:
        samples.drop(s, inplace=True)

    except KeyError:
        pass

#samples['logEC50'] = np.log10(samples.EC50)
#def tval(x, ref=1e3):
#    den = np.std(x) / len(x)

#    return [(i - ref) / den for i in x]
#samples['tval'] = tval(samples.logEC50, ref=3)
#samples['pval'] = [t.sf(np.abs(v), len(samples) - 1) * 2 for v in samples.tval]
#samples[samples.pval > 0.01]
annot = []

for s in samples.iterrows():
    #if s[1].pval > 0.01:
    #    annot.append(np.nan)
    #else:
    annot.append('R' if s[1].EC50 <= 1e3 else 'NR')

samples['annotation'] = annot
samples['sample'] = samples.index
samples.to_csv(os.path.join(data_dir, 'sample_annotated.csv'), index=False)
