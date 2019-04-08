# -*- coding: utf-8 -*-

import os

import pandas as pd
import numpy as np
import itertools as itt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from data_tools.plots import venn, cmap_bkgr

try:
    os.chdir('validation')

except:
    pass


# LOADING DAT-ASS
sdir = '../44_AML_ex_vivo/extreme_Rs_NRs/lasso_results'

predictors = pd.read_csv(os.path.join(sdir, 'predictors.csv'), index_col=0)
predictors[pd.isna(predictors)] = 0
intercepts = pd.read_csv(os.path.join(sdir, 'intercepts.csv'), index_col=0,
                         header=None)
accs = pd.read_csv(os.path.join(sdir, 'accs.csv'), index_col=0, header=None)

# Load the normalized validation dataset
norm = pd.read_csv('data/norm_data.csv', index_col=0) # Raw data all batches
normb = norm[[c for c in norm.columns if 'B' in c]] # Only untreated samples

# Reindexing data and predictors to same phosphosites and order
a, b = map(set, [predictors.index, normb.index])
print map(len, [a, b])

common = list(a.intersection(b))
len(common)

predictors = predictors.loc[common, :]
normb = normb.loc[common, :]

print normb.shape, predictors.shape
#data
normb[pd.isna(normb)] = 0


# Applying models
# ALL the models to the samples
res = pd.DataFrame(index=normb.columns, columns=range(predictors.shape[1]))

for n in range(predictors.shape[1]):
    m = predictors.iloc[:, n].values
    b = intercepts.iloc[n].values[0]
    res[n] = 1 / (1 + np.exp(-((normb.values.T * m).sum(axis=1) + b)))


# Models with ACC>80%
predictors80 = predictors.loc[:, (accs>=0.8).values.flatten()]
res80 = pd.DataFrame(index=normb.columns, columns=range(predictors80.shape[1]))

for n in range(predictors80.shape[1]):
    m = predictors80.iloc[:, n].values
    b = intercepts.iloc[n].values[0]
    res80[n] = 1 / (1 + np.exp(-((normb.values.T * m).sum(axis=1) + b)))

# Predictor's coefficients distribution
# Reloading predictors and data
# Load the normalized validation dataset
norm = pd.read_csv('data/norm_data.csv', index_col=0) # Raw data all batches
normb = norm[[c for c in norm.columns if 'B' in c]] # Only untreated samples

predictors = pd.read_csv(os.path.join(sdir, 'predictors.csv'), index_col=0)

# Subsetting predictors to the ones that had value > 0 in at least one model
predictors = predictors[pd.isna(predictors).sum(axis=1) < 100]

order = predictors.median(axis=1)
order = order.sort_values().index.values

# Let's remove small predictors
order = order[[abs(v) > 0.05 for v in predictors.loc[order, :].median(axis=1)]]

pred_dist = [predictors.loc[r, pd.notna(predictors.loc[r, :])]
             for r in order]
len(pred_dist)
# Validation data subsetted for the predictors:
aux = normb.loc[order, :]

# Plotty plotty!
fig, ax = plt.subplots(nrows=2, ncols=2, #sharey=True, sharex=True,
                       figsize=(20, 15),
                       gridspec_kw={'height_ratios': [1, 3],
                                    'width_ratios': [3, 1]})
fig.subplots_adjust(hspace=0.05)
fig.subplots_adjust(wspace=0.05)
im = ax[1, 0].imshow(aux.T, interpolation='none', aspect='auto', cmap=cmap_bkgr)

ax[0, 0].violinplot(pred_dist, widths=1, showmedians=True)
ax[0, 0].plot([0, len(aux.index) + 1], [0, 0], 'k--', alpha=0.5)

ax[1, 0].set_xticks(range(1, len(aux.index) + 1))
ax[1, 0].set_xticklabels(aux.index, rotation=90, fontsize=8)
ax[1, 0].set_yticks(range(len(aux.columns)))
ax[1, 0].set_yticklabels(aux.columns)

ax[0, 0].set_ylabel('Coefficient values')
ax[1, 0].set_ylabel('Samples')
ax[1, 0].set_xlabel('Predictor phosphosites (ordered by median coef. value)')

rng = range(normb.shape[1])
res = res.loc[aux.columns, :]
res80 = res80.loc[aux.columns, :]
# Same order as colormap
ax[1, 1].violinplot(res, vert=False, showmedians=True, widths=1,
                    positions=rng)
ax[1, 1].violinplot(res80, vert=False, showmedians=True, widths=1,
                    positions=rng)
ax[1, 1].plot([0.5, 0.5], [rng[0] - 1, rng[-1] + 1], 'k--', alpha=0.5)
ax[1, 1].set_xlabel('Prediction probability')
ax[1, 1].legend([Line2D([0], [0], color='C0'), Line2D([0], [0], color='C1')],
                ['ALL', r'$ACC\geq 80%$'], loc=0)

# Fixing the common axes
ax[1, 0].get_shared_x_axes().join(ax[0, 0], ax[1, 0])
ax[1, 0].get_shared_y_axes().join(ax[1, 1], ax[1, 0])

ax[0, 0].set_xlim(0, aux.shape[0]+1)
ax[0, 0].set_xticks([])
ax[1, 1].set_ylim(-1, aux.shape[1])
ax[1, 1].set_yticks([])

fig.tight_layout()
fig.delaxes(ax[0, 1])
fig.colorbar(im, ax=ax[0, 1])
fig.savefig('results/model_coef_validation.pdf')
