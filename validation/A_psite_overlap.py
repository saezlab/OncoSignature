# -*- coding: utf-8 -*-

import os

import pandas as pd
import numpy as np
import itertools as itt
import matplotlib.pyplot as plt

from data_tools.iterables import similarity
from data_tools.plots import venn

try:
    os.chdir('validation')

except:
    pass

#=============================== LOADING DATA ================================#
raw = pd.read_csv('data/raw_all.csv', index_col=0)
rawb = raw[[c for c in raw.columns if 'B' in c]]
#raw = raw[sorted(raw.columns)]

# p-sites with NaN in all samples
(pd.isna(raw).sum(axis=1) == len(raw)).sum()

#============================ PROPORTION OF NaNs =============================#
# Check number of NaNs per sample
size = float(len(rawb))
pct = (pd.isna(rawb).sum() / size).sort_values()

# Make a plot
rng = range(len(pct))
fig, ax = plt.subplots(figsize=(15, 5))

ax.bar(rng, pct.values)
size
ax.set_xticks(rng)
ax.set_xticklabels(pct.index, rotation=90)
ax.set_xlim(rng[0] - 1, rng[-1] + 1)

ax.set_xlabel('Samples')
ax.set_ylabel('Missing values (%)')
fig.tight_layout()
fig.savefig('results/nan_vals.pdf')

#=================== SAMPLE SIMILARITY - MEASURED P-SITES ====================#
psites = dict([(c, set(rawb.index[pd.notna(rawb[c])])) for c in rawb.columns])

sims = pd.DataFrame(index=rawb.columns,
                    columns=rawb.columns)

pos = [2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 46.5, 50.5]
batches = ['Set1', 'Set2', 'Set3', 'Set4', 'Set5', 'Set6', 'Set7', 'Set8',
           'V1', 'V2', 'V3']

# By similarity index
for (name, mode) in [('Jaccard', 'j'), ('Sorensen-Dice', 'sd'), ('Szymkiewicz-Simpson', 'ss')]:
    for (a, b) in itt.product(rawb.columns, repeat=2):
        sims.loc[a, b] = similarity(psites[a], psites[b], mode=mode)

    fig, ax = plt.subplots(figsize=(15, 15))
    im = ax.imshow(sims.values.astype(float), cmap='nipy_spectral', interpolation='none')
    fig.colorbar(im, shrink=0.75)

    ax.set_xticks(pos)
    ax.set_xticklabels(batches, rotation=90)
    ax.set_yticks(pos)
    ax.set_yticklabels(batches)
    ax.set_title(u'%s similarity index of measured p-sites across samples'
                 % name)
    fig.tight_layout()

    fig.savefig('results/%s_sim_psites.pdf' % mode)


    aux = [similarity(psites[a], psites[b], mode=mode) for (a, b)
           in itt.combinations_with_replacement(rawb.columns, 2)]

    fig, ax = plt.subplots()
    ax.hist(aux, bins=50)
    ax.set_xlabel(u'%s similarity index of measured p-sites across samples'
                  % name)
    ax.set_ylabel('Frequency (counts)')
    fig.tight_layout()
    fig.savefig('results/%s_sim_psites_hist.pdf' % mode)

#==================== P-SITE OVERLAP VALIDATION vs. TRAIN ====================#
# Measured p-sites overlap (all training vs. vlaidation batches)
v3 = rawb.columns[-5:].tolist()
v2 = rawb.columns[-8:-5].tolist()
v1 = rawb.columns[-13:-8].tolist()

sets = {'V1': set(), 'V2': set(), 'V3': set(), 'Train': set()}

for (k, v) in psites.items():

    if k in v1:
        sets['V1'].update(v)

    elif k in v2:
        sets['V2'].update(v)

    elif k in v3:
        sets['V3'].update(v)

    else:
        sets['Train'].update(v)

plot = venn(sets.values(), labels=sets.keys(), sizes=True,
            filename='results/psite_overlap_train.pdf')

#============== P-SITE OVERLAP MODELS vs. VALIDATION vs. TRAIN ===============#
sdir = '../44_AML_ex_vivo/extreme_Rs_NRs/lasso_results'
predictors = pd.read_csv(os.path.join(sdir, 'predictors.csv'), index_col=0)
#predictors[pd.isna(predictors)] = 0
predictors.shape

# Subsetting the predictors to the ones that had value > 0 in at least one model
predictors = predictors[pd.isna(predictors).sum(axis=1) < 100]
predictors.shape

sets['Models'] = set(predictors.index)

plot = venn(sets.values(), labels=sets.keys(), sizes=True, filename='results/psite_overlap_train_model.pdf')

# Double-check the p-site overlap of old dataset (train) vs new dataset (train
# + validation)
tdata = pd.read_csv('../44_AML_ex_vivo/data/raw.csv', index_col=0)
#tdatab =  tdata[[c for c in tdata.columns if 'B' in c]]
#len(tdata)
#(pd.isna(tdata).sum(axis=1) == len(tdata)).sum()
#m = []
#for k in tdata.index:
#    if k not in raw.index:
#        m.append(k)
#len(m)
plot = venn([set(raw.index), set(tdata.index)], labels=['NEW', 'OLD'], sizes=True,
            filename='results/oldVSnew.pdf')

#=================== P-SITE OVERLAP MODELS vs. VALIDATION ====================#

del sets['Train']
v = [sets.values()[3]]+sets.values()[:3]
k = [sets.keys()[3]]+sets.keys()[:3]
plot = venn(v, labels=k, sizes=True,
            filename='results/psite_overlap_model.pdf')

#===================== CHECKING OLD vs NEW CORRELATION =======================#
common = list(set(raw.columns).intersection(set(tdata.columns)))

cors = dict()

for c in common:
    aux = raw[c].to_frame().merge(tdata[c].to_frame(), how='outer',
                                  right_index=True, left_index=True,
                                  suffixes=['_NEW', '_OLD'])
    cors[c] = aux.corr().iloc[0, 1]
    aux[pd.isna(aux)] = 0

fig, ax = plt.subplots(figsize=(12, 5))
rng = range(len(cors))
ax.bar(rng, cors.values())
ax.plot([rng[0] - 1, rng[-1] + 1], [1, 1], '--k')
ax.set_xticks(rng)
ax.set_xticklabels(cors.keys(), rotation=90)
ax.set_ylim(.99, 1.005)
ax.set_xlim(rng[0] - 1, rng[-1] + 1)
ax.set_ylabel('Correlation coefficient')
ax.set_title('Correlation among samples from different data files')
fig.tight_layout()
fig
fig.savefig('results/corr_samples.pdf')

old = tdata.copy()
old[pd.isna(old)] = 0
new = raw.copy()
new[pd.isna(new)] = 0

difs = pd.DataFrame([new[c] - old[c] for c in common], index=common).T
difs[pd.isna(difs)] = 0

difs.shape

difs = difs.loc[(difs != 0).sum(axis=1) > 1, :]
difs.shape
len(difs.columns)
difs = difs[tdata.columns]
len(difs.columns)

fig, ax = plt.subplots(figsize=(25, 12))

im = ax.imshow(np.log(difs.values).T, aspect='auto',
               cmap='cool',
               interpolation='none')
ax.set_title(r'$\log_{10}$ differences between data files')
ax.set_xlabel('Non-equal p-sites')
ax.set_ylabel('Samples (batch-ordered)')
ax.set_yticks(range(difs.shape[1]))
ax.set_yticklabels(difs.columns)
#ax.set_xticks(range(difs.shape[0]))
#ax.set_xticklabels(difs.index, size=6, rotation=90)
fig.colorbar(im)
fig.tight_layout()
fig.savefig('results/difs.pdf')
