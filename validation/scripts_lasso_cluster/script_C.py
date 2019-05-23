'''
LASSO Signature - Approach C
============================

Here we just use all the training samples (from the significant Rs/NRs
set of samples) and all p-sites that appear in the validation set (V1
and V2).
'''


from __future__ import division

import os
import time
from collections import Counter

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit as sss

from data_tools.models import Lasso

#==============================================================================#
n = 1000 # NUMBER OF MODELS
tt_ratio = 0.2 # TRAIN/TEST RATIO
njobs = 100 # NUMBER OF CORES
res_dir = 'results/C'
#==============================================================================#

if not os.path.exists(res_dir):
    os.makedirs(res_dir)

print 'Loading data...\n'

val = ['47B', '49B', '50B', '53B', '55B', '57B', '63B', '65B']

data = pd.read_csv('data/norm_data_all.csv', index_col=0)

usecols = [c for c in data.columns if c.split('.')[0].endswith('B')]

data[pd.isna(data)] = 0
x = data[usecols].T

for i, s in enumerate(val):
    if i == 0:
        usesites = set(x.loc[s, x.loc[s, :] != 0].index)
    else:
        usesites.intersection(set(x.loc[s, x.loc[s, :] != 0].index))

annot = pd.read_csv('data/sample_annotated.csv', index_col=5)
annot.index.name = None

y = dict([(str(i[0]) + 'B', 1 if i[1][4] == 'R' else 0)
          for i in annot.iterrows() if pd.notna(i[1][4])])
y = pd.Series(y)

x = x.loc[y.index, usesites]
x.shape
sampler = sss(n_splits=n, test_size=tt_ratio)

predictors = pd.DataFrame(index=x.columns, columns=range(n))
intercepts = pd.Series(index=range(n))
accs = pd.Series(index=range(n))
pprobs = pd.DataFrame(index=x.index, columns=range(n))
trainsets = []
testsets = []

gstart = time.time()

for i, (train, test) in enumerate(sampler.split(x, y)):
    start = time.time()
    msg = 'Computing model #%d' % (i + 1)
    print msg
    print '=' * len(msg) + '\n'

    trainsets.append(list(train))
    testsets.append(list(test))

    model = Lasso(Cs=np.linspace(1e-2, 1e0, 1000), sampler='skf', n_jobs=njobs)
    model.fit_data(x.iloc[train, :], y.iloc[train], silent=True)

    tx = x.iloc[test, :]

    predictors.loc[model.predictors.index, i] = model.predictors.values
    accs[i] = model.score(tx, y.iloc[test])
    pprobs.loc[tx.index, i] = model.predict_proba(tx)[:, 1]
    intercepts[i] = model.intercept_[0]

    print '  Model trained, elapsed time %.4f sec.' % (time.time() - start)
    print '  Number of predictors: %d' % sum(pd.notnull(predictors[i]))
    print '  Accuracy: %.4f\n' % accs[i]

print '\nComputed %d models, elapsed time %.4f sec.' %(n, time.time() - gstart)

predictors.to_csv(os.path.join(res_dir, 'predictors.csv'))
intercepts.to_csv(os.path.join(res_dir, 'intercepts.csv'))
accs.to_csv(os.path.join(res_dir, 'accs.csv'))
pprobs.to_csv(os.path.join(res_dir, 'pprobs.csv'))
np.save(os.path.join(res_dir, 'trainsets.npy'), np.array(trainsets))
np.save(os.path.join(res_dir, 'testsets.npy'), np.array(testsets))
