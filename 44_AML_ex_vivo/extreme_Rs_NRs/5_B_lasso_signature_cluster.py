#! /usr/bin/env python2

import matplotlib
matplotlib.use('Agg')

import time
from collections import Counter

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#from IPython.display import display
from sklearn.model_selection import StratifiedShuffleSplit as sss

from data_tools.models import Lasso
#from data_tools.databases import up_map
#from data_tools.plots import venn

#==============================================================================#
n = 100 # NUMBER OF MODELS
tt_ratio = 0.2 # TRAIN/TEST RATIO
njobs = 100 # NUMBER OF CORES
#==============================================================================#

print 'Loading data...\n'

data = pd.read_csv('data/norm_data.csv', index_col=0)
#display(data.head())

usecols = [c for c in data.columns if c.split('.')[0].endswith('B')]
#print usecols
data[pd.isna(data)] = 0
x = data[usecols].T

#display(x.head())

annot = pd.read_csv('data/sample_annotated.csv', index_col=5)
annot.index.name = None

#display(annot.head())
annot = annot.to_dict()['annotation']

y = pd.Series(index=x.index, dtype=int)

for c in x.index:
    aux = int(c.split('.')[0].strip('B'))

    y[c] = 0 if annot[aux] == 'NR' else 1

#display(y.head())

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

predictors.to_csv('results/predictors.csv')
intercepts.to_csv('results/intercepts.csv')
accs.to_csv('results/accs.csv')
pprobs.to_csv('results/pprobs.csv')
np.save('results/trainsets.npy', np.array(trainsets))
np.save('results/testsets.npy', np.array(testsets))