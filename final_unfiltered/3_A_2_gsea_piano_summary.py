import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

subdirs = [(x[0], [y for y in x[2] if (y.endswith('_rank.txt')
                                       and y.startswith('dist'))])
           for x in os.walk('results/3_gsea')][1:]
#subdirs

aux = dict()

for d, fs in subdirs:
    dname = d.split('/')[-1]
    #aux[dname] = dict()
    for f in fs:
        name = f.split('_')[1]
        df = pd.read_csv(os.path.join(d, f), sep='\t')
        aux['_'.join([dname, name])] = list(df[df.ConsRank == 1].index)

aux
from data_tools.iterables import subsets

print(aux.keys())

ssts = subsets([set(x) for x in aux.values()])

for k, v in ssts.items():
    if (sum(map(int, k)) == 1 or v == set()):
        continue

    else:
        print([list(aux.keys())[i] for i, n in enumerate(k) if int(n)])
        print(v, '\n')
