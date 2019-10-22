import os

import pandas as pd

from data_tools.databases import up_map
from data_tools.iterables import chunk_this

sdir = 'results/2_diff_exp'

files = [f for f in os.listdir(sdir) if f.endswith('_ttop.csv')]

for f in files:
    df = pd.read_csv(os.path.join(sdir, f), index_col=0)
    sig = [a and b for (a, b) in zip(abs(df['logFC']) >= 1,
                                     df['P.Value'] < 0.05)]
    df = df.loc[sig, :]

    ups = list(set([i.split('_')[0] for i in df.index]))

    mapper = dict()

    for ch in chunk_this(ups, 1000):
        aux = up_map(ch)
        mapper.update(aux.set_index('ACC').to_dict()['GENENAME'])

    df.index = ['_'.join([mapper[i.split('_')[0]], i.split('_')[1]])
                for i in df.index]

    df.to_csv(os.path.join(sdir, 'sig_' + f))
