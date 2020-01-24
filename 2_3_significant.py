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
# Significant differential phosphorylation
# ========================================
#
# This script subsets the differential expresssion results and saves
# a table for each contrast containing only the significantly
# differentially expressed p-sites.

import os

import pandas as pd

from data_tools.databases import up_map
from data_tools.iterables import chunk_this

#----------------------------------- INPUT -----------------------------------#
parent_dir = 'results'
#-----------------------------------------------------------------------------#

usedirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir)
           if d.startswith('2_diff_exp')]

for dir_ in usedirs:

    files = [f for f in os.listdir(dir_)
             if (f.endswith('_ttop.csv') and not f.startswith('sig_'))]

    for f in files:
        df = pd.read_csv(os.path.join(dir_, f), index_col=0)
        sig = [a and b for (a, b) in zip(abs(df['logFC']) >= 1,
                                         df['P.Value'] <= 0.05)]
        df = df.loc[sig, :]

        df.index

        ups = list(set([i.split('_')[0] for i in df.index]))

        mapper = dict()

        for ch in chunk_this(ups, 1000):
            aux = up_map(ch)
            mapper.update(aux.set_index('ACC').to_dict()['GENENAME'])

        df.index = ['_'.join([mapper[i.split('_')[0]],
                              '_'.join(i.split('_')[1:])])
                    for i in df.index]

        df.to_csv(os.path.join(dir_, 'sig_' + f))
