
# coding: utf-8

# # GSEA comparison: PIANO vs. MARINA

# Load packages

# In[1]:


import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from data_tools.plots import venn
from data_tools.iterables import unzip_dicts
from data_tools.iterables import subsets

get_ipython().magic(u'matplotlib inline')


# In order to compare the results from the different contrasts, we will use the following common data structure:
# - `<piano|viper>_scores` [`dict`]:
#     - Keys: contrast (e.g.: `'R2_R1'`) [`str`]
#     - Values: `dict`:
#         - Keys: direction (`'up'`/`'dw'`) [`str`]
#         - Values: gene-sets [`set`]

# Load results from `viper`

# In[2]:


subdir = 'results/3_gsea'

viper_df = pd.read_csv(os.path.join(subdir, 'results_gsea_viper.csv'),
                       index_col=0)
contrasts = viper_df.columns.tolist()

viper_scores = {}
for c in contrasts:
    viper_scores[c] = {}
    viper_scores[c]['MARINA up'] = set(viper_df[c].index[viper_df[c] > 0])
    viper_scores[c]['MARINA down'] = set(viper_df[c].index[viper_df[c] < 0])


# Load results from `piano`

# In[3]:


piano_scores = {c:{'PIANO up':set(), 'PIANO down':set()} for c in contrasts}
for c in contrasts:
    files = [f for f in os.listdir(os.path.join(subdir, c)) if f.endswith('_rank.txt')]
    for f in files:
        df = pd.read_table(os.path.join(subdir, c, f),
                           index_col=0)

        if 'up' in f:
            piano_scores[c]['PIANO up'].update(set(df.index))

        elif 'dw' in f:
            piano_scores[c]['PIANO down'].update(set(df.index))

        else:
            pass


# In order to visualize the comparison between algorithms, we will show the Venn diagram of the gene-sets as classified by `piano` and `viper` in each of the contrasts.
# 
# Also, for each contrast, the common gene-sets between algorithms will be shown. This is, the gene-sets contained in the following subsets:
# - MARINA up AND PIANO up
# - MARINA up AND PIANO down
# - MARINA down AND PIANO up
# - MARINA down AND PIANO down
# 
# Let us first define a couple of functions to simplify the code.

# In[4]:


def contains_both(groups):
    return True in ['MARINA' in i for i in groups] and True in ['PIANO' in i for i in groups]

def pprint(c, g, s):
    title = '{} - {} and {}'.format(c, *g)
    print title + '\n' + '-' * len(title) + '\n'
    print '\n'.join(s) + '\n'


# In[5]:


for c in contrasts:
    title = '{} vs. {}'.format(*c.split('_'))
    print title + '\n' + '=' * len(title) + '\n'

    keys, values = unzip_dicts(viper_scores[c], piano_scores[c])
    plot = venn(values, labels=keys, title=title)
    
    # Print the gene-sets common in both algorithms
    for k, v in subsets(values).items():
        groups = [keys[i] for (i, n) in enumerate(k) if n == '1']
        
        if (sum(map(int, k)) == 2 and contains_both(groups)):
            pprint(title, groups, sorted(list(v)))
        
        else:
            continue

