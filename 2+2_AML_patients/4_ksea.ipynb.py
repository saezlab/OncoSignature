
# coding: utf-8

# # Kinase-Substrate Enrichment Analysis
# 
# Considering two drug-treated groups of samples (responders and non-responders), here we will compute the KSEA for the following contrasts between groups:
# - R2 vs. R1: effect of treatment in responders
# - NR2 vs. NR1: effect of treatment in non-responders
# - R1 vs. NR1: difference between responders and non-responders prior to treatment.

# Loading packages

# In[2]:


import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

import kinact
from data_tools.databases import op_kinase_substrate, up_map
from data_tools.strings import join_str_lists
from data_tools.plots import volcano

get_ipython().magic(u'matplotlib inline')


# Load last version of kinase-substrate information from OmniPath.

# In[3]:


#ks_tab = op_kinase_substrate()
#ks_tab.head()


# In[4]:


# Join phosphosite identifiers
#ks_tab['psite'] = join_str_lists(join_str_lists(ks_tab['substrate'],
#                                        ks_tab['residue_type'], sep='_'),
#                                 ks_tab['residue_offset'])
# Map kinase gene names
#kinases = up_map(ks_tab['enzyme'])
#for acc, gname in kinases.values:
#    ks_tab.loc[ks_tab['enzyme'] == acc, 'kinase'] = gname


# Let us build/load the adjacency matrix between kinases and substrates

# In[5]:


#unique_kinases = ks_tab['kinase'].sort_values().unique()
#unique_psites = ks_tab['psite'].sort_values().unique()

#ks_adj = pd.DataFrame(index=unique_psites, columns=unique_kinases)

#for k in unique_kinases:
#    targets = ks_tab.loc[ks_tab['kinase'] == k, 'psite']
#    ks_adj.loc[targets, k] = 1

#ks_adj.to_csv('data/ks_adj.csv')

# Last saved: 20/08/2018 14:17
ks_adj = pd.read_csv('data/ks_adj.csv', index_col=0)


# Load the contrasts' data

# In[6]:


subdir = 'results/2_diff_exp/'

ttop_files = [f for f in os.listdir(subdir) if f.endswith('_ttop.csv')]
ttops = dict([f.replace('_ttop.csv', ''),
              pd.read_csv(subdir + f, index_col=0, usecols=[0, 1])] for f in ttop_files)


# ## KSEA results
# 
# Now, for each of the contrasts a KSEA is computed according to the $\log_2(FC)$ values of the phosphosites and their respective kinases. For each kinase, a score and corresponding p-value are obtained.

# In[7]:


results = {}

for k, v in ttops.items():
    score, pval = kinact.ksea.ksea_mean(v.dropna(), ks_adj, minimum_set_size=1)
    results[k] = pd.DataFrame({'score':score, 'pval':pval})


# In[8]:


pval_thr = 0.1

for k, v in results.items():
    title = 'KSEA {} vs. {}'.format(*k.split('_'))
    print 'Significant kinases of', title

    display(v[v['pval'] <= pval_thr])
    
    plot = volcano(v['score'], -np.log10(v['pval']), thr_fc=1, thr_pval=pval_thr,
                   title=title, filename='results/4_ksea/ksea_%s.pdf' %k)
    ax = plot.gca()
    ax.set_xlabel('KSEA score')
    
    v.to_csv('results/4_ksea/%s_KSEA_results.csv' %k)

