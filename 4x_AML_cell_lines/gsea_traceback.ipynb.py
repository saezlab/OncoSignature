
# coding: utf-8

# In[93]:


import os

import pandas as pd
from IPython.display import display

from data_tools.databases import up_map


# In[94]:


gset_gsymbol = pd.read_csv('data/msigdb_path_to_gsymbol.csv')
gset_gsymbol.columns = ['gset', 'gname']
display(gset_gsymbol.head())
print len(gset_gsymbol)


# In[95]:


psites = pd.read_csv('data/norm_data.csv', index_col=0).index.values
print len(psites)


# In[96]:


mapdf = pd.DataFrame({'psite':psites,
                      'pprot':[i.split('_')[-2] for i in psites]})
display(mapdf.head())


# In[97]:


up_gname = up_map(set(mapdf['pprot']))
up_gname.columns = ['pprot', 'gname']
display(up_gname.head())
print len(up_gname)


# In[98]:


mapdf = pd.merge(mapdf, up_gname, how='left', on='pprot')
display(mapdf.head())
print len(mapdf)


# So, this is the mapping table from measured p-sites to their respective p-proteins and then to their associated gene-sets. From here, one could filter either the enriched gene-sets of the different contrasts or from the differentially expressed p-sites (also depending on the contrast).

# In[102]:


mapdf = pd.merge(mapdf, gset_gsymbol, how='left', on='gname')
display(mapdf.head())
print len(mapdf)


# In[104]:


mapdf.to_csv('results/psite_pprot_gset_table.csv', index=False)

