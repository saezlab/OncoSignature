
# coding: utf-8

# # Biomarkers of interest
# 
# Consider a phosphoproteomic dataset of treated and untreated samples from which some are responding to the treatement and others do not, we define the following groups:
# 
# - $R_1$: Responders before treatment
# - $R_2$: Responders after treatment
# - $NR_1$: Non-responders before treatment
# - $NR_2$: Non-responders after treatment
# 
# Let us define the following phoshposites of interest as:
# 
# 1. **Predictive biomarkers:**
# 
#     1.1. Phospho-sites are downregulated by drug treatment uniquely in R samples when compared to NR samples (where those sites are either not regulated or regulated up after drug treatment) and AT THE SAME time are at HIGHER levels pretreatment in R vs NR samples. They would be very useful to measure pretreatment and use as guidance for predicting efficacy based on pretreatment levels = efficacy-predictive biomarkers
# 
#     1.2. Phospho-sites are upregulated by drug treatment uniquely in R samples when compared to NR samples (where those sites are either not regulated or regulated down after drug treatment) and AT THE SAME time are at LOWER levels pretreatment in R vs NR samples. They would be very useful to measure pretreatment and use as guidance for predicting efficacy based on pretreatment levels = efficacy-predictive biomarkers
# 
# 2. **Excluder biomarkers:**
# 
#     2.1. Phospho-sites are downregulated by drug treatment uniquely in NR samples when compared to R samples (where those sites are either not regulated or regulated up after drug treatment) and AT THE SAME time are at HIGHER levels pretreatment in NR vs R samples. They would be very useful to measure pretreatment and use as guidance for predicting NO efficacy based on pretreatment levels = excluder biomarkers
# 
#     2.1. Phospho-sites are upregulated by drug treatment uniquely in NR samples when compared to R samples (where those sites are either not regulated or regulated down after drug treatment) and AT THE SAME time are at LOWER levels pretreatment in NR vs R samples. They would be very useful to measure pretreatment and use as guidance for predicting NO efficacy based on pretreatment levels = excluder biomarkers
# 
# From these descriptions, there are basically three contrasts to be compared:
# 
# - $R_2$ vs. $R_1$ $\rightarrow$ Effect of treatment in responders
# - $NR_2$ vs. $NR_1$ $\rightarrow$ Effect of treatment in non-responders
# - $R_1$ vs. $NR_1$ $\rightarrow$ Difference between responders and non-responders prior to the treatment
# 
# > **About contrasts and differential expression analysis**
# 
# > In terms of raw intensity, a contrast between samples of two groups, namely $A$ vs. $B$ is basically the ratio between their phosphosite intensities (aka fold-change or $FC$). Hence, for a given phosphosite $p$:
# 
# >$$FC(p) = \frac{p_A}{p_B}$$
# 
# >Such that if $FC(p)>1$ the phosphorilation level of $p_A$ is higher than $p_B$ (usually referred as up-regulation of $p$ in $A$ compared to $B$) and the inverse otherwise (also referred as down-regulation of $p$ in $A$ compared to $B$). It this sense, is very common to work with $\log_2(FC)$ since the result of this transformation yields positive values if $p_A>p_B$ and negative ones otherwise. Also, because $\log_2$ transformation of the raw intensities is a common practice when pre-processing and normalizing these kind of data. This further simplifies the computation of contrasts as simple subtractions, since:
# 
# >$$\log_2(FC(p))=\log_2\left(\frac{p_A}{p_B}\right) = \log_2(p_A) - \log_2(p_B)$$
# 
# Let us begin the analysis by loading some packages and the data from the differential expression analysis.

# ### Importing packages

# In[3]:


import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, Image

from data_tools.databases import up_map
from data_tools.databases import kegg_link
from data_tools.databases import kegg_pathway_mapping
from data_tools.plots import volcano


# In[4]:


get_ipython().magic(u'matplotlib inline')


# ### Loading data
# 
# As defined above, our contrasts of interest are:
# 
# - $R_2$ vs. $R_1$ $\rightarrow$ Effect of treatment in responders
# - $NR_2$ vs. $NR_1$ $\rightarrow$ Effect of treatment in non-responders
# - $R_1$ vs. $NR_1$ $\rightarrow$ Difference between responders and non-responders prior to treatment
# 
# Let us load the output of the differential expression analysis of these contrasts as computed by R's package `limma`. We will only use three of the columns, namely the phosphosite IDs, their $\log_2(FC)$ and the corresponding $p$-value.

# In[5]:


r2r1 = pd.read_csv('results/2_diff_exp/R.2_R.1_ttop.csv',
                   index_col=0, usecols=[0, 1, 4])
nr2nr1 = pd.read_csv('results/2_diff_exp/NR.2_NR.1_ttop.csv',
                     index_col=0, usecols=[0, 1, 4])
r1nr1 = pd.read_csv('results/2_diff_exp/R.1_NR.1_ttop.csv',
                    index_col=0, usecols=[0, 1, 4])

print 'Preview of the data:'
display('R2 vs R1', r2r1.head(),
        'NR2 vs NR1', nr2nr1.head(),
        'R1 vs NR1', r1nr1.head())


# ### Defining the subsets of phosphosites
# 
# We further subdivide the constrasted groups in up and down regulated (even though some cases do not have anything to do with actual biological regulation, this is just a term to express the "directionanilty" of the diferential expression [or phosphorilation in this case] of a contrast). Let us give them arbitrary names for the sake of simplicity as follows:
# 
# - $A\equiv$ Significantly up-regulated phosphosites in $R_2$ vs. $R_1$
# - $B\equiv$ Significantly down-regulated phosphosites in $R_2$ vs. $R_1$
# - $C\equiv$ Significantly up-regulated phosphosites in $NR_2$ vs. $NR_1$
# - $D\equiv$ Significantly down-regulated phosphosites in $NR_2$ vs. $NR_1$
# - $E\equiv$ Significantly up-regulated phosphosites in $R_1$ vs. $NR_1$
# - $F\equiv$ Significantly down-regulated phosphosites in $R_1$ vs. $NR_1$
# 
# Considering each sub-group as a set of phosphosites, we can easily define our biomarkers of interest in as:
# 
# $$
# 1.1\equiv (B\cap E)\backslash D
# $$
# 
# $$
# 1.2\equiv (A\cap F)\backslash C
# $$
# 
# $$
# 2.1\equiv (D\cap F)\backslash B
# $$
# 
# $$
# 2.2\equiv (C\cap E)\backslash A
# $$
# 
# Let us define a function to extract the subsets of a given contrast for any given thresholds of $p$-value and $FC$:

# In[6]:


def sig_subsets(df, pval=0.05, fc=2.0):
    
    p_fc_pairs = zip(df['P.Value'], df['logFC'])
    
    up = df.loc[[p <= pval and f >= np.log2(fc)
                 for (p, f) in p_fc_pairs], :]
    down = df.loc[[p <= pval and f <= -np.log2(fc)
                   for (p, f) in p_fc_pairs], :]
    
    return set(up.index.values), set(down.index.values)


# To avoid repetitive code, let us also define the following modificated volcano plotting function:

# In[7]:


def mod_volcano(df, notin_df, notin_set, notin_name, in_set, in_name,
                pval=0.05, fc=2.0, title=None, figsize=None):

    plot = volcano(df['logFC'], -np.log10(df['P.Value']), legend=False,
                   title=title, thr_pval=pval, thr_fc=fc, figsize=figsize)

    ax = plot.gca()

    notin_sites = [p for p in notin_df.index.values if p not in notin_set]
    
    # Sites of "not in set" but on df positions
    ax.scatter(df.loc[notin_sites, 'logFC'],
               -np.log10(df.loc[notin_sites, 'P.Value']),
               label=r'$\notin %s$' %notin_name, c='C2', s=25)

    # Sites of "in set" but on df positions
    ax.scatter(df.loc[in_set, 'logFC'],
               -np.log10(df.loc[in_set, 'P.Value']),
               label=r'$\in %s$' %in_name, c='C3', s=10)

    ax.legend(loc=0)


# Also a mapper function to link our phosphosite IDs to UniProt ACC, KEGG ID and KEGG pathway:

# In[8]:


def mapper(psites):
    site_up = pd.DataFrame({'PSite':psites,
                            'ACC':[i.split('_')[0] for i in psites]})
    up_keggid = up_map(site_up['ACC'], target='KEGG_ID')
    site_up_keggid = pd.merge(site_up, up_keggid, on='ACC', how='outer')
    
    keggid_keggpath = kegg_link(up_keggid['KEGG_ID'])
    
    if keggid_keggpath.iloc[0, 0] != '':
        keggid_keggpath.columns = ['KEGG_ID', 'KEGG_PATHWAY']

    else:
        keggid_keggpath = pd.DataFrame(columns=['KEGG_ID', 'KEGG_PATHWAY'])
            
    return pd.merge(site_up_keggid,
                    keggid_keggpath,
                    on='KEGG_ID',
                    how='outer')


# In[9]:


colors = {'1.1':'#ff0077',
          '1.2':'#00ffff',
          '2.1':'#ffff00',
          '2.2':'#ff7700'}


# ---

# ## "High" astringency
# 
# Default function values
# 
# $p$-value threshold = $0.05$
# 
# $FC$ threshold = $2^{1.0}$

# In[10]:


pval = 0.05
fc = 2 ** 1.0


# In[11]:


A, B = sig_subsets(r2r1, pval=pval, fc=fc)

C, D = sig_subsets(nr2nr1, pval=pval, fc=fc)

E, F = sig_subsets(r1nr1, pval=pval, fc=fc)


# The number of phosphosites on each of these subsets are as follows:

# In[12]:


print 'A = R2 vs R1 UP: %d\nB = R2 vs R1 DOWN: %d\n' %(len(A), len(B))
print 'C = NR2 vs NR1 UP: %d\nD = NR2 vs NR1 DOWN: %d\n' %(len(C), len(D))
print 'E = R1 vs NR1 UP: %d\nF = R1 vs NR1 DOWN: %d' %(len(E), len(F))


# In[11]:


poi = []


# ### Predictive phosphosites

# #### Case 1.1
# 
# Phospho-sites are downregulated by drug treatment uniquely in R samples when compared to NR samples (where those sites are either not regulated or regulated up after drug treatment) and AT THE SAME time are at HIGHER levels pretreatment in R vs NR samples. They would be very useful to measure pretreatment and use as guidance for predicting efficacy based on pretreatment levels = efficacy-predictive biomarkers
# 
# 
# $$(B\cap E)\backslash D$$

# In[12]:


mod_volcano(r2r1, nr2nr1, D, 'D', E, 'E', title=r'$R_2$ vs. $R_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(B.intersection(E).difference(D))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['1.1'] * len(aux)
    poi.append(aux)


# #### Case 1.2
# 
# Phospho-sites are upregulated by drug treatment uniquely in R samples when compared to NR samples (where those sites are either not regulated or regulated down after drug treatment) and AT THE SAME time are at LOWER levels pretreatment in R vs NR samples. They would be very useful to measure pretreatment and use as guidance for predicting efficacy based on pretreatment levels = efficacy-predictive biomarkers
# 
# $$(A\cap F)\backslash C$$

# In[13]:


mod_volcano(r2r1, nr2nr1, C, 'C', F, 'F', title=r'$R_2$ vs. $R_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(A.intersection(F).difference(C))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['1.2'] * len(aux)
    poi.append(aux)


# ### Excluder phosphosites

# #### Case 2.1
# 
# Phospho-sites are downregulated by drug treatment uniquely in NR samples when compared to R samples (where those sites are either not regulated or regulated up after drug treatment) and AT THE SAME time are at HIGHER levels pretreatment in NR vs R samples. They would be very useful to measure pretreatment and use as guidance for predicting NO efficacy based on pretreatment levels = excluder biomarkers
# 
# $$(D\cap F)\backslash B$$

# In[14]:


mod_volcano(nr2nr1, r2r1, B, 'B', F, 'F', title=r'$NR_2$ vs. $NR_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(D.intersection(F).difference(B))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['2.1'] * len(aux)
    poi.append(aux)


# #### Case 2.2
# 
# Phospho-sites are upregulated by drug treatment uniquely in NR samples when compared to R samples (where those sites are either not regulated or regulated down after drug treatment) and AT THE SAME time are at LOWER levels pretreatment in NR vs R samples. They would be very useful to measure pretreatment and use as guidance for predicting NO efficacy based on pretreatment levels = excluder biomarkers
# 
# $$(C\cap E)\backslash A$$

# In[15]:


mod_volcano(nr2nr1, r2r1, A, 'A', E, 'E', title=r'$NR_2$ vs. $NR_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(C.intersection(E).difference(A))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['2.2'] * len(aux)
    poi.append(aux)


# From the different interesting phosphosites of each case, we mapped their IDs to the related KEGG  pathways (NOTE: not all interesting phosphosites appear in a KEGG pathway, hence we're "losing" or missing some of them).
# 
# Then, we count the number of these that appear on each pathway so that we can show the top pathways that contain the maximum number of our proteins of interest.

# In[16]:


pois = pd.concat(poi)
counts = pd.Series()
for k in pois.groupby('KEGG_PATHWAY'):
    if k[0] != 'path:hsa01100':
        counts[k[0]] = len(k[1])

top5 = counts.sort_values(ascending=False).head(5)
top5


# Color legend:
# 
# - <span style="color:#ff0077">**Case 1.1**</span>
# - <span style="color:#00ffff">**Case 1.2**</span>
# - <span style="color:#ffff00">**Case 2.1**</span>
# - <span style="color:#ff7700">**Case 2.2**</span>

# In[17]:


for pathw in top5.index:
    sub_pois = pois.loc[pois['KEGG_PATHWAY'] == pathw, :]
    
    query = pd.DataFrame({'A':[i.split(':')[1] for i in sub_pois['KEGG_ID'].values],
                          'B':[colors[i] for i in sub_pois['case'].values]})

    print pathw
    display(sub_pois)
    
    pathid = pathw.split(':')[1]
    
    kegg_pathway_mapping(query, pathid)

    display(Image(pathid + '.png'))
    
    os.remove(pathid + '.png')


# ---

# ## "Medium" astringency
# 
# Default function values
# 
# $p$-value threshold = $0.1$
# 
# $FC$ threshold = $2^{0.9}$

# In[18]:


pval = 0.1
fc = 2 ** 0.9


# In[19]:


A, B = sig_subsets(r2r1, pval=pval, fc=fc)

C, D = sig_subsets(nr2nr1, pval=pval, fc=fc)

E, F = sig_subsets(r1nr1, pval=pval, fc=fc)


# The number of phosphosites on each of these subsets are as follows:

# In[20]:


print 'A = R2 vs R1 UP: %d\nB = R2 vs R1 DOWN: %d\n' %(len(A), len(B))
print 'C = NR2 vs NR1 UP: %d\nD = NR2 vs NR1 DOWN: %d\n' %(len(C), len(D))
print 'E = R1 vs NR1 UP: %d\nF = R1 vs NR1 DOWN: %d' %(len(E), len(F))


# In[21]:


poi = []


# ### Predictive phosphosites

# #### Case 1.1
# 
# Phospho-sites are downregulated by drug treatment uniquely in R samples when compared to NR samples (where those sites are either not regulated or regulated up after drug treatment) and AT THE SAME time are at HIGHER levels pretreatment in R vs NR samples. They would be very useful to measure pretreatment and use as guidance for predicting efficacy based on pretreatment levels = efficacy-predictive biomarkers
# 
# 
# $$(B\cap E)\backslash D$$

# In[22]:


mod_volcano(r2r1, nr2nr1, D, 'D', E, 'E', title=r'$R_2$ vs. $R_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(B.intersection(E).difference(D))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['1.1'] * len(aux)
    poi.append(aux)


# #### Case 1.2
# 
# Phospho-sites are upregulated by drug treatment uniquely in R samples when compared to NR samples (where those sites are either not regulated or regulated down after drug treatment) and AT THE SAME time are at LOWER levels pretreatment in R vs NR samples. They would be very useful to measure pretreatment and use as guidance for predicting efficacy based on pretreatment levels = efficacy-predictive biomarkers
# 
# $$(A\cap F)\backslash C$$

# In[23]:


mod_volcano(r2r1, nr2nr1, C, 'C', F, 'F', title=r'$R_2$ vs. $R_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(A.intersection(F).difference(C))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['1.2'] * len(aux)
    poi.append(aux)


# ### Excluder phosphosites

# #### Case 2.1
# 
# Phospho-sites are downregulated by drug treatment uniquely in NR samples when compared to R samples (where those sites are either not regulated or regulated up after drug treatment) and AT THE SAME time are at HIGHER levels pretreatment in NR vs R samples. They would be very useful to measure pretreatment and use as guidance for predicting NO efficacy based on pretreatment levels = excluder biomarkers
# 
# $$(D\cap F)\backslash B$$

# In[24]:


mod_volcano(nr2nr1, r2r1, B, 'B', F, 'F', title=r'$NR_2$ vs. $NR_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(D.intersection(F).difference(B))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['2.1'] * len(aux)
    poi.append(aux)


# #### Case 2.2
# 
# Phospho-sites are upregulated by drug treatment uniquely in NR samples when compared to R samples (where those sites are either not regulated or regulated down after drug treatment) and AT THE SAME time are at LOWER levels pretreatment in NR vs R samples. They would be very useful to measure pretreatment and use as guidance for predicting NO efficacy based on pretreatment levels = excluder biomarkers
# 
# $$(C\cap E)\backslash A$$

# In[25]:


mod_volcano(nr2nr1, r2r1, A, 'A', E, 'E', title=r'$NR_2$ vs. $NR_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(C.intersection(E).difference(A))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['2.2'] * len(aux)
    poi.append(aux)


# From the different interesting phosphosites of each case, we mapped their IDs to the related KEGG  pathways (NOTE: not all interesting phosphosites appear in a KEGG pathway, hence we're "losing" or missing some of them).
# 
# Then, we count the number of these that appear on each pathway so that we can show the top pathways that contain the maximum number of our proteins of interest.

# In[26]:


pois = pd.concat(poi)
counts = pd.Series()
for k in pois.groupby('KEGG_PATHWAY'):
    if k[0] != 'path:hsa01100':
        counts[k[0]] = len(k[1])

top5 = counts.sort_values(ascending=False).head(5)
top5


# Color legend:
# 
# - <span style="color:#ff0077">**Case 1.1**</span>
# - <span style="color:#00ffff">**Case 1.2**</span>
# - <span style="color:#ffff00">**Case 2.1**</span>
# - <span style="color:#ff7700">**Case 2.2**</span>

# In[27]:


for pathw in top5.index:
    sub_pois = pois.loc[pois['KEGG_PATHWAY'] == pathw, :]
    
    query = pd.DataFrame({'A':[i.split(':')[1] for i in sub_pois['KEGG_ID'].values],
                          'B':[colors[i] for i in sub_pois['case'].values]})

    print pathw
    display(sub_pois)
    
    pathid = pathw.split(':')[1]
    
    kegg_pathway_mapping(query, pathid)

    display(Image(pathid + '.png'))
    
    os.remove(pathid + '.png')


# ---

# ## "Low" astringency
# 
# Default function values
# 
# $p$-value threshold = $0.1$
# 
# $FC$ threshold = $2^{0.8}$

# In[28]:


pval = 0.1
fc = 2 ** 0.8


# In[29]:


A, B = sig_subsets(r2r1, pval=pval, fc=fc)

C, D = sig_subsets(nr2nr1, pval=pval, fc=fc)

E, F = sig_subsets(r1nr1, pval=pval, fc=fc)


# The number of phosphosites on each of these subsets are as follows:

# In[30]:


print 'A = R2 vs R1 UP: %d\nB = R2 vs R1 DOWN: %d\n' %(len(A), len(B))
print 'C = NR2 vs NR1 UP: %d\nD = NR2 vs NR1 DOWN: %d\n' %(len(C), len(D))
print 'E = R1 vs NR1 UP: %d\nF = R1 vs NR1 DOWN: %d' %(len(E), len(F))


# In[31]:


poi = []


# ### Predictive phosphosites

# #### Case 1.1
# 
# Phospho-sites are downregulated by drug treatment uniquely in R samples when compared to NR samples (where those sites are either not regulated or regulated up after drug treatment) and AT THE SAME time are at HIGHER levels pretreatment in R vs NR samples. They would be very useful to measure pretreatment and use as guidance for predicting efficacy based on pretreatment levels = efficacy-predictive biomarkers
# 
# 
# $$(B\cap E)\backslash D$$

# In[32]:


mod_volcano(r2r1, nr2nr1, D, 'D', E, 'E', title=r'$R_2$ vs. $R_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(B.intersection(E).difference(D))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['1.1'] * len(aux)
    poi.append(aux)


# #### Case 1.2
# 
# Phospho-sites are upregulated by drug treatment uniquely in R samples when compared to NR samples (where those sites are either not regulated or regulated down after drug treatment) and AT THE SAME time are at LOWER levels pretreatment in R vs NR samples. They would be very useful to measure pretreatment and use as guidance for predicting efficacy based on pretreatment levels = efficacy-predictive biomarkers
# 
# $$(A\cap F)\backslash C$$

# In[33]:


mod_volcano(r2r1, nr2nr1, C, 'C', F, 'F', title=r'$R_2$ vs. $R_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(A.intersection(F).difference(C))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['1.2'] * len(aux)
    poi.append(aux)


# ### Excluder phosphosites

# #### Case 2.1
# 
# Phospho-sites are downregulated by drug treatment uniquely in NR samples when compared to R samples (where those sites are either not regulated or regulated up after drug treatment) and AT THE SAME time are at HIGHER levels pretreatment in NR vs R samples. They would be very useful to measure pretreatment and use as guidance for predicting NO efficacy based on pretreatment levels = excluder biomarkers
# 
# $$(D\cap F)\backslash B$$

# In[34]:


mod_volcano(nr2nr1, r2r1, B, 'B', F, 'F', title=r'$NR_2$ vs. $NR_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(D.intersection(F).difference(B))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['2.1'] * len(aux)
    poi.append(aux)


# #### Case 2.2
# 
# Phospho-sites are upregulated by drug treatment uniquely in NR samples when compared to R samples (where those sites are either not regulated or regulated down after drug treatment) and AT THE SAME time are at LOWER levels pretreatment in NR vs R samples. They would be very useful to measure pretreatment and use as guidance for predicting NO efficacy based on pretreatment levels = excluder biomarkers
# 
# $$(C\cap E)\backslash A$$

# In[35]:


mod_volcano(nr2nr1, r2r1, A, 'A', E, 'E', title=r'$NR_2$ vs. $NR_1$',
            pval=pval, fc=fc, figsize=(9, 7))

selected = list(C.intersection(E).difference(A))
print 'Number of unique phosphosites selected:', len(selected)

if len(selected) > 0:
    aux = mapper(selected)
    display(aux)

    aux.dropna(inplace=True)
    aux['case'] = ['2.2'] * len(aux)
    poi.append(aux)


# From the different interesting phosphosites of each case, we mapped their IDs to the related KEGG  pathways (NOTE: not all interesting phosphosites appear in a KEGG pathway, hence we're "losing" or missing some of them).
# 
# Then, we count the number of these that appear on each pathway so that we can show the top pathways that contain the maximum number of our proteins of interest.

# In[36]:


pois = pd.concat(poi)
counts = pd.Series()
for k in pois.groupby('KEGG_PATHWAY'):
    if k[0] != 'path:hsa01100':
        counts[k[0]] = len(k[1])

top5 = counts.sort_values(ascending=False).head(5)
top5


# Color legend:
# 
# - <span style="color:#ff0077">**Case 1.1**</span>
# - <span style="color:#00ffff">**Case 1.2**</span>
# - <span style="color:#ffff00">**Case 2.1**</span>
# - <span style="color:#ff7700">**Case 2.2**</span>

# In[37]:


for pathw in top5.index:
    sub_pois = pois.loc[pois['KEGG_PATHWAY'] == pathw, :]
    
    query = pd.DataFrame({'A':[i.split(':')[1] for i in sub_pois['KEGG_ID'].values],
                          'B':[colors[i] for i in sub_pois['case'].values]})

    print pathw
    display(sub_pois)
    
    pathid = pathw.split(':')[1]
    
    kegg_pathway_mapping(query, pathid)

    display(Image(pathid + '.png'))
    
    os.remove(pathid + '.png')


# ---
