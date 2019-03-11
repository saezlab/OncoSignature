
# coding: utf-8

# # Regression model for response prediction
# 
# Considering all samples prior to drug treatment (suffixed with **B**) and their annotation as responders (R) and non-responders (NR), we will build a logistic regression model with least absolute shrinkage and selection operator (LASSO) regularization. This will allow us to find which phosphosites are the main predictors of the patient response to the treatment.
# 
# Let us first load some modules.

# In[1]:


import time
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
from sklearn.model_selection import StratifiedShuffleSplit as sss

from data_tools.models import Lasso
from data_tools.databases import up_map
from data_tools.plots import venn


# In[2]:


get_ipython().magic(u'matplotlib inline')


# ## Loading the data
# 
# We will use the previously normalized and batch-corrected data, from which the treated samples will be excluded (those ending with **A**). For convenience we will substitute all NaN's by zeros.

# In[3]:


data = pd.read_csv('data/norm_data.csv', index_col=0)
#display(data.head())

usecols = [c for c in data.columns if c.split('.')[0].endswith('B')]
#print usecols
data[pd.isna(data)] = 0
x = data[usecols].T

#display(x.head())


# Similarly, we will load the annotation as R/NR obtained from the dose-response model and convert it to a binary variable where 0 denotes NR's and 1 correspods to R's.

# In[4]:


annot = pd.read_csv('data/sample_annotated.csv', index_col=0)
annot.index.name = None

#display(annot.head())
annot = annot.to_dict()['condition']

y = pd.Series(index=x.index, dtype=int)
for c in x.index:
    aux = int(c.split('.')[0].strip('B'))
    y[c] = 0 if annot[aux] == 'NR' else 1

#display(y.head())


# ## Training the models
# 
# First of all, let us separate our samples into train and test sets 100 times. This is done with the stratified shuffle split strategy. This is, the samples are separated randomly but keeping the proportion of classes (R's and NR's) among both subsets. The samples are separated in a 4:1 train/test ratio.
# 
# For each train subset, a logistic regression model is computed with LASSO regularization and 10-fold cross-validation (CV) in order to fit the regularization parameter of that model.
# 
# The respective test set is used to predict the classes of those samples and compare to their real annotation, thus obtaining a proper accuracy score of that individual model.

# In[5]:


n = 100
sampler = sss(n_splits=n, test_size=0.2)

predictors = []
trainsets = []
testsets = []
accs = []

gstart = time.time()

for i, (train, test) in enumerate(sampler.split(x, y)):
    start = time.time()
    msg = 'Computing model #%d' % (i + 1)
    print msg
    print '=' * len(msg) + '\n'
    
    trainsets.append(train)
    testsets.append(test)
    
    model = Lasso(Cs=np.linspace(1e-2, 1e0, 1000), sampler='skf', n_jobs=4)
    model.fit_data(x.iloc[train, :], y.iloc[train], silent=True)
    predictors.append(model.predictors)
    accs.append(model.score(x.iloc[test, :], y.iloc[test]))
    
    print '  Model trained, elapsed time %.4f sec.' % (time.time() - start)
    print '  Number of predictors: %d' % len(predictors[-1])
    print '  Accuracy: %.4f\n' % accs[-1]
    
print '\nComputed %d models, elapsed time %.4f sec.' %(n, time.time() - gstart)


# Let's see how the models performed in terms of accuracy:

# In[6]:


accs = pd.Series(accs)

plt.figure()
plt.bar(range(n), accs)
plt.xlabel('Model #')
plt.ylabel('ACC')

plt.figure()
accs.hist()
plt.title('Average ACC = %.3f' % accs.mean())


# In order to find out which phosphosites are better predictors for R/NR classification, let's plot which are the top phosphosites that appeared as non-zero parameters across the different models:

# In[7]:


ps = []

for p in predictors:
    ps.extend(p.index)

pred_counts = pd.Series(Counter(ps)).sort_values(ascending=False)[:25]

fig, ax = plt.subplots()

rng = range(len(pred_counts))

ax.bar(rng, pred_counts.values)
ax.set_xticks(rng)
ax.set_xticklabels(pred_counts.index, rotation=90)
ax.set_xlim(rng[0] - 1, rng[-1] + 1)
ax.set_ylabel('# of models')
fig.tight_layout()


# Mapping these UniProt IDs to gene names:

# In[8]:


ups = [i.split('_')[0] for i in pred_counts.index]
up_map(ups).sort_index()


# ## Best model results

# Now, let us retrieve the best scoring model predictors:

# In[9]:


best = predictors[np.argmax(accs)]

display(up_map([i.split('_')[0] for i in best.index]).sort_index())

fig, ax = plt.subplots()
rng = range(len(best))

ax.bar(rng, best.values, color='k')

ax.set_xlim(rng[0] -1, rng[-1] + 1)
ax.set_xticks(rng)
ax.set_xticklabels(best.index, rotation=90)

fig.tight_layout()


# **Note:** these results indicate potential predictor phosphosites for discriminating R and NR patients prior to treatment, none of the post-treatment phosphoproteomic data has been used in these models.

# ---

# ## Comparison against differential phosphorilation analysis

# Let's now compare our top predictors across all the models (top 25) to the differential expression (phosphorylation) results for the pre-teatment contrast:

# In[10]:


dex_r1nr1 = pd.read_csv('results/2_diff_exp/R1_NR1_ttop.csv', index_col=0)

sig_psites = dex_r1nr1.loc[[(abs(f) >= 1 and p <= 0.05)
                           for (f, p) in zip(dex_r1nr1['logFC'],
                                             dex_r1nr1['P.Value'])], :].index.values


# Now check the intersection of p-sites common to the significantly expressed and the top selected features of the models. NOTE: We consider significantly differentially expressed any p-site with $p\text{-val}\leq 0.05$ and $|\log(FC)| \geq 1$

# In[11]:


top_preds = pred_counts.index

df = up_map([n.split('_')[0] for n in set(top_preds).intersection(sig_psites)])

for p in top_preds:
    i = df.ACC == p.split('_')[0]
    df.loc[i, 'Psite'] = p
    
df


# ## Top models' predictors and overlap with overall models

# Let us now compare the predictors of the subset of models whose accuracy is over 80% and their overlap with the results over all the models.
# 
# First, we select the top models and count how many times each predictor appeared:

# In[12]:


top_model_preds = []
top_models = [p for (i, p) in enumerate(predictors) if accs[i]>=0.8]

print 'There are %d top models (ACC >= 0.8)' % len(top_models)

for p in top_models:
    top_model_preds.extend(p.index)

print 'Number of models each predictor appeared:'
top_mod_pred_counts = pd.Series(Counter(top_model_preds)).sort_values(ascending=False)[:25]
top_mod_pred_counts


# Now we compare out of these top models' predictors, which is their intersection with the overall top predictors

# In[13]:


plot = venn([set(top_mod_pred_counts.index), set(top_preds)],
            labels=['Top models predictors', 'Overall top predictors'])


# ## Ranking predictors against differential expression
# 
# Now we will take our overall models' top 25 predictors and rank them according to their corresponding differential expression (phosphorilation) analyses over the different contrasts, namely: R2vsR1 (effect of drug in responders), NR2vsNR1 (effect of drug in non-responders) and R1vsNR1 (difference between responders and non-responders prior to treatment)

# ### R1vsNR1

# In[14]:


r1nr1_tops = dex_r1nr1.loc[top_preds, :]
r1nr1_tops['abs_logFC'] = abs(r1nr1_tops.logFC)
r1nr1_tops['Significance'] = r1nr1_tops['P.Value'] <= 0.05
r1nr1_tops.sort_values(by='abs_logFC', ascending=False)


# ### R2vsR1

# In[15]:


dex_r2r1 = pd.read_csv('results/2_diff_exp/R2_R1_ttop.csv', index_col=0)

r2r1_tops = dex_r2r1.loc[top_preds, :]
r2r1_tops['abs_logFC'] = abs(r2r1_tops.logFC)
r2r1_tops['Significance'] = r2r1_tops['P.Value'] <= 0.05
r2r1_tops.sort_values(by='abs_logFC', ascending=False)


# ### NR2vsNR1

# In[16]:


dex_nr2nr1 = pd.read_csv('results/2_diff_exp/NR2_NR1_ttop.csv', index_col=0)

nr2nr1_tops = dex_nr2nr1.loc[top_preds, :]
nr2nr1_tops['abs_logFC'] = abs(nr2nr1_tops.logFC)
nr2nr1_tops['Significance'] = nr2nr1_tops['P.Value'] <= 0.05
nr2nr1_tops.sort_values(by='abs_logFC', ascending=False)


# ---
