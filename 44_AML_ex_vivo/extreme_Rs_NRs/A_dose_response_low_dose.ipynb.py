
# coding: utf-8

# # Dose-response analysis of patient samples
# 
# We have 44 samples from patients cultivated *ex vivo* which were treated with different doses of a drug. For each sample and dose, we have four replicates.
# 
# The main idea is to classify these patient samples as responders and non-responders to the treatment. In order to do this, we will fit a scaled Hill curve to each sample. The model is defined as follows:
# 
# $$
# S = m\frac{D^n}{k^n+D^n}
# $$
# 
# Where $S$ is the surviving fraction and $D$ is the dose and $k$, $m$ and $n$ are the function parameters. It can be seen that the model is a mixture between a Hill function and the Michaelis-Menten or Monod kinetics.

# ## Importing packages

# In[1]:


import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import colors
from IPython.display import display
from scipy.optimize import least_squares

from data_tools.models import DoseResponse
from data_tools.iterables import chunk_this


# In[2]:


get_ipython().magic(u'matplotlib inline')

warnings.filterwarnings('ignore')


# ## Loading dose-response data

# In[3]:


df_raw = pd.read_table('data/raw_survival_data.tab')

df_raw.head()


# ## Data normalization
# 
# Let us compute the normalized response (survival fraction) of each sample (relative to $0nM$ dose for that sample).

# In[4]:


rep_cols = [i for i in df_raw.columns if i.startswith('Raw')]
norm_samples = {}

for s in df_raw.groupby('Sample #'):
    sample, data = s

    ref = data.loc[data['slx dose (nM)'] == 0, rep_cols]
    norm = data[rep_cols].values / ref.values
    
    df = pd.DataFrame(norm[2:],
                      index=data['slx dose (nM)'].values[2:],
                      columns=map('R{}'.format, range(1, 5)))
    
    norm_samples[sample] = df


# ---

# > **IMPORTATNT NOTE:**
# 
# > After the normalization and taking into account that the treatment is either not affecting or killing cells, one would not expect survival values much greater than 100% at any point. Nevertheless, if we check all the normalized samples for values higher than 2, we find some outliers:

# In[5]:


for k, v in norm_samples.items():
    if (v.dropna().values > 2).any():
        print 'Sample #%d' %k
        display(v.dropna()[v.dropna().values > 2])


# > On sample #1 we see two different replicates in two different doses whose values are $>20$ ($>2000\%$ cells relative to $0nM$ dose). Taking into account the other replicates and neighboring doses, these two measurements are clear outliers.
# 
# > On sample #33, one could say that the $\sim2.1$ might be feasible, but observing the other replicates, we will also consider this measurement an outlier.
# 
# > In order to improve the fitting of our model, we will discard these three measurements.

# In[6]:


norm_samples[1][norm_samples[1] > 20] = np.nan
norm_samples[33][norm_samples[33] > 2] = np.nan


# ---

# Now, let us plot the normalized dose-response data for each patient. For reference, we also show two dashed lines at 100% survival (black) and 50% survival (red).

# In[7]:


def plot(ax, df):
    ax.errorbar(df.index, df.mean(axis=1), yerr=df.std(axis=1), c='k')
    ax.scatter(np.repeat(df.index.values, df.shape[1]),
               df.values.reshape(-1),
               c=map('C{}'.format, range(df.shape[1])*len(df)))

xmin = min(df_raw['slx dose (nM)'].unique()[2:])
xmax = max(df_raw['slx dose (nM)'].unique()[2:])

n = 6
groups = chunk_this(norm_samples.keys(), n)
m = len(groups)

fig, ax = plt.subplots(m, n, figsize=(3 * n, 3 * m),
                       sharex=True, sharey=True)

# j is vertical index, g is labels of group samples
for j, g in enumerate(groups):
    # i is horizontal index, s is sample id
    for i, s in enumerate(g):
        df = norm_samples[s]

        plot(ax[j, i], df)

        # Plot dashed lines at S=1 and S=0.5 for reference
        ax[j, i].plot([xmin, xmax], [1, 1],
                      '--k', alpha=.5)
        ax[j, i].plot([xmin, xmax], [.5, .5],
                      '--r', alpha=.5)

        ax[j, i].set_title('Sample #%d' %s)

        if j + 1 == m:
            ax[j, i].set_xlabel(r'$D$')

        if i == 0:
            ax[j, i].set_ylabel(r'$S$')

# Remove axis of empty plots in last line
idxs = [i for i in range(n) if len(ax[m - 1, i].get_lines()) == 0]
for i in idxs:
    ax[m - 1, i].axis('off')


fig.legend([Patch(color='C%d' %i) for i in range(df.shape[1])],
           df.columns)
fig.tight_layout()
fig.savefig('results/A_dose_response_low/norm_survival_low.pdf')


# From these plots we can already have an idea of which samples are responding more or less to the drug.

# ## Setting up the model fittings
# 
# In order to fit our model, we will use a non-linear least squares algorithm. Hence, let us define our model.

# In[8]:


def hill(x, k, m, n):
    return m * x ** n / (k ** n + x ** n)


# ## Fitting the models

# In[9]:


fig, ax = plt.subplots(m, n, figsize=(3 * n, 3 * m),
                       sharex=True, sharey=True)
rng = np.linspace(xmin, xmax, 1000)

# Parameter order: k, m, n
x0 = [1e3, 1, -1]
bounds = ([0, 0, -np.inf], [np.inf, np.inf, 0])

fit_params = {}
ec50s = {}

# j is vertical index, g is labels of group samples
for j, g in enumerate(groups):
    # i is horizontal index, s is sample id
    for i, s in enumerate(g):
        xdata = norm_samples[s].index.values[:-1]
        ydata = norm_samples[s].mean(axis=1).values[:-1]

        model = DoseResponse(xdata, ydata, x0=x0, bounds=bounds)

        fit_params[s] = model.params
        ec50s[s] = model.ec()

        ax[j, i].scatter(xdata, ydata, label='Data')
        ax[j, i].plot(rng, hill(rng, *model.params), 'k', label='Fit')
        
        # Plot dashed lines at S=1 and S=0.5 for reference
        ax[j, i].plot([xmin, xmax], [1, 1],
                      '--k', alpha=.5)
        ax[j, i].plot([xmin, xmax], [.5, .5],
                      '--r', alpha=.5)

        ax[j, i].set_title('Sample #%d' %s)
        ax[j, i].text(0, 1.2, r'$k=%.2f$ | $m=%.2f$ | $n=%.2f$'
                      %tuple(model.params), fontsize=9)

        if j + 1 == m:
            ax[j, i].set_xlabel(r'$D$')

        if i == 0:
            ax[j, i].set_ylabel(r'$S$')

# Remove axis of empty plots in last line
idxs = [i for i in range(n) if len(ax[m - 1, i].get_lines()) == 0]
for i in idxs:
    ax[m - 1, i].axis('off')

fig.tight_layout()
fig.savefig('results/A_dose_response_low/fit_survival_low.pdf')


# As can be seen, we could already classify responders and non-responders by only checking the value of $n$. Also known as the Hill coefficient, in the context of dose-response curves, this parameter is many times used as a measure of hypersensitivity as it defines the steepness of the curve.
# 
# In our case, we could say that if $n>-0.05$ the sample could be considered non-responder.
# 
# Nevertheless, when reporting dose-dependent responses, it is many times expressed in the form of $EC_{50}$ (half maximal effective concentration). This is, the concentration ($D$ in our case) at which the response ($S$) is half of the maximal response. Since our normalization is bound to $S\in[0, 1]$ then $S(EC_{50})=0.5$. Following our model definition above, the formula for our $EC_{50}$ is derived as follows:
# 
# $$
# \frac{1}{2}=m\frac{EC_{50}^n}{k^n+EC_{50}^n}\\
# k^n+EC_{50}^n=2mEC_{50}^n\\
# k^n=EC_{50}^n(2m-1)
# $$
# 
# Finally:
# 
# $$
# EC_{50}=\left(\frac{k^n}{2m-1}\right)^{\frac{1}{n}}
# $$

# In[10]:


fig, ax = plt.subplots(figsize=(8, 5))

xmin = min(fit_params.keys())
xmax = max(fit_params.keys())

rng = range(xmin, xmax + 1)

for k, v in fit_params.items():
    
    val = ec50s[k]
    
    if val > 1e25:
        c = 'r'
        
    elif val > 1e4:
        c = 'y'
        
    elif val > 1e3:
        c = 'g'
        
    else:
        c = 'b'
        
    ax.bar(k, val if val != np.inf else 2e96, color=c)

ax.set_xlim(xmin - 1, xmax + 1)

ax.set_yscale('log')
ax.set_ylim(0, 1e5)
ax.set_xticks(rng)
ax.set_xticklabels(fit_params.keys(),
                   rotation=90, ha='center',
                   va='center')

ax.set_xlabel('Sample #')
ax.set_ylabel(r'$EC_{50}$')

ax.legend([Patch(color=i) for i in ['r', 'y', 'g', 'b']],
          [r'$EC_{50}>10^{25}$', r'$10^{25}>EC_{50}>10^{4}$',
           r'$10^{4}>EC_{50}>10^{3}$', r'$EC_{50}<10^{3}$'], loc=0)

fig.tight_layout()
fig.savefig('results/A_dose_response_low/ec50_low.pdf')


# In[11]:


pd.Series(np.log10(ec50s.values())).hist(bins=44)


# The $EC_{50}$ values obtained from the fitted models of each sample have been classified as follows:
# - In red, values over $10^{25}nM$. These are definitely non-responders.
# - In yellow, values between $10^{25}nM$ and $10^{4}nM$. These are responding slightly to the treatment but *could* also be regarded as non-responders.
# - In green, values between $10^{4}nM$ and $10^{3}nM$. These samples are definitely responding to the treatment and are classified as responders.
# - In blue, values below $10^{3}nM$. These are samples which are also responders AND we have experimental proof of their $EC_{50}$.

# In[12]:


sample_class = pd.DataFrame(index=range(1, len(ec50s)), columns=['sample', 'condition'])

for k, v in ec50s.items():
    sample_class.loc[k, 'sample'] = str(k)
    sample_class.loc[k, 'condition'] = 'R' if v <= 1e3 else 'NR'


# In[13]:


sample_class


# In[14]:


sample_class.to_csv('data/sample_annotated.csv', index=False)


# ---

# ## Sample annotation
# 
# Let's check which samples are significantly different from $10^3nM$ in order to consider them for the analysis

# In[15]:


missing = [1, 4, 7, 12, 26]
samples = pd.DataFrame(ec50s.values(), index=ec50s.keys(), columns=['EC50'])

for s in missing:
    samples.drop(s, inplace=True)


# In[16]:


from scipy.stats import t


# In[17]:


samples['logEC50'] = np.log10(samples.EC50)


# In[18]:


def tval(x, ref=1e3):
    den = np.std(x) / len(x)
    
    return [(i - ref) / den for i in x]


# In[19]:


samples['tval'] = tval(samples.logEC50, ref=3)
samples['pval'] = [t.sf(np.abs(v), len(samples) - 1) * 2 for v in samples.tval]


# Let us discard those samples whose $EC_{50}$ is non significantly different from $10^3$ (with $\alpha=0.01$)

# In[20]:


samples[samples.pval > 0.01]


# In[21]:


annot = []

for s in samples.iterrows():
    if s[1].pval > 0.01:
        annot.append(np.nan)
    else:
        annot.append('R' if s[1].EC50 <= 1e3 else 'NR')

samples['annotation'] = annot
samples['sample'] = samples.index
samples.to_csv('data/sample_annotated.csv', index=False)

