import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

os.chdir('/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/validation')

# Loading validation data
data = pd.read_csv('data/norm_data.csv', index_col=0)
data[pd.isna(data)] = 0

usecols = [c for c in data.columns if c.endswith('B')]
data = data[usecols]
data.columns = [c[:-1] for c in data.columns]
#data.head()

# Loading ensemble of models
sdir = '../44_AML_ex_vivo/extreme_Rs_NRs/lasso_results'
predictors = pd.read_csv(os.path.join(sdir, 'predictors.csv'), index_col=0)
predictors[pd.isna(predictors)] = 0
#predictors.head()
intercepts = pd.read_csv(os.path.join(sdir, 'intercepts.csv'), index_col=0,
                         header=None)
#intercepts.head()

accs = pd.read_csv(os.path.join(sdir, 'accs.csv'), index_col=0, header=None)
#accs.head()

# Reindexing data and predictors to same phosphosites and order
a, b = map(set, [predictors.index, data.index])
print map(len, [a, b])

common = list(a.intersection(b))
len(common)

predictors = predictors.loc[common, :]
data = data.loc[common, :]

print data.shape, predictors.shape
#data

# Applying ALL the models to the samples
res = pd.DataFrame(index=data.columns, columns=range(predictors.shape[1]))

for n in range(predictors.shape[1]):
    m = predictors.iloc[:, n].values
    b = intercepts.iloc[n].values[0]
    res[n] = 1 / (1 + np.exp((data.values.T * m).sum(axis=1) + b))

# Applying models with ACC>80%
predictors80 = predictors.loc[:, (accs>=0.8).values.flatten()]
res80 = pd.DataFrame(index=data.columns, columns=range(predictors80.shape[1]))

for n in range(predictors80.shape[1]):
    m = predictors80.iloc[:, n].values
    b = intercepts.iloc[n].values[0]
    res80[n] = 1 / (1 + np.exp((data.values.T * m).sum(axis=1) + b))

# Plotting the results
res['med'] = res.median(axis=1)
res.sort_values(by='med', inplace=True)
res80 = res80.loc[res.index, :]

res.drop('med', axis=1, inplace=True)

rng = range(data.shape[1])
fig, ax = plt.subplots()
ax.violinplot(res, showmedians=True, widths=0.75, positions=rng)
ax.violinplot(res80, showmedians=True, widths=0.75, positions=rng)
ax.plot([rng[0] - 1, rng[-1] + 1], [0.5, 0.5], 'k--', alpha=0.5)

ax.set_xticks(rng)
ax.set_xticklabels(data.columns)
ax.set_xlabel('Sample')
ax.set_ylabel('Prediction probability')
ax.set_xlim([rng[0] - 1, rng[-1] + 1])
ax.legend([Line2D([0], [0], color='C0'), Line2D([0], [0], color='C1')],
          ['ALL', r'$ACC\geq 80%$'], loc='upper left')
fig.tight_layout()
fig

fig.savefig('results/pred_prob.png')
