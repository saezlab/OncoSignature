import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.chdir('/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/44_AML_ex_v'
         'ivo/extreme_Rs_NRs/')
in_dir = 'lasso_results'
out_dir = os.path.join(in_dir, 'acc_80')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
%matplotlib inline

n = 100 # NUMBER OF MODELS
thr = 0.8 # ACC threshold

# Loading results data
#=============================================================================#
# Accuracy measures
accs = pd.read_csv(os.path.join(in_dir, 'accs.csv'), header=None)
accs = accs.iloc[:, 1]
usemodels = accs[accs>=thr].index.tolist()
#accs = accs[usemodels]
print accs.mean()
print accs.std()

fig, ax = plt.subplots()
ax.bar(range(n), accs)
ax.set_xlabel('Model #')
ax.set_ylabel('ACC')
ax.set_xlim(-1, n)
fig.tight_layout()
fig.savefig(os.path.join(in_dir, 'accs.png'))

fig, ax = plt.subplots()
ax.hist(accs.values, bins=7)
ax.set_title('Average ACC = %.3f' % accs.mean())
fig.savefig(os.path.join(in_dir, 'acc_hist.png'))

# Top models (ACC >= thr):
print usemodels
#=============================================================================#
# Intercepts
intercepts = pd.read_csv(os.path.join(in_dir, 'intercepts.csv'), header=None)
intercepts = intercepts.iloc[usemodels, 1]
print intercepts.sum()
#=============================================================================#
# Predictor coefficients
predictors = pd.read_csv(os.path.join(in_dir, 'predictors.csv'), index_col=0)
#predictors.head()
predictors = predictors.loc[:, map(str, usemodels)]

counts = pd.notna(predictors).sum(axis=1)
counts.sort_values(ascending=False, inplace=True)
#counts.head(25)

top_pred = predictors.loc[counts.head(25).index, :]
top_pred['Mean'] = top_pred.mean(axis=1)
top_pred.sort_values(by='Mean', inplace=True)
fpred = []
names = []
for (p, vals) in top_pred.iterrows():
    fpred.append(vals[pd.notna(vals)].tolist())
    names.append(p)

rng = range(1, len(fpred) + 1)

fig, ax = plt.subplots()
ax.boxplot(fpred)
ax.plot([0, len(rng) + 1], [0, 0], 'k--', alpha=0.5)

ax.set_xlim(0, len(rng) + 1)
ax.set_xticks(rng)
ax.set_xticklabels(names, rotation=90)
ax.set_title(r'Top 25 pred. coef. (Models w/ $ACC\geq %.2f$)' %thr)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'box_top_coef.png'))

#fig

# Mean, median and SD of top coefficients
pred_dist = pd.DataFrame(columns=['Mean', 'SD'], index=names)
pred_dist['Mean'] = map(np.mean, fpred)
pred_dist['SD'] = map(np.std, fpred)

rng = range(len(pred_dist))

fig, ax = plt.subplots()
ax.errorbar(rng, pred_dist.Mean, yerr=pred_dist.SD, fmt='_',
            label=r'Mean $\pm$ SD')
ax.scatter(rng, map(np.median, fpred), marker='x', c='C1', label='Median')
ax.plot([-1, len(rng)], [0, 0], 'k--', alpha=0.5)
ax.set_xlim(-1, len(rng))
ax.set_xticks(rng)
ax.set_xticklabels(pred_dist.index, rotation=90)
ax.legend()
ax.set_title(r'Top 25 pred. coef. (Models w/ $ACC\geq %.2f$)' %thr)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'mean_top_coef.png'))


#output = pred_dist
#output.loc['Intercept', :] = [intercepts.mean(), intercepts.std()]
#output.to_csv(os.path.join(out_dir, 'top25_mean_coef.csv'))

#=============================================================================#
# Prediction probabilities (if >=0.5 -> R)
#pprobs = pd.read_csv(os.path.join(sdir, 'pprobs.csv'), index_col=0)
#(pprobs.mean(axis=1)>=0.5).astype(int).sort_index()
#=============================================================================#
# Train/test sample subsets
#testsets = np.load(os.path.join(sdir, 'testsets.npy'))
#trainsets = np.load(os.path.join(sdir, 'trainsets.npy'))
#=============================================================================#
