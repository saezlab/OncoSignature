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
# Biomarkers of interest
# ======================
#
# This script tries to find interesting biomarkers based on the differential
# expression results. Phosphosites of interest are defined as follows:
'''
Efficacy-predictive biomarkers:
===============================
 * Downregulated by the drug in responders while non- or up-regulated in
   the non-responders. At the same time, they are higher pre-treatment in
   responders compared to non-responders.
   {}

 * Upregulated by the drug in responders while non- or down-regulated in
   the non-responders. At the same time, they are lower pre-treatment in
   responders compared to non-responders.
   {}

Excluder biomarkers:
====================
 * Downregulated by the drug in non-responders while non- or up-regulated
   in the responders. At the same time, they are lower pre-treatment in
   responders compared to non-responders (i.e. higher in non-responders).
   {}

 * Upregulated by the drug in non-responders while non- or down-regulated
   in the responders. At the same time, they are higher pre-treatment in
   responders compared to non-responders (i.e. lower in non-responders).
   {}
'''
import os

import pandas as pd

from data_tools.databases import up_map

#----------------------------------- INPUT -----------------------------------#
parent_dir = 'results'
#-----------------------------------------------------------------------------#


def is_significant(df, pval=0.05, logfc=1):
    '''
    Determines if any given pairs of p-values and log fold changes are
    significant.

    :arg pandas.DataFrame df:
        Data frame containing the p-values (first column) and the log
        fold changes (second column). The rest of columns or indexes
        are ignored.
    :arg float pval:
        Optional, 0.05 by default. Sets the threshold for significance
        of p-values.
    :arg float logfc:
        Optional, 1 by default. Sets the threshold for significance of
        log fold changes.

    :return list:
        Contains the logical values for significance of each entry
        (True/False).
    '''

    sig_p = [v <= pval for v in df.iloc[:, 0]]
    sig_fc = [abs(v) >= logfc for v in df.iloc[:, 1]]

    return [p and f for (p, f) in zip(sig_p, sig_fc)]

def translate(ids):
    '''
    Given a set of p-site identifiers, replaces the UniProt AC by its
    corresponding gene symbol. It is assumed that the UniProt AC is the
    first element of the identifier followed by an underscore ('_') and
    whatever else.

    :arg set ids:
        The IDs to translate from UniProt AC to gene symbols.

    :return list:
        The same set of IDs but with gene symbols
    '''

    return [i.replace(i.split('_')[0],
                      up_map([i.split('_')[0]]).iloc[0, 1]) for i in ids]


usedirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir)
           if d.startswith('2_diff_exp')]

for dir_ in usedirs:
    print('Biomarkers for %s'
          % 'cell lines' if dir_.endswith('cl') else 'ex vivo samples')

    # Reading DEX results
    fnames = [fname for fname in os.listdir(dir_)
              if (fname.endswith('_ttop.csv') and not fname.startswith('sig'))]
    ttops = dict((fname.replace('_ttop.csv', ''),
                  pd.read_csv(os.path.join(dir_, fname), index_col=0))
                 for fname in fnames)

    # Extracting the significant p-sites
    sig_sites = dict((k, dict([('up', set()),
                               ('down', set())])) for k in ttops.keys())

    for k, v in ttops.items():
        sig = v.loc[is_significant(v[['P.Value', 'logFC']]), :]
        sig_sites[k]['up'].update(sig.index[sig['logFC'] > 0])
        sig_sites[k]['down'].update(sig.index[sig['logFC'] < 0])

    # Case 1.1
    c11 = (sig_sites['R_TvsR_U']['down']
           - sig_sites['NR_TvsNR_U']['down']
           & sig_sites['R_UvsNR_U']['up'])

    # Case 1.2
    c12 = (sig_sites['R_TvsR_U']['up']
           - sig_sites['NR_TvsNR_U']['up']
           & sig_sites['R_UvsNR_U']['down'])

    # Case 2.1
    c21 = (sig_sites['NR_TvsNR_U']['down']
           - sig_sites['R_TvsR_U']['down']
           & sig_sites['R_UvsNR_U']['down'])

    # Case 2.2
    c22 = (sig_sites['NR_TvsNR_U']['up']
           - sig_sites['R_TvsR_U']['up']
           & sig_sites['R_UvsNR_U']['up'])

    # Printing the results
    if __name__ == '__main__':
        print(__doc__.format(*[', '.join(x) if x else 'None'
                               for x in map(translate, [c11, c12, c21, c22])]))
