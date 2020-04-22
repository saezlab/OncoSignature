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

import numpy as np
import pandas as pd

from data_tools.databases import up_map


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


#----------------------------------- INPUT -----------------------------------#
parent_dir = 'results/perseus'
#-----------------------------------------------------------------------------#
usedirs = ['ExVivo']

for dir_ in usedirs:
    print('Biomarkers for %s' % dir_)
    fulldir = os.path.join(parent_dir, dir_)

    # Loading the results from the differential expression analysis
    dex_files = [f for f in os.listdir(fulldir) if (f.startswith('DEX_')
                 and 'All' not in f)]

    # Extracting the significant p-sites
    sig_sites = dict((k.lstrip('DEX_').rstrip('.txt'), dict([('up', set()),
                               ('down', set())])) for k in dex_files)

    for dexf in dex_files:
        name = dexf.lstrip('DEX_').rstrip('.txt')

        df = pd.read_csv(os.path.join(fulldir, dexf), sep='\t')
        #df.set_index('Fusion', drop=True, inplace=True)
        #df.index.name = None
        #df = df.loc[:, ['Significant', '(-)log10(P)', 'Difference']]
        sig_sites[name]['up'].update([i[1]['Fusion'] for i in df.iterrows()
                                      if (i[1]['Significant'] == '+'
                                          and i[1]['Difference'] > 0)])
        sig_sites[name]['down'].update([i[1]['Fusion'] for i in df.iterrows()
                                        if (i[1]['Significant'] == '+'
                                            and i[1]['Difference'] < 0)])

    # Case 1.1
    c11 = (sig_sites['Responder']['down']
           - sig_sites['NonResponder']['down']
           & sig_sites['RvsNRpretrea']['up'])

    # Case 1.2
    c12 = (sig_sites['Responder']['up']
           - sig_sites['NonResponder']['up']
           & sig_sites['RvsNRpretrea']['down'])

    # Case 2.1
    c21 = (sig_sites['NonResponder']['down']
           - sig_sites['Responder']['down']
           & sig_sites['RvsNRpretrea']['down'])

    # Case 2.2
    c22 = (sig_sites['NonResponder']['up']
           - sig_sites['Responder']['up']
           & sig_sites['RvsNRpretrea']['up'])

    # Printing the results
    if __name__ == '__main__':
        print(__doc__.format(*[', '.join(x) if x else 'None'
                               for x in map(translate, [c11, c12, c21, c22])]))
