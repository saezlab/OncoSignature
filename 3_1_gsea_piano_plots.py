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
# Gene set enrichment analysis plots
# ==================================
#
# This script generates the plots to visualize the results of the GSEA

import os

import pandas as pd
import matplotlib.pyplot as plt

from data_tools.plots import piano_consensus
from data_tools.plots import venn
from data_tools.iterables import subsets

#----------------------------------- INPUT -----------------------------------#
dir_ = 'results/3_gsea'
#-----------------------------------------------------------------------------#

def titlemaker(s):
    '''
    Creates the figure title based on the file name of the PIANO result table.

    :arg str s:
        The file name containing the information about directionality class and
        direction (if applicable). It is assumed to have the format:
        "{1}_{2}_rank.txt" where {1} denotes the directionality class ("dist",
        "mix" or "non") and {2} the direction ("up", "dw" or none)

    :return str:
        The formatted title for the figure
    '''

    ssplit = s.split('_')

    # Check first element of the name for the classification type
    if ssplit[0] == 'dist':
        c = 'Distinct'

    elif ssplit[0] == 'mix':
        c = 'Mixed'

    elif ssplit[0] == 'non':
        c = 'Non'

    else:
        c = None

    # Check second element of the name for directionality
    if ssplit[1] == 'up':
        d = '(up)'

    elif ssplit[1] == 'dw':
        d = '(down)'

    else:
        d = None

    return ('%s-directional %s individual ranks based on p-values' %(c, d) if d
            else '%s-directional individual ranks based on p-values' %c)


supdirs = list(os.walk(dir_))[0][1]

for supdir in supdirs:
    subdirs = [d for d in os.listdir(os.path.join(dir_, supdir))
               if '.' not in d]

    venn_sets = {'NR up':set(),
                 'NR down':set(),
                 'R up':set(),
                 'R down':set()}

    for sdir in subdirs:
        ssdir = os.path.join(dir_, supdir, sdir)

        # Listing all files containing the rank matrices
        files = [f for f in os.listdir(ssdir) if f.endswith('_rank.txt')]

        # Setting up x-limit common for plots on same group
        lim = 200 if supdir == 'ex' else 650

        for f in files:
            df = pd.read_csv(os.path.join(ssdir, f), sep='\t')
            height = len(df) / 5

            fig = piano_consensus(df, figsize=[10, height],
                                  title=titlemaker(f))
            fig.gca().set_xlim(0, lim)
            fig.savefig(os.path.join(ssdir, f.split('.')[0] + '.pdf'))
            plt.close('all')

            if sdir == 'NR_T_NR_U':
                if f == 'dist_up_rank.txt':
                    venn_sets['NR up'].update(df.index.tolist())

                elif f == 'dist_dw_rank.txt':
                    venn_sets['NR down'].update(df.index.tolist())

                else:
                    pass

            elif sdir == 'R_T_R_U':
                if f == 'dist_up_rank.txt':
                    venn_sets['R up'].update(df.index.tolist())

                elif f == 'dist_dw_rank.txt':
                    venn_sets['R down'].update(df.index.tolist())

                else:
                    pass

            else:
                pass


    print(venn_sets.keys())

    aux = subsets(list(venn_sets.values()))

    for k, v in aux.items():
        print([list(venn_sets.keys())[i] for i, n in enumerate(k) if n == '1'])
        print(v)

    venn(list(venn_sets.values()), list(venn_sets.keys()), title='GSEA overlaps',
         filename=os.path.join(dir_, supdir, 'venn_piano.pdf'))
