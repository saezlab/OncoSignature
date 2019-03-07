import os

import pandas as pd

from data_tools.plots import piano_consensus
from data_tools.plots import venn


def titlemaker(s):
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


# Select the directory where to find the piano rank matrices
dir_ = (#'/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/2+2_AML_patients/'
        'results/3_gsea/')

os.chdir(dir_)

subdirs = list(os.walk(os.getcwd()))[0][1]
venn_sets = {'NR up':set(),
             'NR down':set(),
             'R up':set(),
             'R down':set()}

for sdir in subdirs:
    os.chdir(dir_ + sdir + '/')

    # Listing all files containing the rank matrices
    files = [f for f in os.listdir(os.getcwd()) if f.endswith('_rank.txt')]

    # File names are assumed to be {1}_{2}_rank.txt where {1} denotes the
    # directionality class (dist, mix or non) and {2} the direction (up, dw)
    # if applicable
    for f in files:
        df = pd.read_table(f)
        height = len(df) / 5

        fig = piano_consensus(df, figsize=[10, height], title=titlemaker(f),
                              filename=f.split('.')[0] + '.pdf')

        if sdir == 'NR2_NR1':
            if f.endswith('up_rank.txt'):
                venn_sets['NR up'].update(df.index.tolist())

            elif f.endswith('dw_rank.txt'):
                venn_sets['NR down'].update(df.index.tolist())

            else:
                pass

        elif sdir == 'R2_R1':
            if f.endswith('up_rank.txt'):
                venn_sets['R up'].update(df.index.tolist())

            elif f.endswith('dw_rank.txt'):
                venn_sets['R down'].update(df.index.tolist())

            else:
                pass

        else:
            pass

venn(venn_sets.values(), venn_sets.keys(), filename='../venn_piano.pdf')
