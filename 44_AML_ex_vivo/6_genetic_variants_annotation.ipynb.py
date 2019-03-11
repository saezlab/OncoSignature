
# coding: utf-8

# # Comparison between genetic and response annotations
# *Nicolàs Palacio-Escat*

# Here we will compare the annotation similarity between our response annotation (derived from our dose-repsponse data) and the genetic annotation over different genetic features of the samples.
# 
# Since we only have qualitative data available, in order to do this, we will take the different annotation features and compare the subsets of samples corresponding to each category against our subsets of R and NR samples. This similarity can be interpreted as a measure of how that genetic feature is related to the drug response.
# 
# But first, let us load the data and some packages

# In[1]:


import pandas as pd
from IPython.display import display, Image

from data_tools.iterables import similarity


# In[2]:


annot = pd.read_csv('data/genet_annot.csv', index_col=0)


# The annotation and features are as follows:

# In[3]:


annot


# Now we will compare the similarity for each class of each genetic annotation against our response classes.
# 
# This will be computed with three different similarity coefficients, namely:
# 
# - **Jaccard:** Basic similarity metric, basically is the intersection over the union of two sets.
# 
# $$
# s_J(A,B) = \frac{|A\cap B|}{|A\cup B|}
# $$
# 
# - **Sorensen-Dice:** Similar to Jaccard but gives higher weight to similarities between sets.
# 
# $$
# s_{SD}(A,B) = \frac{2|A\cap B|}{|A|+|B|}
# $$
# 
# - **Szymkiewicz–Simpson:** Similarity metric that corrects for differences the difference between set sizes.
# 
# $$
# s_{SS}(A,B) = \frac{|A\cap B|}{\min(|A|,|B|)}
# $$

# In[4]:


for col in annot.columns[1:]:
    string = 'Similarity with %s annotation' % col
    print string
    print '=' * len(string)
    test = annot[col]
    
    # Do not account for samples whose comparing annotation is NaN
    rs = annot[[a and b for (a, b) in zip(annot.ANNOT == 'R',
                                          pd.notna(annot[col]))]].index
    nrs = annot[[a and b for (a, b) in zip(annot.ANNOT == 'NR',
                                           pd.notna(annot[col]))]].index

    for c in set(test):
        if type(c) is not float: # Skip NaN's
            print 'Similarity coefficients for class %s:\n' % c
            print '* Jaccard:'
            print '    - Rs: %.4f' % similarity(rs, test[test == c].index,
                                                mode='j')
            print '    - NRs: %.4f' % similarity(nrs, test[test == c].index,
                                                 mode='j')
            print '* Sorensen-Dice:'
            print '    - Rs: %.4f' % similarity(rs, test[test == c].index,
                                                mode='sd')
            print '    - NRs: %.4f' % similarity(nrs, test[test == c].index,
                                                 mode='sd')
            print '* Szymkiewicz–Simpson:'
            print '    - Rs: %.4f' % similarity(rs, test[test == c].index,
                                                mode='ss')
            print '    - NRs: %.4f' % similarity(nrs, test[test == c].index,
                                                 mode='ss')
            print '\n'


# In order to visualize these results, the following figures have been manually generated in order to show the different similarity indexes across the different annotations.

# In[5]:


display(Image('results/jaccard.png'))
display(Image('results/sorensen_dice.png'))
display(Image('results/szymkiewicz_simpson.png'))


# In general there's no clear genetic annotation that can be used as a proxy of the drug response. Maybe NPM1 (specifically Type A) has some tendency NR's but taking into account the results of the remaining subcategories, it's not very clear if this relation is very significant.
# 
# Similarly, the lack of samples in many subcategories makes taking any robust conclusion impossible.
# 
# In summary, we cannot fully discard as well as we cannot prove the existence of relation between these genomic features and the response to the treatment. The answer to this question should be adressed in a separate study with more samples and maybe a deeper genomic annotation of these (e.g: GWAS).

# ---
