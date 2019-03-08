library(limma)

# Set up CWD and load (normalized) data and target matrix
#setwd('/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/4x_AML_cell_lines/')

data = read.csv('data/norm_data.csv', row.names=1)
targets = read.csv('data/targets.csv')

# Defining contrasts
f = factor(targets$condition, levels=unique(targets$condition))
design = model.matrix(~0 + f)
cont.matrix = makeContrasts(# Effect of treatment in GDM1
                            GDM1.2_GDM1.1=fKPT330_GDM1 - fCtrl_GDM1, 
                            # Effect of treatment in MV411
                            MV411.2_MV411.1=fKPT330_MV411 - fCtrl_MV411,
                            # Effect of treatment in NOMO1
                            NOMO1.2_NOMO1.1=fKPT330_NOMO1 - fCtrl_NOMO1,
                            # Effect of treatment in PL21
                            PL21.2_PL21.1=fKPT330_PL21 - fCtrl_PL21,
                            # Effect of treatment in responders
                            R.2_R.1=(fKPT330_GDM1 + fKPT330_MV411) - (fCtrl_GDM1 + fCtrl_MV411),
                            # Effect pf treatment in non-responders
                            NR.2_NR.1=(fKPT330_NOMO1 + fKPT330_PL21) - (fCtrl_NOMO1 + fCtrl_PL21),
                            # Difference in responders wrt. non-responders prior to treatment
                            R.1_NR.1=(fCtrl_GDM1 + fCtrl_MV411) - (fCtrl_NOMO1 + fCtrl_PL21),
                            levels=design)

# Differential expression analysis
pfit = lmFit(data, design)
pfit = contrasts.fit(pfit, cont.matrix)
pfit = eBayes(pfit)

# Save the results
subdir = 'results/2_diff_exp/'
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

for(c in colnames(cont.matrix)){
    ttop = topTable(pfit, coef=c, adjust='fdr', n=nrow(pfit))
    write.csv(ttop, paste(subdir, paste(c, 'ttop.csv', sep='_'), sep=''))
}
