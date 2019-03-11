library(limma)

# Set up CWD and load (normalized) data and target matrix
#setwd('/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/44_AML_ex_vivo/')

data = read.csv('data/norm_data.csv', row.names=1)
targets = read.csv('data/targets.csv')

# Defining contrasts
f = factor(targets$condition, levels=unique(targets$condition))
design = model.matrix(~0 + f)
cont.matrix = makeContrasts(R2_R1=fR_2 - fR_1, # Effect of treatment in responders
                            NR2_NR1=fNR_2 - fNR_1, # Effect pf treatment in non-responders
                            R1_NR1=fR_1 - fNR_1, # Difference in responders wrt. non-responders prior to treatment
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
