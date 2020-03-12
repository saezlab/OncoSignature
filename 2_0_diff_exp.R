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
# Differential expression analysis
# ================================
#
# This script analyzes the differential expression between treated vs. untreted
# samples (responders and non-responders seperatedly) and between responders
# vs. non-responders before treatment

library(limma)
library(readr)

#----------------------------------- INPUT -----------------------------------#
data_dir <- 'data'
dataf_ex <- 'norm_data_ex.csv'
dataf_cl <- 'norm_data_cl.csv'
annotf <- 'response.csv'

out_dir <- 'results'
#-----------------------------------------------------------------------------#

# Loading the data and annotation
data_ex <- read.csv(paste(data_dir, dataf_ex, sep='/'), row.names=1)
data_cl <- read.csv(paste(data_dir, dataf_cl, sep='/'), row.names=1)
targets <- read_csv(paste(data_dir, annotf, sep='/'))
targets <- as.data.frame(targets[complete.cases(targets$annot), ])

# Removing sample 34 from ex-vivo samples (borderline R, non-significantly
# different from EC50 threshold)

data_ex = data_ex[, !grepl('34', colnames(data_ex))]
targets = targets[!grepl('34', targets$sample), ]

# ============================= Ex-vivo samples ============================= #
# Defining contrasts
subtargets = targets[startsWith(targets$sample, 'EX'), ]
f <- factor(subtargets$condition, levels=unique(subtargets$condition))
design <- model.matrix(~0 + f)
                            # Effect of treatment in responders
cont.matrix <- makeContrasts(R_TvsR_U=fR_T - fR_U,
                             # Effect pf treatment in non-responders
                             NR_TvsNR_U=fNR_T - fNR_U,
                             # Responders vs. non-responders prior to treatment
                             R_UvsNR_U=fR_U - fNR_U,
                             levels=design)

# Differential expression analysis
pfit <- lmFit(data_ex, design)
pfit <- contrasts.fit(pfit, cont.matrix)
pfit <- eBayes(pfit)

# Saving the results
subdir <- paste(out_dir, '2_diff_exp_ex', sep='/')
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

for(c in colnames(cont.matrix)){
    ttop <- topTable(pfit, coef=c, adjust='fdr', n=nrow(pfit))
    write.csv(ttop, paste(subdir, paste(c, 'ttop.csv', sep='_'), sep='/'))
}

# ALL treated vs ALL untreated
aux <- as.character(lapply(subtargets$sample, function(x){
    return(strsplit(x, '_')[[1]][3])
}))
f <- factor(aux, levels=unique(aux))
design <- model.matrix(~0 + f)
# Effect of treatment in responders
cont.matrix <- makeContrasts(TvsU=fT - fU,
                             levels=design)

# Differential expression analysis
pfit <- lmFit(data_ex, design)
pfit <- contrasts.fit(pfit, cont.matrix)
pfit <- eBayes(pfit)

# Saving the results
ttop <- topTable(pfit, coef='TvsU', adjust='fdr', n=nrow(pfit))
write.csv(ttop, paste(subdir, 'TvsU_ttop.csv', sep='/'))

# ============================ Cell line samples ============================ #
# Grouped
# Defining contrasts
subtargets = targets[!startsWith(targets$sample, 'EX'), ]
f <- factor(subtargets$condition, levels=unique(subtargets$condition))
design <- model.matrix(~0 + f)
                             # Effect of treatment in responders
cont.matrix <- makeContrasts(R_TvsR_U=fR_T - fR_U,
                             # Effect pf treatment in non-responders
                             NR_TvsNR_U=fNR_T - fNR_U,
                             # Responders vs. non-responders prior to treatment
                             R_UvsNR_U=fR_U - fNR_U,
                             levels=design)

# Differential expression analysis
pfit <- lmFit(data_cl, design)
pfit <- contrasts.fit(pfit, cont.matrix)
pfit <- eBayes(pfit)

# Saving the results
subdir <- paste(out_dir, '2_diff_exp_cl', sep='/')
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

for(c in colnames(cont.matrix)){
    ttop <- topTable(pfit, coef=c, adjust='fdr', n=nrow(pfit))
    write.csv(ttop, paste(subdir, paste(c, 'ttop.csv', sep='_'), sep='/'))
}
# By individual cell-lines
aux <- gsub('_._', '_', colnames(data_cl))
f <- factor(aux, levels=unique(aux))
design <- model.matrix(~0 + f)
cont.matrix <- makeContrasts(# Effect of treatment in GDM1
                             GDM1_TvsGDM1_U=fGDM1_T - fGDM1_U,
                             # Effect of treatment in MV411
                             MV411_TvsMV411_U=fMV411_T - fMV411_U,
                             # Effect of treatment in PL21
                             PL21_TvsPL21_U=fPL21_T - fPL21_U,
                             # Effect of treatment in NOMO1
                             NOMO1_TvsNOMO1_U=fNOMO1_T - fNOMO1_U,
                             # Effect of treatment in PL21
                             MV411_UvsPL21_U=fMV411_U - fPL21_U,
                             levels=design)

# Differential expression analysis
pfit <- lmFit(data_cl, design)
pfit <- contrasts.fit(pfit, cont.matrix)
pfit <- eBayes(pfit)

# Saving the results
subdir <- paste(out_dir, '2_diff_exp_cl', sep='/')
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

for(c in colnames(cont.matrix)){
    ttop <- topTable(pfit, coef=c, adjust='fdr', n=nrow(pfit))
    write.csv(ttop, paste(subdir, paste(c, 'ttop.csv', sep='_'), sep='/'))
}
