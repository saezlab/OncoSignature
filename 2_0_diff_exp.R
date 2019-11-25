# Copyright (C) 2019 Nicolàs Palacio
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
dataf <- 'norm_data.csv'
annotf <- 'response.csv'

out_dir <- 'results'
#-----------------------------------------------------------------------------#

# Loading the data and annotation
data <- read.csv(paste(data_dir, dataf, sep='/'), row.names=1)
targets <- read_csv(paste(data_dir, annotf, sep='/'))
targets <- as.data.frame(targets[complete.cases(targets$annot), ])

# Defining contrasts
f <- factor(targets$condition, levels=unique(targets$condition))
design <- model.matrix(~0 + f)
                            # Effect of treatment in responders
cont.matrix <- makeContrasts(R_TvsR_U=fR_T - fR_U, 
                             # Effect pf treatment in non-responders
                             NR_TvsNR_U=fNR_T - fNR_U,
                             # Responders vs. non-responders prior to treatment
                             R_UvsNR_U=fR_U - fNR_U,
                             levels=design)

# Differential expression analysis
pfit <- lmFit(data, design)
pfit <- contrasts.fit(pfit, cont.matrix)
pfit <- eBayes(pfit)

# Saving the results
subdir <- paste(out_dir, '2_diff_exp', sep='/')
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

for(c in colnames(cont.matrix)){
    ttop = topTable(pfit, coef=c, adjust='fdr', n=nrow(pfit))
    write.csv(ttop, paste(subdir, paste(c, 'ttop.csv', sep='_'), sep='/'))
}
