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
# Data normalization pipeline
# ===========================
#
# Data cleaning + normalization + batch correction

library(readr)
library(dplyr)
library(vsn)
library(limma)

library(omicToolsTest)

#' Standarizes the column names
#' 
#' @param col [chracter] The column name to standarize
#' 
#' @return [character] The standarized column name
standarize_name <- function(col){
    # Sample comes from ex-vivo experiments
    if(startsWith(col, 'Set')){
        col <- gsub('Set.', 'EX', col)
        treat <- substr(col, nchar(col), nchar(col))
        
        if(treat == 'A'){# Treated
            return(gsub(treat, '_T', col))
        }
        else{# Untreated
            return(gsub(treat, '_U', col))
        }
    }
    # Sample comes from the cell-line experiments
    else{
        spl <- strsplit(col, '_')[[1]]
        rep <- substr(spl[3], nchar(spl[3]), nchar(spl[3]))
        
        if(spl[2] == 'DMSO'){
            treat <- 'U'
        }
        else{
            treat <- 'T'
        }
        return(paste(spl[1], rep, treat, sep='_'))
    }
}

#----------------------------------- INPUT -----------------------------------#
data_dir <- 'data'
dataf <- '20191107_expandSites_table_Pt_CellLines_basicfiltered_LocProb.txt'
annotf <- 'response.csv'

out_dir <- 'results'
#-----------------------------------------------------------------------------#

# Creating output directory
ifelse(!dir.exists(out_dir), dir.create(out_dir, recursive=TRUE), FALSE)

# Set up CWD and load raw data
raw_data <- read_delim(paste(data_dir, fname, sep='/'), sep='\t',
                       escape_double=FALSE, trim_ws=TRUE, na='NaN')
raw_data <- as.data.frame(raw_data)

# Locate the intensity columns within the raw data file
sample_cols <- colnames(raw_data)[(startsWith(colnames(raw_data), 'Set')
                                   & !endsWith(colnames(raw_data), 'NormC'))
                                  | grepl('Rep\\d', colnames(raw_data))]
raw_df <- as.data.frame(raw_data[, sample_cols])

# Simplifying column and row names
colnames(raw_df) <- as.character(lapply(colnames(raw_df), standarize_name))
rownames(raw_df) <- paste(paste(raw_data[, 'T: Protein'],
                                raw_data[, 'C: Amino acid'], sep='_'),
                          raw_data[, 'N: Position'],
                          raw_data[, 'C: Multiplicity'], sep='')

# Assigning responseveness and group annotations
targets <- read_csv(paste(data_dir, annotf, sep='/'))
targets <- as.data.frame(targets[complete.cases(targets), ])

# Removing data columns without annotation
raw_df <- raw_df[, colnames(raw_df) %in% targets$sample]

# Forcing same order of samples in snnotations and data frame
raw_df = raw_df[, targets$sample]

# Consider zeroes as NaN
raw_df[raw_df == 0] <- NA

write.csv(raw_df, paste(data_dir, 'raw.csv', sep='/'))

# Do the plots
subdir <- paste(out_dir, '1_log2', sep='/')
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

magicPlotMaker(log2(raw_df), outpath = subdir, w=15, h=15,
               targets=targets[, c('sample', 'condition')])

# Normalizing and removing batch effect
batches <- targets$batch
df_norm <- as.data.frame(matrix(NA, nrow=nrow(raw_df), ncol=ncol(raw_df)))
colnames(df_norm) <- colnames(raw_df)
rownames(df_norm) <- rownames(raw_df)

for(b in unique(batches)){
    subdf <- raw_df[, targets[targets$batch == b, 'sample']]
    subdf_norm <- as.data.frame(justvsn(as.matrix(subdf)))
    
    df_norm[rownames(subdf_norm), colnames(subdf_norm)] <- subdf_norm
}

df_bcor <- as.data.frame(removeBatchEffect(df_norm, batch=batches))

# Do the plots
subdir <- paste(out_dir, '1_norm', sep='/')
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

magicPlotMaker(df_bcor, outpath = subdir, w=15, h=15,
               targets=targets[, c('sample', 'condition')])

# Save the normalized data
write.csv(df_bcor, paste(data_dir, 'norm_data.csv', sep='/'))

