library(readr)
library(dplyr)
library(vsn)
library(limma)

library(omicToolsTest)

# Set up CWD and load raw data
setwd('/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/validation/')

raw_data = as.data.frame(
            read_delim('data/Phospho(STY)sites_withPtTitles_20181025.txt', '\t',
                       na='NaN'))

# Locate the intensity within the raw data file
sample_cols = colnames(raw_data)[startsWith(colnames(raw_data), ' V')
                               & !endsWith(colnames(raw_data), 'NormC')]
raw_df = as.data.frame(raw_data[, sample_cols])
raw_df$ID = paste0(paste(raw_data[, 142], raw_data[, 118], sep='_'), raw_data[, 137])

# Simplifying column names
colnames(raw_df) = gsub(' V._', '', colnames(raw_df))

# Defining sample conditions
targets = as.data.frame(matrix(NA, nrow=ncol(raw_df) - 1, ncol=2))
colnames(targets) = c('sample', 'condition')
targets[, 'sample'] = colnames(raw_df)[1:ncol(raw_df) - 1]
targets[, 'condition'] = as.character(lapply(colnames(raw_df)[1:ncol(raw_df) - 1],
                                             function(i){substr(i, nchar(i), nchar(i))}))

write.csv(targets, 'data/targets.csv', row.names=F)

# Handling duplicates (sum the raw intensities)
dups_pos = duplicated(raw_df$ID)
dups_id = raw_df[dups_pos, 'ID']
unique_dups = unique(dups_id)

df = raw_df[!(raw_df$ID %in% unique_dups), -ncol(raw_df)]
rownames(df) = raw_df[!(raw_df$ID %in% unique_dups), 'ID']

for(i in unique_dups){
    aux = data.frame(colSums(raw_df[raw_df$ID == i, 1:ncol(df)], na.rm=T))
    colnames(aux) = i
    df = rbind(df, t(aux))
}

# Consider zeroes as NaN
df[df==0] = NA

# Do the plots
subdir = 'results/1_log2/'
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

magicPlotMaker(log2(df), outpath=subdir, targets=targets, w=15, h=15)

# Normalizing and removing batch effect
batches = gsub('_.*', '', gsub(' ', '',  sample_cols))
df_norm = as.data.frame(matrix(NA, nrow=nrow(df), ncol=ncol(df)))
colnames(df_norm) = colnames(df)
rownames(df_norm) = rownames(df)

for(b in unique(batches)){
    subdf = df[, batches == b]
    subdf_norm = as.data.frame(justvsn(as.matrix(subdf)))

    df_norm[rownames(subdf_norm), colnames(subdf_norm)] = subdf_norm
}

df_bcor = as.data.frame(removeBatchEffect(df_norm, batch=batches))

# Do the plots
subdir = 'results/1_norm/'
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

magicPlotMaker(df_bcor, outpath = subdir, targets=targets, w=15, h=15)

# Save the normalized data
write.csv(df_bcor, 'data/norm_data.csv')
