library(readr)
library(dplyr)
library(vsn)

library(omicToolsTest)

# Set up CWD and load raw data
#setwd('/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/4x_AML_cell_lines/')

raw_data = read_delim('data/raw_data.txt', '\t',
					  escape_double=FALSE, trim_ws=TRUE)

# Locate the intensity within the raw data file
raw_df = as.data.frame(raw_data[, 1:32])
raw_df$ID = paste(paste(raw_data$Protein, raw_data$`Amino acid`, sep='_'),
                  raw_data$Position, sep='')

# Simplifying column names
colnames(raw_df) = gsub('Intensity ', '', colnames(raw_df))

# Defining sample conditions
targets = as.data.frame(matrix(NA, nrow=ncol(raw_df) - 1, ncol=2))
colnames(targets) = c('sample', 'condition')
targets[, 'sample'] = colnames(raw_df)[1:ncol(raw_df) - 1]
targets[, 'condition'] = as.character(lapply(colnames(raw_df)[1:ncol(raw_df) - 1],
                                             function(i) gsub('-. ', '_', i)))

write.csv(targets, 'data/targets.csv', row.names=F)

# Handling duplicates (sum the raw intensities)
dups_pos = duplicated(raw_df$ID)
dups_id = raw_df[dups_pos, 'ID']
unique_dups = unique(dups_id)

df = raw_df[!(raw_df$ID %in% unique_dups), -ncol(raw_df)]
rownames(df) = raw_df[!(raw_df$ID %in% unique_dups), 'ID']

for(i in unique_dups){
    aux = data.frame(colSums(raw_df[raw_df$ID == i, 1:ncol(df)]))
    colnames(aux) = i
    df = rbind(df, t(aux))
}

# Consider zeroes as NaN
df[df==0] = NA

# Do the plots
subdir = 'results/1_raw/'
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

magicPlotMaker(df, outpath = subdir, targets=targets)

# Normalize the data
df_norm = as.data.frame(justvsn(as.matrix(df)))

# Do the plots
subdir = 'results/1_norm/'
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

magicPlotMaker(df_norm, outpath = subdir, targets=targets)

# Save the normalized data
write.csv(df_norm, 'data/norm_data.csv')
