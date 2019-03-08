library(biomaRt)
library(viper)

library(omicToolsTest)

# Set up working directory
#setwd('/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/4x_AML_cell_lines/')
subdir = 'results/3_gsea/'
ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

data = read.csv('data/norm_data.csv', stringsAsFactors=F)

# Loading mapping table for UniProtKB to GeneSymbol
upids = gsub('[_].*', '', data$X)

mart = useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')

upid_to_gsymb = getBM(c('uniprotswissprot', 'external_gene_name'),
                      filters='uniprotswissprot', values=upids, mart)
colnames(upid_to_gsymb) = c('UniProtKB', 'GeneSymbol')

all_files = dir('results/2_diff_exp/')
files = all_files[endsWith(all_files, '.csv')]

ttops = list()

for(f in files){
    fname = paste(strsplit(f, '_')[[1]][1:2], collapse='_')
    
    print(paste('Loading contrast:', fname, sep=' '))
    
    # Set-up the subdirectory
    subdir = paste('results/3_gsea', fname, sep='/')
    ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)
    
    # Loading data
    data = as.data.frame(read.csv(paste('results/2_diff_exp', f, sep='/'),
                                  stringsAsFactors=F))

    # Renaming UniProt IDs with GeneSymbol
    for(i in 1:nrow(data)){
        name_site = strsplit(data[i, 1], '_')
        name = name_site[[1]][1]
        site = name_site[[1]][2]
        map_name = upid_to_gsymb[upid_to_gsymb$UniProtKB == name, 2]
        
        if(length(map_name)>1){
            map_name = map_name[1]
            data[i, 1] = paste(map_name, site, sep='_')
        }
        
        else if(length(map_name)==0){
            data[i, 1] = NaN
        }
        
        else{
            data[i, 1] = paste(map_name, site, sep='_')
        }
    }
    
    df = data[data$X != NaN, ]
    
    ttops[[fname]] = df
}

results_raw = runViper(ttops, regulon=pathway_regulon, nperm=10000,
                       mode='signature', minsize=25)
names(results_raw) = names(ttops)

results = makeViperResDf(results_raw)
write.csv(results, 'results/3_gsea/results_gsea_viper.csv', row.names=F)
