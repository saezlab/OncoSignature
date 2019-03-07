library(parallel)
library(piano)
library(ggplot2)
library(biomaRt)
library(dplyr)


run_piano <- function(methods, ttop, gset, ncores = 1){
    # Runs piano in parallel over the methods listed in 'methods'.
    # ARGUMENTS:
    #   * methods [character]: List of gene set statistic methods that
    #     will be computed. These can be: 'fisher', 'stouffer',
    #     'reporter', 'tailStrength', 'wilcoxon', 'mean', 'median',
    #     'sum', 'page', 'maxmean' and/or 'gsea'.
    #   * ttop [data.frame]: Contains the fitted data. It is assumed
    #     that p- and t-values are in columns named 'P.Value' 't'
    #     respectively (as the output of topTable function from limma
    #     package). Row names must correspond to the identifiers used on
    #     the gset.
    #   * gset [GSC]: Output from piano's function 'loadGSC' whose
    #     argument is the mapping table between the entries and the
    #     annotations.
    #   * ncores [numeric]: Number of cores used to run the different
    #     methods in parallel. Default is 1.
    
    ttop = ttop[complete.cases(ttop), ]
    
    res = mclapply(methods, function(m, ttop, gset){
        # Always use t-values except for the methods that do not support it
        if (m == 'fisher' | m == 'stouffer' | m == 'reporter' | m == 'tailStrength'){
            values = ttop$P.Value
        }
        
        else{
            values = ttop$t
        }
        print(paste('- Currently computing method:', m, sep=' '))
        names(values) = rownames(ttop)
        
        dirs = as.data.frame(ttop$logFC)
        row.names(dirs) = rownames(ttop)
        
        return(runGSA(values, gsc=gset, adjMethod='fdr', geneSetStat=m,
                      directions=dirs))
        
    }, ttop, gset, mc.cores=ncores)
    
    names(res) = methods
    
    return(res)
}


# Set up working directory
#setwd('/media/nico/JRC_COMBINE_NICO_BAK/OncoSignature/ORGANIZED/2+2_AML_patients/')

# LOADING GENESETS: -------------------------------------------------------->>>

# - MSigDB genesets (manually downloaded)
msigdb_path_to_gsymbol = read.csv('data/msigdb_path_to_gsymbol.csv',
                                  stringsAsFactors=F)
# - Our IDs
phospho_ids = as.character(read.csv('data/norm_data.csv')$X)
upids = gsub('[_].*', '', phospho_ids) # Extract UniProt IDs
# - Loading Mart for Human
mart = useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')
# - Mapping UniProtIDs to GeneSymbol
upid_to_gsymb = getBM(c('uniprotswissprot', 'external_gene_name'),
                      filters='uniprotswissprot', values=upids, mart)
colnames(upid_to_gsymb) = c('UniProtKB', 'GeneSymbol')
# - Merging UniProt IDs to genesets through GeneSymbol
map_table = left_join(upid_to_gsymb, msigdb_path_to_gsymbol, by='GeneSymbol')
# - Merging phosphosites to genesets through UniProt IDs
upids_psites = data.frame(upids, phospho_ids)
colnames(upids_psites) =  c('UniProtKB', 'PSite')
map_table = left_join(upids_psites, map_table, by='UniProtKB')
# - Cleaning the mapping table
map_table = unique(map_table[, c('PSite', 'GeneSet')])

gset = loadGSC(map_table)

# <<<--------------------------------------------------------------------------

# Which piano methods will be computed
methods = c('fisher', 'stouffer', 'reporter', 'tailStrength', 'wilcoxon',
            'mean', 'median', 'sum', 'page', 'maxmean')

# Find all table top files (contrasts) from which GSEA is to be run
all_files = dir('results/2_diff_exp/')
files = all_files[endsWith(all_files, '_ttop.csv')]

for(f in files){
    name = paste(strsplit(f, '_')[[1]][1:2], collapse='_')
    
    print(paste('Computing GSEA for contrast:', name, sep=' '))
    
    # Set-up the subdirectory
    subdir = paste('results/3_gsea', name, sep='/')
    ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)
    
    # Loading data
    df = as.data.frame(read.csv(paste('results/2_diff_exp', f, sep='/'), row.names=1))
    
    # Running piano
    piano_results = run_piano(methods, df, gset, ncores = 4)
    
    # Saving the environment
    save(list=ls(all.names=T),
         file=paste(subdir, paste(name, 'GSEA.RData', sep='_'), sep='/'),
         envir=.GlobalEnv)
    
    # Consensus heatmap
    pdf(paste(subdir, 'cons_hmap.pdf', sep='/'), width=150, height=200)
    ch = consensusHeatmap(piano_results, cutoff=50, method='median',
                          ncharLabel=50, cellnote='medianPvalue',
                          cex=0.2, plot=T)
    dev.off()
    
    # Consensus scores (distinct-directional)
    #pdf(paste(subdir, 'cs_dist_up.pdf', sep='/'))
    cs_dist_up = piano::consensusScores(piano_results,
                                        class='distinct',
                                        direction='up', plot=F)
    #dev.off()
    
    #pdf(paste(subdir, 'cs_dist_dw.pdf', sep='/'))
    cs_dist_dw = piano::consensusScores(piano_results,
                                        class='distinct',
                                        direction='down', plot=F)
    #dev.off()
    
    write.table(cs_dist_up$rankMat, paste(subdir, 'dist_up_rank.txt', sep='/'),
                sep='\t', quote=F)
    write.table(cs_dist_dw$rankMat, paste(subdir, 'dist_dw_rank.txt', sep='/'),
                sep='\t', quote=F)
    
    write.table(cs_dist_up$pMat, paste(subdir, 'dist_up_pval.txt', sep='/'),
                sep='\t', quote=F)
    write.table(cs_dist_dw$pMat, paste(subdir, 'dist_dw_pval.txt', sep='/'),
                sep='\t', quote=F)
    
    # Consensus scores (mixed-directional)
    #pdf(paste(subdir, 'cs_mix_up.pdf', sep='/'))
    cs_mix_up = piano::consensusScores(piano_results,
                                       class='mixed',
                                       direction='up', plot=F)
    #dev.off()
    
    #pdf(paste(subdir, 'cs_mix_dw.pdf', sep='/'))
    cs_mix_dw = piano::consensusScores(piano_results,
                                       class='mixed',
                                       direction='down', plot=F)
    #dev.off()
    
    write.table(cs_mix_up$rankMat, paste(subdir, 'mix_up_rank.txt', sep='/'),
                sep='\t', quote=F)
    write.table(cs_mix_dw$rankMat, paste(subdir, 'mix_dw_rank.txt', sep='/'),
                sep='\t', quote=F)
    
    write.table(cs_mix_up$pMat, paste(subdir, 'mix_up_pval.txt', sep='/'),
                sep='\t', quote=F)
    write.table(cs_mix_dw$pMat, paste(subdir, 'mix_dw_pval.txt', sep='/'),
                sep='\t', quote=F)
    
    # Consensus scores (non-directional)
    #pdf(paste(subdir, 'cs_non.pdf', sep='/'))
    cs_non = piano::consensusScores(piano_results,
                                    class='non', plot=F)
    #dev.off()
    
    write.table(cs_non$rankMat, paste(subdir, 'non_rank.txt', sep='/'),
                sep='\t', quote=F)
    
    write.table(cs_non$pMat, paste(subdir, 'non_pval.txt', sep='/'),
                sep='\t', quote=F)
}
