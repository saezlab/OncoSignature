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
# Gene Set Enrichment Analysis
# ============================
#
# This script computes the GSEA with 11 different methods and generates a
# consensus score based on ranking

library(snow)
library(piano)
library(ggplot2)
library(biomaRt)
library(dplyr)

#----------------------------------- INPUT -----------------------------------#
res_dir <- 'results_bak/perseus'
out_dir <- 'results_bak/3_gsea'

data_dir <- 'data'
# If exists, reads it, otherwise writes it
mapping_file <- 'psite_to_gset.csv'
# Ignored if above exists
msigdb_file <- 'msigdb_c2_all_v7_2_GS.csv'

# Which piano methods will be computed
methods <- c('fisher',
             'stouffer',
             'reporter',
             'tailStrength',
             'wilcoxon',
             'mean',
             'median',
             'sum')#,
             #'page',
             #'maxmean')
#-----------------------------------------------------------------------------#

#' Runs GSEA in parallel over the requested methods
#'
#' @param methods [character] List of gene set statistic methods that will be
#' computed. These can be: 'fisher', 'stouffer', 'reporter', 'tailStrength',
#' 'wilcoxon', 'mean', 'median', 'sum', 'page', 'maxmean' and/or 'gsea'
#' @param [data.frame] Contains the fitted data. It is assumed that p- and
#' t-values are in columns named 'P.Value' 't' respectively (as the output of
#' topTable function from limma package). Row names must correspond to the
#' identifiers used on the gset
#' @param [GSC] Output from piano's function 'loadGSC' whose argument is the
#' mapping table between the entries and the annotations
#' @param [numeric] Number of cores used to run the different methods in
#' parallel. Default is 1
#'
#' @returns [list] Named list containing the method names and their
#' corresponding results
run_piano <- function(methods, ttop, gset, ncores=1){
    # Remove NaNs
    ttop <- ttop[complete.cases(ttop), ]

    res <- lapply(methods, function(m, ttop, gset){
        # Always use t-values except for the methods that do not support it
        #if (m == 'fisher'
        #    | m == 'stouffer'
        #    | m == 'reporter'
        #    | m == 'tailStrength'){
        values <- 10^(-as.numeric(ttop[, 'X...log10.P.']))
        #}
        #else{
        #    values <- ttop$t
        #}

        names(values) <- rownames(ttop)

        dirs <- as.data.frame(ttop$Difference)
        row.names(dirs) <- rownames(ttop)

        return(runGSA(values, gsc=gset, adjMethod='fdr', geneSetStat=m,
                      directions=dirs, ncpus = ncores))

    }, ttop, gset)

    names(res) <- methods

    return(res)
}

# LOADING GENESETS: -------------------------------------------------------->>>
map_table_path <- paste(data_dir, mapping_file, sep='/')

if (file.exists(map_table_path)){
    map_table <- read.csv(map_table_path, stringsAsFactors=F, row.names = 1)
} else {
    # - MSigDB genesets (manually downloaded)
    msigdb_path_to_gsymbol <- read.csv(paste(data_dir, msigdb_file, sep='/'),
                                       stringsAsFactors=F)
    colnames(msigdb_path_to_gsymbol) <- c('GeneSet', 'GeneSymbol')

    # XXX: OLD, doesn't work
    # - Our IDs
    phospho_ids_ex <- as.character(read.csv(paste(data_dir, 'norm_data_ex.csv',
                                               sep='/'))$X)
    phospho_ids_cl <- as.character(read.csv(paste(data_dir, 'norm_data_cl.csv',
                                               sep='/'))$X)
    phospho_ids <- unique(c(phospho_ids_cl, phospho_ids_ex))
    
    upids <- gsub('[_].*', '', phospho_ids) # Extract UniProt IDs
    # - Loading Mart for Human
    mart <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')
    # - Mapping UniProtIDs to GeneSymbol
    upid_to_gsymb <- getBM(c('uniprotswissprot', 'external_gene_name'),
                           filters='uniprotswissprot', values=upids, mart)
    colnames(upid_to_gsymb) <- c('UniProtKB', 'GeneSymbol')
    # - Merging UniProt IDs to genesets through GeneSymbol
    map_table <- left_join(upid_to_gsymb, msigdb_path_to_gsymbol, by='GeneSymbol')
    # - Merging phosphosites to genesets through UniProt IDs
    upids_psites <- data.frame(upids, phospho_ids)
    colnames(upids_psites) <- c('UniProtKB', 'PSite')
    map_table <- left_join(upids_psites, map_table, by='UniProtKB')
    # - Cleaning the mapping table
    map_table <- unique(map_table[, c('PSite', 'GeneSet')])

    write.csv(map_table, map_table_path)
}

gset <- loadGSC(map_table)
# <<<--------------------------------------------------------------------------
dex_dirs <- list.dirs(res_dir)[2:3]

for (dex_dir in dex_dirs){
    # Find all table top files (contrasts) from which GSEA is to be run
    all_files <- dir(dex_dir)
    files <- all_files[startsWith(all_files, 'DEX')]
    sdir <- as.character(strsplit(dex_dir, '/')[[1]][3])

    for(f in files){
        name <- as.character(strsplit(f, '_')[[1]][2])
        name <- gsub('.txt', '', name)

        print(paste('Computing GSEA for contrast:', name, sep=' '))

        # Set-up the subdirectory
        subdir <- paste(out_dir, sdir, name, sep='/')
        ifelse(!dir.exists(subdir), dir.create(subdir, recursive=T), F)

        # Loading data
        df <- as.data.frame(read.csv(paste(dex_dir, f, sep='/'), sep='\t'))
        rownames(df) <- df[, 'uniprot']
        df <- df[c('Difference', 'X...log10.P.')]
        
        # Running piano
        piano_results <- run_piano(methods, df, gset, ncores = 8)

        # Saving the environment
        save(list=ls(all.names=T),
             file=paste(subdir, paste(name, 'GSEA.RData', sep='_'), sep='/'),
             envir=.GlobalEnv)

        # Consensus heatmap
        pdf(paste(subdir, 'cons_hmap.pdf', sep='/'), width=150, height=200)
        ch <- consensusHeatmap(piano_results, cutoff=50, method='median',
                              ncharLabel=50, cellnote='medianPvalue',
                              cex=0.2, plot=T)
        dev.off()

        # Consensus scores (distinct-directional)
        cs_dist_up <- piano::consensusScores(piano_results,
                                            class='distinct',
                                            direction='up', plot=F)
        cs_dist_dw <- piano::consensusScores(piano_results,
                                            class='distinct',
                                            direction='down', plot=F)

        write.table(cs_dist_up$rankMat, paste(subdir, 'dist_up_rank.txt', sep='/'),
                    sep='\t', quote=F)
        write.table(cs_dist_dw$rankMat, paste(subdir, 'dist_dw_rank.txt', sep='/'),
                    sep='\t', quote=F)
        write.table(cs_dist_up$pMat, paste(subdir, 'dist_up_pval.txt', sep='/'),
                    sep='\t', quote=F)
        write.table(cs_dist_dw$pMat, paste(subdir, 'dist_dw_pval.txt', sep='/'),
                    sep='\t', quote=F)

        # Consensus scores (mixed-directional)
        cs_mix_up <- piano::consensusScores(piano_results,
                                           class='mixed',
                                           direction='up', plot=F)
        cs_mix_dw <- piano::consensusScores(piano_results, class='mixed',
                                            direction='down', plot=F)

        write.table(cs_mix_up$rankMat, paste(subdir, 'mix_up_rank.txt', sep='/'),
                    sep='\t', quote=F)
        write.table(cs_mix_dw$rankMat, paste(subdir, 'mix_dw_rank.txt', sep='/'),
                    sep='\t', quote=F)

        write.table(cs_mix_up$pMat, paste(subdir, 'mix_up_pval.txt', sep='/'),
                    sep='\t', quote=F)
        write.table(cs_mix_dw$pMat, paste(subdir, 'mix_dw_pval.txt', sep='/'),
                    sep='\t', quote=F)

        # Consensus scores (non-directional)
        cs_non <- piano::consensusScores(piano_results,
                                        class='non', plot=F)

        write.table(cs_non$rankMat, paste(subdir, 'non_rank.txt', sep='/'),
                    sep='\t', quote=F)

        write.table(cs_non$pMat, paste(subdir, 'non_pval.txt', sep='/'), sep='\t',
                    quote=F)
    }
}
