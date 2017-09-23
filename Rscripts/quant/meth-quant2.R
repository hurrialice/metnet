library(rtracklayer)
library(readr)
library(dplyr)
library(tibble)
library(GenomicFeatures)
library(GenomicRanges)

rm <- read_table2('RMBase_hg19_all_m6A_site.txt')
rm0 <- rm %>% filter(!is.na(score2)) %>% filter(supportNum > 10) %>%
    dplyr::select(chromosome, modStart, modEnd, modName, strand, 
                  supportNum, pubmedIds, geneName, geneType) %>% 
    dplyr::rename(transcript_id = modName, gene_id = geneName, 
           start = modStart, end = modEnd, seqname = chromosome)
rm0$type <- 'exon'
rmdf <- DataFrame(rm0)
rmgr <- makeGRangesFromDataFrame(rmdf, keep.extra.columns = TRUE)
rm(rm, rm0)
write_rds(rmgr, 'mgr18w.rds')
# read FPKM values.
d <- read_rds('srr_withgtfs.rds')
msites <- read_rds('msites.rds')
home_path <- "/home/qingzhang/meth-qing/stringtie-meths/"

gtf2df <- function(file_to_read){
    gr <- rtracklayer::import(file_to_read, 'gtf')
    df <- as.tibble(mcols(gr)) %>% dplyr::filter(type == 'transcript') %>% 
        dplyr::select( transcript_id, FPKM) %>% dplyr::rename(modName = transcript_id)
}


# initialize with a container
mc <- matrix(nrow = nrow(d), ncol = length(msites), 
                   dimnames = list(d$sra_acc, msites))
for (i in seq(nrow(d))){
    print(i)
    srr_id <- d$sra_acc[i]
    file_to_read <- paste0(home_path, d$wecall[i], '.gtf')
    df <- gtf2df(file_to_read)
    print(df)
    mc[srr_id,] <- df$FPKM[match(colnames(mc), df$modName)]
}

write_rds(mc, 'test_m_raw.rds')

