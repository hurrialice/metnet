# make gtf for subsequent stringtie m6A quantification

library(rtracklayer)
library(GenomicRanges)
library(readr)
library(dplyr)

# make raw rmgr
rm <- read_table2('RMBase_hg19_all_m6A_site.txt')
rm <- rm[-1,]
rm0 <- rm %>% dplyr::select(chromosome, modStart, modEnd, modName, geneName,
                            strand, supportNum )
rm0 <- rename(rm0, transcript_id = modName, gene_id = geneName, 
              start = modStart, end = modEnd, seqname = chromosome)
rm0$type <- 'exon'
rmgr <- makeGRangesFromDataFrame(rm0, keep.extra.columns = TRUE)


# make gtfs
d <- read_rds('srr_dict.rds')
bin_size_pool <- d$bin_size %>% unique
for (i in seq_along(bin_size_pool)){
    bin <- bin_size_pool[i]
    fn <- paste0('RMBase_bin',bin,'.gtf')
    new_gr <- resize(rmgr, bin, fix = 'center')
    rtracklayer::export(new_gr, fn, 'gtf')
}
d$guide_gtf <- paste0('RMBase_bin',d$bin_size,'.gtf')
write_rds(d, 'srr_withgtfs.rds')
