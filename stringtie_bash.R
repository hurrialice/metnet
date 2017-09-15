library(readr)
library(dplyr)
library(tidyr)
dt <- read_tsv("MetDB2_DataSummary_v1.4.tsv")
dt <- dt %>% dplyr::filter(species =='Homo sapiens')
dt_single <- dt %>% dplyr::filter(end_type=="single") %>% dplyr::filter(!source_abbr == "p020")
dt_pair <- dt %>% dplyr::filter(end_type=="pair")

filenames <- read_rds('file-names.rds')
pair_filenames <- paste0(strsplit(dt_pair$sra_acc, split = ',') %>% unlist(), '_accept_hits.bam')
dt_pair$sra_acc <- strsplit(dt_pair$sra_acc, split = ',')
dt_pair$wecall <- mapply(function(a,b,c) paste0(a,"_",b,"_",c),
                          a = dt_pair$source_abbr,
                          b = dt_pair$group_name,
                          c = dt_pair$type)
need <- dt_pair %>% dplyr::select(sra_acc, type, frag_size, wecall)
pair_b2n<- unnest(need)
# first we use stingtie for paired end 

sink("st.sh")
for( i in seq(nrow(pair_b2n))){
    srr <- pair_b2n$sra_acc[i]
    stout <- paste0('/home/qingzhang/meth-qing/test/',pair_b2n$wecall[i],'_',pair_b2n$sra_acc[i],'.gtf')
    file <- paste0("/media/sano/easystore/metdb2_all_bams_hg19/",
                   srr,"/",srr,'_accept_hits.bam')
    
    cat('stringtie -p 8 -G hg19_genes.gtf -o', stout, "-l", srr, file)#,"&")
    cat("\n")
}
cat("echo Done!")
sink()

sink("mergelist.txt")
for( i in seq(nrow(pair_b2n))){
    srr <- pair_b2n$sra_acc[i]
    stout <- paste0('/home/qingzhang/meth-qing/stringtie/',pair_b2n$wecall[i],'_',pair_b2n$sra_acc[i],'.gtf')
    #file <- paste0("/media/sano/easystore/metdb2_all_bams_hg19/",
                  # srr,"/",srr,'_accept_hits.bam')
    
    cat(stout)#,"&")
    cat("\n")
}
sink()


sink("make_count_table.sh")
for (i in seq(nrow(pair_b2n))){
    srr <- pair_b2n$sra_acc[i]
    file <- paste0("/media/sano/easystore/metdb2_all_bams_hg19/",
                   srr,"/",srr,'_accept_hits.bam')
    stout <- paste0('/home/qingzhang/ballgown/',pair_b2n$wecall[i],'_',srr,'.gtf')
    cat("stringtie -e -B -p 8 -G stringtie_merged0.gtf -o ", stout, file, " &")
    cat("\n")
}
sink()

sink("stringtie_abun.sh")
for (i in seq(nrow(pair_b2n))){
    srr <- pair_b2n$sra_acc[i]
    stout <- paste0('/home/qingzhang/meth-qing/test/',pair_b2n$wecall[i],'_',pair_b2n$sra_acc[i],'.gtf')
    file <- paste0("/media/sano/easystore/metdb2_all_bams_hg19/",
                   srr,"/",srr,'_accept_hits.bam')
    stabun <- paste0('/home/qingzhang/meth-qing/abun/',pair_b2n$wecall[i],'_',pair_b2n$sra_acc[i],'.tab')
    cat('stringtie -p 8 -G hg19_genes.gtf -o', stout,"-A", stabun, "-l", srr, file ,"&")
    cat("\n")
}
sink()
####
### to understand what is in the gff file

library(rtracklayer)
sam <- import.gff2('sam.gtf')
