library(readr)
library(dplyr)
library(rtracklayer)


### shell ##########
# shell script for m6A
d <- read_rds('srr_withgtfs.rds')
sink("stringtie-meth.sh")
for( i in seq(nrow(d))){
    srr <- d$sra_acc[i]
    stout <- paste0('/home/qingzhang/meth-qing/stringtie-meths/',d$wecall[i],'.gtf')
    file <- paste0("/media/sano/easystore/metdb2_all_bams_hg19/",
                   srr,"/",srr,'_accept_hits.bam')
    guide <- d$guide_gtf[i]
    cat('stringtie -p 8 -G',guide, '-o', stout, "-l", srr,"-e", file,"&")
    cat("\n")
}
cat("echo Done!")
sink()

# shell for expression pipes
# /home/qingzhang/meth-qing/stringtie-exps-ref_only/abun
# only "_input" is used
input <- d %>% dplyr::filter(grepl('_input', d$wecall))
sink("stringtie-exps.sh")
for( i in seq(nrow(input))){
    srr <- input$sra_acc[i]
    stout <- paste0('/home/qingzhang/meth-qing/stringtie-exps/',input$wecall[i],'.gtf')
    file <- paste0("/media/sano/easystore/metdb2_all_bams_hg19/",
                   srr,"/",srr,'_accept_hits.bam')
    #abun <- paste0('/home/qingzhang/meth-qing/stringtie-exp-ref_only/abun/',srr_dict$wecall[i],".tab")
    cat('stringtie -p 8 -G hg19_genes.gtf -o', stout, "-l", srr, "-e",  file,"&")
    cat("\n")
}
cat("echo Done!")
sink()

# make mergelist
sink("mergelist0918.txt")
for( i in seq(nrow(input))){
    stout <- paste0('/home/qingzhang/meth-qing/stringtie-exps/',input$wecall[i],'.gtf')
    cat(stout)
    cat("\n")
}
sink()

# stringtie --merge
cat('stringtie --merge -p 8 -G hg19_genes.gtf -o stringtie-exps/stringtie_merged_exps.gtf mergelist0918.txt')



# redo stringtie with merged gtf
# give the abundance table!
sink("stringtie-exps-m1.sh")
for( i in seq(nrow(input))){
    srr <- input$sra_acc[i]
    stout <- paste0('/home/qingzhang/meth-qing/stringtie-exps-merged1/',input$wecall[i],'.gtf')
    file <- paste0("/media/sano/easystore/metdb2_all_bams_hg19/",
                   srr,"/",srr,'_accept_hits.bam')
    guide <- "/home/qingzhang/meth-qing/stringtie-exps/stringtie_merged_exps.gtf"
    abun <- paste0('/home/qingzhang/meth-qing/stringtie-exp-merged1/abun/',input$wecall[i],".tab")
    cat('stringtie -p 8 -G',guide, '-o', stout, "-l", srr, "-A",abun, "-e",  file,"&")
    cat("\n")
}
cat("echo Done!")
sink()
# to get the output warnings: 
cat('bash stringtie-exps-m1.sh > m1sh.out')



