library(readr)
library(dplyr)
library(tibble)

sample_tab <- read_table('/home/qingzhang/stringtie-merged/abun/p001_HEK293T_S1_SYSY_input_SRR494613.tab')
gene_pool <- sample_tab$`Gene ID` # tested consistent in all samples


uni_tab <- function(t){
    t <- t %>% dplyr::rename(geneName = `Gene Name`) %>% dplyr::select(geneName, FPKM)
    rep.genes <- unique(t$geneName[duplicated(t$geneName)])
    t_single <- t %>% dplyr::filter(!geneName %in% rep.genes)
    t_multi <- t %>% dplyr::filter(geneName %in% rep.genes)
    tl <- lapply(split(t_multi, t_multi$geneName), function(df){
        geneName <- unique(df$geneName)
        FPKM <- sum(df$FPKM)
        tibble(geneName,FPKM)
    })
    o <- rbind(do.call(rbind, tl), t_single)
    o[order(o$geneName),]
}

d <- read_rds('srr_withgtfs.rds')
d <- d %>% dplyr::filter(grepl("_input", wecall))
home_path <- "/home/qingzhang/stringtie-merged/abun/"
file_names <- paste0(d$wecall, ".tab")


#create a data container
mc <- matrix(nrow = nrow(d), ncol = length(gene_pool), 
             dimnames = list(d$sra_acc, gene_pool))


for (i in seq_along(file_names)){
    print(i)
    file_to_read <- paste0(home_path, file_names[i])
    srr_id <- d$sra_acc[i]
    
    t0 <- read_tsv(file_to_read)
    t1 <- uni_tab(t0)
    
    mc[srr_id,] <- t1$FPKM[match(colnames(mc), t1$geneName)]
}

write_rds(mc, 'test_exp_ens.rds')


