# 
library(readr)
library(preprocessCore)
library(dplyr)
library(igraph)

test_m <- read_rds('selectedM.rds') # made by xiangyu, do not change!!
test_m <- t(as.matrix(test_m))
write_rds(test_m, 'test_m_fin.rds') # direct input to cor!!

test_m <- read_rds('test_m_fin.rds')
bin <- read_rds('bin.rds') # made by xiangyu , nochange to this root file!
bin_rpkm <- read_rds('bin_rpkm.rds') # made by xiangyu, never change this root!
m_dict <- read_rds('msites_BP_fin.rds') 



# add gene id column to bin
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
hggene_ID <- names(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
txdbgene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
siteid <- findOverlaps(bin,txdbgene)
bin$gene_id <- NA
bin$gene_id[queryHits(siteid)] <- txdbgene$gene_id[subjectHits(siteid)]
bindf <- cbind(mcols(bin), mcols(bin_rpkm)) 
# write bindf as 'bin_comp.rds')


bindf <- read_rds('bin_comp.rds')

# filter the selected conditions
filter_conds <- function(df, standard){
    sel.conds <- rownames(standard)
    input.conds <- paste0(sel.conds, '_input')
    out <- df[,c('modstart','gene_id', 'modName',input.conds)]
    colnames(out) <- gsub('_input','', colnames(out))
    out
}
raw_input <- filter_conds(df = bindf, standard = test_m )

log2trans <- function(df){
    nums <- as.matrix(df[,-c(1,2,3)])
    nums <- log2(nums+0.01)
    out <- cbind(df[,c(1,2)], as.data.frame(nums))
    out
}
norm_input <- as.data.frame(log2trans(raw_input))
norm_withid <- norm_input[!is.na(norm_input$gene_id),]



raw_split <- split(norm_withid, norm_withid$gene_id)
del12 <- function(rl){
    ifremove <- vector(length = length(rl))
    for(i in seq_along(rl)){
        #browser()
        df <- rl[[i]]
        #browser()
        if (nrow(df) < 4){
            ifremove[i] <- TRUE
        }
        #browser()
    }
    rl[ifremove] <- NULL
    rl
}
four_plus <- del12(raw_split)

library(reshape2)
overml <- lapply(four_plus, function(df){
    m <- matrix(ncol = nrow(df), nrow = nrow(df))
    for (i in seq(nrow(m))){
        for (j in seq(ncol(m))){
            if (abs(df$modstart[i] - df$modstart[j]) <100){
                m[i,j] <- 1
            }
            else {m[i,j] <- 0}
        }
    }
    m[upper.tri(m, diag = T)] <- 0
    outdf <- melt(m)
    outdf[outdf$value == 1,]
})
library(dplyr)
overml <- lapply(overml, function(df){dplyr::rename(df, ovlp = value)})
# calculate pearson correlation for possibily merged signals:

cor_wgene <- lapply(four_plus, function(df){
    cor <- cor(t(as.matrix(df[-c(1,2)])), method = 'pearson')
    #browser()
    colnames(cor) <- as.character(seq(nrow(df)))
    rownames(cor) <- as.character(seq(nrow(df)))
    cor[upper.tri(cor, diag = T)] <- 0
    outdf <- melt(cor)
    outdf[outdf$value > 0.8,]
})
cor_wgene <- lapply(cor_wgene, function(df){dplyr::rename(df, cor = value)})
comb_bins <- Map(function(o,c){
    df <- left_join(o,c)
    df <- df[complete.cases(df),]
    df
}, c = cor_wgene, o = overml)



maxc <- lapply(comb_bins, function(df){
    edgelist <- df[,c(1,2)]
    g <- graph.data.frame(edgelist, directed = F)
    mc <- lapply(max_cliques(g), as_ids)
    lapply(mc, as.integer)
})
#
trouble_list <- lapply(maxc, function(l){
    as.integer(unique(unlist(l)))
}) # nrow identifier of test_exp dataframe
g <- graph.data.frame(comb_bins$`1000`[,c(1,2)])

safe_test <- Map(function(a,b){
    a[-b,]},
    a = four_plus,
    b = trouble_list)
trouble_test <- Map(function(a,b){
    a[b,]},
    a = four_plus,
    b = trouble_list)

ids <- Map(function(a,b){nrow(a)+nrow(b)},
    a = safe_test,
    b = trouble_test)
ids_true <- lapply(four_plus, nrow)
stopifnot(identical(unlist(ids), unlist(ids))) # checked!
rm(trouble_test)



fp0 <- lapply(four_plus, function(df){df[,-c(1,2)]})
fp <- fp0 %>% lapply(function (x) { rownames(x) <- NULL; x })
ml0 <- ml
fp <- lapply(fp0, unname)
ml <- lapply(ml, function(x) as.list(x$value))

fp_merge <- mapply(ml = ml, fp = fp,
    function (fp, ml) {
        li.subdf <- lapply(ml, function (whichrow) fp[whichrow, ])
        li.df <- lapply(li.subdf, colMeans)
        ans <- do.call(rbind, li.df)
        ans
    }
)

fp_merge <- do.call(DataFrameList, fp_merge)
fp_merge <- fp_merge %>% as.list

stopifnot(unique(lengths(ml) == sapply(fp_merge, nrow))) # checked!
fp_merge1 <- lapply(fp_merge,as.data.frame)
merged_geneids <- Map(function(t,s){rbind(t,s)},
                      t = fp_merge1, 
                      s = lapply(safe_test, function(df){
                          out = df[,-c(1,2)];
                          rownames(out) <- NULL}) )
#
NAfind <- unlist(lapply(merged_geneids, function(df){any(is.na(df))}))

merged_colM <- lapply(merged_geneids, function(df){apply(df, 2, mean)})
table(unlist(lapply(merged_colM, function(df){any(is.na(df))})))

exps <- do.call(rbind, merged_colM)

remain23 <- function(rl){
    ifremove <- vector(length = length(rl))
    for(i in seq_along(rl)){
        #browser()
        df <- rl[[i]]
        #browser()
        if (nrow(df) ==1 | nrow(df) >= 4){
            ifremove[i] <- TRUE
        }
        #browser()
    }
    rl[ifremove] <- NULL
    rl
}
twonthree <- remain23(raw_split)

twonthree <- lapply(twonthree, function(df){df[,-c(1,2)]})
exps_simple <- do.call(rbind,lapply(twonthree, function(df){
    apply(df, 2, mean)
}))
##### mer_exp
mer_exp <- rbind(exps, exps_simple)

exp_sd <- apply(mer_exp, 1, sd)
sdcut <- quantile(exp_sd, 0.25)
id_remain <- names(exp_sd[exp_sd > sdcut])

chose_exp <- mer_exp[id_remain,]


feature_scaling <- function(m){
    sd <- sd(m)
    mean <- mean(m)
    out <- apply(m, c(1,2), function(e){
        (e - mean)/sd
    })
    out
}

fin_exp <- feature_scaling(chose_exp)
test_exp <- t(fin_exp)
write_rds(test_exp, 'test-exp-0824.rds')












###################################

library(gplots)
i <-1
# try with the first gene id
while(TRUE) {
    test1 <- merged_geneids[[i]][, -c(1,2)]
    test1 <- as.matrix(test1)
    heatmap.2(test1, Rowv = F)
    i <- sample(1:length(merged_geneids), 1)
    Sys.sleep(2) }
# how much of this outliers are due to overlapping of signals from other transcripts?
gene_gr <- txdbgene[txdbgene$gene_id %in% names(raw_split),]
siteid <- findOverlaps(gene_gr, gene_gr)
all_over <- as.data.frame(siteid)
remain <- mapply(function(a,b) !a == b, a = all_over$queryHits, b = all_over$subjectHits)
ovlp <-all_over[remain,]
ovlp_ids <- apply(ovlp, c(1,2), function(hits){
    gene_gr$gene_id[match(hits,seq(gene_gr))]
})
gover <- graph.data.frame(ovlp_ids, loops = F)
deg_over <- degree(gover, v = V(gover), loops = F)
deg_over <- sort(deg_over, decreasing = T)
# overall the impact might be minimal. only less than 10% is affected.



# one of the only viable way to reduce size according to what we have:
# choose the variance in top 75%
cal.var <- function(rl){
    out <- list()
    for (i in seq_along(rl)){
        df <- rl[[i]]
        df <- df[,-c(1,2)]
        out[[i]] <- colMeans(df)
    }
    do.call(rbind, out)
}
exp_mean <- cal.var(raw_split)
rownames(exp_mean) <- names(raw_split)

# calculate varience
exp_sd <- apply(exp_mean, 1, sd)
sdcut <- quantile(exp_sd, 0.25)
id_remain <- names(exp_sd[exp_sd>sdcut])

exp_mean_cut <- exp_mean[rownames(exp_mean) %in% id_remain,]
test_exp0 <- t(exp_mean_cut)
write_rds(test_exp0, 'sdcut_testexp.rds')


# make both test_exp and test_m under the same range # feature scaling
feature_scaling <- function(m){
    sd <- sd(m)
    mean <- mean(m)
    out <- apply(m, c(1,2), function(e){
        (e - mean)/sd
    })
    out
}
test_ms <- feature_scaling(test_m)
test_exps <- feature_scaling(test_exp0)


# make adjacency matrix
# this part we need a server!!
my_testm <- test_ms[,m_dict$modName]
write_rds(my_testm, 'test-m0823fs.rds')
write_rds(test_exps, 'test-exp0823fs.rds')

adjmake0 <- function(x, y, quant, pcut){
    cor <- cor(x,y,method = 'spearman')
    rawp <- corPvalueFisher(cor, 40)
    mt <- mt.rawp2adjp(rawp, proc = 'Bonferroni')
    adj <- mt$adjp[,2]
    adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
    cor <- abs(cor)
    diag(cor) <- 0 # here is the difference
    scc_cut <- quantile(cor, quant)
    cor[adjp > pcut | cor < scc_cut] <- 0
    cor[cor > 0] <- 1
    cor
}

adjmake <- function(x, y, quant, pcut){
    cor <- cor(x,y,method = 'spearman')
    rawp <- corPvalueFisher(cor, 40)
    mt <- mt.rawp2adjp(rawp, proc = 'Bonferroni')
    adj <- mt$adjp[,2]
    adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
    cor <- abs(cor)
    scc_cut <- quantile(cor, quant)
    cor[adjp > pcut | cor < scc_cut] <- 0
    cor[cor > 0] <- 1
    cor
}
adj_exp <- adjmake0(x = test_exp, y = test_exp, quant = 0.98, pcut = 0.01)
adj_me <- adjmake(x = test_m, y = test_exp, quant = 0.80, pcut = 0.05)







