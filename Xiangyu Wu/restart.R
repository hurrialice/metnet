library(readr)
library(tibble)
library(topGO)
library(GSEABase)
library(dplyr)
library(igraph)
library(WGCNA)
library(multtest)
library(gProfileR)
library(AnnotationDbi)
library(future)
library(clusterProfiler)
library(preprocessCore)
library(AnnotationDbi)
library(tidyr)
library(GenomicRanges)

# load the m6A site table & transform to GRange format
rm <- read_table2('RMBase_hg19_all_m6A_site.txt')
rm0 <- rm %>% filter(!is.na(score2)) %>% filter(supportNum > 10) %>%
  dplyr::select(chromosome, modStart, modEnd, modName, strand, supportNum, pubmedIds, geneName, geneType) %>% 
  dplyr::rename(transcript_id = modName, gene_id = geneName, start = modStart, end = modEnd, seqname = chromosome)
rm0$type <- 'exon'
rmdf <- DataFrame(rm0)
rmgr <- makeGRangesFromDataFrame(rmdf, keep.extra.columns = TRUE)
write_rds(rmgr, 'mgr18w.rds')

# arrange the samples in manual
test_m <- read_rds('test_m_raw.rds')
homo <- read_tsv('sample_info.tsv')
h <- homo %>% dplyr::select(source_abbr, sra_acc, wecall, type, gse_acc, gsm_acc, group_name)

# separate IP from input
abbr <- h$wecall[match(rownames(test_m), h$sra_acc)]
ifip <- grepl("_ip_", abbr)
test_m_ip <- test_m[ifip,]
test_m_input <- test_m[!(ifip),]

# remove the brain sample in low quality
input_names <- grep('_input', abbr, value = T)[-11]
ip_names <- grep('_ip', abbr, value = T)[-12]
test_m_input <- test_m_input[-which(rownames(test_m_input) == 'SRR456937'),]
test_m_ip <- test_m_ip[-which(rownames(test_m_ip) == 'SRR456936'),]
input_rep_srrs <- list(
  p001 = c('SRR494613', 'SRR494615'),
  p002 = c('SRR456555', 'SRR456556', 'SRR456557'),
  p004_1 = c('SRR903368', 'SRR903369', 'SRR903370'),
  p004_2 = c('SRR903371', 'SRR903372', 'SRR903373'),
  p007_1 = c('SRR847358', 'SRR847359'),
  p007_2 = c('SRR847362', 'SRR847363'),
  p007_3 = c('SRR847366', 'SRR847367'),
  p007_4 = c('SRR847370', 'SRR847371'),
  p007_5 = c('SRR847374', 'SRR847375'),
  p009_1 = c('SRR1182585', 'SRR1182586'),
  p009_2 = c('SRR1182589', 'SRR1182590'),
  p009_3 = c('SRR1182604', 'SRR1182606'),
  p009_4 = c('SRR1182608', 'SRR1182612'),
  p009_5 = c('SRR1182610', 'SRR1182614'),
  p009_6 = c('SRR1182616', 'SRR1182618'),
  p009_7 = c('SRR1182622', 'SRR1182624'),
  p010_1 = c('SRR1035213', 'SRR1035221'),
  p010_2 = c('SRR1035215', 'SRR1035223')
)

ip_rep_srrs <- list(
  p001 = c('SRR494614', 'SRR494616'),
  p002 = c('SRR456551', 'SRR456552', 'SRR456553','SRR456554'),
  p004_1 = c('SRR903374', 'SRR903375', 'SRR903376'),
  p004_2 = c('SRR903377', 'SRR903378', 'SRR903379'),
  p007_1 = c('SRR847360', 'SRR847361'),
  p007_2 = c('SRR847364', 'SRR847365'),
  p007_3 = c('SRR847368', 'SRR847369'),
  p007_4 = c('SRR847372', 'SRR847373'),
  p007_5 = c('SRR847376', 'SRR847377'),
  p009_1 = c('SRR1182582', 'SRR1182583', 'SRR1182584'),
  p009_2 = c('SRR1182587', 'SRR1182588'),
  p009_3 = c('SRR1182603', 'SRR1182605'),
  p009_4 = c('SRR1182607', 'SRR1182611'),
  p009_5 = c('SRR1182609', 'SRR1182613'),
  p009_6 = c('SRR1182615', 'SRR1182617'),
  p009_7 = c('SRR1182621', 'SRR1182623'),
  p010_1 = c('SRR1035214', 'SRR1035222'),
  p010_2 = c('SRR1035216', 'SRR1035224')
)

# name cols and rows of matrix
all_rep_srrs <- unname(c(unlist(input_rep_srrs),unlist(ip_rep_srrs) ))
brain_srrs <- c('SRR456936','SRR456937')
valid_srrs <- h$sra_acc[which(!h$sra_acc %in% brain_srrs)]
single_srrs <- valid_srrs[!valid_srrs %in% all_rep_srrs]
single_conds_name <- paste(h$source_abbr,h$group_name, h$type, sep = '_')[match(single_srrs, h$sra_acc)]
rep_conds_name <- sapply(c(input_rep_srrs, ip_rep_srrs), function(srrs){
  unique(paste(h$source_abbr,h$group_name, h$type, sep = '_')[match(srrs, h$sra_acc)])[1]
}) %>% unname
all_conds <- c(single_conds_name, rep_conds_name)

# merge replicated samples
del.rep <- function(l, orim){
  rep.srrs <- unlist(l)
  all.srrs <- rownames(orim)
  single_srrs <- all.srrs[!(all.srrs %in% rep.srrs)]
  m_single <- orim[single_srrs,]
  rownames(m_single) <- paste(h$source_abbr,h$group_name, sep = '_')[match(rownames(m_single),h$sra_acc)]
  m_rep <- matrix(nrow = length(l), ncol = ncol(orim))
  colnames(m_rep) <- colnames(orim)
  rownames(m_rep) <- sapply(l, function(srrs){
    unique(paste(h$source_abbr,h$group_name, sep = '_')[match(srrs, h$sra_acc)])[1]
  }) %>% unname
  for (i in seq_along(l)){
    srrs_to_merge <- l[[i]]
    m_rep[i,] <- colMeans(orim[srrs_to_merge,])
  }
  rbind(m_single, m_rep)
}

# operations
unim_input <- del.rep(orim = test_m_input,l = input_rep_srrs)
unim_ip <- del.rep(orim = test_m_ip,l = ip_rep_srrs)
stopifnot(identical(rownames(unim_input), rownames(unim_ip)))
stopifnot(identical(colnames(unim_input), colnames(unim_ip)))

# calculate M value ------ methylation level on each methylation site
uni_m <- log2(unim_ip + 0.01) - log2(unim_input + 0.01)
uni_m <- uni_m[match(sort(rownames(uni_m)), rownames(uni_m)),]
write_rds(uni_m, 'test_m.rds')

# heatmap
library(RColorBrewer)
library(gplots)
hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(40)
cols <- palette(brewer.pal(8, "Dark2"))[as.numeric(rownames(uni_m))]
correlation.matrix <- matrix(NA, nrow = nrow(uni_m), ncol=nrow(uni_m))
for (i in 1:nrow(uni_m)){
  for (j in 1:nrow(uni_m)){
    correlation.matrix[i,j]<- cor(uni_m[i,],uni_m[j,])
  }    
}
par(mar=c(7,4,4,7)+0.1)
heatmap.2(correlation.matrix, labCol=rownames(uni_m), trace="none", ColSideColors=cols, srtCol = 20, cexCol = 0.68, col=hmcol, key = FALSE)

# boxplot
library(tidyr)
library(ggplot2)
unim <- gather(as.data.frame(uni_m), sample, M)
staplot <- ggplot(unim, aes(sample, M))
staplot + theme(axis.text.x = element_text(angle = 45)) + geom_boxplot()

# site filtering with mean and SD
sd_and_m_filter <- function(m, qm, qsd){
  site_ms <- apply(m, 2, mean)
  site_sds <- apply(m, 2, sd)
  
  m_cut <- quantile(site_ms, 1 - qm)
  sd_cut <- quantile(site_sds, 1 - qsd)
  
  m_retain <- (site_ms > m_cut)
  sd_retain <- (site_sds > sd_cut)
  
  retain <- mapply(function(a,b) a & b, a = m_retain, b = sd_retain)
  m[,retain]
}

filtered_M <- sd_and_m_filter(uni_m, 0.3, 0.3)

# quantile normalization
m <- preprocessCore::normalize.quantiles(t(filtered_M)) %>% t
colnames(m) <- colnames(filtered_M)
rownames(m) <- rownames(filtered_M)
write_rds(m,'sdm_filtered_testm.rds')

rmgr <- read_rds("mgr18w.rds")
filtered_gr <- rmgr[rmgr$transcript_id %in% colnames(m)]

# add ENID column on the GRange
library(EnsDb.Hsapiens.v75)
g <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(g) <- 'UCSC'
o <- findOverlaps(filtered_gr, g)
filtered_gr$ENID <- NA

dt <- as.data.frame(o) %>% tbl_df()
hits2eid <- tibble(hit = seq(length(g)), eid = g$gene_id)
dt$ensid <- hits2eid$eid[match(dt$subjectHits, hits2eid$hit)]

l <- split(dt, dt$queryHits)
filtered_gr$ENID <- lapply(l, function(df) df$ensid)
mgr <- as.data.frame(filtered_gr) %>% tbl_df()
md_raw <- mgr %>% dplyr::select(transcript_id, ENID) %>% dplyr::rename(modName = transcript_id, ENSEMBL = ENID)
class(md_raw$ENSEMBL) <- 'list'

# remain methylation sites which only are correspond to one gene & build the raw m_dict
for (i in 1:length(md_raw$ENSEMBL)){
  if(length(md_raw$ENSEMBL[[i]]) > 1){
    md_raw$ENSEMBL[[i]] <- NA
  }
}
md_raw$ENSEMBL <- unlist(md_raw$ENSEMBL)
md_raw0 <- md_raw[!is.na(md_raw$ENSEMBL),]

# annotate on each ensembl id
library(org.Hs.eg.db)
library(GO.db)
eid2go <- AnnotationDbi::select(org.Hs.eg.db, keys = md_raw0$ENSEMBL, keytype = 'ENSEMBL', columns = 'GO')
md0 <- left_join(md_raw0, eid2go)

# list of the methylation sites with their corresponding GO terms 
m2go0 <- split(md0, md0$modName)
write_rds(m2go0, 'm_dict_raw.rds')

# the reverse dictionary & count the site number under each EN ID
md_bygene <- split(md_raw0, md_raw0$ENSEMBL)
nsites_bygene <- sapply(md_bygene, function(df) nrow(df))
# select the genes which contain more than one methylation site for further clustering
md_to_retain <- md_bygene[nsites_bygene  == 1] %>% dplyr::bind_rows() %>% 
  dplyr::rename(gene = ENSEMBL, members = modName) %>% 
  mutate(mc = 1L, members = as.list(members)) %>% 
  dplyr::select(gene, mc, members)

md_to_merge <- md_bygene[nsites_bygene > 1]
write_rds(md_to_merge, 'md2merge.rds')
merged.sites <- lapply(md_to_merge, function(df) df$modName %>% unique) %>% unlist %>% unique
testm <- m[,merged.sites]

# to get the correlation adjacency matrix with the methylation sites which need to merge
corm <- cor(testm, method = 'spearman')
corm[corm < 0.3] <- 0 # threshold here!

# Grange of the merged sites
test_mgr <- rmgr[rmgr$transcript_id %in% merged.sites,]
test_mgr <- resize(x = test_mgr, fix = 'center', width = 100)
gr <- tbl_df(as.data.frame(test_mgr)) %>% dplyr::select(seqnames, start, end, strand, transcript_id) %>%  
  dplyr::rename(modName = transcript_id) %>% GRanges()

make_graphs <- function(df, g = gr, cm = corm){
  sites <- unique(df$modName)
  
  # retreive relevant sites
  gr0 <- g[g$modName %in% sites]
  gr0$n <- seq(length(gr0))
  
  # test if overlap in small subgroup
  ov <- findOverlaps(gr0, gr0) %>% as.data.frame()
  same <- mapply(function(a,b) identical(a,b), a = ov[,1], b = ov[,2])
  real_ov <- ov[!same, ]
  if (nrow(real_ov) == 0) {o <- real_ov}
  else {o <- apply(real_ov, c(1,2), 
                   function(n) gr0$modName[match(n, gr0$n)]) %>% as.data.frame}
  
  # retrieve correlation values
  o$cor <- mapply(function(a,b) cm[a,b], a = o$queryHits, b = o$subjectHits )
  o <- o[o$cor > 0,]
  o <- o[!duplicated(o),]
  g <- graph.data.frame(o, directed = F)
  g <- igraph::simplify(g, remove.multiple = T)
}

md_graph <- lapply(md_to_merge, make_graphs)
graph_vs <- sapply(md_graph, function(g) vcount(g))

graphs_to_merge <- md_graph[graph_vs > 0]
graphs_to_retain <- md_graph[graph_vs == 0]

merge_site <- function(g2m){
  members <- cluster_fast_greedy(g2m) %>%  igraph::groups() 
  out <- tibble(mc = names(members), members)
}

clean_cluster <- function(cl){
  o <- dplyr::bind_rows(cl)
  o$gene <- rep(names(cl), sapply(cl, function(df) nrow(df)))
  #o$cid <- paste(o$gene, o$mc, sep = '_')
  o %>% dplyr::select(gene, mc, members) %>% mutate(mc = as.integer(mc))
}

t <- lapply(graphs_to_merge, merge_site) %>% clean_cluster()

clean_unclustered <- function(cl = graphs_to_retain, m2m = md_to_merge){
  m2m <- lapply(m2m, function(df) {
    df$mc <- seq(nrow(df))
    df })
  o <- m2m[names(m2m) %in% names(cl)] %>% dplyr::bind_rows() %>% 
    dplyr::rename(members = modName, gene = ENSEMBL) %>% 
    dplyr::select(gene, mc, members)
  class(o$members) <- 'list'
  o
}

t1 <- clean_unclustered()

t_all <- bind_rows(t, t1, md_to_retain)
t_all$cid <- paste(t_all$gene, t_all$mc, sep = '_')
write_rds(t_all, 'merged_sites.rds')

ms <- read_rds('merged_sites.rds')
md2merge <- read_rds('md2merge.rds')
mdict <- md_raw0
ld <- split(mdict, mdict$ENSEMBL)
ls <- split(ms, ms$gene)
ld_nm <- sapply(ld, function(df) nrow(df)) 
ls_nm <- sapply(ls, function(df) nrow(df))
if_same <- mapply(function(a,b) identical(a,b), a = ld_nm, b = ls_nm) 
p_gid <- names(if_same)[!if_same]
ldp <- ld[!if_same]
lsp <- ls[!if_same]
sites_ldp <- lapply(ldp, function(df) df$modName)
sites_lsp <- lapply(lsp, function(df) df$members %>% unlist %>% unique)
if_same2 <- mapply(function(a,b) identical(length(a), length(b)), a = sites_ldp, b = sites_lsp)
sites_ldp1 <- sites_ldp[!if_same2] 
sites_lsp1 <- sites_lsp[!if_same2] 
missing_sites_by_gene <- mapply(function(d,s) d[ !d %in% s], d = sites_ldp1, s = sites_lsp1)

ls_safe <- ls[!names(ls) %in% names(missing_sites_by_gene)]
ls_tomod <- ls[names(ls) %in% names(missing_sites_by_gene)]
identical(names(ls_tomod), names(missing_sites_by_gene))

ls_modfin <- Map(function(sites, df){
  out_sites <- sapply(sites, as.list) %>% unname
  new_sites <- c(df$members, out_sites)
  g <- df$gene %>% unique
  o <- tibble(gene = g, mc = seq(new_sites), members = new_sites)
  o$cid <- paste(o$gene, o$mc, sep = '_')
  o
} ,
sites = missing_sites_by_gene,
df = ls_tomod)

total_l <- c(ls_safe, ls_modfin)
merged_sites_mod <- bind_rows(total_l)
write_rds(merged_sites_mod,'merged_sites.rds')

m2m <- function(test_m = m, mg = merged_sites_mod){
  c <- matrix(ncol = nrow(mg), nrow = nrow(test_m),
              dimnames = list(rownames(test_m), mg$cid))
  for (i in seq(mg$members)){
    sites <- mg$members[i] %>% unlist
    #browser()
    if (length(sites) == 1){
      c[,i] <- test_m[,sites]
    }
    else{
      c[,i] <- rowMeans(test_m[,sites])
    }
  }
  c
}

new_testm <- m2m()
apply(new_testm, 2, mean)
write_rds(new_testm, "new_testm.rds")

merged_sites <- read_rds("merged_sites.rds")
new_testm <- read_rds("new_testm.rds")

new_testm <- new_testm[match(sort(rownames(new_testm)),rownames(new_testm)),]

# ENID2GO
ENID2GO_BP <- inverseList(annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "ensembl"))
ENID2GO_CC <- inverseList(annFUN.org("CC", mapping = "org.Hs.eg.db", ID = "ensembl"))
ENID2GO_MF <- inverseList(annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "ensembl"))
save(ENID2GO_BP,ENID2GO_CC,ENID2GO_MF,file = "ENID2GO.Rdata")
load("ENID2GO.Rdata")

merged_sites_unlist = merged_sites[rep(1:nrow(merged_sites),elementNROWS(merged_sites$members)),]
merged_sites_unlist$members_unlist = unlist(merged_sites$members)
merged_sites <- merged_sites_unlist[-3]
merged_sites <- arrange(merged_sites, as.numeric(gsub("m6A_site_", "", merged_sites$members_unlist)))
merged_sites <- merged_sites[-2]
t <- split(merged_sites,merged_sites$cid)
modName <- sapply(t, function(a) a[,3])
names(modName) <- names(t)
cid <- names(t)

ensg <- lapply(split(merged_sites,merged_sites$cid),function(x) unique(x$gene))
a <- NULL
for (i in 1:length(ensg)){
  a <- append(a, ensg[[i]])
}
ENSG <- as.tibble(cbind(names(ensg),a))
colnames(ENSG) <- c("cid","gene")
write_rds(ENSG, "ENSG.rds")

ENID <- ENSG[match(cid,ENSG$cid),]$gene
m_dict <- cbind(cid,ENID)
m_dict <- as.tibble(m_dict)
m_dict$modName <- modName

# build co-methylation network (key parameters: quant: quantile percentage; pcut: p-adjust cutoff)
adjmake <- function(x, y, quant, pcut){
  cor <- cor(x, y, method = 'pearson')
  rawp <- corPvalueFisher(cor, 38)
  mt <- mt.rawp2adjp(rawp, proc = 'BH')
  #p.adjust(rawp, method = "fdr", n = length(rawp))
  adj <- mt$adjp[,2]
  adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
  cor <- abs(cor)
  scc_cut <- quantile(cor, quant)
  cor[adjp > pcut | cor < scc_cut] <- 0
  cor[cor > 0] <- 1
  diag(cor) <- 0
  cor
}
adj_meth <- adjmake(x = new_testm, y = new_testm, quant = 0.80, pcut = 0.01)
write_rds(adj_meth,"adj_meth.rds")
adj_meth <- read_rds("adj_meth.rds")

# get the igraph format of the co-methylation network
gme <- graph.adjacency(adj_meth, mode = 'undirected')
write_rds(gme,"gme.rds")

# degree distribution of co-methylation network
png(filename = 'degree-distribution.png', width = 700, height = 550, units = 'px')

#dd <- degree_distribution(gme)
dd <- degree_distribution(gme)
plot(log10(seq(length(dd))),log10(dd), 
     xlab = 'log10(degree)', ylab = 'log10(frequency)', 
     main = 'degree distribution of co-methylation network', pch = 20, col = 'red')
avg.degree <- mean(degree(gme, V(gme)))
abline(v = avg.degree, lwd = 1, lty = 2)
dev.off()

# ================================================================================
# ============================== HUB-BASED METHOD ================================
# ================================================================================
# find neighbors
sites <- colnames(new_testm)
find_nb <- function(namelist = sites, g = gme){
  l = list()
  for (i in 1:length(sites)){
    site <- sites[i]
    l[[i]] <- as_ids(neighbors(g, site))
  }
  l
}
m_dict$nb <- find_nb(m_dict$cid)
m_dict$nb_ENID <- lapply(m_dict$nb, function(x) unique(ENSG[match(x, ENSG$cid),]$gene))
write_rds(m_dict,"m_dict.rds")

# only select the gene which have more than 9 neighbors
m_dict_trim <- dplyr::filter(m_dict, lengths(nb_ENID) > 9)

# m_dict_BP_trim$entrez <- lapply(m_dict_BP_trim$nb_ENID, function(x)
# bitr(x, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE))
# to get the corresponding ENTREZ ID on each neighbor
m_dict_trim$entrez <- lapply(m_dict_trim$nb_ENID, function(a)
  unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = a, columns ="ENTREZID", keytype="ENSEMBL")$ENTREZID)))
write_rds(m_dict_trim,"m_dict_trim.rds")

# create the basic m_dict on 3 various GO aspects
m_dict_BP_trim <- m_dict_trim
m_dict_BP_trim$GO_exact_BP <- ENID2GO_BP[match(m_dict_BP_trim$ENID, names(ENID2GO_BP))]
fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
makeslim <- function(golist , slim. = slim){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'BP'))
      if (!all(unique(temp$Count) == 0)){
        ans[[i]] <- rownames(temp[temp$Count > 1,])
      }
      else {ans[[i]] <- character()}
    }
    else {ans[[i]] <- character()}
  }
  ans
}
m_dict_BP_trim$GO_slim <- makeslim(m_dict_BP_trim$GO_exact_BP)
write_rds(m_dict_BP_trim, 'm_dict_BP_trim.rds')

m_dict_CC_trim <- m_dict_trim
m_dict_CC_trim$GO_exact_CC <- ENID2GO_CC[match(m_dict_CC_trim$ENID, names(ENID2GO_CC))]
makeslim <- function(golist , slim. = slim){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'CC'))
      if (!all(unique(temp$Count) == 0)){
        ans[[i]] <- rownames(temp[temp$Count > 1,])
      }
      else {ans[[i]] <- character()}
    }
    else {ans[[i]] <- character()}
  }
  ans
}
m_dict_CC_trim$GO_slim <- makeslim(m_dict_CC_trim$GO_exact_CC)
write_rds(m_dict_CC_trim, 'm_dict_CC_trim.rds')

m_dict_MF_trim <- m_dict_trim
m_dict_MF_trim$GO_exact_MF <- ENID2GO_MF[match(m_dict_MF_trim$ENID, names(ENID2GO_MF))]
makeslim <- function(golist , slim. = slim){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'MF'))
      if (!all(unique(temp$Count) == 0)){
        ans[[i]] <- rownames(temp[temp$Count > 1,])
      }
      else {ans[[i]] <- character()}
    }
    else {ans[[i]] <- character()}
  }
  ans
}
m_dict_MF_trim$GO_slim <- makeslim(m_dict_MF_trim$GO_exact_MF)
write_rds(m_dict_MF_trim, 'm_dict_MF_trim.rds')

# ===================================================================
# topGO annotation

m_dict_BP_trim <- read_rds("m_dict_BP_trim.rds")
test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'fisher test')
load("geneID2GO.RData")

have <- names(geneID2GO_BP)
annGO <- function(ids, geneID2GO. = geneID2GO_BP){
  geneList <- factor(as.integer(have %in% ids))
  names(geneList) <- have
  GOdata <- new("topGOdata", ontology = 'BP', 
                allGenes = geneList, annot = annFUN.gene2GO,
                gene2GO = geneID2GO.)
  resultFisher <- getSigGroups(GOdata, test.stat)
  t <- sort(score(resultFisher))
  t <- t[t < 0.1] 
  t
}
m_dict_BP_trim$GO_predict <- mapply(annGO, ids = m_dict_BP_trim$entrez)
write_rds(m_dict_BP_trim, "dict_trim_BP.rds")

m_dict_CC_trim <- read_rds("m_dict_CC_trim.rds")
test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'fisher test')
load("geneID2GO.RData")

have <- names(geneID2GO_CC)
annGO <- function(ids, geneID2GO. = geneID2GO_CC){
  geneList <- factor(as.integer(have %in% ids))
  names(geneList) <- have
  GOdata <- new("topGOdata", ontology = 'CC', 
                allGenes = geneList, annot = annFUN.gene2GO,
                gene2GO = geneID2GO.)
  resultFisher <- getSigGroups(GOdata, test.stat)
  t <- sort(score(resultFisher))
  t <- t[t < 0.1] 
  t
}

m_dict_CC_trim$GO_predict <- mapply(annGO, ids = m_dict_CC_trim$entrez)
write_rds(m_dict_CC_trim, "dict_trim_CC.rds")

m_dict_MF_trim <- read_rds("m_dict_MF_trim.rds")
test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'fisher test')
load("geneID2GO.RData")

have <- names(geneID2GO_MF)
annGO <- function(ids, geneID2GO. = geneID2GO_MF){
  geneList <- factor(as.integer(have %in% ids))
  names(geneList) <- have
  GOdata <- new("topGOdata", ontology = 'MF', 
                allGenes = geneList, annot = annFUN.gene2GO,
                gene2GO = geneID2GO.)
  resultFisher <- getSigGroups(GOdata, test.stat)
  t <- sort(score(resultFisher))
  t <- t[t < 0.1] 
  t
}
m_dict_MF_trim$GO_predict <- mapply(annGO, ids = m_dict_MF_trim$entrez)
write_rds(m_dict_MF_trim, "dict_trim_MF.rds")

# =====================================================================
# clusterProfiler annotation

background <- unique(unlist(m_dict_trim$entrez))
# BP annotation
i <- 0
GO_BP_predict <- function(a){
  i <<- i + 1; print(i)
  GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
  #GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, universe = background, readable = TRUE)
  if (nrow(GO_pre) == 0){
    GO_pre <- NULL
  }
  else {
    GO_pre <- GO_pre[,c(1:6,9)]
  }
}
BP <- lapply(m_dict_BP_trim$nb_ENID[1:10], GO_BP_predict)
write_rds(BP, "BP_term.rds")

# CC annotation
i <- 0
GO_CC_predict <- function(a){
  i <<- i + 1; print(i)
  GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")
  #GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "CC", pvalueCutoff = 0.01, universe = background, readable = TRUE)
  if (nrow(GO_pre) == 0){
    GO_pre <- NULL
  }
  else {
    GO_pre <- GO_pre[,c(1:6,9)]
  }
}
CC <- lapply(m_dict_CC_trim$nb_ENID, GO_CC_predict)
write_rds(CC, "CC_term.rds")

# MF annotation
i <- 0
GO_MF_predict <- function(a){
  i <<- i + 1; print(i)
  GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
  #GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "MF", pvalueCutoff = 0.01, universe = background, readable = TRUE)
  if (nrow(GO_pre) == 0){
    GO_pre <- NULL
  }
  else {
    GO_pre <- GO_pre[,c(1:6,9)]
  }
}
MF <- lapply(m_dict_MF_trim$nb_ENID, GO_MF_predict)
write_rds(MF, "MF_term.rds")

m_dict_trim <- read_rds("m_dict_trim.rds")
background <- unique(unlist(m_dict_trim$entrez))
i <- 0
GO_ALL_predict <- function(a){
  i <<- i + 1; print(i)
  GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", universe = background, ont = "ALL")
  #GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, universe = background, readable = TRUE)
  if (nrow(GO_pre) == 0){
    GO_pre <- NULL
  }
  else {
    GO_pre <- GO_pre[,c(1:6,9)]
  }
}
GO_term <- lapply(m_dict_trim$nb_ENID, GO_ALL_predict)
write_rds(GO_term, "GO_all_term.rds")
# ==========================================================================================================================
# pvalue cut
p_vals = sort(c(10^(-1), 10^(-1.5), 10^(-2), 10^(-2.5),
                10^(-3), 10^(-3.5), 10^(-4), 10^(-4.5),
                10^(-5), 10^(-5.5), 10^(-6), 10^(-6.5)))

cutpval <- function(lg,p){
  lapply(lg,function(n){
    names(n[n < p])
  })
}

fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
enslim <- function(go){
  apply(go, c(1,2), function(char){
    char <- unname(unlist(char, recursive = FALSE))
    if (length(char)>0){
      mygo <- GOCollection(char)
      temp <- goSlim(mygo, slim, 'MF')
      if (!all(unique(temp$Count) == 0)){
        #browser()
        temp <- temp[temp$Count>0,]
        temp <- temp[with(temp, order(-Count, -Percent)),]
        rownames(temp)
      }
      else {character()}
    }
    else {character()}
  })
}

sel.slim <- function(m, nt, selGO){
  l <- split(m, col(m))
  names(l) <- colnames(m)
  l <- lapply(l, function(ll){
    names(ll) <- rownames(selGO)
    ll})
  lcut <- lapply(l, function(sl){
    lapply(sl, function(s){
      ans <- s[1:nt]
      ans[!is.na(ans)]
    })
  })
  lcut
}

nt_vec <- c(2,4,6,8,10,12,14,16,18,20)
#nt_vec <- seq(10)

make_pa2 <- function(pval = p_vals, nts = nt_vec, eval){
  sel.GO <- mapply(cutpval, p = pval, MoreArgs = list(lg = eval$GO_predict))
  colnames(sel.GO) <- as.character(pval)
  rownames(sel.GO) <- eval$cid
  pa2 <- mapply(sel.slim, nt = nts, 
                MoreArgs = list(m = enslim(sel.GO),
                                selGO = sel.GO))
  colnames(pa2) <- as.character(nts)
  pa2
}

m_dict_BP_trim <- read_rds("dict_trim_BP.rds")
m_dict_CC_trim <- read_rds("dict_trim_CC.rds")
m_dict_MF_trim <- read_rds("dict_trim_MF.rds")
pa2_true_BP <- make_pa2(eval = m_dict_BP_trim)
pa2_true_CC <- make_pa2(eval = m_dict_CC_trim)
pa2_true_MF <- make_pa2(eval = m_dict_MF_trim)
write_rds(pa2_true_BP, "pa2_true_BP.rds")
write_rds(pa2_true_CC, "pa2_true_CC.rds")
write_rds(pa2_true_MF, "pa2_true_MF.rds")
# dotplot
#dotplot(enrichGO(gene = m_dict_BP_trim$entrez[[1]], 
#                 OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, universe = background, readable = TRUE))
# view the enriched GO induced graph
#a <- enrichGO(gene = m_dict_BP_trim$entrez[[1]], OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, readable = TRUE)
#enrichMap(a)

# permutated network
tt <- 100*length(E(gme))
gran <- rewire(gme, keeping_degseq(niter = tt))
m_dict_ran <- m_dict

find_nb <- function(namelist = sites, g = gran){
  l = list()
  for (i in 1:length(sites)){
    site <- sites[i]
    l[[i]] <- as_ids(neighbors(g, site))
  }
  l
}
m_dict_ran$nb <- find_nb(m_dict_ran$cid)
m_dict_ran$nb_ENID <- lapply(m_dict_ran$nb, function(x) unique(ENSG[match(x, ENSG$cid),]$gene))
write_rds(m_dict_ran,"m_dict_ran.rds")

m_dict_ran$nb_ENID_real <- m_dict$nb_ENID 
m_dict_ran_trim <- dplyr::filter(m_dict_ran, lengths(nb_ENID_real) > 9)
m_dict_ran_trim$nb_ENID_real <- NULL

m_dict_ran_trim$entrez <- lapply(m_dict_ran_trim$nb_ENID, function(a)
  unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = a, columns ="ENTREZID", keytype="ENSEMBL")$ENTREZID)))
write_rds(m_dict_ran_trim,"m_dict_ran_trim.rds")

m_dict_ran_BP_trim <- m_dict_ran_trim
m_dict_ran_BP_trim$GO_exact_BP <- ENID2GO_BP[match(m_dict_ran_BP_trim$ENID, names(ENID2GO_BP))]
fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
makeslim <- function(golist , slim. = slim){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'BP'))
      if (!all(unique(temp$Count) == 0)){
        ans[[i]] <- rownames(temp[temp$Count > 1,])
      }
      else {ans[[i]] <- character()}
    }
    else {ans[[i]] <- character()}
  }
  ans
}
m_dict_ran_BP_trim$GO_slim <- makeslim(m_dict_ran_BP_trim$GO_exact_BP)
write_rds(m_dict_ran_BP_trim, 'm_dict_ran_BP_trim.rds')

m_dict_ran_CC_trim <- m_dict_ran_trim
m_dict_ran_CC_trim$GO_exact_CC <- ENID2GO_CC[match(m_dict_ran_CC_trim$ENID, names(ENID2GO_CC))]
makeslim <- function(golist , slim. = slim){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'CC'))
      if (!all(unique(temp$Count) == 0)){
        ans[[i]] <- rownames(temp[temp$Count > 1,])
      }
      else {ans[[i]] <- character()}
    }
    else {ans[[i]] <- character()}
  }
  ans
}
m_dict_ran_CC_trim$GO_slim <- makeslim(m_dict_ran_CC_trim$GO_exact_CC)
write_rds(m_dict_ran_CC_trim, 'm_dict_ran_CC_trim.rds')

m_dict_ran_MF_trim <- m_dict_ran_trim
m_dict_ran_MF_trim$GO_exact_MF <- ENID2GO_MF[match(m_dict_ran_MF_trim$ENID, names(ENID2GO_MF))]
makeslim <- function(golist , slim. = slim){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'MF'))
      if (!all(unique(temp$Count) == 0)){
        ans[[i]] <- rownames(temp[temp$Count > 1,])
      }
      else {ans[[i]] <- character()}
    }
    else {ans[[i]] <- character()}
  }
  ans
}
m_dict_ran_MF_trim$GO_slim <- makeslim(m_dict_ran_MF_trim$GO_exact_MF)
write_rds(m_dict_ran_MF_trim, 'm_dict_ran_MF_trim.rds')

#=====================================================================

pa2_true_BP <- read_rds("pa2_true_BP.rds")
pa2_true_CC <- read_rds("pa2_true_CC.rds")
pa2_true_MF <- read_rds("pa2_true_MF.rds")
pa2_ran_BP <- read_rds("pa2_ran_BP.rds")
pa2_ran_CC <- read_rds("pa2_ran_CC.rds")
pa2_ran_MF <- read_rds("pa2_ran_MF.rds")

el <- function(l, eval.){
  e0 <- eval. %>% dplyr::select(cid, GO_slim)
  l <- l[[1]]
  e0$predict <- l[match(e0$cid, names(l))]
  e0$hit <- mapply(function(a,b) intersect(a,b), a = e0$GO_slim, b = e0$predict)
  e0 <- e0[lengths(e0$GO_slim)>0 & lengths(e0$predict)>0,]
  pre <- sum(lengths(e0$hit)) / sum(lengths(e0$GO_slim))
  spe <- sum(lengths(e0$hit)) / sum(lengths(e0$predict))
  out = list(e0, pre, spe)
  names(out) <- c('alldat', 'precision', 'specificity')
  out
}

pa2all_true_BP <- apply(pa2_true_BP, c(1,2), el, eval. = m_dict_BP_trim)
pa2all_true_CC <- apply(pa2_true_CC, c(1,2), el, eval. = m_dict_CC_trim)
pa2all_true_MF <- apply(pa2_true_MF, c(1,2), el, eval. = m_dict_MF_trim)
pa2all_ran_BP <- apply(pa2_ran_BP, c(1,2), el, eval. = m_dict_ran_BP_trim)
pa2all_ran_CC <- apply(pa2_ran_CC, c(1,2), el, eval. = m_dict_ran_CC_trim)
pa2all_ran_MF <- apply(pa2_ran_MF, c(1,2), el, eval. = m_dict_ran_MF_trim)

# Pre & Spe comparison
save(m_dict_BP_trim, m_dict_ran_BP_trim, 
     pa2all_ran_BP, pa2all_true_BP, 
     pa2_true_BP, pa2_ran_BP, file = 'hub-BP.RData')
save(m_dict_CC_trim, m_dict_ran_CC_trim, 
     pa2all_ran_CC, pa2all_true_CC, 
     pa2_true_CC, pa2_ran_CC, file = 'hub-CC.RData')
save(m_dict_MF_trim, m_dict_ran_MF_trim, 
     pa2all_ran_MF, pa2all_true_MF, 
     pa2_true_MF, pa2_ran_MF, file = 'hub-MF.RData')

better_BP <- function(r = pa2all_ran_BP, t = pa2all_true_BP){
  tspe <- apply(t, c(1,2), function(l3){l3[[1]]$specificity})
  tpre <- apply(t, c(1,2), function(l3){l3[[1]]$precision})
  rspe <- apply(r, c(1,2), function(l3){l3[[1]]$specificity})
  rpre <- apply(r, c(1,2), function(l3){l3[[1]]$precision})
  comp <- list(true = list(tpre, tspe),
               random = list(rpre, rspe),
               inc = list(p = tpre - rpre,
                          s = tspe - rspe))
  comp
}
better_CC <- function(r = pa2all_ran_CC, t = pa2all_true_CC){
  tspe <- apply(t, c(1,2), function(l3){l3[[1]]$specificity})
  tpre <- apply(t, c(1,2), function(l3){l3[[1]]$precision})
  rspe <- apply(r, c(1,2), function(l3){l3[[1]]$specificity})
  rpre <- apply(r, c(1,2), function(l3){l3[[1]]$precision})
  comp <- list(true = list(tpre, tspe),
               random = list(rpre, rspe),
               inc = list(p = tpre - rpre,
                          s = tspe - rspe))
  comp
}
better_MF <- function(r = pa2all_ran_MF, t = pa2all_true_MF){
  tspe <- apply(t, c(1,2), function(l3){l3[[1]]$specificity})
  tpre <- apply(t, c(1,2), function(l3){l3[[1]]$precision})
  rspe <- apply(r, c(1,2), function(l3){l3[[1]]$specificity})
  rpre <- apply(r, c(1,2), function(l3){l3[[1]]$precision})
  comp <- list(true = list(tpre, tspe),
               random = list(rpre, rspe),
               inc = list(p = tpre - rpre,
                          s = tspe - rspe))
  comp
}

comp_BP <- better_BP()
comp_CC <- better_CC()
comp_MF <- better_MF()
write_rds(comp_BP, 'compare-hub-BP.rds')
write_rds(comp_CC, 'compare-hub-CC.rds')
write_rds(comp_MF, 'compare-hub-MF.rds')

# ====================================================
i <- 0
GO_BP_predict <- function(a){
  i <<- i + 1; print(i)
  GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
  #GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, universe = background, readable = TRUE)
  if (nrow(GO_pre) == 0){
    GO_pre <- NULL
  }
  else {
    GO_pre <- GO_pre[,c(1:6,9)]
  }
}
BP <- lapply(m_dict_ran_BP_trim$nb_ENID, GO_BP_predict)

#==========gprofiler--------Proxy Error==========
#j <- 0
#GO_gprofiler <- lapply(m_dict_trim$nb_ENID[1:50], function(a){
#  j <<- j + 1; print(j)
#  gprofiler(a, organism = "hsapiens", src_filter = "GO")})
#================================================

#KEGG analysis
i <- 0
KEGG_predict <- function(a){
  i <<- i + 1; print(i)
  KEGG_pre <- enrichKEGG(gene = a, organism = "hsa")
  if (nrow(KEGG_pre) == 0){
    KEGG_pre <- NULL
  }
}
KEGG <- lapply(m_dict_trim$entrez[1:200], KEGG_predict)
write_rds(KEGG, "KEGG_term.rds")

bitr_kegg(geneID = m_dict_trim$entrez[1], fromType = "ncbi-geneid", toType = "Path", organism = "hsa")
# browse the KEGG pathway on KEGG official database with the enrich genes highlighted
browseKEGG(KEGG[[253]], "hsa00561")


# =============================================================================
# ========================== MODULE-BASED METHOD ==============================
# =============================================================================
# create the ABC-format file for MCL clustering
adj <- exportNetworkToCytoscape(adj_meth)$edgeData
write_tsv(adj,"adj.tsv")
adj <- read_tsv("adj.tsv")
adj <- adj[,1:3]
adj <- as.matrix(adj)
for (i in 1:nrow(adj)){
  adj [i,] <- paste0(adj[i,1]," ",adj[i,2]," ",adj[i,3])
}
adj <- as.data.frame(adj[,1])
write_tsv(adj,"adjfile.tsv")

# ================== mcl clustering (in bash shell) =====================
# mcl adjfile.tsv -I 1.8 --abc -o out.adjfile1.8
# mcl adjfile.tsv -I 1 --abc -o out.adjfile1
# =======================================================================

adj1.8 <- as.matrix(read_csv("out.adjfile1.8", col_names = FALSE, na = c("NA")))
adj1.4 <- as.matrix(read_csv("out.adjfile1.4", col_names = FALSE, na = c("NA")))

mclmodule <- NULL
for (i in 1:length(adj1.4)){
  mclmodule[[i]] <- strsplit(adj1.4[i],"\t")
}

mcl_module <- NULL
for (i in 1:length(mclmodule)){
  if (elementNROWS(mclmodule[[i]]) > 9){
    mcl_module <- append(mcl_module, mclmodule[[i]])
  }
}

cluster <- NULL
meth <- NULL
for (j in 1:length(mcl_module)){
  cluster <- append(cluster, rep(j,length(mcl_module[[j]])))
  meth <- append(meth, mcl_module[[j]])
}

# create the m_dict for module-based 
module_num <- c(1:length(unique(cluster)))
module_size <- unlist(lapply(mcl_module, function(a) length(a)))
module_dict <- as.tibble(cbind(module_num, module_size))
module_dict$cid <- mcl_module
module_dict$ENID <- lapply(module_dict$cid, function(x) unique(ENSG[match(x, ENSG$cid),]$gene))
module_dict$entrez <- lapply(module_dict$ENID, function(a) 
  unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = a, columns ="ENTREZID", keytype="ENSEMBL")$ENTREZID)))
module_dict_trim <- dplyr::filter(module_dict, lengths(entrez) > 3)
module_background <- unique(unlist(module_dict_trim$entrez))

#enrichGO(gene = module_dict_trim$entrez[[1]], 
#         OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, universe = module_background, readable = TRUE)

i <- 0
GO_predict <- function(a){
  i <<- i + 1; print(i)
  GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "BP", readable = TRUE)
  if (nrow(GO_pre) == 0){
    GO_pre <- NULL
  }
  else {
    GO_pre <- GO_pre[,c(1:6,9)]
  }
}
module_GO <- lapply(module_dict_trim$entrez, GO_predict)
write_rds(module_GO, "module_GO_BP.rds")

# KEGG analysis
i <- 0
KEGG_predict <- function(a){
  i <<- i + 1; print(i)
  KEGG_pre <- enrichKEGG(gene = a, organism = "hsa")
}
KEGG_module <- lapply(module_dict_trim$entrez, KEGG_predict)
write_rds(KEGG_module, "KEGG_module_term.rds")