library(readr)
library(igraph)
library(dplyr)
library(igraph)
library(tibble)
library(GenomicRanges)
library(BiocParallel)
register(SerialParam())
m <- read_rds('sdm_filtered_testm.rds')
rmgr <- read_rds('mgr18w.rds')
sel.gr <- rmgr[rmgr$transcript_id %in% colnames(m)]


library(EnsDb.Hsapiens.v75)
g <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(g) <- 'UCSC' 
o <- findOverlaps(sel.gr, g)
sel.gr$ENS_id <- NA



dt <- as.data.frame(o) %>% tbl_df()
hits2eid <- tibble(hit = seq(length(g)), eid = g$gene_id)
dt$ensid <- hits2eid$eid[match(dt$subjectHits, hits2eid$hit)]
l <- split(dt, dt$queryHits)
sel.gr$ENS_id <- lapply(l, function(df){
  df$ensid
})
rm(hits2eid, l, dt, g)


md_raw <- as.data.frame(sel.gr) %>% tbl_df() %>% 
  dplyr::select(transcript_id, ENS_id) %>% 
  dplyr::rename(modName = transcript_id, ENSEMBL = ENS_id)
class(md_raw$ENSEMBL) <- 'list'

# a flattened version of m_dict
md_raw0 <- tidyr::unnest(md_raw)
md_bygene <- split(md_raw0, md_raw0$ENSEMBL)
nsites_bygene <- sapply(md_bygene, function(df) nrow(df))

md_to_retain <- md_bygene[nsites_bygene  == 1] %>% dplyr::bind_rows() %>% 
  dplyr::rename(gene = ENSEMBL, members = modName) %>% 
  mutate(mc = 1L, members = as.list(members)) %>% 
  dplyr::select(gene, mc, members)

md2merge <- md_bygene[nsites_bygene > 1]
rm(nsites_bygene)



# make a pc-friendly testset
# md2merge <- head(md_to_merge)
relevant.sites <- lapply(md2merge, function(df) df$modName %>% unique) %>% 
  unlist %>% unique
testm <- m[,relevant.sites]
corm <- cor(testm, method = 'spearman')
corm[corm < 0.3] <- 0 # threshold here!

# prepare grange

test_mgr <- rmgr[rmgr$transcript_id %in% relevant.sites]
rm(rmgr)
test_mgr <- resize(x = test_mgr, fix = 'center', width = 100)
gr <- as.data.frame(test_mgr) %>% tbl_df() %>% 
  dplyr::select(seqnames, start, end, strand, transcript_id) %>%  
  dplyr::rename(modName = transcript_id) %>% GRanges()


make_graphs <- function(df, g = gr, cm = corm){
  #browser()
  
  sites <- df$modName %>%  unique()
  # retreive relevant sites
  gr0 <- g[g$modName %in% sites]
  gr0$n <- seq(length(gr0))
  
  # test if overlap in small subgroup
  ov <- findOverlaps(gr0, gr0) %>% as.data.frame()
  #browser()
  same <- mapply(function(a,b) identical(a,b), a = ov[,1], b = ov[,2])
  real_ov <- ov[!same, ]
  
  if (nrow(real_ov) == 0) {o <- real_ov}
  else {o <- apply(real_ov, c(1,2), 
                   function(n) gr0$modName[match(n, gr0$n)]) %>% as.data.frame}
  #real_ov <- apply(real_ov, c(1,2), function(n) gr0$modName[match(n, gr0$n)])
  #o <- real_ov[!duplicated(real_ov),]
  
  #browser()
  # retrieve correlation values
  o$cor <- mapply(function(a,b) cm[a,b], a = o$queryHits, b = o$subjectHits )
  o <- o[o$cor > 0,]
  #browser()
  o <- o[!duplicated(o),]
  g <- graph.data.frame(o, directed = F)
  g <- igraph::simplify(g, remove.multiple = T)
  
}

md_graph <- lapply(md2merge, make_graphs)
graph_vs <- sapply(md_graph, function(g) vcount(g))
# may reduce computation load
graphs_to_merge <- md_graph[graph_vs > 0] 
graphs_to_retain <- md_graph[graph_vs == 0] #redundant!

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

clean_unclustered <- function(cl = graphs_to_retain, m2m = md2merge){
  m2m <- lapply(m2m, function(df) {
    df$mc <- seq(nrow(df))
    df })
  o <- m2m[names(m2m) %in% names(cl)] %>% dplyr::bind_rows() %>% 
     dplyr::rename(members = modName, gene=ENSEMBL) %>% 
    dplyr::select(gene, mc, members)
  class(o$members) <- 'list'
  o
}


t <- lapply(graphs_to_merge, merge_site) %>% clean_cluster()
t1 <- clean_unclustered()
t_all <- bind_rows(t, t1, md_to_retain)
t_all$cid <- paste(t_all$gene, t_all$mc, sep = '_')

write_rds(t_all, 'merged_sites.rds')