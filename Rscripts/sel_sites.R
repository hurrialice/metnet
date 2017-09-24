library(readr)
library(igraph)
m <- read_rds('sdm_filtered_testm.rds')


md_to_merge <- read_rds('md2merge.rds')
md2merge <- head(md_to_merge)
relevant.sites <- lapply(md2merge, function(df) df$modName %>% unique) %>% unlist %>% unique

test_sites<- unique(c(sample(colnames(m), 40), relevant.sites))
testm <- m[,test_sites]
corm <- cor(testm, method = 'spearman')
corm[corm < 0.2] <- 0 # threshold here!

# prepare grange
rmgr <- read_rds('mgr18w.rds')
test_mgr <- rmgr[rmgr$transcript_id %in% test_sites]
rm(rmgr)
test_mgr <- resize(x = test_mgr, fix = 'center', width = 50)
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
  g <- simplify(g, remove.multiple = T)
  
}

md_graph <- lapply(md2merge, make_graphs)
graph_vs <- sapply(md_graph, function(g) vcount(g))
# may reduce computation load
graphs_to_merge <- md_graph[graph_vs > 0] 
graphs_to_retain <- md_graph[graph_vs == 0] #redundant!

merge_site <- function(g2m){
  members <- cluster_fast_greedy(g2m) %>%  groups() 
  out <- tibble(mc = names(members), members)
}

lapply(graphs_to_merge, merge_site)
