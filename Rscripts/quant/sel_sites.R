# objective:
# select sites that  - 
# 1. scale to N(0,1)
# 2. top 25% variance
# 3. top 25% mean
# 4. greedy clutsering for sites on same gene!


# ideal, returns 10000 sites with 3000 features.



library(readr)
library(dplyr)
look <- function(m) m[1:10, 1:10]

# scale
uni_m <- read_rds('test_m0921.rds')
nor_m <- preprocessCore::normalize.quantiles(t(uni_m)) %>% t
s.m  <- t(scale(t(nor_m)))
s.m <- apply(nor_m, c(1,2), function(x) (x - mean(x))/sd(x) )
attr(s.m,'dimnames') <- attr(uni_m, 'dimnames')


# select by sd and mean
sdnm_filter <- function(m, qm, qsd){
    
    site_ms <- apply(m, 2, mean)
    site_sds <- apply(m, 2, sd)
    
    m_cut <- quantile(site_ms, 1 - qm)
    sd_cut <- quantile(site_sds, 1 - qsd)
    
    m_retain <- (site_ms > m_cut)
    sd_retain <- (site_sds > sd_cut)
    
    retain <- mapply(function(a,b) a & b, a = m_retain, b = sd_retain)
    m[,retain]
    
}

s1m <- sdnm_filter(s.m, 0.3, 0.3)
m <- preprocessCore::normalize.quantiles(t(s1m)) %>% t

selsites <- colnames(m)
attr(m,'dimnames') <- attr(s1m, 'dimnames')
write_rds(m,'sdm_filtered_testm.rds')



## load msite origin info

rmgr <- read_rds('mgr18w.rds') # from meth-quant2.R
sel.gr <- rmgr[rmgr$transcript_id %in% selsites]


library(EnsDb.Hsapiens.v75)
g <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(g) <- 'UCSC' 
o <- findOverlaps(sel.gr, g)
sel.gr$ENS_id <- NA
rm(rmgr,g)


dt <- as.data.frame(o) %>% tbl_df()
hits2eid <- tibble(hit = seq(length(g)), eid = g$gene_id)
dt$ensid <- hits2eid$eid[match(dt$subjectHits, hits2eid$hit)]
l <- split(dt, dt$queryHits)
sel.gr$ENS_id <- lapply(l, function(df){
    df$ensid
})
rm(hits2eid, l)



mgr <- as.data.frame(sel.gr) %>% tbl_df()
md_raw <- mgr %>% dplyr::select(transcript_id, ENS_id) %>% dplyr::rename(modName = transcript_id, ENSEMBL = ENS_id)
class(md_raw$ENSEMBL) <- 'list'

# a flattened version of m_dict
md_raw0 <- tidyr::unnest(md_raw)


# refer to GO, argublly meaningful ?
library(org.Hs.eg.db)
library(GO.db)
eid2go <- AnnotationDbi::select(org.Hs.eg.db, keys = md_raw0$ENSEMBL, 
                      keytype = 'ENSEMBL', columns = 'GO' )
md0 <- left_join(md_raw0, eid2go)
m2go0 <- split(md0, md0$modName)
write_rds(m2go0, 'm_dict_raw.rds') # a putative methylation site dict


# a reverse dictionary 
md_bygene <- split(md_raw0, md_raw0$ENSEMBL)
nsites_bygene <- sapply(md_bygene, function(df) nrow(df))
md_to_merge <- md_bygene[nsites_bygene > 1]
rm(nsites_bygene)

# calculate correlation among the sites
m <- read_rds('sdm_filtered_testm.rds')
write_rds(md_to_merge, 'md2merge.rds')
