
# 'merged_sites_by_graph' contains a lot problems,
# this script is to first resolve those missing sites and then 
# get test_m object from merged features.




library(dplyr)
library(readr)
ms <- read_rds('merged_sites.rds')
m <- read_rds('sdm_filtered_testm.rds')
md2merge <- read_rds('md2merge.rds')
look <- function(m) m[1:10, 1:10]

mdict <- read_rds('m_dict_raw.rds') # this is consistent with m object
ld <- split(mdict, mdict$ENSEMBL)
ls <- split(ms, ms$gene)
ld_nm <- sapply(ld, function(df) nrow(df)) 
ls_nm <- sapply(ls, function(df) nrow(df))
if_same1 <- mapply(function(a,b) identical(a,b), a = ld_nm, b = ls_nm) 
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
identical(names(ld_tomod), names(missing_sites_by_gene))

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



# 
m2m <- function(test_m = m, mg = merged_sites_mod, md = mdict){
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
write_rds(new_testm, 'new_testm.rds')
