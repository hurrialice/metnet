library(readr)
library(org.Hs.eg.db)
te<- read_rds('test_exp_f.rds')
gene_names <- colnames(te)
map_ens <- AnnotationDbi::select(og, keys = gene_names, keytype = "SYMBOL", columns = "ENSEMBL" )

no_ens_sym <- map_ens$SYMBOL[is.na(map_ens$ENSEMBL)] %>% unique()


mod_te <- function(m){
  # cut those genes without ensg id
  m0 <- m[, -which(colnames(m) %in% no_ens_sym) ]
  new_colname <- map_ens$EMSEMBL[match(colnames(m0), map_ens$SYMBOL)]
  colnames(m0) <- new_colname
  m0
}

test_exp_new <- mod_te(te)
write_rds(test_exp_new, 'test_exp_fin.rds')
