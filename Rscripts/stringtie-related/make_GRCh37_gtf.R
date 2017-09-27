# Objective: make a uniform gene model in gtf
# for merging sites via greedy clustering, and GO analysis

library(EnsDb.Hsapiens.v75)
library(rtracklayer)
e <- exons(EnsDb.Hsapiens.v75,columns = c("tx_id", "gene_id"))
colnames(mcols(e))[1] <- 'transcript_id'
e$type <- 'exon'

cds <- cdsBy(EnsDb.Hsapiens.v75, columns =c ('tx_id', 'gene_id'))
cds <- unlist(cds)
cds$type <- 'cds'
colnames(mcols(cds))[1] <- 'transcript_id'
cds$exon_rank <- NULL

pre_guide <- do.call("c", list(e,cds)) 
rm(e, cds)

t <- keepStandardChromosomes(pre_guide, pruning.mode = 'coarse')
newStyle <- mapSeqlevels(seqlevels(t),"UCSC")
t <- renameSeqlevels(t, newStyle)
export(t, 'ens_v75.gtf')
