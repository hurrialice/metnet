# from "test_exp_raw" file.
# first - identify those problematic samples and delete 
look <- function(m) m[1:10, 1:10]
library(readr)
library(dplyr)
source("quant_assist.R")


# load raw file SRRX*SYM
exp_raw <- read_rds('test_exp_raw.rds')
m_raw <- read_rds('test_m_raw.rds')

# load the smaple dictionary
homo <- read_rds("srr_withgtfs.rds")
h <- homo %>% dplyr::select(wecall, sra_acc)
h$source_abbr <- mapply(function(a,b) gsub(paste0("_",a), "", b), a = h$sra_acc, b = h$wecall)

# split into three matrices
m_ip_raw <- m_raw[h$sra_acc[grepl("_ip_", h$wecall)],]
m_input_raw <- m_raw[h$sra_acc[grepl("_input_", h$wecall)],]

# identify the samples to be deleted, e.g. humanBrain
srrs_to_delete <- c('SRR456936','SRR456937')
delete.from.matrix <- function(orim, del.rows){
        orim[ -which(rownames(orim) %in% del.rows),]
}
exp_del <- delete.from.matrix(exp_raw, srrs_to_delete)
m_ip_del <- delete.from.matrix(m_ip_raw, srrs_to_delete)
m_input_del <- delete.from.matrix(m_input_raw, srrs_to_delete)

### for display and check by eyes (only differ by 'BR' ids, presumably bio-reps)
# check_abbr(input_rep_srrs)
# check_abbr(ip_rep_srrs)

exp_merge <- merge_reps(input_rep_srrs, exp_del, "_input")
ip_merge <- merge_reps(ip_rep_srrs, m_ip_del, "_ip")
input_merge <- merge_reps(input_rep_srrs, m_input_del, "_input")

# calculate the methylation level
meth_merge <- log2((ip_merge + 0.01) / (input_merge + 0.01))
Note <- "here the three stored matrices have gone thorough stringtie quantification, reading from gtf files and merging the techreps and bioreps. However these data has not been normalised. Start your journey here!"
save(exp_merge, ip_merge, input_merge, Note, file = "raw_merge.RData")
