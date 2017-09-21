# objective: 
# clean up replicates among srrs.

library(readr)
library(dplyr)
library(tidyr)
test_m <- read_rds('test_m_raw.rds')
test_exp <- read_rds('test_exp_raw.rds')
homo <- read_tsv('sample_info.tsv')
h <- homo %>% dplyr::select(source_abbr, sra_acc, wecall, type, gse_acc, gsm_acc, group_name)

# separate IP from input
abbr <- h$wecall[match(rownames(test_m), h$sra_acc)]
ifip <- grepl("_ip_", abbr)
test_m_ip <- test_m[ifip,]
test_m_input <- test_m[!(ifip),]
rm(test_m)

# remove human brain sample
input_names <- grep('_input', abbr, value = T)[-11] #identified by raw eye!
ip_names <- grep('_ip', abbr, value = T)[-12]
test_m_input <- test_m_input[-which(rownames(test_m_input) == 'SRR456937'),]
test_m_ip <- test_m_ip[-which(rownames(test_m_ip) == 'SRR456936'),]
test_exp <- test_exp[-which(rownames(test_exp) == 'SRR456937'),]

# manual curated! sad!
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

# name col and rows
all_rep_srrs <- unname(c(unlist(input_rep_srrs),unlist(ip_rep_srrs) ))
brain_srrs <- c('SRR456936','SRR456937')
valid_srrs <- h$sra_acc[which(!h$sra_acc %in% brain_srrs)]
single_srrs <- valid_srrs[!valid_srrs %in% all_rep_srrs]
single_conds_name <- paste(h$source_abbr,h$group_name, h$type, sep = '_')[match(single_srrs, h$sra_acc)]
rep_conds_name <- sapply(c(input_rep_srrs, ip_rep_srrs), function(srrs){
    unique(paste(h$source_abbr,h$group_name, h$type, sep = '_')[match(srrs, h$sra_acc)])[1]
}) %>% unname
all_conds <- c(single_conds_name, rep_conds_name)
any(duplicated(all_conds))


# function to merge replicated conditions
del.rep <- function(l, orim){
    rep.srrs <- unlist(l)
    all.srrs <- rownames(orim)
    single_srrs <- all.srrs[!(all.srrs %in% rep.srrs)]
    m_single <- orim[single_srrs,]
    rownames(m_single) <- paste(h$source_abbr,h$group_name, sep = '_')[match(rownames(m_single),h$sra_acc)]
    
    #browser()
    # init
    m_rep <- matrix(nrow = length(l), ncol = ncol(orim))
    colnames(m_rep) <- colnames(orim)
    rownames(m_rep) <- sapply(l, function(srrs){
        unique(paste(h$source_abbr,h$group_name, sep = '_')[match(srrs, h$sra_acc)])[1]
    }) %>% unname
    
    for (i in seq_along(l)){
        srrs_to_merge <- l[[i]]
        #new_name <- paste(srrs_to_merge, collapse = "_") #double check by name
        m_rep[i,] <- colMeans(orim[srrs_to_merge,])
    }
    
    rbind(m_single, m_rep)
}

# operations
unim_input <- del.rep(orim = test_m_input,l = input_rep_srrs)
unim_ip <- del.rep(orim = test_m_ip,l = ip_rep_srrs)
uni_exp <- del.rep(orim = test_exp, l = input_rep_srrs)
stopifnot(identical(rownames(unim_input), rownames(unim_ip)))
stopifnot(identical(colnames(unim_input), colnames(unim_ip)))
uni_m <- log2(unim_ip+0.01) - log2(unim_input+0.01)

write_rds(uni_m, 'test_m0921.rds')
write_rds(uni_exp, 'test_exp0921.rds')
