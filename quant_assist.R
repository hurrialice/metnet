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


## merge the reps (first include the )
merge_reps <- function(srr_list, m, pat){
        
        # init container
        mrep <- list()
        
        # first process the reps
        for (i in seq_along(srr_list)){
                reps <- srr_list[[i]]
                mrep[[i]] <- colMeans(m[reps, ])
        }
        # browser()
        # merge the reps into a single matrix
        mrep_m <- do.call(rbind, mrep)
        rownames(mrep_m) <- sapply(check_abbr(srr_list), function(chars) chars[1])
        colnames(mrep_m) <- colnames(m)
        
        # for those lines that do not need to be merged,
        single_srrs <- rownames(m)[! rownames(m) %in% unique(unlist(srr_list))]
        msin_m <- m[single_srrs,]
        rownames(msin_m) <- h$source_abbr[match(rownames(msin_m), h$sra_acc)]
        
        # merge the single and reps
        out_m <- rbind(msin_m, mrep_m)
        o <- out_m[order(rownames(out_m)),]
        rownames(o) <- gsub(pat, "", rownames(o))
        return(o)
        
}



## check if the corresponding source abbr are identical
check_abbr <- function(srr_list){
        possi_abbrs <- list()
        for (i in seq_along(srr_list)){
                reps <- srr_list[[i]]
                possi_abbrs[[i]] <- unique(h$source_abbr[h$sra_acc %in% reps])
        }
        return(possi_abbrs)
}


