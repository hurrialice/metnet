# this version abandoned the use of future and uses a smaller 
# simulated data set (1/100 in scale), friendly to run on PC.
library(tibble)
library(igraph)
library(readr)
library(topGO)
library(GSEABase)
library(dplyr)
library(future)
library(Matrix)
future::plan(multicore)
options(future.globals.maxSize= 5368709120)
# load ground-truth dictionaries
m_dict. <- read_rds('msites_BP_all.rds')
load('geneID2GO.RData')
e_dict. <- read_rds('e_dict.rds')
have  <- names(geneID2GO_BP) # entrez ids with GO annotation
paste0('our work starts at ', Sys.time())
# however we need some outliers: entrez ids without GO annotation
# this step is used to prevent unprecendented result... when run topGO

edict_have <- e_dict.$exp[lengths(e_dict.$GO_slim) > 0]
edict_none <- e_dict.$exp[lengths(e_dict.$GO_slim) == 0]

geneids <- c(sample(edict_have, size = 8000), sample(edict_none, size = 2000))
geneids <- sample(geneids, size = length(geneids))
msites <- sample(m_dict.$modName, size = 7000)

####
### 1. toy data generator #######
gmm <- static.power.law.game(no.of.nodes = 7000, no.of.edges = 100000, 
                             exponent.in = -1, exponent.out = 3, multiple = FALSE, 
                             finite.size.correction = TRUE)
gee <- static.power.law.game(no.of.nodes = 10000, no.of.edges = 400000, 
                             exponent.in = -1, exponent.out = 3, multiple = FALSE, 
                             finite.size.correction = TRUE)
gme <- sample_bipartite(n1 = 7000, n2 = 10000, type = 'gnm', m = 300000, directed = FALSE)
#give names to the nodes of graphs
V(gme)$name <- c(msites,geneids)
V(gee)$name <- geneids
V(gmm)$name <- msites
gu <- igraph::union(gme,gee, gmm)
paste0('get GU at ', Sys.time())



## 2. ranwalker ##########
# function for random walk



M <- Matrix(as_adj(gu)) # with name
writeMM(M, 'sparseM0911.txt')
paste(sum(M), 'edges')

print('getM')
pc0 <- apply(M, 2, function(c){
    if (sum(c) == 0){
        c
    }
    else{ c/sum(c)}
})
pc0 <- Matrix(pc0)
print('pc0 ok')
rw <- function(thres = 10^(-10), back, p0m){
    p0 <- p0m
    p <- p0m
    i <- 1
    plast <- Matrix(0, ncol = ncol(p0), nrow = nrow(p0))
    while (max(abs(p-plast)) > thres){
        plast <- p
        p <- (1-back) * (p0m %*% p) + back * p0 
        print(paste0('result did not converge until ', i, ' ', Sys.time()))
        i <- i +1
    }
    p
}

rw_mat <- rw(back = 0.7, p0m = pc0)

paste0('random walk finished at ', Sys.time())

# find the likely nodes to be stepped.
trace_nodes <- function(m, top){
    outm <- apply(m,2, function(v){
        names(sort(v)[1:top])
    })
    df <- data.frame(node = colnames(outm))
    #browser()
    df$rw10 <- split(outm, col(outm))
    as_tibble(df)
}

rw_trace <- trace_nodes(rw_mat, 10)
print('made rw_trace')

rm(rw_mat, pc0, M)






### 3. eval_raw ######
# first only retain the selected raw dict with GO slim! (reduce the computation load)
e_have <- e_dict.[match(geneids, e_dict.$exp), ] %>% dplyr::filter(lengths(GO_slim) >=1)
m_have <- m_dict.[match(msites, m_dict.$modName),] %>% dplyr::filter(lengths(GO_slim) >=1)
rm(e_dict., m_dict.)



print('start make eval')

make_eval_raw <- function(g, md, ed, rwt){
    md0 <- md %>% dplyr::rename(node = modName, GO_exact = GO_exact_BP) %>% dplyr::select(node, GO_exact, GO_slim)
    ed0 <- ed %>% dplyr::rename(node = exp) %>% dplyr::select(node, GO_exact, GO_slim)
    raw <- rbind(md0, ed0)
    step1nb <- future({lapply(sapply(X = raw$node, FUN = igraph::ego,
                                     order = 1, mindist = 1, graph = g), as_ids)})
    step2nb <- future({lapply(sapply(X = raw$node, FUN = igraph::ego,
                                     order = 2, mindist = 1, graph = g), as_ids)})
    print('made step1nb and step2nb futures'); Sys.time()
    
    
    degf <- future({igraph::degree(g)})
    betf <- future({estimate_betweenness(g, directed = FALSE, cutoff = 1000)})
    
    
    evf <- future({evcent(g)$vector})
    prf <- future({page_rank(g)$vector})
    clf <- future({closeness(g)})
    pwf <- future({power_centrality(g)})
    
    print('made evf prf clf, pwf futures')
    
    deg <- value(degf)
    print(paste('get deg at '), Sys.time())
    bet <- value(betf)
    print(paste('get bet at '), Sys.time())
    ev <- value(evf)
    print(paste('get ev at '), Sys.time())
    pr <- value(prf)
    print(paste('get pr at '), Sys.time())
    cl <- value(clf)
    print(paste('get cl at '), Sys.time())
    pw <- value(pwf)
    print(paste('get pw at '), Sys.time())
    
    
    raw$step1nb <- value(step1nb)
    print(paste0('step1nb finished at ', Sys.time()))
    raw$step2nb <- value(step2nb)
    print(paste0('step2nb finished at ', Sys.time()))
    raw$rw10 <- rwt$rw10[match(raw$node, rwt$node)]
    print(paste0('rw appended at ', Sys.time()))
    raw$deg <- deg[match(raw$node,names(deg))]
    raw$between <- bet[match(raw$node, names(bet))]
    raw$evcent <- ev[match(raw$node,names(ev))]
    print(paste0('deg and bet and evcent finished at ', Sys.time()))
    raw$page_rank <- pr[match(raw$node,names(pr))]
    print(paste0('page_rank finished at ', Sys.time()))
    raw$closeness <- cl[match(raw$node,names(cl))]
    print(paste0('closeness finished at ', Sys.time()))
    raw$power <- pw[match(raw$node,names(pw))]
    print(paste0('power finished at ', Sys.time()))
    raw
}

eval_raw <- make_eval_raw(gu, md = m_have, ed = e_have, rwt = rw_trace)
write_rds(eval_raw, 'eraw-test1.rds')
paste0('all finished at ', Sys.time())


##### 4. med filters ####



# subfunc-1: transfer a mixed node vec to their matching geneids
any2exp <- function(mixchar, yes = have, md = m_have, ed = e_have){
    mids <- grep('m6A', mixchar)
    meth_char <- mixchar[mids]
    exp_char <- mixchar[-mids]
    
    if (length(mids) == 0){
        exp_char <- mixchar
    }
    
    meth2exp <- md$gene_id[match(meth_char,md$modName)]
    m2e <- meth2exp[!is.na(meth2exp)]
    unique(c(m2e, exp_char))
} # input char, out char



# 1. filter entire eval-raw to methyl-based evidence only
meth_filter <- function(eval, yes = have){
    
    eval0 <- eval %>% dplyr::select(-node, -power, - evcent, -GO_exact, -GO_slim,
                                    -page_rank, -closeness, -deg, -betweeness)
    
    only_meth <- as_tibble(apply(eval0, c(1,2), function(cl){
        char_vec <- cl[[1]]
        grep( 'm6A', char_vec, value = TRUE)
    }))
    
    m2et_raw <- as_tibble(apply(only_meth, c(1,2), function(cl){
        char <- cl[[1]]
        any2exp(char)}))
    
    base_eval <- eval %>% dplyr::select(node, GO_slim)
    out0 <- bind_cols(base_eval, m2et_raw)
    out0
} # eval alike



# 2. filter entire eval-meth to expressions only
exp_filter <- function(eval, yes = have){
    
    eval0 <- eval %>% dplyr::select(-node, -power, - evcent, -GO_exact, -GO_slim,
                                    -page_rank, -closeness, -deg, -betweeness)
    
    only_exp <- as_tibble(apply(eval0, c(1,2), function(cl){
        char_vec <- cl[[1]]
        mids <- grep( 'm6A', char_vec)
        char_vec[-mids]
    }))
    
    #browser()
    e2et_raw <- as_tibble(apply(only_exp, c(1,2), function(cl){
        char <- cl[[1]]
        any2exp(char)}))
    
    
    #browser()
    base_eval <- eval %>% dplyr::select(node, GO_slim)
    out0 <- bind_cols(base_eval, e2et_raw)
    out0
} # out = eval alike




# 3. double-kill 
dk_filter <- function(eval, yes = have){
    eval0 <- eval %>% dplyr::select(-node, -power, - evcent, -GO_exact, -GO_slim,
                                    -page_rank, -closeness, -deg, -betweeness)
    e2et_raw <- as_tibble(apply(eval0, c(1,2), function(cl){
        char <- cl[[1]]
        any2exp(char)}))
    base_eval <- eval %>% dplyr::select(node, GO_slim)
    out0 <- bind_cols(base_eval, e2et_raw)
} 



# operaitons : real data ##
eval_meth <- eval_raw %>% dplyr::filter(
    node %in% grep('m6A', eval_raw$node, value = TRUE)
)
### input eval_meth and output 3 tibbles for evaluaitons.
mvg <- meth_filter(eval= eval_meth, yes = have)
evg <- exp_filter(eval = eval_meth, yes = have)
dvg <- dk_filter(eval = eval_meth, yes = have)


# opreations : positive control
eval_exps <- eval_raw %>% dplyr::filter(
   ! node %in% grep('m6A', eval_raw$node, value = TRUE)
)
### basically to select the exps using the same principle
mvgpc <- meth_filter(eval= eval_exps, yes = have)
evgpc <- exp_filter(eval = eval_exps, yes = have)
dvgpc <- dk_filter(eval = eval_exps, yes = have)


###  5. GO get ############

# using topGO package
test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'fisher test')
annGO <- function(gl, geneID2GO. = geneID2GO_BP ){
    f <- list()
    for (i in seq_along(gl)){
        f[[i]] <- future({
            ids <- gl[[i]]
            geneList <- factor(as.integer(have %in% ids))
            if (length(levels(geneList)) == 1){character()}
            else {names(geneList) <- have
            GOdata <- new("topGOdata", ontology = 'BP', 
                          allGenes = geneList, annot = annFUN.gene2GO,
                          gene2GO = geneID2GO.)
            resultFisher <- getSigGroups(GOdata, test.stat)
            t <- sort(score(resultFisher))
            t <- t[t < 0.1] 
            t}
        }) 
    }
    f
    v <- lapply(f, FUN = value)
    v
} # tested robust to null predictors

# input cvg evg and dvg to make table with GO terms
xvg2GO <- function(t){
    base_eval <- t %>% dplyr::select(node, GO_slim)
    exps_eval <- t %>% dplyr::select(-node, -GO_slim)
    base_eval$step1nb <- annGO(gl = exps_eval$step1nb)
    base_eval$step2nb <- annGO(gl = exps_eval$step2nb)
    base_eval$rw10 <- annGO(gl = exps_eval$rw10)
    base_eval
}

# cut the result by p_val
p_vals = sort(c(10^(-1), 10^(-1.5), 10^(-2), 10^(-2.5), 
                10^(-3), 10^(-3.5), 10^(-4), 10^(-4.5), 
                10^(-5), 10^(-5.5), 10^(-6), 10^(-6.5)))

cutpval <- function(t, p){
    GO_eval <- t %>% dplyr::select(-node, -GO_slim)
    m  <- apply(GO_eval, c(1,2), function(a){
        a <- a[[1]]
        names(a)[a < p]
    })
    m
}

# a wrapper for all funcitons in this section
sel_go  <- function(vg, p = p_vals){
    o <- Map(cutpval, p = p_vals, 
        MoreArgs = list(t = xvg2GO(t = vg)))
    names(o) <- as.character(p)
    o
}

## operations ##

ev_go <- sel_go(evg)
mv_go <- sel_go(mvg)
dv_go <- sel_go(dvg)

### positive control

ev_go_pc <- sel_go(evgpc)
mv_go_pc <- sel_go(mvgpc)
dv_go_pc <- sel_go(dvgpc)



#### 6. slim ##########

fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)

# identify any matching go slim terms that have count bigger than 1
enslim <- function(got){
    gl <- split(got, col(got))
    names(gl) <- colnames(got)
    golist <- unlist(gl, recursive = F)
    
    ansf <- list()
    
    for (i in seq_along(golist)){
        ansf[[i]] <- future({
            gos <- golist[[i]]
            
            
            if (length(gos)>0){
                mygo <- GOCollection(gos)
                temp <- goSlim(mygo, slim, 'BP')
                if (!all(unique(temp$Count) == 0)){
                    temp <- temp[temp$Count>0,]
                    temp <- temp[with(temp, order(-Count, -Percent)),]
                    rownames(temp)
                }
                else { character()}
            }
            else {character()}
        })
    }
    
    
    ans <- lapply(ansf, value)
    out <- matrix(data = ans, ncol = ncol(got), byrow = F)
    colnames(out) <- colnames(got)
    out
}

# apply the same operation as enslim, but the object is a list
enslim.l <- function(l){
    out <- list()
    for (i in seq_along(l)){
        df <- l[[i]]
        if (length(df) == 0){out[[i]] <- character()}
        else{out[[i]] <- enslim(df) }
    }
    names(out) <- names(l)
    out
}

#
sel.slim <- function(l,seln){
    l1 <- lapply(l, function(mat){
        attr0 <- attributes(mat)
        ans <- lapply(mat, function (chr) {
            t <- chr[seq(seln)]
            t <- t[!is.na(t)]
        })
        attributes(ans) <- attr0
        ans
    })
    l1
}

nt_vec <- seq(10)
sel.slim2 <- function(lm, n, p = p_vals){
    nl  <- Map(sel.slim, seln = n, 
               MoreArgs = list(l = lm))
    #browser()
    sl <- unlist(nl, recursive = F)
    out <- matrix(sl, ncol = length(n), byrow = F)
    colnames(out) <- as.character(n)
    rownames(out) <- as.character(p)
    out
}

## operations ###

pa2_mv <- sel.slim2(lm = enslim.l(mv_go), n = nt_vec)
pa2_ev <- sel.slim2(lm = enslim.l(ev_go), n = nt_vec)
pa2_dv <- sel.slim2(lm = enslim.l(dv_go), n = nt_vec)


pa2_mv_pc <- sel.slim2(lm = enslim.l(mv_go_pc), n = nt_vec)
pa2_ev_pc <- sel.slim2(lm = enslim.l(mv_go_pc), n = nt_vec)
pa2_dv_pc <- sel.slim2(lm = enslim.l(mv_go_pc), n = nt_vec)



#### 7. performance #########

# make a matrix of matrix for ground truth and p2(pa2)
eval_align_pred <- function(p2, t){
    ind <- (apply(p2, c(1,2), lengths) >0 )
    # make a filter for rows
    stay_row <- apply(ind, 1, function(lg){all(lg)}) 
    p2_new <- p2[stay_row,] # the basis of evaluation
    
    t0 <- t %>% dplyr::select(node, GO_slim) # the ground_truth for evulation
    # to each element in the matrix
    apply(p2_new, c(1,2), function(a){
        a <- a[[1]]
        a <- as.tibble(a)
        #browser()
        as.tibble(cbind(t0, a))
    })
}

# input re, pattern, out list of 3 elements: alldat, 
# precision and specificity
eval_pat <- function(re_align, pattern){
    rt_names <- c('node', 'GO_slim', pattern)
    alldat <- apply(re_align, c(1,2), function(m){
        m <- m[[1]]
        m <- m[,rt_names]
        m$inter <- mapply(function(a,b){intersect(a,b)}, 
                          a = m$GO_slim, b = m[,pattern][[1]])
        
        l <- lengths(m[, pattern][[1]])
        stay_row <- (l > 0)
        m[stay_row,]
    })
    #browser()
    precision <- apply(alldat, c(1,2), function(m){
        m <- m[[1]]
        #browser()
        sum(lengths(m$inter)) / sum(lengths(m$GO_slim))
    })
    specificity <- apply(alldat, c(1,2), function(m){
        m <- m[[1]]
        sum(lengths(m$inter)) / sum(lengths(m[, pattern][[1]]))
    })
    out <- list(alldat, precision, specificity)
    names(out) <- c('alldata', 'precision', 'specificity')
    out
}

comp_pat <- function(pa2, go){
    re <- eval_align_pred(p2 = pa2, t = go)
    re_step1nb <- eval_pat(re_align = re, pattern = 'step1nb')
    re_step2nb <- eval_pat(re_align = re, pattern = 'step2nb')
    re_rw10 <- eval_pat(re_align = re, pattern = 'rw10')
    o <- list(re_step1nb, re_step2nb, re_rw10)
    names(o) <- c('step1nb','step2nb','rw10')
    o
}
# operations ##

m_comp <- comp_pat(pa2 = pa2_mv, go = mv_go)
e_comp <- comp_pat(pa2 = pa2_ev, go = ev_go)
d_comp <- comp_pat(pa2 = pa2_dv, go = dv_go)


m_comp_pc <- comp_pat(pa2 = pa2_mv_pc, go = mv_go_pc)
e_comp_pc <- comp_pat(pa2 = pa2_ev_pc, go = ev_go_pc)
d_comp_pc <- comp_pat(pa2 = pa2_dv_pc, go = dv_go_pc)





