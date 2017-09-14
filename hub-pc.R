## now we presumed that we have finished network construction 
library(igraph)
library(readr)
library(tibble)
library(dplyr)
library(topGO)
library(GSEABase)
library(BiocParallel)
load('geneID2GO.RData')
adj_exp <- read_rds('adj-coexp0824.rds') #same par as NAR paper
geneids <- colnames(adj_exp)
e_dict. <- read_rds('e_dict.rds')
e_dict. <- e_dict.[e_dict.$exp %in% geneids,] 
have <- names(geneID2GO_BP)
# tested: rewired plot gives same idseq as original.


########## sub-function ##############################
# from graph to a small dict containing the neighbors (as char)
graph2dict <- function(g, e_dict = e_dict.){
    idseq <- as_ids(V(g))
    e_dict <- e_dict[match(idseq,e_dict$exp),]
    e_dict$nb <- mapply(as_ids, 
                        mapply(neighbors, v = e_dict$exp, 
                               MoreArgs = list(graph = g)))
    e_dict$hs <- hub_score(g)$vector[match(e_dict$exp, 
                                           names(hub_score(g)$vector))]
    ifhave <- mapply(function(cv){any(cv %in% have)}, cv = e_dict$nb)
    eval <- e_dict[ifhave,]
    eval <- eval[with(eval, order(-hs)), ][1:1000,] #to minimize workload
    eval
}
# note here the sequence of eval is determined by its hub score, 
# independent on idseq, which is the colnames of adj_exp.

## annotate with enriched GO terms
i <- 0
test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'fisher test')
annGO <- function(ids, geneID2GO. = geneID2GO_BP){
    i <<- i+1; print(i)
    geneList <- factor(as.integer(have %in% ids))
    names(geneList) <- have
    GOdata <- new("topGOdata", ontology = 'BP', 
                allGenes = geneList, annot = annFUN.gene2GO,
                gene2GO = geneID2GO.)
    resultFisher <- getSigGroups(GOdata, test.stat)
    t <- sort(score(resultFisher))
    t <- t[t < 0.1]
    t
}

p_vals = sort(c(10^(-1), 10^(-1.5), 10^(-2), 10^(-2.5), 
                10^(-3), 10^(-3.5), 10^(-4), 10^(-4.5), 
                10^(-5), 10^(-5.5), 10^(-6), 10^(-6.5)))
cutpval <- function(lg,p){
    lapply(lg,function(n){
        names(n[n < p])
    })
}

### make slim from enriched GO terms
fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
enslim <- function(go){
    apply(go, c(1,2), function(char){
        char <- unname(unlist(char, recursive = FALSE))
        if (length(char)>0){
            mygo <- GOCollection(char)
            temp <- goSlim(mygo, slim, 'BP')
            if (!all(unique(temp$Count) == 0)){
                temp <- temp[temp$Count>0,]
                temp <- temp[with(temp, order(-Count, -Percent)),]
                rownames(temp)
            }
            else {character()}
        }
        else {character()}
    })
}
nt_vec <- seq(10)
sel.slim <- function(m, nt, selGO){
    l <- split(m, col(m))
    names(l) <- colnames(m)
    l <- lapply(l, function(ll){
        names(ll) <- rownames(selGO)
        ll})
    lcut <- lapply(l, function(sl){
        lapply(sl, function(s){
            ans <- s[1:nt]
            ans[!is.na(ans)]
        })
    })
    lcut
}
el <- function(l, eval.){
    e0 <- eval. %>% dplyr::select(exp, GO_slim)
    l <- l[[1]]
    e0$predict <- l[match(e0$exp, names(l))] #aux!
    e0$hit <- mapply(function(a,b) intersect(a,b), a = e0$GO_slim, b = e0$predict)
    e0 <- e0[lengths(e0$GO_slim)>0 & lengths(e0$predict)>0,]
    pre <- sum(lengths(e0$hit))/sum(lengths(e0$GO_slim))
    spe <- sum(lengths(e0$hit)) / sum(lengths(e0$predict))
    out = list(e0, pre, spe)
    names(out) <- c('alldat', 'precision', 'specificity')
    out
}

########### main functions ##########
### true
# main function -1: make_eval
make_eval <- function(g){
    eval <- graph2dict(g)
    i <- 0
    eval$GO_predict <- mapply(annGO, ids = eval$nb)
    eval
}
# main function -2: make pa2
make_pa2 <- function(pval = p_vals, nts = nt_vec, eval){
    sel.GO <- mapply(cutpval, p = pval, MoreArgs = list(lg = eval$GO_predict))
    colnames(sel.GO) <- as.character(pval)
    rownames(sel.GO) <- eval$exp
    pa2 <- mapply(sel.slim, nt = nts, 
                  MoreArgs = list(m = enslim(sel.GO),
                                  selGO = sel.GO))
    colnames(pa2) <- as.character(nts)
    pa2
}
# main function -3: test if superior?
better <- function(r = pa2all_ran, t = pa2all_true){
    tspe <- apply(t, c(1,2), function(l3){l3[[1]]$specificity})
    tpre <- apply(t, c(1,2), function(l3){l3[[1]]$precision})
    rspe <- apply(r, c(1,2), function(l3){l3[[1]]$specificity})
    rpre <- apply(r, c(1,2), function(l3){l3[[1]]$precision})
    comp <- list(true = list(tpre, tspe),
                 random = list(rpre, rspe),
                 inc = list(p = tpre - rpre,
                            s = tspe - rspe))
    comp
}
############# operations ###############
gee <- graph.adjacency(adj_exp, mode = 'undirected')
tt <- 100*length(E(gee))
gran <- rewire(gee, keeping_degseq(niter = tt))
eval_true <- make_eval(gee)
eval_ran <- make_eval(gran)
pa2_true <- make_pa2(eval = eval_true)
pa2_ran <- make_pa2(eval = eval_ran)
pa2all_true <- apply(pa2_true, c(1,2), el, eval. = eval_true)
pa2all_ran <- apply(pa2_ran, c(1,2), el, eval. = eval_ran)

save(eval_true, eval_ran, 
     pa2all_ran, pa2all_true, 
     pa2_true, pa2_ran, file = 'hubpc-tr-0824.RData')

# evaluate how the performance superior to random
comp <- better()
write_rds(comp, 'compare-hubpc0824.rds')
