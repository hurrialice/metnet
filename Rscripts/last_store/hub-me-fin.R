## now we presumed that we have finished network construction 
library(igraph)
library(readr)
library(tibble)
library(dplyr)
library(topGO)
library(GSEABase)
library(BiRewire)
load('geneID2GO.RData')
adj_me <- read_rds('adj-me0824-85.rds')
m_dict. <- read_rds('msites_BP_fin.rds') # note the first column name is modName
have <- names(geneID2GO_BP)


############# subfunctions ###########################
graph2dict <- function(g, m_dict = m_dict.){
    m_dict <- m_dict[match(as_ids(V(g))[!V(g)$type],m_dict$modName),] #rearrange to match to exact site
    m_dict$nb <- mapply(as_ids, mapply(neighbors, v = m_dict$modName, MoreArgs = list(graph = g)))
    m_dict$hs <- hub_score(g)$vector[match(m_dict$modName, names(hub_score(g)$vector))]
    ifhave <- mapply(function(cv){any(cv %in% have)}, cv = m_dict$nb)
    eval <- m_dict[ifhave,] # no number selection
    eval
}
test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'fisher test')
i <- 0
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
sel.slim <- function(m, nt = nt_vec, selGO){
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
nt_vec <- seq(10)
el <- function(l, eval. ){
    e0 <- eval. %>% dplyr::select(modName, GO_slim)
    l <- l[[1]]
    e0$predict <- l[match(e0$modName, names(l))] #aux!
    e0$hit <- mapply(function(a,b) intersect(a,b), a = e0$GO_slim, b = e0$predict)
    e0 <- e0[lengths(e0$GO_slim)>0 & lengths(e0$predict)>0,]
    pre <- sum(lengths(e0$hit))/sum(lengths(e0$GO_slim))
    spe <- sum(lengths(e0$hit)) / sum(lengths(e0$predict))
    out = list(e0, pre, spe)
    names(out) <- c('alldat', 'precision', 'specificity')
    out
}

#################### main function #######################
# main function -1: make eval (dataframe)
make_eval <- function(g){
    eval <- graph2dict(g)
    eval$GO_predict <- mapply(annGO, ids = eval$nb)
    eval
}
# main function -2: make pa2
make_pa2 <- function(pval = p_vals, nts = nt_vec, eval){
    sel.GO <- mapply(cutpval, p = pval, MoreArgs = list(lg = eval$GO_predict))
    colnames(sel.GO) <- as.character(pval)
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

######## operations ################
gme <- graph.incidence(adj_me)
adj_ran <- birewire.rewire.bipartite(adj_me); print('2')
gran <- graph.incidence(adj_ran)
eval_true <- make_eval(gme)
eval_ran <- make_eval(gran)
pa2_true <- make_pa2(eval = eval_true)
pa2_ran <- make_pa2(eval = eval_ran)
pa2all_true <- apply(pa2_true, c(1,2), el, eval. = eval_true)
pa2all_ran <- apply(pa2_ran, c(1,2), el, eval. = eval_ran)

save(eval_true, eval_ran, 
     pa2all_ran, pa2all_true, 
     pa2_true, pa2_ran, file = 'hubme-tr-0824.RData')

# evaluate how the performance superior to random
comp <- better()
write_rds(comp, 'compare-hubme-0824.rds')

