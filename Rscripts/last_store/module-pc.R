# further questions # 
# 1. how well is the module assignment agrees with actual clustering result?
# 2. much can be 'definately' clustered into clusters?
library(readr)
library(igraph)
library(GSEABase)
library(topGO)
library(dplyr)
library(tibble)
library(WGCNA)
library(tidyr)
test_exp <- read_rds('test-exp-0824.rds')
geneids <- colnames(test_exp)
load('geneID2GO.RData')
e_dict. <- read_rds('e_dict.rds')
e_dict. <- e_dict.[e_dict.$exp %in% geneids,] 
have <- names(geneID2GO_BP) #ubiquitous
# use the inherent sequence for e_dict.

################## sub-functions ################
impexp <- function (path) {
    s <- readLines(path)
    spaceline <- which(s == "")
    stopifnot(length(spaceline) == 2)
    needed <- s[(spaceline[2] + 1L):length(s)]
    file <- tempfile()
    writeLines(needed, file)
    read_tsv(file)
}
split_node <- function(char) {
    a <- strsplit(char, ", ")
    lapply(a, as.character)
}
mcode_filter <- function(oric){
    l <- oric$`Node IDs`
    retain <- sapply(l, function(chars){sum(chars %in% have) > 2})
    oric <- oric[retain,]
    oric$Cluster <- seq(nrow(oric))
    oric
}

cal_ME <- function(test_exp. = test_exp, e2c){
    mod_exp <- test_exp.[,e2c$exp]
    MEs <- moduleEigengenes(mod_exp, colors = e2c$cluster)$eigengenes
    colnames(MEs) <- gsub('ME','',colnames(MEs))
    MEs
    }
ass_mod <- function(cME){
    apply(cME, 1, function(r){
        colnames(cME)[which(r == max(r))[1]] # here only one assignment is made
    })
} # output a named char vector
del_self <- function(one, c){
    if (one %in% c){c[!c %in% one]}
    else{c}
}
oric_alike <- function(e2c){
    l <- split(e2c, e2c$cluster)
    cids <- lapply(l, function(df){unname(unlist(df$exp))})
    cids
}
ana_modass <- function(e2c, cME, gm, e_dict = e_dict.){
    # first derieve information from cME
    # proved: range and sd show strong positive correlation!
    e <- tibble(exp = names(gm), module = gm,
                    range = apply(cME, 1, function(nvec){max(nvec) - min(nvec)}), 
                    sd = apply(cME, 1, sd)) 
    e <- e[with(e, order(-range, -sd)),]
    oric <- oric_alike(e2c)
    e$cids <- oric[match(e$module,names(oric))]
    # does this assigned cluster contain the exp itself?
    # used for performance interpretation (later make a cumulative plot)
    e$ifcontain <- mapply(function(a,b) a %in% b, a = e$exp, b = e$cids)
    e$newc <- mapply(del_self, one = e$exp, c = e$cids)
    e <- left_join(e, e_dict)
    e[1:1000,] # choose a small subset for testing (be deleted when on server)
} 
# gives a full version tibble for topGO, 
#similar to 'eval' object in hub-based predictions.

# GO enrichment analysis
test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'fisher test')
annGO <- function(ids, geneID2GO. = geneID2GO_BP){
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


# make slim
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
nt_vec = seq(10)
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


####################### middle function ##########################
make_exp2cluster <- function(mcode_path){
    o <- impexp(mcode_path)
    o$`Node IDs` <- split_node(o$`Node IDs`)
    o <- mcode_filter(o)
    exp2cluster <- o[,c(1,5)] %>% unnest %>% as.data.frame()
    colnames(exp2cluster) <-c('cluster','exp')
    exp2cluster
}
ran_exp2cluster <- function(ori_e){
    new <- ori_e
    cpool <- unique(new$cluster)
    new$cluster <- mapply(function(c){sample(cpool[!cpool %in% c], size = 1)}, 
                          c = ori_e$cluster)
    new
}


########### main function ###########
# main function -3: make eval
# eval ubi-def: all-in-one information quiry set
make_eval <- function(e2c, cME, gm){
    eval <- ana_modass(e2c, cME, gm)
    eval$GO_predict <- mapply(annGO, ids = eval$newc)
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
######################   operations   ###############
e2c_true <- make_exp2cluster('mcode0825.txt')
e2c_ran <- ran_exp2cluster(e2c_true)

ME_true <- cal_ME(e2c = e2c_true)
ME_ran <- cal_ME(e2c = e2c_ran)

cME_true <- cor(test_exp, ME_true, method = 'spearman')
cME_ran <- cor(test_exp, ME_ran, method = 'spearman')

gm_true <- ass_mod(cME = cME_true)
gm_ran <- ass_mod(cME = cME_ran)


eval_true <- make_eval(e2c = e2c_true, cME = cME_true, gm = gm_true)
eval_ran <- make_eval(e2c = e2c_ran, cME = cME_ran, gm = gm_ran)

pa2_true <-  make_pa2(eval = eval_true)
pa2_ran <- make_pa2(eval = eval_ran)

pa2all_true <- apply(pa2_true, c(1,2), el, eval. = eval_true)
pa2all_ran <- apply(pa2_ran, c(1,2), el, eval. = eval_ran)

save(eval_true, eval_ran, 
     pa2all_ran, pa2all_true, 
     pa2_true, pa2_ran, file = 'modulepc-tr-0824.RData')

comp <- better()
write_rds(comp, 'compare-modulepc-0824.rds')

