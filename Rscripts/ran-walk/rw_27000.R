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
gee <- static.power.law.game(no.of.nodes = 20000, no.of.edges = 1600000, 
                             exponent.in = -1, exponent.out = 3, multiple = FALSE, 
                             finite.size.correction = TRUE)
gme <- sample_bipartite(n1 = 7000, n2 = 20000, type = 'gnm', m = 300000, directed = FALSE)
#give names to the nodes of graphs
V(gme)$name <- c(msites,geneids)
V(gee)$name <- geneids
V(gmm)$name <- msites
gu <- igraph::union(gme,gee, gmm)
paste0('get GU at ', Sys.time())
##########################


## 2. make random walks ##########
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
### 3. eval_raw object ######
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
    #pwf <- future({power_centrality(g)})
    
    print('made evf prf clf, pwf futures')
    
    deg <- value(degf)
    print(paste('get deg',Sys.time()))
    bet <- value(betf)
    print(paste('get bet',Sys.time()))
    ev <- value(evf)
    print(paste('get ev',Sys.time()))
    pr <- value(prf)
    print(paste('get pr',Sys.time()))
    cl <- value(clf)
    print(paste('get cl',Sys.time()))
    #pw <- value(pwf)
    #print('get pw'); Sys.time()
    
    
    
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
    raw
}

eval_raw <- make_eval_raw(gu, md = m_have, ed = e_have, rwt = rw_trace)
write_rds(eval_raw, 'eraw-test1.rds')
paste0('all finished at ', Sys.time())