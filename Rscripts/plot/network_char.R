library(readr)
library(multtest)
library(WGCNA)
library(preprocessCore)
test_m <- read_rds('test-m0810.rds')
test_exp <- read_rds('test-exp0810.rds')
mininfo <- function(m){m[1:10,1:10]}
library(GenomicFeatures)
m_info <- read_rds('msites_BP_fin.rds')
bin_rpkm <- read_rds('bin_rpkm.rds')
bin <- read_rds('bin.rds')


calm <- function(bin. = bin, bin_rpkm. = bin_rpkm, par, test_m. = test_m){
    msite_all <-  bin.$modName
    ip_names <- sort(grep('_ip', colnames(mcols(bin_rpkm.)), value = TRUE))
    input_names <- sort(grep('_input', colnames(mcols(bin_rpkm.)), value = TRUE))
    stopifnot(identical(gsub('_ip','', ip_names), gsub('_input','', input_names)))
    ip <- as.matrix(mcols(bin_rpkm.)[,ip_names])
    input <- as.matrix(mcols(bin_rpkm.)[,input_names])
    M_vals <- t(log2(ip+0.01) - (1-par)*log2(input+0.01) ) #### here to modify this formular
    rownames(M_vals) <- gsub('_ip','', rownames(M_vals))
    colnames(M_vals) <- msite_all
    msites <- colnames(test_m.)
    conds <- rownames(test_m.)
    selMval <- M_vals[conds, msites]
    selMval
}
par_vec <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1)
mlist <- Map(calm, par = par_vec) 

m <- t(normalize.quantiles(t(test_m)))
colnames(m) <- colnames(test_m)
rownames(m) <- rownames(test_m)
test_m <- m
write_rds(test_m, 'test_m0814.rds')


# make self-correlation 
makeself <- function(test_m){
    geneids <- m_info$gene_id
    have <- colnames(test_exp) %in% geneids
    self <- test_exp[, which(have)]
    a <- as.character(m_info$gene_id[match(colnames(test_m), m_info$modName)])
    self <- self[,a]
    self
}
self <- makeself(test_m)

# make the plot, no background
palette(rainbow(20))
cor_true <- diag(cor(mlist[[1]], self, method = 'spearman'))


png(filename = 'density-dynamics-par.png',width = 700, height = 550, units = 'px')
plot(density(cor_true), 
     main = 'density dynamics with respect to par', xlab = 'spearman correlation')
text(0.7, 1.7, 'par = 0')
abline(v = 0, lty = 2)
for (i in (seq_along(mlist)-1)){
    c <- diag(cor(mlist[[i+1]], self, method = 'spearman'))
    lines(density(c), col = i+1)
    text(-0.9, 1.5-0.05*i, sprintf(paste0('par = ',par_vec[i+1])), col = i+1, pos = 4)
    rm(c)
}
dev.off()

######
# permutate conditions
m <- mlist[[1]]
png(file = "conds-permutate.png", width = 700, height = 550, units = "px")
plot(c(-1.2,1.2),c(0,3),type = "n", 
     xlab = 'Spearman Correlation', 
     ylab = 'density', main = 'spearman correlation for condition permutation')
for (i in 1:100){
    rownames(m) <- sample(rownames(m), size = nrow(m))
    m <- m[order(rownames(m)), ]
    cor_me <- diag(cor(self, m, method = 'spearman'))
    lines(density(cor_me), col = 'grey')
}
abline(v = 0, lwd = 1, lty = 2)
lines(density(cor_true), col = 'red', lwd = 2)
dev.off()

# permutate sites
png(file = "sites-permutate.png", width = 700, height = 550, units = "px")
plot(c(-1.2,1.2),c(0,2),type = "n", 
     xlab = 'Spearman Correlation', 
     ylab = 'density', main = 'spearman correlation for sites permutation')
others <- test_exp
set.seed(1)
for (i in 1:100){
    inds <- sample(seq(ncol(others)), size = 2597)
    a <- others[, inds]
    cor_others <- cor(test_m, a, method = 'spearman')
    lines(density(cor_others), col = 'grey')
} 
abline(v = 0, lwd = 1, lty = 2)
lines(density(cor_true), col = 'red', lwd = 2)
dev.off()
print('done')


# make the network co-exp and exp-meth
adjmake0 <- function(x, y, quant, pcut){
    cor <- cor(x,y,method = 'spearman')
    rawp <- corPvalueFisher(cor, 40)
    mt <- mt.rawp2adjp(rawp, proc = 'Bonferroni')
    adj <- mt$adjp[,2]
    adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
    cor <- abs(cor)
    scc_cut <- quantile(cor, quant)
    cor[adjp > pcut | cor < scc_cut] <- 0
    cor[cor > 0] <- 1
    diag(cor) <- 0
    cor
}
adjmake <- function(x, y, quant, pcut){
    cor <- cor(x,y,method = 'spearman')
    rawp <- corPvalueFisher(cor, 40)
    mt <- mt.rawp2adjp(rawp, proc = 'Bonferroni')
    adj <- mt$adjp[,2]
    adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
    cor <- abs(cor)
    scc_cut <- quantile(cor, quant)
    cor[adjp > pcut | cor < scc_cut] <- 0
    cor[cor > 0] <- 1
    cor
}
adj_exp <- adjmake0(x = test_exp, y = test_exp, quant = 0.98, pcut = 0.01)
adj_me <- adjmake(x = test_m, y = test_exp, quant = 0.80, pcut = 0.05)
write_rds(adj_exp, 'adj-exp-stren.rds')
write_rds(adj_me, 'adj-me-0815.rds')

# degree distribution
library(igraph)
gee <- graph.adjacency(adj_exp, mode = 'undirected')
gme <- graph.incidence(adj_me)
# 1. degree (distribution and averag

d <- sort(degree(gee), decreasing = TRUE)
d <- d[d>0]
d2 <- sort(degree(gme), decreasing = TRUE)
d2 <- d2[d2 > 0]
write_rds(d, 'most-degree-nodes-exp.rds')
write_rds(d2, 'most-degree-meth-exp.rds')


# co-expression network
png(filename = 'degree-distribution-stren-coexp.png', width = 700, height = 550, units = 'px')
dd <- degree_distribution(gee)

plot.new()
x <- seq(length(dd))
y <- dd + 0.001
plot(y, xlim=c(0.001,100), ylim=c(0.001, 1),log="xy", xlab="degree",ylab="frequency")
ticks <- seq(0, 2.5, by = 0.5)
labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(0.01, 0.1, 1, 10, 100, 1000), labels=labels)





plot(log10(seq(length(dd))),log10(dd+0.001), 
     xlab = 'log10(degree)', ylab = 'log10(frequency)', 
     main = 'degree distribution of co-expression network', pch = 20, col = 'red')
avg.degree <- mean(degree(gee, V(gee)))
abline(v = avg.degree, lwd = 1, lty = 2)
dev.off()
# methylation-expression network
png(filename = 'degree-distribution-meth2exp.png', width = 700, height = 550, units = 'px')
ddme <- degree_distribution(gme)
plot(log10(seq(length(ddme))),log10(ddme+0.001), 
     xlab = 'log10(degree)', ylab = 'log10(frequency)', 
     main = 'degree distribution of methylation-expression network', pch = 20, col = 'red')
avg.degree <- mean(degree(gme, V(gme)))
abline(v = avg.degree, lwd = 1, lty = 2)
dev.off()


# 2. degree distribution of average clustering coefficient
t <- 100*length(E(gee))
gr_exp <- rewire(gee, keeping_degseq(niter = t))
gr_exp
cal_acc <- function(g){
    cc <- transitivity(g, type = 'local', isolates = 'zero')
    d <- degree(g, loops = FALSE)
    df <- data.frame(v = names(d), deg = unname(d), lcc = cc)
    acc <- unlist(lapply(split(df, df$deg), function(sdf){
        mean(sdf$lcc)
    }))
    acc
}

ee_acc <- cal_acc(gee)
print('ee_acc')
ree_acc <- cal_acc(gr_exp)
print('ree_acc')
gl.cc <- transitivity(gee, type = 'global')
print(gl.cc)
gl.rcc <- transitivity(gr_exp, type = 'global')
print(gl.rcc)


png(filename = 'avg-lcc-distribution.png', width = 700, height = 550, units = 'px')
plot(as.numeric(names(ee_acc)),ee_acc, pch = 20, 
     xlab = 'degree', ylab = 'Average clustering coefficient',
     main = 'Average clustering coefficient distribution')
points(as.numeric(names(ree_acc)),ree_acc, pch = 20, col = 'red')
abline(h = gl.cc, lwd = 1, lty = 2)
abline(h = gl.rcc, lwd = 1, lty = 2, col = 'red')
dev.off()
