library(readr)
library(gplots)
eval_CC <- read_rds('CCcompare.rds')
eval_MF <- read_rds('MFcompare.rds')

inc.CC.pre <- eval_CC$improve[[1]]
inc.CC.spe <- eval_MF$improve[[2]]

CC_pre_heat <- heatmap.2(inc.CC.pre, Rowv = FALSE)
CC_spe_heat <- heatmap.2(inc.CC.spe, Rowv = FALSE)

library(scatterplot3d)
pop <- function(nl, rrm, nrm){
    lapply(nl, function(l){
        lapply(l, function(df){
            df[-rrm, -nrm]
        })
    })
}
mod.eval.CC <- pop(eval_CC, rrm = c(1,2,3), nrm = c(1,2))
mod.eval.MF <- pop(eval_MF, rrm = c(1,2,3), nrm = c(1,2))

prep43d <- function(mod_df){
    spe <- as.data.frame(as.table(mod_df))
    spe <- as.data.frame(apply(spe, 2, as.numeric))
    colnames(spe) <- c('-lg(Pval)','Slim Cutoff','Performance')
    spe$`-lg(Pval)` <- (-1)*log10(spe$`-lg(Pval)`)
    spe
}
real.pre.CC <- prep43d(mod.eval.CC$real[[1]])
real.spe.CC <- prep43d(mod.eval.CC$real[[2]])
random.pre.CC <- prep43d(mod.eval.CC$random[[1]])
random.spe.CC <- prep43d(mod.eval.CC$random[[2]])

real.pre.MF <- prep43d(mod.eval.MF$real[[1]])
real.spe.MF <- prep43d(mod.eval.MF$real[[2]])
random.pre.MF <- prep43d(mod.eval.MF$random[[1]])
random.spe.MF <- prep43d(mod.eval.MF$random[[2]])

# make the plot for s3d_CC_spe
s3d_CC_spe <- scatterplot3d(real.spe.CC, type = 'h', color = 'blue',
                     angle = 120, scale.y = 0.8, 
                     pch = '', main = 'CC prediction precision')
s3d_CC_spe$points3d(random.spe.CC,  type = 'h', col = 'white', pch = 1)
s3d_CC_spe$points3d(random.spe.CC,  type = 'h', col = 'white', pch = 1)
s3d_CC_spe$points3d(random.spe.CC,  type = 'p', col = 'white', pch = 16)
s3d_CC_spe$points3d(real.spe.CC, type = 'p', col = 'blue', pch = 16)

s3d_CC_spe <- scatterplot3d(real.spe.CC, type = 'h', color = 'blue',
                            angle = 120, scale.y = 0.8, 
                            pch = '', main = 'CC prediction precision')

# make the plot for s3d_CC_pre
s3d_CC_pre <- scatterplot3d(real.pre.CC, type = 'h', color = 'blue',
                            angle = 120, scale.y = 0.8, 
                            pch = '', main = 'CC prediction precision')
s3d_CC_pre$points3d(random.pre.CC,  type = 'h', col = 'white', pch = 1)
s3d_CC_pre$points3d(random.pre.CC,  type = 'h', col = 'white', pch = 1)
s3d_CC_pre$points3d(random.pre.CC,  type = 'p', col = 'white', pch = 16)
s3d_CC_pre$points3d(real.pre.CC, type = 'p', col = 'blue', pch = 16)


# make the plot for s3d_MF_pre
s3d_MF_pre <- scatterplot3d(real.pre.MF, type = 'h', color = 'blue',
                            angle = 120, scale.y = 0.8, 
                            pch = '', main = 'CC prediction precision')
s3d_MF_pre$points3d(random.pre.MF,  type = 'h', col = 'white', pch = 1)
s3d_MF_pre$points3d(random.pre.MF,  type = 'h', col = 'white', pch = 1)
s3d_MF_pre$points3d(random.pre.MF,  type = 'p', col = 'white', pch = 16)
s3d_MF_pre$points3d(real.pre.MF, type = 'p', col = 'blue', pch = 16)

# make the plot for s3d_MF_spe
s3d_MF_spe <- scatterplot3d(real.spe.MF, type = 'h', color = 'blue',
                            angle = 120, scale.y = 0.8, 
                            pch = '', main = 'CC prediction precision')
s3d_MF_spe$points3d(random.spe.MF,  type = 'h', col = 'white', pch = 1)
s3d_MF_spe$points3d(random.spe.MF,  type = 'h', col = 'white', pch = 1)
s3d_MF_spe$points3d(random.spe.MF,  type = 'p', col = 'white', pch = 16)
s3d_MF_spe$points3d(real.spe.MF, type = 'p', col = 'blue', pch = 16)


######################################
############# hub based ##############
library(readr)
library(gplots)
eval_hub <- read_rds('hub-compare.rds')

inc.hub.pre <- eval_hub$improve[[1]]
inc.hub.spe <- eval_hub$improve[[2]]

hub_pre_heat <- heatmap.2(inc.hub.pre, Rowv = FALSE)
hub_spe_heat <- heatmap.2(inc.hub.spe, Rowv = FALSE)

library(scatterplot3d)
library(reshape2)
pop <- function(nl, rrm, nrm){
    lapply(nl, function(l){
        lapply(l, function(df){
            df[-rrm, -nrm]
        })
    })
}
mod.eval.hub <- pop(eval_hub, rrm = c(1), nrm = c(1))
prep43d <- function(mod_df){
    spe <- as.data.frame(melt(mod_df))
    spe <- as.data.frame(apply(spe, 2, as.numeric))
    colnames(spe) <- c('-lg(Pval)','Slim Cutoff','Performance')
    spe$`-lg(Pval)` <- (-1)*log10(spe$`-lg(Pval)`)
    spe
}
real.pre.hub <- prep43d(mod.eval.hub$real[[1]])
real.spe.hub <- prep43d(mod.eval.hub$real[[2]])
random.pre.hub <- prep43d(mod.eval.hub$random[[1]])
random.spe.hub <- prep43d(mod.eval.hub$random[[2]])
s3d_hub_spe <- scatterplot3d(real.spe.hub, type = 'h', color = 'blue',
                            angle = 120, scale.y = 0.8, 
                            pch = '', main = 'hub-based prediction specificity')
s3d_hub_spe$points3d(random.spe.hub,  type = 'h', col = 'white', pch = 1)
s3d_hub_spe$points3d(random.spe.hub,  type = 'h', col = 'white', pch = 1)
s3d_hub_spe$points3d(random.spe.hub,  type = 'p', col = 'white', pch = 16)
s3d_hub_spe$points3d(real.spe.hub, type = 'p', col = 'blue', pch = 16)
s3d_hub_spe$points3d(random.spe.hub, type = 'p', col = 'red', pch = 16)


s3d_hub_spe <- scatterplot3d(real.pre.hub, type = 'h', color = 'blue',
                             angle = 120, scale.y = 0.8, 
                             pch = '', main = 'hub-based prediction precision')
s3d_hub_spe$points3d(random.pre.hub,  type = 'h', col = 'white', pch = 1)
s3d_hub_spe$points3d(random.pre.hub,  type = 'h', col = 'white', pch = 1)
s3d_hub_spe$points3d(random.pre.hub,  type = 'p', col = 'white', pch = 16)
s3d_hub_spe$points3d(real.pre.hub, type = 'p', col = 'blue', pch = 16)
s3d_hub_spe$points3d(random.pre.hub, type = 'p', col = 'red', pch = 16)





################## about BP #########################
eval.BP <- read_rds('compare-hubme-0817.rds')
inc.BP.pre <- eval.BP$inc[[1]]
inc.BP.spe <- eval.BP$inc[[2]]
eval.BP$random <- lapply(eval.BP$random, function(m) {
    colnames(m) <- colnames(inc.BP.pre)
    m})


library(gplots)
BP_pre_heat <- heatmap.2(inc.BP.pre, Rowv = FALSE)
BP_spe_heat <- heatmap.2(inc.BP.spe, Rowv = FALSE)

library(scatterplot3d)
pop <- function(nl, rrm, nrm){
    lapply(nl, function(l){
        lapply(l, function(df){
            df[-rrm, -nrm]
        })
    })
}
mod.eval.BP <- pop(eval.BP, rrm = c(1,2,3,4,5,6), nrm = c(1,2))


prep43d <- function(mod_df){
    spe <- as.data.frame(as.table(mod_df))
    spe <- as.data.frame(apply(spe, 2, as.numeric))
    colnames(spe) <- c('-lg(Pval)','Slim Cutoff','Performance')
    spe$`-lg(Pval)` <- (-1)*log10(spe$`-lg(Pval)`)
    spe
}
real.pre.BP <- prep43d(mod.eval.BP$real[[1]])
real.spe.BP <- prep43d(mod.eval.BP$real[[2]])
random.pre.BP <- prep43d(mod.eval.BP$random[[1]])
random.spe.BP <- prep43d(mod.eval.BP$random[[2]])

s3d_CC_spe <- scatterplot3d(real.spe.BP, type = 'h', color = 'blue',
                            angle = 120, scale.y = 0.8, 
                            pch = '', main = 'BP prediction specificity if mod.mcode')
s3d_CC_spe$points3d(random.spe.BP,  type = 'h', col = 'white', pch = 1)
s3d_CC_spe$points3d(random.spe.BP,  type = 'h', col = 'white', pch = 1)
s3d_CC_spe$points3d(random.spe.BP,  type = 'p', col = 'white', pch = 16)
s3d_CC_spe$points3d(real.spe.BP, type = 'p', col = 'blue', pch = 16)
s3d_CC_spe$points3d(random.spe.BP, type = 'p', col = 'red', pch = 16)


s3d_CC_spe <- scatterplot3d(real.pre.BP, type = 'h', color = 'blue',
                            angle = 120, scale.y = 0.8, 
                            pch = '', main = 'BP prediction precision if mod.mcode')
s3d_CC_spe$points3d(random.pre.BP,  type = 'h', col = 'white', pch = 1)
s3d_CC_spe$points3d(random.pre.BP,  type = 'h', col = 'white', pch = 1)
s3d_CC_spe$points3d(random.pre.BP,  type = 'p', col = 'white', pch = 16)
s3d_CC_spe$points3d(real.pre.BP, type = 'p', col = 'blue', pch = 16)
s3d_CC_spe$points3d(random.pre.BP, type = 'p', col = 'red', pch = 16)



###### why the result become not so good> #######
### hypo: msite most correlate with the cluster that contains itself.
posML <-read_rds('Ml_epos.rds')


library(readr)
test_m <- read_rds('cor-m0724.rds')
test_exp <- read_rds('test-exp-all.rds')

library(reshape2)
m <-melt(test_m)
colnames(m) <- c('condition','site', 'M_val')
m_info_BP <- read_rds('msites_BP_fin.rds')
m$gene_id <- m_info_BP$gene_id[match(m$site, m_info_BP$modName)]
e <- melt(test_exp)
colnames(e) <- c('condition','gene_id', 'E_val')
library(dplyr)
mids <- m_info_BP$gene_id
e <-e[e$gene_id %in% mids,]

m_con <- split(m, m$condition)
e_con <- split(e, e$condition)
allexp <- Map(function(a,b){
    b$M_val <- a$M_val[match(b$gene_id, a$gene_id)]
    b
}, a = m_con, b = e_con)

essen <- lapply(allexp, function(con1){
    n <- con1 %>% dplyr::select(E_val, M_val) %>% filter(E_val < 100)
})

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
ll<-lapply(
    seq(1,40), 
    function(i) qplot(data=data.frame(essen[[i]]), E_val, M_val, alpha = 0.1)
)
do.call(grid.arrange,ll)

#example
first <- list(a = 1, b = 2, c = 3)
second <- list(a = 2, b = 3, c = 4)
Map(c, first, second)

##############

library(readr)
test_m <- read_rds('cor-m0724.rds')
test_exp <- read_rds('test-exp-all.rds')


m_info_BP <- read_rds('msites_BP_fin.rds')
library(reshape2)
library(dplyr)
pcc <- cor(test_m, test_exp, method = 'pearson')
df <- data.frame(site  = rownames(pcc))
df$geneid <- m_info_BP$gene_id[match(df$site, m_info_BP$modName)]

mpcc <- melt(pcc)
colnames(mpcc) <- c('site', 'geneid','cor')
test <- left_join(df, mpcc)



df$cor <- Map(function(a,b) {
    mpcc %>% dplyr::filter(mpcc$site ==a & mpcc$geneid == b)
}, a  = df$site, b = df$geneid)

##################
h <- read_rds('hub-random-network0805.rds')
hub_compare <- read_rds('hub-compare.rds')

ran_spe <- h$rnspe
ran_pre <- h$rnpre
true_spe <- hub_compare$real[[2]]
true_pre <- hub_compare$real[[1]]
library(gplots)
inc_pre <- true_pre - ran_pre
inc_spe <- true_spe - ran_spe
heatmap.2(inc_pre, Rowv = F, Colv = F)
heatmap.2(inc_spe, Rowv = F, Colv = F)
dimnames(inc_pre)

####################
### alternative ####
tp <- true_pre[12,]
ts <- true_spe[12,]
rp <- ran_pre[12,]
rs <- ran_spe[12,]    

random <- melt(rp)
random$specificity <- rs
random$id <- 'random'  
colnames(random) <- c('precision', 'specificity','id')


exact <- melt(tp)
exact$specificity <- ts
exact$id <- 'true'
colnames(exact)[1] <- 'precision'

all <- rbind(random, exact)
########################################################
library(ggplot2)
theme_set(theme_bw())
pdf(file = 'hub-based-performance-mod.pdf', width = 8, height = 6)
gg <- ggplot(all, aes(x=specificity, y=precision)) + 
    geom_point(aes(col=id), size = 5,alpha = 0.6) + 
    xlim(c(0.07, 0.4)) + 
    ylim(c(0.05, 0.65)) + 
    theme(legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
    labs(subtitle="Specificity vs. Precision", 
         y="Precision", 
         x="Specificity", 
         title="Hub-based performance improvement")
plot(gg)
dev.off()





#######################################################
library(reshape2)
rs <- melt(ran_spe)
colnames(rs)[3] <- 'specificity'
rp <- melt(ran_pre)
colnames(rp)[3] <- 'precision'
library(dplyr)
random <-left_join(rp,rs)
random$id <- 'random'


ts <- melt(true_spe)
colnames(ts)[3] <- 'specificity'
tp <- melt(true_pre)
colnames(tp)[3] <- 'precision'
library(dplyr)
exact <-left_join(tp,ts)
exact$id <- 'true'

alldat <- rbind(exact, random)
alldat$Var1 <- as.numeric(alldat$Var1)
alldat$Var1 <- log10(alldat$Var1)
colnames(alldat)[1] <- 'lg[Pval]'

options(scipen=999)  # turn-off scientific notation like 1e+48
library(ggplot2)
theme_set(theme_bw())  # pre-set the bw theme.

# Scatterplot
pdf(file = 'hub-based-performance.pdf', width = 10, height = 7)
gg <- ggplot(alldat, aes(x=specificity, y=precision)) + 
    geom_point(aes(col=id, size=`lg[Pval]`), alpha = 0.3) + 
    xlim(c(0.07, 0.5)) + 
    ylim(c(0.05, 0.65)) + 
    theme(legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
    labs(subtitle="Specificity vs. Precision", 
         y="Precision", 
         x="Specificity", 
         title="Hub-based performance improvement")
plot(gg)
dev.off()





#######
mod_compare <- read_rds('BP-new-compare.rds')
heatmap.2(mod_compare$improve[[1]])
heatmap.2(mod_compare$improve[[1]], Colv = F, Rowv = F)
heatmap.2(mod_compare$improve[[2]], Colv = F, Rowv = F)
ts <- mod_compare$real[[2]][-c(1,2),]
rs <- mod_compare$random[[2]][-c(1,2),]
tp <- mod_compare$real[[1]][-c(1,2),]
rp <- mod_compare$random[[1]][-c(1,2),]
a <- melt(ts)
b <- melt(rs)
a$id <- 'true'
b$id <- 'random'
dat1 <- rbind(a,b)
colnames(dat1)[3] <- 'specificity'
c <- melt(tp)
d <- melt(rp)
c$id <- 'true'
d$id <- 'random'
dat2 <- rbind(c,d)
colnames(dat2)[3] <- 'precision'
alldat2 <- left_join(dat1, dat2)
alldat2$`lg[Pval]` <- log10(as.numeric(alldat2$Var1))


options(scipen=999)  # turn-off scientific notation like 1e+48
library(ggplot2)
theme_set(theme_bw())  # pre-set the bw theme.

# Scatterplot
pdf(file = 'module-based-performance.pdf', width = 10, height = 7)
gg <- ggplot(alldat2, aes(x=specificity, y=precision)) + 
    geom_point(aes(col=id, size=`lg[Pval]`), alpha = 0.5) + 
    xlim(c(0.07, 0.22)) + 
    ylim(c(0.1, 0.37)) + 
    theme(legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
    labs(subtitle="Specificity vs. Precision", 
         y="Precision", 
         x="Specificity", 
         title="Hub-based performance improvement")
plot(gg)
dev.off()
#####################
library(readr)
pa2all <- read_rds('pa2all0731.rds')
pa2alldat <- apply(pa2all, c(1,2), function(l3){unlist(l3, recursive = FALSE)$alldat})
test <- pa2alldat[1:2,1:2]
pvec <- rownames(test)
mapply(function(df, par1, par2) {
    df$par1 = par1_vec
    df$par2 = par2_vec
    df
    }
    df = test,
    par1 = 
    )


seq <- matrix(c(list = letters[3]), 4)

###########################
ts <- mod_compare$real[[2]][12,]
rs <- mod_compare$random[[2]][12,]
tp <- mod_compare$real[[1]][12,]
rp <- mod_compare$random[[1]][12,]
random <- melt(rp)
random$specificity <- rs
random$id <- 'random'  
colnames(random) <- c('precision', 'specificity','id')


exact <- melt(tp)
exact$specificity <- ts
exact$id <- 'true'
colnames(exact)[1] <- 'precision'

all <- rbind(random, exact)
###########################
pdf(file = 'module-based-performance-mod.pdf', width = 8, height = 6)
gg <- ggplot(all, aes(x=specificity, y=precision)) + 
    geom_point(aes(col=id), alpha = 0.5, size = 5) + 
    xlim(c(0.1, 0.18)) + 
    ylim(c(0.1, 0.37)) + 
    theme(legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
    labs(subtitle="Specificity vs. Precision", 
         y="Precision", 
         x="Specificity", 
         title="Hub-based performance improvement")
plot(gg)
dev.off()

###########################




library(tidyr)
library(dplyr)
library(reshape2)
p <- do.call(rbind.fill, test0)

flatten <- function(ml){
    n <- ml
    for(i in seq(nrow(ml))){
        for (j in seq(ncol(ml))){
            df <- ml[i,j][[1]]
            df$pa1 = rownames(ml)[i]
            df$pa2 = colnames(ml)[j]
            n[i,j] <- list(df)
        }
    }
    n
}

t <- flatten(pa2alldat)
all <- do.call(rbind, t)
