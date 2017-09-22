# objective
# filter genes, so that:
# 1. log2 transformed
# 2. standardize to N(0,1)
# 3. top 90% variance & top 90% mean

# ideal size: 2w

library(readr)
library(dplyr)
uni_exp <- read_rds('test_exp0921.rds')
# shortcut for visualizing large matrix
look <- function(m) m[1:10, 1:10]

# log2 transform
expm <- log2(uni_exp + 0.01)

# scale to N(0,1)
m <- mean(expm)
sd <- sd(expm)
s.expm <- apply(expm, c(1,2), function(x) (x - m)/sd )

# fliter required sd and mean cutoffs
sdnm_filter <- function(m, qm, qsd){
    # calculate
    site_ms <- apply(m, 2, mean)
    site_sds <- apply(m, 2, sd)
    
    m_cut <- quantile(site_ms, 1 - qm)
    sd_cut <- quantile(site_sds, 1 - qsd)
    
    m_retain <- (site_ms > m_cut)
    sd_retain <- (site_sds > sd_cut)
    
    retain <- mapply(function(a,b) a & b, a = m_retain, b = sd_retain)
    m[,retain]
}

test_m <- sdnm_filter(s.expm, 0.9 ,  0.9)
write_rds(test_m, 'test_m0922.rds')
