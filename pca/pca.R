setwd("/mnt/d/PROJECTS/preterm_language/pca")

library(readxl)
library(tidyverse)
library(factoextra)
# library(FactoMineR)
library(corrplot)
library(psych)


## Load table
# from final version of vars
dat_fin <- read.csv("../data/id_vars_fin.csv")
row.names(dat_fin) <- dat_fin$AP_ID

# column as variables
vars <- c("WISC_VCI_CS", "WISC_PR_CS", "WISC_WM_CS", "WISC_PS_CS", "srs-rrb", "srs-sci",
        "bayley22_cog_comp", "bayley22_language_comp", "bayley22_motor_comp", "parca22_cognitive", "parca22_language",
        "wppsi4_verb_compr_raw", "wppsi4_visuo_sp_raw", "wppsi4_fluid_res_raw", "wppsi4_working_mem_raw", "wppsi4_proc_speed_raw", "srs4_rrb_raw", "srs4_sci_raw",
        "MC_COUNT_TOTAL_FAILS_nooffails")

dat <- tb_merge[,vars]


## PCA
dat <- sapply(dat, as.numeric)
res <- prcomp(na.omit(dat), center = T, scale = T) ##?

summary(res)

# PCA permutation
# https://github.com/lucyvanes/preterm-outcomes/blob/main/1_PCA/1_PCA.R
#===========================================
#    Run PCA with permutation testing
#===========================================
# for each permutation, shuffle rows in each column, re-compute PCA
# function for the permutation testing taken from:
# http://bioops.info/2015/01/permutation-pca/


sign.pc<-function(x,R=5000,s=10, cor=T,...){
  pc.out<-princomp(x,cor=cor,...)  # run PCA
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s] # proportion of variance for each PC
  pve.perm<-matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    x.perm<-apply(x,2,sample)
    pc.perm.out<-princomp(x.perm,cor=cor,...)
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s]
  }
  pval<-apply(t(pve.perm)>pve,1,sum)/R # calculate p-values
  return(list(pve=pve,pval=pval))
}

pca_sign <- sign.pc(as.data.frame(scale(na.omit(dat), center=T, scale=T)),cor=T)
pca_sign # the first two PCs are significant


## p-val
# somthing wrong
## pmat <- corr.test(res$x, adjust = "none")$p

library(Hmisc)

pmat <- rcorr(res$x, na.omit(dat), type = c("pearson"))$P %>% round(., digits = 4)
pmat <- pmat[-(1:ncol(res$x)), -(ncol(na.omit(dat)):ncol(pmat))]

rmat <- rcorr(res$x, na.omit(dat), type = c("pearson"))$r %>% round(., digits = 4)
rmat <- rmat[-(1:ncol(res$x)), -(ncol(na.omit(dat)):ncol(rmat))]



# viz
png("scree.png", width = 700, height = 580)
fviz_eig(res)
dev.off()


png("scree_allpc.png", width = 700, height = 580)
fviz_eig(res, ncp = 19)
dev.off()



png("corr.png", width = 580, height = 580)
corrplot(res$rotation, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20')
dev.off()


png("corr_new.png", width = 580, height = 580)
corrplot(rmat, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20')
dev.off()