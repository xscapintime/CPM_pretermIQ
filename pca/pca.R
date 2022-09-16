setwd("/mnt/d/PROJECTS/preterm_language/pca")
rm(list = ls())

library(readxl)
library(tidyverse)
library(factoextra)
# library(FactoMineR)
library(corrplot)
library(psych)
library(caret)
library(RColorBrewer)


## Load table
# from final version of vars
dat_fin <- read.csv("../data/id_vars_fin.csv")
dat_fin <- dat_fin %>% filter(!is.na(AP_ID))
row.names(dat_fin) <- dat_fin$AP_ID

# column as variables
# vars <- c("WISC_VCI_CS", "WISC_PR_CS", "WISC_WM_CS", "WISC_PS_CS", "srs-rrb", "srs-sci",
#         "bayley22_cog_comp", "bayley22_language_comp", "bayley22_motor_comp", "parca22_cognitive", "parca22_language",
#         "wppsi4_verb_compr_raw", "wppsi4_visuo_sp_raw", "wppsi4_fluid_res_raw", "wppsi4_working_mem_raw", "wppsi4_proc_speed_raw", "srs4_rrb_raw", "srs4_sci_raw",
#         "MC_COUNT_TOTAL_FAILS_nooffails")


dat <- dat_fin[,-c(1:7)]
dim(dat)

# change colnames
colnames(dat) <- c("WISC VC 8yo", "WISC PR 8yo", "WISC WM 8yo", "WISC PS 8yo", "SRS TOTAL 8yo",
        "Bayley Cog 22mo", "Bayley Lang 22mo", "Bayley Motor 22mo", "PARCA Cog 22mo", "PARCA Lang 22mo",
        "WPPSI VC 4yo", "WPPSI VS 4yo", "WPPSI FR 4yo", "WPPSI WM 4yo", "WPPSI PS 4yo", "SRS TOTAL 4yo",
        "MCHAT Failed 8yo")


dat <- dat[,c("WISC VC 8yo", "WISC PR 8yo", "WISC WM 8yo", "WISC PS 8yo", "SRS TOTAL 8yo", "MCHAT Failed 8yo",
        "WPPSI VC 4yo", "WPPSI VS 4yo", "WPPSI FR 4yo", "WPPSI WM 4yo", "WPPSI PS 4yo", "SRS TOTAL 4yo",
        "Bayley Cog 22mo", "Bayley Lang 22mo", "Bayley Motor 22mo", "PARCA Cog 22mo", "PARCA Lang 22mo")] 



## match with subj with FC
fc_files <- list.files("../data/fc_individual/", ".txt$")
subj <- gsub(".fc.txt", "", fc_files)

dat <- dat[subj,]
dim(dat)


## remove vars that NA > 30% subj
## check first
for (i in 1:(ncol(dat))) {
  if (is.na(dat[,i]) %>% sum() / nrow(dat) > 0.3) {
      print(paste0(colnames(dat)[i], ": ", round(is.na(dat[,i]) %>% sum() / nrow(dat), 2)))
  }

}  # all of vars have > 30% NA

# [1] "WISC VC 8yo: 0.36"
# [1] "WISC PR 8yo: 0.36"
# [1] "WISC WM 8yo: 0.36"
# [1] "WISC PS 8yo: 0.37"
# [1] "SRS TOTAL 8yo: 0.39"
# [1] "MCHAT Failed 8yo: 0.38"
# [1] "WPPSI VC 4yo: 0.38"
# [1] "WPPSI VS 4yo: 0.38"
# [1] "WPPSI FR 4yo: 0.38"
# [1] "WPPSI WM 4yo: 0.37"
# [1] "WPPSI PS 4yo: 0.38"
# [1] "SRS TOTAL 4yo: 0.39"
# [1] "Bayley Cog 22mo: 0.38"
# [1] "Bayley Lang 22mo: 0.38"
# [1] "Bayley Motor 22mo: 0.38"
# [1] "PARCA Cog 22mo: 0.39"
# [1] "PARCA Lang 22mo: 0.39"

## check subjs that with a lot NAs
for (i in 1:(nrow(dat))) {
  if (is.na(dat[i,]) %>% sum() != 0) {
    print(paste0(row.names(dat)[i], ": ", is.na(dat[i,]) %>% sum()))
  }
}

# [1] "AP004: 2"
# [1] "AP005: 2"
# [1] "AP006: 17"
# [1] "AP009: 17"
# [1] "AP010: 1"
# [1] "AP011: 1"
# [1] "AP012: 1"
# [1] "AP016: 17"
# [1] "AP019: 17"
# [1] "AP020: 1"
# [1] "AP023: 17"
# [1] "AP025: 17"
# [1] "AP026: 17"
# [1] "AP027: 17"
# [1] "AP030: 17"
# [1] "AP033: 17"
# [1] "AP034: 17"
# [1] "AP037: 17"
# [1] "AP045: 17"
# [1] "AP046: 17"
# [1] "AP061: 17"
# [1] "AP070: 17"
# [1] "AP071: 17"
# [1] "AP073: 17"
# [1] "AP075: 17"
# [1] "AP077: 17"
# [1] "AP079: 17"
# [1] "AP081: 17"
# [1] "AP082: 17"
# [1] "AP084: 17"
# [1] "AP085: 17"
# [1] "AP089: 17"
# [1] "AP090: 17"
# [1] "AP092: 17"
# [1] "AP094: 17"
# [1] "AP095: 17"
# [1] "AP096: 17"
# [1] "AP097: 17"
# [1] "AP098: 17"
# [1] "AP102: 17"
# [1] "AP109: 17"
# [1] "AP120: 17"
# [1] "AP121: 1"
# [1] "AP122: 17"
# [1] "AP123: 17"
# [1] "AP129: 17"
# [1] "BIPP002: 1"
# [1] "BIPP003: 1"
# [1] "BIPP008: 1"
# [1] "BIPP010: 17"
# [1] "BIPP012: 17"
# [1] "BIPP013: 17"
# [1] "BIPP014: 17"
# [1] "BIPP017: 6"
# [1] "BIPP021: 1"
# [1] "BIPP022: 1"
# [1] "BIPP023: 1"
# [1] "BIPP027: 12"


dat_keep <- dat[apply(is.na(dat), 1, sum) <= 6,]

# check var NA%
for (i in 1:(ncol(dat_keep))) {
  # if (is.na(dat_keep[,i]) %>% sum() / nrow(dat_keep) > 0.3) {
    print(paste0(colnames(dat_keep)[i], ": ", round(is.na(dat_keep[,i]) %>% sum() / nrow(dat_keep), 2)))
  # }

}
dim(dat_keep)
# [1] 75 17


## PCA original wmethod
## remove all NA
dat <- sapply(dat, as.numeric)
na.omit(dat) %>% dim
# [1] 61 17

row.names(dat) <- subj
res <- prcomp(na.omit(dat), center = T, scale = T) ## no NA allowed, 89 subjects left

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
pca_sign # the first 3 PCs are significant

# export PCs
# write.csv(res$x, file = "../data/var19_pca.csv", quote = F)
write.csv(res$x, file = "../data/var17_pca.csv", quote = F)


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



png("corr_rotation.png", width = 580, height = 580)
corrplot(res$rotation, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20',
         col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         tl.col = "black")
dev.off()


png("corr_rmat.png", width = 580, height = 580)
corrplot(rmat, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20',
         col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         tl.col = "black")
dev.off()



## PCA after KNN imputation
# KNN imputation
imp <- preProcess(dat_keep, method = "knnImpute", k = 5)
demo_imp <- predict(imp, dat_keep)
dim(demo_imp)
# [1] 75 17

### PCA
dat <- sapply(demo_imp, as.numeric)
row.names(dat) <- row.names(demo_imp)
res <- prcomp(na.omit(dat), center = T, scale = T) ## no NA allowed, 89 subjects left

summary(res)

# PCA permutation
# https://github.com/lucyvanes/preterm-outcomes/blob/main/1_PCA/1_PCA.R
#===========================================
#    Run PCA with permutation testing
#===========================================
# for each permutation, shuffle rows in each column, re-compute PCA
# function for the permutation testing taken from:
# http://bioops.info/2015/01/permutation-pca/


pca_sign <- sign.pc(as.data.frame(scale(na.omit(dat), center=T, scale=T)),cor=T)
pca_sign # the first 2 PCs are significant

# export PCs
# write.csv(res$x, file = "../data/var19_imp_pca.csv", quote = F)
write.csv(res$x, file = "../data/var17_imp_pca.csv", quote = F)



## p-val
# somthing wrong
## pmat <- corr.test(res$x, adjust = "none")$p

library(Hmisc)

pmat <- rcorr(res$x, na.omit(dat), type = c("pearson"))$P %>% round(., digits = 4)
pmat <- pmat[-(1:ncol(res$x)), -(ncol(na.omit(dat)):ncol(pmat))]

rmat <- rcorr(res$x, na.omit(dat), type = c("pearson"))$r %>% round(., digits = 4)
rmat <- rmat[-(1:ncol(res$x)), -(ncol(na.omit(dat)):ncol(rmat))]



# viz
png("scree_imp.png", width = 700, height = 580)
fviz_eig(res)
dev.off()


png("scree_allpc_imp.png", width = 700, height = 580)
fviz_eig(res, ncp = 19)
dev.off()



png("corr_rotation_imp.png", width = 580, height = 580)
corrplot(res$rotation, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20',
         col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         tl.col = "black")
dev.off()


png("corr_rmat_imp.png", width = 580, height = 580)
corrplot(rmat, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20',
         col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         tl.col = "black")
dev.off()