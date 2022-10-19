rm(list = ls())
setwd("/mnt/d/PROJECTS/preterm_language/pca_with_externaltest")
# and only 8 yo vars

library(readxl)
library(tidyverse)
library(factoextra)
# library(FactoMineR)
library(corrplot)
library(psych)
library(VIM)
library(RColorBrewer)

## function
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



## Load table
# from final version of vars
dat_fin <- read.csv("../data/id_vars_fin.csv")
dim(dat_fin)

dat_fin <- dat_fin %>% filter(!is.na(AP_ID))
row.names(dat_fin) <- dat_fin$AP_ID

# only keep PT
dat_fin <- dat_fin %>% filter(group == "PT")
dim(dat_fin)
# [1] 126  28

# tidy colnams
colnames(dat_fin) <- gsub(".", " ", colnames(dat_fin), fixed = TRUE)

## only keep 8 yo data
dat_fin <- dat_fin[,c(1,2,3,4,7,8,9,10,11,13,14)] # srs rrb and sci


# column as variables
# vars <- c("WISC_VCI_CS", "WISC_PR_CS", "WISC_WM_CS", "WISC_PS_CS", "srs-rrb", "srs-sci",
#         "bayley22_cog_comp", "bayley22_language_comp", "bayley22_motor_comp", "parca22_cognitive", "parca22_language",
#         "wppsi4_verb_compr_raw", "wppsi4_visuo_sp_raw", "wppsi4_fluid_res_raw", "wppsi4_working_mem_raw", "wppsi4_proc_speed_raw", "srs4_rrb_raw", "srs4_sci_raw",
#         "MC_COUNT_TOTAL_FAILS_nooffails")


data <- dat_fin[,-c(1:5)]
dim(data)
# [1] 126 6
colnames(data)                          
# [1] "WISC VC 8yo"
# [2] "WISC PR 8yo"
# [3] "WISC WM 8yo"
# [4] "WISC PS 8yo"
# [5] "SRS RRB 8yo"
# [6] "SRS SCI 8yo"


## match with subj with FC
fc_files <- list.files("../data/fc_individual/", ".txt$")
subj <- gsub(".fc.txt", "", fc_files)

dat_fc <- data[rownames(data) %in% subj,]
dim(dat_fc)
# [1] 82  6

## remove vars that NA > 30% subj
## check first
for (i in 1:(ncol(dat_fc))) {
  n_na <- is.na(dat_fc[,i]) %>% sum()
  if (n_na / nrow(dat_fc) > 0) {
      print(paste0(colnames(dat_fc)[i], ": ", n_na , ", ", round(is.na(dat_fc[,i]) %>% sum() / nrow(dat_fc), 2)))
  }

}

# [1] "WISC PS 8yo: 1, 0.01"
# [1] "SRS Total 8yo: 4, 0.05"


## check subjs that with a lot NAs
for (i in 1:(nrow(dat_fc))) {
  if (is.na(dat_fc[i,]) %>% sum() != 0) {
    print(paste0(row.names(dat_fc)[i], ": ", is.na(dat_fc[i,]) %>% sum()))
  }
}

# [1] "AP010: 2"
# [1] "AP011: 2"
# [1] "BIPP002: 2"
# [1] "BIPP003: 2"
# [1] "BIPP021: 1"


dat_keep <- dat_fc[apply(is.na(dat_fc), 1, sum) <= 6,]

# check var NA%
for (i in 1:(ncol(dat_keep))) {
  # if (is.na(dat_keep[,i]) %>% sum() / nrow(dat_keep) > 0.3) {
    print(paste0(colnames(dat_keep)[i], ": ", round(is.na(dat_keep[,i]) %>% sum() / nrow(dat_keep), 2)))
  # }

}
dim(dat_keep)
# [1] 82 6


## PCA original wmethod
## remove all NA
dat <- sapply(data, as.numeric)
na.omit(dat) %>% dim
# [1] 113 6

row.names(dat) <- row.names(data)
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
pca_sign # the first 2 PC is significant

# export PCs
# write.csv(res$x, file = "../data/var19_pca.csv", quote = F)
write.csv(res$x, file = "../data/var6_8yo_pca.csv", quote = F)


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


# png("scree_allpc.png", width = 700, height = 580)
# fviz_eig(res, ncp = 19)
# dev.off()



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
# KNN imputation, subjs have FC
## impute
check_na <- cbind(
   lapply(
     lapply(dat_keep, is.na)
     , sum)
)

names(check_na[check_na[,1] != 0,])

knn <- kNN(dat_keep, variable = names(check_na[check_na[,1] != 0,]), k=3, imp_var=F)
row.names(knn) <- row.names(dat_keep)
dim(knn)
# [1] 82  6

# export imputed vars
write.csv(knn, file = "../data/vars_8yo_82subj_imputed.csv", quote = F)


### PCA
dat <- sapply(knn, as.numeric)
row.names(dat) <- row.names(knn)
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
pca_sign # the first 2 PC is significant

# export PCs
# write.csv(res$x, file = "../data/var19_imp_pca.csv", quote = F)
write.csv(res$x, file = "../data/var6_8yo_sbj82_imp_pca.csv", quote = F)



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


# png("scree_allpc_imp.png", width = 700, height = 580)
# fviz_eig(res, ncp = 19)
# dev.off()



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




# KNN imputation, all subjs w/wo FC
## impute
check_na <- cbind(
   lapply(
     lapply(data, is.na)
     , sum)
)

names(check_na[check_na[,1] != 0,])

knn <- kNN(data, variable = names(check_na[check_na[,1] != 0,]), k=3, imp_var=F)
row.names(knn) <- row.names(data)
dim(knn)
# [1] 126   6

# export imputed vars
# write.csv(knn, file = "../data/vars_116subj_imputed.csv", quote = F)
write.csv(knn, file = "../data/vars_8yo_126subj_imputed.csv", quote = F)


### PCA
dat <- sapply(knn, as.numeric)
row.names(dat) <- row.names(knn)
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
pca_sign # the first 2 PC is significant

# export PCs
# write.csv(res$x, file = "../data/var19_imp_pca.csv", quote = F)
# write.csv(res$x, file = "../data/var17_sbj116_imp_pca.csv", quote = F)
write.csv(res$x, file = "../data/vars_8yo_sbj126_imp_pca.csv", quote = F)



## p-val
# somthing wrong
## pmat <- corr.test(res$x, adjust = "none")$p

library(Hmisc)

pmat <- rcorr(res$x, na.omit(dat), type = c("pearson"))$P %>% round(., digits = 4)
pmat <- pmat[-(1:ncol(res$x)), -(ncol(na.omit(dat)):ncol(pmat))]

rmat <- rcorr(res$x, na.omit(dat), type = c("pearson"))$r %>% round(., digits = 4)
rmat <- rmat[-(1:ncol(res$x)), -(ncol(na.omit(dat)):ncol(rmat))]



# viz
# png("scree_imp_sbj116.png", width = 700, height = 580)
png("scree_imp_sbj126.png", width = 700, height = 580)
fviz_eig(res)
dev.off()


# png("scree_allpc_imp_sbj116.png", width = 700, height = 580)
# png("scree_allpc_imp_sbj126.png", width = 700, height = 580)
# fviz_eig(res, ncp = 19)
# dev.off()



# png("corr_rotation_imp_sbj116.png", width = 580, height = 580)
png("corr_rotation_imp_sbj126.png", width = 580, height = 580)
corrplot(res$rotation, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20',
         col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         tl.col = "black")
dev.off()


# png("corr_rmat_imp_sbj116.png", width = 580, height = 580)
png("corr_rmat_imp_sbj126.png", width = 580, height = 580)
corrplot(rmat, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20',
         col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         tl.col = "black")
dev.off()