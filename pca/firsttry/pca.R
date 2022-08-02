setwd("/mnt/d/PROJECTS/preterm_language/pca")

library(readxl)
library(factoextra)
# library(FactoMineR)
library(corrplot)
library(psych)


## Load table
ep_tb <- read_excel("ePrime_BIPP_master_file_GEORGE_LHv4 correct tmcq.xlsx")


## filter
# -999 and -998
ep_tb[ep_tb == -999] <- NA
ep_tb[ep_tb == -998] <- NA

# column as variables
# bayley22_language_comp, parca22_language, wppsi4_verb_compr_raw, wisc8_vci_cs
# no sex
dat <- ep_tb[,c("bayley22_language_comp", "parca22_language", "wppsi4_verb_compr_raw", "wisc8_vci_cs")]


## PCA
res <- prcomp(na.omit(dat), center = TRUE, scale = T) ##?

## p-val
pmat <- corr.test(res$rotation)$p

# viz
png("scree_4v.png",  width = 700, height = 580)
fviz_eig(res)
dev.off()


png("corr_4v.png", width = 580, height = 580)
corrplot(res$rotation, p.mat = pmat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20')
dev.off()



