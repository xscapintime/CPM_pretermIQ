setwd("/mnt/d/PROJECTS/preterm_language/pca")

library(readxl)
library(tidyverse)
library(factoextra)
# library(FactoMineR)
library(corrplot)
library(psych)


## Load table
tb1 <- read_excel("../data/AP_marking_JULY 2022.xlsx", sheet = "Master")
tb1 <- tb1[-1, ]
tb2 <- read_excel("../data/ePrime_BIPP_master_file_GEORGE_LHv4 correct tmcq.xlsx")
tb3 <- read_excel("../data/mchat_at22m_EP_ids.xlsx")


## merge
tmp <- merge(tb1, tb2, by.x = "AP_ID", by.y = "AP_id")
tb_merge <- merge(tmp, tb3, by.x = "Eprime_ID", by.y = "id")


## filter
# -999 and -998
tb_merge[tb_merge == -999] <- NA
tb_merge[tb_merge == -998] <- NA

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



