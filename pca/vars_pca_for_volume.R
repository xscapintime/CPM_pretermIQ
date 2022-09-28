setwd("/mnt/d/PROJECTS/preterm_language/pca")
rm(list = ls())

library(readxl)
library(tidyverse)
library(VIM)

## merge description with pca
# load merged vars
dat_fin <- read.csv("../data/id_vars_fin.csv")

# load neonatal sickness
neo_sick <- read_excel("../data/Dana's neonatal sickness FINAL.xlsx", sheet = "Sheet2")
dim(neo_sick)
str(neo_sick)

vars_sick <- left_join(
  dat_fin,
  neo_sick,
  by = c("Eprime_ID" = "ID"),
  copy = T,
  suffix = c(".var", ".sick"),
  keep = F,
  na_matches = "never"
)

vars_sick <- vars_sick[, -c(27,26)]
colnames(vars_sick)[25] <- "Neonatal Sickness"

# load IMD
imd <- read_excel("../data/251 EAP and ePrime cases 27 9 2019 new.xlsx")
imd <- imd[,c("ID", "Y1_IMD_SCORE")]

vars_sick_imd <- left_join(
  vars_sick,
  imd,
  by = c("Eprime_ID" = "ID"),
  copy = T,
  suffix = c(".var", ".imd"),
  keep = F,
  na_matches = "never"
)

colnames(vars_sick_imd)[26] <- "IMD Score"


# filter 
vars_sick_imd <- vars_sick_imd %>% filter(!is.na(AP_ID))
row.names(vars_sick_imd) <- vars_sick_imd$AP_ID

# only keep PT
vars_sick_imd <- vars_sick_imd %>% filter(group == "PT")
dim(vars_sick_imd)
# [1] 126  24

# tidy colnams
colnames(vars_sick_imd) <- gsub(".", " ", colnames(vars_sick_imd), fixed = TRUE)


# load fd mean
fd <- read.csv("../data/fd_mean_max.txt", sep = " ")
fd <- fd[, c("Eprime_ID", "fd.m", "fd.max")]

vars_sick_imd_fd <- right_join(
  vars_sick_imd,
  fd,
  by = c("Eprime_ID" = "Eprime_ID"),
  copy = T,
  suffix = c(".var", ".fd"),
  keep = F,
  na_matches = "never"
)


## for 75 pt subj with fc
# load pca 75 subj
pca <- read.csv("../data/var17_sbj75_imp_pca.csv")

## merge vars and pca
fin75 <- right_join(
  vars_sick_imd_fd,
  pca,
  by = c("AP_ID" = "X"),
  copy = T,
  suffix = c(".var", ".pca"),
  keep = F,
  na_matches = "never"
)

# export
write.csv(fin75, file = "../data/pt_sbj75withfc_newvars_fd_pca_merged.csv", quote = F)


## impute
check_na <- cbind(
   lapply(
     lapply(fin75, is.na)
     , sum)
)

names(check_na[check_na[,1] != 0,])

knn75 <- kNN(fin75, variable=names(check_na[check_na[,1] != 0,]), k=3, imp_var=F)

# export
write.csv(knn75, file = "../data/pt_sbj75withfc_newvars_imp_fd_pca_merged.csv", quote = F)



#### would not work as have to merge with fd
# for all 126 pt subj
# load pca 126 subj
# pca <- read.csv("../data/var17_sbj116_imp_pca.csv")
pca <- read.csv("../data/var17_sbj126_imp_pca.csv")

## merge vars and pca
fin126 <- full_join(
  vars_sick_imd,
  pca,
  by = c("AP_ID" = "X"),
  copy = T,
  suffix = c(".var", ".pca"),
  keep = F,
  na_matches = "never"
)

# export
# write.csv(fin116, file = "../data/pt_sbj116_newvars_pca_merged.csv", quote = F)
write.csv(fin126, file = "../data/pt_sbj126_newvars_pca_merged.csv", quote = F)


## impute
check_na <- cbind(
   lapply(
     lapply(fin126, is.na)
     , sum)
)

names(check_na[check_na[,1] != 0,])

knn126 <- kNN(fin126, variable=names(check_na[check_na[,1] != 0,]), k=3, imp_var=F)

# export
# write.csv(knn116, file = "../data/pt_sbj116_newvars_imp_pca_merged.csv", quote = F)
write.csv(knn126, file = "../data/pt_sbj126_newvars_imp_pca_merged.csv", quote = F)