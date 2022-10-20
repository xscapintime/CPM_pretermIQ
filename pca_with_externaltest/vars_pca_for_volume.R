setwd("/mnt/d/PROJECTS/preterm_language/pca_with_externaltest")
rm(list = ls())

library(readxl)
library(tidyverse)
library(VIM)

## merge description with pca
# load merged vars
dat_fin <- read.csv("../data/id_vars_fin.csv")

## only keep 8 yo data
dat_fin <- dat_fin[,c(1,2,3,4,7,8,9,10,11,13,14)] # srs rrb and sci

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

vars_sick <- vars_sick[, -c(13,14)]
colnames(vars_sick)[12] <- "Neonatal Sickness"

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

colnames(vars_sick_imd)[13] <- "IMD Score"


# filter 
vars_sick_imd <- vars_sick_imd %>% filter(!is.na(AP_ID))
row.names(vars_sick_imd) <- vars_sick_imd$AP_ID

# only keep PT
vars_sick_imd <- vars_sick_imd %>% filter(group == "PT")
dim(vars_sick_imd)
# [1] 126  13

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


## for 82 pt subj with fc
# load pca 82 subj
pca <- read.csv("../data/var6_8yo_sbj82_imp_pca.csv")

## merge vars and pca
fin82 <- right_join(
  vars_sick_imd_fd,
  pca,
  by = c("AP_ID" = "X"),
  copy = T,
  suffix = c(".var", ".pca"),
  keep = F,
  na_matches = "never"
)

# export
write.csv(fin82, file = "../data/pt_sbj82withfc_newvars_fd_pca_merged.csv", quote = F)


## impute
check_na <- cbind(
   lapply(
     lapply(fin82, is.na)
     , sum)
)

names(check_na[check_na[,1] != 0,])

knn82 <- kNN(fin82, variable=names(check_na[check_na[,1] != 0,]), k=3, imp_var=F)

# export
write.csv(knn82, file = "../data/pt_sbj82withfc_newvars_imp_fd_pca_merged.csv", quote = F)


# temp <- fin82[,-c(5,6,12,seq(15,28))]

# check_na_t <- cbind(
#    lapply(
#      lapply(temp, is.na)
#      , sum)
# )

# names(check_na_t[check_na_t[,1] != 0,])
# knnt <- kNN(temp, variable=names(check_na_t[check_na_t[,1] != 0,]), k=3, imp_var=F)
