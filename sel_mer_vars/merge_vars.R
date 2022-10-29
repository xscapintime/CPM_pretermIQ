# setwd("/mnt/d/PROJECTS/preterm_language/sel_mer_vars")
rm(list = ls())

library(readxl)
library(tidyverse)

## Load table
tb1 <- read_excel("../data/AP_marking_17_aug21 sept 28 22.xlsx", sheet = "Master") #AP_marking_JULY 2022.xlsx
tb1 <- tb1[-1, ]
#tb1 <- tb1 %>% filter(!is.na(AP_ID) & !is.na(Eprime_ID))
#tb1$Eprime_ID <- as.numeric(tb1$Eprime_ID) %>% as.character()
tb1 <- tb1 %>% filter(grepl("AP|BIPP", AP_ID))

vars <- c("AP_ID", "Sex (1=m, 2=f)", "WISC_FULL_CS", "WISC_VCI_CS", "WISC_PR_CS", "WISC_WM_CS",
            "WISC_PS_CS", "srs-rrb", "srs-sci")
dat <- tb1[,vars]

## mark group
dat$group <- ifelse(tb1$GA > 38, "FT", "PT")
dim(dat)
# [1] 188 10

# based on excel color
dat[dat$AP_ID %in% c("BIPP047", "BIPP048", "BIPP050", "BIPP051"),]$group <- "FT"
dat[is.na(dat$group),]$group <- "PT"


table(dat$group)
#  FT  PT 
#  56 132

# -999 and -998
# tb_merge[tb_merge == -999] <- NA
# tb_merge[tb_merge == -998] <- NA


## merge duplicated
dup_ap <- na.omit(dat$AP_ID)[na.omit(dat$AP_ID) %>% duplicated()]
dup_ap
# no duplicates



############################ no need #####################
#dup_ep <- na.omit(dat$Eprime_ID)[na.omit(dat$Eprime_ID) %>% duplicated()]
#dup_ep
# no duplicates

# check overlap
# (dat %>% filter(AP_ID %in% dup_ap))$Eprime_ID %>% na.omit() %in% dup_ep
# no overlap

# library(purrrlyr)
# library(zoo)
# # https://stackoverflow.com/a/40046888/14498100
# # dat %>% filter(AP_ID %in% dup_ap) %>% group_by(AP_ID) %>%  by_slice(function(x) { 
# #     na.locf(na.locf(x, na.rm=F), fromLast=T, na.rm=F) }, 
# #     .collate = "rows") %>% distinct()

# dup_clean <- dat %>% filter(Eprime_ID %in% dup_ep) %>% group_by(Eprime_ID) %>%  by_slice(function(x) { 
#     na.locf(na.locf(x, na.rm=F), fromLast=T, na.rm=F) }, 
#     .collate = "rows") %>% distinct() %>% relocate(Eprime_ID, .after = AP_ID)

# # remove and bind
# dat_fin <- rbind(dat %>% filter(!Eprime_ID %in% dup_ep & !AP_ID %in% dup_ap), dup_clean) %>% arrange(AP_ID, Eprime_ID)

#######################################################




# sex to 1=m, 2=f, merge 2 sex column
# dat$sex <- ifelse(dat$s == "Male", 1, 2)
dat$`Sex (1=m, 2=f)` <- dat$`Sex (1=m, 2=f)` %>% as.numeric()
# dat <- dat %>% mutate(sex_mer = coalesce(`Sex (1=m, 2=f)`, sex))
colnames(dat)[2] <- "sex"


# check data types
str(dat)
dat$WISC_FULL_CS <- dat$WISC_FULL_CS %>% as.numeric()
dat$WISC_VCI_CS <- dat$WISC_VCI_CS %>% as.numeric()
dat$WISC_PR_CS <- dat$WISC_PR_CS %>% as.numeric()
dat$WISC_WM_CS <- dat$WISC_WM_CS %>% as.numeric()
dat$WISC_PS_CS <- dat$WISC_PS_CS %>% as.numeric()


## ">90" --> 90
dat$`srs-rrb`[!grepl("^\\d+", dat$`srs-rrb`)][!dat$`srs-rrb`[!grepl("^\\d+$", dat$`srs-rrb`)] %>% is.na()] <- "90"
dat$`srs-sci`[!grepl("^\\d+", dat$`srs-sci`)][!dat$`srs-sci`[!grepl("^\\d+$", dat$`srs-sci`)] %>% is.na()] <- "90"
# dat$`srs-tot`[!grepl("^\\d+", dat$`srs-tot`)][!dat$`srs-tot`[!grepl("^\\d+$", dat$`srs-tot`)] %>% is.na()] <- "90"

dat[, c("srs-rrb", "srs-sci")] <- sapply(dat[, c("srs-rrb", "srs-sci")], as.numeric)


# check data types
summary(dat)
str(dat)


# export
write.table(dat, file = "../data/vars_wisc_srs_8yo.csv", sep = ",",
              quote = F, row.names = F, col.names = T)
