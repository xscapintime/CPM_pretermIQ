# setwd("/mnt/d/PROJECTS/preterm_language/pca")
rm(list = ls())

library(readxl)
library(tidyverse)

## Load table
tb1 <- read_excel("../data/AP_marking_17_aug21 sept 28 22.xlsx", sheet = "Master") #AP_marking_JULY 2022.xlsx
tb1 <- tb1[-1, ]
tb1 <- tb1 %>% filter(!is.na(AP_ID) & !is.na(Eprime_ID))
tb1$Eprime_ID <- as.numeric(tb1$Eprime_ID) %>% as.character()
tb1 <- tb1 %>% filter(grepl("AP|BIPP", AP_ID))

tb2 <- read_excel("../data/ePrime_BIPP_master_file_GEORGE_LHv4 correct tmcq.xlsx")
tb2$EP_id <- as.character(tb2$EP_id)




# tmp <- merge(tb1, tb2, by.x = "AP_ID", by.y = "AP_id")
# tb_merge <- merge(tmp, tb3, by.x = "Eprime_ID", by.y = "id")


tb_merge <- full_join(
  tb1,
  tb2,
  by = c("AP_ID" = "AP_id", "Eprime_ID" = "EP_id"),
  copy = T,
  suffix = c(".1", ".2"),
  keep = F,
  na_matches = "never"
)

# -999 and -998
tb_merge[tb_merge == -999] <- NA
tb_merge[tb_merge == -998] <- NA


# add srs total and seperate values
vars <- c("AP_ID", "Eprime_ID", "group", "Sex (1=m, 2=f)", "sex", "age8",
        "WISC_FULL_CS", "wisc8_full_cs")


dat <- tb_merge[,vars]

## merge duplicated
dup_ap <- na.omit(dat$AP_ID)[na.omit(dat$AP_ID) %>% duplicated()]
dup_ap
# [1] "BIPP015" "BIPP016"

dup_ep <- na.omit(dat$Eprime_ID)[na.omit(dat$Eprime_ID) %>% duplicated()]
dup_ep
#  [1] "2233" "6324" "5534" "9431" "5832" "6594" "9458"
#  [8] "7153" "6415" "6652" "7224" "9429" "7448" "7654"
# [15] "7519" "9503" "7699" "9406" "7427" "9444" "9511"
# [22] "7439" "9664" "7207" "6470" "6487" "7238" "7376"
# [29] "7275" "9412" "9465" "7369" "7687"

# check overlap
(dat %>% filter(AP_ID %in% dup_ap))$Eprime_ID %>% na.omit() %in% dup_ep



library(purrrlyr)
library(zoo)
# https://stackoverflow.com/a/40046888/14498100
# dat %>% filter(AP_ID %in% dup_ap) %>% group_by(AP_ID) %>%  by_slice(function(x) { 
#     na.locf(na.locf(x, na.rm=F), fromLast=T, na.rm=F) }, 
#     .collate = "rows") %>% distinct()

dup_clean <- dat %>% filter(Eprime_ID %in% dup_ep) %>% group_by(Eprime_ID) %>%  by_slice(function(x) { 
    na.locf(na.locf(x, na.rm=F), fromLast=T, na.rm=F) }, 
    .collate = "rows") %>% distinct() %>% relocate(Eprime_ID, .after = AP_ID)

# remove and bind
dat_fin <- rbind(dat %>% filter(!Eprime_ID %in% dup_ep & !AP_ID %in% dup_ap), dup_clean) %>% arrange(AP_ID, Eprime_ID)

# sex to 1=m, 2=f, merge 2 sex column
dat_fin$sex <- ifelse(dat_fin$sex == "Male", 1, 2)
dat_fin$`Sex (1=m, 2=f)` <- dat_fin$`Sex (1=m, 2=f)` %>% as.numeric()

dat_fin <- dat_fin %>% mutate(sex_mer = coalesce(`Sex (1=m, 2=f)`, sex))


# merge 2 wisc
dat_fin$WISC_FULL_CS <- dat_fin$WISC_FULL_CS %>% as.numeric()
dat_fin <- dat_fin %>% mutate(wisc_full = coalesce(WISC_FULL_CS, wisc8_full_cs))

# check data types
summary(dat_fin)
str(dat_fin)

colnames(dat_fin)[4] <- "sex_AP"

# export
write.table(dat_fin, file = "../data/fulliq_8yo.csv", sep = ",",
              quote = F, row.names = F, col.names = T)
