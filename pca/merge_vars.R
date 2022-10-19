setwd("/mnt/d/PROJECTS/preterm_language/pca")
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


tb3 <- read_excel("../data/mchat_at22m_EP_ids.xlsx")
tb3$id <- as.character(tb3$id)
tb3$MC_COUNT_TOTAL_FAILS_nooffails <- as.numeric(tb3$MC_COUNT_TOTAL_FAILS_nooffails)


# tmp <- merge(tb1, tb2, by.x = "AP_ID", by.y = "AP_id")
# tb_merge <- merge(tmp, tb3, by.x = "Eprime_ID", by.y = "id")


tmp <- full_join(
  tb1,
  tb2,
  by = c("AP_ID" = "AP_id", "Eprime_ID" = "EP_id"),
  copy = T,
  suffix = c(".1", ".2"),
  keep = F,
  na_matches = "never"
)

tb_merge <- left_join(
  tmp,
  tb3,
  by = c("Eprime_ID" = "id"),
  copy = T,
  suffix = c(".12", ".3"),
  keep = F,
  na_matches = "never"
)


# -999 and -998
tb_merge[tb_merge == -999] <- NA
tb_merge[tb_merge == -998] <- NA


# add srs total and seperate values
vars <- c("AP_ID", "Eprime_ID", "group", "sex", "age22", "age4", "age8",
        "WISC_VCI_CS", "WISC_PR_CS", "WISC_WM_CS", "WISC_PS_CS", "srs-tot", "srs-rrb", "srs-sci", ## AP_marking_JULY, 8 yo
        "bayley22_cog_comp", "bayley22_language_comp", "bayley22_motor_comp", "parca22_cognitive", "parca22_language",
        "wppsi4_verb_compr_raw", "wppsi4_visuo_sp_raw", "wppsi4_fluid_res_raw", "wppsi4_working_mem_raw", "wppsi4_proc_speed_raw",
        "srs4_sum_subscales", "srs4_rrb_raw", "srs4_sci_raw",
        "MC_COUNT_TOTAL_FAILS_nooffails")


dat <- tb_merge[,vars]

## merge duplicated
dup_ap <- na.omit(dat$AP_ID)[na.omit(dat$AP_ID) %>% duplicated()]
dup_ap
# [1] "BIPP015" "BIPP016"

dup_ep <- na.omit(dat$Eprime_ID)[na.omit(dat$Eprime_ID) %>% duplicated()]
dup_ep
#  [1] 6324 9431 5832 6594 7153 6415
#  [7] 6652 7224 9429 7448 7654 7519
# [13] 7699 9406 7427 9444 7439 7207
# [19] 6470 6487 7238 7275 7687

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

# sex to 1=m, 2=f
dat_fin$sex <- ifelse(dat_fin$sex == "Male", 1, 2)


# check data types
summary(dat_fin)
# sapply(dat_fin[, c("WISC_VCI_CS", "WISC_PR_CS", "WISC_WM_CS", "WISC_PS_CS",
#                     "srs-rrb", "srs-sci")], as.numeric)
dat_fin$WISC_VCI_CS <- dat_fin$WISC_VCI_CS %>% as.numeric()
dat_fin$WISC_PR_CS <- dat_fin$WISC_PR_CS %>% as.numeric()
dat_fin$WISC_WM_CS <- dat_fin$WISC_WM_CS %>% as.numeric()

## "na" --> NA
dat_fin$WISC_PS_CS[!grepl("^\\d+$", dat_fin$WISC_PS_CS)][!dat_fin$WISC_PS_CS[!grepl("^\\d+$", dat_fin$WISC_PS_CS)] %>% is.na()] <- NA
dat_fin$WISC_PS_CS <- dat_fin$WISC_PS_CS %>% as.numeric()


## ">90" --> 90
dat_fin$`srs-rrb`[!grepl("^\\d+", dat_fin$`srs-rrb`)][!dat_fin$`srs-rrb`[!grepl("^\\d+$", dat_fin$`srs-rrb`)] %>% is.na()] <- "90"
dat_fin$`srs-sci`[!grepl("^\\d+", dat_fin$`srs-sci`)][!dat_fin$`srs-sci`[!grepl("^\\d+$", dat_fin$`srs-sci`)] %>% is.na()] <- "90"
dat_fin$`srs-tot`[!grepl("^\\d+", dat_fin$`srs-tot`)][!dat_fin$`srs-tot`[!grepl("^\\d+$", dat_fin$`srs-tot`)] %>% is.na()] <- "90"

dat_fin[, c("srs-tot", "srs-rrb", "srs-sci")] <- sapply(dat_fin[, c("srs-tot", "srs-rrb", "srs-sci")], as.numeric)


str(dat_fin)


# change colnames
colnames(dat_fin)[8:28] <- c("WISC VC 8yo", "WISC PR 8yo", "WISC WM 8yo", "WISC PS 8yo",
        "SRS Total 8yo", "SRS RRB 8yo", "SRS SCI 8yo",
        "Bayley Cog 22mo", "Bayley Lang 22mo", "Bayley Motor 22mo", "PARCA Cog 22mo", "PARCA Lang 22mo",
        "WPPSI VC 4yo", "WPPSI VS 4yo", "WPPSI FR 4yo", "WPPSI WM 4yo", "WPPSI PS 4yo",
        "SRS Total 4yo", "SRS RRB 4yo", "SRS SCI 4yo",
        "MCHAT Fails 22mo")


dat_fin <- dat_fin[,c(colnames(dat)[1:7], "WISC VC 8yo", "WISC PR 8yo", "WISC WM 8yo", "WISC PS 8yo",
        "SRS Total 8yo", "SRS RRB 8yo", "SRS SCI 8yo",
        "WPPSI VC 4yo", "WPPSI VS 4yo", "WPPSI FR 4yo", "WPPSI WM 4yo", "WPPSI PS 4yo",
        "SRS Total 4yo", "SRS RRB 4yo", "SRS SCI 4yo",
        "Bayley Cog 22mo", "Bayley Lang 22mo", "Bayley Motor 22mo", "PARCA Cog 22mo", "PARCA Lang 22mo", "MCHAT Fails 22mo")] 


# export
write.table(dat_fin, file = "../data/id_vars_fin.csv", sep = ",",
              quote = F, row.names = F, col.names = T)
