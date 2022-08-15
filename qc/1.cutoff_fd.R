### ============ filtering subjects ================= ###
### delete subjects with outliner mean/max FD value
rm(list = ls())
setwd("/mnt/d/PROJECTS/preterm_language/qc")

#### load libraries  ####
library(dplyr)
# plots
library(readxl)
library(seqinr)       # col2alpha
library(ggplot2)
library(hexbin)
library(RColorBrewer) # brewer.pal
# stats / modelling
library(matrixStats)  # colMaxs
library(psych)        # fisherz
library(nlme)         # lme
library(mgcv)         # gam(m)
library(vows)         # bs = "bq"
library(AICcmodavg)   # ?
library(MuMIn)        # r.squaredGLMM
#functions
source("rp.main.R")


#### path to save plots ####
plot.path = "./plots/"

# #### demographic information - age, sex etc. (as vectors of length ns) ####
# ID list
# id <- read_excel("../data/For Liyang - IDs resting state preprocessed.xlsx") ## AP ID
# stats <- read_excel("../data/ePrime_BIPP_master_file_GEORGE_LHv4 correct tmcq.xlsx")

## updated sheets
## can add vars later
# id <-  read_excel("../data/LH scanning database AP EAP BIPP.xlsx", sheet = "AP and BIPP QC_final", skip = 5)
# id <- id$AP_ID[grepl("AP|BIPP", id$AP_ID)]

# stats <- read_excel("../data/AP_marking_JULY 2022.xlsx", sheet = "Master", skip=1)

# stats <- stats %>% filter(AP_ID %in% id)



#### framewise displacement - get from fmriprep .tsv file ####
fd.path = "../data/framwise_displacement/"

## get id
files <- list.files(fd.path, pattern = "_framwise_displacement.tsv")
id <- sub("sub-", "", files) %>% sub("_framwise_displacement.tsv", "", .)

#### set nt and ns ####
nt = 400 #number of time series
ns = length(id) #number of subjects


fd <- matrix(NA, nrow=nt-1, ncol=ns)
for (i in 1:ns) {
  if (i%%1==0) print(i)
    if (file.exists(paste0(fd.path, "sub-",id[i],"_framwise_displacement.tsv"))) {
      fd[,i] = matrix(as.numeric(unlist(read.table(paste(fd.path,"sub-",id[i],"_framwise_displacement.tsv",sep = "")))))
 }
}

dim(fd) # col = subjects
colnames(fd) <- id

# omit na column
fd <- t(na.omit(t(fd)))
dim(fd)
# [1] 399 107

#fd = array(NA,dim=c(nt-1,ns))
#for (i in 1:ns) {
#  fd[,i] = rowSums(abs(diff(mot[,4:6,i])))+(50*(pi/180))*rowSums(abs(diff(mot[,1:3,i])))
#}
fd.m = colMeans(fd)   # mean fd for each subject
fd.max = colMaxs(fd)  # max fd

# ns = dim(fd)[2] # number of subjects
# nt = dim(fd)[1] + 1 ## 400 in this case


#### histograms before FD exclusions (limits require manual adjustment) ####
# mean fd
mean_cutoff <- 0.15
pdf(paste0(plot.path,"hist_mean_fd_cutoff_", mean_cutoff, ".pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg="white")
hist(fd.m,30, xlab=expression(paste(mu, " FD (mm)",sep="")), ylab="frequency",main="",col="grey90",xlim=c(0,0.2))
abline(v=mean_cutoff,lty=2,col="grey") #lines for cutoffs - for exclusion 
dev.off()

# max fd
max_cutoff <- 1
pdf(paste0(plot.path,"hist_max_fd_cutoff_",max_cutoff, ".pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg="white")
hist(fd.max,30, xlab=expression(paste("max FD (mm)",sep="")), ylab="frequency",main="",col="grey90",xlim=c(0,1.5))
abline(v=max_cutoff,lty=2,col="grey") #lines for cutoffs - for exclusion 
dev.off()

#### exclude FD max and means beyond cut-off ####
#find Ids and exclude.. 
# add vars later
# FD_df <- stats %>% select(AP_id, sex, group, ga, birth_weight, bayley22_language_comp, parca22_language,AP_id,
#                     wppsi4_verb_compr_raw, wisc8_vci_cs) %>%
#                     cbind(., fd.m, fd.max)
# rownames(FD_df) <- FD_df$AP_id

FD_df <- cbind(fd.m, fd.max) %>% as.data.frame()

toexclu <- unique(c(which(fd.m >= mean_cutoff), which(fd.max >= max_cutoff)))
toexclu                       
# [1] 10 23

#so now we exlclude IDs
#we exlclude 2 IDs - fd exclusions

#### data frame after FD
FD_df_after_exclusions <- FD_df[-toexclu,]
dim(FD_df_after_exclusions)
# [1] 105   2

#### histograms AFTER FD exclusions (limits require manual adjustment) ####
# mean fd
pdf(paste0(plot.path,"hist_mean_fd_after_FD_", mean_cutoff, "co_exclusion.pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg="white")
hist(FD_df_after_exclusions$fd.m,30, xlab=expression(paste(mu, " FD (mm)",sep="")), ylab="frequency",main="",col="grey90",xlim=c(0,0.15))
abline(v=mean_cutoff, lty=2,col="grey") #lines for cutoffs - for exclusion 
dev.off()

# max fd
pdf(paste0(plot.path,"hist_max_fd_after_FD_", max_cutoff, "co_exclusions.pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg="white")
hist(FD_df_after_exclusions$fd.max,30, xlab=expression(paste("max FD (mm)",sep="")), ylab="frequency",main="",col="grey90",xlim=c(0,1))
abline(v=max_cutoff, lty=2,col="grey") #lines for cutoffs - for exclusion 
dev.off()


## save 
save(FD_df_after_exclusions, file = "FD_df_after_exclusions.RData")
save.image("1.cutoff_fd.RData")