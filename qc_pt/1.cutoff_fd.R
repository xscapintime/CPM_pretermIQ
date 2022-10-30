### ============ filtering subjects ================= ###
### delete subjects with outliner mean/max FD value
rm(list = ls())
setwd("/mnt/d/PROJECTS/preterm_language/qc")

#### load libraries  ####
library(dplyr)
# library(readxl)
# library(seqinr)       # col2alpha
# library(ggplot2)
# library(hexbin)
# library(RColorBrewer) # brewer.pal
# stats / modelling
library(matrixStats)  # colMaxs
# library(psych)        # fisherz
# library(nlme)         # lme
# library(mgcv)         # gam(m)
# library(vows)         # bs = "bq"
# library(AICcmodavg)   # ?
# library(MuMIn)        # r.squaredGLMM
# #functions
# source("rp.main.R")


#### path to save plots ####
plot.path <- "./plots/"

# #### demographic information - age, sex etc. (as vectors of length ns) ####
#### merged var from differnet spreadsheets ####
#stats <- read.csv("../data/id_vars_fin.csv")
stats <- read.csv("../data/vars_wisc_srs_8yo.csv")
stats <- stats[!is.na(stats$AP_ID),]
row.names(stats) <- stats$AP_ID

# only keep PT
stats <- stats %>% filter(group == "PT")
dim(stats)
# [1] 132 10
######## 132 PT subjects with outcomes/variables (NAs included) ######


#### framewise displacement - get from fmriprep .tsv file ####
### to check which subjects have maching fmri data ###
fd.path <- "../data/framwise_displacement/"

## get id
all_files <- list.files(fd.path, pattern = "_framwise_displacement.tsv")

# only keep PT
id <- sub("sub-", "", all_files) %>% sub("_framwise_displacement.tsv", "", .)
id <- id[id %in% stats$AP_ID]

files <- paste0("sub-", id, "_framwise_displacement.tsv")
length(files)
# [1] 100
###### 100 PT subjects with fmri data ######


#### set nt and ns ####
# not used
nt <- 400 #number of time series
ns <- length(files) #number of subjects

# fd <- matrix(NA, nrow=nt-1, ncol=ns)
# for (i in 1:ns) {
#   if (i%%1==0) print(i)
#     if (file.exists(paste0(fd.path, "sub-",id[i],"_framwise_displacement.tsv"))) {
#       fd[,i] = matrix(as.numeric(unlist(read.table(paste(fd.path,"sub-",id[i],"_framwise_displacement.tsv",sep = "")))))
#  }
# }

# check subjects with incorrect number of framewise diepalcement
# exclude them
data <- lapply(files, function(file) {
  read.csv(paste0(fd.path, file), stringsAsFactors = FALSE, sep = " ", header = F)})

wronglen <- unlist(lapply(data, nrow)) == 399

## which subj
id[!wronglen]
# [1] "AP049" "AP067"
data_keep <- data[wronglen]

# merge all time series to one dataframe
fd <- Reduce(cbind, data_keep)


dim(fd) # col = subjects
# [1] 399  98
### 2 subjects were filtered ###

colnames(fd) <- id[wronglen]


# omit na column
# no NA
fd <- t(na.omit(t(fd)))
dim(fd)
# [1] 399  98


#fd = array(NA,dim=c(nt-1,ns))
#for (i in 1:ns) {
#  fd[,i] = rowSums(abs(diff(mot[,4:6,i])))+(50*(pi/180))*rowSums(abs(diff(mot[,1:3,i])))
#}
fd.m <- colMeans(fd)   # mean fd for each subject
fd.max <- colMaxs(fd)  # max fd
names(fd.max) <- colnames(fd)

fivenum(fd.m)
fivenum(fd.max)

# ns = dim(fd)[2] # number of subjects
# nt = dim(fd)[1] + 1 ## 400 in this case


#### histograms before FD exclusions (limits require manual adjustment) ####
# mean fd
mean_cutoff <- 0.11
pdf(paste0(plot.path,"hist_mean_fd_cutoff_", mean_cutoff, ".pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg="white")
hist(fd.m,30, xlab=expression(paste(mu, " FD (mm)",sep="")), ylab="frequency",main="",col="grey90",xlim=c(0,0.2))
abline(v=mean_cutoff,lty=2,col="grey") #lines for cutoffs - for exclusion 
dev.off()

# max fd
max_cutoff <- 0.6
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
# [1]  9 10 37 42 53 54 65 12 87

row.names(FD_df)[toexclu]
# [1] "AP014"   "AP015"   "AP060"  
# [4] "AP069"   "AP111"   "AP112"  
# [7] "AP131"   "AP018"   "BIPP033"


#so now we exlclude IDs
#we exlclude 9 IDs - fd exclusions

#### data frame after FD
FD_df_after_exclusions <- FD_df[-toexclu,]
dim(FD_df_after_exclusions)
# [1] 89  2

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


FD_df_after_exclusions <- merge(stats, FD_df_after_exclusions, by = 0 )
row.names(FD_df_after_exclusions) <- FD_df_after_exclusions$AP_ID
FD_df_after_exclusions <- FD_df_after_exclusions[, -1]


## save 
save(FD_df_after_exclusions, file = "FD_df_after_exclusions.RData")
#save.image("1.cutoff_fd.RData")