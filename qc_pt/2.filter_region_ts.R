#### =========== load time series =========== ####
### filtering subjects out by their mising regions
rm(list = ls())
setwd("/mnt/d/PROJECTS/preterm_language/qc_pt")

library(dplyr)
library(VIM)
# plots
# library(readxl)
# library(seqinr)       # col2alpha
# library(ggplot2)
# library(hexbin)
# library(RColorBrewer) # brewer.pal
# # stats / modelling
# library(matrixStats)  # colMaxs
# library(psych)        # fisherz
# library(nlme)         # lme
# library(mgcv)         # gam(m)
# library(vows)         # bs = 'bq'
# library(AICcmodavg)   # ?
# library(MuMIn)        # r.squaredGLMM
library(stringr)
# #functions
# source('rp.main.R')

#### path to save plots ####
plot.path <- "./plots/"

## load FD_df_after_exclusions
load("FD_df_after_exclusions.RData")

nt <- 400 #number of time series
nroi <- 376 # number of regions
ts.path <- "../data/timeseries_byglasserROIs/" # path to time series

all_files <- list.files(ts.path, pattern = "_timeseries_sept29.txt")
id_to_use <- str_split(all_files, "_", simplify = T)[, 1] %>% gsub("sub-", "", .) %in% FD_df_after_exclusions$AP_ID
files <- all_files[id_to_use]
length(files)
# 89

# check subjects with incorrect dimensions of time series
# exclude them
data <- lapply(files, function(file) {
  read.table(paste0(ts.path, file), stringsAsFactors = FALSE, sep = "", header = T)[, -c(1:2)]})
names(data) <- str_split(files, "_", simplify = T)[, 1]

wrong_dimen <- (unlist(lapply(data, nrow)) == 400 & unlist(lapply(data, ncol)) == 376)

# which subjects have incorrect dimensions
names(data)[!wrong_dimen]
# [1] "sub-AP068" "sub-AP099" "sub-AP132"


# check which regions are missing in these subjects
colninx <- paste0("Mean_", seq(1,376))
for (id in paste0("sub-", c("AP068", "AP099", "AP132"))) {
    ms_r <- setdiff(colninx, colnames(data[[id]])) %>% gsub("Mean_","",.) %>% str_c(., collapse = " ")
    print(paste0(id, ": missing region ", ms_r))
}
# [1] "sub-AP068: missing region 30 135 136 142 252 256 370 371 372 373 374 375 376"
# [1] "sub-AP099: missing region 136"
# [1] "sub-AP132: missing region 316"


# so keep AP099 and AP132 as they missed Hippocampus regions
#### regions 136 (left hippocampus), and 316 (right hippocampus) ####
# delete region 136 and 316 for all the subjects
ts_nohippo <- lapply(data, function(ts_df) {
  ts_df[,!(colnames(ts_df) %in% c("Mean_136", "Mean_316"))]
})

# delete AP068
ts_nohippo_exsbj <- ts_nohippo[names(ts_nohippo)!= "sub-AP068"]
# how many subjects left
length(ts_nohippo_exsbj)
# [1] 88
# dimensions of all TS, all have 374 regions
unlist(lapply(ts_nohippo_exsbj, ncol) == 374) %>% table()
# TRUE 
#  88 


## convert to 3d array ##
dim_order <- c(nt, ncol(ts_nohippo_exsbj[[1]]), length(ts_nohippo_exsbj))
# ts_final <- lm2a(ts_nohippo_exsbj, dim.order = dim_order, dimlab.list = "subjects")
ts_final<- array(NA, dim = dim_order)
for (i in 1:length(ts_nohippo_exsbj)) {
  ts_final[,,i] <- array(unlist(ts_nohippo_exsbj[[i]]), dim = dim_order[-3])
}

dim(ts_final)
# [1] 400 374  88


#### create final df ####
id_final <- gsub("sub-", "", names(ts_nohippo_exsbj))

#df_final in the correct ID order as the ts_final:
df_final <- FD_df_after_exclusions[id_final,]
dim(df_final)
# [1] 88  12


#### histograms AFTER FD and missing region subjects exclusions (limits require manual adjustment)  ####
# mean fd
mean_cutoff <- 0.11
pdf(paste0(plot.path,"hist_mean_fd_after_FD_", mean_cutoff ,"co_exclusion_FINAL.pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(df_final$fd.m,30, xlab=expression(paste(mu, ' FD (mm)',sep='')), ylab='frequency',main='',col='grey90',xlim=c(0,0.15))
abline(v=mean_cutoff,lty=2,col='grey') #lines for cutoffs - for exclusion 
dev.off()

# max fd
max_cutoff <- 0.6
pdf(paste0(plot.path,"hist_max_fd_after_FD_", max_cutoff ,"co_exclusion_FINAL.pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(df_final$fd.max,30, xlab=expression(paste('max FD (mm)',sep='')), ylab='frequency',main='',col='grey90',xlim=c(0,1))
abline(v=max_cutoff,lty=2,col='grey') #lines for cutoffs - for exclusion 
dev.off()


## save obejct for later scripts
save(ts_final, file = "ts_final.RData")
save(df_final, file = "df_final.RData")
write.table(df_final, "../data/pt_vars_fd.txt", quote = F)
#save.image("2.filter_region_ts.RData")


## KNN imputation to fill the NAs
# check which subjs has NAs
for (i in 1:(nrow(df_final))) {
  if (is.na(df_final[i,]) %>% sum() != 0) {
    print(paste0(row.names(df_final)[i], ": ", is.na(df_final[i,]) %>% sum()))
  }
}
# [1] "AP010: 2"
# [1] "AP011: 2"
# [1] "BIPP002: 2"
# [1] "BIPP003: 2"
# [1] "BIPP021: 2"

check_na <- cbind(
   lapply(
     lapply(df_final, is.na)
     , sum)
)
names(check_na[check_na[,1] != 0,])  
# [1] "WISC_FULL_CS" "WISC_PS_CS"   "srs.rrb"     
# [4] "srs.sci" 

# KNN
knn <- kNN(df_final, variable = names(check_na[check_na[,1] != 0,]), k=5, imp_var=F)
row.names(knn) <- row.names(df_final)
dim(knn)

# export
write.table(knn, "../data/pt_vars_fd_imp.txt", quote = F)




# # verify dimensions of text files
# ns <- nrow(FD_df_after_exclusions) #new number of subjects
# id_clean_FD <- row.names(FD_df_after_exclusions)
# for (i in 1:ns) {
#   if (i%%10==0) print(i)
#   temp = read.table(paste(ts.path,'sub-',id_clean_FD[i],'_glasserROI_timeseries_sept29.txt',sep=''), header=TRUE)[,-c(1:2)]
#   if (!(dim(temp)[1]==nt & dim(temp)[2]==nroi)) print(paste('dimensions for subj. ',toString(id_clean_FD[i]),' incorrect',sep=''))
# }
# rm(temp)

# [1] 10
# [1] 20
# [1] 30
# [1] "dimensions for subj. AP068 incorrect"
# [1] "dimensions for subj. AP099 incorrect"
# [1] 40
# [1] 50
# [1] "dimensions for subj. AP132 incorrect"
# [1] 60
# [1] 70


# load time series
# ts = array(NA,dim=c(nt,nroi,ns)) #3d array of all the time series
# for (i in 1:ns) {
#   if (i%%10==0) print(i)
#   ts[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_clean_FD[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi))
# }


#first exclude those with missing regions 
# #### load time series - excluding the 3 with a missing region
# # find the missing regions
# colninx <- paste0("Mean_", seq(1,376))
# for (id in c("AP068", "AP099", "AP132")) {
#     tmp <- read.table(paste(ts.path,'sub-', id,'_glasserROI_timeseries_sept29.txt',sep=''), header=TRUE)[,-c(1:2)]
#     ms_r <- setdiff(colninx, colnames(tmp)) %>% gsub("Mean_","",.) %>% str_c(., collapse = " ")
#     print(paste0(id, ": missing region ", ms_r))

# }

# [1] "AP068: missing region 30 135 136 142 252 256 370 371 372 373 374 375 376"
# [1] "AP099: missing region 136"
# [1] "AP132: missing region 316"



# #### time series - only of subjects with no NA regions ####
# id_clean_FD_clean_ts <- id_clean_FD[-which(id_clean_FD %in% c("AP068", "AP099", "AP132"))]
# id_clean_FD_clean_ts <- as.character(id_clean_FD_clean_ts)
# ns_no_NA = length(id_clean_FD_clean_ts)

# ## TS 3D array excluded 
# ts_no_NA_regions = array(NA,dim=c(nt,nroi,ns_no_NA)) #3d array of all the time series
# for (i in 1:ns_no_NA) {
#   if (i%%10==0) print(i)
#   ts_no_NA_regions[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_clean_FD_clean_ts[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi))
# }

# # we will exclude Hippocampus regions from glasser atlas in these participants 
# # as we already have a hippocampus parcellation from the subcortical parcellation
# #### regions 136 (left hippocampus), and 316 (right hippocampus) ####
# ts_no_NA_regions_noHipp <- ts_no_NA_regions[,-c(136, 316),] #exclude the glasser hippocampus regions
# dim(ts_no_NA_regions_noHipp)              
# # [1] 400 374 73
# # [1] 400 374  80

# #now, for those with regions missing, we will exclude them completely (n=3),
# #those with a hippocampus region missing, we will keep them, (n=4)
# #and remove the other hippocampus region
# #### load time series - excluding those with a missing region that isn't hippocampus ####


# #### ts for those with one hipp region missing: ####
# id_with_hippNAs_ts <- id_clean_FD[which(id_clean_FD %in% c("AP099", "AP132"))] # AP068 has too many regions missing
# id_with_hippNAs_ts <- as.character(id_with_hippNAs_ts)
# ns_hippNAs = length(id_with_hippNAs_ts)
# nroi_hippNA = 375 # they each have ONE hipp region missing (none have both missing)

# ts_hippNAs = array(NA,dim=c(nt,nroi_hippNA,ns_hippNAs)) #3d array of all the time series
# for (i in 1:ns_hippNAs) {
#   if (i%%10==0) print(i)
#   ts_hippNAs[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_with_hippNAs_ts[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi_hippNA))
# }

# #exclude the other hipp region:
# which(id_with_hippNAs_ts=="AP099") #1 - missing region 136
# which(id_with_hippNAs_ts=="AP132") #2 - missing region 316

# ts_hippNAs_1 <- ts_hippNAs[,-c(316), c(1)] #exclude right (316) - left (136) already missing
# ts_hippNAs_2 <- ts_hippNAs[,-c(135), c(2)]  #exclude left (136) -  right (316) already missing, 135 is becuz now the order and total n are changed


# #then merge the no_NA_noHipp ts with the NA_noHipp ts
# library(abind)
# ts.final <- abind(ts_hippNAs_1, ts_hippNAs_2, ts_no_NA_regions_noHipp)
# ts.final %>% dim
# [1] 400 374  75
## regions:376-2, subjects:76-3+2
# [1] 400 374  82
## subjects + 7

