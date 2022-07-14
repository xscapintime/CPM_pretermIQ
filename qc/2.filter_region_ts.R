#### =========== load time series =========== ####
### filtering subjects out by their mising regions
setwd("/mnt/d/PROJECTS/preterm_language/qc")

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
library(vows)         # bs = 'bq'
library(AICcmodavg)   # ?
library(MuMIn)        # r.squaredGLMM
library(stringr)
#functions
source('rp.main.R')

#### path to save plots ####
plot.path = "./plots/"

## load FD_df_after_exclusions
load("FD_df_after_exclusions.RData")

nt= 400 #number of time series
nroi = 376 # number of regions
ts.path = "./timeseries_byglasserROIs/" # path to time series


# verify dimensions of text files
ns = nrow(FD_df_after_exclusions) #new number of subjects
id_clean_FD <- as.character(FD_df_after_exclusions$AP_id)
for (i in 1:ns) {
  if (i%%10==0) print(i)
  temp = read.table(paste(ts.path,'sub-',id_clean_FD[i],'_glasserROI_timeseries_sept29.txt',sep=''), header=TRUE)[,-c(1:2)]
  if (!(dim(temp)[1]==nt & dim(temp)[2]==nroi)) print(paste('dimensions for subj. ',toString(id_clean_FD[i]),' incorrect',sep=''))
}
rm(temp)

# [1] 10
# [1] 20
# [1] 30
# [1] 40
# [1] "dimensions for subj. AP068 incorrect"
# [1] 50
# [1] 60
# [1] "dimensions for subj. AP099 incorrect"
# [1] 70
# [1] 80
# [1] 90
# [1] "dimensions for subj. AP132 incorrect"

# load time series
# ts = array(NA,dim=c(nt,nroi,ns)) #3d array of all the time series
# for (i in 1:ns) {
#   if (i%%10==0) print(i)
#   ts[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_clean_FD[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi))
# }


#first exclude those with missing regions 
#### load time series - excluding the 3 with a missing region
# find the missing regions
colninx <- paste0("Mean_", seq(1,376))
for (id in c("AP068", "AP099", "AP132")) {
    tmp <- read.table(paste(ts.path,'sub-', id,'_glasserROI_timeseries_sept29.txt',sep=''), header=TRUE)[,-c(1:2)]
    ms_r <- setdiff(colninx, colnames(tmp)) %>% gsub("Mean_","",.) %>% str_c(., collapse = " ")
    print(paste0(id, ": region ", ms_r))

}

# [1] "AP068: region 30 135 136 142 252 256 370 371 372 373 374 375 376"
# [1] "AP099: region 136"
# [1] "AP132: region 316"



#### time series - only of subjects with no NA regions ####
id_clean_FD_clean_ts <- id_clean_FD[- which(id_clean_FD %in% c("AP068", "AP099", "AP132"))]
id_clean_FD_clean_ts <- as.character(id_clean_FD_clean_ts)
ns_no_NA = length(id_clean_FD_clean_ts)

## TS 3D array excluded 
ts_no_NA_regions = array(NA,dim=c(nt,nroi,ns_no_NA)) #3d array of all the time series
for (i in 1:ns_no_NA) {
  if (i%%10==0) print(i)
  ts_no_NA_regions[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_clean_FD_clean_ts[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi))
}

# we will exclude Hippocampus regions from glasser atlas in these participants 
# as we already have a hippocampus parcellation from the subcortical parcellation
#### regions 136 (left hippocampus), and 316 (right hippocampus) ####
ts_no_NA_regions_noHipp <- ts_no_NA_regions[,-c(136, 316),] #exclude the glasser hippocampus regions


#now, for those with regions missing, we will exclude them completely (n=3),
#those with a hippocampus region missing, we will keep them, (n=4)
#and remove the other hippocampus region 
#### load time series - excluding those with a missing region that isn't hippocampus ####

# previously checked
# [1] "AP068: region 30 135 136 142 252 256 370 371 372 373 374 375 376"
# [1] "AP099: region 136"
# [1] "AP132: region 316"

#### ts for those with one hipp region missing: ####
id_with_hippNAs_ts <- id_clean_FD[which(id_clean_FD %in% c("AP099", "AP132"))] # AP068 has too many regions missing
id_with_hippNAs_ts <- as.character(id_with_hippNAs_ts)
ns_hippNAs = length(id_with_hippNAs_ts)
nroi_hippNA = 375 # they each have ONE hipp region missing (none have both missing)

ts_hippNAs = array(NA,dim=c(nt,nroi_hippNA,ns_hippNAs)) #3d array of all the time series
for (i in 1:ns_hippNAs) {
  if (i%%10==0) print(i)
  ts_hippNAs[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_with_hippNAs_ts[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi_hippNA))
}

#exclude the other hipp region:
which(id_with_hippNAs_ts=="AP099") #1 - missing region 136
which(id_with_hippNAs_ts=="AP132") #2 - missing region 316

ts_hippNAs_1 <- ts_hippNAs[,-c(136), c(1)] #exclude left (136) - right (316) already missing
ts_hippNAs_2 <- ts_hippNAs[,-c(315), c(2)]  #exclude right(316) - left (136) already missing, 315 is becuz now the order and total n are changed


#then merge the no_NA_noHipp ts with the NA_noHipp ts
library(abind)
ts.final <- abind(ts_hippNAs_1, ts_hippNAs_2, ts_no_NA_regions_noHipp)
ts.final %>% dim                     
# [1] 400 374  91 
## regions:376-2, subjects:107-15-3+2

#### create final df ####
id_final <- unlist(list(id_with_hippNAs_ts[1], id_with_hippNAs_ts[2],
                  id_clean_FD_clean_ts))

#df_final in the correct ID order as the ts_final:
FD_df_after_all_exclusions <- FD_df_after_exclusions[-which(id_clean_FD == "AP068"),] #final inclusions but wrong order
dim(FD_df_after_all_exclusions)
# [1] 91 11

df_final <- FD_df_after_all_exclusions %>%
  arrange(factor(FD_df_after_all_exclusions$AP_id, levels = id_final))


#### histograms AFTER FD and missing region subjects exclusions (limits require manual adjustment)  ####
# mean fd
mean_cutoff <- 0.1
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

save(df_final, file = "df_final.RData")
