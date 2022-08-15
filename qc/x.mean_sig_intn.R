#### THIS IS FOR LATER: mean signal intensities ####
### === no need for this for now === ###
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

### n of regions, ts, subjects...
nroi <- 376 # number of regions
ns_no_NA <- length(id_clean_FD_clean_ts)



#### 1.  now get average activation of each region for each subject
mean_intensity_per_region_per_subj <- array(NA, dim=c(nroi, ns_no_NA))
for (i in 1:ns_no_NA) {
  if (i%%10==0) print(i)
  mean_intensity_per_region_per_subj[,i] <- colMeans(ts_no_NA_regions[,,i]) 
}
mean_intensity_per_region_per_subj

#### 2. now get z-scores for each subject across regions
z_intensity_per_region_per_subj <- scale(mean_intensity_per_region_per_subj) ### NOT SURE 

#### 3. for each region: which subjects have regions with mean intensitiy < 2 S.D away from the mean (i.e. z < -1.96)
z_intensity_per_region_per_subj_2SD_away <- list(NA, length=nroi)
for (i in 1:nroi) {
  #if (i%%10==0) print(i)
  if(!sum(z_intensity_per_region_per_subj[i,] <=-1.96)==0)
  print(i)
  #z_intensity_per_region_per_subj_2SD_away[[i]] <- which(z_intensity_per_region_per_subj[i,] <=-1.96)
}
z_intensity_per_region_per_subj_2SD_away


#which regions for each subject have signal < 2SD 
z_intensity_per_region_per_subj_2SD_away_subjects <- list(NA, length=ns_no_NA)
for (i in 1:ns_no_NA) {
  if (i%%10==0) print(i)
  z_intensity_per_region_per_subj_2SD_away_subjects[[i]] <- which(z_intensity_per_region_per_subj[,i] <=-2.33)
}
z_intensity_per_region_per_subj_2SD_away_subjects

#print out regions with a scan that had LOW signal intensity
for (i in 1:nroi) {
  if(!sum(z_intensity_per_region_per_subj[i,] <=-2.33)==0)
    print(i)
}

#print out regions with a scan that had LOW signal intensity in > 5% of subjects
#can tweak the cut-off ; Z < -1.96 or -2.58, or Z <-2.81
(5*ns_no_NA)/100 #5% of subjs = 10.15= ~10 or 11 
(10*ns_no_NA)/100 #10% of subjs = 20.3 = ~20 or 21

for (i in 1:nroi) {
  if(sum(z_intensity_per_region_per_subj[i,] <=-2.33)>20)
    print(i)
}
