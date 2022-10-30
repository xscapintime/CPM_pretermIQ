setwd("/mnt/d/PROJECTS/preterm_language/qc")

#### install libraries ####
# plots
# install.packages("readxl")
# install.packages("seqinr")       # col2alpha
# install.packages("ggplot2")
# install.packages("hexbin")
# install.packages("RColorBrewer") # brewer.pal
# # stats / modelling
# install.packages("matrixStats")  # colMaxs
# install.packages("psych")        # fisherz
# install.packages("nlme")         # lme
# install.packages("mgcv")         # gam(m)
# install.packages("vows")         # bs = 'bq'
# install.packages("AICcmodavg")   # ?
# install.packages("MuMIn")        # r.squaredGLMM
# install.packages("dplyr")


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
library(vows)         # bs = 'bq'
library(AICcmodavg)   # ?
library(MuMIn)        # r.squaredGLMM
#functions
source('rp.main.R')


#### path to save plots ####
plot.path = './plots/'

# #### demographic information - age, sex etc. (as vectors of length ns) ####
# ID list
id <- read_excel('For Liyang - IDs resting state preprocessed.xlsx')
id <- id$IDs

#### set nt and ns ####
nt= 400 #number of time series
ns= length(id) #number of subjects - 220

#### framewise displacement - get from fmriprep .tsv file ####
fd.path = "./framwise_displacement/"

fd=matrix(NA, nrow=nt-1, ncol=ns)
for (i in 1:ns) {
  if (i%%1==0) print(i)
    fd[,i] = matrix(as.numeric(unlist(read.table(paste(fd.path,'sub-',id[i],'_framwise_displacement.tsv',sep = '')))))
}

dim(fd) # col = subjects
# colnames(fd) <- id


#fd = array(NA,dim=c(nt-1,ns))
#for (i in 1:ns) {
#  fd[,i] = rowSums(abs(diff(mot[,4:6,i])))+(50*(pi/180))*rowSums(abs(diff(mot[,1:3,i])))
#}
fd.m = colMeans(fd)   # mean fd for each subject
fd.max = colMaxs(fd)  # max fd
ns = dim(fd)[2] # number of subjects
nt = dim(fd)[1] + 1 ## 400 in this case


#### histograms before FD exclusions (limits require manual adjustment) ####
# mean fd
mean_cutoff <- 0.13
pdf(paste0(plot.path,"hist_mean_fd_cutoff_", mean_cutoff, ".pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(fd.m,30, xlab=expression(paste(mu, ' FD (mm)',sep='')), ylab='frequency',main='',col='grey90',xlim=c(0,0.8))
abline(v=mean_cutoff,lty=2,col='grey') #lines for cutoffs - for exclusion 
dev.off()

# max fd
max_cutoff <- 1
pdf(paste0(plot.path,"hist_max_fd_cutoff_",max_cutoff, ".pdf"),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(fd.max,30, xlab=expression(paste('max FD (mm)',sep='')), ylab='frequency',main='',col='grey90',xlim=c(0,17))
abline(v=max_cutoff,lty=2,col='grey') #lines for cutoffs - for exclusion 
dev.off()

#### exclude FD max and means beyond cut-off ####
#find Ids and exclude.. 
FD_df <- as.data.frame(cbind(id)) #, age, male, fd.m, fd.max, group, surveyno))
rownames(FD_df) <- FD_df$id

which(fd.m >= 0.13) # 23
which(fd.max >= 1) # 11 23

#so now we exlclude IDs 4  29  60  63  94 130 171 181 199 208
#we exlclude 10 IDs - fd exclusions
fd_exclusions <- c(4,  29,  60,  63,  94, 130, 171, 181, 199, 208)
FD_df[fd_exclusions,]

#repeated subjects exclusions:
repeated_subjects <- c("UCCHILDHB061", "UCCHILDG19", "UCCHILDHB106","EUCHILDG02" )
FD_df[repeated_subjects,]
which(FD_df==c("UCCHILDHB061")) #86
which(FD_df==c("EUCHILDG02")) #1
repeated_subjects_exclusions <- c(1,86) #chose to exclude the scan with higher FD mean -since qualitatively both seem fine

#### data frame after FD and repeated subject exclusions ####
FD_df_after_exclusions <- FD_df[-c(cbind(fd_exclusions,repeated_subjects_exclusions)),]
FD_df_after_exclusions$fd.max <- as.numeric(as.character(FD_df_after_exclusions$fd.max))
FD_df_after_exclusions$fd.m <- as.numeric(as.character(FD_df_after_exclusions$fd.m))
FD_df_after_exclusions$age <- as.numeric(as.character(FD_df_after_exclusions$age))
FD_df_after_exclusions$id <- as.character(FD_df_after_exclusions$id)
FD_df_after_exclusions$surveyno <- as.character(FD_df_after_exclusions$surveyno)


#### histograms AFTER FD exclusions (limits require manual adjustment) ####
# mean fd
pdf(paste(plot.path,'hist_mean_fd_after_FD_0.4co_exclusion.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(FD_df_after_exclusions$fd.m,30, xlab=expression(paste(mu, ' FD (mm)',sep='')), ylab='frequency',main='',col='grey90',xlim=c(0,0.5))
abline(v=0.4,lty=2,col='grey') #lines for cutoffs - for exclusion 
dev.off()
# max fd
pdf(paste(plot.path,'hist_max_fd_after_FD_4co_exclusions.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(FD_df_after_exclusions$fd.max,30, xlab=expression(paste('max FD (mm)',sep='')), ylab='frequency',main='',col='grey90',xlim=c(0,5))
abline(v=4,lty=2,col='grey') #lines for cutoffs - for exclusion 
dev.off()

#### load time series ####
nroi = 376 # number of regions
ts.path = "~/OneDrive/BIPP/Projects/Functional connectivity/timeseries_byglasserROIs/" # path to time series

# verify dimensions of text files
ns = nrow(FD_df_after_exclusions) #new number of subjects
id_clean_FD <- as.character(FD_df_after_exclusions$id)
for (i in 1:ns) {
  if (i%%10==0) print(i)
  temp = read.table(paste(ts.path,'sub-',id_clean_FD[i],'_glasserROI_timeseries_sept29.txt',sep=''), header=TRUE)[,-c(1:2)]
  if (!(dim(temp)[1]==nt & dim(temp)[2]==nroi)) print(paste('dimensions for subj. ',toString(id_clean_FD[i]),' incorrect',sep=''))
}
rm(temp)

# load time series
ts= array(NA,dim=c(nt,nroi,ns)) #3d array of all the time series
for (i in 1:ns) {
  if (i%%10==0) print(i)
  ts[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_clean_FD[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi))
}

#### low mean intensity exclusions - to use LATER in a more conservative analysis ####
# will first run a liberal analysis - including all - then check with the conservative analysis - excluding those with low signal intensity
# from previous paper: "30 cortical regions (near frontal and temporal poles) were excluded due to low regional mean signal, defined by a low Z score of mean signal intensity (Z <  1.96) in at least one scan "
# mean z-score ....
#exclude those at < 2SD
#first exclude those with missing regions 
#### load time series - excluding the 7 with a missing region
which(id_clean_FD=="UCCHILDHB009") #37 - missing region 136
which(id_clean_FD=="UCCHILDHB074") #88 - missing region 74
which(id_clean_FD=="UCCHILDHB076") #90 - missing region 75
which(id_clean_FD=="UCCHILDHB093") #106 - missing region 316
which(id_clean_FD=="UCCHILDJB011") #130 - missing region 316
which(id_clean_FD=="UCCHILDJB074") #182 - missing region 136
which(id_clean_FD=="UCCHILDJB083") #190 - missing region 230

#### time series - only of subjects with no NA regions ####
id_clean_FD_clean_ts <- id_clean_FD[-c(37, 88, 90, 106, 130, 182, 190)]
id_clean_FD_clean_ts <- as.character(id_clean_FD_clean_ts)
ns_no_NA = length(id_clean_FD_clean_ts)
ts_no_NA_regions = array(NA,dim=c(nt,nroi,ns_no_NA)) #3d array of all the time series
for (i in 1:ns_no_NA) {
  if (i%%10==0) print(i)
  ts_no_NA_regions[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_clean_FD_clean_ts[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi))
}

# we will exclude Hippocampus regions from glasser atlas in these participants 
# as we already have a hippocampus parcellation from the subcortical parcellation
# regions 136 (left hippocampus), and 316 (right hippocampus)
ts_no_NA_regions_noHipp <- ts_no_NA_regions[,-c(136, 316),] #exclude the glasser hippocampus regions

#now, for those with regions missing, we will exclude them completely (n=3),
#those with a hippocampus region missing, we will keep them, (n=4)
#and remove the other hippocampus region 
#### load time series - excluding those with a missing region that isn't hippocampus ####
which(id_clean_FD=="UCCHILDHB009") #37 - missing region 136
which(id_clean_FD=="UCCHILDHB074") #88 - missing region 74 - exclude completely
which(id_clean_FD=="UCCHILDHB076") #90 - missing region 75 - exclude completely
which(id_clean_FD=="UCCHILDHB093") #106 - missing region 316
which(id_clean_FD=="UCCHILDJB011") #130 - missing region 316
which(id_clean_FD=="UCCHILDJB074") #182 - missing region 136
which(id_clean_FD=="UCCHILDJB083") #190 - missing region 230 - exclude completely

#### ts for those with one hipp region missing: ####
id_with_hippNAs_ts <- id_clean_FD[c(37, 106, 130, 182)]
id_with_hippNAs_ts <- as.character(id_with_hippNAs_ts)
ns_hippNAs = length(id_with_hippNAs_ts)
nroi_hippNA = 375 # they each have ONE hipp region missing (none have both missing)
ts_hippNAs = array(NA,dim=c(nt,nroi_hippNA,ns_hippNAs)) #3d array of all the time series
for (i in 1:ns_hippNAs) {
  if (i%%10==0) print(i)
  ts_hippNAs[,,i] = array(as.numeric(as.character(unlist(read.delim(paste(ts.path,'sub-',id_with_hippNAs_ts[i],'_glasserROI_timeseries_sept29.txt',sep=''))[,-c(1:2)]))),dim=c(nt,nroi_hippNA))
}

#exclude the other hipp region:
which(id_with_hippNAs_ts=="UCCHILDHB009") #1 - missing region 136
which(id_with_hippNAs_ts=="UCCHILDHB093") #2 - missing region 316
which(id_with_hippNAs_ts=="UCCHILDJB011") #3 - missing region 316
which(id_with_hippNAs_ts=="UCCHILDJB074") #4 - missing region 136

ts_hippNAs_2.3 <- ts_hippNAs[,-c(136), c(2,3)] #exclude left (136) - right (316) already missing
ts_hippNAs_1.4 <- ts_hippNAs[,-c(315), c(1,4)]  #exclude right(316) - left (136) already missing 

#then merge the no_NA_noHipp ts with the NA_noHipp ts
library(abind)
ts.final <- abind(ts_hippNAs_1.4, ts_hippNAs_2.3, ts_no_NA_regions_noHipp)

#### create final df ####
id_final <- unlist(list(id_with_hippNAs_ts[1], id_with_hippNAs_ts[4], 
                  id_with_hippNAs_ts[2], id_with_hippNAs_ts[3],
                  id_clean_FD_clean_ts))

which(FD_df_after_exclusions$id=="UCCHILDHB074") #88
which(FD_df_after_exclusions$id=="UCCHILDHB076") #90
which(FD_df_after_exclusions$id=="UCCHILDJB083") #190
FD_df_after_all_exclusions <- FD_df_after_exclusions[-c(88,90,190),] #final inclusions but wrong order
#df_final in the correct ID order as the ts_final:
df_final <-   FD_df_after_all_exclusions %>%
  arrange(factor(FD_df_after_all_exclusions$id, levels = id_final))

#### histograms AFTER FD and missing region subjects exclusions (limits require manual adjustment)  ####
# mean fd
pdf(paste(plot.path,'hist_mean_fd_after_FD_0.4co_exclusion_FINAL.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(df_final$fd.m,30, xlab=expression(paste(mu, ' FD (mm)',sep='')), ylab='frequency',main='',col='grey90',xlim=c(0,0.5))
abline(v=0.4,lty=2,col='grey') #lines for cutoffs - for exclusion 
dev.off()
# max fd
pdf(paste(plot.path,'hist_max_fd_after_FD_4co_exclusions_FINAL.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
hist(df_final$fd.max,30, xlab=expression(paste('max FD (mm)',sep='')), ylab='frequency',main='',col='grey90',xlim=c(0,5))
abline(v=4,lty=2,col='grey') #lines for cutoffs - for exclusion 
dev.off()

#### THIS IS FOR LATER: mean signal intensities ####
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


#### loop to create fc ####
nroi.f = dim(ts.final)[2]  #new number of rois -previous nroi (initial) = nroi
ns.f = dim(ts.final)[3] 
fc.raw = array(NA,dim=c(nroi.f,nroi.f,ns.f))  # raw FC - before fisher transformation
for (i in 1:ns.f) {
  if (i%%10==0) print(i)
  fc.raw[,,i] = cor(ts.final[,,i])
}

# fisher transform correlations
fc = fisherz(fc.raw)

# set matrix diagonals to NaN
for (n in 1:ns.f) {
  fc[,,n][seq(from=1,to=(nroi.f^2),by=(nroi.f+1))] = NaN
}

#### subjs 76: "UCCHILDHB050" and 106:"UCCHILDHB090" have almost competely missing brains in rsFMRI ####
nasubjs.in.fc.conn.raw = list(NA)

for (i in 1:nroi.f) {
  print(i)
  for (j in 1:nroi.f) {
    nasubjs.in.fc.conn.raw[[i]]=which(is.na(fc[i,j,]))
  }
}

#### visualise FC distributions (raw / fisher Z) ####
#can detect outliers here

#set up density plot parameters
triup = upper.tri(matrix(nrow=nroi.f,ncol=nroi.f))

dens.n = 500                                            # number of points in distribution
dens.i = -1.5; dens.f = 2.5                                 # initial and final values for distributions - determine by looking at the range on fc
abl = c(-1.5,0.5,1.5,2.5);                                      # horizontal lines on plots - also depends on your range
dens.x = seq(from=dens.i,to=dens.f,length.out = dens.n) # x-axis points at which distribution is sampled
dens.y.raw = dens.y = matrix(NA, nrow=ns.f, ncol=dens.n)  # initialise variable for storing values of distributions

# densities of correlation distribution
for (i in 1:ns.f) {
  if (i%%10 == 0) print(i)
  dens.y.raw[i,] = density(fc.raw[,,i][triup],n=dens.n,from=dens.i,to=dens.f, na.rm = TRUE)$y # raw
  dens.y[i,] = density(fc[,,i][triup],n=dens.n,from=dens.i,to=dens.f, na.rm = TRUE)$y         # fisher transformed
}

# identify ranges of values, to fix axes
summary(c(dens.y.raw,dens.y)) # to set ylim for the plot:
age_final = df_final$age
survey_final=df_final$surveyno
fd.m_final=df_final$fd.m
group_final=df_final$group
male_final = df_final$male

# plots
#creates a page for each subject
pdf(paste(plot.path,'ind_fc_dist.pdf',sep=''),width=6,height=5)
for (n in 1:ns.f) {
  if (n%%10 == 0) print(n)
  par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 1.65, cex.axis = 1.5, cex.main = 1.5, font.main = 1, bg='white')
  plot(0, type='n', ylab='density', xlab='correlation', main=bquote(ID~.(survey_final[n])*'  |  age = '*.(signif(age_final[n],3))*' y  |  '*mu*'FD = '*.(round(fd.m_final[n],2))*' mm'), ylim=c(0,3), xlim=c(dens.i,dens.f)) 
  abline(v=abl,col='grey',lty=2)
  lines(dens.x,dens.y.raw[n,],col='steelblue',lwd=2)
  lines(dens.x,dens.y[n,],col='darkorange',lwd=2)
  legend('topright',c('raw','fisher Z'),col=c('steelblue','darkorange'),lty=1,lwd=3,bty='n')
}
dev.off()

# coordinates of ROI centroids (dimension nroi x 3 (= x, y, z))
#pos = read.table("/Users/lailahadaya/OneDrive/BIPP/Projects/Functional connectivity/MMP_w_FreesSurfer_subcort_coords_MNI.txt")
#pos = read.table("~/OneDrive/BIPP/Projects/Functional connectivity/MMP_w_FreesSurfer_subcort_coords_MNI.txt")
pos = MMP_MNI_376_UCCHILD_coords <- read.table("~/OneDrive/BIPP/Projects/Functional connectivity/MMP_MNI_376_UCCHILD_coords.txt", quote="\"", comment.char="")

# matrix of euclidean distance between regions
dist = matrix(nrow=nroi.f,ncol=nroi.f)
for (i in 1:nroi.f) {; for (j in 1:nroi.f) {; dist[i,j] = sqrt( sum( (pos[i,]-pos[j,])^2 ) ); }; }

# s = sort(age,decreasing=F,index.return=T)
# ix = s$ix

# indices of upper triangular part of matrix
triup = upper.tri(matrix(nrow=nroi.f,ncol=nroi.f))
# mean FC
fc.m_final = vector(length=ns.f)
for (n in 1:ns.f) fc.m_final[n] = mean(fc[,,n][triup],na.rm=T)

#### plot violin and boxplots of mean FD and max FD by group - run t-tests ####
#http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
library(ggplot2)
# Basic violin plot
pdf(paste(plot.path,'violinplot_fdm_bygroup.pdf',sep=''),width=6,height=5)
violin_fd.m <- ggplot(df_final, aes(x=group, y=fd.m), trim=FALSE) + 
  geom_violin()
violin_fd.m + geom_boxplot(width=0.1) #add boxplot 
dev.off()

pdf(paste(plot.path,'violinplot_fdmax_bygroup.pdf',sep=''),width=6,height=5)
violin_fd.max <- ggplot(df_final, aes(x=group, y=fd.max), trim=FALSE) + 
  geom_violin()
violin_fd.max + geom_boxplot(width=0.1) #add boxplot 
dev.off()

#### assumptions ####
library(lawstat)
levene.test(df_final$fd.m, df_final$group, location = c("mean"))
levene.test(df_final$fd.max, df_final$group, location = c("mean"))

shapiro.test(df_final$fd.m)
shapiro.test(df_final$fd.max)

#### Boxplots  ####
library(ggstatsplot)
pdf(paste(plot.path,'fdm_ttest_bygroup.pdf',sep=''),width=6,height=5)
ggstatsplot::ggbetweenstats(
  plot.type = "box",
  data = df_final,
  x = group,
  xlab="Group",
  ylab="fd.m",
  y = fd.m,
  title = "Mean fd in FT vs VPT",
  type = "nonparametric",
  pairwise.comparisons = TRUE,
  pairwise.annotation = "p.value",
  pairwise.display = "significant",
  p.adjust.method="fdr",
  var.equal=TRUE,
  notch=FALSE,
  messages = FALSE,
  ggstatsplot.layer = FALSE,
)
dev.off()

pdf(paste(plot.path,'fdmax_ttest_bygroup.pdf',sep=''),width=6,height=5)
ggstatsplot::ggbetweenstats(
  plot.type = "box",
  data = df_final,
  x = group,
  xlab="Group",
  ylab="fd.max",
  y = fd.max,
  title = "fd max in FT vs VPT",
  type = "nonparametric",
  pairwise.comparisons = TRUE,
  pairwise.annotation = "p.value",
  pairwise.display = "significant",
  p.adjust.method="fdr",
  var.equal=TRUE,
  notch=FALSE,
  messages = FALSE,
  ggstatsplot.layer = FALSE,
)
dev.off()

#### violin plots: ####
pdf(paste(plot.path,'violin_fdm_ttest_bygroup.pdf',sep=''),width=6,height=5)
ggstatsplot::ggbetweenstats(
  plot.type = "violin",
  data = df_final,
  x = group,
  xlab="Group",
  ylab="fd.m",
  y = fd.m,
  title = "mean FD in FT vs VPT",
  type = "nonparametric",
  pairwise.comparisons = TRUE,
  pairwise.annotation = "p.value",
  pairwise.display = "significant",
  p.adjust.method="fdr",
  var.equal=TRUE,
  notch=FALSE,
  messages = FALSE,
  ggstatsplot.layer = FALSE,
)
dev.off()

pdf(paste(plot.path,'violin_fdmax_ttest_bygroup.pdf',sep=''),width=6,height=5)
ggstatsplot::ggbetweenstats(
  plot.type = "violin",
  data = df_final,
  x = group,
  xlab="Group",
  ylab="fd.max",
  y = fd.max,
  title = "fd max in FT vs VPT",
  type = "nonparametric",
  pairwise.comparisons = TRUE,
  pairwise.annotation = "p.value",
  pairwise.display = "significant",
  p.adjust.method="fdr",
  var.equal=TRUE,
  notch=FALSE,
  messages = FALSE,
  ggstatsplot.layer = FALSE,
)
dev.off()



#### mean FC as a function of mean FD ####
#relationship between FD and mean FC
#here we want to plot the inxn by group 
#### plot interaction 
#https://philippmasur.de/2018/11/26/visualizing-interaction-effects/
library(sjPlot)
library(sjmisc)
library(ggplot2)
theme_set(theme_sjplot())
# Setting up the building blocks
basic_plot <- ggplot(df_final ,
                     aes(x = fd.m_final,
                         y = fc.m_final,
                         color = group_final)) +
  theme_bw() +
  labs(x = "fd.m",
       y = "fc.m",
       color = "group")
# Colored scatterplot
basic_plot +
  geom_point()
# Colored scatterplot and regression lines
pdf(paste(plot.path,'fdm_fcm_inxn_bygroup.pdf',sep=''),width=6,height=5)
basic_plot +
  geom_point(alpha = .3, 
             size = .9) +
  geom_smooth(method = "lm") 
dev.off()

#l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
l = lm(y ~ x1 + x2 + x1*x2, data=df)
newdat = expand.grid(x2='male',x1=range(fd.m))
#pred.m = predictSE.lme(l, newdata=expand.grid(x2='male',x1=range(fd.m)), se.fit=T, level=0) #predicting the error bars - CIs for the lme 
#pred.f = predictSE.lme(l, newdata=expand.grid(x2='female',x1=range(fd.m)), se.fit=T, level=0)
#pred = list(); pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
pdf(paste(plot.path,'fd_m_fc_raw_m_lme.pdf',sep=''),width=6,height=5)
ggplot(df, aes(x=x1, y=y)) +
  labs(x = expression(paste(mu, ' FD (mm)',sep='')), y = expression(paste(mu, ' FC',sep='')), title = bquote(r^2*' = '*.(toString(signif(r.squaredGLMM(l)[1,1],2)))*', '*p[FD]*' = '*.(toString(signif(summary(l)$tTable[2,5],2))))) +
  geom_point(size=1, colour = 'grey') +
  #geom_path(aes(group = id), colour = 'grey') + # spaghetti plot
  geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
  geom_line(data=newdat, aes(y=pred$fit),size=1) +
  theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#plot the interaction by group
# look at whether movement is related to mean FC - for each group seperately then see if there's an interaction
df_final_with_fcmean <- cbind(df_final,fc.m_final)
pdf(paste(plot.path,'fd_m_fc_m_lm_CONTROLS.pdf',sep=''),width=6,height=5)
ggplot(df_final_with_fcmean[df_final_with_fcmean$group==0,], aes(x=fd.m, y=fc.m_final)) +
  theme_bw() +
  geom_point(alpha = .3, 
             size = .9) +
  geom_smooth(method='lm', formula= y~x)+
  ggtitle("FD FC in full term controls")
dev.off()

pdf(paste(plot.path,'fd_m_fc_m_lm_PRETERMS.pdf',sep=''),width=6,height=5)
ggplot(df_final_with_fcmean[df_final_with_fcmean$group==1,], aes(x=fd.m, y=fc.m_final)) +
  theme_bw() +
  geom_point(alpha = .3, 
             size = .9) +
  geom_smooth(method='lm', formula= y~x)+
  ggtitle("FD FC in very preterms")
dev.off()

library(Hmisc)
#cor between FD and FC in controls
rcorr(df_final_with_fcmean[df_final_with_fcmean$group==0,"fd.m"], df_final_with_fcmean[df_final_with_fcmean$group==0,"fc.m_final"], type = "spearman")
#cor between FD and FC in preterms
rcorr(df_final_with_fcmean[df_final_with_fcmean$group==1,"fd.m"], df_final_with_fcmean[df_final_with_fcmean$group==1,"fc.m_final"], type = "spearman")

### another way of plotting this:
#plot(fd.m,fc.m,xlab='...',ylab='...')
fd.pred = seq(from=min(fd.m),to=max(fd.m),length.out=100)
pred = predict(l,newdata=data.frame('x1'=fd.pred),interval='confidence')
#polygon(c(rev(fd.pred), fd.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
#lines(fd.pred,pred[,1],col='black',lwd=3);

#inxn
summary(lm(fc.m_final~fd.m_final + group_final + fd.m_final*group_final))



#### "Satterthwaite plot" - edge-wise corr(FC,FD) (across subjects) as a function of distance ####
# set-up colorbar for hexbin
color.bar <- function(lut, alpha, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  lut = sapply(lut, col2alpha, alpha)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, at=ticks, labels=ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# correlate connectivity and motion (fd) for each edge
conn.raw.vs.mot = array(NA,dim=c(nroi.f,nroi.f))
for (i in 1:nroi.f) {
  print(i)
  for (j in 1:nroi.f) {
    conn.raw.vs.mot[i,j] = cor(fc[i,j,],fd.m_final,method='pearson')
  }
}

##hexbin
# set-up
bin = hexbin(dist[triup], conn.raw.vs.mot[triup], xbins=55)
my_colors = colorRampPalette(rev(brewer.pal(11,'Spectral')))
# plot
pdf(paste(plot.path,'satt_plot_raw_hexbin.pdf',sep=''),width=6,height=5)
par(mar=c(1, 1, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
x = dist[triup]; l = lm(conn.raw.vs.mot[triup]~x); #abline(l,col='orange',lwd=1.5)
spear = cor.test(dist[triup],conn.raw.vs.mot[triup],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
pl = plot(bin, colramp=my_colors , legend=F,xlab='distance (mm)', ylab='corr. of FC and FD (r)',main=rp.main.sp(sp.rho,sp.p,2)) 
hexVP.abline(pl$plot.vp,l,col = 'white', lwd = 5)
hexVP.abline(pl$plot.vp,l,col = gray(.6), lwd = 2)
dev.off()

# plot colorbar
pdf(paste(plot.path,'satt_plot_raw_hexbin_colorbar.pdf',sep=''),width=1.5,height=6)
color.bar(rev(brewer.pal(11,'Spectral')), alpha = 0.8, min = min(bin@count), max =  max(bin@count), nticks = 6, ticks = rev(seq(min(bin@count), max(bin@count), len=3)))
dev.off()

# Sat plot for each group
controls_in_finaldf <- which(df_final$group==0)
VPTs_in_finaldf <- which(df_final$group==1)

fc_controls <- fc[,,c(controls_in_finaldf)]
fc_VPTs <-fc[,,c(VPTs_in_finaldf)]

#### correlate connectivity and motion (fd) for each edge ####
#for controls:
conn.raw.vs.mot_CONTROLS = array(NA,dim=c(nroi.f,nroi.f))
for (i in 1:nroi.f) {
  print(i)
  for (j in 1:nroi.f) {
    conn.raw.vs.mot_CONTROLS[i,j] = cor(fc_controls[i,j,],fd.m_final[controls_in_finaldf],method='pearson')
  }
}
#for VPTs:
conn.raw.vs.mot_VPTS = array(NA,dim=c(nroi.f,nroi.f))
for (i in 1:nroi.f) {
  print(i)
  for (j in 1:nroi.f) {
    conn.raw.vs.mot_VPTS[i,j] = cor(fc_VPTs[i,j,],fd.m_final[VPTs_in_finaldf],method='pearson')
  }
}

##hexbin
# set-up for controls:
bin = hexbin(dist[triup], conn.raw.vs.mot_CONTROLS[triup], xbins=55)
my_colors = colorRampPalette(rev(brewer.pal(11,'Spectral')))
# plot for controls
pdf(paste(plot.path,'satt_plot_raw_hexbin_CONTROLS.pdf',sep=''),width=6,height=5)
par(mar=c(1, 1, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
x = dist[triup]; l = lm(conn.raw.vs.mot_CONTROLS[triup]~x); #abline(l,col='orange',lwd=1.5)
spear = cor.test(dist[triup],conn.raw.vs.mot_CONTROLS[triup],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
pl = plot(bin, colramp=my_colors , legend=F,xlab='distance (mm)', ylab='CONTROLS corr. of FC and FD (r)',main=rp.main.sp(sp.rho,sp.p,2)) 
hexVP.abline(pl$plot.vp,l,col = 'white', lwd = 5)
hexVP.abline(pl$plot.vp,l,col = gray(.6), lwd = 2)
dev.off()

# set-up for VPTs:
bin = hexbin(dist[triup], conn.raw.vs.mot_VPTS[triup], xbins=55)
my_colors = colorRampPalette(rev(brewer.pal(11,'Spectral')))
# plot for VPTs
pdf(paste(plot.path,'satt_plot_raw_hexbin_VPTs.pdf',sep=''),width=6,height=5)
par(mar=c(1, 1, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
x = dist[triup]; l = lm(conn.raw.vs.mot_VPTS[triup]~x); #abline(l,col='orange',lwd=1.5)
spear = cor.test(dist[triup],conn.raw.vs.mot_VPTS[triup],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
pl = plot(bin, colramp=my_colors , legend=F,xlab='distance (mm)', ylab='VPT corr. of FC and FD (r)',main=rp.main.sp(sp.rho,sp.p,2)) 
hexVP.abline(pl$plot.vp,l,col = 'white', lwd = 5)
hexVP.abline(pl$plot.vp,l,col = gray(.6), lwd = 2)
dev.off()

## you dont need this.. 
#mean FD VS age 
# lme (random effect intercept)
#df = data.frame(x1 = age, x2 = male, y = fd.m, id = id)
#l = lme(y ~ x1 + x2, random = ~ 1|id, data = df)
#newdat = expand.grid(x2='male',x1=range(age))
#pred.m = predictSE.lme(l, newdata=expand.grid(x2='male',x1=range(age)), se.fit=T, level=0)
#pred.f = predictSE.lme(l, newdata=expand.grid(x2='female',x1=range(age)), se.fit=T, level=0)
#pred = list(); pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
#pdf(paste(plot.path,'age_fd_m_lme.pdf',sep=''),width=6,height=5)
#ggplot(df, aes(x=x1, y=y)) +
 # labs(x = 'age (years)', y = expression(paste(mu,' FD (mm)',sep='')), title = bquote(r^2*' = '*.(toString(signif(r.squaredGLMM(l)[1,1],2)))*', '*p[age]*' = '*.(toString(signif(summary(l)$tTable[2,5],2))))) +
 # geom_point(size=1, colour = 'grey') +
 # geom_path(aes(group = id), colour = 'grey') + # spaghetti plot
 # geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
 # geom_line(data=newdat, aes(y=pred$fit),size=1) +
 # theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
 # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()


# # MMP region ID's
# nroi.subc = 16  # number of subcortical regions
# nroi.cort = 360 # number of cortical regions
# mmp.subc = vector(length=nroi); mmp.subc[1:nroi.subc] = 1; mmp.subc[(nroi.subc+1):nroi] = 0
# mmp.cort = vector(length=nroi); mmp.cort[1:nroi.subc] = 0; mmp.cort[(nroi.subc+1):nroi] = 1

# # region names
#Glasser suppl material: https://static-content.springer.com/esm/art%3A10.1038%2Fnature18933/MediaObjects/41586_2016_BFnature18933_MOESM330_ESM.pdf
#make sure you exclude excluded regions 
nm = read.csv('~/R/fMRI/HCP.fsaverage_names.txt',header=F,as.is = array(1,dim=nroi))
 
nm.simple = vector(length = nroi)
nm.simple[1:8] = c('L thalamus','L caudate','L putamen','L pallidum','L hippocampus','L amygdala','L accumbens','L diencephalon')
nm.simple[9:16] = c('R thalamus','R caudate','R putamen','R pallidum','R hippocampus','R amygdala','R accumbens','R diencephalon')
for (n in 17:196) nm.simple[n] = paste('L ',strsplit(nm[n,1],'_')[[1]][2]) # left
for (n in 197:376) nm.simple[n] = paste('R ',strsplit(nm[n,1],'_')[[1]][2]) # right
 
# nm.all = nm.simple
# nm = nm.simple[hcp.keep.id]

# # save variables
# save('age','male','id','fd.m','fc','pos','nm','nm.all','dist',file='~/Desktop/fc_mmp.RData')

# # motion as function of age and sex, including interaction
# 
# # lme
# df = data.frame(x1 = age, x2 = male, y = fd.m, id = id)
# l = lme(y ~ x1 + x2 + x1:x2, random = ~ 1|id, data = df)
# 
# # gamm
# df = data.frame(x1 = age, x2 = male, y = fd.m, id = id); x1.pred = age.pred
# g = gamm(y ~ x2 + s(x1, by = x2, bs = 'bq', k = 10), random=list(id=~1), method = 'REML', data = df) # no main age effect, age(male) + age(female)
# 
# # gamm with sex as an ordered factor
# df = data.frame(x1 = age, x2 = ordered(male), y = fd.m, id = id); x1.pred = age.pred
# g = gamm(y ~ x2 + s(x1, by = x2, bs = 'bq', k = 10), random=list(id=~1), method = 'REML', data = df) # no main age effect, age(male) + age(female)
