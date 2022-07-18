#### ======= loop to create fc ========= ####
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


# load ts.final, df_final
load("ts.final.RData")
load("df_final.RData")

# create fc
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

## save object
save(fc, file = "fc.RData")
## save as individual matrices for python CPM



#### subjs 76: "UCCHILDHB050" and 106:"UCCHILDHB090" have almost competely missing brains in rsFMRI ####
# nasubjs.in.fc.conn.raw = list(NA)

# for (i in 1:nroi.f) {
#   print(i)
#   for (j in 1:nroi.f) {
#     nasubjs.in.fc.conn.raw[[i]]=which(is.na(fc[i,j,]))
#   }
# }
# not sure what it does


#### visualise FC distributions (raw / fisher Z) ####
#can detect outliers here

#set up density plot parameters
triup = upper.tri(matrix(nrow=nroi.f,ncol=nroi.f))

dens.n = 500   # number of points in distribution
dens.i = -1.5; dens.f = 2.5      # initial and final values for distributions - determine by looking at the range on fc
abl = c(-1.5,0.5,1.5,2.5);       # horizontal lines on plots - also depends on your range
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

df_final$sex <- ifelse(df_final$sex == "Male", 1, 0)

# age_final = df_final$age
id_final=df_final$AP_id
age_final = rep(8, times=nrow(df_final)) ## need another spreadsheet?
fd.m_final=df_final$fd.m
group_final=df_final$group
male_final = df_final$male


# plots
#creates a page for each subject
pdf(paste(plot.path,'ind_fc_dist.pdf',sep=''),width=6,height=5)
for (n in 1:ns.f) {
  if (n%%10 == 0) print(n)
  par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 1.65, cex.axis = 1.5, cex.main = 1.5, font.main = 1, bg='white')
  plot(0, type='n', ylab='density', xlab='correlation', main=bquote(ID~.(id_final[n])*'  |  age = '*.(signif(age_final[n],3))*' y  |  '*mu*'FD = '*.(round(fd.m_final[n],2))*' mm'), ylim=c(0,3), xlim=c(dens.i,dens.f)) 
  abline(v=abl,col='grey',lty=2)
  lines(dens.x,dens.y.raw[n,],col='steelblue',lwd=2)
  lines(dens.x,dens.y[n,],col='darkorange',lwd=2)
  legend('topright',c('raw','fisher Z'),col=c('steelblue','darkorange'),lty=1,lwd=3,bty='n')
}
dev.off()


# coordinates of ROI centroids (dimension nroi x 3 (= x, y, z))
pos <- read.table("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", quote="\"", comment.char="")


# matrix of euclidean distance between regions
dist = matrix(nrow=nroi.f,ncol=nroi.f)
for (i in 1:nroi.f) {
    for (j in 1:nroi.f) {
        dist[i,j] = sqrt(sum((pos[i,]-pos[j,])^2))
    } 
}

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
  geom_violin() + geom_boxplot(width=0.1) #add boxplot
violin_fd.m
dev.off()

pdf(paste(plot.path,'violinplot_fdmax_bygroup.pdf',sep=''),width=6,height=5)
violin_fd.max <- ggplot(df_final, aes(x=group, y=fd.max), trim=FALSE) + 
  geom_violin() + geom_boxplot(width=0.1) #add boxplot
violin_fd.max
dev.off()


#### assumptions ####
library(lawstat)
## check if our data sets fulfill the homogeneity of variance assumption before we perform the t-test or Analysis of Variance (ANOVA).
levene.test(df_final$fd.m, df_final$group, location = c("mean")) # p-value = 0.9407
levene.test(df_final$fd.max, df_final$group, location = c("mean")) # p-value = 0.1452
# PT and FT have equal variance?


## check if a continuous variable follows a normal distribution.
shapiro.test(df_final$fd.m) # p-value = 0.05158, may be normally distributed
shapiro.test(df_final$fd.max) # p-value = 1.227e-05, not normally distributed


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

