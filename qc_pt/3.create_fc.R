#### ======= loop to create fc ========= ####
rm(list = ls())
setwd("/mnt/d/PROJECTS/preterm_language/qc_pt")

library(dplyr)
# plots
# library(readxl)
# library(seqinr)       # col2alpha
# library(ggplot2)
# library(hexbin)
# library(RColorBrewer) # brewer.pal
# # stats / modelling
# library(matrixStats)  # colMaxs
library(psych)        # fisherz
# library(nlme)         # lme
# library(mgcv)         # gam(m)
# library(vows)         # bs = 'bq'
# library(AICcmodavg)   # ?
# library(MuMIn)        # r.squaredGLMM
# library(stringr)
#functions
source('rp.main.R')

#### path to save plots ####
plot.path <- "./plots/"


# load ts_final, df_final
load("ts_final.RData")
load("df_final.RData")


# create fc
nroi.f <- dim(ts_final)[2]  ## 374, regions #new number of rois -previous nroi (initial) = nroi
ns.f <- dim(ts_final)[3] ## 88, subjects

fc.raw = array(NA,dim=c(nroi.f,nroi.f,ns.f))  # raw FC - before fisher transformation
for (i in 1:ns.f) {
  if (i%%10==0) print(i)
  fc.raw[,,i] = cor(ts_final[,,i])
}
dim(fc.raw)
# [1] 374 374 88

# fisher transform correlations
fc <- fisherz(fc.raw)
dim(fc)
# [1] 374 374 88


## save as individual matrices for python CPM
# diagonals can't be NA
for (n in 1:ns.f) {
  fc[,,n][seq(from=1,to=(nroi.f^2),by=(nroi.f+1))] <- 3
}

# export
for (i in 1:ns.f) {
  write.table(fc[,,i], file = paste0("../data/fc_individual_pt/", row.names(df_final)[i], ".fc.txt"), sep = " ",
              quote = F, row.names = F, col.names = F)
}


## for in R
# set matrix diagonals to NaN
for (n in 1:ns.f) {
  fc[,,n][seq(from=1,to=(nroi.f^2),by=(nroi.f+1))] <- NaN
}

## save object
save(fc, file = "fc.RData")


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
# indices of upper triangular part of matrix
triup <- upper.tri(matrix(nrow=nroi.f,ncol=nroi.f))

# mean FC
fc.m_final <- vector(length=ns.f)
for (n in 1:ns.f) fc.m_final[n] <- mean(fc[,,n][triup],na.rm=T)
save(fc.m_final, file = "fc.m_final.RData")