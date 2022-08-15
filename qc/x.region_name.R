### ============== region names ============== ###
## ??? ###
rm(list = ls())
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


nroi = 376


#Glasser suppl material: https://static-content.springer.com/esm/art%3A10.1038%2Fnature18933/MediaObjects/41586_2016_BFnature18933_MOESM330_ESM.pdf
#make sure you exclude excluded regions 
# nm = read.csv('~/R/fMRI/HCP.fsaverage_names.txt',header=F,as.is = array(1,dim=nroi))
nm <- read_excel("../data/UCCHILD 374 region labels and numbers.xlsx")


# nm.simple = vector(length = nroi)
# nm.simple[1:8] = c('L thalamus','L caudate','L putamen','L pallidum','L hippocampus','L amygdala','L accumbens','L diencephalon')
# nm.simple[9:16] = c('R thalamus','R caudate','R putamen','R pallidum','R hippocampus','R amygdala','R accumbens','R diencephalon')
# for (n in 17:196) nm.simple[n] = paste('L ',strsplit(nm[n,1],'_')[[1]][2]) # left
# for (n in 197:376) nm.simple[n] = paste('R ',strsplit(nm[n,1],'_')[[1]][2]) # right
 