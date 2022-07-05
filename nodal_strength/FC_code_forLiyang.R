### == Install Packages ====
# install.packages("readxl")
# install.packages("openxlsx")
# install.packages("purrr")
# install.packages("ggplot2")
# install.packages("gdata")
# install.packages("diceR")
# install.packages("fields")
# install.packages("RANN")

### == Load libraries ====
library(dplyr)
library(readxl)
library(seqinr)       # col2alpha
library(ggplot2)
library(hexbin)
library(RColorBrewer) # brewer.pal
library(matrixStats)  # colMaxs
library(psych)        # fisherz
library(nlme)         # lme
library(mgcv)         # gam(m)
library(vows)         # bs = 'bq'
library(AICcmodavg)   # ?
library(MuMIn)        # r.squaredGLMM
library(openxlsx)
library(data.table)
library(purrr)
library(gdata)
library(diceR)
library(alluvial)
library(haven)
library(naniar)
library(Hmisc)
library(onewaytests)
library(lawstat)
library(data.table)
library(coin)
library(rcompanion)
library(fields)
library(caret)
library(RANN)
library(scales)
library(Matrix)

### == Load functions ====
source('rp.main.R')
#source('plot.smooth.mat.R') # function to plot smooth matrices with cells of arbitrary size # based on image.plot
#source('plot.unsmooth.mat.R') # function to plot unsmooth matrices with cells of arbitrary size # based on image.plot
source('node.str.R') # function to calculate node strength, ignoring self-connections (on the diagonal)

### == Load fake dataframe ====
load("fake_dataset.RData")
dim(fake_dataset)

### ==== FC - Calculate global connectivity ====
# averaged the upper tri per subject and got a mean FC per subject
load('../cpm/restOpen.rda')
fc <- restOpen
dim(fc) # 268 * 268 * 144
ns.f <- dim(fc)[3] # participants
nroi.f <- dim(fc)[2] # connectivity matrice

mean_upper_tri <- vector(length=ns.f)
for (n in 1:ns.f) { ## for each participant
  # mean_upper_tri[n] <- mean(fc[,,n][triup]) # no triup
  mat <- fc[,,n]
  # mat[upper.tri(mat)] <- 0
  # mean_upper_tri[n] <- mean(mat)
  mean_upper_tri[n] <- mean(triu(mat)) # triu from Matrix works better
}

## check diff.
# fc[,,n][upper.tri(fc[,,n], diag=T)] [fc[,,n][upper.tri(fc[,,n], diag=T)] != 0] %>% mean()
# triu(fc[,,n])[triu(fc[,,n]) != 0] %>% mean()


# linear model checking for an effect of group, correcting for age, sex and motion
summary(lm(mean_upper_tri ~ fake_dataset$Group + fake_dataset$Age + fake_dataset$Gender + fake_dataset$fd.m))

### ==== FC - regional strength differences between the groups ====
# only keep positive connections: turn negative correlations to zeros
fc_pos <- fc
fc_pos[which(fc < 0)] <- 0
## Regional measures - calculate node strength
str = apply(fc_pos,3,node.str) # strength

# linear models and extraction of paremeters for each regional strength
str.p = str.t = vector(length = nroi.f)

for (n in 1:nroi.f) {
  if (n%%10 == 0) print(n)
  # str
  l = lm(str[n,] ~ fake_dataset$Group + fake_dataset$Age + fake_dataset$Gender + fake_dataset$fd.m);
  str.p[n] = summary(l)$coefficients[2,4] # Pr(>|t|) 
  str.t[n] = summary(l)$coefficients[2,3] # t value
  # ...
}

#which is significant before FDR
which(str.p <0.05)
# [1]  87 114 268

# FDR correct p-values
str.p.fdr = p.adjust(str.p, method = 'fdr')

#which is significant after FDR
which(str.p.fdr <0.05)
# none


hist(str.p)
hist(str.p.fdr)
hist(str.t)
