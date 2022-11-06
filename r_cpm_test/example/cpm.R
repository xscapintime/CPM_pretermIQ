# install.packages("NetworkToolbox")
setwd("/mnt/d/PROJECTS/preterm_language/cpm")

library(NetworkToolbox)
library(tidyverse)


# bstat, Behavioral statistic for each participant with neural data (a vector)
behav <- behavOpen

# neuralarray
# https://drive.google.com/file/d/1ugwi7nRrlHQYuGPzEB4wYzsizFrIMvKR/view

load("restOpen.rda")
class(restOpen)
dim(restOpen) # 3D array

## Run cpmIV
res <- cpmIV(neuralarray = restOpen, bstat = behav,
                model = "linear", method = "mean", cores = 2)

save(res, file = "cpm_test.Rdata")
#save.image(file="testCPM.RData")

write.table(res$posMask, file = "pos_behavopen_restopen.csv", quote = F, sep = ",",
            row.names = F, col.names = F)
write.table(res$negMask, file = "neg_behavopen_restopen.csv", quote = F, sep = ",",
            row.names = F, col.names = F)

# permutation
system.time(
  perm_res <- cpmIVperm(neuralarray = restOpen, bstat = behav,
             model = "linear", method = "mean", cores = 2, iter = 20)
)

perm_res
# Positive Prediction Negative Prediction
# p-value                 0.6                 0.15

## cpmEV
# split dataset to training and testing
sample <- sample(c(TRUE, FALSE), dim(restOpen)[3], replace=TRUE, prob=c(0.7,0.3))
train  <- restOpen[,,sample]
test   <- restOpen[,,!sample]

dim(train)
dim(test)

train_b <- behav[sample]
test_b <- behav[!sample]


# cpmEV
system.time(
  res_ev <- cpmEV(train_na = train, train_b = train_b, valid_na = test, valid_b = test_b,
                thresh = .01, overlap = T)
)
## overlap=F doesn't work, overlap=T --> leave-one-out
# r in cpmEV(train_na = train, train_b = train_b, valid_na = test, valid_b = test_b,  : 
# object 'train_mats_tmp' not found





# Plot cpmIV results
load("cpm_test.Rdata")

png("cpmtest_behave-restopen.png", width = 700, height = 580)
cpmPlot(res, visual.nets = T)
dev.off()
## Error in cpmPlot(res) : could not find function "cpmPlot" # solved, newest version on github
