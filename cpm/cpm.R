# install.packages("NetworkToolbox")

library(NetworkToolbox)


# bstat, Behavioral statistic for each participant with neural data (a vector)
behav <- behavOpen

# neuralarray
# https://drive.google.com/file/d/1ugwi7nRrlHQYuGPzEB4wYzsizFrIMvKR/view

load("restOpen.rda")
class(restOpen)
dim(restOpen) # 3D array

# Run cpmIV
res <- cpmIV(neuralarray = restOpen, bstat = behav,
                model = "linear", method = "mean", cores = 2)

save(res, file = "cpm_test.Rdata")
#save.image(file="testCPM.RData")

write.table(res$posMask, file = "pos_behavopen_restopen.csv", quote = F, sep = ",",
            row.names = F, col.names = F)
write.table(res$negMask, file = "neg_behavopen_restopen.csv", quote = F, sep = ",",
            row.names = F, col.names = F)


# Plot cpmIV results
load("cpm_test.Rdata")

png("cpmtest_behave-restopen.png", width = 700, height = 580)
cpmPlot(res, visual.nets = T)
dev.off()
## Error in cpmPlot(res) : could not find function "cpmPlot" # solved, newest version on github
