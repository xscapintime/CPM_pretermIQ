library(NetworkToolbox)

### == Load fake dataframe ====
load("fake_dataset.RData")
dim(fake_dataset)

### ==== FC - Calculate global connectivity ====
# averaged the upper tri per subject and got a mean FC per subject
load('../cpm/restOpen.rda')
fc <- restOpen
dim(fc) # 268 * 268 * 144


# use stats above as the bstat in cpm ## not sure
# cpm
# Run cpmIV
res <- cpmIV(neuralarray = fc, bstat = fake_dataset$fd.m,
                model = "linear", method = "mean", cores = 2)

#               r     p   mae
# positive -0.040 0.633 0.054
# negative -0.123 0.141 0.058
#           rmse
# positive 0.070
# negative 0.071
save(res, file = "cpm_fakefd.Rdata")

write.table(res$posMask, file = "pos_fake_restopen.csv", quote = F, sep = ",",
            row.names = F, col.names = F)
write.table(res$negMask, file = "neg_fake_restopen.csv", quote = F, sep = ",",
            row.names = F, col.names = F)


# Plot cpmIV results
load("cpm_fakefd.Rdata")

png("cpm_fakefd-restopen.png", width = 700, height = 580)
cpmPlot(res, visual.nets = T)
dev.off()

