#### ========== mean FC as a function of mean FD ============== ####
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

#### path to save plots ####
plot.path = "./plots/"

# load ts.final, df_final, fc.m_final
load("ts.final.RData")
load("df_final.RData")
load("fc.RData")
load("fc.m_final.RData")

df_final$sex <- ifelse(df_final$sex == "Male", 1, 0)
df_final$group <- ifelse(df_final$group == "PT", 1, 0)

nroi.f = dim(ts.final)[2]  #new number of rois -previous nroi (initial) = nroi
ns.f = dim(ts.final)[3] 



#relationship between FD and mean FC
#here we want to plot the inxn by group
#### plot interaction
#https://philippmasur.de/2018/11/26/visualizing-interaction-effects/
library(sjPlot)
library(sjmisc)
library(ggplot2)
theme_set(theme_sjplot())
# Setting up the building blocks
## in scrpit 4
# id_final=df_final$AP_id
# age_final = rep(8, times=nrow(df_final)) ## need another spreadsheet?
# fd.m_final=df_final$fd.m
# group_final=df_final$group
# male_final = df_final$male

pdf(paste(plot.path,'fdm_fcm_inxn_bygroup.pdf',sep=''),width=6,height=5)
basic_plot <- ggplot(df_final,
                     aes(x = fd.m,
                         y = fc.m_final,
                         color = group)) +
  theme_bw() +
  labs(x = "fd.m",
       y = "fc.m",
       color = "group") +
  geom_point() +  # Colored scatterplot
  geom_point(alpha = .3, 
             size = .9) + 
  geom_smooth(method = "lm") # Colored scatterplot and regression lines
basic_plot
dev.off()


####### not sure what's this for ################
#l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
# l = lm(y ~ x1 + x2 + x1*x2, data=df)
# newdat = expand.grid(x2='male',x1=range(fd.m))
# #pred.m = predictSE.lme(l, newdata=expand.grid(x2='male',x1=range(fd.m)), se.fit=T, level=0) #predicting the error bars - CIs for the lme 
# #pred.f = predictSE.lme(l, newdata=expand.grid(x2='female',x1=range(fd.m)), se.fit=T, level=0)
# #pred = list(); pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
# pdf(paste(plot.path,'fd_m_fc_raw_m_lme.pdf',sep=''),width=6,height=5)
# ggplot(df, aes(x=x1, y=y)) +
#   labs(x = expression(paste(mu, ' FD (mm)',sep='')), y = expression(paste(mu, ' FC',sep='')), title = bquote(r^2*' = '*.(toString(signif(r.squaredGLMM(l)[1,1],2)))*', '*p[FD]*' = '*.(toString(signif(summary(l)$tTable[2,5],2))))) +
#   geom_point(size=1, colour = 'grey') +
#   #geom_path(aes(group = id), colour = 'grey') + # spaghetti plot
#   geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
#   geom_line(data=newdat, aes(y=pred$fit),size=1) +
#   theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()
#############################################

l = lm(fc.m_final~df_final$fd.m + df_final$group + df_final$fd.m*df_final$group) ## ?



#plot the interaction by group
# look at whether movement is related to mean FC - for each group seperately then see if there's an interaction
df_final_with_fcmean <- cbind(df_final, fc.m_final)

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


################################################
### another way of plotting this:
#plot(fd.m,fc.m,xlab='...',ylab='...')
fd.pred = seq(from=min(df_final$fd.m),to=max(df_final$fd.m), length.out=100)
pred = predict(l,newdata=data.frame('x1'=fd.pred), interval='confidence') ## ?
#polygon(c(rev(fd.pred), fd.pred), c(rev(pred[ ,3]), pred[ ,2]), col = col2alpha('grey',alpha=0.5), border = NA)
#lines(fd.pred,pred[,1],col='black',lwd=3);

#inxn
summary(lm(fc.m_final~df_final$fd.m + df_final$group + df_final$fd.m*df_final$group))
################################################



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
    conn.raw.vs.mot[i,j] = cor(fc[i,j,], df_final$fd.m, method='pearson')
  }
}

##hexbin
# set-up
# coordinates of ROI centroids (dimension nroi x 3 (= x, y, z))
pos <- read.table("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", quote="\"", comment.char="")

triup = upper.tri(matrix(nrow=nroi.f,ncol=nroi.f))

# matrix of euclidean distance between regions
dist = matrix(nrow=nroi.f,ncol=nroi.f)
for (i in 1:nroi.f) {
    for (j in 1:nroi.f) {
        dist[i,j] = sqrt(sum((pos[i,]-pos[j,])^2))
    } 
}
save(dist, file = "dist.RData")

## bin and color
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
    conn.raw.vs.mot_CONTROLS[i,j] = cor(fc_controls[i,j,], df_final$fd.m[controls_in_finaldf],method='pearson')
  }
}
#for VPTs:
conn.raw.vs.mot_VPTS = array(NA,dim=c(nroi.f,nroi.f))
for (i in 1:nroi.f) {
  print(i)
  for (j in 1:nroi.f) {
    conn.raw.vs.mot_VPTS[i,j] = cor(fc_VPTs[i,j,], df_final$fd.m[VPTs_in_finaldf],method='pearson')
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
