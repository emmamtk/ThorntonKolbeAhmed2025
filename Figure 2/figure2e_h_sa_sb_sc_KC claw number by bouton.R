#KC claw number per bouton analysis.
#originally written 3/8/2024
#load packages------
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(patchwork)
theme_set(theme_classic())
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
#KC info and claw info only for olfactory KCs (those in the main calyx)
kc.info <- read.csv('KC_info.csv')
claw.info <- read.csv('claw_info.csv')
#pn and bouton info
pn.info <- read.csv('PN_info.csv')
bout.info <- read.csv('bouton_info_KCconnectivity.csv')[,1:8]
#calculate number of claws on each bouton
for(i in 1:nrow(bout.info)){
  bclaws <- filter(claw.info, pn.boutId == bout.info$pn.boutId[i])
  #number of claws on that bouton
  bout.info$nclaw[i] <- nrow(bclaws)
}
#filter out boutons with 0 claws
bout.info <- bout.info[bout.info$nclaw>0,]


#plot number of claws per bouton distribution (2e)-----------
ggplot(bout.info, aes(x=nclaw))+geom_histogram(binwidth = 1)+geom_vline(xintercept = quantile(bout.info$nclaw))+
  geom_vline(xintercept = bout.info$nclaw[bout.info$pn.glomerulus=='DC1'], color='green')+
  geom_vline(xintercept = bout.info$nclaw[bout.info$pn.glomerulus=='DM3'], color='cadetblue')


ggsave('clawperbouton_histogram_DC1_values.pdf',
       width = 6, height = 4,units = 'in')
#plot number of claws by size of bouton (2e)----------
#shapiro test for normality
shapiro.test(bout.info$nclaw)
shapiro.test(bout.info$nsyn)
#none are normal will use spearman correlation formula to calculate values
test <- cor.test(bout.info$nclaw,bout.info$nsyn, method = 'spearman')
ggplot(bout.info,aes(nsyn,nclaw))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  geom_point(data = filter(bout.info, pn.glomerulus == 'DC1'), aes(x=nsyn, y = nclaw), color = 'green')+
  geom_point(data = filter(bout.info, pn.glomerulus == 'DM3'), aes(x=nsyn, y = nclaw), color = 'cadetblue')+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('clawperbouton_byboutonsize_withDC1DM3.pdf',
       width = 6, height = 6,units = 'in')


#bouton size along different anatomical axes (S2a)--------
#shapiro test for normality
shapiro.test(bout.info$center.x)
shapiro.test(bout.info$center.y)
shapiro.test(bout.info$center.z)
shapiro.test(bout.info$nsyn)
#none are normal will use spearman correlation formula to calculate values
test <- cor.test(bout.info$center.x,bout.info$nsyn, method = 'spearman')
ggplot(bout.info,aes(center.x,nsyn))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('boutonsizebybouton_xlocation.pdf',
       width = 6, height = 6,units = 'in')

test <- cor.test(bout.info$center.y,bout.info$nsyn, method = 'spearman')
ggplot(bout.info,aes(center.y,nsyn))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('boutonsizebybouton_ylocation.pdf',
       width = 6, height = 6,units = 'in')

test <- cor.test(bout.info$center.z,bout.info$nsyn, method = 'spearman')
ggplot(bout.info,aes(center.z,nsyn))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('boutonsizebybouton_zlocation.pdf',
       width = 6, height = 6,units = 'in')


#number of claws per bouton by bouton position (S2b)-------------
#plot number of claws by bouton ML position
mlmid <- median(bout.info$center.x)

shapiro.test(abs(mlmid-bout.info$center.x))
shapiro.test(bout.info$nclaw)

test <- cor.test(abs(mlmid-bout.info$center.x),bout.info$nclaw, method = 'spearman')
ggplot(bout.info,aes(abs(mlmid-bout.info$center.x),nclaw))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('clawperbouton_byboutonX.pdf',
       width = 6, height = 6,units = 'in')

#plot number of claws by bouton AP position
test <- cor.test(bout.info$center.y,bout.info$nclaw, method = 'spearman')
ggplot(bout.info,aes(center.y,nclaw))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('clawperbouton_byboutonY.pdf', width = 6, height = 6, units = 'in')
#plot number of claws by bouton DV position
test <- cor.test(bout.info$center.z,bout.info$nclaw, method = 'spearman')
ggplot(bout.info,aes(center.z,nclaw))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('clawperbouton_byboutonZ.pdf', width = 6, height = 6, units = 'in')
#plot boutons colored by number of claws (S2c)--------------------
#adding a column to assign bins to each claw per bouton value
bout.info <- bout.info%>% mutate(new_bin = cut(nclaw, breaks = c(seq(0,50,5),85)))
#bins are 5 apart between 0 and 50 there is another bin for values between 50 and 85
#making a palette for the binned data wtih 10 values 11th value will be black for boutons with more than 50 claws
pal <- colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7"))(10)
#plot ML DV
ggplot(bout.info, aes(x= center.x, y = center.z, color = new_bin))+geom_point(shape = 20)+
  scale_color_manual(values=c(pal,'black'))
ggsave('clawperbouton_ML_DV.pdf', width = 6, height = 4, units = 'in')
#plot AP DV
ggplot(bout.info, aes(x= center.y, y = center.z, color = new_bin))+geom_point(shape = 20)+
  scale_color_manual(values=c(pal,'black'))
ggsave('clawperbouton_AP_DV.pdf', width = 6, height = 4, units = 'in')

#########fig 2h CLAW SPREAD###############
#for each kc calculate maximum spread on each axis---------
for(i in 1:nrow(kc.info)){
  kcclaws <- filter(claw.info, kc.bodyId == kc.info$kc.bodyId[i])
  kc.info$x.spread[i] <- abs(max(kcclaws$center.x)-min(kcclaws$center.x))
  kc.info$y.spread[i] <- abs(max(kcclaws$center.y)-min(kcclaws$center.y))
  kc.info$z.spread[i] <- abs(max(kcclaws$center.z)-min(kcclaws$center.z))
}
#plot ML (x) spread -----------
#size of M-L axis
mlscale <- max(bout.info$center.x)-min(bout.info$center.x)
ggplot()+
  #ML quarters
  geom_vline(xintercept=c(mlscale*.25,mlscale*.5,mlscale*.75,mlscale),color = 'red')+
  #all KCs
  stat_bin(data = filter(kc.info, nclaw>2), aes(x=x.spread),geom = 'line',binwidth = 100,color = 'black')+
  geom_vline(xintercept = median(kc.info$x.spread), color = 'black')+
  #gamma KCs
  stat_bin(data = filter(kc.info, kc.largetype=='gamma',nclaw>2), aes(x=x.spread),geom = 'line',binwidth = 100,color = 'darkblue')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'gamma')$x.spread), color = 'darkblue')+
  #alpha'beta' KCs
  stat_bin(data = filter(kc.info, kc.largetype=='alphaprimebetaprime',nclaw>2), aes(x=x.spread),geom = 'line',binwidth = 100,color = 'green')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'alphaprimebetaprime')$x.spread), color = 'green')+
  #alpha beta KCs
  stat_bin(data = filter(kc.info, kc.largetype=='alphabeta',nclaw>2), aes(x=x.spread),geom = 'line',binwidth = 100,color = 'cadetblue')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'alphabeta')$x.spread), color = 'cadetblue')+
  ggtitle('blue gamma, green abprime lblue ab')+
  xlim(0,8000)+ylim(0,110)
ggsave('ML_x_clawspread_bytype.pdf',
       width = 4, height = 3,units = 'in')
#plot AP (y) spread ----------------
apscale <- max(bout.info$center.y)-min(bout.info$center.y)
ggplot()+
  #AP quarters
  geom_vline(xintercept=c(apscale*.25,apscale*.5,apscale*.75,apscale),color = 'red')+
  #all KCs
  stat_bin(data = filter(kc.info, nclaw>2), aes(x=y.spread),geom = 'line',binwidth = 100,color = 'black')+
  geom_vline(xintercept = median(kc.info$y.spread), color = 'black')+
  #gamma KCs
  stat_bin(data = filter(kc.info, kc.largetype=='gamma',nclaw>2), aes(x=y.spread),geom = 'line',binwidth = 100,color = 'darkblue')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'gamma')$y.spread), color = 'darkblue')+
  #alpha'beta' KCs
  stat_bin(data = filter(kc.info, kc.largetype=='alphaprimebetaprime',nclaw>2), aes(x=y.spread),geom = 'line',binwidth = 100,color = 'green')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'alphaprimebetaprime')$y.spread), color = 'green')+
  #alpha beta KCs
  stat_bin(data = filter(kc.info, kc.largetype=='alphabeta',nclaw>2), aes(x=y.spread),geom = 'line',binwidth = 100,color = 'cadetblue')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'alphabeta')$y.spread), color = 'cadetblue')+
  ggtitle('blue gamma, green abprime lblue ab')+
  xlim(0,8000)+ylim(0,110)
ggsave('AP_y_clawspread_bytype.pdf',
       width = 4, height = 3,units = 'in')
#plot DV spread (z)--------------
dvscale <- max(bout.info$center.z)-min(bout.info$center.z)
ggplot()+
  #ML quarters
  geom_vline(xintercept=c(dvscale*.25,dvscale*.5,dvscale*.75,dvscale),color = 'red')+
  #all KCs
  stat_bin(data = filter(kc.info, nclaw>2), aes(x=z.spread),geom = 'line',binwidth = 100,color = 'black')+
  geom_vline(xintercept = median(kc.info$z.spread), color = 'black')+
  #gamma KCs
  stat_bin(data = filter(kc.info, kc.largetype=='gamma',nclaw>2), aes(x=z.spread),geom = 'line',binwidth = 100,color = 'darkblue')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'gamma')$z.spread), color = 'darkblue')+
  #alpha'beta' KCs
  stat_bin(data = filter(kc.info, kc.largetype=='alphaprimebetaprime',nclaw>2), aes(x=z.spread),geom = 'line',binwidth = 100,color = 'green')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'alphaprimebetaprime')$z.spread), color = 'green')+
  #alpha beta KCs
  stat_bin(data = filter(kc.info, kc.largetype=='alphabeta',nclaw>2), aes(x=z.spread),geom = 'line',binwidth = 100,color = 'cadetblue')+
  geom_vline(xintercept = median(filter(kc.info, kc.largetype == 'alphabeta')$z.spread), color = 'cadetblue')+
  ggtitle('blue gamma, green abprime lblue ab')+
  xlim(0,8000)+ylim(0,110)
ggsave('DV_z_clawspread_bytype.pdf',
       width = 4, height = 3,units = 'in')


#plot claw number vs claw spread----------------------
ggplot(filter(kc.info, nclaw>2, kc.largetype != 'unknown'), 
       aes(x= as.factor(nclaw), y = x.spread))+
  geom_jitter()+geom_boxplot()+facet_grid(~kc.largetype)+
  ylim(0,6500)
ggsave('ML_x_clawspread_byclawnumber.pdf',
       width = 10, height = 4,units = 'in')

ggplot(filter(kc.info, nclaw>2, kc.largetype != 'unknown'), 
       aes(x= as.factor(nclaw), y = y.spread))+
  geom_jitter()+geom_boxplot()+facet_grid(~kc.largetype)+
  ylim(0,6500)
ggsave('AP_y_clawspread_byclawnumber.pdf',
       width = 10, height = 4,units = 'in')

ggplot(filter(kc.info, nclaw>2, kc.largetype != 'unknown'), 
       aes(x= as.factor(nclaw), y = z.spread))+
  geom_jitter()+geom_boxplot()+facet_grid(~kc.largetype)+
  ylim(0,6500)
ggsave('DV_z_clawspread_byclawnumber.pdf',
       width = 10, height = 4,units = 'in')

