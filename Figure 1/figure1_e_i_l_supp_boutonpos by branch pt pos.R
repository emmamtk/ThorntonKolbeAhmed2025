#bouton branch point analyses
#originally written 1/17/2024
#load packages------
library(dplyr)
theme_set(theme_classic())
#set working directory and load files------
setwd("C:/Users/emmamtk/Dropbox (University of Michigan)/AhmedThorntonKolbe2022/DataFiles/CSVs for paper")
pn.info <- read.csv('PN_info.csv')
bout.info <- read.csv('bouton_info_collateral_bps.csv')

#calculate x,y,z differences betwen each bouton and its branch point (collateral length proxy)--------
bout.info$dx <- abs(bout.info$center.x- bout.info$bp.x)
bout.info$dy <- abs(bout.info$center.y- bout.info$bp.y)
bout.info$dz <- abs(bout.info$center.z- bout.info$bp.z)
bout.info$dbp <- sqrt((bout.info$dx^2)+(bout.info$dy^2)+(bout.info$dz^2))

#shapiro test for normality
shapiro.test(bout.info$center.x)
shapiro.test(bout.info$center.y)
shapiro.test(bout.info$center.z)
shapiro.test(bout.info$dbp)
shapiro.test(bout.info$coll.length)
#non are normal will use spearman correlation formula to calculate values

#bouton position on various axes by linear distance from bp------
test <- cor.test(bout.info$dbp,bout.info$center.x, method = 'spearman')
ggplot(bout.info,aes(dbp,center.x))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_distancefrombp_ML.pdf',
       width = 6, height = 6,units = 'in')
test <- cor.test(bout.info$dbp,bout.info$center.y, method = 'spearman')
ggplot(bout.info,aes(dbp,center.y))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_distancefrombp_AP.pdf',
       width = 6, height = 6,units = 'in')
test <- cor.test(bout.info$dbp,bout.info$center.z, method = 'spearman')
ggplot(bout.info,aes(dbp,center.z))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_distancefrombp_DV.pdf',
       width = 6, height = 6,units = 'in')
#bouton position on various axes by collateral length------
test <- cor.test(bout.info$coll.length,bout.info$center.x, method = 'spearman')
ggplot(bout.info,aes(coll.length,center.x))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_collateral_length_ML.pdf',
       width = 6, height = 6,units = 'in')
test <- cor.test(bout.info$coll.length,bout.info$center.y, method = 'spearman')
ggplot(bout.info,aes(coll.length,center.y))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_collateral_length_AP.pdf',
       width = 6, height = 6,units = 'in')
test <- cor.test(bout.info$coll.length,bout.info$center.z, method = 'spearman')
ggplot(bout.info,aes(coll.length,center.z))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_collateral_length_DV.pdf',
       width = 6, height = 6,units = 'in')
#bouton position on various axes by collateral length remove outlier------
test <- cor.test(bout.info$coll.length[bout.info$coll.length!=max(bout.info$coll.length)],
                 bout.info$center.x[bout.info$coll.length!=max(bout.info$coll.length)], method = 'spearman')
ggplot(bout.info[bout.info$coll.length!=max(bout.info$coll.length),],aes(coll.length,center.x))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_collateral_length_ML_nooutlier.pdf',
       width = 6, height = 6,units = 'in')
test <- cor.test(bout.info$coll.length[bout.info$coll.length!=max(bout.info$coll.length)],
                 bout.info$center.y[bout.info$coll.length!=max(bout.info$coll.length)], method = 'spearman')
ggplot(bout.info[bout.info$coll.length!=max(bout.info$coll.length),],aes(coll.length,center.y))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_collateral_length_AP_nooutlier.pdf',
       width = 6, height = 6,units = 'in')
test <- cor.test(bout.info$coll.length[bout.info$coll.length!=max(bout.info$coll.length)],
                 bout.info$center.z[bout.info$coll.length!=max(bout.info$coll.length)], method = 'spearman')
ggplot(bout.info[bout.info$coll.length!=max(bout.info$coll.length),],aes(coll.length,center.z))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_collateral_length_DV_nooutlier.pdf',
       width = 6, height = 6,units = 'in')

#bouton position and branchpoint position------
#based on ML position of bp
test <- cor.test(bout.info$bp.x,bout.info$center.x, method = 'spearman')
ggplot(bout.info,aes(bp.x,center.x))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_bpx_ml.pdf',
       width = 6, height = 6,units = 'in')

test <- cor.test(bout.info$bp.x,bout.info$center.y, method = 'spearman')
ggplot(bout.info,aes(bp.x,center.y))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_bpx_AP.pdf',
       width = 6, height = 6,units = 'in')

test <- cor.test(bout.info$bp.x,bout.info$center.z, method = 'spearman')
ggplot(bout.info,aes(bp.x,center.z))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_bpx_DV.pdf',
       width = 6, height = 6,units = 'in')

#based on AP position of bp
test <- cor.test(bout.info$bp.y,bout.info$center.x, method = 'spearman')
ggplot(bout.info,aes(bp.y,center.x))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_bpy_ml.pdf',
       width = 6, height = 6,units = 'in')
test <- cor.test(bout.info$bp.y,bout.info$center.y, method = 'spearman')
ggplot(bout.info,aes(bp.y,center.y))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_bpy_AP.pdf',
       width = 6, height = 6,units = 'in')
test <- cor.test(bout.info$bp.y,bout.info$center.z, method = 'spearman')
ggplot(bout.info,aes(bp.y,center.z))+
  stat_smooth(method = 'lm', col = 'red', se = T, level =0.95)+
  geom_point(shape =16)+
  labs(title = paste('S= ',test$statistic, 'p=', test$p.value, 'rho=', round(test$estimate,5 )))
ggsave('olfactorybouton_bpy_DV.pdf',
       width = 6, height = 6,units = 'in')

