#bouton branch point analyses
#originally written 1/17/2024
#load packages------
library(dplyr)
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
pn.info <- read.csv('PN_info.csv')
bout.info <- read.csv('bouton_info_KCconnectivity.csv')[,1:12]

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
#non are normal will use spearman correlation formula to calculate values

#bouton position on various axes by collateral length------
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

