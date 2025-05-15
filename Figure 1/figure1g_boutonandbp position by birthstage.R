#PN bouton location patterns
#originally written 11/28/2023
#loading packages------
library(dplyr)
library(ggplot2)
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
pn.info <- read.csv('PN_info.csv')
bout.info <- read.csv('bouton_info_collateral_bps.csv')
#add lineage information to bout.info
bout.info <- left_join(bout.info, pn.info[,c(1,4,6)])
#Jefferis 2007 Fig 1- boutons by NB origin and birth stage-----
#axons seperate in the iACT by birth time and lineage- embryonic vs larval and by parent NB
#add a variable that is birth time and lineage combined
bout.info$pn.birth.lineage <- paste(bout.info$pn.lineage, bout.info$birthstage, sep = ' ')
#make this variable a factor with an order
bout.info$pn.birth.lineage <- factor(bout.info$pn.birth.lineage, 
                                     levels = c('adPN embryonic','adPN early larval','adPN late larval','adPN NA',
                                                'lPN early larval', 'ilPN NA', 'lvPN NA'))
#color key
ggplot(bout.info, aes(x= center.x, y=center.z, color = pn.birth.lineage))+geom_point(shape = 16)+
  theme_classic()+xlab('medial to lateral')+ylab('ventral to dorsal')+
  xlim(-21600,-13200)+ylim(-17100,-12300)
ggsave('olfactorybouton_bybirthlineage_ml_dv_colorkey.pdf',
       width = 6, height = 6,units = 'in')
#ml dv 
ggplot(bout.info, aes(x= center.x, y=center.z, color = pn.birth.lineage))+geom_point(shape = 16)+
  theme_classic()+theme(legend.position = '')+xlab('medial to lateral')+ylab('ventral to dorsal')+
  xlim(-21600,-13200)+ylim(-17100,-12300)
ggsave('olfactorybouton_bybirthlineage_ml_dv.pdf',
       width = 6, height = 6,units = 'in')
#ml ap
ggplot(bout.info, aes(x= center.x, y=center.y, color = pn.birth.lineage))+geom_point(shape = 16)+
  theme_classic()+theme(legend.position = '')+xlab('medial to lateral')+ylab('posterior to anterior')+
  xlim(-21600,-13200)+ylim(-5400,400)
ggsave('olfactorybouton_bybirthlineage_ml_ap.pdf',
       width = 6, height = 6,units = 'in')
#ap dv
ggplot(bout.info, aes(x= center.y, y=center.z, color = pn.birth.lineage))+geom_point(shape = 16)+
  theme_classic()+theme(legend.position = '')+xlab('posterior to anterior')+ylab('ventral to dorsal')+
  xlim(-5400,400)+ylim(-17100,-12300)
ggsave('olfactorybouton_bybirthlineage_ap_dv.pdf',
       width = 6, height = 6,units = 'in')

#histogram of boutons x coordinate (medial to lateral)
#set y min and max based on the tallest histogram so that all histograms for all axes could be on the same axes
ymin <- 0
ymax <- 30
#set bin minimum and maximum values based on minimum and maximum values for bouton x coordinates
min(bout.info$center.x)
bmin <- -21600
max(bout.info$center.x)
bmax <- -13200
#plots histogram of bouton x coordinates separating by birthstage and lineage
birthlin <- levels(bout.info$pn.birth.lineage)
for(i in birthlin){
  ggplot(filter(bout.info, pn.birth.lineage == i), aes(x= center.x))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(bout.info, pn.birth.lineage == i)$center.x))+
    theme_classic()+ylim(ymin,ymax)+xlab(paste(i,'x location'))
  
  ggsave(paste(i,'_xlocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}
#histogram of bouton centers y coordinates (anterior to posterior)
#set bin minimum and maximum values based on minimum and maximum values for claw center y coordinates
min(bout.info$center.y)
bmin <- -5400
max(bout.info$center.y)
bmax <- 400

#plots histogram of bouton y coordinates separating by birthstage and lineage
for(i in birthlin){
  ggplot(filter(bout.info, pn.birth.lineage == i), aes(x= center.y))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(bout.info, pn.birth.lineage == i)$center.y))+
    theme_classic()+ylim(ymin,ymax)+xlab(paste(i,'y location'))
  
  ggsave(paste(i,'_ylocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}
#histogram of bouton center z coordinates (dorsal to ventral)
#set bin minimum and maximum values based on minimum and maximum values for claw center z coordinates
min(bout.info$center.z)
bmin <- -17100
max(bout.info$center.z)
bmax <- -12300
#plots histogram of bouton z coordinates separating by birthstage and lineage
for(i in birthlin){
  ggplot(filter(bout.info, pn.birth.lineage == i), aes(x= center.z))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(bout.info, pn.birth.lineage == i)$center.z))+
    theme_classic()+ylim(ymin,ymax)+xlab(paste(i,'z location'))
  
  ggsave(paste(i,'_zlocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}


#make branch point dataframe-----
bp.info <- unique(bout.info[,c(2,3,9:12,16:18)])
#plot locations of PN branch points ap dv
ggplot(bp.info, aes(x= bp.y, y=bp.z, color = pn.birth.lineage))+geom_point(shape = 16)+
  theme_classic()+theme(legend.position = '')+xlab('posterior to anterior')+ylab('ventral to dorsal')+
  xlim(-5400,400)+ylim(-18000,-12300)
ggsave('branchpoint_bybirthlineage_ap_dv.pdf',
       width = 6, height = 6,units = 'in')
#plots histogram of bouton y coordinates separating by birthstage and lineage
for(i in birthlin){
  ggplot(filter(bp.info, pn.birth.lineage == i), aes(x= bp.y))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(bp.info, pn.birth.lineage == i)$bp.y))+
    theme_classic()+ylim(ymin,ymax)+xlab(paste(i,'y location'))
  
  ggsave(paste(i,'_ylocation_bp.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}
