#Claw location by clone and type
#originally written 9/7/2023
#edited for clarity 11/14/2023
#load packages------
library(ggplot2)
library(dplyr)
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
#KC info and claw info only for olfactory KCs (those in the main calyx)
kc.info <- read.csv('KC_info.csv')
claw.info <- read.csv('claw_info.csv')
#counts of KCs of each clone and type
#gamma
nrow(filter(kc.info, nb.cluster=='a',kc.largetype=='gamma'))
nrow(filter(kc.info, nb.cluster=='b',kc.largetype=='gamma'))
nrow(filter(kc.info, nb.cluster=='c',kc.largetype=='gamma'))
nrow(filter(kc.info, nb.cluster=='d',kc.largetype=='gamma'))
#alphaprimebetaprime
nrow(filter(kc.info, nb.cluster=='a',kc.largetype=='alphaprimebetaprime'))
nrow(filter(kc.info, nb.cluster=='b',kc.largetype=='alphaprimebetaprime'))
nrow(filter(kc.info, nb.cluster=='c',kc.largetype=='alphaprimebetaprime'))
nrow(filter(kc.info, nb.cluster=='d',kc.largetype=='alphaprimebetaprime'))
#alphabeta
nrow(filter(kc.info, nb.cluster=='a',kc.largetype=='alphabeta'))
nrow(filter(kc.info, nb.cluster=='b',kc.largetype=='alphabeta'))
nrow(filter(kc.info, nb.cluster=='c',kc.largetype=='alphabeta'))
nrow(filter(kc.info, nb.cluster=='d',kc.largetype=='alphabeta'))




#plot location of claw centers in 2d colored by clonal origin----------
ggplot(claw.info[!is.na(claw.info$nb.cluster),], aes(x= center.x, y=center.z, color = nb.cluster))+geom_point(shape = 1)+
  theme_classic()+xlab('medial to lateral')+ylab('ventral to dorsal')+
  xlim(-21600,-13200)+ylim(-17100,-12300)
ggsave('claw location by clone key.pdf',
       width = 3, height = 3,units = 'in')

ggplot(claw.info[!is.na(claw.info$nb.cluster),], aes(x= center.x, y=center.z, color = nb.cluster))+geom_point(shape = 1)+
  theme_classic()+theme(legend.position = '')+xlab('medial to lateral')+ylab('ventral to dorsal')+
  xlim(-21600,-13200)+ylim(-17100,-12300)
ggsave('filtered_rotated_clawlocation_byclone_ml_dv.pdf',
       width = 6, height = 6,units = 'in')

ggplot(claw.info[!is.na(claw.info$nb.cluster),], aes(x= center.x, y=center.y, color = nb.cluster))+geom_point(shape = 1)+
  theme_classic()+theme(legend.position = '')+xlab('medial to lateral')+ylab('posterior to anterior')+
  xlim(-21600,-13200)+ylim(-5400,400)
ggsave('filtered_rotated_clawlocation_byclone_ml_ap.pdf',
       width = 6, height = 6,units = 'in')

ggplot(claw.info[!is.na(claw.info$nb.cluster),], aes(x= center.y, y=center.z, color = nb.cluster))+geom_point(shape = 1)+
  theme_classic()+theme(legend.position = '')+xlab('posterior to anterior')+ylab('ventral to dorsal')+
  xlim(-5400,400)+ylim(-17100,-12300)
ggsave('filtered_rotated_clawlocation_byclone_ap_dv.pdf',
       width = 6, height = 6,units = 'in')
#plot location of claw centers in 2d colored by type-------------
#medial lateral vs ventral dorsal with key
ggplot(filter(claw.info[!is.na(claw.info$kc.largetype),], kc.largetype != 'unknown'), 
       aes(x= center.x, y=center.z, color =kc.largetype))+
  geom_point(shape = 1)+
  theme_classic()+xlab('medial to lateral')+ylab('ventral to dorsal')+
  xlim(-21600,-13200)+ylim(-17100,-12300)
ggsave('claw location by type key.pdf',
       width = 3, height = 3,units = 'in')

#medial lateral vs ventral dorsal
ggplot(filter(claw.info[!is.na(claw.info$kc.largetype),], kc.largetype != 'unknown'), 
       aes(x= center.x, y=center.z, color =kc.largetype))+
  geom_point(shape = 1)+
  theme_classic()+theme(legend.position = '')+xlab('medial to lateral')+ylab('ventral to dorsal')+
  xlim(-21600,-13200)+ylim(-17100,-12300)
ggsave('filtered_rotated_clawlocation_bytype_ml_dv.pdf',
       width = 6, height = 6,units = 'in')

#medial lateral vs posterior to anterior
ggplot(filter(claw.info[!is.na(claw.info$kc.largetype),], kc.largetype != 'unknown'), 
       aes(x= center.x, y=center.y, color = kc.largetype))+
  geom_point(shape = 1)+
  theme_classic()+theme(legend.position = '')+xlab('medial to lateral')+ylab('posterior to anterior')+
  xlim(-21600,-13200)+ylim(-5400,400)
ggsave('filtered_rotated_clawlocation_bytype_ml_ap.pdf',
       width = 6, height = 6,units = 'in')

#posterior to anterior vs ventral dorsal
ggplot(filter(claw.info[!is.na(claw.info$kc.largetype),], kc.largetype != 'unknown'), 
       aes(x= center.y, y=center.z, color = kc.largetype))+
  geom_point(shape = 1)+
  theme_classic()+theme(legend.position = '')+xlab('posterior to anterior')+ylab('ventral to dorsal')+
  xlim(-5400,400)+ylim(-17100,-12300)
ggsave('filtered_rotated_clawlocation_bytype_ap_dv.pdf',
       width = 6, height = 6,units = 'in')


#histogram of claw centers x coordinate (medial to lateral)--------------
#set y min and max based on the tallest histogram so that all histograms for all axes could be on the same axes
ymin <- 0
ymax <- 410
#set bin minimum and maximum values based on minimum and maximum values for claw center x coordinates
min(claw.info$center.x)
bmin <- -21600
max(claw.info$center.x)
bmax <- -13200
#plots histogram of claw center x coordinates
ggplot(claw.info, aes(x= center.x))+
  stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
  geom_vline(xintercept = median(claw.info$center.x, col = 'red'))+
  theme_classic()+ylim(ymin,ymax)
ggsave('all claws x location.pdf',
       width = 6, height = 3,units = 'in')

#plots histogram of claw center x coordinates separating claws by clonal origin
clones <- c('a','b','c','d')
for(i in clones){
  ggplot(filter(claw.info, nb.cluster == i), aes(x= center.x))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(claw.info, nb.cluster == i)$center.x), col = 'red')+
    theme_classic()+ylim(ymin,ymax)+xlab(paste('clone',i,'x location'))
  
  ggsave(paste('clone',i,'_xlocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}

#plots histogram of claw center x coordinates separating claws by type 
types = c('gamma','abprime','ab')
for(i in types){
  ggplot(filter(claw.info,kc.largetype == i), aes(x= center.x))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(claw.info, kc.largetype == i)$center.x), col = 'red')+
    theme_classic()+ylim(ymin,ymax)+xlab(paste(i,'x location'))
  
  ggsave(paste(i,'_xlocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}

#plots histogram of claw center x coordinates separating claws by clone and type 
for(i in clones){
  for(j in types){
    ggplot(filter(claw.info,kc.largetype == j, nb.cluster == i), aes(x= center.x))+
      stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
      geom_vline(xintercept = median(filter(claw.info, kc.largetype == j, nb.cluster == i)$center.x), col = 'red')+
      theme_classic()+ylim(ymin,ymax)+xlab(paste('clone',i,j,'x location'))
    
    ggsave(paste('clone',i,'_',j,'_xlocation.pdf', sep = ''),
           width = 6, height = 3,units = 'in')
    
  }
}


#histogram of claw centers y coordinate (anterior to posterior)---------
#set bin minimum and maximum values based on minimum and maximum values for claw center y coordinates
min(claw.info$center.y)
bmin <- -5400
max(claw.info$center.y)
bmax <- 400

#plots histogram of claw center y coordinates 
ggplot(claw.info, aes(x= center.y))+
  stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
  geom_vline(xintercept = median(claw.info$center.y, col = 'red'))+
  theme_classic()+ylim(ymin,ymax)
ggsave('all claws y location.pdf',
       width = 6, height = 3,units = 'in')

#plots histogram of claw center y coordinates separating claws by clonal origin
clones <- c('a','b','c','d')
for(i in clones){
  ggplot(filter(claw.info, nb.cluster == i), aes(x= center.y))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(claw.info, nb.cluster == i)$center.y), col = 'red')+
    theme_classic()+ylim(ymin,ymax)+xlab(paste('clone',i,'y location'))
  
  ggsave(paste('clone',i,'_ylocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}

#plots histogram of claw center y coordinates separating claws by type 
types = c('gamma','abprime','ab')
for(i in types){
  ggplot(filter(claw.info,kc.largetype == i), aes(x= center.y))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(claw.info, kc.largetype == i)$center.y), col = 'red')+
    theme_classic()+ylim(ymin,ymax)+xlab(paste(i,'y location'))
  
  ggsave(paste(i,'_ylocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}

#plots histogram of claw center y coordinates separating by clone and type
for(i in clones){
  for(j in types){
    ggplot(filter(claw.info,kc.largetype == j, nb.cluster == i), aes(x= center.y))+
      stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
      geom_vline(xintercept = median(filter(claw.info, kc.largetype == j, nb.cluster == i)$center.y), col = 'red')+
      theme_classic()+ylim(ymin,ymax)+xlab(paste('clone',i,j,'y location'))
    
    ggsave(paste('clone',i,'_',j,'_ylocation.pdf', sep = ''),
           width = 6, height = 3,units = 'in')
    
  }
}
#histogram of claw centers z coordinate (ventral to dorsal-------
#set bin minimum and maximum values based on minimum and maximum values for claw center z coordinates
min(claw.info$center.z)
bmin <- -17100
max(claw.info$center.z)
bmax <- -12300

#plots histogram of claw center z coordinates
ggplot(claw.info, aes(x= center.z))+
  stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
  geom_vline(xintercept = median(claw.info$center.z, col = 'red'))+
  theme_classic()+ylim(ymin,ymax)
ggsave('all claws z location.pdf',
       width = 6, height = 3,units = 'in')

#plots histogram of claw center z coordinates separating claws by clonal origin
clones <- c('a','b','c','d')
for(i in clones){
  ggplot(filter(claw.info, nb.cluster == i), aes(x= center.z))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(claw.info, nb.cluster == i)$center.z), col = 'red')+
    theme_classic()+ylim(ymin,ymax)+xlab(paste('clone',i,'z location'))
  
  ggsave(paste('clone',i,'_zlocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}

#plots histogram of claw center z coordinates separating claws by type
types = c('gamma','abprime','ab')
for(i in types){
  ggplot(filter(claw.info,kc.largetype == i), aes(x= center.z))+
    stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
    geom_vline(xintercept = median(filter(claw.info, kc.largetype == i)$center.z), col = 'red')+
    theme_classic()+ylim(ymin,ymax)+xlab(paste(i,'z location'))
  
  ggsave(paste(i,'_zlocation.pdf', sep = ''),
         width = 6, height = 3,units = 'in')
}

#plots histogram of claw center z coordinates separating claws by clonal origin and type
for(i in clones){
  for(j in types){
    ggplot(filter(claw.info,kc.largetype == j, nb.cluster == i), aes(x= center.z))+
      stat_bin(geom = 'step', breaks = seq(bmin,bmax,100))+
      geom_vline(xintercept = median(filter(claw.info, kc.largetype == j, nb.cluster == i)$center.z), col = 'red')+
      theme_classic()+ylim(ymin,ymax)+xlab(paste('clone',i,j,'z location'))
    
    ggsave(paste('clone',i,'_',j,'_zlocation.pdf', sep = ''),
           width = 6, height = 3,units = 'in')
    
  }
}

