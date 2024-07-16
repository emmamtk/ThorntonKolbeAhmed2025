#KCs connected to each bouton vary with where a bouton is in space
#originally written 9/14/2023
#edited for clarity 11/14/2023
#load packages-----
library(dplyr)
library(natverse)
library(rlist)
library(ggplot2)
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
pn.skel <- list.load('ca_pn_skeletons_rotated.rds')#pn skeletons
ca <- list.load('calyxmesh_rotated.rds')#calyx mesh
apdv<- as.matrix.data.frame(read.csv('anteriorposterior_DV_POV.csv')[,2:5])#saved POV for calyx anterior posterior vs dorsal ventral
mldv<- as.matrix.data.frame(read.csv('mediallateral_DV_POV.csv')[,2:5]) #saved POV for calyx medial lateral vs dorsal ventral
claw.info <- read.csv('claw_info.csv')
bout.info <- read.csv('bouton_info_KCconnectivity.csv')[,1:8]
#get total number of claws per bouton, and claws separated by type and clone
for(i in 1:nrow(bout.info)){
  bclaws <- filter(claw.info, pn.boutId == bout.info$pn.boutId[i])
  bout.info$totalclaw[i] <- nrow(bclaws)
  bout.info$gamma[i] <- sum(bclaws$kc.largetype%in%'gamma')
  bout.info$abprime[i] <-  sum(bclaws$kc.largetype%in%'alphaprimebetaprime')
  bout.info$ab[i] <-  sum(bclaws$kc.largetype%in%'alphabeta')
  bout.info$clonea[i] <- sum(bclaws$nb.cluster%in%'a')
  bout.info$cloneb[i] <- sum(bclaws$nb.cluster%in%'b')
  bout.info$clonec[i] <- sum(bclaws$nb.cluster%in%'c')
  bout.info$cloned[i] <- sum(bclaws$nb.cluster%in%'d')
}
#filter out boutons that don't synapse with any olfactory KCs
bout.info <- filter(bout.info, totalclaw>0)


#example PN 1: DC1----------------------
#pn body Id for DC1 pn
dc1pn <- '1640594274'
#get number of claws from each type and each clone for each bouton from DC1
dc1 <- filter(bout.info, pn.bodyId == dc1pn)
#calculate type and clone totals so I can see which boutons have claws where all KC types and clonal origins are assigned
dc1$typetotal <- dc1$gamma+dc1$abprime+dc1$ab
dc1$clonetotal <- dc1$clonea+dc1$cloneb+dc1$clonec+dc1$cloned
dc1 <- dc1[,c(1:8,16,17,9:15)]
#choose two representative boutons
dc1bout1 <- '1640594274_4'
dc1bout2 <- '1640594274_9'

#example PN 2: DM3 (shares KCs with DC1)-------
#pn body id for DM3 PN
dm3pn <- '755518957'
#representative bouton (shares KCs with dc1bout1)
dm3bout3 <- '755518957_2'

#get number of claws from each type and each clone for each bouton from DM3-
dm3 <- filter(bout.info, pn.bodyId == dm3pn)
dm3$typetotal <- dm3$gamma+dm3$abprime+dm3$ab
dm3$clonetotal <- dm3$clonea+dm3$cloneb+dm3$clonec+dm3$cloned


#calyx with both PNs and their boutons---------------------
#boutons are made up of claw centers to give better dimension
#selected boutons are colored in red
#calyx mesh
plot3d(ca, alpha = 0.05, add = T)
#dc1 PN
plot3d(pn.skel[names(pn.skel)%in%dc1pn], col = 'grey25', lwd = 5)
#dm3 PN
plot3d(pn.skel[names(pn.skel)%in%dm3pn], col = 'grey45', lwd = 5)
#dc1 PN boutons
spheres3d(x= claw.info$center.x[claw.info$pn.boutId%in%dc1$pn.boutId[!dc1$pn.boutId%in%c(dc1bout1, dc1bout2)]],
       y= claw.info$center.y[claw.info$pn.boutId%in%dc1$pn.boutId[!dc1$pn.boutId%in%c(dc1bout1, dc1bout2)]],
       z= claw.info$center.z[claw.info$pn.boutId%in%dc1$pn.boutId[!dc1$pn.boutId%in%c(dc1bout1, dc1bout2)]],
       col = 'grey25',
       radius = 150,
       alpha = 0.7,
       lit = F,
       add = T)
#dc1 selected PN boutons
spheres3d(x= claw.info$center.x[claw.info$pn.boutId%in%dc1$pn.boutId[dc1$pn.boutId%in%c(dc1bout1, dc1bout2)]],
       y= claw.info$center.y[claw.info$pn.boutId%in%dc1$pn.boutId[dc1$pn.boutId%in%c(dc1bout1, dc1bout2)]],
       z= claw.info$center.z[claw.info$pn.boutId%in%dc1$pn.boutId[dc1$pn.boutId%in%c(dc1bout1, dc1bout2)]],
       col = 'red',
       radius = 150,
       alpha = 0.7,
       lit = F,
       add = T)
#dm3 PN boutons
spheres3d(x= claw.info$center.x[claw.info$pn.boutId%in%dm3$pn.boutId[!dm3$pn.boutId%in%dm3bout3]],
       y= claw.info$center.y[claw.info$pn.boutId%in%dm3$pn.boutId[!dm3$pn.boutId%in%dm3bout3]],
       z= claw.info$center.z[claw.info$pn.boutId%in%dm3$pn.boutId[!dm3$pn.boutId%in%dm3bout3]],
       col = 'grey45',
       radius = 150,
       alpha = 0.7,
       lit = F,
       add = T)
#dm3 selected PN bouton
spheres3d(x= claw.info$center.x[claw.info$pn.boutId%in%dm3$pn.boutId[dm3$pn.boutId%in%dm3bout3]],
       y= claw.info$center.y[claw.info$pn.boutId%in%dm3$pn.boutId[dm3$pn.boutId%in%dm3bout3]],
       z= claw.info$center.z[claw.info$pn.boutId%in%dm3$pn.boutId[dm3$pn.boutId%in%dm3bout3]],
       col = 'red',
       radius = 150,
       alpha = 0.7,
       lit = F,
       add = T)

#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('dc1_dm3_ap.png', fmt = 'png')

#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('dc1_dm3_ml.png', fmt = 'png')
#DM3 PN and boutons visible dc1 and calyx invisible-------------------
#calyx and dc1 are invisible so that the POV is the same and images can be overlayed later
plot3d(ca, alpha = 0, add = T)
plot3d(pn.skel[names(pn.skel)%in%dc1pn], col = 'black', lwd = 5)
plot3d(pn.skel[names(pn.skel)%in%dm3pn], col = 'grey45', lwd = 5,alpha = 0)
spheres3d(x= claw.info$center.x[claw.info$pn.bodyId == dc1pn],
       y= claw.info$center.y[claw.info$pn.bodyId == dc1pn],
       z= claw.info$center.z[claw.info$pn.bodyId == dc1pn],
       col = 'black',
       radius = 150,
       alpha = 0.7,
       lit = F,
       add = T)
spheres3d(x= claw.info$center.x[claw.info$pn.bodyId == dm3pn],
       y= claw.info$center.y[claw.info$pn.bodyId == dm3pn],
       z= claw.info$center.z[claw.info$pn.bodyId == dm3pn],
       col = 'grey45',
       radius = 150,
       alpha = 0,
       lit = F,
       add = T)

#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('dc1_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('dc1_ml.png', fmt = 'png')
#DC1 PN and boutons visible dm3 and calyx invisible-----------
#calyx and dm3 are invisible so that the POV is the same and images can be overlayed later
plot3d(ca, alpha = 0, add = T)
plot3d(pn.skel[names(pn.skel)%in%dc1pn], col = 'grey25', lwd = 5, alpha = 0)
plot3d(pn.skel[names(pn.skel)%in%dm3pn], col = 'black', lwd = 5)
plot3d(x= claw.info$center.x[claw.info$pn.bodyId == dc1pn],
       y= claw.info$center.y[claw.info$pn.bodyId == dc1pn],
       z= claw.info$center.z[claw.info$pn.bodyId == dc1pn],
       col = 'grey25',
       type= 'p',
       size = 20,
       add = T, alpha= 0)
spheres3d(x= claw.info$center.x[claw.info$pn.bodyId == dm3pn],
       y= claw.info$center.y[claw.info$pn.bodyId == dm3pn],
       z= claw.info$center.z[claw.info$pn.bodyId == dm3pn],
       color = 'black',
       radius = 150,
       alpha = 0.7,
       lit = F,
       add = T)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('dm3_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('dm3_ml.png', fmt = 'png')



#plot percent of a boutons partners that are of each KC type by that boutons coordinate location with line of best fit----
#info for selected boutons
selectbouts <- bout.info[bout.info$pn.boutId%in%c(dc1bout1,dc1bout2,dm3bout3),]
  

#percent of a boutons partners that are gamma kcs by boutons y (anterior posterior) coordinate
ggplot()+
  geom_point(data =bout.info, aes(x=center.y, y = gamma/totalclaw))+geom_smooth(data =bout.info, aes(x=center.y, y = gamma/totalclaw))+
  geom_point(data =selectbouts, aes(x=center.y, y = gamma/totalclaw),color='red')+
  theme_classic()
ggsave('boutons_gammashare_ap_withbouts.pdf',
       width = 6,
       units = 'in')
#percent of a boutons partners that are abprime kcs by boutons y (anterior posterior) coordinate
ggplot()+
  geom_point(data=bout.info, aes(x=center.y, y = abprime/totalclaw))+
  geom_smooth(data=bout.info, aes(x=center.y, y = abprime/totalclaw))+
  geom_point(data =selectbouts, aes(x=center.y, y = abprime/totalclaw),color='red')+
  theme_classic()
ggsave('boutons_abprimeshare_ap_withbouts.pdf',
       width = 6,
       units = 'in')
#percent of a boutons partners that are alpha beta kcs by boutons y (anterior posterior) coordinate
ggplot(bout.info, aes(x=center.y, y = ab/totalclaw))+
  geom_point(data = bout.info, aes(x=center.y, y = ab/totalclaw))+
  geom_smooth(data = bout.info, aes(x=center.y, y = ab/totalclaw))+
  geom_point(data =selectbouts, aes(x=center.y, y = ab/totalclaw),color='red')+
  theme_classic()
ggsave('boutons_abshare_ap_withbouts.pdf',
       width = 6,
       units = 'in')

#percent of a boutons partners that are clone a kcs by boutons x (medial lateral) coordinate
ggplot()+
  geom_point(data= bout.info, aes(x=center.x, y = clonea/totalclaw))+
  geom_smooth(data= bout.info, aes(x=center.x, y = clonea/totalclaw))+
  geom_point(data =selectbouts, aes(x=center.x, y = clonea/totalclaw),color='red')+
  theme_classic()+ylim(-.25,1)
ggsave('boutons_cloneashare_ml_withbouts.pdf',
       width = 6,
       units = 'in')
#percent of a boutons partners that are clone b kcs by boutons x (medial lateral) coordinate
ggplot()+
  geom_point(data= bout.info, aes(x=center.x, y = cloneb/totalclaw))+
  geom_smooth(data= bout.info, aes(x=center.x, y = cloneb/totalclaw))+
  geom_point(data =selectbouts, aes(x=center.x, y = cloneb/totalclaw),color='red')+
  theme_classic()+ylim(-.25,1)
ggsave('boutons_clonebshare_ml_withbouts.pdf',
       width = 6,
       units = 'in')
#percent of a boutons partners that are clone c kcs by boutons x (medial lateral) coordinate
ggplot()+
  geom_point(data= bout.info, aes(x=center.x, y = clonec/totalclaw))+
  geom_smooth(data= bout.info, aes(x=center.x, y = clonec/totalclaw))+
  geom_point(data =selectbouts, aes(x=center.x, y = clonec/totalclaw),color='red')+
  theme_classic()+ylim(-.25,1)
ggsave('boutons_clonecshare_ml_withbouts.pdf',
       width = 6,
       units = 'in')
#percent of a boutons partners that are clone d kcs by boutons x (medial lateral) coordinate
ggplot()+
  geom_point(data= bout.info, aes(x=center.x, y = cloned/totalclaw))+
  geom_smooth(data= bout.info, aes(x=center.x, y = cloned/totalclaw))+
  geom_point(data =selectbouts, aes(x=center.x, y = cloned/totalclaw),color='red')+
  theme_classic()+ylim(-.25,1)
ggsave('boutons_clonedshare_ml_withbouts.pdf',
       width = 6,
       units = 'in')

