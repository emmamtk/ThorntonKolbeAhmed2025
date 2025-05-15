#branch point location analyses
#originally written 1/31/24
#edited 3/29/2024
#load packages------
library(natverse)
library(rlist)
library(dplyr)
library(ggbeeswarm)
theme_set(theme_classic())
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
cohensd <- function(x,y){
  #absolute value of difference in group means
  diffinmeans <- abs(mean(x)-mean(y))
  #pooled standard deviation- root mean square of the standard deviations
  pooledsd <- sqrt(((sd(x)^2)+(sd(y)^2))/2)
  #cohens d
  return(diffinmeans/pooledsd)
}
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
pn.info <- read.csv('PN_info.csv')
bout.info <- read.csv('bouton_info_collateral_bps.csv')
ca.pn <- list.load('ca_pn_skeletons_rotated.RDS')
mldv <- as.matrix(read.csv('mediallateral_DV_POV.csv')[,2:5])
#split ml axis into 4 equal parts-----------------
#largest ml coordinate(most lateral)
max(bout.info$bp.x)
#smallest ml coordinate (most medial)
min(bout.info$bp.x)
#bin width (differences between max and min/4)
bw <- (max(bout.info$bp.x)-min(bout.info$bp.x))/4

#make list of PNs with branch points ordered in L to M order-----
pnbranchpts <- list()
for(i in 1:nrow(pn.info)){
  pn <- pn.info$pn.bodyId[i]
  bps <- unique(bout.info$bp.x[bout.info$pn.bodyId== pn])
  bps <- bps[order(bps)]
  pnbranchpts[[i]] <- bps
  names(pnbranchpts)[i] <- pn
  pn.info$nbp[i] <- length(bps)
}

#add maximum spread of branch point from the center of the calyx--------
for(i in 1:nrow(pn.info)){
  ca.center <- median(unique(bout.info$bp.x))
  pn.info$maxdist[i] <- max(abs(bout.info$bp.x[bout.info$pn.bodyId == pn.info$pn.bodyId[i]]-ca.center))
}

#MAIN FIGURE: plot maximum branch point distance from center by number of branch points---------
#calculate number of branch points for each neuron
for(i in 1:nrow(pn.info)){
  pn.info$nbp[i] <- length(unique(filter(bout.info, pn.bodyId == pn.info$pn.bodyId[i])$bp))
}
#add factor to combine all neurons with 5 or more branch points
pn.info$nbpfactor <- pn.info$nbp
pn.info$nbpfactor[pn.info$nbp>=5] <- 'greater or equal to 5'
pn.info$nbpfactor <- as.factor(pn.info$nbpfactor)
#plot
ggplot(pn.info, aes(x=nbpfactor, y = maxdist, label = pn.glomerulus))+
  #quarters of the range of bp locations
  geom_hline(yintercept =c(bw,2*bw,3*bw), col = 'grey')+
  #data- each point is a neuron with the ML centroid of its bps and the ML distance between them
  geom_boxplot()+geom_beeswarm()+
  #labels points with glomerulus name
  geom_text()+
  ylab('maximum ML distance from calyx center')+xlab('number of branch points per neuron')

ggsave('branchpoint_MLmaxdist_allneurons.pdf',
       width = 20, height = 8,units = 'in') 

#how many neurons and glomeruli in each group for pn.info
a='greater or equal to 5'
length(unique(pn.info$pn.bodyId[pn.info$nbpfactor==a]))
length(unique(pn.info$pn.glomerulus[pn.info$nbpfactor==a]))

#plot neurons with 2 or fewer branch points and more than 2 branch points in different colors---------
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$nbp<=2]], col='brown2', lwd=5)
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$nbp>2]], col='cadetblue', lwd=5)
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('allbp.png',
             fmt = 'png')
#plotting neuron groups from other papers--------
#Zheng community PNs----------
#how many peripheral PNs are in the Zheng 'community'? 5/10
community <- c('DM2','DP1m','VM2','DL2v','DM3','DM4','DM1','VM3','VA2','VA4')
#how many community PNs in the hemibrain?
length(pn.info$pn.bodyId[pn.info$pn.glomerulus%in%community])
#how many of these have more than 2 bps and thus are peripheral? 6/15
sum(pn.info$pn.bodyId[pn.info$nbp>2]%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%community])
#how many of these have more than 2 bps and thus are peripheral? 6/15
sum(pn.info$pn.bodyId[pn.info$nbp<=2]%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%community])
#plot Zheng community PNs
clear3d()
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%community]], col='black', lwd=5)
plot3d(ca, add = T, alpha = 0.1)
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('zheng_ca.png',
             fmt = 'png')
#plot Zheng community PNs with peripheral PNs
clear3d()
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%community]], col='black', lwd=5)
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$nbp>2]], col='cadetblue', lwd=3)
plot3d(ca, add = T, alpha = 0.1)
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('zheng_peripheral_ca.png',
             fmt = 'png')

#plot Zheng community PNs with centeral PNs
clear3d()
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%community]], col='black', lwd=5)
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$nbp<=2]], col='brown2', lwd=3)
plot3d(ca, add = T, alpha = 0.1)
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('zheng_central_ca.png',
             fmt = 'png')

#tanaka peripheral PNs-------------
tanakap <- c('DM1','VA4','VC1','VM2')
#how many tanaka peripheral PNs?
length(pn.info$pn.bodyId[pn.info$pn.glomerulus%in%tanakap])
#how many of these have more than 2 bps and thus are peripheral? 6/15
sum(pn.info$pn.bodyId[pn.info$nbp>2]%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%tanakap])

#plot tanaka peripheral PNs
clear3d()
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%tanakap]], col='black', lwd=5)
plot3d(ca, add = T, alpha = 0.1)
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('tanaka_ca.png',
             fmt = 'png')
#plot tanaka peripheral PNs with pheripheral PNs
clear3d()
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%tanakap]], col='black', lwd=5)
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$nbp>2]], col='cadetblue', lwd=3)
plot3d(ca, add = T, alpha = 0.1)
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('tanaka_peripheral_ca.png',
             fmt = 'png')
#Plot tanaka peripheral PNs with central PNs
clear3d()
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$nbp<=2]], col='brown2', lwd=3)
plot3d(ca.pn[names(ca.pn)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus%in%tanakap]], col='black', lwd=5)
plot3d(ca, add = T, alpha = 0.1)
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('tanaka_central_ca.png',
             fmt = 'png')



#bp number by peripheral spread plot with tanaka and zheng community labeled------
ggplot(pn.info, aes(x=nbpfactor, y = maxdist, label = pn.glomerulus))+
  #quarters of the range of bp locations
  geom_hline(yintercept =c(bw,2*bw,3*bw), col = 'grey')+
  #data- each point is a neuron with the ML centroid of its bps and the ML distance between them
  geom_boxplot()+geom_beeswarm()+
  geom_beeswarm(data = pn.info[pn.info$pn.glomerulus%in%community,], aes(x=nbpfactor,y = maxdist), color = 'red')+
  geom_beeswarm(data = pn.info[pn.info$pn.glomerulus%in%tanakap,], aes(x=nbpfactor,y = maxdist), color = 'blue')+
  #labels points with glomerulus name
  geom_text()+
  ylab('maximum ML distance from calyx center')+xlab('number of branch points per neuron')

ggsave('branchpoint_MLmaxdist_allneurons_zheng_tanaka.pdf',
       width = 20, height = 8,units = 'in') 
