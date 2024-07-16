#PN and KC skeleton snapshots
#originally written 12/4/2023
#load packages------
library(plotly)
library(natverse)
library(neuprintr)
library(dplyr)
library(rlist)
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/Data/connectome analysis redo/Figure Code")
full.pn <- list.load('pn_skeletons_rotated.rds')
full.kc <- list.load('kc_skeletons_rotated.rds')
ca.pn <- list.load('ca_pn_skeletons_rotated.rds')
ca.kc <- list.load('ca_kc_skeletons_rotated.rds')
ur.ca.pn <- list.load('ca_pn_skeletons.rds')
ur.ca.kc <- list.load('kc_skeletons_calyx.rds')
pn.info <- read.csv('PN_info.csv')
kc.info <- read.csv('KC_info.csv')
claw.info <- read.csv('claw_info.csv')
bout.info <- read.csv('bouton_info_KCconnectivity.csv')[1:12]
ca <- list.load('calyxmesh_rotated.rds')#calyx mesh
apdv<- as.matrix.data.frame(read.csv('anteriorposterior_DV_POV.csv')[,2:5])#saved POV for calyx anterior posterior vs dorsal ventral
mldv<- as.matrix.data.frame(read.csv('mediallateral_DV_POV.csv')[,2:5]) #saved POV for calyx medial lateral vs dorsal ventral
fullapdv <- as.matrix.data.frame(read.csv('fullskeleton_anteriorposterior_DV_POV.csv')[,2:5])
fullmldv <- as.matrix.data.frame(read.csv('fullskeleton_mediallateral_DV_POV.csv')[,2:5])
starryview <- as.matrix.data.frame(read.csv('calyx_starrynight_POV.csv')[,2:5])
#filter skeletons to only ones being included in the analyses
ca.pn <- ca.pn[names(ca.pn)%in%pn.info$pn.bodyId]
ca.kc <- ca.kc[names(ca.kc)%in%kc.info$kc.bodyId]
full.pn <- full.pn[names(full.pn)%in%pn.info$pn.bodyId]
full.kc <- full.kc[names(full.kc)%in%kc.info$kc.bodyId]
#setting POVS for calyx-------
#apdv <- par3d()$userMatrix
#write.csv(apdv, 'anteriorposterior_DV_POV.csv')
#mldv <- par3d()$userMatrix
#write.csv(mldv, 'mediallateral_DV_POV.csv')
#setting POVs for full skeletons-----
#full skeleton ML DV view
#fullmldv <- par3d()$userMatrix
#write.csv(fullmldv, 'fullskeleton_mediallateral_DV_POV.csv')
#full skeleton AP DV view
#fullapdv <- par3d()$userMatrix
#write.csv(fullapdv, 'fullskeleton_anteriorposterior_DV_POV.csv')
#color palettes------
#pns- 51 different colors
pn.pal <- c('yellow','wheat', 'violetred', 'violet','turquoise','tomato','thistle','tan','springgreen',
         'slategrey','slateblue','skyblue','sienna','seagreen','sandybrown','salmon','grey','royalblue','rosybrown',
         'red','purple','powderblue','plum','pink','peru','palevioletred','paleturquoise','palegreen','orchid','orangered',
         'orange','moccasin','mediumpurple','mediumorchid','maroon','limegreen','lightsteelblue','lightseagreen',
         'lightpink','lavender','lightcoral','khaki','hotpink','indianred','green','gold','chocolate','cadetblue',
         'brown','aquamarine','darkmagenta')
#kcs- 3 different colors
#gamma, alpha'beta', alpha beta
kc.pal <- c('#020F5C','#517802','#4787ED')

#PN soma-----
#add soma point to pn.info
for(i in 1:length(pn.info$pn.bodyId)){
  pn.info$soma.pt[i] <- full.pn[names(full.pn)%in%pn.info$pn.bodyId[i]][[1]]$soma
}
pn.info$soma.pt <- unlist(pn.info$soma.pt)
#add x,y,z coordinates for the point specified as soma
for(i in 1:length(pn.info$pn.bodyId)){
  pn.info$soma.x[i] <- full.pn[names(full.pn)%in%pn.info$pn.bodyId[i]][[1]]$d[pn.info$soma.pt[i],3]
  pn.info$soma.y[i] <- full.pn[names(full.pn)%in%pn.info$pn.bodyId[i]][[1]]$d[pn.info$soma.pt[i],4]
  pn.info$soma.z[i] <- full.pn[names(full.pn)%in%pn.info$pn.bodyId[i]][[1]]$d[pn.info$soma.pt[i],5]
}

#kc soma------
#add soma point to kc.info
for(i in 1:length(kc.info$kc.bodyId)){
  kc.info$soma.pt[i] <- full.kc[names(full.kc)%in%kc.info$kc.bodyId[i]][[1]]$soma
}

kc.info$soma.pt <- unlist(kc.info$soma.pt)
#add x,y,z coordinates for the point specified as soma
for(i in 1:length(kc.info$kc.bodyId)){
  kc.info$soma.x[i] <- full.kc[names(full.kc)%in%kc.info$kc.bodyId[i]][[1]]$d[kc.info$soma.pt[i],3]
  kc.info$soma.y[i] <- full.kc[names(full.kc)%in%kc.info$kc.bodyId[i]][[1]]$d[kc.info$soma.pt[i],4]
  kc.info$soma.z[i] <- full.kc[names(full.kc)%in%kc.info$kc.bodyId[i]][[1]]$d[kc.info$soma.pt[i],5]
}
#full KC skeletons colored by type----
plot3d(full.kc[kc.info$kc.bodyId[kc.info$kc.type=='kcg']], 
       col = kc.pal[1])
spheres3d(x= kc.info$soma.x[kc.info$kc.type=='kcg'],
          y= kc.info$soma.y[kc.info$kc.type=='kcg'],
          z= kc.info$soma.z[kc.info$kc.type=='kcg'],
          col = kc.pal[1],
          radius = 375,
          alpha = 0.7,
          lit = F,
          add = T)
plot3d(full.kc[kc.info$kc.bodyId[kc.info$kc.largetype=='alphaprimebetaprime']], 
       col = kc.pal[2])
spheres3d(x= kc.info$soma.x[kc.info$kc.largetype=='alphaprimebetaprime'],
          y= kc.info$soma.y[kc.info$kc.largetype=='alphaprimebetaprime'],
          z= kc.info$soma.z[kc.info$kc.largetype=='alphaprimebetaprime'],
          col = kc.pal[2],
          radius = 375,
          alpha = 0.7,
          lit = F,
          add = T)
plot3d(full.kc[kc.info$kc.bodyId[kc.info$kc.largetype=='alphabeta']], 
       col = kc.pal[3])
spheres3d(x= kc.info$soma.x[kc.info$kc.largetype=='alphabeta'],
          y= kc.info$soma.y[kc.info$kc.largetype=='alphabeta'],
          z= kc.info$soma.z[kc.info$kc.largetype=='alphabeta'],
          col = kc.pal[3],
          radius = 375,
          alpha = 0.7,
          lit = F,
          add = T)
plot3d(full.pn, alpha = 0)
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view without axes
view3d(userMatrix = fullapdv, zoom =.5)
rgl.snapshot('fullkc_ap.png', fmt = 'png')
#medial lateral view without axes
view3d(userMatrix = fullmldv, zoom =.5)
rgl.snapshot('fullkc_ml.png', fmt = 'png')

#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view with axes
view3d(userMatrix = fullapdv, zoom =.5)
rgl.snapshot('fullkc_ap_axes.png', fmt = 'png')
#medial lateral view with axes
view3d(userMatrix = fullmldv, zoom =.5)
rgl.snapshot('fullkc_ml_axes.png', fmt = 'png')

#full PN skeletons colored by type----
gloms <- unique(pn.info$pn.glomerulus)
for(i in 1:length(gloms)){
  pns <- pn.info$pn.bodyId[pn.info$pn.glomerulus == gloms[i]]
  plot3d(full.pn[names(full.pn)%in%pns], col=pn.pal[i])
  spheres3d(x= pn.info$soma.x[pn.info$pn.glomerulus == gloms[i]],
         y= pn.info$soma.y[pn.info$pn.glomerulus == gloms[i]],
         z= pn.info$soma.z[pn.info$pn.glomerulus == gloms[i]],
         col = pn.pal [i],
         radius = 850,
         alpha = 0.0,
         lit = F,
         add = T)
}
#invisible KCs---
plot3d(full.kc, alpha = 0)
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view without axes
view3d(userMatrix = fullapdv, zoom =.5)
rgl.snapshot('fullpn_ap.png', fmt = 'png')
#medial lateral view without axes
view3d(userMatrix = fullmldv, zoom =.5)
rgl.snapshot('fullpn_ml.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view with axes
view3d(userMatrix = fullapdv, zoom =.5)
rgl.snapshot('fullpn_ap_axes.png', fmt = 'png')
#medial lateral view with axes
view3d(userMatrix = fullmldv, zoom =.5)
rgl.snapshot('fullpn_ml_axes.png', fmt = 'png')
#example single full PN skeletons--------
plot3d(full.pn[names(full.pn)%in%'5813039235'])
#invisible PNs
plot3d(full.pn, alpha = 0)
#invisible KCs
plot3d(full.kc, alpha = 0)
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#medial lateral view without axes
view3d(userMatrix = fullmldv, zoom =.5)
rgl.snapshot('example_fullpn_ml.png', fmt = 'png')
#full PN and KC skeletons colored by types----
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view without axes
view3d(userMatrix = fullapdv, zoom =.5)
rgl.snapshot('fullpnkc_ap.png', fmt = 'png')
#medial lateral view without axes
view3d(userMatrix = fullmldv, zoom =.5)
rgl.snapshot('fullpnkc_ml.png', fmt = 'png')

#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view with axes
view3d(userMatrix = fullapdv, zoom =.5)
rgl.snapshot('fullpnkc_ap_axes.png', fmt = 'png')
#medial lateral view with axes
view3d(userMatrix = fullmldv, zoom =.5)
rgl.snapshot('fullpnkc_ml_axes.png', fmt = 'png')
#calyx KCs colored by type PNs invisible so angle is the same-----
plot3d(ca.kc[names(ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.type=='kcg']], 
       col = kc.pal[1], add = T)
spheres3d(x= kc.info$soma.x[kc.info$kc.type=='kcg'],
          y= kc.info$soma.y[kc.info$kc.type=='kcg'],
          z= kc.info$soma.z[kc.info$kc.type=='kcg'],
          col = kc.pal[1],
          radius = 300,
          alpha = 0.7,
          lit = F,
          add = T)
plot3d(ca.kc[names(ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.largetype=='alphaprimebetaprime']], 
       col = kc.pal[2], add = T)
spheres3d(x= kc.info$soma.x[kc.info$kc.largetype=='alphaprimebetaprime'],
          y= kc.info$soma.y[kc.info$kc.largetype=='alphaprimebetaprime'],
          z= kc.info$soma.z[kc.info$kc.largetype=='alphaprimebetaprime'],
          col = kc.pal[2],
          radius = 300,
          alpha = 0.7,
          lit = F,
          add = T)
plot3d(ca.kc[names(ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.largetype=='alphabeta']], 
       col = kc.pal[3], add = T)
spheres3d(x= kc.info$soma.x[kc.info$kc.largetype=='alphabeta'],
          y= kc.info$soma.y[kc.info$kc.largetype=='alphabeta'],
          z= kc.info$soma.z[kc.info$kc.largetype=='alphabeta'],
          col = kc.pal[3],
          radius = 300,
          alpha = 0.7,
          lit = F,
          add = T)

plot3d(ca.pn[names(ca.pn)], alpha=0)
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.4)
rgl.snapshot('ca_kc_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.4)
rgl.snapshot('ca_kc_ml.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.4)
rgl.snapshot('ca_kc_ap_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.4)
rgl.snapshot('ca_kc_ml_axes.png', fmt = 'png')

#calyx PNs colored by type-----
clear3d()
gloms <- unique(pn.info$pn.glomerulus)
for(i in 1:length(gloms)){
  pns <- pn.info$pn.bodyId[pn.info$pn.glomerulus == gloms[i]]
  plot3d(ca.pn[names(ca.pn)%in%pns], col=pn.pal[i], lwd=2)
  spheres3d(x= claw.info$center.x[claw.info$pn.glomerulus == gloms[i]],
         y= claw.info$center.y[claw.info$pn.glomerulus == gloms[i]],
         z= claw.info$center.z[claw.info$pn.glomerulus == gloms[i]],
         col = pn.pal [i],
         radius = 150,
         lit= F,
         alpha = 0.5)
}
#invisible KCs
plot3d(ca.kc[names(ca.kc)], alpha = 0)
spheres3d(x= kc.info$soma.x,
          y= kc.info$soma.y,
          z= kc.info$soma.z,
          radius = 300,
          alpha = 0,
          add = T)

#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.4)
rgl.snapshot('ca_pn_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.4)
rgl.snapshot('ca_pn_ml.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.4)
rgl.snapshot('ca_pn_ap_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.4)
rgl.snapshot('ca_pn_ml_axes.png', fmt = 'png')


#counts of different features----
#olfactory PNs-105
nrow(pn.info)
#olfactory KCs-1728
nrow(kc.info)
#clone A KCs 450
nrow(filter(kc.info, nb.cluster =='a'))
#clone B KCs 383
nrow(filter(kc.info, nb.cluster =='b'))
#clone C KCs- 424
nrow(filter(kc.info, nb.cluster =='c'))
#Clone D KCs- 415
nrow(filter(kc.info, nb.cluster =='d'))
#gamma KCs 605
nrow(filter(kc.info, kc.largetype =='gamma'))
#a'b' KCs-246
nrow(filter(kc.info, kc.largetype =='alphaprimebetaprime'))
#ab KCs-799
nrow(filter(kc.info, kc.largetype =='alphabeta'))
#Clone A g-182, a'b'-60, ab-208
nrow(filter(kc.info, nb.cluster == 'a', kc.largetype =='gamma'))
nrow(filter(kc.info, nb.cluster == 'a', kc.largetype =='alphaprimebetaprime'))
nrow(filter(kc.info, nb.cluster == 'a', kc.largetype =='alphabeta'))
#Clone B g-113, a'b'-64, ab-206
nrow(filter(kc.info, nb.cluster == 'b', kc.largetype =='gamma'))
nrow(filter(kc.info, nb.cluster == 'b', kc.largetype =='alphaprimebetaprime'))
nrow(filter(kc.info, nb.cluster == 'b', kc.largetype =='alphabeta'))
#Clone C g-158, a'b'-61, ab-204
nrow(filter(kc.info, nb.cluster == 'c', kc.largetype =='gamma'))
nrow(filter(kc.info, nb.cluster == 'c', kc.largetype =='alphaprimebetaprime'))
nrow(filter(kc.info, nb.cluster == 'c', kc.largetype =='alphabeta'))
#Clone D g-151, a'b'-61, ab-181
nrow(filter(kc.info, nb.cluster == 'd', kc.largetype =='gamma'))
nrow(filter(kc.info, nb.cluster == 'd', kc.largetype =='alphaprimebetaprime'))
nrow(filter(kc.info, nb.cluster == 'd', kc.largetype =='alphabeta'))
#total KC claws-10,195
nrow(claw.info)
#total boutons-442
nrow(bout.info)
#average presynaptic sites per bouton- 56/ 50
mean(bout.info$nsyn)
median(bout.info$nsyn)
#average claws per bouton 24/22
mean(bout.info$nkc)
median(bout.info$nkc)
#average post synaptic sites per claw-16.45/16
mean(claw.info$npostsyn)
median(claw.info$npostsyn)
#invisible calyx PNs with one PN visible----
chosenpn <- '5813039315'
sisterpn <- '754534424'
strangerpn <- '5813039235'

plot3d(ca.pn[names(ca.pn)%in%chosenpn], col='black', lwd=4)
plot3d(ca.pn[names(ca.pn)%in%sisterpn], col='grey34', lwd=4)
plot3d(ca.pn[names(ca.pn)%in%strangerpn], col='red', lwd=4)
plot3d(ca.pn[names(ca.pn)], alpha = 0)
plot3d(ca.kc, alpha=0)
spheres3d(x= claw.info$center.x[claw.info$pn.bodyId == chosenpn],
       y= claw.info$center.y[claw.info$pn.bodyId == chosenpn],
       z= claw.info$center.z[claw.info$pn.bodyId == chosenpn],
       col = 'black',
       radius = 150,
       lit= F,
       alpha = 0.5,
       add = T)
spheres3d(x= claw.info$center.x[claw.info$pn.bodyId == sisterpn],
          y= claw.info$center.y[claw.info$pn.bodyId == sisterpn],
          z= claw.info$center.z[claw.info$pn.bodyId == sisterpn],
          col = 'grey34',
          radius = 150,
          lit= F,
          alpha = 0.5,
          add = T)
spheres3d(x= claw.info$center.x[claw.info$pn.bodyId == strangerpn],
          y= claw.info$center.y[claw.info$pn.bodyId == strangerpn],
          z= claw.info$center.z[claw.info$pn.bodyId == strangerpn],
          col = 'red',
          radius = 150,
          lit= F,
          alpha = 0.5,
          add = T)
plot3d(ca.kc[names(ca.kc)%in%'754547386'], col='grey20', lwd=2)
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.4)
rgl.snapshot('ca_cognate_da1_dm6_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.4)
rgl.snapshot('ca_cognate_da1_dm6_ml.png', fmt = 'png')
#stary calyx view
view3d(userMatrix = starryview, zoom =.3)
rgl.snapshot('ca_pnkc_starry.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.4)
rgl.snapshot('ca_cognate_da1_dm6_ap_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.4)
rgl.snapshot('ca_cognate_da1_dm6_ml_axes.png', fmt = 'png')

#starry calyx-------
starrykc <- c('#121020','#141c3a','#324c5c','#7ea4b0','#4c6394')
plot3d(ca.kc, 
       col = rep(starrykc,369), add = T, lwd =2)
spheres3d(x= kc.info$soma.x,
          y= kc.info$soma.y,
          z= kc.info$soma.z,
          col = rep(starrykc, 369),
          radius = 250,
          alpha = 0.7,
          lit = F,
          add = T)
starrypn <- c('#cdd27e','#b4ad63','#fddf2c','#f4fdad')
for(i in 1:length(gloms)){
  pns <- pn.info$pn.bodyId[pn.info$pn.glomerulus == gloms[i]]
  plot3d(ca.pn[names(ca.pn)%in%pns], col=rep(starrypn,11)[i], lwd=2)
  spheres3d(x= claw.info$center.x[claw.info$pn.glomerulus == gloms[i]],
            y= claw.info$center.y[claw.info$pn.glomerulus == gloms[i]],
            z= claw.info$center.z[claw.info$pn.glomerulus == gloms[i]],
            col = rep(starrypn,15)[i],
            radius = 75,
            lit= F,
            alpha = 1)
}
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.3)
rgl.snapshot('ca_starry_ap_nosoma.png', fmt = 'png')
#best starry night view
view3d(userMatrix = starryview, zoom =.4)
rgl.snapshot('ca_starry_starry_soma.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.3)
rgl.snapshot('ca_starry_ap_axes_nosoma.png', fmt = 'png')
#anterior posterior view
view3d(userMatrix =starryview, zoom =.3)
rgl.snapshot('ca_starry_starry_axes_nosoma.png', fmt = 'png')
#defining best starry night view
#starryview <- par3d()$userMatrix
#write.csv(starryview, 'calyx_starrynight_POV.csv')

#rotated and unrotated calyx--------------
#unrotated KCs------------------
plot3d(ur.ca.kc[names(ur.ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.type=='kcg']], 
       col = kc.pal[1], add = T)

plot3d(ur.ca.kc[names(ur.ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.largetype=='alphaprimebetaprime']], 
       col = kc.pal[2], add = T)

plot3d(ur.ca.kc[names(ur.ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.largetype=='alphabeta']], 
       col = kc.pal[3], add = T)

plot3d(ur.ca.pn[names(ur.ca.pn)], alpha=0)
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('ur_ca_kc_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('ur_ca_kc_ml.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('ur_ca_kc_ap_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('ur_ca_kc_ml_axes.png', fmt = 'png')

#unrotated PNs------------
clear3d()
gloms <- unique(pn.info$pn.glomerulus)
for(i in 1:length(gloms)){
  pns <- pn.info$pn.bodyId[pn.info$pn.glomerulus == gloms[i]]
  plot3d(ur.ca.pn[names(ur.ca.pn)%in%pns], col=pn.pal[i], lwd=2)
}
#invisible KCs
plot3d(ur.ca.kc[names(ur.ca.kc)], alpha = 0)

#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('ur_ca_pn_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('ur_ca_pn_ml.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('ur_ca_pn_ap_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('ur_ca_pn_ml_axes.png', fmt = 'png')


#rotated KCs------------------
plot3d(ca.kc[names(ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.type=='kcg']], 
       col = kc.pal[1], add = T)

plot3d(ca.kc[names(ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.largetype=='alphaprimebetaprime']], 
       col = kc.pal[2], add = T)

plot3d(ca.kc[names(ca.kc)%in%kc.info$kc.bodyId[kc.info$kc.largetype=='alphabeta']], 
       col = kc.pal[3], add = T)

plot3d(ca.pn[names(ca.pn)], alpha=0)
#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('ca_kc_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('ca_kc_ml.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('ca_kc_ap_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('ca_kc_ml_axes.png', fmt = 'png')

#unrotated PNs------------
clear3d()
gloms <- unique(pn.info$pn.glomerulus)
for(i in 1:length(gloms)){
  pns <- pn.info$pn.bodyId[pn.info$pn.glomerulus == gloms[i]]
  plot3d(ca.pn[names(ca.pn)%in%pns], col=pn.pal[i], lwd=2)
}
#invisible KCs
plot3d(ca.kc[names(ca.kc)], alpha = 0)

#invisible axes
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('ca_pn_ap.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('ca_pn_ml.png', fmt = 'png')
#visible axes
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.6)
rgl.snapshot('ca_pn_ap_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.6)
rgl.snapshot('ca_pn_ml_axes.png', fmt = 'png')


#calyx with tanaka peripheral PNs labeled
#calyx with food fovea PNs labeled