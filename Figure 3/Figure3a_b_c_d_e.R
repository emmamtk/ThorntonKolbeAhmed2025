#Random Sampling Figure Graphs
#originally written 2/26/2024
#load packages------
library(dplyr)
library(rlist)
library(ggplot2)
library(natverse)
library(plotly)
library(ggridges)
library(ggbeeswarm)
library(FSA)
theme_set(theme_classic())
#load data------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper/10k random model conditional input analysis mean sd zscore")
ana.list <- list.load('coninp_analysis.RDS')
pn.skel <- list.load('ca_pn_skeletons_rotated.RDS')

setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
pn.info <- read.csv('PN_info.csv')
claw.info <- read.csv('claw_info.csv')
bout.info <- read.csv('bouton_info_KCconnectivity.csv')

ca <- list.load('calyxmesh_rotated.rds')#calyx mesh
apdv<- as.matrix.data.frame(read.csv('anteriorposterior_DV_POV.csv')[,2:5])#saved POV for calyx anterior posterior vs dorsal ventral
mldv<- as.matrix.data.frame(read.csv('mediallateral_DV_POV.csv')[,2:5]) #saved POV for calyx medial lateral vs dorsal ventral
############Figure 3 A####
#plot boutons centers colored by glomerulus type-----------
pn.pal <- c('yellow','wheat', 'violetred', 'violet','turquoise','tomato','thistle','tan','springgreen',
            'slategrey','slateblue','skyblue','sienna','seagreen','sandybrown','salmon','grey','royalblue','rosybrown',
            'red','purple','powderblue','plum','pink','peru','palevioletred','paleturquoise','palegreen','orchid','orangered',
            'orange','moccasin','mediumpurple','mediumorchid','maroon','limegreen','lightsteelblue','lightseagreen',
            'lightpink','lavender','lightcoral','khaki','hotpink','indianred','green','gold','chocolate','cadetblue',
            'brown','aquamarine','darkmagenta')

gloms <- unique(pn.info$pn.glomerulus)
#without calyx and axes
for(i in 1:length(gloms)){
  spheres3d(x= bout.info$center.x[bout.info$pn.glomerulus == gloms[i]],
            y= bout.info$center.y[bout.info$pn.glomerulus == gloms[i]],
            z= bout.info$center.z[bout.info$pn.glomerulus == gloms[i]],
            col = pn.pal [i],
            radius = 200,
            alpha = 1,
            lit = F,
            add = T)
}
plot3d(ca, alpha = 0, add= T)
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.5)
rgl.snapshot('allboutons_AP.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.5)
rgl.snapshot('allboutons_ML.png', fmt = 'png')
#with calyx and axes
clear3d()
for(i in 1:length(gloms)){
  spheres3d(x= bout.info$center.x[bout.info$pn.glomerulus == gloms[i]],
            y= bout.info$center.y[bout.info$pn.glomerulus == gloms[i]],
            z= bout.info$center.z[bout.info$pn.glomerulus == gloms[i]],
            col = pn.pal [i],
            radius = 200,
            alpha = 1,
            lit = F,
            add = T)
}
plot3d(ca, alpha = 0.1, add= T)
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.5)
rgl.snapshot('allboutons_AP_ca_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.5)
rgl.snapshot('allboutons_ML_ca_axes.png', fmt = 'png')

#plot bouton centers included in each KC group model colored by glomerulus type--------

#dont need to include calyx or axes because it will be the same as the all boutons one and I will trace the calyx anyway
#for kc types
for(j in c('gamma','alphaprimebetaprime','alphabeta')){
  jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, kc.largetype == j)$pn.boutId),]
  for(i in 1:length(gloms)){
    spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == gloms[i]],
              y= jbout.info$center.y[jbout.info$pn.glomerulus == gloms[i]],
              z= jbout.info$center.z[jbout.info$pn.glomerulus == gloms[i]],
              col = pn.pal [i],
              radius = 200,
              alpha = 1,
              lit = F,
              add = T)
  }
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,'_boutons_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,'_boutons_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
}
#for nb clones
for(j in c('a','b','c','d')){
  jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, nb.cluster == j)$pn.boutId),]
  for(i in 1:length(gloms)){
    spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == gloms[i]],
              y= jbout.info$center.y[jbout.info$pn.glomerulus == gloms[i]],
              z= jbout.info$center.z[jbout.info$pn.glomerulus == gloms[i]],
              col = pn.pal [i],
              radius = 200,
              alpha = 1,
              lit = F,
              add = T)
  }
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,'_boutons_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,'_boutons_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
}
#for kc clone and types
for(k in c('gamma','alphaprimebetaprime','alphabeta')){
  for(j in c('a','b','c','d')){
    jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, nb.cluster == j, kc.largetype == k)$pn.boutId),]
    for(i in 1:length(gloms)){
      spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == gloms[i]],
                y= jbout.info$center.y[jbout.info$pn.glomerulus == gloms[i]],
                z= jbout.info$center.z[jbout.info$pn.glomerulus == gloms[i]],
                col = pn.pal [i],
                radius = 200,
                alpha = 1,
                lit = F,
                add = T)
    }
    #without calyx and axes
    plot3d(ca, alpha = 0, add= T)
    axes3d(edges = 'bbox',alpha=0)
    #anterior posterior view
    view3d(userMatrix = apdv, zoom =.5)
    rgl.snapshot(paste(k,j,'_boutons_AP.png',sep=''), fmt = 'png')
    #medial lateral view
    view3d(userMatrix = mldv, zoom =.5)
    rgl.snapshot(paste(k,j,'_boutons_ML.png',sep=''), fmt = 'png')
    #clear before next KC group
    clear3d()
  }
}


#plot KCs claw distribution for each KC group---------
#without calyx and axes
spheres3d(x= claw.info$center.x,
            y= claw.info$center.y,
            z= claw.info$center.z,
            radius = 25,
            alpha = 1,
            lit = F,
            add = T)
plot3d(ca, alpha = 0, add= T)
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.5)
rgl.snapshot('allclaws_AP.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.5)
rgl.snapshot('allclaws_ML.png', fmt = 'png')

#with calyx and axes
clear3d()
spheres3d(x= claw.info$center.x,
          y= claw.info$center.y,
          z= claw.info$center.z,
          radius = 25,
          alpha = 1,
          lit = F,
          add = T)
plot3d(ca, alpha = 0.1, add= T)
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.5)
rgl.snapshot('allclaws_AP_ca_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.5)
rgl.snapshot('allclaws_ML_ca_axes.png', fmt = 'png')


#for kc types
for(j in c('gamma','alphaprimebetaprime','alphabeta')){
  spheres3d(x= claw.info$center.x[claw.info$kc.largetype==j],
            y= claw.info$center.y[claw.info$kc.largetype==j],
            z= claw.info$center.z[claw.info$kc.largetype==j],
            radius = 25,
            alpha = 1,
            lit = F,
            add = T)
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,'_claws_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,'_claws_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
}
#for clones
for(j in c('a','b','c','d')){
  spheres3d(x= claw.info$center.x[claw.info$nb.cluster==j],
            y= claw.info$center.y[claw.info$nb.cluster==j],
            z= claw.info$center.z[claw.info$nb.cluster==j],
            radius = 25,
            alpha = 1,
            lit = F,
            add = T)
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,'_claws_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,'_claws_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
}
#for clones
for(k in c('gamma','alphaprimebetaprime','alphabeta')){
  for(j in c('a','b','c','d')){
    spheres3d(x= claw.info$center.x[claw.info$nb.cluster==j&claw.info$kc.largetype==k],
              y= claw.info$center.y[claw.info$nb.cluster==j&claw.info$kc.largetype==k],
              z= claw.info$center.z[claw.info$nb.cluster==j&claw.info$kc.largetype==k],
              radius = 25,
              alpha = 1,
              lit = F,
              add = T)
    #without calyx and axes
    plot3d(ca, alpha = 0, add= T)
    axes3d(edges = 'bbox',alpha=0)
    #anterior posterior view
    view3d(userMatrix = apdv, zoom =.5)
    rgl.snapshot(paste(j,k,'_claws_AP.png',sep=''), fmt = 'png')
    #medial lateral view
    view3d(userMatrix = mldv, zoom =.5)
    rgl.snapshot(paste(j,k,'_claws_ML.png',sep=''), fmt = 'png')
    #clear before next KC group
    clear3d()
  }
}

################Figure 3 B####
#histogram of z scores by type-----
allz <- data.frame(zscores=as.vector(ana.list$all$zscore))
gammaz <- data.frame(zscores=as.vector(ana.list$gamma$zscore))
abprimez <- data.frame(zscores=as.vector(ana.list$abprime$zscore))
abz <- data.frame(zscores=as.vector(ana.list$ab$zscore))

ggplot()+theme_classic()+
  geom_vline(xintercept=seq(-5,5,1), color= 'grey')+
  geom_histogram(data= allz, aes(x=zscores), binwidth = 0.5, color = 'black')+
  geom_vline(xintercept = median(allz$zscores))+
  geom_histogram(data = gammaz, aes(x= zscores), binwidth = 0.5, color = 'darkblue')+
  geom_vline(xintercept = median(gammaz$zscores),color='darkblue')+
  geom_histogram(data = abprimez, aes(x= zscores), binwidth = 0.5, color = 'green')+
  geom_vline(xintercept = median(abprimez$zscores),color='green')+
  geom_histogram(data = abz, aes(x= zscores), binwidth = 0.5, color = 'cadetblue')+
  geom_vline(xintercept = median(abz$zscores),color='cadetblue')+
  ylim(0,1500)+xlim(-5.5,5.5)+
  ggtitle('blue gamma, green abprime lblue ab')
ggsave('zscoresbykctype_5z.pdf',
       width = 6, height = 3,units = 'in')

#histogram of z scores by clone------
cloneaz <- data.frame(zscores=as.vector(ana.list$clonea$zscore))
clonebz <- data.frame(zscores=as.vector(ana.list$cloneb$zscore))
clonecz <- data.frame(zscores=as.vector(ana.list$clonec$zscore))
clonedz <- data.frame(zscores=as.vector(ana.list$cloned$zscore))
ggplot()+theme_classic()+
  geom_vline(xintercept=seq(-5,5,1), color= 'grey')+
  geom_histogram(data= allz, aes(x=zscores), binwidth = 0.5)+
  geom_vline(xintercept = median(allz$zscores))+
  geom_histogram(data= cloneaz, aes(x=zscores), binwidth = 0.5, color = 'red')+
  geom_vline(xintercept = median(cloneaz$zscores), color = 'red')+
  geom_histogram(data= clonebz, aes(x=zscores), binwidth = 0.5, color = 'orange')+
  geom_vline(xintercept = median(clonebz$zscores), color = 'orange')+
  geom_histogram(data= clonecz, aes(x=zscores), binwidth = 0.5, color = 'gold')+
  geom_vline(xintercept = median(clonecz$zscores), color = 'gold')+
  geom_histogram(data= clonedz, aes(x=zscores),binwidth = 0.5, color = 'green')+
  geom_vline(xintercept = median(clonedz$zscores),color = 'green')+
  ylim(0,1500)+xlim(-5.5,5.5)+
  ggtitle('red a, ora b, gold c, green d')
ggsave('zscoresbykcclone_5z.pdf',
       width = 6, height = 3,units = 'in')
#histogram of z scores by clone and type-----
#clone A
cloneagammaz <- data.frame(zscores=as.vector(ana.list$cloneagamma$zscore))
cloneaabprimez <- data.frame(zscores=as.vector(ana.list$cloneaabprime$zscore))
cloneaabz <- data.frame(zscores=as.vector(ana.list$cloneaab$zscore))
ggplot()+theme_classic()+
  geom_vline(xintercept=seq(-5,5,1), color= 'grey')+
  geom_histogram(data= allz, aes(x=zscores), binwidth = 0.5, color = 'black')+
  geom_vline(xintercept = median(allz$zscores))+
  geom_histogram(data= cloneagammaz, aes(x=zscores), binwidth = 0.5, color = 'darkblue')+
  geom_vline(xintercept = median(cloneagammaz$zscores),color='darkblue')+
  geom_histogram(data= cloneaabprimez, aes(x=zscores), binwidth = 0.5, color = 'green')+
  geom_vline(xintercept = median(cloneaabprimez$zscores),color='green')+
  geom_histogram(data= cloneaabz, aes(x=zscores), binwidth = 0.5, color = 'cadetblue')+
  geom_vline(xintercept = median(cloneaabz$zscores),color='cadetblue')+
  ylim(0,1500)+xlim(-5.5,5.5)+
  ggtitle('clone A types')
ggsave('zscoresbytype_cloneA_5z.pdf',
       width = 6, height = 3,units = 'in')

#############Figure 3 C D E####
#identify most over convergent pair of glomeruli and the least convergent pair of glomeruli that has 1 glomeruli in common------
#z scores for all KC random bouton model
allz <- ana.list$all$zscore
#most over convergent pair that is not self self comparison
arrayInd(which(ana.list$all$zscore == allz[order(allz)][2600]),
         dim(ana.list$all$zscore))
rownames(ana.list$all$zscore)[28]
colnames(ana.list$all$zscore)[23]
#DP1l and VA2 z score is 
ana.list$all$zscore['DP1l','VA2']
#lowest z score that DP1l has with another glomeruli
min(ana.list$all$zscore['DP1l',])
#lowest z score that VA2 has with another glomeruli
min(ana.list$all$zscore['VA2',])
#we will go with VA2 partner
arrayInd(which(ana.list$all$zscore == min(ana.list$all$zscore['VA2',])),
         dim(ana.list$all$zscore))
colnames(ana.list$all$zscore)[3]
#it is DA2
#plot histograms of number of kcs cosampling the two pairs of glomeruli for each KC group-----------
#VA2 and DP1l histograms of shared KCs------
#plot histogram of expected numbers of shared KCs with actual for all KC random model

mod <- unlist(lapply(list.load('10k_allkcs_rbout_condinp.rds'), '[','DP1l','VA2'))
mod <- as.data.frame(mod)
ob <- read.csv('allkcsconinp.csv', row.names = 1)['DP1l','VA2']
#plot
ggplot(mod,aes(x=mod))+geom_histogram(binwidth = 1)+
  geom_vline(aes(xintercept=ob),color = 'red')+
  theme_classic()+ylim(0,2500)+xlim(-1,65)+
  xlab('allkcs DP1l VA2')
ggsave('cinput_DP1l_VA2_allkcs.pdf', units = 'in', width =3)
mean(mod$mod)
#plot for KC clone groups
for(c in c('A','B','C','D')){
  #plot histogram of expected number of shared KCs
  mod <- unlist(lapply(list.load(paste('10k_clone',c,'kcs_rbout_condinp.rds',sep = '')), '[','DP1l','VA2'))
  mod <- as.data.frame(mod)
  ob <- read.csv(paste('clone',c,'coninp.csv', sep = ''), row.names = 1)['DP1l','VA2']
  #plot
  ggplot(mod,aes(x=mod))+geom_histogram(binwidth = 1)+
    geom_vline(aes(xintercept=ob),color = 'red')+
    theme_classic()+ylim(0,2500)+xlim(-1,65)+
    xlab(paste('clone',c,'DP1l VA2'))
  ggsave(paste('cinput_DP1l_VA2_clone',c,'.pdf',sep = ''), units = 'in', width =3)
}
#plot for KC type groups
for(c in c('gamma','abprime','alphabeta')){
  #plot histogram of expected number of shared KCs
  mod <- unlist(lapply(list.load(paste('10k_',c,'kcs_rbout_condinp.rds',sep = '')), '[','DP1l','VA2'))
  mod <- as.data.frame(mod)
  ob <- read.csv(paste(c,'coninp.csv', sep = ''), row.names = 1)['DP1l','VA2']
  #plot
  ggplot(mod,aes(x=mod))+geom_histogram(binwidth = 1)+
    geom_vline(aes(xintercept=ob),color = 'red')+
    theme_classic()+ylim(0,2500)+xlim(-1,65)+
    xlab(paste(c,'DP1l VA2'))
  ggsave(paste('cinput_DP1l_VA2_',c,'.pdf',sep = ''), units = 'in', width =3)
}
#plot for KC clone and type groups
for(c in c('A','B','C','D')){
  for(t in c('gamma','abprime','alphabeta')){
    #plot histogram of expected number of shared KCs
    mod <- unlist(lapply(list.load(paste('10k_clone',c,t,'kcs_rbout_condinp.rds',sep = '')), '[','DP1l','VA2'))
    mod <- as.data.frame(mod)
    ob <- read.csv(paste('clone',c,t,'coninp.csv', sep = ''), row.names = 1)['DP1l','VA2']
    min(mod)
    max(mod)
    ob
    #plot
    ggplot(mod,aes(x=mod))+geom_histogram(binwidth = 1)+
      geom_vline(aes(xintercept=ob),color = 'red')+
      theme_classic()+ylim(0,2500)+xlim(-1,65)+
      xlab(paste('clone',c,t,'DP1l VA2'))
    ggsave(paste('cinput_DP1l_VA2_clone',c,t,'.pdf',sep = ''), units = 'in', width =3)
  }
}
#VA2 and DA2 histograms of shared KCs----------
#plot histogram of expected numbers of shared KCs with actual for all KC random model

mod <- unlist(lapply(list.load('10k_allkcs_rbout_condinp.rds'), '[','DA2','VA2'))
mod <- as.data.frame(mod)
ob <- read.csv('allkcsconinp.csv', row.names = 1)['DA2','VA2']
#plot
ggplot(mod,aes(x=mod))+geom_histogram(binwidth = 1)+
  geom_vline(aes(xintercept=ob),color = 'red')+
  theme_classic()+ylim(0,2500)+xlim(-1,65)+
  xlab('allkcs DA2 VA2')
ggsave('cinput_DA2_VA2_allkcs.pdf', units = 'in', width =3)
mean(mod$mod)
#plot for KC clone groups
for(c in c('A','B','C','D')){
  #plot histogram of expected number of shared KCs
  mod <- unlist(lapply(list.load(paste('10k_clone',c,'kcs_rbout_condinp.rds',sep = '')), '[','DA2','VA2'))
  mod <- as.data.frame(mod)
  ob <- read.csv(paste('clone',c,'coninp.csv', sep = ''), row.names = 1)['DA2','VA2']
  #plot
  ggplot(mod,aes(x=mod))+geom_histogram(binwidth = 1)+
    geom_vline(aes(xintercept=ob),color = 'red')+
    theme_classic()+ylim(0,2500)+xlim(-1,65)+
    xlab(paste('clone',c,'DA2 VA2'))
  ggsave(paste('cinput_DA2_VA2_clone',c,'.pdf',sep = ''), units = 'in', width =3)
}
#plot for KC type groups
for(c in c('gamma','abprime','alphabeta')){
  #plot histogram of expected number of shared KCs
  mod <- unlist(lapply(list.load(paste('10k_',c,'kcs_rbout_condinp.rds',sep = '')), '[','DA2','VA2'))
  mod <- as.data.frame(mod)
  ob <- read.csv(paste(c,'coninp.csv', sep = ''), row.names = 1)['DA2','VA2']
  #plot
  ggplot(mod,aes(x=mod))+geom_histogram(binwidth = 1)+
    geom_vline(aes(xintercept=ob),color = 'red')+
    theme_classic()+ylim(0,2500)+xlim(-1,65)+
    xlab(paste(c,'DA2 VA2'))
  ggsave(paste('cinput_DA2_VA2_',c,'.pdf',sep = ''), units = 'in', width =3)
}
#plot for KC clone and type groups
for(c in c('A','B','C','D')){
  for(t in c('gamma','abprime','alphabeta')){
    #plot histogram of expected number of shared KCs
    mod <- unlist(lapply(list.load(paste('10k_clone',c,t,'kcs_rbout_condinp.rds',sep = '')), '[','DA2','VA2'))
    mod <- as.data.frame(mod)
    ob <- read.csv(paste('clone',c,t,'coninp.csv', sep = ''), row.names = 1)['DA2','VA2']
    min(mod)
    max(mod)
    ob
    #plot
    ggplot(mod,aes(x=mod))+geom_histogram(binwidth = 1)+
      geom_vline(aes(xintercept=ob),color = 'red')+
      theme_classic()+ylim(0,2500)+xlim(-1,65)+
      xlab(paste('clone',c,t,'DA2 VA2'))
    ggsave(paste('cinput_DA2_VA2_clone',c,t,'.pdf',sep = ''), units = 'in', width =3)
  }
}
#overconvergent neuron skeletons VA2 and DPl1--------
#identify KCs that cosample VA2 and DP1l
commonkcs <- unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DP1l'])[
  unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])%in%
    unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DP1l'])]
#VA2
plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='VA2']], col = 'red',lwd = 2)
spheres3d(x= bout.info$center.x[bout.info$pn.glomerulus == 'VA2'],
          y= bout.info$center.y[bout.info$pn.glomerulus == 'VA2'],
          z= bout.info$center.z[bout.info$pn.glomerulus == 'VA2'],
          col = 'red',
          radius = 250,
          lit= F,
          alpha = 0.5)
shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'VA2'],
          y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'VA2'],
          z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'VA2'],
          col = 'brown4',
          radius = 100,
          lit= F,
          alpha = 1)
#DP1l
plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='DP1l']],col = 'black',lwd = 2)
spheres3d(x= bout.info$center.x[bout.info$pn.glomerulus == 'DP1l'],
          y= bout.info$center.y[bout.info$pn.glomerulus == 'DP1l'],
          z= bout.info$center.z[bout.info$pn.glomerulus == 'DP1l'],
          col = 'black',
          radius = 250,
          lit= F,
          alpha = 0.5)
shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'DP1l'],
          y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'DP1l'],
          z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'DP1l'],
          col = 'black',
          radius = 100,
          lit= F,
          alpha = 1)
#with calyx and axes
plot3d(ca, alpha = 0.1, add= T)
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.5)
rgl.snapshot('VA2_DP1l_skeletons_AP_ca_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.5)
rgl.snapshot('VA2_DP1l_skeletons_ML_ca_axes.png', fmt = 'png')
#clear before next KC group
clear3d()
#without calyx and axes
plot3d(ca, alpha = 0, add= T)
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.5)
rgl.snapshot('VA2_DP1l_skeletons_AP.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.5)
rgl.snapshot('VA2_DP1l_skeletons_ML.png', fmt = 'png')
#clear before next KC group
clear3d()
#underconvergent skeletons VA2 and DA2-------
#VA2
plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='VA2']], col = 'red',lwd = 2)
spheres3d(x= bout.info$center.x[bout.info$pn.glomerulus == 'VA2'],
          y= bout.info$center.y[bout.info$pn.glomerulus == 'VA2'],
          z= bout.info$center.z[bout.info$pn.glomerulus == 'VA2'],
          col = 'red',
          radius = 250,
          lit= F,
          alpha = 0.5)
shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'VA2'],
          y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'VA2'],
          z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'VA2'],
          col = 'brown4',
          radius = 100,
          lit= F,
          alpha = 1)
#DA2
commonkcs <- unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])[
  unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DA2'])%in%
    unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])]

plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='DA2']],col = 'black', lwd = 2)
spheres3d(x= bout.info$center.x[bout.info$pn.glomerulus == 'DA2'],
          y= bout.info$center.y[bout.info$pn.glomerulus == 'DA2'],
          z= bout.info$center.z[bout.info$pn.glomerulus == 'DA2'],
          col = 'black',
          radius = 250,
          lit= F,
          alpha = 0.5)
shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'DA2'],
          y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'DA2'],
          z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'DA2'],
          col = 'black',
          radius = 100,
          lit= F,
          alpha = 1)
#with calyx and axes
plot3d(ca, alpha = 0.1, add= T)
axes3d(edges = 'bbox',alpha=1)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.5)
rgl.snapshot('VA2_DA2_skeletons_AP_ca_axes.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.5)
rgl.snapshot('VA2_DA2_skeletons_ML_ca_axes.png', fmt = 'png')
#clear before next KC group
clear3d()
#without calyx and axes
plot3d(ca, alpha = 0, add= T)
axes3d(edges = 'bbox',alpha=0)
#anterior posterior view
view3d(userMatrix = apdv, zoom =.5)
rgl.snapshot('VA2_DA2_skeletons_AP.png', fmt = 'png')
#medial lateral view
view3d(userMatrix = mldv, zoom =.5)
rgl.snapshot('VA2_DA2_skeletons_ML.png', fmt = 'png')
#clear before next KC group
clear3d()

#plot skeletons of two pairs of glomeruli with boutons included in KC group random model shown- color co-sampling claws------
#types VA2 and DP1l
for(j in c('gamma','alphaprimebetaprime','alphabeta')){
  jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, kc.largetype == j)$pn.boutId),]
  
  commonkcs <- unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DP1l'])[
    unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])%in%
      unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DP1l'])]
  commonkcs <- commonkcs[commonkcs%in%claw.info$kc.bodyId[claw.info$kc.largetype==j]]
  #VA2
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='VA2']], col = 'red',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'VA2'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'VA2'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'VA2'],
            col = 'red',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'VA2'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'VA2'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'VA2'],
            col = 'brown4',
            radius = 100,
            lit= F,
            alpha = 1)
  #DP1l
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='DP1l']],col = 'black',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'DP1l'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'DP1l'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'DP1l'],
            col = 'black',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'DP1l'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'DP1l'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'DP1l'],
            col = 'black',
            radius = 100,
            lit= F,
            alpha = 1)
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,'_VA2_DP1l_skeleton_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,'_VA2_DP1l_skeleton_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
}
#types VA2 and DA2
for(j in c('gamma','alphaprimebetaprime','alphabeta')){
  jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, kc.largetype == j)$pn.boutId),]
  
  commonkcs <- unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DA2'])[
    unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])%in%
      unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DA2'])]
  commonkcs <- commonkcs[commonkcs%in%claw.info$kc.bodyId[claw.info$kc.largetype==j]]
  #VA2
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='VA2']], col = 'red',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'VA2'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'VA2'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'VA2'],
            col = 'red',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'VA2'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'VA2'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'VA2'],
            col = 'brown4',
            radius = 100,
            lit= F,
            alpha = 1)
  #DA2
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='DA2']],col = 'black',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'DA2'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'DA2'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'DA2'],
            col = 'black',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'DA2'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'DA2'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'DA2'],
            col = 'black',
            radius = 100,
            lit= F,
            alpha = 1)
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,'_VA2_DA2_skeleton_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,'_VA2_DA2_skeleton_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
}
#clones VA2 and DP1l
for(j in c('a','b','c','d')){
  jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, nb.cluster == j)$pn.boutId),]
  
  commonkcs <- unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DP1l'])[
    unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])%in%
      unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DP1l'])]
  commonkcs <- commonkcs[commonkcs%in%claw.info$kc.bodyId[claw.info$nb.cluster==j]]
  #VA2
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='VA2']], col = 'red',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'VA2'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'VA2'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'VA2'],
            col = 'red',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'VA2'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'VA2'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'VA2'],
            col = 'brown4',
            radius = 100,
            lit= F,
            alpha = 1)
  #DP1l
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='DP1l']],col = 'black',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'DP1l'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'DP1l'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'DP1l'],
            col = 'black',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'DP1l'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'DP1l'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'DP1l'],
            col = 'black',
            radius = 100,
            lit= F,
            alpha = 1)
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,'_VA2_DP1l_skeleton_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,'_VA2_DP1l_skeleton_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
}
#clones VA2 and DA2
for(j in c('a','b','c','d')){
  jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, nb.cluster == j)$pn.boutId),]
  
  commonkcs <- unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DA2'])[
    unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])%in%
      unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DA2'])]
  commonkcs <- commonkcs[commonkcs%in%claw.info$kc.bodyId[claw.info$nb.cluster==j]]
  #VA2
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='VA2']], col = 'red',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'VA2'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'VA2'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'VA2'],
            col = 'red',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'VA2'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'VA2'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'VA2'],
            col = 'brown4',
            radius = 100,
            lit= F,
            alpha = 1)
  #DA2
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='DA2']],col = 'black',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'DA2'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'DA2'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'DA2'],
            col = 'black',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'DA2'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'DA2'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'DA2'],
            col = 'black',
            radius = 100,
            lit= F,
            alpha = 1)
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,'_VA2_DA2_skeleton_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,'_VA2_DA2_skeleton_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
}
#clones and types VA2 and DP1l
for(j in c('a','b','c','d')){
  for(k in c('gamma','alphaprimebetaprime','alphabeta')){
  jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, nb.cluster == j, kc.largetype == k)$pn.boutId),]
  
  commonkcs <- unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DP1l'])[
    unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])%in%
      unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DP1l'])]
  commonkcs <- commonkcs[commonkcs%in%claw.info$kc.bodyId[claw.info$nb.cluster==j]]
  commonkcs <- commonkcs[commonkcs%in%claw.info$kc.bodyId[claw.info$kc.largetype==k]]
  #VA2
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='VA2']], col = 'red',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'VA2'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'VA2'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'VA2'],
            col = 'red',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'VA2'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'VA2'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'VA2'],
            col = 'brown4',
            radius = 100,
            lit= F,
            alpha = 1)
  #DP1l
  plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='DP1l']],col = 'black',lwd = 2)
  spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'DP1l'],
            y= jbout.info$center.y[jbout.info$pn.glomerulus == 'DP1l'],
            z= jbout.info$center.z[jbout.info$pn.glomerulus == 'DP1l'],
            col = 'black',
            radius = 250,
            lit= F,
            alpha = 0.5)
  shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
  spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'DP1l'],
            y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'DP1l'],
            z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'DP1l'],
            col = 'black',
            radius = 100,
            lit= F,
            alpha = 1)
  #without calyx and axes
  plot3d(ca, alpha = 0, add= T)
  axes3d(edges = 'bbox',alpha=0)
  #anterior posterior view
  view3d(userMatrix = apdv, zoom =.5)
  rgl.snapshot(paste(j,k,'_VA2_DP1l_skeleton_AP.png',sep=''), fmt = 'png')
  #medial lateral view
  view3d(userMatrix = mldv, zoom =.5)
  rgl.snapshot(paste(j,k,'_VA2_DP1l_skeleton_ML.png',sep=''), fmt = 'png')
  #clear before next KC group
  clear3d()
  }
}
#clones and types VA2 and DA2
for(j in c('a','b','c','d')){
  for(k in c('gamma','alphaprimebetaprime','alphabeta')){
    jbout.info <- bout.info[bout.info$pn.boutId%in%unique(filter(claw.info, nb.cluster == j, kc.largetype == k)$pn.boutId),]
    
    commonkcs <- unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DA2'])[
      unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='VA2'])%in%
        unique(claw.info$kc.bodyId[claw.info$pn.glomerulus=='DA2'])]
    commonkcs <- commonkcs[commonkcs%in%claw.info$kc.bodyId[claw.info$nb.cluster==j]]
    commonkcs <- commonkcs[commonkcs%in%claw.info$kc.bodyId[claw.info$kc.largetype==k]]
    #VA2
    plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='VA2']], col = 'red',lwd = 2)
    spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'VA2'],
              y= jbout.info$center.y[jbout.info$pn.glomerulus == 'VA2'],
              z= jbout.info$center.z[jbout.info$pn.glomerulus == 'VA2'],
              col = 'red',
              radius = 250,
              lit= F,
              alpha = 0.5)
    shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
    spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'VA2'],
              y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'VA2'],
              z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'VA2'],
              col = 'brown4',
              radius = 100,
              lit= F,
              alpha = 1)
    #DA2
    plot3d(pn.skel[names(pn.skel)%in%pn.info$pn.bodyId[pn.info$pn.glomerulus=='DA2']],col = 'black',lwd = 2)
    spheres3d(x= jbout.info$center.x[jbout.info$pn.glomerulus == 'DA2'],
              y= jbout.info$center.y[jbout.info$pn.glomerulus == 'DA2'],
              z= jbout.info$center.z[jbout.info$pn.glomerulus == 'DA2'],
              col = 'black',
              radius = 250,
              lit= F,
              alpha = 0.5)
    shareclaw.info <- claw.info[claw.info$kc.bodyId%in%commonkcs,]
    spheres3d(x= shareclaw.info$center.x[shareclaw.info$pn.glomerulus == 'DA2'],
              y= shareclaw.info$center.y[shareclaw.info$pn.glomerulus == 'DA2'],
              z= shareclaw.info$center.z[shareclaw.info$pn.glomerulus == 'DA2'],
              col = 'black',
              radius = 100,
              lit= F,
              alpha = 1)
    #without calyx and axes
    plot3d(ca, alpha = 0, add= T)
    axes3d(edges = 'bbox',alpha=0)
    #anterior posterior view
    view3d(userMatrix = apdv, zoom =.5)
    rgl.snapshot(paste(j,k,'_VA2_DA2_skeleton_AP.png',sep=''), fmt = 'png')
    #medial lateral view
    view3d(userMatrix = mldv, zoom =.5)
    rgl.snapshot(paste(j,k,'_VA2_DA2_skeleton_ML.png',sep=''), fmt = 'png')
    #clear before next KC group
    clear3d()
  }
}

