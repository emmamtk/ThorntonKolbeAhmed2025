#Claw location by clone and type
#originally written 9/21/2023
#edited for clarity 11/14/2023
#loading packages------
library(dplyr)
library(ggplot2)
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
pn.info <- read.csv('PN_info.csv')
bout.info <- read.csv('bouton_info_KCconnectivity.csv')
kc.info <- read.csv('KC_info.csv')
claw.info <- read.csv('claw_info.csv')
#make dataframe that is glomeruli as rows and connectivity to each KC type as columns
#get list of unique glomeruli
glomeruli <- data.frame(pn.glomerulus=unique(pn.info$pn.glomerulus))
#add number of boutons per glomerulus and order glomeruli by this
glomeruli$nbout=0
for (i in 1:nrow(glomeruli)){
  glomeruli$nbout[i] <- nrow(unique(filter(bout.info, pn.glomerulus == glomeruli$pn.glomerulus[i])))
}
glomeruli$pn.glomerulus <- factor(glomeruli$pn.glomerulus, 
                                  levels=glomeruli$pn.glomerulus[order(glomeruli$nbout)] )
#get number of claws of each KC type attached to boutons of this glomerulus
for(i in 1:nrow(glomeruli)){
  gclaw <- filter(claw.info, pn.glomerulus == glomeruli$pn.glomerulus[i])
  glomeruli$totalclaw[i] <- nrow(gclaw)
  glomeruli$gamma.claw[i] <- nrow(filter(gclaw, kc.largetype == 'gamma'))
  glomeruli$abprime.claw[i]<- nrow(filter(gclaw, kc.largetype == 'alphaprimebetaprime'))
  glomeruli$ab.claw[i]<- nrow(filter(gclaw, kc.largetype == 'alphabeta'))
  glomeruli$clonea.claw[i]<- nrow(filter(gclaw, nb.cluster == 'a'))
  glomeruli$cloneb.claw[i]<- nrow(filter(gclaw, nb.cluster == 'b'))
  glomeruli$clonec.claw[i]<- nrow(filter(gclaw, nb.cluster == 'c'))
  glomeruli$cloned.claw[i]<- nrow(filter(gclaw, nb.cluster == 'd'))
}

#data frame: claw to glomeruli proportion by type--------
#this data frame will have the proportion of claws of each kc type receiving input from each of the glomeruli
#basically: nclaw of kc type receiving input from glomeruli a/ total claws from kc type
#if there are no claws of that KC type receiving input it will assign the proportion as 0 to avoid dividing by 0
c2g.prop <- glomeruli
for(i in 3:10){
  if(sum(c2g.prop[,i]) == 0) c2g.prop[,i] <- 0
  if(sum(c2g.prop[,i]) > 0) c2g.prop[,i] <- c2g.prop[,i]/sum(c2g.prop[,i]) 
}


#data frame: convert above to a form that I can make a dot plot from--------
#start data frame with info from all kcs
c2g.prop.dot <- cbind(c2g.prop[,1:2], #info about each glomerulus
                      c2g.prop[,3], #proportion of claws from all KCs that receive inputs from each glomerulus
                      c2g.prop[,3]/c2g.prop[,3], #fold change in proportion of claws that receive inputs from each glomerulus from proportion for all kcs
                      colnames(c2g.prop)[3]) #column name of column 3 which is the group of Kenyon cells that column represents
colnames(c2g.prop.dot)[3:5] <- c('claw.glom.sample','fold.prop.change','kctype')
#add data only from other kc groups
for(i in 4:10){
  x <- cbind(c2g.prop[,1:2],#info about each glomerulus
             c2g.prop[,i], #proportion of claws from i group of Kenyon cells that receive inputs from each glomerulus
             c2g.prop[,i]/c2g.prop[,3], #fold change in proportion of claws from i kcs that receive inputs from each glomerulus from proportion for all kcs
             colnames(c2g.prop)[i]) #column name of i which is the group of Kenyon cells that column represents
  colnames(x)[3:5] <- c('claw.glom.sample','fold.prop.change','kctype')
  c2g.prop.dot <- rbind(c2g.prop.dot,x)
}

#adding a column to assign bins to each fold change value so that it can be colored correctly in the dot plot 
#bins are .1 apart between -.1 and 2.5 there is another bin for values between 2.5 and 10
c2g.prop.dot<- c2g.prop.dot%>% mutate(new_bin = cut(fold.prop.change, breaks = c(seq(-.1,2.5,.1),10)))
#making a palette for the binned data
#seashell color represents a fold change of 1 (no change) there are 11 values from .1 to 1
pal <- colorRampPalette(c('skyblue3','seashell'))
pal(11)
#seashell color represents a fold change of 1 (no change) there are 16 values from 1 to 2.5
pal2 <- colorRampPalette(c('seashell','firebrick3'))
pal2(16)
#fold change of 0 is represented by skyblue4 and foldchange of anything above 2.5 is represented by firebrick red
#this allows more of a dynamic color range where most of the data is
scalecolor <- c('skyblue4', "#6CA6CD", "#7AADD0", "#89B5D3", "#98BDD6","#A6C5DA", "#B5CDDD", "#C4D5E0", "#D2DDE4", "#E1E5E7", "#F0EDEA", 
  "#FFF5EE", "#FBE7E0", "#F8D9D3", "#F5CBC6", "#F1BDB8", "#EEB0AB", "#EBA29E", "#E79490", "#E48683", "#E17876" ,
  "#DD6B68", "#DA5D5B", "#D74F4E", "#D34140", "#D03333", "firebrick4")


#add ordering to kc types so that the dot plot columns will be in the order I want them to be in
c2g.prop.dot$kctype <- factor(c2g.prop.dot$kctype, levels = c('totalclaw','gamma.claw','abprime.claw','ab.claw',
                                                              'clonea.claw','cloneb.claw','clonec.claw','cloned.claw'))
#make dot plot where dot position is determined by KC type (column) and glomeruli (row)
#dot size is the proportion of claws of that kc type getting input from boutons of that glomerulus
#dot color is fold change in that proportion from that of all kcs
ggplot(c2g.prop.dot,
       aes(x=kctype, y = pn.glomerulus, size = claw.glom.sample,color = new_bin))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = scalecolor)
ggsave('claw on glomerular proportion by KC type .pdf', 
       device = 'pdf', width =12 ,height=8, units='in')
