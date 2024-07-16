#within and between type comparisons of bouton locations
#originally written 1/24/2024
#load packages------
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
bout.info <- read.csv('bouton_info_KCconnectivity.csv')[,1:12]
#list of PNs with boutons ordered Medial to Lateral that the next for loop can pull from----------
pnbouts <- list()
for(i in 1:nrow(pn.info)){
  pn <- pn.info$pn.bodyId[i]
  boutons <- bout.info$pn.boutId[bout.info$pn.bodyId == pn][
    order(bout.info$center.x[bout.info$pn.bodyId == pn])]
  pnbouts[[i]] <- boutons
  names(pnbouts)[i] <- pn
}
#pairwise distances between boutons on different types of PNs with the same number of boutons------

pwbout <- vector(length = 0)
for(n in c(1:9)){
  #PNs with the given number of boutons
  pns <- filter(pn.info, nbout ==n)
  #boutons of pns
  nbouts <- pnbouts[names(pnbouts)%in%pns$pn.bodyId]
  #glomeruli represented in the group
  gloms <- unique(pns$pn.glomerulus)
  #iterate through glomeruli, comparing each PN to PNs of different glomeruli avoiding redundant comparisons
  for(a in 1:(length(gloms)-1)){
    #pns in glom a
    glomapn <- pns$pn.bodyId[pns$pn.glomerulus%in%gloms[a]]
    #pns not in glom a that haven't been compared yet
    otherpn <- pns$pn.bodyId[pns$pn.glomerulus%in%gloms[(a+1):length(gloms)]]
    #compute pairwise distances between each PN belonging to glomA and otherpns
    if(length(nbouts) >1){  
      for(i in 1:length(glomapn)){
        for(j in 1:length(otherpn)){
          xyboutdist <- vector(length = 0)
          for(k in 1:n){
            boutx <- filter(bout.info, pn.boutId == nbouts[names(nbouts)%in%glomapn[i]][[1]][k])
            bouty <- filter(bout.info, pn.boutId == nbouts[names(nbouts)%in%otherpn[j]][[1]][[k]])
            xyboutdist <- c(xyboutdist, 
                            sqrt(((boutx$center.x-bouty$center.x)^2)+
                                   ((boutx$center.y-bouty$center.y)^2)+
                                   ((boutx$center.z-bouty$center.z)^2))
            )
          }
          pwbout <- c(pwbout,xyboutdist)
        }
      }
    }
  }
  
}

#filter out pn types with only 1 neuron----------
for(i in 1:length(unique(pn.info$pn.glomerulus))){
  pn.info$nneuron[pn.info$pn.glomerulus==unique(pn.info$pn.glomerulus)[i]] <- 
    nrow(filter(pn.info, pn.glomerulus == unique(pn.info$pn.glomerulus)[i]))
}
pn.info <- pn.info[pn.info$nneuron!=1,]
bout.info <- bout.info[bout.info$pn.bodyId%in%pn.info$pn.bodyId,]

#make glomeruli dataframe to store info about PNs of the same type----------
glomeruli <-unique(pn.info[,c(3:6,10)])
#list of pairwise distances between cognate boutons on PNs of the same type---------
glomdeltabout <- list()
for(g in 1:nrow(glomeruli)){
  #defines glomerulus of interest
  glom <- glomeruli$pn.glomerulus[g]
  #defines PNs that innervate the glomerulus of interest
  glom.pns <- pn.info$pn.bodyId[pn.info$pn.glomerulus==glom]
  #list of bouton ids for each PN ordered by ML position
  glom.bouts <- pnbouts[names(pnbouts)%in%glom.pns]
  #determine what the most common number of boutons is for PNs of this type
  glomeruli$nbout[g] <- getmode(lengths(glom.bouts))
  #limit list of PNs to those with the most common number of boutons
  glom.bouts <- glom.bouts[lengths(glom.bouts)%in%glomeruli$nbout[g]]
  #record how many neurons have the most common number of boutons
  glomeruli$neuronsbouts[g] <- length(glom.bouts)
  #compute pairwise distances between cognate branch points
  pwboutdist <- vector(length = 0)
  if(length(glom.bouts) >1){
    for(p in 1:(length(glom.bouts)-1)){
      for(q in (p+1):length(glom.bouts)){
        pqboutdist <- vector(length = 0)
        for(i in 1:glomeruli$nbout[g]){
          boutx <- filter(bout.info, pn.boutId == glom.bouts[[p]][[i]])
          bouty <- filter(bout.info, pn.boutId == glom.bouts[[q]][[i]])
          pqboutdist <- c(pqboutdist, 
                          sqrt(((boutx$center.x-bouty$center.x)^2)+
                                 ((boutx$center.y-bouty$center.y)^2)+
                                 ((boutx$center.z-bouty$center.z)^2))
          )
        }
        pwboutdist <- c(pwboutdist,pqboutdist)
      }
    }
  }
  if(length(glom.bouts)==1){
    pwboutdist <- c(pwboutdist,NA)
  }
  #store pairwise distances between cognate branch points
  glomdeltabout[[g]] <- pwboutdist
  #name list entry after the PN type/ glomerulus
  names(glomdeltabout)[g] <- glomeruli$pn.glomerulus[g]
}


#plot pairwise distances between cognate PNs and other PNs----
#get vector of pairwise distances between cognate branch points on sister PNs
pwsisterbouts <- vector(length = 0)
for(i in 1:length(glomdeltabout)){
  pwsisterbouts <- c(pwsisterbouts, glomdeltabout[[i]])
}
pwsisterbouts <- pwsisterbouts[!is.na(pwsisterbouts)]
#select 98 random stranger distances to plot
pwbout.rand <- unique(pwbout[round(runif(98,1,2263))])
#make data into a plottable data frame
plotbouts <- as.data.frame(rbind(cbind(pwsisterbouts, 'sister'),cbind(pwbout.rand, 'stranger')))
colnames(plotbouts) <- c('distance', 'type')
#shapiro test for normality
shapiro.test(pwsisterbouts)#data not normal
shapiro.test(pwbout.rand)#data is normal but will do rank sum test still
#wilcox rank sum test 
test <- wilcox.test(pwsisterbouts, pwbout.rand)
m1 <- round(median(pwsisterbouts),2)
m2 <- round(median(pwbout.rand),2)
#cohensD for effect size
cd <- cohensd(pwbout.rand,pwsisterbouts)
#plot data
ggplot(plotbouts, aes(x= as.factor(type), y = as.numeric(distance)))+geom_beeswarm()+geom_boxplot()+
  ggtitle(paste("wilcox test p =",round(test$p.value,15),'medians',m1, ' ', m2))+
  theme(plot.title = element_text(size = 5))+
  ylab('distance between cognate boutons')+
  ylim(0,4600)
ggsave('cognatebouton_distances.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')
