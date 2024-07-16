#branch point location analyses
#originally written 1/31/24
#edited 3/29/2024
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
#pw distances between different types of PNs with the same number of branch points------
#for comparisons sake compute pairwise distances between each PN and a PN
#with the same number of branch points but different glomerulus
#for each neuron compute the pw distances btw the ML position of branch points with all PNs that aren't of the same type
#list of PNs with the X coordinates of their branch points
pwdist <- vector(length = 0)
for(n in 1:6){
  #PNs with the given number of branch points
  pns <- filter(pn.info, nbp ==n)
  #branch points of pns
  pnbps <- pnbranchpts[names(pnbranchpts)%in%pns$pn.bodyId]
  #glomeruli represented in the group
  gloms <- unique(pns$pn.glomerulus)
  #iterate through glomeruli, comparing each PN to PNs of different glomeruli avoiding redundant comparisons
  for(a in 1:(length(gloms)-1)){
    #pns in glom a
    glomapn <- pns$pn.bodyId[pns$pn.glomerulus%in%gloms[a]]
    #pns not in glom a that haven't been compared yet
    otherpn <- pns$pn.bodyId[pns$pn.glomerulus%in%gloms[(a+1):length(gloms)]]
    #compute pairwise distances between each PN belonging to glomA and otherpns
    for(i in 1:length(glomapn)){
      for(j in 1:length(otherpn)){
        pwdist <- c(pwdist,
                    abs(pnbps[names(pnbps)%in%glomapn[i]][[1]]-pnbps[names(pnbps)%in%otherpn[j]][[1]])
        )
      }
    }
    
  }
  
}
#filter out pn types with only 1 neuron--------------
for(i in 1:length(unique(pn.info$pn.glomerulus))){
  pn.info$nneuron[pn.info$pn.glomerulus==unique(pn.info$pn.glomerulus)[i]] <- 
    nrow(filter(pn.info, pn.glomerulus == unique(pn.info$pn.glomerulus)[i]))
}
pn.info <- pn.info[pn.info$nneuron!=1,]
bout.info <- bout.info[bout.info$pn.bodyId%in%pn.info$pn.bodyId,]
#make glomeruli dataframe to store comparisons in
glomeruli <-unique(pn.info[,c(3,6)])
#Pairwise distances between cognate BPs on PNs------
#for each glomerulus calculate pairwise differences in ml position of cognate branch points
#list that will hold the name of each glomerulus/ PN type and the vector of pw differences in ML location of 
#cognate branch points
glomdeltabp <- list()
for(g in 1:nrow(glomeruli)){
  #defines glomerulus of interest
  glom <- glomeruli$pn.glomerulus[g]
  #defines PNs that innervate the glomerulus of interest
  glom.pns <- pn.info$pn.bodyId[pn.info$pn.glomerulus==glom]
  #list that will hold x coordinates of branchpoints on PNs of this type
  pn.branchpt <- pnbranchpts[names(pnbranchpts)%in%glom.pns]
  #determine what the most common number of branch points is for PNs of this type
  glomeruli$nbp[g] <- getmode(lengths(pn.branchpt))
  #limit list of branch points to those with the most common number of branch points
  pn.branchpt <- pn.branchpt[lengths(pn.branchpt)%in%glomeruli$nbp[g]]
  #record how many neurons have the most common number of branch points
  glomeruli$neuronsbp[g] <- length(pn.branchpt)
  #compute pairwise distances between cognate branch points
  pwbpdist <- vector(length = 0)
  if(length(pn.branchpt) >1){
    for(p in 1:(length(pn.branchpt)-1)){
      for(q in (p+1):length(pn.branchpt)){
        pwbpdist <- c(pwbpdist,abs(pn.branchpt[[p]]-pn.branchpt[[q]]))
      }
    }
  }
  if(length(pn.branchpt)==1){
    pwbpdist <- c(pwbpdist,NA)
  }
  #store pairwise distances between cognate branch points
  glomdeltabp[[g]] <- pwbpdist
  #name list entry after the PN type/ glomerulus
  names(glomdeltabp)[g] <- glomeruli$pn.glomerulus[g]
}
#plot pairwise distances between cognate PNs and other PNs----
#get vector of pairwise distances between cognate branch points on sister PNs
pwsisterdist <- vector(length = 0)
for(i in 1:length(glomdeltabp)){
  pwsisterdist <- c(pwsisterdist, glomdeltabp[[i]])
}
pwsisterdist <- pwsisterdist[!is.na(pwsisterdist)]
#select 88 random stranger distances to plot
pwdist.rand <-unique(pwdist[round(runif(88,1,2818))])
#make data into a plottable data frame
plotdist <- as.data.frame(rbind(cbind(pwsisterdist, 'sister'),cbind(pwdist.rand, 'stranger')))
colnames(plotdist) <- c('distance', 'type')
#shapiro test for normality
shapiro.test(pwsisterdist)
shapiro.test(pwdist.rand)#data is nor normally distributed will do a wilcox rank sum test
#wilcox rank sum test 
test <- wilcox.test(pwsisterdist, pwdist.rand)
m1 <- round(median(pwsisterdist),2)
m2 <- round(median(pwdist.rand),2)
#cohensD for effect size
cd <- cohensd(pwsisterdist,pwdist.rand)
#plot data
ggplot(plotdist, aes(x= type, y = as.numeric(distance)))+geom_beeswarm()+
  geom_boxplot()+
  theme_classic()+
  ggtitle(paste("wilcox test p =",round(test$p.value,6),'medians',m1, ' ', m2,
                ' es:', round(cd, 4)))+
  theme(plot.title = element_text(size = 5))+
  ylab('cognate bp ml difference')+
  ylim(0,4600)
ggsave('cognatebp_ml_differences.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')




