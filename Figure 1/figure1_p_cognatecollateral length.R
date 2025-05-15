#collateral length analyses
#load packages------
library(dplyr)
library(ggbeeswarm)
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
theme_set(theme_classic())
#set working directory and load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
pn.info <- read.csv('PN_info.csv')
bout.info <- read.csv('bouton_info_collateral_bps.csv')
#calculate bp to bouton distance (collateral length proxy)------------
bout.info$dx <- abs(bout.info$center.x- bout.info$bp.x)
bout.info$dy <- abs(bout.info$center.y- bout.info$bp.y)
bout.info$dz <- abs(bout.info$center.z- bout.info$bp.z)
bout.info$dbp <- sqrt((bout.info$dx^2)+(bout.info$dy^2)+(bout.info$dz^2))
#get ML position of branch points for each neuron ordered ML
pnbranchpts <- list()
for(i in 1:nrow(pn.info)){
  pn <- pn.info$pn.bodyId[i]
  bps <- unique(bout.info$bp.x[bout.info$pn.bodyId== pn])
  bps <- bps[order(bps)]
  pnbranchpts[[i]] <- bps
  names(pnbranchpts)[i] <- pn
  pn.info$nbp[i] <- length(bps)
}
#get collateral lengths ordered ML for each PN
pncollaterals <- list()
for(i in 1:nrow(pn.info)){
  pn <- pn.info$pn.bodyId[i]
  collateral <- bout.info$coll.length[bout.info$pn.bodyId== pn]
  bps <- bout.info$bp.x[bout.info$pn.bodyId== pn]
  collateral <- collateral[order(bps)]
  pncollaterals[[i]] <- collateral
  names(pncollaterals)[i] <- pn
}

#pw differences in collateral length between different types of PNs with the same number of collaterals------
#get number of collaterals (really bouton number) for each PN
for(i in 1:nrow(pn.info)){
  pn.info$ncoll[i] <- length(pncollaterals[names(pncollaterals) == pn.info$pn.bodyId[i]][[1]])
}

pwdiff <- vector(length = 0)
for(n in 1:7){
  #PNs with the given number of branch points
  pns <- filter(pn.info, ncoll ==n)
  #branch points of pns
  pn.coll <- pncollaterals[names(pncollaterals)%in%pns$pn.bodyId]
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
        pwdiff <- c(pwdiff,
                    abs(pn.coll[names(pn.coll)%in%glomapn[i]][[1]]-pn.coll[names(pn.coll)%in%otherpn[j]][[1]])
        )
      }
    }
    
  }
  
}
#compare cognate collateral lengths--------------------
#make glomeruli dataframe to store comparisons in
glomeruli <-data.frame(pn.glomerulus = unique(pn.info[,3]))
#filter out pn types with only 1 neuron
for(i in 1:length(unique(pn.info$pn.glomerulus))){
  pn.info$nneuron[pn.info$pn.glomerulus==unique(pn.info$pn.glomerulus)[i]] <- 
    nrow(filter(pn.info, pn.glomerulus == unique(pn.info$pn.glomerulus)[i]))
}
pn.info <- pn.info[pn.info$nneuron!=1,]
bout.info <- bout.info[bout.info$pn.bodyId%in%pn.info$pn.bodyId,]

#Pairwise differences in collateral length------
#for each glomerulus calculate pairwise differences in collateral length of cognate boutons
#list that will hold the name of each glomerulus/ PN type and the vector of collateral lengths ordered ml
glomdeltacol <- list()
for(g in 1:nrow(glomeruli)){
  #defines glomerulus of interest
  glom <- glomeruli$pn.glomerulus[g]
  #defines PNs that innervate the glomerulus of interest
  glom.pns <- pn.info$pn.bodyId[pn.info$pn.glomerulus==glom]
  #list that will hold collateral lengths on PNs of this type
  pn.coll <- pncollaterals[names(pncollaterals)%in%glom.pns]
  #determine what the most common number of boutons is for PNs of this type
  glomeruli$ncoll[g] <- getmode(lengths(pn.coll))
  #limit list of branch points to those with the most common number of branch points
  pn.coll <- pn.coll[lengths(pn.coll)%in%glomeruli$ncoll[g]]
  #record how many neurons have the most common number of branch points
  glomeruli$neuronscoll[g] <- length(pn.coll)
  #compute pairwise distances between cognate branch points
  pwdbp <- vector(length = 0)
  if(length(pn.coll) >1){
    for(p in 1:(length(pn.coll)-1)){
      for(q in (p+1):length(pn.coll)){
        pwdbp <- c(pwdbp,abs(pn.coll[[p]]-pn.coll[[q]]))
      }
    }
  }
  if(length(pn.coll)==1){
    pwdbp <- c(pwdbp,NA)
  }
  #store pairwise distances between cognate branch points
  glomdeltacol[[g]] <- pwdbp
  #name list entry after the PN type/ glomerulus
  names(glomdeltacol)[g] <- glomeruli$pn.glomerulus[g]
}

#plot pairwise distances between cognate PNs and other PNs----
#get vector of pairwise distances between cognate branch points on sister PNs
pwsistercoll <- vector(length = 0)
for(i in 1:length(glomdeltacol)){
  pwsistercoll <- c(pwsistercoll, glomdeltacol[[i]])
}
pwsistercoll <- pwsistercoll[!is.na(pwsistercoll)]
#select 98 random stranger distances to plot
pwcoll.rand <- unique(pwdiff[round(runif(105,1,2206))])
pwcoll.rand <- pwcoll.rand[1:98]
#make data into a plottable data frame
plotcoll <- as.data.frame(rbind(cbind(pwsistercoll, 'sister'),cbind(pwcoll.rand, 'stranger')))
colnames(plotcoll) <- c('distance', 'type')
#shapiro test for normality
shapiro.test(pwsistercoll)
shapiro.test(pwcoll.rand)#data is nor normally distributed will do a wilcox rank sum test
#wilcox rank sum test 
test <- wilcox.test(pwsistercoll, pwcoll.rand)
m1 <- round(median(pwsistercoll),2)
m2 <- round(median(pwcoll.rand),2)
#cohensD for effect size
cd <- cohensd(pwsistercoll,pwcoll.rand)
#plot data
ggplot(plotcoll, aes(x= type, y = as.numeric(distance)))+geom_beeswarm()+geom_boxplot()+
  ggtitle(paste("wilcox test p =",round(test$p.value,10),'medians',m1, ' ', m2, ' es:',round(cd,4)))+
  theme(plot.title = element_text(size = 5))+
  ylab('difference in cognate collateral length')+
  ylim(0,4600)
ggsave('cognatecollateral_differences.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

