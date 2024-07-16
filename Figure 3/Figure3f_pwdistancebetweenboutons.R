#pairwise distances between boutons
#originally written 07/02/2024
#based on code from 02/26/2024
#load packages-------
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
cohensd <- function(x,y){
  #absolute value of difference in group means
  diffinmeans <- abs(mean(x)-mean(y))
  #pooled standard deviation- root mean square of the standard deviations
  pooledsd <- sqrt(((sd(x)^2)+(sd(y)^2))/2)
  #cohens d
  return(diffinmeans/pooledsd)
}
#load files------
setwd("C:/Users/emmamtk/Desktop/CSVs for paper")
bout.info <- read.csv('bouton_info_KCconnectivity.csv')

#pw distances between boutons--------
#make data frame with pairwise distances between each bouton and the number of KCs they have in common
#get list of every unique pair of boutons
pwbout <- matrix(ncol = 2, nrow=0)
for(i in 1:435){
  pwbout <- rbind(pwbout, 
                  cbind(bout.info$pn.boutId[i],bout.info$pn.boutId[(i+1):nrow(bout.info)]))
}
#make pw bout a dataframe with col names
pwbout <- data.frame(pwbout)
colnames(pwbout) <- c('bout1','bout2')

#populate pairwise distances between boutons and number of shared KCs
for(i in 1:nrow(pwbout)){
  bout.1 <- filter(bout.info, pn.boutId == pwbout$bout1[i])
  bout.2 <- filter(bout.info, pn.boutId == pwbout$bout2[i])
  pwbout$x.dist[i] <- bout.1$center.x-bout.2$center.x
  pwbout$y.dist[i] <- bout.1$center.y-bout.2$center.y
  pwbout$z.dist[i] <- bout.1$center.z-bout.2$center.z
  pwbout$dist[i] <- sqrt((pwbout$x.dist[i]^2) +
                           (pwbout$y.dist[i]^2)+
                           (pwbout$z.dist[i]^2))
  pwbout$ncommonkc[i] <- sum(((bout.1[,13:1740]>1)+(bout.2[,13:1740]>1))==2)
}
write.csv(pwbout, 'pw_bouton_distances.csv', row.names=FALSE)
pwbout <- read.csv('pw_bouton_distances.csv')

#add column for plotting called kc group that is 0 or any number of shared KCs
pwbout$kcgroup <- 'any number of shared KCs'
pwbout$kcgroup[pwbout$ncommonkc==0]<- '0'

#what percent of bouton pairs share no KCs
nrow(filter(pwbout, kcgroup == 0))/nrow(pwbout)
#what percent of bouton pairs share 1 or more KCs
nrow(filter(pwbout, kcgroup == 'any number of shared KCs'))/nrow(pwbout)
#plot pw distances with the same number of points in each group---------
#number of pairs that share 1 or more kcs
nrow(filter(pwbout, kcgroup == 'any number of shared KCs'))
#number of pairs that share 0 kcs
nrow(filter(pwbout, kcgroup == 0))
#generate 1000 random numbers between 1 and 75,331
randn <- unique(round(runif(n=1000,min=1,max=75331)))
randn <- unique(c(randn,unique(round(runif(n=10,min=1,max=75331)))))
randn <- randn[1:1000]
#generate 1000 random numbers between 1 and 19,499
randn2 <- unique(round(runif(n=1000,min=1,max=19499)))
randn2 <- unique(c(randn2,unique(round(runif(n=10,min=1,max=19499)))))
randn2 <- randn2[1:1000]
#select 1000 pairs of boutons that share no KCs
plot.pwbout <- filter(pwbout, kcgroup=='0')[randn,]
#select 1000 pairs of boutons that share 1 or more kcs to the dataframe
plot.pwbout <- rbind(plot.pwbout, filter(pwbout, kcgroup=='any number of shared KCs')[randn2,])
write.csv(plot.pwbout, 'pairwise_bouton_distances_1ksample.csv')
plot.pwbout <- read.csv('pairwise_bouton_distances_1ksample.csv')
#plot for total distance------
#shapiro test for normality
shapiro.test(filter(plot.pwbout, kcgroup == '0')$dist)
shapiro.test(filter(plot.pwbout, kcgroup == 'any number of shared KCs')$dist)
#data is not normal will use wilcox rank sum test
#wilcox rank sum test 
test <- wilcox.test(filter(plot.pwbout, kcgroup ==0)$dist,
                    filter(plot.pwbout, kcgroup =='any number of shared KCs')$dist,)
m1 <- round(median(filter(plot.pwbout, kcgroup ==0)$dist),2)
m2 <- round(median(filter(plot.pwbout, kcgroup =='any number of shared KCs')$dist),2)
#effect size as measured by cohen's D
cd <- cohensd(filter(plot.pwbout, kcgroup==0)$dist,
              filter(plot.pwbout, kcgroup=='any number of shared KCs')$dist)
#plot data
ggplot(plot.pwbout, aes(x= kcgroup, y = as.numeric(dist)))+geom_beeswarm()+geom_boxplot()+
  ggtitle(paste("wilcox test p =",test$p.value,'medians',m1, ' ', m2, 'es:',cd))+
  theme(plot.title = element_text(size = 5))+
  ylab('pairwise bouton dists')+
  ylim(0,8000)
ggsave('pwboutondists_1000sample.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

#plot for x dist ML-----
#convert x.dist to absolute value of x.dist
plot.pwbout$x.dist <- abs(plot.pwbout$x.dist)
#shapiro test for normality
shapiro.test(filter(plot.pwbout, kcgroup == 0)$x.dist)
shapiro.test(filter(plot.pwbout, kcgroup == 'any number of shared KCs')$x.dist)
#data is not normal will use wilcox rank sum test
#wilcox rank sum test 
test <- wilcox.test(filter(plot.pwbout, kcgroup ==0)$x.dist,
                    filter(plot.pwbout, kcgroup =='any number of shared KCs')$x.dist,)
m1 <- round(median(filter(plot.pwbout, kcgroup ==0)$x.dist),2)
m2 <- round(median(filter(plot.pwbout, kcgroup =='any number of shared KCs')$x.dist),2)
#effect size as measured by cohen's D
cd <- cohensd(filter(plot.pwbout, kcgroup==0)$x.dist,
              filter(plot.pwbout, kcgroup=='any number of shared KCs')$x.dist)
#plot data
ggplot(plot.pwbout, aes(x= kcgroup, y = as.numeric(x.dist)))+geom_beeswarm()+geom_boxplot()+
  ggtitle(paste("wilcox test p =",test$p.value,'medians',m1, ' ', m2, 'es:',cd))+
  theme(plot.title = element_text(size = 5))+
  ylab('pairwise bouton x.dists')+
  ylim(0,8000)
ggsave('pwboutonxdists_1000sample.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')
#plot for y dist AP------
#convert y.dist to absolute value of y.dist
plot.pwbout$y.dist <- abs(plot.pwbout$y.dist)
#shapiro test for normality
shapiro.test(filter(plot.pwbout, kcgroup == 0)$y.dist)
shapiro.test(filter(plot.pwbout, kcgroup == 'any number of shared KCs')$y.dist)
#data is not normal will use wilcox rank sum test
#wilcox rank sum test 
test <- wilcox.test(filter(plot.pwbout, kcgroup ==0)$y.dist,
                    filter(plot.pwbout, kcgroup =='any number of shared KCs')$y.dist,)
m1 <- round(median(filter(plot.pwbout, kcgroup ==0)$y.dist),2)
m2 <- round(median(filter(plot.pwbout, kcgroup =='any number of shared KCs')$y.dist),2)
#effect size as measured by cohen's D
cd <- cohensd(filter(plot.pwbout, kcgroup==0)$y.dist,
              filter(plot.pwbout, kcgroup=='any number of shared KCs')$y.dist)
#plot data
ggplot(plot.pwbout, aes(x= kcgroup, y = as.numeric(y.dist)))+geom_beeswarm()+geom_boxplot()+
  ggtitle(paste("wilcox test p =",test$p.value,'medians',m1, ' ', m2, 'es:',cd))+
  theme(plot.title = element_text(size = 5))+
  ylab('pairwise bouton y.dists')+
  ylim(0,8000)
ggsave('pwboutonydists_1000sample.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')
#plot for z dist DV------
#convert z.dist to absolute value of z.dist
plot.pwbout$z.dist <- abs(plot.pwbout$z.dist)
#shapiro test for normality
shapiro.test(filter(plot.pwbout, kcgroup == 0)$z.dist)
shapiro.test(filter(plot.pwbout, kcgroup == 'any number of shared KCs')$z.dist)
#data is not normal will use wilcox rank sum test
#wilcox rank sum test 
test <- wilcox.test(filter(plot.pwbout, kcgroup ==0)$z.dist,
                    filter(plot.pwbout, kcgroup =='any number of shared KCs')$z.dist,)
m1 <- round(median(filter(plot.pwbout, kcgroup ==0)$z.dist),2)
m2 <- round(median(filter(plot.pwbout, kcgroup =='any number of shared KCs')$z.dist),2)
#effect size as measured by cohen's D
cd <- cohensd(filter(plot.pwbout, kcgroup==0)$z.dist,
              filter(plot.pwbout, kcgroup=='any number of shared KCs')$z.dist)
#plot data
ggplot(plot.pwbout, aes(x= kcgroup, y = as.numeric(z.dist)))+geom_beeswarm()+geom_boxplot()+
  ggtitle(paste("wilcox test p =",test$p.value,'medians',m1, ' ', m2,'es=',cd))+
  theme(plot.title = element_text(size = 5))+
  ylab('pairwise bouton z.dists')+
  ylim(0,8000)
ggsave('pwboutonzdists_1000sample.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')
