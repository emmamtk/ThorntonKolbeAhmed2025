setwd("C:/Users/emmamtk/Desktop/CSVs for paper/Figure 6")
#load packages
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
cohensd <- function(x,y){
  #absolute value of difference in group means
  diffinmeans <- abs(mean(x)-mean(y))
  #pooled standard deviation- root mean square of the standard deviations
  pooledsd <- sqrt(((sd(x)^2)+(sd(y)^2))/2)
  #cohens d
  return(diffinmeans/pooledsd)
}
#experiments
#20220511- VT033008 ablation with MZ19 label, empty attp40 control for DTA. anti GFP, anti dsRed, anti brp
#20220819- VT033008 ablation with 91G04 label, empty attp40 control for DTA. anti GFP, anti dsRed, anti brp

#load data
nc82sum <- read.csv('nc82_ca_volume_summary.csv')
#plot nc82 measurements------------

#shapiro test for normality
shapiro.test(filter(nc82sum,genotype=='caryP')$norm.nc82) 
shapiro.test(filter(nc82sum,genotype=='DTA')$norm.nc82) #not normal will use wilcox test
#wilcox ranks sum test for non normal data
test <- wilcox.test(filter(nc82sum,genotype=='caryP')$norm.nc82,filter(nc82sum,genotype=='DTA')$norm.nc82)
m1 <- median(filter(nc82sum,genotype=='caryP')$norm.nc82)
m2 <- median(filter(nc82sum,genotype=='DTA')$norm.nc82)
#cohensd effect size
cd <- cohensd(filter(nc82sum, genotype == 'caryP')$norm.nc82,
              filter(nc82sum, genotype == 'DTA')$norm.nc82)
ggplot(nc82sum, aes(x= genotype, y=norm.nc82))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,1.7)+
  theme_classic()+
  ylab('normalized nc82 intensity')+
  xlab(element_blank())+
  ggtitle(paste("wilcox p =",round(test$p.value,6),                 
                subtitle = paste('medians',m1, ' ', m2)))

ggsave('normalized nc82.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')





#plot nc82 volume measurements--------
#shapiro test for normality
shapiro.test(filter(nc82sum, genotype=='caryP')$volume) 
shapiro.test(filter(nc82sum, genotype=='DTA')$volume) #both normal can use a t test
#t test for normal data
test <- t.test(filter(nc82sum, genotype=='caryP')$volume,filter(nc82sum, genotype=='DTA')$volume)
m1 <- median(filter(nc82sum, genotype=='caryP')$volume)
m2 <- median(filter(nc82sum, genotype=='DTA')$volume)
#cohens d effect size
cd <- cohensd(filter(nc82sum, genotype == 'caryP')$volume,
        filter(nc82sum, genotype == 'DTA')$volume)

ggplot(nc82sum, aes(x= genotype, y=volume))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,20000)+
  theme_classic()+
  ylab('volume')+
  xlab(element_blank())+
  ggtitle(paste("t p =",round(test$p.value,6),                 
                subtitle = paste('med',round(m1,1), ' ', round(m2,1), 'es', cd)))
ggsave('nc82_calyxvolume.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

