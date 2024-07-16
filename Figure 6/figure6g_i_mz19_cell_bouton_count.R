#Ablation Figure Collected Data
#started 04-26-2023
setwd("C:/Users/emmamtk/Dropbox (University of Michigan)/AhmedThorntonKolbe2022/CSVs for paper/Figure 6")
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
#20211201- VT033008 ablation with MZ19 label, luc RNAi control for DTA. anti GFP, anti dsRed, anti Chat
#20220511- VT033008 ablation with MZ19 label, empty attp40 control for DTA. anti GFP, anti dsRed, anti brp
#load data
mz19 <- read.csv('mz19_vt08_summary_data.csv')
#MZ19 bouton counts------
#shapiro test for normality
shapiro.test(filter(mz19,label=='control')$mz19.bouton.count) 
shapiro.test(filter(mz19,label=='DTA')$mz19.bouton.count) #both normal can use a t test
#t test for normal data
test <- t.test(filter(mz19,label=='control')$mz19.bouton.count,filter(mz19,label=='DTA')$mz19.bouton.count)
m1 <- median(filter(mz19,label=='control')$mz19.bouton.count)
m2 <- median(filter(mz19,label=='DTA')$mz19.bouton.count)
#cohens d for effect size
cd <- cohensd(filter(mz19, label=='control')$mz19.bouton.count,
              filter(mz19, label=='DTA')$mz19.bouton.count)

ggplot(mz19, aes(x= label, y=mz19.bouton.count))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,85)+
  theme_classic()+
  ylab('mz19 bouton count')+
  xlab(element_blank())+
  ggtitle(paste("t test p =",round(test$p.value,6),                 
                subtitle = paste('medians',m1, ' ', m2, 'es',cd)))

ggsave('mz19_boutoncount.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')  

#MZ19 cell count-------------
#mz19 cell count graph
shapiro.test(filter(mz19,label=='control')$mz19.cell.count) 
shapiro.test(filter(mz19,label=='DTA')$mz19.cell.count) #both normal can use a t test
#t test for normal data
test <- t.test(filter(mz19,label=='control')$mz19.cell.count,
               filter(mz19,label=='DTA')$mz19.cell.count)
m1 <- median(filter(mz19,label=='control')$mz19.cell.count, na.rm = T)
m2 <- median(filter(mz19,label=='DTA')$mz19.cell.count,na.rm = T)
#cohens d for effect size
cd <- cohensd(filter(mz19, label == 'control')$mz19.cell.count[!is.na(filter(mz19, label == 'control')$mz19.cell.count)],
        filter(mz19, label == 'DTA')$mz19.cell.count[!is.na(filter(mz19, label == 'DTA')$mz19.cell.count)])

ggplot(mz19, aes(x= label, y=mz19.cell.count))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,20)+
  theme_classic()+
  ylab('mz19 bouton per cell')+
  xlab(element_blank())+
  ggtitle(paste("t test p =",round(test$p.value,6),                 
                subtitle = paste('medians',m1, ' ', m2, 'es',cd)))

ggsave('mz19_cellcount.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')


