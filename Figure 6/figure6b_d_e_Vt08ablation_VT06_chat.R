#20230308 Analysis 
#03 28 2023
setwd("C:/Users/emmamtk/Dropbox (University of Michigan)/AhmedThorntonKolbe2022/CSVs for paper/Figure 6")
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
#loading data
sumdata <- read.csv('20230308_new_summary_data.csv')
###bouton count plots####
#all/chat boutons scatter plot------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP')$chat.bout) 
shapiro.test(filter(sumdata, genotype=='DTA')$chat.bout) #both normal can use a t test
#t test for normal data
test <- t.test(filter(sumdata, genotype == 'caryP')$chat.bout,filter(sumdata, genotype == 'DTA')$chat.bout)
m1 <- median(filter(sumdata, genotype == 'caryP')$chat.bout)
m2 <- median(filter(sumdata, genotype == 'DTA')$chat.bout)
#cohen's d for effect size
cd <- cohensd(filter(sumdata, genotype == 'caryP')$chat.bout,
              filter(sumdata, genotype == 'DTA')$chat.bout)

ggplot(sumdata, aes(x= genotype, y=chat.bout))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,600)+
  theme_classic()+
  ylab('all ChAT+ boutons')+
  xlab(element_blank())+
  ggtitle(paste("t test p =",round(test$p.value,6),
                subtitle = paste('medians',m1, ' ', m2, 'es:',cd)))

ggsave('chat_bouton.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')
 
#vt033006 boutons scatter plot-----------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP')$vt033006.bout) 
shapiro.test(filter(sumdata, genotype=='DTA')$vt033006.bout) #both normal can use a t test
#t test for normal data
test <- t.test(filter(sumdata, genotype == 'caryP')$vt033006.bout,filter(sumdata, genotype == 'DTA')$vt033006.bout)
m1 <- median(filter(sumdata, genotype == 'caryP')$vt033006.bout)
m2 <- median(filter(sumdata, genotype=='DTA')$vt033006.bout)
#cohens d for effect size
cd <- cohensd(filter(sumdata, genotype == 'caryP')$vt033006.bout,
              filter(sumdata, genotype == 'DTA')$vt033006.bout)
ggplot(sumdata, aes(x= genotype, y=vt033006.bout))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  ylim(0,600)+
  theme_classic()+
  ylab('all vt033006+ boutons')+
  xlab(element_blank())+
  ggtitle(paste("t test p =",round(test$p.value,6),                 
                subtitle = paste('medians',m1, ' ', m2, 'es',cd)))

ggsave('vt033006_bouton.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')  


#vt033008 boutons scatter plot-----------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP')$vt033008.bout) 
shapiro.test(filter(sumdata, genotype=='DTA')$vt033008.bout) #both normal can use a t test
#t test for normal data
test <- t.test(filter(sumdata, genotype == 'caryP')$vt033008.bout,filter(sumdata, genotype == 'DTA')$vt033008.bout)
m1 <- median(filter(sumdata, genotype == 'caryP')$vt033008.bout)
m2 <- median(filter(sumdata, genotype=='DTA')$vt033008.bout)
#cohens d for effect size
cd <- cohensd(filter(sumdata, genotype == 'caryP')$vt033008.bout,
              filter(sumdata, genotype == 'DTA')$vt033008.bout)

ggplot(sumdata, aes(x= genotype, y=vt033008.bout))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,600)+
  theme_classic()+
  ylab('all vt033008+ boutons')+
  xlab(element_blank())+
  ggtitle(paste("t test p =",round(test$p.value,10),                 
                subtitle = paste('medians',m1, ' ', m2, 'es',cd)))

ggsave('vt033008_bouton.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')
  


#vt033006 only plot-----------------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP')$vt033006only.bout) #not normal will use wilcox rank test
shapiro.test(filter(sumdata, genotype=='DTA')$vt033006only.bout)
#wilcox rank test for non normal data
test <- wilcox.test(vt033006only.bout~genotype, data= sumdata, paired=FALSE, exact=FALSE, conf.int=TRUE)
m1 <- median(filter(sumdata, genotype == 'caryP')$vt033006only.bout)
m2 <- median(filter(sumdata, genotype=='DTA')$vt033006only.bout)
#cohens d for effect size
cd <- cohensd(filter(sumdata, genotype == 'caryP')$vt033006only.bout,
              filter(sumdata, genotype == 'DTA')$vt033006only.bout)

ggplot(sumdata, aes(x= genotype, y=vt033006only.bout))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,600)+
  theme_classic()+
  ylab('VT033006+ only boutons')+
  xlab(element_blank())+
  ggtitle(paste("wrst p =",round(test$p.value,6),
                subtitle = paste('medians',m1, ' ', m2, 'es',cd)))

ggsave('vt033006only_bouton.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

#survivor boutons (chat and vt033006 not vt033008 boutons)-------------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype == 'caryP')$chatonly.bout+
               filter(sumdata, genotype == 'caryP')$vt033006only.bout)
shapiro.test(filter(sumdata, genotype == 'DTA')$chatonly.bout+
               filter(sumdata, genotype == 'DTA')$vt033006only.bout)
#t test normal data
test <-  t.test(filter(sumdata, genotype == 'caryP')$chatonly.bout+
                  filter(sumdata, genotype=='caryP')$vt033006only.bout,
                filter(sumdata, genotype == 'DTA')$chatonly.bout+
                  filter(sumdata, genotype=='DTA')$vt033006only.bout)

m1 <- median(filter(sumdata, genotype == 'caryP')$chatonly.bout+
               filter(sumdata, genotype=='caryP')$vt033006only.bout)
m2 <- median(filter(sumdata, genotype == 'DTA')$chatonly.bout+filter(sumdata, genotype=='DTA')$vt033006only.bout)

#cohens d for effect size
cd <- cohensd(filter(sumdata, genotype == 'caryP')$chatonly.bout+
                filter(sumdata, genotype=='caryP')$vt033006only.bout,
              filter(sumdata, genotype == 'DTA')$chatonly.bout+
                filter(sumdata, genotype=='DTA')$vt033006only.bout)

ggplot(sumdata, aes(x= genotype, y=chatonly.bout+vt033006only.bout))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,600)+
  theme_classic()+
  ylab('ChAT only and VT033006 only boutons')+
  xlab(element_blank())+
  ggtitle(paste("t test p =",round(test$p.value,6),                 
                subtitle = paste('medians',m1, ' ', m2, 'es', cd)))

ggsave('survivor_boutons.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')


#chat only boutons---------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype == 'caryP')$chatonly.bout)
shapiro.test(filter(sumdata, genotype == 'DTA')$chatonly.bout)
#t test normal data
test <-  t.test(filter(sumdata, genotype == 'caryP')$chatonly.bout,
                filter(sumdata, genotype == 'DTA')$chatonly.bout)

m1 <- median(filter(sumdata, genotype == 'caryP')$chatonly.bout)
m2 <- median( filter(sumdata, genotype == 'DTA')$chatonly.bout)
#cohens d for effect size
cd <- cohensd(filter(sumdata, genotype == 'caryP')$chatonly.bout,
              filter(sumdata, genotype == 'DTA')$chatonly.bout)

ggplot(sumdata, aes(x= genotype, y=chatonly.bout))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,600)+
  theme_classic()+
  ylab('ChAT only boutons')+
  xlab(element_blank())+
  ggtitle(paste("t test p =",round(test$p.value,6),                 
                subtitle = paste('medians',m1, ' ', m2, 'es',cd)))

ggsave('chatonly_boutons.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')



#stacked bargraph plot for bouton counts----------------
#making data frame for plotting stacked bar graphs for bouton counts
b06only <- sumdata[,c(1,2,9,15)]
b06only$type <- 'vt033006'
b08only <- sumdata[,c(1,2,10,16)]
b08only$type <- 'vt033008'
b06and08 <- sumdata[,c(1,2,11,17)]
b06and08$type <- 'vt033006 and vt033008'
bchat <- sumdata[,c(1,2,8)]
bchat$n.soma <- 0
bchat$type <- 'chat only'

colnames(b06and08) <- c('Label','genotype', 'n.bout','n.soma','type')
colnames(b06only) <- c('Label','genotype', 'n.bout','n.soma','type')
colnames(b08only) <- c('Label','genotype', 'n.bout','n.soma','type')
colnames(bchat) <- c('Label','genotype', 'n.bout','n.soma','type')
boutcounts <- rbind(bchat, b06only, b06and08, b08only)
rm(b06and08)
rm(b08only)
rm(b06only)
rm(bchat)
write.csv(boutcounts,'20230308_boutoncounts.csv', row.names = FALSE)

boutcounts <- read.csv('20230308_boutoncounts.csv')

boutcounts$Label <- factor(boutcounts$Label, 
                           levels= unique(c(boutcounts$Label[boutcounts$genotype=='caryP'],boutcounts$Label[boutcounts$genotype=='DTA'])))
boutcounts$type <- factor(boutcounts$type, 
                          levels = c('chat only' ,'vt033008','vt033006 and vt033008' ,'vt033006'))



ggplot(boutcounts, 
       aes(fill=type, y=n.bout, x=Label, group= genotype)) + 
  geom_bar(position="stack", stat="identity")+
  ylab('number of calyx boutons')+ 
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x= element_blank(),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values=c('#f76f73','#ffb703','#07c4c5','#027fdc'))
ggsave('bouton_counts.pdf',
       device = 'pdf',
       width = 10,height =6 ,units = 'in')

ggplot(boutcounts, 
       aes(fill=type, y=n.soma, x=Label, group= genotype)) + 
  geom_bar(position="stack", stat="identity")+
  ylab('number of PN soma')+ 
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x= element_blank(),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values=c('#f76f73','#ffb703','#07c4c5','#027fdc'))
ggsave('cell_counts.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')






###soma count plots####
#vt033006 only soma plot--------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP', Label!=15)$vt033006only.soma) 
shapiro.test(filter(sumdata, genotype=='DTA')$vt033006only.soma)#not normal use wilcox rank sum test
#wilcox rank test for non normal data
test <- wilcox.test(vt033006only.soma~genotype, data= sumdata, paired=FALSE, exact=FALSE, conf.int=TRUE)
m1 <- median(filter(sumdata, genotype=='caryP', Label!=15)$vt033006only.soma)
m2 <- median(filter(sumdata, genotype=='DTA')$vt033006only.soma)
#cohens d for effect size
cd <- cohensd(filter(sumdata, genotype=='caryP')$vt033006only.soma,
                filter(sumdata, genotype=='DTA')$vt033006only.soma)

ggplot(filter(sumdata, Label!=15), aes(x= genotype, y=vt033006only.soma))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,215)+
  theme_classic()+
  ylab('VT033006 only soma count')+
  xlab(element_blank())+
  ggtitle(paste("wrst p =",round(test$p.value,6),
                subtitle= paste('median ',m1, ' ', m2 , 'es',cd)))

ggsave('vt033006only_soma.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')


#vt033008 only soma plot-------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP', Label!=15)$vt033008only.soma) 
shapiro.test(filter(sumdata, genotype=='DTA')$vt033008only.soma)#not normal use wilcox rank sum test
#wilcox rank test for non normal data
test <- wilcox.test(vt033008only.soma~genotype, data= sumdata, paired=FALSE, exact=FALSE, conf.int=TRUE)
m1 <- median(filter(sumdata, genotype=='caryP', Label!=15)$vt033008only.soma)
m2 <- median(filter(sumdata, genotype=='DTA')$vt033008only.soma)
#cohens d for effect size
cd <- cohensd(filter(sumdata, genotype=='caryP')$vt033008only.soma,
              filter(sumdata, genotype=='DTA')$vt033008only.soma)

ggplot(filter(sumdata, Label!=15), aes(x= genotype, y=vt033008only.soma))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,215)+
  theme_classic()+
  ylab('vt033008 only soma count')+
  xlab(element_blank())+
  ggtitle(paste("wrst p =",round(test$p.value,6),
                subtitle= paste('median ',m1, ' ', m2 )))

ggsave('vt033008only_soma.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

#vt033006 soma plot-------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP', Label!=15)$vt033006.soma) #not normal will use wilcox rank test
shapiro.test(filter(sumdata, genotype=='DTA')$vt033006.soma)
#wilcox rank test for non normal data
test <- wilcox.test(vt033006.soma~genotype, data= filter(sumdata,Label!=15), paired=FALSE, exact=FALSE, conf.int=TRUE)
m1 <- median(filter(sumdata, genotype=='caryP', Label!=15)$vt033006.soma)
m2 <- median(filter(sumdata, genotype=='DTA')$vt033006.soma)

ggplot(filter(sumdata,Label!=15), aes(x= genotype, y=vt033006.soma))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,215)+
  theme_classic()+
  ylab('VT033006 soma count')+
  xlab(element_blank())+
  ggtitle(paste("wrst p =",round(test$p.value,6),
                subtitle= paste('median ',m1, ' ', m2 )))

ggsave('vt033006_soma.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

#vt033008 soma plot-----------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP', Label!=15)$vt033008.soma) #not normal will use wilcox rank test
shapiro.test(filter(sumdata, genotype=='DTA')$vt033008.soma)
#wilcox rank test for non normal data
test <- wilcox.test(vt033008.soma~genotype, data= filter(sumdata,Label!=15), paired=FALSE, exact=FALSE, conf.int=TRUE)
m1 <- median(filter(sumdata, genotype=='caryP', Label!=15)$vt033008.soma)
m2 <- median(filter(sumdata, genotype=='DTA')$vt033008.soma)


ggplot(filter(sumdata,Label!=15), aes(x= genotype, y=vt033008.soma))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,215)+
  theme_classic()+
  ylab('vt033008 soma count')+
  xlab(element_blank())+
  ggtitle(paste("wrst p =",round(test$p.value,6),
                subtitle= paste('median ',m1, ' ', m2 )))

ggsave('vt033008_soma.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')


#overlap soma plot--------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP', Label!=15)$overlap.soma) 
shapiro.test(filter(sumdata, genotype=='DTA')$overlap.soma)#not normal use wilcox rank sum test
#wilcox rank test for non normal data
test <- wilcox.test(overlap.soma~genotype, data= sumdata, paired=FALSE, exact=FALSE, conf.int=TRUE)
m1 <- median(filter(sumdata, genotype=='caryP', Label!=15)$overlap.soma)
m2 <- median(filter(sumdata, genotype=='DTA')$overlap.soma)
#cohens d for effect size
cd <- cohensd(filter(sumdata, genotype=='caryP')$overlap.soma,
              filter(sumdata, genotype=='DTA')$overlap.soma)
ggplot(filter(sumdata, Label!=15), aes(x= genotype, y=overlap.soma))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  geom_boxplot()+
  ylim(0,215)+
  theme_classic()+
  ylab('overlap only soma count')+
  xlab(element_blank())+
  ggtitle(paste("wrst p =",round(test$p.value,6),
                subtitle= paste('median ',m1, ' ', m2 )))

ggsave('overlap_soma.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

#all soma plot---------------
#shapiro test for normality
shapiro.test(filter(sumdata, genotype=='caryP', Label!=15)$all.soma) 
shapiro.test(filter(sumdata, genotype=='DTA')$all.soma)#not normal use wilcox rank sum test
#wilcox rank test for non normal data
test <- wilcox.test(all.soma~genotype, data= sumdata, paired=FALSE, exact=FALSE, conf.int=TRUE)
m1 <- median(filter(sumdata, genotype=='caryP', Label!=15)$all.soma)
m2 <- median(filter(sumdata, genotype=='DTA')$all.soma)
ggplot(filter(sumdata, Label!=15), aes(x= genotype, y=all.soma))+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_beeswarm(size=2)+
  ylim(0,215)+
  theme_classic()+
  ylab('all soma count')+
  xlab(element_blank())+
  ggtitle(paste("wrst p =",round(test$p.value,6),
                subtitle= paste('median ',m1, ' ', m2 )))

ggsave('all_soma.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')
