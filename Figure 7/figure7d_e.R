#20240521 FIJI measure ROI output processing
#originally written 05/31/2024
#load packages----
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(PairedData)
cohensd <- function(x,y){
  #absolute value of difference in group means
  diffinmeans <- abs(mean(x)-mean(y))
  #pooled standard deviation- root mean square of the standard deviations
  pooledsd <- sqrt(((sd(x)^2)+(sd(y)^2))/2)
  #cohens d
  return(diffinmeans/pooledsd)
}
#load data----
data <- read.csv('dpr4OE_summarydata.csv')
#cumulative distribution by sample---------- 
ggplot(filter(data, roi == 'mz19'), aes(x = kc.norm.f)) + 
  stat_ecdf(aes(group = interaction(genotype,sampleid,roi), color = genotype), geom = 'smooth') + 
  xlim(0,3) + 
  theme_classic()  + 
  scale_color_manual(values = c('#73CE4F','#C3C8C1')) +
  labs(x = "Normalized fluorescence", y = 'Probability') + 
  ggtitle('mz19 Cumulative Distribution Graphs by Sample') + 
  theme(plot.title = element_text(size = 10))
ggsave('mz19cumulativedistribution.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

ggplot(filter(data, roi == 'nonmz19'), aes(x = kc.norm.f)) + 
  stat_ecdf(aes(group = interaction(genotype,sampleid), color = genotype), geom = 'smooth') + 
  xlim(0,3) + 
  theme_classic()  + 
  scale_color_manual(values = c('#73CE4F','#C3C8C1')) +
  labs(x = "Normalized fluorescence", y = 'Probability') + 
  ggtitle('nonmz19 Cumulative Distribution Graphs by Sample') + 
  theme(plot.title = element_text(size = 10))
ggsave('nonmz19cumulativedistribution.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

#mean normalized fluorescence by hemisphere-------
datahemi <- unique(data[,c(1,2,4,6,7,17)])
for(i in 1:nrow(datahemi)){
  datahemi$kc.norm.f[i] <- mean(filter(data, sampleid == datahemi$sampleid[i],roi ==datahemi$roi[i])$kc.norm.f)
  datahemi$boutoncount[i] <- length(filter(data, sampleid == datahemi$sampleid[i],roi ==datahemi$roi[i])$kc.norm.f)
}

#shapiro test for normality
shapiro.test(filter(datahemi, roi== 'mz19',genotype == 'dpr4')$kc.norm.f)
shapiro.test(filter(datahemi, roi== 'mz19',genotype == 'empty')$kc.norm.f)
#data normal will use t test
test <- t.test(filter(datahemi, roi== 'mz19',genotype == 'dpr4')$kc.norm.f,
                    filter(datahemi, roi== 'mz19',genotype == 'empty')$kc.norm.f)
m1 <- median(filter(datahemi, roi== 'mz19',genotype == 'dpr4')$kc.norm.f)
m2 <- median(filter(datahemi, roi== 'mz19',genotype == 'empty')$kc.norm.f)

#cohen's D effect size
cd <- cohensd(filter(datahemi, roi == 'mz19', genotype == 'dpr4')$kc.norm.f,
        filter(datahemi, roi == 'mz19', genotype == 'empty')$kc.norm.f)

ggplot(filter(datahemi, roi == 'mz19'), aes(x= genotype, y = kc.norm.f))+
  geom_beeswarm(size=2)+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_boxplot()+
  theme_classic()+
  ylab('gamma kc in mz19 normalized fluorescence')+
  xlab(element_blank())+
  ylim(0,2)+
  ggtitle(paste("t test p =",round(test$p.value,5),                 
                subtitle = paste('medians',round(m1,4), ' ', round(m2,4), 'es',cd)))

ggsave('mz19_normalizedKCbyhemisphere.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')


#shapiro test for normality
shapiro.test(filter(datahemi, roi== 'nonmz19',genotype == 'dpr4')$kc.norm.f)
shapiro.test(filter(datahemi, roi== 'nonmz19',genotype == 'empty')$kc.norm.f)
#data normal will use t test
test <- t.test(filter(datahemi, roi== 'nonmz19',genotype == 'dpr4')$kc.norm.f,
               filter(datahemi, roi== 'nonmz19',genotype == 'empty')$kc.norm.f)
m1 <- median(filter(datahemi, roi== 'nonmz19',genotype == 'dpr4')$kc.norm.f)
m2 <- median(filter(datahemi, roi== 'nonmz19',genotype == 'empty')$kc.norm.f)
#cohen's D effect size
cd <- cohensd(filter(datahemi, roi == 'nonmz19', genotype == 'dpr4')$kc.norm.f,
        filter(datahemi, roi == 'nonmz19', genotype == 'empty')$kc.norm.f)

ggplot(filter(datahemi, roi == 'nonmz19'), aes(x= genotype, y = kc.norm.f))+
  geom_beeswarm(size=2)+
  stat_summary(fun = median, geom = "crossbar", width=.4, color = 'dark grey')+
  geom_boxplot()+
  theme_classic()+
  ylab('gamma kc in non mz19 normalized fluorescence')+
  xlab(element_blank())+
  ylim(0,2)+
  ggtitle(paste("t test p =",round(test$p.value,4),                 
                subtitle = paste('medians',round(m1,4), ' ', round(m2,4),'es',cd)))

ggsave('nonmz19_normalizedKCbyhemisphere.pdf',
       device = 'pdf',
       width = 4,height =6 ,units = 'in')

#mean normalized fluorescence by hemisphere- paired samples--------

#stats for paired wilcoxon signed rank test for paired samples
dp4test <- wilcox.test(filter(datahemi, roi == 'mz19', genotype == 'dpr4')$kc.norm.f,
            filter(datahemi, roi == 'nonmz19',genotype == 'dpr4')$kc.norm.f, paired = T)
emptytest <- wilcox.test(filter(datahemi, roi == 'mz19', genotype == 'empty')$kc.norm.f,
            filter(datahemi, roi == 'nonmz19',genotype == 'empty')$kc.norm.f, paired = T)
#make factors so that data plots in the right order
datahemi$roi <- factor(datahemi$roi, levels = c('nonmz19','mz19'))
datahemi$genotype <- factor(datahemi$genotype, levels = c('empty', 'dpr4'))
#plot paired data
ggplot(datahemi, aes(x=roi, y = kc.norm.f))+geom_beeswarm()+geom_line(aes(group = sampleid))+
  facet_wrap(~genotype)+theme_classic()+ylim(0,1.8)+
  ggtitle(paste('dpr4 paired p = ',round(dp4test$p.value,5),'con paired p =',round(emptytest$p.value,5)))
ggsave('paired_normalizedKCbyhemisphere.pdf',
       device = 'pdf',
       width = 8,height =6 ,units = 'in')
