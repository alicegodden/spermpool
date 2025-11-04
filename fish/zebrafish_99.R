# Title: 99th percentile analysis zebrafish
# Author: Dr. Daniel Marcu


#load all the packages 

library(dplyr)
library(ggplot2)
library(ggeasy)
library(cowplot)







#quantiles


duration = chromosome1$LRT_CO      # the eruption durations 
q1 <- quantile(duration, c(.999)) 
duration = chromosome2$LRT_CO     # the eruption durations 
q2 <- quantile(duration, c(.999)) 
duration = chromosome3$LRT_CO     # the eruption durations 
q3 <- quantile(duration, c(.999)) 
duration = chromosome4$LRT_CO     # the eruption durations 
q4 <- quantile(duration, c(.999)) 
duration = chromosome5$LRT_CO     # the eruption durations 
q5 <- quantile(duration, c(.999)) 
duration = chromosome6$LRT_CO     # the eruption durations 
q6 <- quantile(duration, c(.999)) 
duration = chromosome7$LRT_CO     # the eruption durations 
q7 <- quantile(duration, c(.999)) 
duration = chromosome8$LRT_CO     # the eruption durations 
q8 <- quantile(duration, c(.999)) 
duration = chromosome9$LRT_CO     # the eruption durations 
q9 <- quantile(duration, c(.999)) 
duration = chromosome10$LRT_CO     # the eruption durations 
q10 <- quantile(duration, c(.999)) 
duration = chromosome11$LRT_CO     # the eruption durations 
q11 <- quantile(duration, c(.999)) 
duration = chromosome12$LRT_CO     # the eruption durations 
q12 <- quantile(duration, c(.999)) 
duration = chromosome13$LRT_CO     # the eruption durations 
q13 <- quantile(duration, c(.999)) 
duration = chromosome14$LRT_CO     # the eruption durations 
q14 <- quantile(duration, c(.999)) 
duration = chromosome15$LRT_CO     # the eruption durations 
q15 <- quantile(duration, c(.999)) 
duration = chromosome16$LRT_CO     # the eruption durations 
q16 <- quantile(duration, c(.999)) 
duration = chromosome17$LRT_CO     # the eruption durations 
q17 <- quantile(duration, c(.999)) 
duration = chromosome18$LRT_CO     # the eruption durations 
q18 <- quantile(duration, c(.999)) 
duration = chromosome19$LRT_CO     # the eruption durations 
q19 <- quantile(duration, c(.999)) 
duration = chromosome20$LRT_CO     # the eruption durations 
q20 <- quantile(duration, c(.999)) 
duration = chromosome21$LRT_CO     # the eruption durations 
q21 <- quantile(duration, c(.999)) 
duration = chromosome22$LRT_CO     # the eruption durations 
q22 <- quantile(duration, c(.999)) 
duration = chromosome22$LRT_CO     # the eruption durations 
q22 <- quantile(duration, c(.999)) 
duration = chromosome23$LRT_CO     # the eruption durations 
q23 <- quantile(duration, c(.999)) 
duration = chromosome24$LRT_CO     # the eruption durations 
q24 <- quantile(duration, c(.999)) 
duration = chromosome25$LRT_CO     # the eruption durations 
q25 <- quantile(duration, c(.999)) 

(q1+q2+q3+q4+q5+q6+q7+q8+q9+q10+q11+q12+q13+q14+q15+q16+q17+q18+q19+q20+q21+q22+q23+q24+q25)/25

quant<- c(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18,q19,q20,q21,q22,q23,q24,q25)



T0_T4_hetPoolLikelihoods <- Male31

#convert all the relevant collums from character to numeric 

T0_T4_hetPoolLikelihoods$LRT_CO <- as.numeric(as.character(T0_T4_hetPoolLikelihoods$LRT_CO))

T0_T4_hetPoolLikelihoods$LRT_UC <- as.numeric(as.character(T0_T4_hetPoolLikelihoods$LRT_UC))

T0_T4_hetPoolLikelihoods$LRT_FU <- as.numeric(as.character(T0_T4_hetPoolLikelihoods$LRT_FU))

T0_T4_hetPoolLikelihoods$qual <- as.numeric(as.character(T0_T4_hetPoolLikelihoods$qual))

T0_T4_hetPoolLikelihoods$Maf_F <- as.numeric(as.character(T0_T4_hetPoolLikelihoods$Maf_F))

T0_T4_hetPoolLikelihoods$Maf_U <- as.numeric(as.character(T0_T4_hetPoolLikelihoods$Maf_U))

T0_T4_hetPoolLikelihoods$Maf_C <- as.numeric(as.character(T0_T4_hetPoolLikelihoods$Maf_C))

T0_T4_hetPoolLikelihoods$Maf_O <- as.numeric(as.character(T0_T4_hetPoolLikelihoods$Maf_O))

T0_T4_hetPoolLikelihoods$Tcov <- as.numeric(T0_T4_hetPoolLikelihoods$Tcov)

T0_T4_hetPoolLikelihoods$pos <- as.numeric(T0_T4_hetPoolLikelihoods$pos)


#filter for max coverage (xcoverage) x 2 =  386.21 for D1T0_T4

DT0_T4 <- T0_T4_hetPoolLikelihoods[(T0_T4_hetPoolLikelihoods[,3]<543.4),]

DT0_T4 <- na.omit(DT0_T4 )

mean(DT0_T4$Tcov)


#male 31 543.4

#male 32 532.850
#male 33  555.120

#Donor 11 200.684475*2 = 398.88
#Donor 8 206.69 *2 = 413.38
#Donor 13 174.9 *2 = 349.89
#Donor 12 176.31 *2 = 352.62

#Donor1 261.17*2=522.34
#Donor2 221.03*2=442.06



#D1T0_T4 = (95.2936+97.8141)x2 = 386.21
#D1T0_T24 = (95.2936+77.70)x2 = 345.98

#D2T0_T2 = (75.401+77.2314)*2 = 305.2648
#D2T0_T4 = (75.401+96.1164)*2 = 343.0348
#D2T0_T8 = (75.401+68.310)*2 = 287.422
#D2T0_T12 = (75.401+77.2314)*2 = 305.2648
#D2T0_T24 = (75.401+71.664)*2 = 294.13
#D2T0_T48 = (75.401+68.2987)*2 = 287.3994



DT0_T4_filter1 <- DT0_T4[(DT0_T4[,18]>0.5),]

DT0_T4_filter <- DT0_T4_filter1[(DT0_T4_filter1[,11]<0.01),]



write.csv(DT0_T4_filter,'/Users/work/OneDrive - University of East Anglia/human_project/result/Likelihood_ratio_results/overlap_all/Male12_LRTs.csv')


#Ensure the Tcov, LRT_CO, and pos are numerical 
#Must run each chromosome at a time

# 15.2814 T2  11.20498
# 15.12895 T4 11.06666
#  15.07233 T8 11.13373
# 15.48539 T12 11.36361 
# 18.40472  T24 11.3188
# 20.87731 T48 11.47106 

#mean 11.25981

# D1 T24 17.65599 - 99th quuantile 11.42835  - 10.47786 
# D1 T4 16.82565  - 11.31484   - 10.43837 99.9

# 11.37/10.46 and 11.30 for donor 1 and 2

#  Male 8-  21.34524 + 12.81651=17.08   99 quantile
# Male 11 - (12.19478+11.64023)/2 = 11.91751
#  male 12 - 19.20963 (19.20963+11.5195)/2 = 15.36457
#  male 13 - 11.75276  (11.75276+11.40678)/2 = 11.6

#male 31
#male 32 11.54181, 99.9 37.87657


D1T24 -hetPool


#Chromosome 1

chromosome1 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007112.7"),]

highlight_chr1 <- chromosome1 %>% 
  filter(LRT_CO>=13.99598  )

ch1_ <- chromosome1 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr1, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 1") +
  ggeasy::easy_center_title()
ch1 <- ch1_ + geom_hline(yintercept=13.99598  , linetype="dashed", color = "red")
ch1 



#Chromosome 2

chromosome2 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007113.7"),]

highlight_chr2 <- chromosome2 %>% 
  filter(LRT_CO>=13.99598   )


ch2_ <- chromosome2 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr2, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 2") +
  ggeasy::easy_center_title()
ch2 <- ch2_ + geom_hline(yintercept=13.99598     , linetype="dashed", color = "red")
ch2

#Chromosome 3

chromosome3 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007114.7"),]

highlight_chr3 <- chromosome3 %>% 
  filter(LRT_CO>=13.99598 )


ch3_ <- chromosome3 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr3, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 3") +
  ggeasy::easy_center_title()
ch3 <- ch3_ + geom_hline(yintercept=13.99598   , linetype="dashed", color = "red")

ch3

#Chromosome 4

chromosome4 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007115.7"),]

highlight_chr4 <- chromosome4 %>% 
  filter(LRT_CO>=13.99598   )


ch4_ <- chromosome4 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr4, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 4") +
  ggeasy::easy_center_title()
ch4 <- ch4_ + geom_hline(yintercept=13.99598 , linetype="dashed", color = "red")
ch4
#Chromosome 5

chromosome5 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007116.7"),]

highlight_chr5 <- chromosome5 %>% 
  filter(LRT_CO>=13.99598)


ch5_ <- chromosome5 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr5, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 5") +
  ggeasy::easy_center_title()
ch5 <- ch5_ + geom_hline(yintercept=13.99598    , linetype="dashed", color = "red")
ch5



#Chromosome 6

chromosome6 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007117.7"),]

highlight_chr6 <- chromosome6 %>% 
  filter(LRT_CO>=13.99598 )


ch6_ <- chromosome6 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr6, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 6") +
  ggeasy::easy_center_title()

ch6 <- ch6_ + geom_hline(yintercept=13.99598  , linetype="dashed", color = "red")
ch6


#Chromosome 7

chromosome7 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007118.7"),]

highlight_chr7 <- chromosome7 %>% 
  filter(LRT_CO>=13.99598)


ch7_ <- chromosome7 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr7, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 7") +
  ggeasy::easy_center_title()
ch7 <- ch7_ + geom_hline(yintercept=13.99598 , linetype="dashed", color = "red")
ch7

#Chromosome 8

chromosome8 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007119.7"),]

highlight_chr8 <- chromosome8 %>% 
  filter(LRT_CO>=13.99598 )


ch8_ <- chromosome8 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr8, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 8") +
  ggeasy::easy_center_title()

ch8 <- ch8_ + geom_hline(yintercept=13.99598   , linetype="dashed", color = "red")
ch8






#Chromosome 9

chromosome9 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007120.7"),]

highlight_chr9 <- chromosome9 %>% 
  filter(LRT_CO>=13.99598 )


ch9_ <- chromosome9 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr9, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 9") +
  ggeasy::easy_center_title()

ch9 <- ch9_+ geom_hline(yintercept=13.99598, linetype="dashed", color = "red")
ch9
#Chromosome 10

chromosome10 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007121.7"),]

highlight_chr10 <- chromosome10 %>% 
  filter(LRT_CO>=13.99598)


ch10_ <- chromosome10 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr10, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 10") +
  ggeasy::easy_center_title()

ch10 <- ch10_+ geom_hline(yintercept=13.99598   , linetype="dashed", color = "red")
ch10



#Chromosome 11

chromosome11 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007122.7"),]

highlight_chr11 <- chromosome11 %>% 
  filter(LRT_CO>=13.99598 )


ch11_ <- chromosome11 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr11, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 11") +
  ggeasy::easy_center_title()
ch11 <- ch11_+ geom_hline(yintercept=13.99598    , linetype="dashed", color = "red")
ch11
#Chromosome 12

chromosome12 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007123.7"),]

highlight_chr12 <- chromosome12 %>% 
  filter(LRT_CO>=13.99598 )


ch12_ <- chromosome12 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr12, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 12") +
  ggeasy::easy_center_title()
ch12 <- ch12_+ geom_hline(yintercept=13.99598    , linetype="dashed", color = "red")
ch12



#Chromosome 13

chromosome13 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007124.7"),]

highlight_chr13 <- chromosome13 %>% 
  filter(LRT_CO>=13.99598)


ch13_ <- chromosome13 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr13, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 13") +
  ggeasy::easy_center_title()

ch13 <- ch13_+ geom_hline(yintercept=13.99598  , linetype="dashed", color = "red")
ch13
#Chromosome 14

chromosome14 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007125.7"),]

highlight_chr14 <- chromosome14 %>% 
  filter(LRT_CO>=13.99598 )


ch14_ <- chromosome14 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr14, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 14") +
  ggeasy::easy_center_title()
ch14 <- ch14_+ geom_hline(yintercept=13.99598  , linetype="dashed", color = "red")
ch14
#Chromosome 15

chromosome15 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007126.7"),]

highlight_chr15 <- chromosome15 %>% 
  filter(LRT_CO>=13.99598  )


ch15_ <- chromosome15 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr15, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 15") +
  ggeasy::easy_center_title()

ch15 <- ch15_+ geom_hline(yintercept=13.99598   , linetype="dashed", color = "red")
ch15

#Chromosome 16

chromosome16 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007127.7"),]

highlight_chr16 <- chromosome16 %>% 
  filter(LRT_CO>=13.99598)


ch16_ <- chromosome16 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr16, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 16") +
  ggeasy::easy_center_title()
ch16 <- ch16_+ geom_hline(yintercept=13.99598     , linetype="dashed", color = "red")
ch16

#Chromosome 17

chromosome17 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007128.7"),]

highlight_chr17 <- chromosome17 %>% 
  filter(LRT_CO>=13.99598 )


ch17_ <- chromosome17 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr17, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 17") +
  ggeasy::easy_center_title()

ch17 <- ch17_+ geom_hline(yintercept=13.99598    , linetype="dashed", color = "red")
ch17


#Chromosome 18

chromosome18 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007129.7"),]

highlight_chr18 <- chromosome18 %>% 
  filter(LRT_CO>=13.99598  )


ch18_ <- chromosome18 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr18, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 18") +
  ggeasy::easy_center_title()
ch18 <- ch18_+ geom_hline(yintercept=13.99598   , linetype="dashed", color = "red")
ch18
#Chromosome 19

chromosome19 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007130.7"),]

highlight_chr19 <- chromosome19 %>% 
  filter(LRT_CO>=13.99598)

ch19_ <- chromosome19 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr19, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 19") +
  ggeasy::easy_center_title()
ch19 <- ch19_+ geom_hline(yintercept=13.99598  , linetype="dashed", color = "red")
ch19
#Chromosome 20

chromosome20 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007131.7"),]

highlight_chr20 <- chromosome20 %>% 
  filter(LRT_CO>=13.99598  )

ch20_ <- chromosome20 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr20, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 20") +
  ggeasy::easy_center_title()
ch20 <- ch20_+ geom_hline(yintercept=13.99598   , linetype="dashed", color = "red")
ch20


#Chromosome 21 
chromosome21 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007132.7"),]

highlight_chr21 <- chromosome21 %>% 
  filter(LRT_CO>=13.99598)


ch21_ <- chromosome21 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr21, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 21") +
  ggeasy::easy_center_title()

ch21 <- ch21_+ geom_hline(yintercept=13.99598   , linetype="dashed", color = "red")
ch21
#Chromosome 22 
chromosome22 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007133.7"),]

highlight_chr22 <- chromosome22 %>% 
  filter(LRT_CO>=13.99598  )


ch22_ <- chromosome22 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr22, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 22") +
  ggeasy::easy_center_title()
ch22 <- ch22_+ geom_hline(yintercept=13.99598    , linetype="dashed", color = "red")
ch22
#Chromosome X 
chromosome23 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007134.7"),]

highlight_chr23 <- chromosome23 %>% 
  filter(LRT_CO>=13.99598 )


ch23_ <- chromosome23 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr23, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 23") +
  ggeasy::easy_center_title()
ch23 <- ch23_ + geom_hline(yintercept=13.99598 , linetype="dashed", color = "red")
ch23
#Chromosome 24 
chromosome24 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007135.7"),]

highlight_chr24 <- chromosome24 %>% 
  filter(LRT_CO>=13.99598  )


ch24_ <- chromosome24 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr24, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 24") +
  ggeasy::easy_center_title()

ch24 <- ch24_ + geom_hline(yintercept=13.99598 , linetype="dashed", color = "red" )

ch24 


#Chromosome 25 
chromosome25 <- DT0_T4_filter[ DT0_T4_filter$ref %in% c("NC_007136.7"),]

highlight_chr25 <- chromosome25 %>% 
  filter(LRT_CO>=13.99598  )


ch25_ <- chromosome25 %>% 
  ggplot(aes(x=pos,y=LRT_CO)) + 
  geom_point(alpha=0.1) +  
  geom_point() + geom_point(color = 'grey') + 
  geom_point(data=highlight_chr25, 
             aes(x=pos,y=LRT_CO), 
             color='red',
             size=1 ) +
  theme_bw()+
  xlab("Position (bp) ")+
  labs(y=expression("LRT"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom") +   ggtitle("Chromosome 25") +
  ggeasy::easy_center_title()

ch25 <- ch25_+ geom_hline(yintercept=13.99598 , linetype="dashed", color = "red" )

ch25 


plot <- plot_grid(ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9, ch10, ch11, ch12, ch13, ch14, ch15, ch16, ch17, ch18, ch19, ch20, ch21, ch22, ch23, ch24, ch25)

plot +  ggtitle("Male 31, C vs O ") + ggeasy::easy_center_title() 



final1 = rbind(chromosome1,chromosome2,chromosome3,chromosome4,chromosome5,chromosome6,chromosome7,chromosome8,chromosome9,chromosome10,chromosome11,chromosome12,chromosome13,chromosome14,chromosome15,chromosome16,chromosome17,chromosome18,chromosome19,chromosome20,chromosome21,chromosome22,chromosomeX,chromosomeY)

male31_sig = rbind(highlight_chr1,highlight_chr2,highlight_chr3,highlight_chr4,highlight_chr5,highlight_chr6,highlight_chr7,highlight_chr8,highlight_chr9,highlight_chr10,highlight_chr11,highlight_chr12,highlight_chr13,highlight_chr14,highlight_chr15,highlight_chr16,highlight_chr17,highlight_chr18,highlight_chr19,highlight_chr20,highlight_chr21,highlight_chr22,highlight_chr23,highlight_chr24, highlight_chr25)



write.csv(male31_sig,'/Users/work/OneDrive - University of East Anglia/human_project/result/Likelihood_ratio_results/male31_sig.csv')



----
  
  
  
library(GenomicRanges)
library(TREGEL) # load package



library(ggplot2)

install.packages("plyr")                           # Install plyr package
library("plyr")     


install.packages("dplyr")                          # Install dplyr package
library("dplyr")   

# Load the loci file and then substitute a numeric chromsome rather than a NC_0000...)
Donor <- male31_sig %>%
  mutate(chr = factor(ref, levels = c("NC_007112.7", "NC_007113.7", "NC_007114.7", "NC_007115.7", "NC_007116.7", "NC_007117.7", "NC_007118.7", "NC_007119.7", "NC_007120.7", "NC_007121.7", "NC_007122.7", "NC_007123.7", "NC_007124.7", "NC_007125.7", "NC_007126.7","NC_007127.7", "NC_007128.7", "NC_007129.7", "NC_007130.7", "NC_007131.7", "NC_007132.7", "NC_007133.7", "NC_007134.7", "NC_007135.7", "NC_007136.7"), 
                      labels = c("1", "2", "3", "4", "5", "6", "7", "8","9","10","11", "12","13", "14", "15", "16","17","18","19","20","21","22","23","24", "25")))
head(Donor)

#include only the columns needed
Donor_Loci <- subset(Donor , select = c(chr, pos))
head(Donor_Loci)
tail(Donor_Loci)





# change name of the columns 
Donor1_Loci <- Donor_Loci %>% 
  rename(
    start = pos
  )



Donor1_Loci<- plyr::rename(Donor_Loci, c("pos" = "start"))


head(Donor1_Loci)

Donor1_Loci$start <- as.numeric(Donor1_Loci$start)


# add a End collumn (this will be the same as the start - TREGEL requeires a end column)
Donor1_Loci$end <- Donor1_Loci$start + 10000

head(Donor1_Loci)

na.omit(Donor1_Loci)
# create a unique ID for each row (TREGEL requires this)

Donor1_Loci$ID <- sample(1:nrow(Donor1_Loci), nrow(Donor1_Loci), F)


head(Donor1_Loci)

Donor1_Loci$ID <- as.numeric(Donor1_Loci$ID)
Donor1_Loci$chr <- as.numeric(Donor1_Loci$chr)

# make Grange files for the data
Male31_gene <- makeGRangesFromDataFrame(Donor1_Loci, keep.extra.columns=TRUE)

head(Male31_gene)
# save the GRange file
saveRDS(Male31_gene, file = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/TREGEL/extdata/query/Loci.RDS")



# load the gene annotation file 

Gene_annotation1 <- plyr::rename(gene_annotation_zebrafish, c("Gene.stable.ID" = "Gene", "Gene.start..bp." = "Start", "Gene.end..bp." = "End", "Gene.type" = "Element", "Chromosome.scaffold.name" = "Chromosome" ))



# create a unique ID for each row (TREGEL requires this)

Gene_annotation1$ID <- sample(1:nrow(Gene_annotation1), nrow(Gene_annotation1), F)
head(Gene_annotation1)


Gene_annotation1$Start <- as.numeric(Gene_annotation1$Start)
Gene_annotation1$End <- as.numeric(Gene_annotation1$End)
Gene_annotation1$Chromosome <- as.numeric(Gene_annotation1$Chromosome)
Gene_annotation1$ID <- as.numeric(Gene_annotation1$ID)


head(Gene_annotation1)

# change name of the columns 

Gene_annotation1<- plyr::rename(Gene_annotation1, c("Chromosome" = "chr"))



head(Gene_annotation1)

# make Grange files for the data
Gene_annotation2 <- makeGRangesFromDataFrame(Gene_annotation1, keep.extra.columns=TRUE)

head(Gene_annotation2)




# save the GRange file
saveRDS(Gene_annotation2, file = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/TREGEL/extdata/subject/Gene_annotation.RDS")



myQueryFolder <- file.path(system.file('extdata', package = 'TREGEL'),"query") #  this is loci file (make sure it is in the right location)
mySubjectFolder <- file.path(system.file('extdata', package = 'TREGEL'),"subject") #  this is the gene annotation file (make sure it is in the right location)

myTREGEL <- TREGELDataSetFromDir(queryFolder=myQueryFolder,subjectFolder=mySubjectFolder)


# check that where the files are 
S4Vectors::metadata(myTREGEL)

# read the dataset
myTREGEL <- loadAnnotations(myTREGEL)

myTREGEL

# check what files are there 
names(myTREGEL)

# find overlap ovelap 
myTREGEL <- fOverlaps(myTREGEL)

head(S4Vectors::metadata(myTREGEL)$detailDT)

overlap <- (S4Vectors::metadata(myTREGEL)$detailDT)

overlap








# change name of the collumns


overlap2<- plyr::rename(myTREGEL, c("ID" = "subjID"))



overlap2 <- myTREGEL %>% 
  rename(
    
    ID = subjID)




head(overlap2)

# merge the gene annotation to loci file to add the extra info
Overlap <- merge(Gene_annotation2,overlap2, by = "ID")

head(Overlap)


# remove duplicates based on ID
Overlap1 <- Overlap[!duplicated(Overlap$Gene), ]
Overlap1


#save overlap file
write.csv(Overlap1,'/Users/work/OneDrive - University of East Anglia/human_project/result/Likelihood_ratio_results/overlap_all/Male8_overlap.csv')

# check how many elements there are 

unique(D2T48_gene[c("Element")])


Overlap1 <- D2T48_gene

# count each elemenet 
nrow(Overlap1[ Overlap1$Element == "protein_coding", ])

nrow(Overlap1[ Overlap1$Element == "lncRNA", ])

nrow(Overlap1[ Overlap1$Element == "unprocessed_pseudogene", ])

nrow(Overlap1[ Overlap1$Element == "transcribed_unprocessed_pseudogene", ])

nrow(Overlap1[ Overlap1$Element == "transcribed_processed_pseudogene", ])

nrow(Overlap1[ Overlap1$Element == "unitary_pseudogene", ])

nrow(Overlap1[ Overlap1$Element == "transcribed_unitary_pseudogene", ])

nrow(Overlap1[ Overlap1$Element == "processed_pseudogene", ])

nrow(Overlap1[ Overlap1$Element == "TEC", ])

nrow(Overlap1[ Overlap1$Element == "miRNA", ])

nrow(Overlap1[ Overlap1$Element == "IG_V_pseudogene", ])

nrow(Overlap1[ Overlap1$Element == "polymorphic_pseudogene", ])

nrow(Overlap1[ Overlap1$Element == "snoRNA", ])

# create data frame w/ the information above
df <- data.frame(Element=c("protein_coding", 
                           "lncRNA", 
                           "unprocessed_pseudogene",
                           "transcribed_unprocessed_pseudogene",
                           "transcribed_unprocessed_pseudogene", 
                           "unprocessed_pseudogene",
                           "processed_pseudogene",
                           "TEC"),
                 Occurrence=c(364, 201, 13, 13, 3, 13, 5, 1))

df

# create graph with the information from above
p<-ggplot(data=df, aes(x=Occurrence, y=Element)) +
  geom_bar(stat="identity") + theme_minimal() + labs(title = "D1T4")
p







# change name of the collumns

Gene_annotation1 <- Gene_ann_zebrafish %>% 
  rename(
    
    Gene = Gene.stable.ID,
    Start = Gene.start..bp.,
    Transcript = Transcript.stable.ID,
    End = Gene.end..bp.,
    Element = Gene.type,
    Chromosome = Chromosome.scaffold.name)

# create a unique ID for each row (TREGEL requires this)

Gene_annotation1$ID <- sample(1:nrow(Gene_annotation1), nrow(Gene_annotation1), F)
head(Gene_annotation1)

# make Grange files for the data
Gene_annotation2 <- makeGRangesFromDataFrame(Gene_annotation1, keep.extra.columns=TRUE)

head(Gene_annotation2)




# save the GRange file
saveRDS(Gene_annotation2, file = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/TREGEL/extdata/subject/Gene_annotation.RDS")

