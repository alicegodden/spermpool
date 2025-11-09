# Title: Phenotypic models- human
# Author: Dr. Daniel Marcu

#Phenotypic data on human sperm


library(tidyverse)
library(readxl)
library(dplyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(Rmisc)
library(DHARMa)
library(car)
library(ggbeeswarm)
library(Hmisc)


Phenotypic_test_result_V4<- read_excel("Library/CloudStorage/OneDrive-UniversityofEastAnglia/human_project/phenotypic_data/Phenotypic_test_result_V4.xlsx")


phenodata <- Phenotypic_test_result_V4
# Assuming your data frame is named df
# Subset the data for Time values 3 and 5 (selected T4 and T24)
phenodata <- Phenotypic_test_result_V4[Phenotypic_test_result_V4$Time %in% c(0, 1, 3, 5), ]

# Calculate the mean of VT_y for each Time value
mean_VT_y <- tapply(subset_data$VT_y, subset_data$Time, mean)


# Assuming your data frame is named Phenotypic_test_result_V4
# Subset the data for Time values 3 and 5
phenodata <- Phenotypic_test_result_V4[Phenotypic_test_result_V4$Time %in% c(0, 1, 2, 4, 6), ]

# Calculate the mean of VT_y for each Time value excluding NA values
mean_VT_y <- tapply(subset_data$AB_y, subset_data$Time, mean, na.rm = TRUE)

subset_data <- Phenotypic_test_result_V4[Phenotypic_test_result_V4$Time %in% c(0, 1, 2, 4, 6), ]

# Calculate the mean of VT_y for each Time value excluding NA values
mean_VT_y <- tapply(subset_data$TB_n, subset_data$Time, mean, na.rm = TRUE)

# Calculate the standard error of the mean for each Time value
se_VT_y <- tapply(subset_data$TB_n, subset_data$Time, function(x) sd(x, na.rm = TRUE) / sqrt(length(x)))

sd_VT_y <- tapply(subset_data$TB_n, subset_data$Time, sd, na.rm = TRUE)

ci_VT_y <- tapply(subset_data$TB_n, subset_data$Time, function(x) t.test(x, na.rm = TRUE)$conf.int)

# Now, 'se_VT_y' contains the standard error for each Time value




# Print the results
cat("Mean of VT_y for Time 3:", mean_VT_y[3], "\n")
cat("Mean of VT_y for Time 5:", mean_VT_y[5], "\n")



str(phenodata)


data1<- phenodata[ phenodata$Sample_ID != "GD164", ]

data2<- data1[ data1$Sample_ID != "GD44", ]

data3<- data2[ data2$Sample_ID != "GD028", ]

data4<- data3[ data3$Sample_ID != "MR96", ]

data5<- data4[ data4$Sample_ID != "MR70", ]

phenodata <- data5
str(phenodata)


phenodata$Sample_ID<-as.factor(phenodata$Sample_ID)

discard2<-which(phenodata$Time==3)
phenodata1<-phenodata[-discard2,]

discard3<-which(phenodata1$Time==0)
phenodata2<-phenodata1[-discard3,]

discard4<-which(phenodata2$Time==1)
phenodata3<-phenodata2[-discard4,]

discard5<-which(phenodata3$Time==2)
phenodata4<-phenodata3[-discard5,]



phenodata$Sample_ID<-as.factor(phenodata$Sample_ID)
phenodata$Time<-as.factor(phenodata$Time)

phenodata$VT_y <- round(phenodata$VT_y)
phenodata$VT_n <- round(phenodata$VT_n)

selected_data <- phenodata[, c("Time", "Sample_ID", "VT_y", "VT_n")]

subset_data <- na.omit(selected_data)


model1 <- glmer(cbind(VT_y, VT_n) ~ Time + (1 | Sample_ID), data = phenodata, family = binomial)
library(lmerTest)
install.packages("lmerTest")
# Obtain the summary of the model
model_summary <- summary(model1)


# Fit the GLMM
model1 <- glmer(cbind(VT_y, VT_n) ~ Time + (1 | Sample_ID), data = phenodata, family = binomial)


summary_with_df <- lmerTest::summary(model1, ddf = "Kenward-Roger")

summary(model1, ddf = "Kenward-Roger")

summary_with_df(model1)

# Use the lmerTest::summary function to obtain p-values and degrees of freedom
anova_with_df <- anova(model1)
# Print the summary with p-values and degrees of freedom
print(summary_with_df)




# Extract the degrees of freedom for fixed effects
df_fixed <- model_summary$coefficients$df[,"t value"]

# Print the degrees of freedom for fixed effects
print(df_fixed)



summary(model1)
anova(model1)

# Use the Anova function from the car package for Type III analysis
Anova_model1 <- Anova(model1, type = "III")

# Print the ANOVA table with p-values and degrees of freedom
print(Anova_model1)


str(phenodata)


phenodata$SB_y <- round(phenodata$SB_y)
phenodata$SB_n <- round(phenodata$SB_n)
selected_data <- phenodata[, c("Time", "Sample_ID", "SB_y", "SB_n")]

subset_data <- na.omit(selected_data)


#First model testing for difference in AB for times 1,4 and 6:
model1<-glmer(cbind(SB_y,SB_n)~Time+(1|Sample_ID), data=subset_data, family=binomial)
summary(model1)



anova(model1)
confint(model1)

#Example model for velocity parameter:
selected_data <- phenodata[, c("Time", "Sample_ID", "TB_y", "TB_n")]

filtered_data <- selected_data %>%
  filter(Time %in% c(4, 6))

#First model testing for difference in TB for times 1,4 and 6:
model2<-glmer(cbind(TB_y,TB_n)~Time+(1|Sample_ID), data=filtered_data, family=binomial)
summary(model2)
# Obtain confidence intervals for fixed effects
ci_fixed <- confint(model2)

# Print the confidence intervals
print(ci_fixed)

anova(model2)

confint(modelname)



#Example model for velocity parameter:
selected_data <- phenodata[, c("Time", "Sample_ID", "ALH")]

filtered_data <- selected_data %>%
  filter(Time %in% c(2,4, 6))

# View the filtered data
print(filtered_data)


model3<-lmer(ALH~Time+(1|Sample_ID), data=filtered_data)
summary(model3)

# Obtain confidence intervals for fixed effects
ci_fixed <- confint(model3)

# Print the confidence intervals
print(ci_fixed)

# Use the lmerTest::anova function to obtain p-values and degrees of freedom
anova_with_df <- anova(model3)

# Print the ANOVA table with p-values and degrees of freedom
print(anova_with_df)
anova(model3)
confint(model3)





Donor4<- SW_Parameter[ SW_Parameter$Time != "2", ]
Donor<- Donor[ Donor$time != "1", ]
Donor<- Donor[ Donor$time != "2", ]
unique(Donor[c("time")])




SW_Parameter

SW_Parameter

model3<-lmer(SurvivalSW~Time+(1|Sample_ID), data=Donor4)
summary(model3)
confint(model3)




df <- ll2$VCL

hist(df)




summary(model3)
confint(model3)

anova(model3)


model2<-glmer(VCL~Time+(1|Sample_ID), data=subset_data, family=poisson)
summary(model2)

files <- list.files(path = "/phenotypic_data/phenotypic_kinematics2", pattern = "*.xls", full.names = F)

# get the file names
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x)-4)
}
substrRight(files[1],16)

# load all CASA files
ll <- dir("/Users/work/OneDrive - University of East Anglia/human_project/phenotypic_data/phenotypic_kinematics/", full.names=T) %>% map_dfr(read_excel, skip=9, .id="sample")
str(ll) # missing the names though :( 
table(ll$sample)


ll$sample <- plyr::mapvalues(ll$sample, from=c(1:length(unique(ll$sample))), to=substrRight(files, 16))
ll$sample <- as.factor(ll$sample)
table(ll$sample)
ll$maleTreat <- as.factor(substr(ll$sample, start=1, stop=5))


ll1<- ll[ ll$maleTreat != "T48_G", ]
ll2<- ll1[ ll1$maleTreat != "T4_NA", ]
ll3<- ll2[ ll2$maleTreat != "T24_N", ]


unique(ll3[c("maleTreat")])

ll3 <- ll1

data1<- ll1[ ll1$sample != "GD164", ]

data2<- data1[ data1$sample != "GD44", ]

data3<- data2[ data2$sample != "GD028", ]

data4<- data3[ data3$sample != "MR96", ]

data5<- data4[ data4$sample != "MR70", ]


Donor <- data5 %>%
  mutate(time = factor(maleTreat, levels = c("BSA_G", "bsa_G", "SW_GD", "T24_G", "T4_GD", "T4_MR", "WASH_", "T24_N", "T4_NA"  ), 
                       labels = c("0", "0", "2","6", "4", "4", "1", "5", "3")))
head(Donor)
unique(Donor[c("time")])
Donor$time <- factor(Donor$time, levels = c("0", "1", "2", "3", "4", "5", "6"))

unique(Donor[c("time")])
Donor1 <- Donor


Donor<- Donor[ Donor$time != "0", ]
Donor<- Donor[ Donor$time != "1", ]
Donor<- Donor[ Donor$time != "2", ]
Donor<- Donor[ Donor$time != "4", ]
Donor<- Donor[ Donor$time != "6", ]
unique(Donor[c("time")])


Donor$VCL<-as.numeric(Donor$VCL)

Donor$time<-as.numeric(Donor$time)

VCL   VSL   VAP    LIN    STR    WOB   ALH   BCF
#Example model for velocity parameter:

model3<-lmer(BCF~time+(1|sample), data=Donor)
summary(model3)
confint(model3)

# Fit the lmer model to the data
fit <- lmer(BCF ~ time + (1|sample), data=Donor)

r <- cor(Donor$time, Donor$VCL)
p <- cor.test(Donor$time, Donor$VCL)$p.value

ggplot(Donor, aes(x=time, y=VCL)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  xlab("time") +
  ylab("VCL") +
  ggtitle(paste("Correlation between time and BCF (r =", round(r, 2), ", p =", round(p, 2), ")"))

#Model for count data (sperm concentration):
Donor_VCL <- na.omit(Donor1$VCL)

VCL1<- ggplot(data = Donor1 ,aes(x = time, y = VCL, fill = time))+
  scale_fill_viridis_d( option = "G")+
  geom_violin(alpha=0.6, position = position_dodge(width = .05),size=0.01,color=NA) +
  stat_summary(fun = mean, geom = 'pointrange', size = 0.45, color = 'black')+ stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3, size = 0.5, color = 'black')+
  ggbeeswarm::geom_quasirandom(shape = 21,size=0.1, dodge.width = .009, color="black",alpha=.009,show.legend = F)+
  theme_minimal()+
  ylab(  c("Curvilinear velocity (VCL μm/s) ")  )  +
  xlab(  c("Sperm Subpopulations")  )  +
  guides(fill = guide_legend(override.aes = list(alpha = 5.5,color="black"))) 

VCL_plot <- VCL1 + scale_x_discrete(labels=c("SL (Raw)", "SL (Washed)", "SL(Swim up)","LL (Four hour)", "LL (Twenty-four hour)")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(angle=45, hjust=0.5))  
VCL_plot


VSL1<- ggplot(data = Donor1 ,aes(x = time, y = VSL, fill = time))+
  scale_fill_viridis_d( option = "G")+
  geom_violin(alpha=0.6, position = position_dodge(width = .05),size=0.01,color=NA) +
  stat_summary(fun = mean, geom = 'pointrange', size = 0.45, color = 'black')+ stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3, size = 0.5, color = 'black')+
  ggbeeswarm::geom_quasirandom(shape = 21,size=0.1, dodge.width = .009, color="black",alpha=.009,show.legend = F)+
  theme_minimal()+
  ylab(  c("Straight line velocity (VSL μm/s) ")  )  +
  xlab(  c("Sperm Subpopulations")  )  +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5,color="black"))) 
VSL_plot <- VSL1 + scale_x_discrete(labels=c("SL (Raw)", "SL (Washed)", "SL(Swim up)","LL (Four hour)", "LL (Twenty-four hour)")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(angle=45, hjust=0.5))  
VSL_plot


WOB1<- ggplot(data = Donor1 ,aes(x = time, y = WOB, fill = time))+
  scale_fill_viridis_d( option = "G")+
  geom_violin(alpha=0.6, position = position_dodge(width = .05),size=0.01,color=NA) +
  stat_summary(fun = mean, geom = 'pointrange', size = 0.45, color = 'black')+ stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3, size = 0.5, color = 'black')+
  ggbeeswarm::geom_quasirandom(shape = 21,size=0.1, dodge.width = .009, color="black",alpha=.009,show.legend = F)+
  theme_minimal()+
  ylab(  c("Wobble (%) ")  )  +
  xlab(  c("Sperm Subpopulations")  )  +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5,color="black"))) 
WOB_plot <- WOB1 + scale_x_discrete(labels=c("SL (Raw)", "SL (Washed)", "SL(Swim up)","LL (Four hour)", "LL (Twenty-four hour)")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(angle=45, hjust=0.5))  
WOB_plot



VAP1<- ggplot(data = Donor1 ,aes(x = time, y = VAP, fill = time))+
  scale_fill_viridis_d( option = "G")+
  geom_violin(alpha=0.6, position = position_dodge(width = .05),size=0.01,color=NA) +
  stat_summary(fun = mean, geom = 'pointrange', size = 0.45, color = 'black')+ stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3, size = 0.5, color = 'black')+
  ggbeeswarm::geom_quasirandom(shape = 21,size=0.1, dodge.width = .009, color="black",alpha=.009,show.legend = F)+
  theme_minimal()+
  ylab(  c("Average path velocity (μm/s) ")  )  +
  xlab(  c("Sperm Subpopulations")  )  +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5,color="black"))) 
VAP_plot <- VAP1 + scale_x_discrete(labels=c("SL (Raw)", "SL (Washed)","SL(Swim up)", "LL (Four hour)", "LL (Twenty-four hour)")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(angle=45, hjust=0.5))  
VAP_plot


STR1<- ggplot(data = Donor1 ,aes(x = time, y = STR, fill = time))+
  scale_fill_viridis_d( option = "G")+
  geom_violin(alpha=0.6, position = position_dodge(width = .05),size=0.01,color=NA) +
  stat_summary(fun = mean, geom = 'pointrange', size = 0.45, color = 'black')+ stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3, size = 0.5, color = 'black')+
  ggbeeswarm::geom_quasirandom(shape = 21,size=0.1, dodge.width = .009, color="black",alpha=.009,show.legend = F)+
  theme_minimal()+
  ylab(  c("Straightness (μm/s) ")  )  +
  xlab(  c("Sperm Subpopulations")  )  +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5,color="black"))) 
STR_plot <- STR1 + scale_x_discrete(labels=c("SL (Raw)", "SL (Washed)", "SL(Swim up)","LL (Four hour)", "LL (Twenty-four hour)")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(angle=45, hjust=0.5))  
STR_plot

LIN1<- ggplot(data = Donor1 ,aes(x = time, y = LIN, fill = time))+
  scale_fill_viridis_d( option = "G")+
  geom_violin(alpha=0.6, position = position_dodge(width = .05),size=0.01,color=NA) +
  stat_summary(fun = mean, geom = 'pointrange', size = 0.45, color = 'black')+ stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3, size = 0.5, color = 'black')+
  ggbeeswarm::geom_quasirandom(shape = 21,size=0.1, dodge.width = .009, color="black",alpha=.009,show.legend = F)+
  theme_minimal()+
  ylab(  c("Linear coefficient (%) ")  )  +
  xlab(  c("Sperm Subpopulations")  )  +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5,color="black"))) 
LIN_plot <- LIN1 + scale_x_discrete(labels=c("SL (Raw)", "SL (Washed)", "SL(Swim up)","LL (Four hour)", "LL (Twenty-four hour)")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(angle=45, hjust=0.5))  
LIN_plot

BCF1<- ggplot(data = Donor1 ,aes(x = time, y = BCF, fill = time))+
  scale_fill_viridis_d( option = "G")+
  geom_violin(alpha=0.6, position = position_dodge(width = .05),size=0.01,color=NA) +
  stat_summary(fun = mean, geom = 'pointrange', size = 0.45, color = 'black')+ stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3, size = 0.5, color = 'black')+
  ggbeeswarm::geom_quasirandom(shape = 21,size=0.1, dodge.width = .009, color="black",alpha=.009,show.legend = F)+
  theme_minimal()+
  ylab(  c("Beat cross frequency (Hz) ")  )  +
  xlab(  c("Sperm Subpopulations")  )  +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5,color="black"))) 
BCF_plot <- BCF1 + scale_x_discrete(labels=c("SL (Raw)", "SL (Washed)","SL(Swim up)", "LL (Four hour)", "LL (Twenty-four hour)")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(angle=45, hjust=0.5))  
BCF_plot



AHL1<- ggplot(data = Donor1 ,aes(x = time, y = ALH, fill = time))+
  scale_fill_viridis_d( option = "G")+
  geom_violin(alpha=0.6, position = position_dodge(width = .05),size=0.01,color=NA) +
  stat_summary(fun = mean, geom = 'pointrange', size = 0.45, color = 'black')+ stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3, size = 0.5, color = 'black')+
  ggbeeswarm::geom_quasirandom(shape = 21,size=0.1, dodge.width = .009, color="black",alpha=.009,show.legend = F)+
  theme_minimal()+
  ylab(  c("Lateral head displacement (μm) ")  )  +
  xlab(  c("Sperm Subpopulations")  )  +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5,color="black"))) 
AHL_plot <- AHL1 + scale_x_discrete(labels=c("SL (Raw)", "SL (Washed)", "SL(Swim up)", "LL (Four hour)", "LL (Twenty-four hour)")) +
  theme(legend.position="none") + theme(axis.text.x=element_text(angle=45, hjust=0.5))  
AHL_plot



all_plot2<- multiplot(VCL_plot, VSL_plot, VAP_plot, LIN_plot, STR_plot, WOB_plot, AHL_plot, BCF_plot,  cols=4)
all_plot2




model2<-glmer(VCL~time+(1|sample), data=Donor, family=poisson)
summary(model2)



anova(model2)
confint(model2)

#Example model for velocity parameter:

model3<-lmer(VCL~Time+(1|Sample_ID), data=Dataname)
summary(model3)
anova(model3)

