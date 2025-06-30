# Author : Sara Irish & Alice Godden
# Title: GLMM analyses on TE counts that have been filtered for high coverage confident calls, and the experimental groups retain TEs not in the control sample

library(dplyr)
library(tidyr)
library(lubridate)
library(survival)
library(survminer)
library(ggplot2)
library(dabestr)
library(lme4)
library(coxme)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(car)


# read in your data
df <- read.csv("zebrafish_TE_counts.csv")
df
head (df)
str(df)


#repro_long is the name of the new longform dataset

repro_long <- df %>% 
  
  pivot_longer(
    
    cols = c(`Total`,`DNA`, `LINE`, `SINE`, `LTR`, `SATELLITTE`, `RC`), #list all TE columns in dataset here RC for fish
    
    names_to = "TE_type", #this creates a column named TE_type, containing Total and TE by family
    
    values_to = "value") #this creates a column named value which contains the offspring numbers produced each day



head(repro_long)


# plot the TE data

repro_plot <-ggplot(data=repro_long, aes(x=TE_type, y=value, group=Time, color=Time)) +
  
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5))+#adds raw data to plot
  
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1.2, width=0.0, position = position_dodge(0.5)) +#calculates and adds error bars
  
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) + #calculates and adds mean points
  
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1.2, position = position_dodge(0.5)) + #adds lines bbetween mean points
  
  theme_classic()+#theme that eliminates gridlines
  
  labs(y="TE count", x="TE", title ="Zebrafish- TE counts- unique to treatment group") + #axis labels
  
  theme(
    # Make plot title bold and adjust size
    plot.title = element_text(size = 18, face = "bold"), # Added face = "bold"
    
    # Make y-axis title bold and adjust size
    axis.title.y = element_text(size = 16, face = "bold"), # Added face = "bold"
    
    # Make x-axis title bold and adjust size
    axis.title.x = element_text(size = 16, face = "bold"), # Added face = "bold"
    
    # Make axis tick text bold and adjust size
    axis.text = element_text(size = 14, face = "bold"), # Added face = "bold" for both x and y axis text
    
    # Make legend text bold and adjust size
    legend.text = element_text(size = 14, face = "bold"), # Added face = "bold"
    
    # Make legend title bold (if you have one, e.g., using labs(color = "Time Point"))
    legend.title = element_text(size = 14, face = "bold"), # Often forgotten!
    
    legend.position = c(0.8, 0.9) #x.y coordinates for where you want the legend to appear
  )
# generate plot
repro_plot


# Save the plot to a file
ggsave(
  filename = "zfish_TE_count_plot.png",   # The name of your file and its extension
  plot = repro_plot,            # The ggplot object you want to save
  width = 8,                    # Width in inches (adjust as needed)
  height = 6,                   # Height in inches (adjust as needed)
  dpi = 600                     # Resolution (dots per inch) for raster images
)

# now make a histogram
hist <- hist(repro_long$value)
ggsave(
  filename = "zfish_TE_count_hist.png",   # The name of your file and its extension
  plot = repro_plot,            # The ggplot object you want to save
  width = 8,                    # Width in inches (adjust as needed)
  height = 6,                   # Height in inches (adjust as needed)
  dpi = 600                     # Resolution (dots per inch) for raster images
)

#models
repro_long$Time <- as.factor(repro_long$Time)
##Then, change the variables names to numbers
levels(repro_long$Time) <- list('1' = 'C', '2' = 'O') # for MC-'1' = 'raw', '2' = 'central', '3' ='outer'
levels(repro_long$Time)

# running poisson model
rep_m1 <- glmer(value ~ Time + (1 | MaleID), family = poisson, data = repro_long)
summary(rep_m1)


# try a negative binomial if there is overdispersion

# Fit a Negative Binomial GLMM
rep_m_negbin <- glmer.nb(value ~ Time + (1 | MaleID), data = repro_long)
summary(rep_m_negbin)

# analyse if there is a significant change in TE by family/group
df$Time <- as.factor(df$Time)
levels(df$Time) <- list('1' = 'C', '2' = 'O') #, '3' ='outer') 1 is control grouup/ref group
levels(df$Time)


# running poisson model
rep_m1 <- glmer(Total ~ Time + (1 | MaleID), family = poisson, data = df)
summary(rep_m1)


# try a negative binomial if there is overdispersion

# Fit a Negative Binomial GLMM
rep_m_negbin_total <- glmer.nb(Total ~ Time + (1 | MaleID), data = df)
summary(rep_m_negbin_total)

rep_m_negbin_DNA <- glmer.nb(DNA ~ Time + (1 | MaleID), data = df)
summary(rep_m_negbin_DNA)

rep_m_negbin_LINE <- glmer.nb(LINE ~ Time + (1 | MaleID), data = df)
summary(rep_m_negbin_LINE)

rep_m_negbin_SINE <- glmer.nb(SINE ~ Time + (1 | MaleID), data = df)
summary(rep_m_negbin_SINE)

rep_m_negbin_LTR <- glmer.nb(LTR ~ Time + (1 | MaleID), data = df)
summary(rep_m_negbin_LTR)

rep_m_negbin_SATELLITTE <- glmer.nb(SATELLITTE ~ Time + (1 | MaleID), data = df)
summary(rep_m_negbin_SATELLITTE)


rep_m_negbin_RC <- glmer.nb(RC ~ Time + (1 | MaleID), data = df)
summary(rep_m_negbin_RC)
# fit poisson
# Change these Negative Binomial GLMMs to Poisson GLMMs

# Fit a Poisson GLMM for Total
rep_m_poisson_total <- glmer(Total ~ Time + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_total)

# Fit a Poisson GLMM for DNA
rep_m_poisson_DNA <- glmer(DNA ~ Time + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_DNA)

# Fit a Poisson GLMM for LINE
rep_m_poisson_LINE <- glmer(LINE ~ Time + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_LINE)

# Fit a Poisson GLMM for SINE
rep_m_poisson_SINE <- glmer(SINE ~ Time + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_SINE)

# Fit a Poisson GLMM for LTR
rep_m_poisson_LTR <- glmer(LTR ~ Time + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_LTR)

# Fit a Poisson GLMM for SATELLITTE
rep_m_poisson_SATELLITTE <- glmer(SATELLITTE ~ Time + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_SATELLITTE)


rep_m_poisson_RC <- glmer(RC ~ Time + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_RC)

# running to see central vs outer 2
# --- 1. Set 'Time2' as the reference level for the 'Time' factor ---
# Assuming your Time levels are '1', '2', '3'
df$Time <- relevel(df$Time, ref = "1")

#then re-run the above
