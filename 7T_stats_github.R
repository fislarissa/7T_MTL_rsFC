##########################################################################################################################

#Laminar 7T functional connectivity of the medial temporal lobe reflects patterns 
#of tau vulnerability and compensation in older adults

#Analysis of Z03/ B04 data - 7 Tesla resting-state fMRI with cognitive and AD markers
#Larissa Fischer - Multimodal Neuroimaging Lab, DZNE Magdeburg, 2025

#########################################################################################################################


#Import packages ################################################################
library(readxl)
library(writexl)
library(tidyr)
library(dplyr)
library(car)
library(lmtest)
library(ggplot2)
library(sjPlot)
library(sjmisc)
library(Hmisc)
library(car)
library(MuMIn)
library(sjstats)
library(ppcor)
library(effsize)
library(lubridate)
library(emmeans)
library(corrplot)
library(purrr)
library(broom)


#set working directory:
setwd("/Users/your/path")

#load data
load("/.../data.RData")
#...

#functions#######################################################################

plot_theme <- function() {
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )
}

summarise_group <- function(df, split_apoe = FALSE) {
  if(split_apoe){
    grp <- c("APOE4_carrier")
  } else {
    grp <- NULL
  }
  
  df %>%
    group_by(across(all_of(grp))) %>%
    summarise(
      n = n(),
      n_female = sum(sex == "f", na.rm = TRUE),
      n_APOE4 = sum(APOE4_carrier == "1", na.rm = TRUE),
      
      mean_age = mean(age, na.rm = TRUE),
      sd_age   = sd(age, na.rm = TRUE),
      min_age  = min(age, na.rm = TRUE),
      max_age  = max(age, na.rm = TRUE),
      
      mean_edu = mean(education, na.rm = TRUE),
      sd_edu   = sd(education, na.rm = TRUE),
      
      mean_AB_4240_ratio = mean(AB_4240_ratio, na.rm = TRUE),
      sd_AB_4240_ratio   = sd(AB_4240_ratio, na.rm = TRUE),
      mean_pTau217_fujirebio = mean(pTau217_fujirebio, na.rm = TRUE),
      sd_pTau217_fujirebio   = sd(pTau217_fujirebio, na.rm = TRUE),
      mean_AT_Term = mean(AT_Term, na.rm = TRUE),
      sd_AT_Term   = sd(AT_Term, na.rm = TRUE),
      mean_GFAP = mean(GFAP, na.rm = TRUE),
      sd_GFAP   = sd(GFAP, na.rm = TRUE),
      
      mean_AB_4240_ratio_tp2 = mean(AB_4240_ratio_tp2, na.rm = TRUE),
      sd_AB_4240_ratio_tp2   = sd(AB_4240_ratio_tp2, na.rm = TRUE),
      mean_pTau217_fujirebio_tp2 = mean(pTau217_fujirebio_tp2, na.rm = TRUE),
      sd_pTau217_fujirebio_tp2   = sd(pTau217_fujirebio_tp2, na.rm = TRUE),
      mean_AT_Term_tp2 = mean(AT_Term_tp2, na.rm = TRUE),
      sd_AT_Term_tp2   = sd(AT_Term_tp2, na.rm = TRUE),
      mean_GFAP_tp2 = mean(GFAP_tp2, na.rm = TRUE),
      sd_GFAP_tp2   = sd(GFAP_tp2, na.rm = TRUE),
      
      mean_metaROI = mean(metaROI_DVR, na.rm = TRUE),
      sd_metaROI   = sd(metaROI_DVR, na.rm = TRUE),
      
      mean_RSC = mean(RSC_DVR, na.rm = TRUE),
      sd_RSC   = sd(RSC_DVR, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_cognition <- function(df, split_apoe = FALSE) {
  if(split_apoe){
    grp <- c("APOE4_carrier")
  } else {
    grp <- NULL
  }
  
  df %>%
    group_by(across(all_of(grp))) %>%
    summarise(
      n = n(),
      
      mean_VLMT = mean(VLMT_recog, na.rm = TRUE),
      sd_VLMT   = sd(VLMT_recog, na.rm = TRUE),
      
      mean_RCFT = mean(RCFT_DR, na.rm = TRUE),
      sd_RCFT   = sd(RCFT_DR, na.rm = TRUE),
      
      mean_WMS = mean(WMS_DR, na.rm = TRUE),
      sd_WMS   = sd(WMS_DR, na.rm = TRUE),
      
      mean_VLMT_DR = mean(VLMT_DR, na.rm = TRUE),
      sd_VLMT_DR   = sd(VLMT_DR, na.rm = TRUE),
      
      mean_EM_composite = mean(EM_composite, na.rm = TRUE),
      sd_EM_composite   = sd(EM_composite, na.rm = TRUE),
      
      mean_SDMT = mean(SDMT, na.rm = TRUE),
      sd_SDMT   = sd(SDMT, na.rm = TRUE),
      
      mean_CAG = mean(CAG, na.rm = TRUE),
      sd_CAG   = sd(CAG, na.rm = TRUE),
      .groups = "drop"
    )
}



#demographics/ basics#############################################################
################################################################################


# without APOE split
demographics_table <- summarise_group(data_current, split_apoe = FALSE)
cog_table <- summarise_cognition(data_current, split_apoe = FALSE)
# with APOE split
demographics_table_APOE <- summarise_group(data_current, split_apoe = TRUE)
cog_table_APOE <- summarise_cognition(data_current, split_apoe = TRUE)


# time between measurements Z03
summary_time_diffs <- data_current %>% 
  summarise(
    across(
      .cols = c(diff_days_MRI_PET, diff_days_MRI_blood_I, diff_days_MRI_cog_I, diff_days_MRI_cog_II, 
                diff_days_cog_I_cog_II, diff_days_blood_I_blood_II),
      .fns = list(
        n   = ~sum(!is.na(.)),
        mean = ~mean(., na.rm = TRUE),
        sd   = ~sd(., na.rm = TRUE),
        min  = ~min(., na.rm = TRUE),
        max  = ~max(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )


#check hippocampal volume
lm_hc_vol <- lm(Sum_HC_bilat_TIVcorr ~ scale(metaROI_DVR)+ scale(AT_Term)+scale(GFAP)+APOE4_carrier+ scale(age) + sex +scale(education), data = data_current)
summary(lm_hc_vol)
shapiro.test(rstandard(lm_hc_vol)) 
bptest(lm_hc_vol) 
vif(lm_hc_vol)
tab_model(lm_hc_vol, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "lm_hc_vol.doc")

# check AD pathology
lm_apoe_tau_MTL <- lm(metaROI_DVR ~ APOE4_carrier + age + sex + education, data = data_current)
summary(lm_apoe_tau_MTL)
shapiro.test(rstandard(lm_apoe_tau_MTL)) 
bptest(lm_apoe_tau_MTL) 
vif(lm_apoe_tau_MTL)
tab_model(lm_apoe_tau_MTL, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "lm_apoe_tau_MTL.doc")

lm_apoe_tau_RSC <- lm(RSC_DVR ~ APOE4_carrier + age + sex + education, data = data_current)
summary(lm_apoe_tau_RSC)
shapiro.test(rstandard(lm_apoe_tau_RSC)) 
bptest(lm_apoe_tau_RSC) 
vif(lm_apoe_tau_RSC)
tab_model(lm_apoe_tau_RSC, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "lm_apoe_tau_RSC.doc")

lm_apoe_plasma <- lm(AT_Term ~ APOE4_carrier + age + sex + education, data = data_current)
summary(lm_apoe_plasma)
shapiro.test(rstandard(lm_apoe_plasma)) 
bptest(lm_apoe_plasma) 
vif(lm_apoe_plasma)
tab_model(lm_apoe_plasma, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "lm_apoe_plasma.doc")


# change in cognition
summary(data_current$EM_composite_change)
shapiro.test(data_current$EM_composite_change)
hist(data$EM_composite_change)
qqnorm(data$EM_composite_change); qqline(data$EM_composite_change)
t.test(data_current$EM_composite_change, mu = 0)

# change in AT-term
summary(data_current$AT_Term_change)
shapiro.test(data_current$AT_Term_change)
hist(data$AT_Term_change)
qqnorm(data$AT_Term_change); qqline(data$AT_Term_change)
t.test(data_current$AT_Term_change, mu = 0)

cor.test(data_current$AT_Term, data_current$AT_Term_tp2)
model <- lm(AT_Term_tp2 ~ AT_Term + age + APOE4_carrier + sex + education, data = data_current)
summary(model)
shapiro.test(rstandard(model)) 
bptest(model) 
vif(model)




#95% CI p value cluster
p <- 0.048 #change accordingly
M <- 10000  #permutationen
z <- 1.96
SE <- sqrt((p * (1 - p)) / M)
CI_lower <- p - z * SE
CI_upper <- p + z * SE
CI_lower <- round(CI_lower, 3)
CI_upper <- round(CI_upper, 3)
cat("95%-CI: [", CI_lower,",",CI_upper, "]\n")


#cohort age######################################################################
#FC cohort age###################################################################

#H1.1 age

 #1
plot1_1_1 <- ggplot(data_7T_rs_A_T_,aes(y=DG_LH_HB.BA36_LH,x=age))+
  geom_point()+  geom_smooth(method = "lm", color = "#0054ffff") +
  plot_theme()+labs(
    x = "Age in years",
    y = "rsFC left BA36 - left DG")
plot1_1_1


model1_1_1 <- lm(DG_LH_HB.BA36_LH ~ age+APOE4_carrier + sex + education, data = data_7T_rs_A_T_) #+ WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model1_1_1) 
shapiro.test(rstandard(model1_1_1)) 
bptest(model1_1_1) 
vif(model1_1_1)
tab_model(model1_1_1, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_1_1.doc")

#2
plot1_1_2 <- ggplot(data_7T_rs_A_T_,aes(y=SUB_RH_HB_deep.BA36_LH,x=age))+
  geom_point()+geom_smooth(method="lm")+plot_theme()+labs(
    x = "age",
    y = "FC right deep subiculum - left BA36")
plot1_1_2 
model1_1_2 <- lm(SUB_RH_HB_deep.BA36_LH ~ age+APOE4_carrier + sex + education, data = data_7T_rs_A_T_) #+ WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model1_1_2)
shapiro.test(rstandard(model1_1_2))
bptest(model1_1_2) 
vif(model1_1_2)
tab_model(model1_1_2, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_1_2.doc")


#graph age
plot2_1_1 <- ggplot(data_7T_rs_A_T_,aes(y=CA3_LH_HH_LocalEfficiency,x=age))+
  geom_point()+geom_smooth(method="lm")+plot_theme()+labs(
    x = "age",
    y = "left CA3 head Local Efficiency")
plot2_1_1 
model2_1_1 <- lm(CA3_LH_HH_LocalEfficiency ~ age+APOE4_carrier + sex + education, data = data_7T_rs_A_T_) #+ WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model2_1_1)
shapiro.test(rstandard(model2_1_1)) 
bptest(model2_1_1) 
vif(model2_1_1)
tab_model(model2_1_1, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_1.doc")

model2_1_2 <- lm(Network_LocalEfficiency ~ age+APOE4_carrier + sex + education, data = data_7T_rs_A_T_) #+ WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model2_1_2)
shapiro.test(rstandard(model2_1_2)) 
bptest(model2_1_2) 
vif(model2_1_2)
tab_model(model2_1_2, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_2.doc")



# apoe
# nothing significant


# age x apoe
#plot 
plot1_1_3 <- ggplot(data_7T_rs_A_T_,aes(y=CA3_LH_HB_deep.CA3_LH_HB_superficial,x=age, color=APOE4_carrier))+
  geom_point()+geom_smooth(method="lm")+plot_theme()+labs(
    x = "Age in years",
    y = "rsFC left CA3 superficial - deep",
    color = expression(italic("APOE4")~group))+
  scale_color_manual(values = c("0" = "#0054ffff", "1" = "#c21759ff"),
                     labels = c("0" = "Non-Carrier", "1" = "Carrier"))
plot1_1_3   


model1_1_3 <- lm(CA3_LH_HB_deep.CA3_LH_HB_superficial ~ age*APOE4_carrier + sex + education, data = data_7T_rs_A_T_) #+ WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model1_1_3) #sanity check: same estimate, T-value etc as in CONN
shapiro.test(rstandard(model1_1_3)) 
bptest(model1_1_3) 
vif(model1_1_3)   
tab_model(model1_1_3, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_1_3.doc")

#Cohens f squared
# full and reduced models
mdl_full <- lm(CA3_LH_HB_deep.CA3_LH_HB_superficial ~ age*APOE4_carrier + sex + education, data = data_7T_rs_A_T_)
mdl_reduced <- lm(CA3_LH_HB_deep.CA3_LH_HB_superficial ~ age+APOE4_carrier + sex + education, data = data_7T_rs_A_T_)
# R²
R2_full <- r.squaredGLMM(mdl_full)[1]   # marginal R²
R2_red  <- r.squaredGLMM(mdl_reduced)[1]
# Delta R²
deltaR2 <- R2_full - R2_red
# f²
f2 <- deltaR2 / (1 - R2_full)
R2_full; R2_red; deltaR2; f2

#subgroups APOE
data_apoe4_pos <- subset(data_7T_rs_A_T_, APOE4_carrier == 1)
data_apoe4_neg <- subset(data_7T_rs_A_T_, APOE4_carrier == 0)

model_apoe4_pos <- lm(CA3_LH_HB_deep.CA3_LH_HB_superficial ~ age + sex + education,data = data_apoe4_pos)
summary(model_apoe4_pos)
vif(model_apoe4_pos)
bptest(model_apoe4_pos)
shapiro.test(rstandard(model_apoe4_pos))
tab_model(model_apoe4_pos,df.method = "satterthwaite",show.stat = TRUE,show.se = TRUE,show.std = TRUE,file = "model_1_1_3_apoe4_carriers.doc")

model_apoe4_neg <- lm(CA3_LH_HB_deep.CA3_LH_HB_superficial ~ age + sex + education, data = data_apoe4_neg)
summary(model_apoe4_neg)
vif(model_apoe4_neg)
bptest(model_apoe4_neg)
shapiro.test(rstandard(model_apoe4_neg))
tab_model( model_apoe4_neg,df.method = "satterthwaite",show.stat = TRUE, show.se = TRUE,show.std = TRUE,file = "model_1_1_3_apoe4_noncarriers.doc")




#FC - cognition cohort age########################################################

#3.1.lower rsfc - lower performance
# EM COMPOSITE
# age

model3_1_1 <- lm(EM_composite ~ scale(DG_LH_HB.BA36_LH)+scale(age)+ APOE4_carrier + sex + scale(education), data = data_7T_rs_A_T_cog)
summary(model3_1_1) 
shapiro.test(rstandard(model3_1_1)) 
bptest(model3_1_1) 
vif(model3_1_1)
tab_model(model3_1_1, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model3_1_1.doc")

# age x apoe
plot3_1_2 <- ggplot(data_7T_rs_A_T_cog,aes(y=CA3_LH_HB_deep.CA3_LH_HB_superficial,x=EM_composite))+
  geom_point()+geom_smooth(method="lm", color = "#c21759ff")+plot_theme()+labs(
    x = "Episodic memory performance (scaled)",
    y = "rsFC left CA3 superficial - deep")
plot3_1_2 

model3_1_2 <- lm(EM_composite ~ scale(CA3_LH_HB_deep.CA3_LH_HB_superficial)+scale(age)+APOE4_carrier+ sex + scale(education), data = data_7T_rs_A_T_cog)
summary(model3_1_2) 
shapiro.test(rstandard(model3_1_2)) 
bptest(model3_1_2) 
vif(model3_1_2)
tab_model(model3_1_2, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model3_1_2_with_apoe.doc")


#FC - longit cognition cohort age###############################################

#age

data_plot <- data_7T_rs_A_T_cog %>%
  mutate(age_group = ifelse(age <= 70, "≤70", ">70"))
plot3_1_1_change <- 
  ggplot(data_plot, aes(x = DG_LH_HB.BA36_LH, y = EM_composite_change, color = age_group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(
    x = "DG–BA36 FC",
    y = "EM composite change",
    color = "Age group")
plot3_1_1_change

model3_1_1_change <- lm(EM_composite_change ~ scale(DG_LH_HB.BA36_LH)*scale(age)+APOE4_carrier+ sex + scale(education), data = data_7T_rs_A_T_cog)
summary(model3_1_1_change) #+diff_days_cog_I_cog_II no effect
shapiro.test(rstandard(model3_1_1_change)) 
bptest(model3_1_1_change) 
vif(model3_1_1_change)
tab_model(model3_1_1_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model3_1_1_change.doc")

#age x APOE
model3_1_2_change <- lm(EM_composite_change ~ scale(CA3_LH_HB_deep.CA3_LH_HB_superficial)+scale(age)+ APOE4_carrier+ sex + scale(education), data = data_7T_rs_A_T_cog)
summary(model3_1_2_change) #+diff_days_cog_I_cog_II no effect
shapiro.test(rstandard(model3_1_2_change)) 
bptest(model3_1_2_change) 
vif(model3_1_2_change)


#cohort blood AT term############################################################
#FC cohort blood AT term##########################################################

#H1.2

plot1_2_1 <- ggplot(data_7T_rs_blood,aes(y=CA1_LH_HH.BA36_RH,x=AT_Term))+
  geom_point()+geom_smooth(method="lm", color = "#c21759ff")+plot_theme()+labs(
    x = "AT-Term at baseline",
    y = "rsFC left CA1 head - right BA36")
plot1_2_1 

model1_2_1 <- lm(CA1_LH_HH.BA36_RH ~ AT_Term+age + sex + education, data = data_7T_rs_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model1_2_1)
shapiro.test(rstandard(model1_2_1)) 
bptest(model1_2_1) 
vif(model1_2_1)
tab_model(model1_2_1, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_2_1.doc")



#FC cohort longit blood AT term##########################################################

#take connection from above:

#CHANGE
#H1.2
plot1_2_1_change <- ggplot(data_7T_rs_blood,aes(y=CA1_LH_HH.BA36_RH,x=AT_Term_change))+
  geom_point()+geom_smooth(method="lm")+plot_theme()+labs(
    x = "AT-Term change",
    y = "rsFC left CA1 head - right BA36")
plot1_2_1_change
model1_2_1_change <- lm(CA1_LH_HH.BA36_RH ~ AT_Term_change+age + sex + education, data = data_7T_rs_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model1_2_1_change)
shapiro.test(rstandard(model1_2_1_change)) 
bptest(model1_2_1_change) 
vif(model1_2_1_change)
tab_model(model1_2_1_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_2_1_change.doc")


#graph:

plot1_2_6_change <- ggplot(data_7T_rs_blood,aes(y=Network_GlobalEfficiency,x=AT_Term_change))+
  geom_point()+geom_smooth(method="lm", color = "#0054ffff")+plot_theme()+labs(
    x = "Change in AT-Term",
    y = "Network global efficiency")
plot1_2_6_change

model2_1_6_change <- lm(Network_GlobalEfficiency ~ AT_Term_change+age + sex + education, data = data_7T_rs_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model2_1_6_change)
shapiro.test(rstandard(model2_1_6_change)) 
bptest(model2_1_6_change) 
vif(model2_1_6_change)
tab_model(model2_1_6_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_6_change.doc")

model2_1_7_change <- lm(Tail_LH_GlobalEfficiency ~ AT_Term_change+age + sex + education, data = data_7T_rs_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model2_1_7_change)
shapiro.test(rstandard(model2_1_7_change)) 
bptest(model2_1_7_change) 
vif(model2_1_7_change)
tab_model(model2_1_7_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_7_change.doc")

model2_1_8_change <- lm(DG_LH_HB_GlobalEfficiency ~ AT_Term_change+age + sex + education, data = data_7T_rs_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model2_1_8_change)
shapiro.test(rstandard(model2_1_8_change)) 
bptest(model2_1_8_change) 
vif(model2_1_8_change)
tab_model(model2_1_8_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_8_change.doc")

model2_1_9_change <- lm(CA1_LH_HB_deep_GlobalEfficiency ~ AT_Term_change+age + sex + education, data = data_7T_rs_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model2_1_9_change)
shapiro.test(rstandard(model2_1_9_change)) 
bptest(model2_1_9_change) 
vif(model2_1_9_change)
tab_model(model2_1_9_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_9_change.doc")

model2_1_10_change <- lm(BA35_LH_deep_GlobalEfficiency ~ AT_Term_change+age + sex + education, data = data_7T_rs_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model2_1_10_change)
shapiro.test(rstandard(model2_1_10_change)) 
bptest(model2_1_10_change) 
vif(model2_1_10_change)
tab_model(model2_1_10_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_10_change.doc")

model2_1_11_change <- lm(BA35_RH_superficial_AveragePathLength ~ AT_Term_change+age + sex + education, data = data_7T_rs_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model2_1_11_change)
shapiro.test(rstandard(model2_1_11_change)) 
bptest(model2_1_11_change) 
vif(model2_1_11_change)
tab_model(model2_1_11_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_11_change.doc")



#FC - cognition cohort blood######################################################

#H3.2
#MTL
#base 

plot3_2_0_base <- ggplot(data_7T_rs_blood_cog,aes(x=AT_Term,y=EM_composite))+
  geom_point()+geom_smooth(method="lm")
plot3_2_0_base
model3_2_0_base <- lm(EM_composite ~ scale(AT_Term)+ CA1_LH_HH.BA36_RH+scale(age)+ sex + scale(education), data = data_7T_rs_blood_cog)
summary(model3_2_0_base)
shapiro.test(rstandard(model3_2_0_base)) 
bptest(model3_2_0_base) 
vif(model3_2_0_base)
tab_model(model3_2_0_base, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model3_2_0_base.doc")
#no base, no weak evidence
# no interaction FC * APOE or GFAP


#CAG
model3_2_0_base_2 <- lm(CAG ~ scale(AT_Term)+ sex + scale(education), data = data_7T_rs_blood_cog)
summary(model3_2_0_base_2)  #no base, no weak evidence 
shapiro.test(rstandard(model3_2_0_base_2)) 
bptest(model3_2_0_base_2) 
vif(model3_2_0_base_2)

#SDMT
model3_2_0_base_3 <- lm(SDMT ~ scale(AT_Term)+ scale(age)+ sex + scale(education), data = data_7T_rs_blood_cog)
summary(model3_2_0_base_3)  #no base, no weak evidence
shapiro.test(rstandard(model3_2_0_base_3)) 
bptest(model3_2_0_base_3) 
vif(model3_2_0_base_3)

#FC - longit cognition cohort blood######################################################


model3_2_0_change_base <- lm(EM_composite_change ~ scale(AT_Term)+ scale(age)+ sex + scale(education), data = data_7T_rs_blood_cog)
summary(model3_2_0_change_base)
shapiro.test(rstandard(model3_2_0_change_base)) 
bptest(model3_2_0_change_base) 
vif(model3_2_0_change_base)
# no base, no weak evidence

model3_2_0_change_change <- lm(EM_composite_change ~ scale(AT_Term_change)+ scale(age)+ sex + scale(education), data = data_7T_rs_blood_cog)
summary(model3_2_0_change_change)
shapiro.test(rstandard(model3_2_0_change_change)) 
bptest(model3_2_0_change_change) 
vif(model3_2_0_change_change)
# no base, no weak evidence


#cohort PET MTL###################################################################
#FC cohort PET MTL################################################################

#H1.3

#tau
#graph
model2_1_3 <- lm(ERC_RH_superficial_LocalEfficiency ~ metaROI_DVR + age + sex + education, data = data_7T_rs_PET_MTL)
summary(model2_1_3)
shapiro.test(rstandard(model2_1_3)) 
bptest(model2_1_3)
vif(model2_1_3)
tab_model(model2_1_3, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_3.doc")


#tau x GFAP

#strongest connection
#median split
data_1_3_1 <- data_7T_rs_PET_MTL %>%
  mutate(
    GFAP_z = scale(GFAP)[, 1],
    GFAP_group = ifelse(
      GFAP_z <= median(GFAP_z, na.rm = TRUE),
      "Low GFAP",
      "High GFAP"))
plot1_3_1 <- ggplot(
  data_1_3_1,
  aes(x = metaROI_DVR, y = SUB_LH_HH.CA1_LH_HH, colour = GFAP_group)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  labs(
    x = "Medial temporal lobe tau burden",
    y = "rsFC left CA1 head - left subiculum head",
    colour = "GFAP:") +
  scale_colour_manual(
    values = c("Low GFAP" = "#0054ffff", "High GFAP" = "#c21759ff")
  ) +plot_theme() + theme(legend.position = "top")
plot1_3_1

model1_3_1 <- lm(SUB_LH_HH.CA1_LH_HH ~  metaROI_DVR*GFAP + age + sex + education, data = data_7T_rs_PET_MTL) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model1_3_1) 
shapiro.test(rstandard(model1_3_1)) 
bptest(model1_3_1) 
vif(model1_3_1) 
tab_model(model1_3_1, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_3_1.doc")

#low and high gfap
slopes <- emtrends(model1_3_1, ~ GFAP, var = "metaROI_DVR",
                   at = list(GFAP = quantile(data_7T_rs_PET_MTL$GFAP,
                                             probs = c(.05, .5, .95),
                                             na.rm = TRUE)))
summary(slopes, infer = TRUE)



#other
model1_3_2 <- lm(CA1_LH_HH.BA36_LH ~ metaROI_DVR*GFAP + age + sex + education, data = data_7T_rs_PET_MTL)
summary(model1_3_2) 
shapiro.test(rstandard(model1_3_2))  
bptest(model1_3_2)
vif(model1_3_2)
tab_model(model1_3_2, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_3_2.doc")

model1_3_3 <- lm(PHC_LH_superficial.ERC_LH_superficial ~ metaROI_DVR*GFAP + age + sex + education, data = data_7T_rs_PET_MTL)
summary(model1_3_3) 
shapiro.test(rstandard(model1_3_3)) 
bptest(model1_3_3)
vif(model1_3_3)
tab_model(model1_3_3, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_3_3.doc")

model1_3_4 <- lm(SUB_LH_HH.Tail_LH ~ metaROI_DVR*GFAP + age + sex + education, data = data_7T_rs_PET_MTL)
summary(model1_3_4) 
shapiro.test(rstandard(model1_3_4)) 
bptest(model1_3_4)
vif(model1_3_4)
tab_model(model1_3_4, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_3_4.doc")

model1_3_5 <- lm(ERC_LH_superficial.BA36_LH ~ metaROI_DVR*GFAP + age + sex + education, data = data_7T_rs_PET_MTL)
summary(model1_3_5) 
shapiro.test(rstandard(model1_3_5)) 
bptest(model1_3_5)
vif(model1_3_5)
tab_model(model1_3_5, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_3_5.doc")



#FC - cognition cohort PET MTL####################################################

#H3.2: 
#If higher rsFC is related to cognitive reserve processes, rsFC will moderate 
#the relationship between pathology burden and memory performance such that in 
#older adults with higher rsFC strength, the association of higher pathology 
#and lower memory performance will be weaker.

#MTL
#base

plot3_2_1_base <- ggplot(data_7T_rs_PET_MTL_cog,aes(x=metaROI_DVR,y=EM_composite))+
  geom_point()+geom_smooth(method="lm")
plot3_2_1_base
model3_2_1_base <- lm(EM_composite ~ scale(metaROI_DVR)+ scale(age)+ sex + scale(education), data = data_7T_rs_PET_MTL_cog)
summary(model3_2_1_base)
tab_model(model3_2_1_base, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model3_2_1_base.doc")

#moderation
#plot median split
data_plot <- data_7T_rs_PET_MTL_cog %>% 
  mutate(
    rsFC_z = scale(SUB_LH_HH.CA1_LH_HH)[,1],
    rsFC_group = ifelse(rsFC_z <= median(rsFC_z), "Low rsFC", "High rsFC"))
plot3_2_1 <- ggplot(data_plot, aes(x = metaROI_DVR, y = EM_composite, colour = rsFC_group)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  labs(
    x = "Medial temporal lobe tau burden",
    y = "Episodic memory performance (scaled)",
    colour = "Left CA1 head - left subiculum head:"
  ) +
  scale_colour_manual(values = c("Low rsFC" = "#0054ffff", "High rsFC" = "#c21759ff")) +
  plot_theme() +
  theme(legend.position = "top")
plot3_2_1


model3_2_1 <- lm(EM_composite ~ scale(SUB_LH_HH.CA1_LH_HH)*scale(metaROI_DVR)+ scale(age) + sex + scale(education), data = data_7T_rs_PET_MTL_cog)
summary(model3_2_1)#+WM_DVR #+Sum_HC_bilat_TIVcorr
shapiro.test(rstandard(model3_2_1)) 
bptest(model3_2_1) 
vif(model3_2_1)
tab_model(model3_2_1, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model3_2_1.doc")


#Cohens f squared
# full and reduced models
mdl_full <- lm(EM_composite ~ scale(SUB_LH_HH.CA1_LH_HH)*scale(metaROI_DVR)+ scale(age) + sex + scale(education), data = data_7T_rs_PET_MTL_cog)
mdl_reduced <- lm(EM_composite ~ scale(SUB_LH_HH.CA1_LH_HH)+scale(metaROI_DVR)+ scale(age) + sex + scale(education), data = data_7T_rs_PET_MTL_cog)
# R²
R2_full <- r.squaredGLMM(mdl_full)[1]   # marginal R²
R2_red  <- r.squaredGLMM(mdl_reduced)[1]
# Delta R²
deltaR2 <- R2_full - R2_red
# f²
f2 <- deltaR2 / (1 - R2_full)
R2_full; R2_red; deltaR2; f2


#FC subgroups (median split from plot above)

model_low <- lm(
  EM_composite ~ scale(metaROI_DVR) + scale(age) + sex + scale(education),
  data = data_plot %>% filter(rsFC_group == "Low FC"))
summary(model_low)
tab_model(model_low, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model_3_2_1_low.doc")


model_high <- lm(
  EM_composite ~ scale(metaROI_DVR) + scale(age) + sex + scale(education),
  data = data_plot %>% filter(rsFC_group == "High FC"))
summary(model_high)
tab_model(model_high, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model_3_2_1_high.doc")



#CAG
model3_2_1_base_2 <- lm(CAG ~ scale(metaROI_DVR)+ sex + scale(education), data = data_7T_rs_PET_MTL_cog)
summary(model3_2_1_base_2) #no base, no weak evidence
shapiro.test(rstandard(model3_2_1_base_2)) 
bptest(model3_2_1_base_2)
vif(model3_2_1_base_2)

#SDMT
model3_2_1_base_3 <- lm(SDMT ~ scale(metaROI_DVR)+ scale(age)+ sex + scale(education), data = data_7T_rs_PET_MTL_cog)
summary(model3_2_1_base_3) #no base, no weak evidence
shapiro.test(rstandard(model3_2_1_base_3)) 
bptest(model3_2_1_base_3)
vif(model3_2_1_base_3)


#FC - longit cognition cohort PET MTL####################################################

#base
plot3_2_1_change_base <- ggplot(data_7T_rs_PET_MTL_cog,aes(x=metaROI_DVR,y=EM_composite_change))+
  geom_point()+geom_smooth(method="lm")
plot3_2_1_change_base
model3_2_1_change_base <- lm(EM_composite_change ~ scale(metaROI_DVR)+ scale(age)+ sex + scale(education), data = data_7T_rs_PET_MTL_cog)
summary(model3_2_1_change_base)
shapiro.test(rstandard(model3_2_1_change_base)) 
bptest(model3_2_1_change_base)
vif(model3_2_1_change_base)

#FC as predictor
plot3_2_1_change <- ggplot(data_7T_rs_PET_MTL_cog,aes(x=SUB_LH_HH.CA1_LH_HH,y=EM_composite_change))+
  geom_point()+geom_smooth(method="lm")+plot_theme()+
  labs(
    x = "rsFC left CA1 head - left subiculum head",
    y = "Change in episodic memory performance (scaled)")
plot3_2_1_change

model3_2_1_change <- lm(EM_composite_change ~scale(metaROI_DVR)+ scale(SUB_LH_HH.CA1_LH_HH)+APOE4_carrier+ scale(age) + sex + scale(education), data = data_7T_rs_PET_MTL_cog)
summary(model3_2_1_change)
shapiro.test(rstandard(model3_2_1_change)) 
bptest(model3_2_1_change)
vif(model3_2_1_change)
tab_model(model3_2_1_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model3_2_1_change.doc")



#cohort PET RSC###################################################################
#FC cohort PET RSC################################################################

#H1.4

#tau
#strongest connection
plot1_4_1 <- ggplot(data_7T_rs_PET_RSC,aes(y=RSC_superficial_LH.Tail_RH,x=RSC_DVR))+
  geom_point()+geom_smooth(method="lm", color = "#c21759ff")+plot_theme()+
  labs(
    x = "Retrosplenial cortex tau burden",
    y = "rsFC right CA1 tail - left superficial RSC")
plot1_4_1


model1_4_1 <- lm(RSC_superficial_LH.Tail_RH ~ RSC_DVR + age + sex + education, data = data_7T_rs_PET_RSC)
summary(model1_4_1) ##+WM_DVR #+Sum_HC_bilat_TIVcorr
shapiro.test(rstandard(model1_4_1)) 
bptest(model1_4_1) 
vif(model1_4_1)  
tab_model(model1_4_1, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_4_1.doc")

#other
plot1_4_2 <- ggplot(data_7T_rs_PET_RSC,aes(y=SUB_LH_HB_deep.RSC_deep_LH,x=RSC_DVR))+
  geom_point()+geom_smooth(method="lm")+plot_theme()
plot1_4_2
model1_4_2 <- lm(SUB_LH_HB_deep.RSC_deep_LH ~ RSC_DVR + age + sex + education, data = data_7T_rs_PET_RSC)
summary(model1_4_2) 
shapiro.test(rstandard(model1_4_2)) 
bptest(model1_4_2)
vif(model1_4_2)
tab_model(model1_4_2, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_4_2.doc")


plot1_4_3 <- ggplot(data_7T_rs_PET_RSC,aes(y=RSC_deep_LH.Tail_RH,x=RSC_DVR))+
  geom_point()+geom_smooth(method="lm")+plot_theme()
plot1_4_3
model1_4_3 <- lm(RSC_deep_LH.Tail_RH ~ RSC_DVR + age + sex + education, data = data_7T_rs_PET_RSC)
summary(model1_4_3) 
shapiro.test(rstandard(model1_4_3)) 
bptest(model1_4_3)
vif(model1_4_3)
tab_model(model1_4_3, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model1_4_3.doc")



# graph
plot2_1_4 <- ggplot(data_7T_rs_PET_RSC,aes(y=RSC_superficial_LH_GlobalEfficiency,x=RSC_DVR))+
  geom_point()+geom_smooth(method="lm")+plot_theme()
plot2_1_4
model2_1_4 <- lm(RSC_superficial_LH_GlobalEfficiency ~ RSC_DVR + age + sex + education, data = data_7T_rs_PET_RSC)
summary(model2_1_4) ##+WM_DVR #+Sum_HC_bilat_TIVcorr
shapiro.test(rstandard(model2_1_4)) 
bptest(model2_1_4) 
vif(model2_1_4)  
tab_model(model2_1_4, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_4.doc")

plot2_1_5 <- ggplot(data_7T_rs_PET_RSC,aes(y=RSC_superficial_LH_AveragePathLength,x=RSC_DVR))+
  geom_point()+geom_smooth(method="lm")+plot_theme()
plot2_1_5
model2_1_5 <- lm(RSC_superficial_LH_AveragePathLength ~ RSC_DVR + age + sex + education, data = data_7T_rs_PET_RSC)
summary(model2_1_5) ##+WM_DVR #+Sum_HC_bilat_TIVcorr
shapiro.test(rstandard(model2_1_5)) 
bptest(model2_1_5) 
vif(model2_1_5)  
tab_model(model2_1_5, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model2_1_5.doc")



#FC - cognition cohort PET RSC####################################################

#H3.2: 
#If higher rsFC is related to cognitive reserve processes, rsFC will moderate 
#the relationship between pathology burden and memory performance such that in 
#older adults with higher rsFC strength, the association of higher pathology 
#and lower memory performance will be weaker.

#RSC
#median split
data_plot <- data_7T_rs_PET_RSC %>%
  mutate(
    rsFC_z = scale(RSC_superficial_LH.Tail_RH)[,1], 
    rsFC_group = ifelse(rsFC_z <= median(rsFC_z), "Low rsFC", "High rsFC"))
plot3_2_2 <- ggplot(data_plot, aes(x = RSC_DVR, y = EM_composite, colour = rsFC_group)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  labs(
    x = "RSC Tau DVR",
    y = "EM composite",
    colour = "rsFC (CA1 (Tail) – Superficial RSC)"
  ) +
  scale_colour_manual(values = c("Low rsFC" = "#1f78b4", "High rsFC" = "#e31a1c")) +
  plot_theme() +
  theme(legend.position = "top")
plot3_2_2

#base
model3_2_2_base <- lm(EM_composite ~ scale(RSC_DVR)+ scale(age)+ sex + scale(education), data = data_7T_rs_PET_RSC_cog)
summary(model3_2_2_base)
shapiro.test(rstandard(model3_2_2_base)) 
bptest(model3_2_2_base) 
vif(model3_2_2_base)  
tab_model(model3_2_2_base, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model3_2_2_base.doc")
#no base, no weak evidence

#CAG
model3_2_2_base_2 <- lm(CAG ~ scale(RSC_DVR)+ sex + scale(education), data = data_7T_rs_PET_RSC_cog)
summary(model3_2_2_base_2) #no base, no weak evidence
shapiro.test(rstandard(model3_2_2_base_2)) 
bptest(model3_2_2_base_2) 
vif(model3_2_2_base_2)  

#SDMT
model3_2_2_base_3 <- lm(SDMT ~ scale(RSC_DVR)+ scale(age)+ sex + scale(education), data = data_7T_rs_PET_RSC_cog)
summary(model3_2_2_base_3) #no base, no weak evidence
shapiro.test(rstandard(model3_2_2_base_3)) 
bptest(model3_2_2_base_3) 
vif(model3_2_2_base_3)  


#FC - longit cognition cohort PET RSC####################################################

model3_2_2_change <- lm(EM_composite_change ~ scale(RSC_DVR)+ scale(age)+ sex + scale(education), data = data_7T_rs_PET_RSC_cog)
summary(model3_2_2_change)
shapiro.test(rstandard(model3_2_2_change)) 
bptest(model3_2_2_change) 
vif(model3_2_2_change)  


#################################################################################
#additional analyses###############################################################

#additional results conn cvr###############################################################

#AT term change
model_cvr_1_change <- lm(SUB_RH_HB_superficial.CA1_LH_HH ~ AT_Term_change+age + sex + education, data = data_7T_cvr_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model_cvr_1_change)
shapiro.test(rstandard(model_cvr_1_change)) 
bptest(model_cvr_1_change) 
vif(model_cvr_1_change) 
tab_model(model_cvr_1_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model_cvr_1_change.doc")

model_cvr_2_change <- lm(DG_RH_HB.CA1_LH_HH ~ AT_Term_change+age + sex + education, data = data_7T_cvr_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model_cvr_2_change)
shapiro.test(rstandard(model_cvr_2_change)) 
bptest(model_cvr_2_change) 
vif(model_cvr_2_change) 
tab_model(model_cvr_2_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model_cvr_2_change.doc")

model_cvr_3_change <- lm(DG_RH_HH.DG_RH_HB ~ AT_Term_change+age + sex + education, data = data_7T_cvr_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model_cvr_3_change)
shapiro.test(rstandard(model_cvr_3_change)) 
bptest(model_cvr_3_change) 
vif(model_cvr_3_change) 
tab_model(model_cvr_3_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model_cvr_3_change.doc")

#graph
model_cvr_graph_change <- lm(RSC_superficial_LH_GlobalEfficiency ~ AT_Term_change+age + sex + education, data = data_7T_cvr_blood) #+WM_DVR #+Sum_HC_bilat_TIVcorr
summary(model_cvr_graph_change)
shapiro.test(rstandard(model_cvr_graph_change)) 
bptest(model_cvr_graph_change) 
vif(model_cvr_graph_change) 
tab_model(model_cvr_graph_change, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model_cvr_7_graph_change.doc")



#MTL tau*GFAP
model_cvr_4 <- lm(SUB_LH_HB_deep.Tail_LH ~ metaROI_DVR*GFAP + age + sex + education, data = data_7T_cvr_PET_MTL) 
summary(model_cvr_4)
shapiro.test(rstandard(model_cvr_4)) 
bptest(model_cvr_4) 
vif(model_cvr_4) 
tab_model(model_cvr_4, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model_cvr_4.doc")

model_cvr_5 <- lm(SUB_LH_HB_superficial.Tail_LH ~ metaROI_DVR*GFAP + age + sex + education, data = data_7T_cvr_PET_MTL)
summary(model_cvr_5)  
shapiro.test(rstandard(model_cvr_5)) 
bptest(model_cvr_5) 
vif(model_cvr_5) 
tab_model(model_cvr_5, df.method = "satterthwaite", show.stat = TRUE, show.se = TRUE, show.std = TRUE, file = "model_cvr_5.doc")



#sanity checks FC################################################################

# more FC within HC than ERC-HC, more ipsilateral thann contralateral
# rs and cvr separately (use respective data frame)

#specific sup ERC - sup CA3
fc_ERC_CA3 <- data_7T_rs_A_T_ %>%
  dplyr::select(
    CA3_LH_HB_superficial.ERC_LH_superficial,
    CA3_RH_HB_superficial.ERC_RH_superficial,
    CA3_RH_HB_superficial.ERC_LH_superficial,
    CA3_LH_HB_superficial.ERC_RH_superficial
  ) %>%
  rename(
    ipsi_LH = CA3_LH_HB_superficial.ERC_LH_superficial,
    ipsi_RH = CA3_RH_HB_superficial.ERC_RH_superficial,
    contra_LH_RH = CA3_RH_HB_superficial.ERC_LH_superficial,
    contra_RH_LH = CA3_LH_HB_superficial.ERC_RH_superficial
  ) %>%
  mutate(
    ERC_CA3_ipsi    = rowMeans(cbind(ipsi_LH, ipsi_RH), na.rm = TRUE),
    ERC_CA3_contra  = rowMeans(cbind(contra_LH_RH, contra_RH_LH), na.rm = TRUE))
summary(fc_ERC_CA3)


#specific sup ERC - sup CA1
fc_ERC_CA1 <- data_7T_rs_A_T_ %>%
  dplyr::select(
    CA1_LH_HB_superficial.ERC_LH_superficial,
    CA1_RH_HB_superficial.ERC_RH_superficial,
    CA1_RH_HB_superficial.ERC_LH_superficial,
    CA1_LH_HB_superficial.ERC_RH_superficial
  ) %>%
  rename(
    ipsi_LH = CA1_LH_HB_superficial.ERC_LH_superficial,
    ipsi_RH = CA1_RH_HB_superficial.ERC_RH_superficial,
    contra_LH_RH = CA1_RH_HB_superficial.ERC_LH_superficial,
    contra_RH_LH = CA1_LH_HB_superficial.ERC_RH_superficial
  ) %>%
  mutate(
    ERC_CA1_ipsi   = rowMeans(cbind(ipsi_LH, ipsi_RH), na.rm = TRUE),
    ERC_CA1_contra = rowMeans(cbind(contra_LH_RH, contra_RH_LH), na.rm = TRUE))
summary(fc_ERC_CA1)


#specific sup ERC - sup SUB
fc_ERC_SUB <- data_7T_rs_A_T_ %>%
  dplyr::select(
    SUB_LH_HB_superficial.ERC_LH_superficial,
    SUB_RH_HB_superficial.ERC_RH_superficial,
    SUB_RH_HB_superficial.ERC_LH_superficial,
    SUB_LH_HB_superficial.ERC_RH_superficial
  ) %>%
  rename(
    ipsi_LH = SUB_LH_HB_superficial.ERC_LH_superficial,
    ipsi_RH = SUB_RH_HB_superficial.ERC_RH_superficial,
    contra_LH_RH = SUB_RH_HB_superficial.ERC_LH_superficial,
    contra_RH_LH = SUB_LH_HB_superficial.ERC_RH_superficial
  ) %>%
  mutate(
    ERC_SUB_ipsi    = rowMeans(cbind(ipsi_LH, ipsi_RH), na.rm = TRUE),
    ERC_SUB_contra  = rowMeans(cbind(contra_LH_RH, contra_RH_LH), na.rm = TRUE))
summary(fc_ERC_SUB)


# within HC: ca3 - ca1 
fc_CA3_CA1 <- data_7T_rs_A_T_ %>%
  dplyr::select(
    CA3_LH_HB_deep.CA1_LH_HB_deep,
    CA3_RH_HB_deep.CA1_RH_HB_deep,
    CA3_LH_HB_deep.CA1_RH_HB_deep,
    CA3_RH_HB_deep.CA1_LH_HB_deep
  ) %>%
  rename(
    ipsi_LH = CA3_LH_HB_deep.CA1_LH_HB_deep,
    ipsi_RH = CA3_RH_HB_deep.CA1_RH_HB_deep,
    contra_LH_RH = CA3_LH_HB_deep.CA1_RH_HB_deep,
    contra_RH_LH = CA3_RH_HB_deep.CA1_LH_HB_deep
  ) %>%
  mutate(
    CA3_CA1_ipsi   = rowMeans(cbind(ipsi_LH, ipsi_RH), na.rm = TRUE),
    CA3_CA1_contra = rowMeans(cbind(contra_LH_RH, contra_RH_LH), na.rm = TRUE))
summary(fc_CA3_CA1)

#within HC: ca1 - sub
fc_CA1_SUB <- data_7T_rs_A_T_ %>%
  dplyr::select(
    SUB_LH_HB_deep.CA1_LH_HB_deep,
    SUB_RH_HB_deep.CA1_RH_HB_deep,
    SUB_LH_HB_deep.CA1_RH_HB_deep,
    SUB_RH_HB_deep.CA1_LH_HB_deep
  ) %>%
  rename(
    ipsi_LH = SUB_LH_HB_deep.CA1_LH_HB_deep,
    ipsi_RH = SUB_RH_HB_deep.CA1_RH_HB_deep,
    contra_LH_RH = SUB_LH_HB_deep.CA1_RH_HB_deep,
    contra_RH_LH = SUB_RH_HB_deep.CA1_LH_HB_deep
  ) %>%
  mutate(
    CA1_SUB_ipsi   = rowMeans(cbind(ipsi_LH, ipsi_RH), na.rm = TRUE),
    CA1_SUB_contra = rowMeans(cbind(contra_LH_RH, contra_RH_LH), na.rm = TRUE))
summary(fc_CA1_SUB)



# plot mean for ipsi and contra
fc_mean <- tibble(
  connection = c("ERC_CA3", "ERC_CA1", "ERC_SUB", "CA3_CA1", "CA1_SUB"),
  ipsi_mean   = c(
    mean(fc_ERC_CA3$ERC_CA3_ipsi, na.rm = TRUE),
    mean(fc_ERC_CA1$ERC_CA1_ipsi, na.rm = TRUE),
    mean(fc_ERC_SUB$ERC_SUB_ipsi, na.rm = TRUE),
    mean(fc_CA3_CA1$CA3_CA1_ipsi, na.rm = TRUE),
    mean(fc_CA1_SUB$CA1_SUB_ipsi, na.rm = TRUE)
  ),
  contra_mean = c(
    mean(fc_ERC_CA3$ERC_CA3_contra, na.rm = TRUE),
    mean(fc_ERC_CA1$ERC_CA1_contra, na.rm = TRUE),
    mean(fc_ERC_SUB$ERC_SUB_contra, na.rm = TRUE),
    mean(fc_CA3_CA1$CA3_CA1_contra, na.rm = TRUE),
    mean(fc_CA1_SUB$CA1_SUB_contra, na.rm = TRUE)))

# pivot longer for plotting
fc_plot <- fc_mean %>%
  pivot_longer(cols = c(ipsi_mean, contra_mean),
               names_to = "hemisphere",
               values_to = "mean_FC")
fc_plot$connection <- factor(
  fc_plot$connection,
  levels = c("ERC_CA3", "ERC_CA1", "ERC_SUB", "CA3_CA1", "CA1_SUB"))


# plot
ggplot(fc_plot, aes(x = connection, y = hemisphere, fill = mean_FC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_FC, 2)), color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels = c(
    CA1_SUB = "CA1 - subiculum rsFC",
    CA3_CA1 = "CA3 - CA1 rsFC",
    ERC_CA1 = "ERC - CA1 rsFC",
    ERC_CA3 = "ERC - CA3 rsFC",
    ERC_SUB = "ERC - subiculum rsFC"
  )) +
  scale_y_discrete(labels = c(
    ipsi_mean = "FC ipsilateral",
    contra_mean = "FC contralateral"
  )) +
  theme_classic() +
  labs(x = NULL, y = NULL, fill = "mean rsFC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# plot all ipsy boxplots
fc_long <- bind_cols(
  dplyr::select(fc_ERC_CA3, ERC_CA3_ipsi),
  dplyr::select(fc_ERC_CA1, ERC_CA1_ipsi),
  dplyr::select(fc_ERC_SUB, ERC_SUB_ipsi),
  dplyr::select(fc_CA3_CA1, CA3_CA1_ipsi),
  dplyr::select(fc_CA1_SUB, CA1_SUB_ipsi)
) %>%
  pivot_longer(
    cols = dplyr::everything(),
    names_to = "connection",
    values_to = "r")

ggplot(fc_long, aes(x = connection, y = r)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Functional connectivity (r)")



#for all
# put all datasets into a named list
fc_list <- list(
  ERC_CA3 = fc_ERC_CA3,
  ERC_CA1 = fc_ERC_CA1,
  ERC_SUB = fc_ERC_SUB,
  CA3_CA1 = fc_CA3_CA1,
  CA1_SUB = fc_CA1_SUB)

results <- map_df(names(fc_list), function(name) {
  df <- fc_list[[name]] %>%
    mutate(
      ipsi_mean   = rowMeans(cbind(ipsi_LH, ipsi_RH), na.rm = TRUE),
      contra_mean = rowMeans(cbind(contra_LH_RH, contra_RH_LH), na.rm = TRUE))
  test <- t.test(df$ipsi_mean, df$contra_mean, paired = TRUE)
  tibble(
    connection = name,
    mean_ipsi = mean(df$ipsi_mean, na.rm = TRUE),
    mean_contra = mean(df$contra_mean, na.rm = TRUE),
    median_ipsi = median(df$ipsi_mean, na.rm = TRUE),
    median_contra = median(df$contra_mean, na.rm = TRUE),
    mean_diff = mean(df$ipsi_mean - df$contra_mean, na.rm = TRUE),
    median_diff = median(df$ipsi_mean - df$contra_mean, na.rm = TRUE),
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2],
    t = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
})


# multiple testing correction (e.g. FDR)
results <- results %>%
  mutate(p_FDR = p.adjust(p, method = "fdr"))
results


