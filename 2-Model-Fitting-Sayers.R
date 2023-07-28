# Chris Sayers
# Last updated: July 27, 2023

# Script designed to create and compare linear mixed-effects models
# Requires objects from 1-Data-Wrangling-Sayers.R

library(tidyverse)
library(ggpubr)
library(MuMIn)
library(glmmTMB)
library(performance)
library(car)
library(DHARMa)
library(lawstat)
library(emmeans)
library(multcomp)
library(multcompView)

# creating a data frame for modeling purposes
HgSamples <- CollectiveData %>%
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Hg_Concentration") %>% 
  # only including full species names
  filter(!(str_detect(Species_Latin_Name, " sp."))) %>%
  # removing NA values in key variables so the models can run
  filter(!is.na(Hg_Concentration), !is.na(Trophic_Niche),
         !is.na(Primary_Habitat), !is.na(Season),
         !is.na(Migratory_Status), !is.na(Mining_Present_Yes_No),
         !is.na(Site_Name), !is.na(Banding_Station_Name),
         !is.na(Family), !is.na(Species_Latin_Name)) %>%
  # Creating new column with natural-log transformed THg concentrations
  mutate(lHg_Concentration = log(Hg_Concentration)) %>% 
  # Creating new columns for temporal tests
  mutate(MonthYear = make_date(Year, Month),
         Date = make_date(Year, Month, Day),
         Julian_Date = yday(Date),
         Year = as.numeric(Year))

# calculating how many levels are in each variable
unique(HgSamples$Tissue_Type)
unique(HgSamples$Trophic_Niche)
unique(HgSamples$Primary_Habitat)
unique(HgSamples$Mining_Present_Yes_No)
unique(HgSamples$Season)
unique(HgSamples$Migratory_Status)
unique(HgSamples$Site_Name)
unique(HgSamples$Banding_Station_Name)
unique(HgSamples$Family)
unique(HgSamples$Species_Latin_Name)
unique(HgSamples$Year)

#----------------------- DATA VISUALIZATION ------------------------------------
# much of the strategy below is from Zurr et al. (2010)  https://doi.org/10.1111/j.2041-210X.2009.00001.x

# OUTLIERS & NORMALITY OF RESPONSE VARIABLE ----------
ggdensity(HgSamples$Hg_Concentration, xlab = "Hg Concentration (µg/g)") # some high outliers 

ggqqplot(HgSamples$Hg_Concentration, ylab = "Hg Concentration (µg/g)") # tails stray from far from normal
shapiro.test(HgSamples$Hg_Concentration) # W = 0.29137, p-value < 2.2e-16, not normal
# log-transformation may be necessary for linear models

# log-transformed data diagnostics
ggdensity(HgSamples$lHg_Concentration, xlab = "ln[Hg Concentration (µg/g)]")
ggqqplot(HgSamples$lHg_Concentration) # much better than before
shapiro.test(HgSamples$lHg_Concentration) # W = 0.98232, p-value = 5.134e-16, not normal

# Cleavland dotplot of raw data by Family
ggplot(HgSamples, aes(x = Hg_Concentration, y = Family, color = Tissue_Type)) +
  geom_point() +
  labs(x = "Hg Concentration (µg/g ww)", y = "Family")

# Cleavland dotplot of natural-log transformed data, no apparent outliers anymore
ggplot(HgSamples, aes(x = log(Hg_Concentration), y = Family, color = Tissue_Type)) +
  geom_point() +
  labs(x = "ln(Hg) Concentration (µg/g ww)", y = "Family")


# FUNCTIONAL TRAIT MODEL -------------------------------------------------------

# global functional trait model
fullmodel <- glmmTMB(log(Hg_Concentration) ~ 
                     
                     # Hg toxicokinetics
                     Tissue_Type +
                    
                     # functional traits
                     Trophic_Niche +
                     Primary_Habitat +
                     Migratory_Status +
                    
                     # local conditions
                     Mining_Present_Yes_No +
                  
                     # crossed random effects
                     (1 | Site_Name/Banding_Station_Name) +
                     (1 | Family/Species_Common_Name/Band_Num) +
                     (1 | Year),
                   data = HgSamples, family = "gaussian", REML = F)

summary(fullmodel)
as.data.frame(confint(fullmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(fullmodel)
car::Anova(fullmodel, type = 3)

# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(fullmodel)) # VERY close to 0

#library(lmerTest)
#library(ggResidpanel)
#fullmodel1 <- lmer(log(Hg_Concentration) ~ Tissue_Type + Trophic_Niche + 
#                    Primary_Habitat + Migratory_Status + Mining_Present_Yes_No +
#                    (1 | Site_Name/Banding_Station_Name) + (1 | Family/Species_Common_Name/Band_Num) + (1 | Year),
#                  data = HgSamples, REML = F)
#
#resid_panel(fullmodel1, plots = "all", type = NA, bins = 30,
#            smoother = T, qqline = T, qqbands = T, scale = 1,
#            theme = "bw", axis.text.size = 10, title.text.size = 12,
#            title.opt = TRUE, nrow = NULL)
## residual plots look great

simulateResiduals(fullmodel, plot = T, refit = F, use.u = T)
shapiro.test(residuals(fullmodel)) # W = 0.97427, p-value < 2.2e-16, not normal
# residual plots look okay

# Checking for data points with high leverage
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
inf <- influence_mixed(fullmodel)
infIndexPlot(inf)

# Checking for normality of random effects
rand <- as.data.frame(ranef(fullmodel)) %>% 
  filter(grpvar %in% c("Site_Name"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.90097, p-value = 0.00204, we are not normal

rand <- as.data.frame(ranef(fullmodel)) %>%
  filter(grpvar %in% c("Banding_Station_Name:Site_Name"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.92365, p-value = 0.001487, we are not normal

rand <- as.data.frame(ranef(fullmodel)) %>% 
  filter(grpvar %in% c("Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.94334, p-value = 0.02, we are not normal

rand <- as.data.frame(ranef(fullmodel)) %>% 
  filter(grpvar %in% c("Species_Common_Name:Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.9918, p-value = 0.07674, we are normal

rand <- as.data.frame(ranef(fullmodel)) %>% 
  filter(grpvar %in% c("Band_Num:Species_Common_Name:Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.96609, p-value < 2.2e-16, we are not normal

rand <- as.data.frame(ranef(fullmodel)) %>% 
  filter(grpvar %in% c("Year"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.90177, p-value = 0.1414, we are normal

# Checking for autocorrelation/independence
acf(HgSamples$Hg_Concentration) # raw data is autocorrelated
acf(residuals(fullmodel)) # random effects variable corrects for this
runs.test(residuals(fullmodel)) # we do not have autocorrelated data


# CANDIDATE MODEL SET -----------------------------------------------

fullmodel <- glmmTMB(log(Hg_Concentration) ~ 
                       
                       # Hg toxicokinetics
                       Tissue_Type +
                       
                       # functional traits
                       Trophic_Niche +
                       Primary_Habitat +
                       Migratory_Status +
                       
                       # local conditions
                       Mining_Present_Yes_No +
                       
                       # crossed random effects
                       (1 | Site_Name/Banding_Station_Name) +
                       (1 | Family/Species_Common_Name/Band_Num) +
                       (1 | Year),
                     data = HgSamples, family = "gaussian", REML = F)

summary(fullmodel)
as.data.frame(confint(fullmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(fullmodel)
car::Anova(fullmodel, type = 3)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out1 <- MuMIn::dredge(fullmodel, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out1)
options(na.action = "na.omit")
write.csv(d.out1, "Outputs/functional-trait-model-selection.csv")

# figuring out which REs are significant in our top model structure
m0 <- glmmTMB(log(Hg_Concentration) ~ Tissue_Type + Trophic_Niche + Primary_Habitat + Mining_Present_Yes_No, data = HgSamples, family = "gaussian", REML = F)
m1 <- glmmTMB(log(Hg_Concentration) ~ (1 | Site_Name/Banding_Station_Name) + Tissue_Type + Trophic_Niche + Primary_Habitat + Mining_Present_Yes_No, data = HgSamples, family = "gaussian", REML = F)
m2 <- glmmTMB(log(Hg_Concentration) ~ (1 | Family/Species_Common_Name/Band_Num) + Tissue_Type + Trophic_Niche + Primary_Habitat + Mining_Present_Yes_No, data = HgSamples, family = "gaussian", REML = F)
m3 <- glmmTMB(log(Hg_Concentration) ~ (1 | Year) + Tissue_Type + Trophic_Niche + Primary_Habitat  + Mining_Present_Yes_No, data = HgSamples, family = "gaussian", REML = F)
anova(m0, m1, m2, m3) # all have p < 0.05

# 1st place model by a long-shot
topmodel <- glmmTMB(log(Hg_Concentration) ~ Tissue_Type + Trophic_Niche + 
                   Primary_Habitat + Mining_Present_Yes_No +
                   (1 | Site_Name/Banding_Station_Name) + (1 | Family/Species_Common_Name/Band_Num) + (1 | Year),
                 data = HgSamples, family = "gaussian", REML = F)

summary(topmodel)
as.data.frame(confint(topmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
# RE estimates are SDs
# back-transformed ASGM-present estimate = 3.762148
exp(1.32499)
performance::r2(topmodel)
car::Anova(topmodel, type = 3)

# computing post-hoc comparisons to determine significant differences among the modeled means
emmeans(topmodel, "Tissue_Type", type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topmodel, "Trophic_Niche", type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topmodel, "Primary_Habitat", type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topmodel, "Mining_Present_Yes_No", type = "response") %>% 
  cld(Letter = "abcdefg")

boxplot(log(Hg_Concentration) ~ Tissue_Type, HgSamples)
boxplot(log(Hg_Concentration) ~ Trophic_Niche, HgSamples)
boxplot(log(Hg_Concentration) ~ Primary_Habitat, HgSamples)
boxplot(log(Hg_Concentration) ~ Mining_Present_Yes_No, HgSamples)
boxplot(log(Hg_Concentration) ~ Year, HgSamples)


# TEMPORAL MODEL -----------------------------------------------

BloodHgSamples <- HgSamples %>% 
  filter(Tissue_Type == "Blood_Hg_ppm")

# global temporal model
temporalmodel <- glmmTMB(log(Hg_Concentration) ~ 
                        
                        # temporal variation   
                        Season +
                        
                        # crossed random effects
                        (1 | Site_Name/Banding_Station_Name) +
                        (1 | Family/Species_Common_Name/Band_Num) +
                        (1 | Year),
                      
                      data = BloodHgSamples, family = "gaussian", REML = F)

summary(temporalmodel)
as.data.frame(confint(temporalmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(temporalmodel)
car::Anova(temporalmodel, type = 3)


# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(temporalmodel)) # VERY close to 0

#library(lmerTest)
#library(ggResidpanel)
#temporalmodel1 <- lmer(log(Hg_Concentration) ~ Season + (1 | Site_Name/Banding_Station_Name)
#                         + (1 | Family/Species_Common_Name/Band_Num) + (1 | Year),
#                         data = BloodHgSamples, REML = F)
#
#resid_panel(temporalmodel1, plots = "all", type = NA, bins = 30,
#            smoother = T, qqline = T, qqbands = T, scale = 1,
#            theme = "bw", axis.text.size = 10, title.text.size = 12,
#            title.opt = TRUE, nrow = NULL)
## residual plots look great

simulateResiduals(temporalmodel, plot = T, refit = F, use.u = T)
shapiro.test(residuals(temporalmodel)) # W = 0.97751, p-value = 9.174e-11, not normal
# residual plots look okay

# Checking for data points with high leverage
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
inf <- influence_mixed(temporalmodel)
infIndexPlot(inf)

# Checking for normality of random effects
rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Site_Name"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.98091, p-value = 0.8716, we are normal

rand <- as.data.frame(ranef(temporalmodel)) %>%
  filter(grpvar %in% c("Banding_Station_Name:Site_Name"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.97825, p-value = 0.7014, we are normal

rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.9463, p-value = 0.07334, we are normal

rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Species_Common_Name:Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.98615, p-value = 0.1281, we are normal

rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Band_Num:Species_Common_Name:Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.97703, p-value = 8.22e-11, we are not normal

rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Year"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.94921, p-value = 0.6254, we are normal

# Checking for autocorrelation/independence
acf(BloodHgSamples$Hg_Concentration) # raw data is autocorrelated
acf(residuals(temporalmodel)) # random effects variable corrects for this
runs.test(residuals(temporalmodel)) # we do not have autocorrelated data


# CANDIDATE MODEL SET -----------------------------------------------

temporalmodel <- glmmTMB(log(Hg_Concentration) ~ 
                           
                           # temporal variation   
                           Season +
                           
                           # crossed random effects
                           (1 | Site_Name/Banding_Station_Name) +
                           (1 | Family/Species_Common_Name/Band_Num) +
                           (1 | Year),
                         
                         data = BloodHgSamples, family = "gaussian", REML = F)

summary(temporalmodel)
as.data.frame(confint(temporalmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(temporalmodel)
car::Anova(temporalmodel, type = 3)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out2 <- MuMIn::dredge(temporalmodel, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out2)
options(na.action = "na.omit")
write.csv(d.out2, "Outputs/temporal-model-selection.csv")

# figuring out which REs are significant in our top model structure
m0 <- glmmTMB(log(Hg_Concentration) ~ Season, data = HgSamples, family = "gaussian", REML = F)
m1 <- glmmTMB(log(Hg_Concentration) ~ (1 | Site_Name/Banding_Station_Name) + Season, data = HgSamples, family = "gaussian", REML = F)
m2 <- glmmTMB(log(Hg_Concentration) ~ (1 | Family/Species_Common_Name) + Season, data = HgSamples, family = "gaussian", REML = F)
m3 <- glmmTMB(log(Hg_Concentration) ~ (1 | Year) + Season, data = HgSamples, family = "gaussian", REML = F)
anova(m0, m1, m2, m3) # all have p < 0.05

# 1st place model by a long-shot
temporalmodel <- glmmTMB(log(Hg_Concentration) ~ Season + (1 | Site_Name/Banding_Station_Name)
                         + (1 | Family/Species_Common_Name/Band_Num) + (1 | Year),
                         data = BloodHgSamples, family = "gaussian", REML = F)

summary(temporalmodel)
as.data.frame(confint(temporalmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
# back-transformed wet season estimate = 0.763 = dry season is 31% higher Hg
exp(-0.2703)
performance::r2(temporalmodel)
car::Anova(temporalmodel, type = 3)


# computing post-hoc comparisons to determine significant differences among the modeled means
emmeans(temporalmodel, "Season", type = "response") %>% 
  cld(Letter = "abcdefg")
