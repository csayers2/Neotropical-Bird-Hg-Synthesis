# Chris Sayers
# Last updated: July 30, 2023

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
  #filter(!(str_detect(Species_Common_Name, " sp."))) %>% 
  # removing NA values in key variables so the models can run
  filter(!is.na(Hg_Concentration), !is.na(Trophic_Niche),
         !is.na(Primary_Habitat), !is.na(Season),
         !is.na(Migratory_Status), !is.na(Mining_Present_Yes_No),
         !is.na(Site_Name), !is.na(Banding_Station_Name),
         !is.na(Family), !is.na(Species_Latin_Name), !is.na(Band_Num)) %>%
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
length(unique(HgSamples$Band_Num))
unique(HgSamples$Year)

#----------------------- DATA VISUALIZATION ------------------------------------
# much of the strategy below is from Zurr et al. (2010)  https://doi.org/10.1111/j.2041-210X.2009.00001.x

# OUTLIERS & NORMALITY OF RESPONSE VARIABLE ----------
ggdensity(HgSamples$Hg_Concentration, xlab = "Hg Concentration (µg/g)") # some high outliers 

ggqqplot(HgSamples$Hg_Concentration, ylab = "Hg Concentration (µg/g)") # tails stray from far from normal
shapiro.test(HgSamples$Hg_Concentration) # W = 0.28904, p-value < 2.2e-16, not normal
# log-transformation will be necessary for linear models

# log-transformed data diagnostics
ggdensity(HgSamples$lHg_Concentration, xlab = "ln[Hg Concentration (µg/g)]")
ggqqplot(HgSamples$lHg_Concentration) # much better than before
shapiro.test(HgSamples$lHg_Concentration) # W = 0.98234, p-value = 2.547e-16, not normal

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
FTmodel <- glmmTMB(log(Hg_Concentration) ~ 
                     
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
                     (1 | Family/Species_Latin_Name/Band_Num) +
                     (1 | Year),
                   
                   data = HgSamples, family = "gaussian", REML = F)

summary(FTmodel)
as.data.frame(confint(FTmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(FTmodel)
car::Anova(FTmodel, type = 3)


# CHECKING MODEL ASSUMPTIONS -------------------------------------
# Checking for homogeneity of variance & normality of residuals
mean(residuals(FTmodel)) # VERY close to 0

#library(lmerTest)
#library(ggResidpanel)
#FTmodel1 <- lmer(log(Hg_Concentration) ~ Tissue_Type + Trophic_Niche + 
#                    Primary_Habitat + Migratory_Status + Mining_Present_Yes_No +
#                    (1 | Site_Name/Banding_Station_Name) + (1 | Family/Species_Common_Name/Band_Num) + (1 | Year),
#                  data = HgSamples, REML = F)
#
#resid_panel(FTmodel1, plots = "all", type = NA, bins = 30,
#            smoother = T, qqline = T, qqbands = T, scale = 1,
#            theme = "bw", axis.text.size = 10, title.text.size = 12,
#            title.opt = TRUE, nrow = NULL)
## residual plots look great

simulateResiduals(FTmodel, plot = T, refit = F, use.u = T)
shapiro.test(residuals(FTmodel)) # W = 0.97453, p-value < 2.2e-16, not normal
# residual plots look okay

# Checking for data points with high leverage
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
inf <- influence_mixed(FTmodel)
infIndexPlot(inf)

# Checking for normality of random effects
rand <- as.data.frame(ranef(FTmodel)) %>% 
  filter(grpvar %in% c("Site_Name"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.91255, p-value = 0.003976, not normal

rand <- as.data.frame(ranef(FTmodel)) %>%
  filter(grpvar %in% c("Banding_Station_Name:Site_Name"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.93048, p-value = 0.002543, not normal

rand <- as.data.frame(ranef(FTmodel)) %>% 
  filter(grpvar %in% c("Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.94058, p-value = 0.01296, not normal

rand <- as.data.frame(ranef(FTmodel)) %>% 
  filter(grpvar %in% c("Species_Latin_Name:Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.99215, p-value = 0.0872, normal

rand <- as.data.frame(ranef(FTmodel)) %>% 
  filter(grpvar %in% c("Band_Num:Species_Latin_Name:Family"))
ggqqplot(rand$condval) # tails stray from normal
shapiro.test(rand$condval) # W = 0.96647, p-value < 2.2e-16, not normal

rand <- as.data.frame(ranef(FTmodel)) %>% 
  filter(grpvar %in% c("Year"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.91285, p-value = 0.2006, normal

# Checking for autocorrelation/independence
acf(HgSamples$Hg_Concentration) # raw data is autocorrelated
acf(residuals(FTmodel)) # random effects variable corrects for this
runs.test(residuals(FTmodel)) # we do not have autocorrelated data


# CANDIDATE MODEL SET -----------------------------------------------

FTmodel <- glmmTMB(log(Hg_Concentration) ~ 
                       
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
                       (1 | Family/Species_Latin_Name/Band_Num) +
                       (1 | Year),
                   
                     data = HgSamples, family = "gaussian", REML = F)

summary(FTmodel)
as.data.frame(confint(FTmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(FTmodel)
car::Anova(FTmodel, type = 3)

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out1 <- MuMIn::dredge(FTmodel, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out1)
options(na.action = "na.omit")
write.csv(d.out1, "Outputs/functional-trait-model-selection.csv")

# figuring out which REs are significant in our top model structure
m0 <- glmmTMB(log(Hg_Concentration) ~ Tissue_Type + Trophic_Niche + Primary_Habitat + Mining_Present_Yes_No, data = HgSamples, family = "gaussian", REML = F)
m1 <- glmmTMB(log(Hg_Concentration) ~ (1 | Site_Name/Banding_Station_Name) + Tissue_Type + Trophic_Niche + Primary_Habitat + Mining_Present_Yes_No, data = HgSamples, family = "gaussian", REML = F)
m2 <- glmmTMB(log(Hg_Concentration) ~ (1 | Family/Species_Latin_Name/Band_Num) + Tissue_Type + Trophic_Niche + Primary_Habitat + Mining_Present_Yes_No, data = HgSamples, family = "gaussian", REML = F)
m3 <- glmmTMB(log(Hg_Concentration) ~ (1 | Year) + Tissue_Type + Trophic_Niche + Primary_Habitat  + Mining_Present_Yes_No, data = HgSamples, family = "gaussian", REML = F)
anova(m0, m1, m2, m3) # all have p < 0.05

# 1st place model by a long-shot
topFTmodel <- glmmTMB(log(Hg_Concentration) ~ Tissue_Type + Trophic_Niche + 
                        Primary_Habitat + Mining_Present_Yes_No + 
                        (1 | Site_Name/Banding_Station_Name) +
                        (1 | Family/Species_Latin_Name/Band_Num) + (1 | Year),
                      data = HgSamples, family = "gaussian", REML = F)

summary(topFTmodel)
as.data.frame(confint(topFTmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
# RE estimates are SDs
# back-transformed ASGM-present estimate = 3.78
performance::r2(topFTmodel)
car::Anova(topFTmodel, type = 3)

# computing post-hoc comparisons to determine significant differences among the modeled means
emmeans(topFTmodel, ~ Tissue_Type, type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topFTmodel, ~ Trophic_Niche, type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topFTmodel, ~ Primary_Habitat, type = "response") %>% 
  cld(Letter = "abcdefg")

emmeans(topFTmodel, ~ Mining_Present_Yes_No, type = "response") %>% 
  cld(Letter = "abcdefg")

boxplot(log(Hg_Concentration) ~ Tissue_Type, HgSamples)
boxplot(log(Hg_Concentration) ~ Trophic_Niche, HgSamples)
boxplot(log(Hg_Concentration) ~ Primary_Habitat, HgSamples)
boxplot(log(Hg_Concentration) ~ Mining_Present_Yes_No, HgSamples)
boxplot(log(Hg_Concentration) ~ Year, HgSamples)


# TEMPORAL MODEL -----------------------------------------------

BloodHgSamples <- HgSamples %>% 
  filter(Tissue_Type == "Blood_Hg_ppm", !is.na(Season),
         Migratory_Status == "Resident")

# global temporal model
temporalmodel <- glmmTMB(log(Hg_Concentration) ~ 
                        
                        # temporal variation   
                        Season*Trophic_Niche +
                        
                        # crossed random effects
                        (1 | Site_Name/Banding_Station_Name) +
                        (1 | Family/Species_Latin_Name/Band_Num) +
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
shapiro.test(residuals(temporalmodel)) # W = 0.96664, p-value = 2.821e-10, not normal
# residual plots look okay

# Checking for data points with high leverage
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
inf <- influence_mixed(temporalmodel)
infIndexPlot(inf)

# Checking for normality of random effects
rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Site_Name"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.98724, p-value = 0.978, normal

rand <- as.data.frame(ranef(temporalmodel)) %>%
  filter(grpvar %in% c("Banding_Station_Name:Site_Name"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.96075, p-value = 0.2557, normal

rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.88479, p-value = 0.001581, not normal

rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Species_Latin_Name:Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.98291, p-value = 0.1511, not normal

rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Band_Num:Species_Latin_Name:Family"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.97017, p-value = 1.82e-09, not normal

rand <- as.data.frame(ranef(temporalmodel)) %>% 
  filter(grpvar %in% c("Year"))
ggqqplot(rand$condval) # few tail stragglers, but the rest looks okay
shapiro.test(rand$condval) # W = 0.94322, p-value = 0.5409, normal

# Checking for autocorrelation/independence
acf(BloodHgSamples$Hg_Concentration) # raw data is autocorrelated
acf(residuals(temporalmodel)) # random effects variable corrects for this
runs.test(residuals(temporalmodel)) # we do not have autocorrelated data


# CANDIDATE MODEL SET -----------------------------------------------

# refactor reference level
#BloodHgSamples$Trophic_Niche <- as.factor(BloodHgSamples$Trophic_Niche)
#BloodHgSamples$Trophic_Niche <- relevel(BloodHgSamples$Trophic_Niche, ref = "Aquatic predator")

# global temporal model
temporalmodel <- glmmTMB(log(Hg_Concentration) ~ 
                           
                           # temporal variation   
                           Season*Trophic_Niche +
                           
                           # crossed random effects
                           (1 | Site_Name/Banding_Station_Name) +
                           (1 | Family/Species_Latin_Name/Band_Num) +
                           (1 | Year),
                         
                         data = BloodHgSamples, family = "gaussian", REML = F)

summary(temporalmodel)
as.data.frame(confint(temporalmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(temporalmodel)
car::Anova(temporalmodel, type = 3)

boxplot(log(Hg_Concentration) ~ Season, BloodHgSamples)
t.test(log(Hg_Concentration) ~ Season, BloodHgSamples)
exp(-3.793337) # dry season Hg estimate
exp(-3.342712) # wet season Hg estimate
0.03534098/0.02252033

options(na.action = "na.fail")
# computes marginal and conditional R^2
d.out2 <- MuMIn::dredge(temporalmodel, extra = list("Rsq" = function(x){performance::r2(x)}))
View(d.out2)
options(na.action = "na.omit")
write.csv(d.out2, "Outputs/temporal-model-selection.csv")

# figuring out which REs are significant in our top model structure
m0 <- glmmTMB(log(Hg_Concentration) ~ Season, data = HgSamples, family = "gaussian", REML = F)
m1 <- glmmTMB(log(Hg_Concentration) ~ (1 | Site_Name/Banding_Station_Name) + Season, data = HgSamples, family = "gaussian", REML = F)
m2 <- glmmTMB(log(Hg_Concentration) ~ (1 | Family/Species_Latin_Name/Band_Num) + Season, data = HgSamples, family = "gaussian", REML = F)
m3 <- glmmTMB(log(Hg_Concentration) ~ (1 | Year) + Season, data = HgSamples, family = "gaussian", REML = F)
anova(m0, m1, m2, m3) # all have p < 0.05

# 1st place model by a long-shot
temporalmodel <- glmmTMB(log(Hg_Concentration) ~ Season + Trophic_Niche +
                           (1 | Site_Name/Banding_Station_Name) +
                           (1 | Family/Species_Common_Name/Band_Num) + (1 | Year),
                         data = BloodHgSamples, family = "gaussian", REML = F)

summary(temporalmodel)
as.data.frame(confint(temporalmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
# back-transformed wet season estimate = 0.763 = dry season is 31% higher Hg
performance::r2(temporalmodel)
car::Anova(temporalmodel, type = 3)


# computing post-hoc comparisons to determine significant differences among the modeled means
emmeans(temporalmodel, ~ Season*Trophic_Niche, type = "response") %>% 
  cld(Letter = "abcdefg")



# temporal model to determine aggregate seasonal effect across trophic niches
temporalmodel <- glmmTMB(log(Hg_Concentration) ~ Season + 
                           (1 | Site_Name/Banding_Station_Name) +
                           (1 | Family/Species_Latin_Name/Band_Num) +
                           (1 | Year),
                         data = BloodHgSamples, family = "gaussian", REML = F)

summary(temporalmodel)
as.data.frame(confint(temporalmodel)) %>% 
  mutate(Estimate = exp(Estimate), `2.5 %` = exp(`2.5 %`), `97.5 %` = exp(`97.5 %`))
performance::r2(temporalmodel)
car::Anova(temporalmodel, type = 3)



BloodHgSamples <- HgSamples %>% 
  filter(Tissue_Type == "Blood_Hg_ppm")

# calculating tissue sample sizes for y axis
ss <- BloodHgSamples %>%
  group_by(Trophic_Niche) %>%
  summarize(n = n()) %>%
  mutate(Trophic_Niche.s = str_c(Trophic_Niche, "\n(n = ", n, ")"))

# calculating predicted Hg values with top model structure
# type = "random" gives prediction intervals rather than confidence intervals
pr <- ggpredict(temporalmodel, terms = c("Trophic_Niche", "Season"),
                type = "random", back.transform = T) %>%
  rename(Trophic_Niche = x, Season = group) %>%
  left_join(ss, by = "Trophic_Niche") %>% 
  # this function is critical to order the facets properly
  transform(Season = factor(Season, levels = c("Wet", "Dry"),
                            labels = c("Wet", "Dry"))) %>% 
  # ordering the levels based on maximum predicted mean across tissue types
  group_by(Trophic_Niche) %>% 
  mutate(max_predicted = max(predicted)) %>% 
  #group_by(Season) %>% 
  #mutate(spacing = min(conf.low)/2)
  group_by(Season) %>% 
  mutate(spacing = max(conf.high) * 2)

# create final data frame with raw data and predicted means to plot
df <- BloodHgSamples %>% 
  # this function is critical to order the facets properly
  transform(Season = factor(Season, levels = c("Wet", "Dry"),
                            labels = c("Wet", "Dry"))) %>%  
  full_join(pr, by = c("Trophic_Niche", "Season"))

ggplot() +
  geom_point(data = df, mapping = aes(x = Hg_Concentration, y = reorder(Trophic_Niche.s, max_predicted),
                                      fill = Season),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
             color = "#E41A1C", alpha = 0.4, size = 2) +
  geom_pointrange(data = df, aes(x = predicted, y = reorder(Trophic_Niche.s, max_predicted),
                                 xmin = conf.low, xmax = conf.high, shape = Season),
                  position = position_dodge(0.6), size = 0.6, linewidth = 0.6) +
  labs(x = "Predicted THg (µg/g)", y = "Trophic niche") +
  scale_x_continuous(expand = c(0.1, 0),
                     trans = "log",
                     breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 75),
                     labels = c("0", "0.001", "0.01", "0.1", "1", "10", "75")) + 
  scale_shape_manual(limits = c("Dry", "Wet"), values = c(17, 16)) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        aspect.ratio = 1)
