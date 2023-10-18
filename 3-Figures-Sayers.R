# Chris Sayers
# Last updated: July 30, 2023

# Script designed to produce graphs featured in publication and annex
# Requires objects from 1-Data-Wrangling-Sayers.R and 2-Model-Fitting-Sayers.R

library(tidyverse)
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(scales)
library(glmmTMB)
library(ggeffects)
library(multcomp)
library(ggtext)
library(glue)
library(ggview)

# resolve namespace conflicts and creating necessary functions
select <- dplyr::select
"%nin%" <- Negate("%in%")

# 1st place models from model selection
topFTmodel <- glmmTMB(log(Hg_Concentration) ~ Tissue_Type + Trophic_Niche + 
                        Primary_Habitat + Mining_Present_Yes_No + 
                        (1 | Site_Name/Banding_Station_Name) +
                        (1 | Family/Species_Latin_Name/Band_Num) + (1 | Year),
                      data = HgSamples, family = "gaussian", REML = F)

temporalmodel <- glmmTMB(log(Hg_Concentration) ~ Season*Trophic_Niche +
                           (1 | Site_Name/Banding_Station_Name) +
                           (1 | Family/Species_Common_Name/Band_Num) + (1 | Year),
                         data = BloodHgSamples, family = "gaussian", REML = F)
# PREDICTED TROPHIC NICHE x ASGM ----------------------------------------------

# calculating tissue sample sizes for y axis
ss <- HgSamples %>%
  group_by(Trophic_Niche, Tissue_Type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Tissue_Type, values_from = n, values_fill = 0) %>%
  select(Trophic_Niche, n_Blood = Blood_Hg_ppm, n_Body = Body_Hg_ppm, n_Tail = Tail_Hg_ppm) %>% 
  mutate(Trophic_Niche.s = str_c(Trophic_Niche, "\n(n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"),
         n_Total = sum(n_Blood, n_Body, n_Tail))

# calculating predicted Hg values with top model structure
# type = "random" gives prediction intervals rather than confidence intervals
pr <- ggpredict(topFTmodel, terms = c("Trophic_Niche", "Tissue_Type", "Mining_Present_Yes_No"),
                type = "random", back.transform = T) %>%
  rename(Trophic_Niche = x, Tissue_Type = group, Mining_Present_Yes_No = facet) %>%
  left_join(ss, by = "Trophic_Niche") %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("No", "Yes"),
                                           labels = c("ASGM absent", "ASGM present"))) %>% 
  # adding post-hoc comparisons
  mutate(Tukey = ifelse(Trophic_Niche %in% c("Terrestrial vertivore", "Invertivore"), "a", 
                        ifelse(Trophic_Niche %in% c("Aquatic predator"), "ab",
                               ifelse(Trophic_Niche %in% c("Nectarivore"), "abc",
                                      ifelse(Trophic_Niche %in% c("Omnivore"), "b", "c"))))) %>%
  # only pasting the tukey levels once in the plot
  #mutate(Tukey = if_else(Tissue_Type == "Whole blood", Tukey, NA_character_)) %>%
  mutate(Tukey = if_else(Tissue_Type == "Tail feather", Tukey, NA_character_)) %>%
  # ordering the levels based on maximum predicted mean across tissue types
  group_by(Trophic_Niche) %>% 
  mutate(max_predicted = max(predicted)) %>% 
  #group_by(Tissue_Type) %>% 
  #mutate(spacing = min(conf.low)/2)
  group_by(Tissue_Type) %>% 
  mutate(spacing = max(conf.high) * 2) %>% 
  # adjusting lower bound of confidence interval
  mutate(conf.low = if_else(conf.low < 0.001, 0.001, conf.low))

# create final data frame with raw data and predicted means to plot
df <- HgSamples %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("No", "Yes"),
                                           labels = c("ASGM absent", "ASGM present"))) %>% 
  full_join(pr, by = c("Trophic_Niche", "Tissue_Type", "Mining_Present_Yes_No"))


ggplot() +
  geom_point(data = df, mapping = aes(x = Hg_Concentration, y = reorder(Trophic_Niche.s, max_predicted),
                                      color = Tissue_Type, fill = Mining_Present_Yes_No),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
             alpha = 0.4, size = 2) +
  geom_pointrange(data = df, aes(x = predicted, y = reorder(Trophic_Niche.s, max_predicted),
                                 xmin = conf.low, xmax = conf.high, shape = Mining_Present_Yes_No),
                  position = position_dodge(0.6), size = 0.6, linewidth = 0.6) +
  geom_text(data = df, mapping = aes(x = spacing, y = Trophic_Niche.s, label = Tukey), size = 5) +
  labs(x = "Predicted THg (µg/g)", y = "Trophic niche") +
  scale_x_continuous(expand = c(0.1, 0),
                     trans = "log",
                     breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 75),
                     labels = c("0", "0.001", "0.01", "0.1", "1", "10", "75")) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), guide = "none") +
  scale_shape_manual(limits = c("ASGM present", "ASGM absent"), values = c(17, 16)) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
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

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 5)

ggsave("Publication-Figures/Fig1_Predicted_TrophicxASGM_AllTissues.jpg", dpi = 1200, width = 15, height = 5)
#ggsave("Publication-Figures/Fig1_Predicted_TrophicxASGM_AllTissues.tiff", dpi = 600, width = 15, height = 5)


# PREDICTED PRIMARY HABITAT x ASGM ------------------------------------------

# calculating tissue sample sizes for y axis
ss <- HgSamples %>%
  group_by(Primary_Habitat, Tissue_Type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Tissue_Type, values_from = n, values_fill = 0) %>%
  select(Primary_Habitat, n_Blood = Blood_Hg_ppm, n_Body = Body_Hg_ppm, n_Tail = Tail_Hg_ppm) %>% 
  mutate(Primary_Habitat.s = str_c(Primary_Habitat, "\n(n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"),
         n_Total = sum(n_Blood, n_Body, n_Tail))

# calculating predicted Hg values with top model structure
# type = "random" gives prediction intervals rather than confidence intervals
pr <- ggpredict(topFTmodel, terms = c("Primary_Habitat", "Tissue_Type", "Mining_Present_Yes_No"),
                type = "random", back.transform = T) %>%
  rename(Primary_Habitat = x, Tissue_Type = group, Mining_Present_Yes_No = facet) %>%
  left_join(ss, by = "Primary_Habitat") %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("No", "Yes"),
                                           labels = c("ASGM absent", "ASGM present"))) %>% 
  # adding post-hoc comparisons
  mutate(Tukey = ifelse(Primary_Habitat %in% c("Aquatic"), "a", 
                        ifelse(Primary_Habitat %in% c("Lowland deciduous forest", "Grassland/scrub",
                                                      "Lowland evergreen forest"), "ab", "b"))) %>% 
  # only pasting the tukey levels once in the plot
  #mutate(Tukey = if_else(Tissue_Type == "Whole blood", Tukey, NA_character_)) %>%
  mutate(Tukey = if_else(Tissue_Type == "Tail feather", Tukey, NA_character_)) %>%
  # ordering the levels based on maximum predicted mean across tissue types
  group_by(Primary_Habitat) %>% 
  mutate(max_predicted = max(predicted)) %>% 
  #group_by(Tissue_Type) %>% 
  #mutate(spacing = min(conf.low)/2)
  group_by(Tissue_Type) %>% 
  mutate(spacing = max(conf.high) * 3) %>% 
  # adjusting lower bound of confidence interval
  mutate(conf.low = if_else(conf.low < 0.001, 0.001, conf.low))

# create final data frame with raw data and predicted means to plot
df <- HgSamples %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("No", "Yes"),
                                           labels = c("ASGM absent", "ASGM present"))) %>% 
  full_join(pr, by = c("Primary_Habitat", "Tissue_Type", "Mining_Present_Yes_No"))

ggplot() +
  geom_point(data = df, mapping = aes(x = Hg_Concentration, y = reorder(Primary_Habitat.s, max_predicted),
                                      color = Tissue_Type, fill = Mining_Present_Yes_No),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
             alpha = 0.4, size = 2) +
  geom_pointrange(data = df, aes(x = predicted, y = reorder(Primary_Habitat.s, max_predicted),
                      xmin = conf.low, xmax = conf.high, shape = Mining_Present_Yes_No),
                  position = position_dodge(0.6), size = 0.6, linewidth = 0.6) +
  geom_text(data = df, mapping = aes(x = spacing, y = Primary_Habitat.s, label = Tukey), size = 5) +
  labs(x = "Predicted THg (µg/g)", y = "Primary habitat") +
  scale_x_continuous(expand = c(0.1, 0),
                     trans = "log",
                     breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 75),
                     labels = c("0", "0.001", "0.01", "0.1", "1", "10", "75")) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), guide = "none") +
  scale_shape_manual(limits = c("ASGM present", "ASGM absent"), values = c(17, 16)) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        aspect.ratio = 1)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 5)

ggsave("Publication-Figures/Fig2_Predicted_HabitatxASGM_AllTissues.jpg", dpi = 1200, width = 15, height = 5)
#ggsave("Publication-Figures/Fig2_Predicted_HabitatxASGM_AllTissues.tiff", dpi = 600, width = 15, height = 5)


# PREDICTED FAMILY x ASGM -----------------------------------------------------

# calculating tissue sample sizes for y axis
ss <- HgSamples %>%
  group_by(Family, Tissue_Type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Tissue_Type, values_from = n, values_fill = 0) %>%
  select(Family, n_Blood = Blood_Hg_ppm, n_Body = Body_Hg_ppm, n_Tail = Tail_Hg_ppm) %>% 
  mutate(Family.s = str_c(Family, "\n(n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"),
         n_Total = sum(n_Blood, n_Body, n_Tail))

# calculating predicted Hg values with top model structure
# type = "random" gives prediction intervals rather than confidence intervals
pr <- ggpredict(topFTmodel, terms = c("Family", "Tissue_Type", "Mining_Present_Yes_No"),
                type = "random", back.transform = T) %>%
  rename(Family = x, Tissue_Type = group, Mining_Present_Yes_No = facet) %>%
  left_join(ss, by = "Family") %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("No", "Yes"),
                                           labels = c("ASGM absent", "ASGM present"))) %>% 
  # ordering the levels based on maximum predicted mean across tissue types
  group_by(Family) %>% 
  mutate(max_predicted = max(predicted)) %>% 
  #group_by(Tissue_Type) %>% 
  #mutate(spacing = min(conf.low)/2)
  group_by(Tissue_Type) %>% 
  mutate(spacing = max(conf.high) * 2) %>% 
  # adjusting lower bound of confidence interval
  mutate(conf.low = if_else(conf.low < 0.001, 0.001, conf.low))

# create final data frame with raw data and predicted means to plot
df <- HgSamples %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("No", "Yes"),
                                           labels = c("ASGM absent", "ASGM present"))) %>% 
  full_join(pr, by = c("Family", "Tissue_Type", "Mining_Present_Yes_No")) %>% 
  filter(n_Total > 24)

ggplot() +
  geom_point(data = df, mapping = aes(x = Hg_Concentration, y = reorder(Family.s, max_predicted),
                                      color = Tissue_Type, fill = Mining_Present_Yes_No),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
             alpha = 0.4, size = 2) +
  geom_pointrange(data = df, aes(x = predicted, y = reorder(Family.s, max_predicted),
                                 xmin = conf.low, xmax = conf.high, shape = Mining_Present_Yes_No),
                  position = position_dodge(0.6), size = 0.6, linewidth = 0.6) +
  labs(x = "Predicted THg (µg/g)", y = "Family (n ≥ 25)") +
  scale_x_continuous(expand = c(0.1, 0),
                     trans = "log",
                     breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 75),
                     labels = c("0", "0.001", "0.01", "0.1", "1", "10", "75")) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), guide = "none") +
  scale_shape_manual(limits = c("ASGM present", "ASGM absent"), values = c(17, 16)) +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  #theme_classic() +
  facet_grid(~ Tissue_Type, scales = "free") +
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
        aspect.ratio = 3)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 10)

ggsave("Publication-Figures/Fig3_Predicted_FamilyxASGM_AllTissues.jpg", dpi = 1200, width = 15, height = 10)
#ggsave("Publication-Figures/Fig3_Predicted_FamilyxASGM_AllTissues.tiff", dpi = 600, width = 15, height = 10)


# PREDICTED SITE -----------------------------------------------------

# calculating tissue sample sizes for y axis
ss <- HgSamples %>%
  mutate(Site_Name2 = str_c(Site_Name, ", ", Country)) %>%
  # denoting ASGM presence
  mutate(Site_Name2 = if_else(Mining_Present_Yes_No == "Yes", str_c("* ", Site_Name2), Site_Name2)) %>% 
  group_by(Site_Name, Site_Name2, Tissue_Type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Tissue_Type, values_from = n, values_fill = 0) %>%
  select(Site_Name, Site_Name2, n_Blood = Blood_Hg_ppm, n_Body = Body_Hg_ppm, n_Tail = Tail_Hg_ppm) %>%
  mutate(Site_Name.s = str_c(Site_Name2, "\n(n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"), 
         n_Total = sum(n_Blood, n_Body, n_Tail))

# calculating predicted Hg values with top model structure
# type = "random" gives prediction intervals rather than confidence intervals
pr <- ggpredict(topFTmodel, terms = c("Site_Name", "Tissue_Type"),
                type = "random", back.transform = T) %>%
  rename(Site_Name = x, Tissue_Type = group) %>%
  left_join(ss, by = "Site_Name") %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  # ordering the levels based on maximum predicted mean across tissue types
  group_by(Site_Name) %>% 
  mutate(max_predicted = max(predicted)) %>% 
  #group_by(Tissue_Type) %>% 
  #mutate(spacing = min(conf.low)/2)
  group_by(Tissue_Type) %>% 
  mutate(spacing = max(conf.high) * 2) %>% 
  # adjusting lower bound of confidence interval
  mutate(conf.low = if_else(conf.low < 0.001, 0.001, conf.low))

# create final data frame with raw data and predicted means to plot
df <- HgSamples %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  full_join(pr, by = c("Site_Name", "Tissue_Type")) %>% 
  filter(n_Total > 24)

ggplot() +
  geom_point(data = df, mapping = aes(x = Hg_Concentration, y = reorder(Site_Name.s, max_predicted),
                                      color = Tissue_Type),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
             alpha = 0.4, size = 2) +
  geom_pointrange(data = df, aes(x = predicted, y = reorder(Site_Name.s, max_predicted),
                                 xmin = conf.low, xmax = conf.high),
                  position = position_dodge(0.6), size = 0.6, linewidth = 0.6) +
  labs(x = "Predicted THg (µg/g)", y = "Site (n ≥ 25)") +
  scale_x_continuous(expand = c(0.1, 0),
                     trans = "log",
                     breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 75),
                     labels = c("0", "0.001", "0.01", "0.1", "1", "10", "75")) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), guide = "none") +
  scale_fill_discrete(guide = "none") +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
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
        aspect.ratio = 3)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 10)

ggsave("Publication-Figures/Predicted_Site_AllTissues_Hist.jpg", dpi = 1200, width = 15, height = 10)
#ggsave("Publication-Figures/Predicted_Site_AllTissues_Hist.tiff", dpi = 600, width = 15, height = 10)


# PREDICTED TROPHIC NICHE x SEASON ---------------------------------------------

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
  mutate(spacing = max(conf.high) * 2) %>% 
  # adjusting lower bound of confidence interval
  mutate(conf.low = if_else(conf.low < 0.001, 0.001, conf.low))

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

ggview(device = "jpeg", units = "in", dpi = 1200, width = 10, height = 6)

ggsave("Publication-Figures/Predicted_TrophicxSeason_Blood_Hist.jpg", dpi = 1200, width = 7, height = 5)

# TROPHIC NICHE VARIATION ---------------------------------------------

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Trophic_Niche), !is.na(Concentration)) %>%
  group_by(Trophic_Niche, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration)) %>%
  pivot_wider(names_from = Tissue_Type, values_from = c(n, mean, sd), values_fill = list(n = 0, mean = NA, sd = NA)) %>%
  select(Trophic_Niche, n_Blood = n_Blood_Hg_ppm, n_Body = n_Body_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         mean_Blood = mean_Blood_Hg_ppm, mean_Body = mean_Body_Hg_ppm, mean_Tail = mean_Tail_Hg_ppm,
         sd_Blood = sd_Blood_Hg_ppm, sd_Body = sd_Body_Hg_ppm, sd_Tail = sd_Tail_Hg_ppm) %>% 
  mutate(Trophic_Niche = str_c(Trophic_Niche, "\n(n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"))

p1 <- df %>%
  select(Trophic_Niche, mean = mean_Blood, sd = sd_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Trophic_Niche, mean = mean_Body, sd = sd_Body, n_Body = n_Body) %>%
  mutate(Tissue_Type = "Body feather")
p3 <- df %>%
  select(Trophic_Niche, mean = mean_Tail, sd = sd_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3) %>%
  filter(!is.na(mean)) %>%
  group_by(Trophic_Niche) %>% 
  mutate(sum = sum(n_Blood, n_Body, n_Tail, na.rm = T)) %>% 
  # creating a system to better rank the y axis
  #group_by(Trophic_Niche) %>% 
  #mutate(global_mean = mean(mean)) %>%
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Trophic_Niche) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather")))

ggplot(final, mapping = aes(x = mean, y = reorder(Trophic_Niche, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15, linewidth = 0.6) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Trophic niche") +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        aspect.ratio = 1)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 5)

ggsave("Publication-Figures/TrophicNiche_AllTissues_Hist.jpg", dpi = 1200, width = 15, height = 5)


# SPECIES HABITAT VARIATION ---------------------------------------------

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Primary_Habitat), !is.na(Concentration)) %>%
  group_by(Primary_Habitat, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration)) %>%
  pivot_wider(names_from = Tissue_Type, values_from = c(n, mean, sd), values_fill = list(n = 0, mean = NA, sd = NA)) %>%
  select(Primary_Habitat, n_Blood = n_Blood_Hg_ppm, n_Body = n_Body_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         mean_Blood = mean_Blood_Hg_ppm, mean_Body = mean_Body_Hg_ppm, mean_Tail = mean_Tail_Hg_ppm,
         sd_Blood = sd_Blood_Hg_ppm, sd_Body = sd_Body_Hg_ppm, sd_Tail = sd_Tail_Hg_ppm) %>% 
  mutate(Primary_Habitat = str_c(Primary_Habitat, "\n(n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"))

p1 <- df %>%
  select(Primary_Habitat, mean = mean_Blood, sd = sd_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Primary_Habitat, mean = mean_Body, sd = sd_Body, n_Body = n_Body) %>%
  mutate(Tissue_Type = "Body feather")
p3 <- df %>%
  select(Primary_Habitat, mean = mean_Tail, sd = sd_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3) %>%
  filter(!is.na(mean)) %>%
  group_by(Primary_Habitat) %>% 
  mutate(sum = sum(n_Blood, n_Body, n_Tail, na.rm = T)) %>%
  # creating a system to better rank the y axis
  #group_by(Primary_Habitat) %>% 
  #mutate(global_mean = mean(mean)) %>%
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Primary_Habitat) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather")))

ggplot(final, mapping = aes(x = mean, y = reorder(Primary_Habitat, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15, linewidth = 0.6) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Primary habitat") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        aspect.ratio = 1)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 5)

ggsave("Publication-Figures/SpeciesHabitat_AllTissues_Hist.jpg", dpi = 1200, width = 15, height = 5)


# ORDER VARIATION ---------------------------------------------------------

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Order) | Order == "", !is.na(Concentration)) %>%
  group_by(Order, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration)) %>%
  pivot_wider(names_from = Tissue_Type, values_from = c(n, mean, sd), values_fill = list(n = 0, mean = NA, sd = NA)) %>%
  select(Order, n_Blood = n_Blood_Hg_ppm, n_Body = n_Body_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         mean_Blood = mean_Blood_Hg_ppm, mean_Body = mean_Body_Hg_ppm, mean_Tail = mean_Tail_Hg_ppm,
         sd_Blood = sd_Blood_Hg_ppm, sd_Body = sd_Body_Hg_ppm, sd_Tail = sd_Tail_Hg_ppm) %>% 
  mutate(Order = str_c(Order, " (n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"))

p1 <- df %>%
  select(Order, mean = mean_Blood, sd = sd_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Order, mean = mean_Body, sd = sd_Body, n_Body = n_Body) %>%
  mutate(Tissue_Type = "Body feather")
p3 <- df %>%
  select(Order, mean = mean_Tail, sd = sd_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3) %>%
  filter(!is.na(mean)) %>%
  group_by(Order) %>% 
  mutate(sum = sum(n_Blood, n_Body, n_Tail, na.rm = T)) %>% 
  # creating a system to better rank the y axis
  #group_by(Order) %>% 
  #mutate(global_mean = mean(mean)) %>%
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Order) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather")))

ggplot(final, mapping = aes(x = mean, y = reorder(Order, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15, linewidth = 0.6) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Order") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"), 
        aspect.ratio = 1)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 5)

ggsave("Publication-Figures/Order_AllTissues_Hist.jpg", dpi = 1200, width = 15, height = 5)


# FAMILY VARIATION --------------------------------------------------------

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Family), !is.na(Concentration)) %>%
  group_by(Family, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration)) %>%
  pivot_wider(names_from = Tissue_Type, values_from = c(n, mean, sd), values_fill = list(n = 0, mean = NA, sd = NA)) %>%
  select(Family, n_Blood = n_Blood_Hg_ppm, n_Body = n_Body_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         mean_Blood = mean_Blood_Hg_ppm, mean_Body = mean_Body_Hg_ppm, mean_Tail = mean_Tail_Hg_ppm,
         sd_Blood = sd_Blood_Hg_ppm, sd_Body = sd_Body_Hg_ppm, sd_Tail = sd_Tail_Hg_ppm) %>% 
  mutate(Family = str_c(Family, " (n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"))

p1 <- df %>%
  select(Family, mean = mean_Blood, sd = sd_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Family, mean = mean_Body, sd = sd_Body, n_Body = n_Body) %>%
  mutate(Tissue_Type = "Body feather")
p3 <- df %>%
  select(Family, mean = mean_Tail, sd = sd_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3) %>%
  filter(!is.na(mean)) %>%
  group_by(Family) %>% 
  mutate(sum = sum(n_Blood, n_Body, n_Tail, na.rm = T)) %>% 
  # creating a system to better rank the y axis
  #group_by(Family) %>% 
  #mutate(global_mean = mean(mean)) %>%
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Family) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather")))

ggplot(final, mapping = aes(x = mean, y = reorder(Family, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.2, linewidth = 0.6) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Family") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        aspect.ratio = 3)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 10)

ggsave("Publication-Figures/Family_AllTissues_Hist.jpg", dpi = 1200, width = 15, height = 10)


# SPECIES VARIATION --------------------------------------------------------

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Species_Latin_Name), !is.na(Concentration)) %>%
  group_by(Species_Latin_Name, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration)) %>%
  pivot_wider(names_from = Tissue_Type, values_from = c(n, mean, sd), values_fill = list(n = 0, mean = NA, sd = NA)) %>%
  select(Species_Latin_Name, n_Blood = n_Blood_Hg_ppm, n_Body = n_Body_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         mean_Blood = mean_Blood_Hg_ppm, mean_Body = mean_Body_Hg_ppm, mean_Tail = mean_Tail_Hg_ppm,
         sd_Blood = sd_Blood_Hg_ppm, sd_Body = sd_Body_Hg_ppm, sd_Tail = sd_Tail_Hg_ppm) %>% 
  mutate(Sample_Size = str_c("n = ", n_Blood, ", ", n_Body, ", ", n_Tail)) %>%
  mutate(Species_Latin_Name = glue("*{Species_Latin_Name}* ({Sample_Size})"))

p1 <- df %>%
  select(Species_Latin_Name, mean = mean_Blood, sd = sd_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Species_Latin_Name, mean = mean_Body, sd = sd_Body, n_Body = n_Body) %>%
  mutate(Tissue_Type = "Body feather")
p3 <- df %>%
  select(Species_Latin_Name, mean = mean_Tail, sd = sd_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3) %>%
  filter(!is.na(mean)) %>%
  group_by(Species_Latin_Name) %>% 
  mutate(sum = sum(n_Blood, n_Body, n_Tail, na.rm = T)) %>% 
  filter(sum > 9) %>%
  # creating a system to better rank the y axis
  #group_by(Species_Latin_Name) %>% 
  #mutate(global_mean = mean(mean)) %>%
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Species_Latin_Name) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather")))

ggplot(final, mapping = aes(x = mean, y = reorder(Species_Latin_Name, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15, linewidth = 0.6) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Species (n ≥ 10)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_markdown(), # this is the key step for italic species names
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        aspect.ratio = 3)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 10)

ggsave("Publication-Figures/Species_AllTissues_Hist.jpg", dpi = 1200, width = 15, height = 10)


# SITE VARIATION ---------------------------------------------------------

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  filter(Site_Name %nin% c("Unknown")) %>% 
  mutate(Site_Name = str_c(Site_Name, ", ", Country)) %>%
  mutate(Site_Name = if_else(Mining_Present_Yes_No == "Yes", str_c("* ", Site_Name), Site_Name)) %>% 
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Site_Name), !is.na(Concentration)) %>%
  group_by(Site_Name, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration)) %>%
  pivot_wider(names_from = Tissue_Type, values_from = c(n, mean, sd), values_fill = list(n = 0, mean = NA, sd = NA)) %>%
  select(Site_Name, n_Blood = n_Blood_Hg_ppm, n_Body = n_Body_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         mean_Blood = mean_Blood_Hg_ppm, mean_Body = mean_Body_Hg_ppm, mean_Tail = mean_Tail_Hg_ppm,
         sd_Blood = sd_Blood_Hg_ppm, sd_Body = sd_Body_Hg_ppm, sd_Tail = sd_Tail_Hg_ppm) %>% 
  mutate(Site_Name = str_c(Site_Name, " (n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"),
         n_Total = sum(n_Blood, n_Body, n_Tail, na.rm = T)) %>% 
  filter(n_Total > 0)

p1 <- df %>%
  select(Site_Name, mean = mean_Blood, sd = sd_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Site_Name, mean = mean_Body, sd = sd_Body, n_Body = n_Body) %>%
  mutate(Tissue_Type = "Body feather")
p3 <- df %>%
  select(Site_Name, mean = mean_Tail, sd = sd_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3) %>%
  filter(!is.na(mean)) %>%
  # creating a system to better rank the y axis
  #group_by(Site_Name) %>% 
  #mutate(global_mean = mean(mean)) %>%
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Site_Name) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather"))) 

ggplot(final, mapping = aes(x = mean, y = reorder(Site_Name, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15, linewidth = 0.6) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Site") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_classic(base_size = 14) +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        aspect.ratio = 3)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 15, height = 10)

ggsave("Publication-Figures/Site_AllTissues_Hist.jpg", dpi = 1200, width = 15, height = 10)


# RISK ASSESSMENT ---------------------------------------------------------

# Risk graph for blood
bloodrisk <- CollectiveData %>%
  filter(!is.na(Blood_Hg_ppm), !is.na(Species_Latin_Name)) %>%
  group_by(Species_Latin_Name) %>%
  mutate(Sample_Size = str_c("n = ", n()),
         Species_Latin_Name = glue("*{Species_Latin_Name}* ({Sample_Size})"),
         Risk_Levels = cut(Blood_Hg_ppm, breaks = c(0, 0.7, 1.2, 1.7, 2.2, Inf),
                           labels = c("< 0.7 µg/g ww", "≥ 0.7 µg/g ww", "≥ 1.2 µg/g ww", "≥ 1.7 µg/g ww", "≥ 2.2 µg/g ww")),
         # reordering risk categories
         Risk_Levels = factor(Risk_Levels, levels = c("≥ 2.2 µg/g ww","≥ 1.7 µg/g ww", "≥ 1.2 µg/g ww", "≥ 0.7 µg/g ww", "< 0.7 µg/g ww"))) %>%
  filter(n() > 4) %>%
  # Trying to create some sort of weighted average for the ordering of the y axis - this was sooo painstaking
  mutate(Ext = ifelse(Risk_Levels == "≥ 2.2 µg/g ww", 1, 0),
         Prop_Ext = sum(Ext)/n(),
         High = ifelse(Risk_Levels == "≥ 1.7 µg/g ww", 1, 0),
         Prop_High = sum(High)/n(),
         Med = ifelse(Risk_Levels == "≥ 1.2 µg/g ww", 1, 0),
         Prop_Med = sum(Med)/n(),
         Low = ifelse(Risk_Levels == "≥ 0.7 µg/g ww", 1, 0),
         Prop_Low = sum(Low)/n(),
         None = ifelse(Risk_Levels == "< 0.7 µg/g ww", 1, 0),
         Prop_None = sum(None)/n(),
         Score = sum(Ext, High, Med, Low, None)) %>%
  ggplot(mapping = aes(y = reorder(reorder(reorder(reorder(
    reorder(reorder(reorder(Species_Latin_Name, desc(Species_Latin_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightgray")) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)",
       fill = "Whole blood\nrisk categories") +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_markdown(), # this is the key step for italic species names
        legend.title = element_text(face = "bold"),
        #legend.position = "none",
        aspect.ratio = 1)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 12, height = 8)


# Risk graph for body feathers
bodyrisk <- CollectiveData %>%
  filter(!is.na(Species_Latin_Name), !is.na(Body_Hg_ppm)) %>%
  group_by(Species_Latin_Name) %>%
  mutate(Sample_Size = str_c("n = ", n()),
         Species_Latin_Name = glue("*{Species_Latin_Name}* ({Sample_Size})"),
         Risk_Levels = cut(Body_Hg_ppm, breaks = c(0, 2.4, 3.4, 4.5, 5.3, Inf),
                           labels = c("< 2.4 µg/g fw", "≥ 2.4 µg/g fw", "≥ 3.4 µg/g fw", "≥ 4.5 µg/g fw", "≥ 5.3 µg/g fw")),
         # reordering risk categories
         Risk_Levels = factor(Risk_Levels, levels = c("≥ 5.3 µg/g fw", "≥ 4.5 µg/g fw",
                                                      "≥ 3.4 µg/g fw", "≥ 2.4 µg/g fw", "< 2.4 µg/g fw"))) %>%
  filter(n() > 4) %>% 
  # Trying to create some sort of weighted average for the ordering of the y axis - this was sooo painstaking
  mutate(Ext = ifelse(Risk_Levels == "≥ 5.3 µg/g fw", 1, 0),
         Prop_Ext = sum(Ext)/n(),
         High = ifelse(Risk_Levels == "≥ 4.5 µg/g fw", 1, 0),
         Prop_High = sum(High)/n(),
         Med = ifelse(Risk_Levels == "≥ 3.4 µg/g fw", 1, 0),
         Prop_Med = sum(Med)/n(),
         Low = ifelse(Risk_Levels == "≥ 2.4 µg/g fw", 1, 0),
         Prop_Low = sum(Low)/n(),
         None = ifelse(Risk_Levels == "< 2.4 µg/g fw", 1, 0),
         Prop_None = sum(None)/n(),
         Score = sum(Ext, High, Med, Low, None)) %>% 
  ggplot(mapping = aes(y = reorder(reorder(reorder(reorder(
    reorder(reorder(reorder(Species_Latin_Name, desc(Species_Latin_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightgray")) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)", fill = "Body feather\nrisk categories") +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_markdown(), # this is the key step for italic species names
        legend.title = element_text(face = "bold"),
        #legend.position = "none",
        aspect.ratio = 1)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 12, height = 8)


# Risk graph for rectrices
tailrisk <- CollectiveData %>%
  filter(!is.na(Species_Latin_Name), !is.na(Tail_Hg_ppm)) %>%
  group_by(Species_Latin_Name) %>%
  mutate(Sample_Size = str_c("n = ", n()),
         Species_Latin_Name = glue("*{Species_Latin_Name}* ({Sample_Size})"),
         Risk_Levels = cut(Tail_Hg_ppm, breaks = c(0, 3, 4.7, 6.4, 7.7, Inf),
                           labels = c("< 3 µg/g fw", "≥ 3 µg/g fw", "≥ 4.7 µg/g fw", "≥ 6.4 µg/g fw", "≥ 7.7 µg/g fw")),
         # reordering risk categories
         Risk_Levels = factor(Risk_Levels, levels = c("≥ 7.7 µg/g fw", "≥ 6.4 µg/g fw", "≥ 4.7 µg/g fw", "≥ 3 µg/g fw", "< 3 µg/g fw"))) %>%
  filter(n() > 4) %>% 
  # Trying to create some sort of weighted average for the ordering of the y axis - this was sooo painstaking
  mutate(Ext = ifelse(Risk_Levels == "≥ 7.7 µg/g fw", 1, 0),
         Prop_Ext = sum(Ext)/n(),
         High = ifelse(Risk_Levels == "≥ 6.4 µg/g fw", 1, 0),
         Prop_High = sum(High)/n(),
         Med = ifelse(Risk_Levels == "≥ 4.7 µg/g fw", 1, 0),
         Prop_Med = sum(Med)/n(),
         Low = ifelse(Risk_Levels == "≥ 3 µg/g fw", 1, 0),
         Prop_Low = sum(Low)/n(),
         None = ifelse(Risk_Levels == "< 3 µg/g fw", 1, 0),
         Prop_None = sum(None)/n(),
         Score = sum(Ext, High, Med, Low, None)) %>% 
  ggplot(mapping = aes(y = reorder(reorder(reorder(reorder(
    reorder(reorder(reorder(Species_Latin_Name, desc(Species_Latin_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightgray")) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)", fill = "Tail feather\nrisk categories") +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_markdown(), # this is the key step for italic species names
        legend.title = element_text(face = "bold"),
        #legend.position = "none",
        aspect.ratio = 1)

ggview(device = "jpeg", units = "in", dpi = 1200, width = 12, height = 8)


ggarrange(bloodrisk, bodyrisk, tailrisk, labels = c("a)", "b)", "c)"), nrow = 3, ncol = 1)
ggview(device = "jpeg", units = "in", dpi = 1200, width = 12, height = 24)

ggsave("Publication-Figures/Fig4_Facet_Risk_AllTissues.jpg", dpi = 1200, width = 12, height = 24)
#ggsave("Publication-Figures/Fig4_Facet_Risk_AllTissues.tiff", dpi = 600, width = 12, height = 24)


# HG MAPS ------------------------------------------------------------------
library(gcookbook)
library(sp)
library(raster)

# Google satellite imagery as a background
library(ggmap)
#register_google(key = "Insert personal key here", write = TRUE)

# For google map, you have to give the center of the window you are looking at.
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid
# get the map info
satellitemap <- get_map(location = c(lon = -83.7, lat = 3.1), maptype = "satellite", source = "google", zoom = 4)

# Blood map with points as circles
MapData <- CollectiveData %>%
  filter(!is.na(Blood_Hg_ppm)) %>% 
  arrange(Blood_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

# setting breaks for the color ramp
my_breaks = c(0.001, 0.01, 0.1, 1, 3)

bloodmap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat, color = Blood_Hg_ppm, size = Blood_Hg_ppm)) +
  #scale_color_viridis_c() +
  scale_color_gradient(high = "red", low = "yellow", trans = "log", breaks = my_breaks,
                       labels = my_breaks) +
  scale_size(range = c(10, 2), trans = "log") + # setting the point size
  #labs(x = "Longitude", y = "Latitude", color = "Whole blood\nTHg (µg/g ww)") + 
  labs(color = "Whole blood\nTHg (µg/g ww)") +
  theme_classic(base_size = 14) +
  theme(#axis.title.x = element_text(face = "bold"),
    #axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1) +
  guides(size = "none") # removing size legend 

ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 6)


# Body feather map with points
MapData <- CollectiveData %>%
  filter(!is.na(Body_Hg_ppm)) %>% 
  arrange(Body_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

# setting breaks for the color ramp
my_breaks = c(0.001, 0.01, 0.1, 1, 10, 60)

bodymap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Banding_Station_Long, y = Banding_Station_Lat,
                                           color = Body_Hg_ppm, size = Body_Hg_ppm)) +
  scale_color_gradient(high = "red", low = "yellow", trans = "log", breaks = my_breaks,
                       labels = my_breaks) +
  scale_size(range = c(10,2), trans = "log") + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "Body feather\nTHg (µg/g fw)") + 
  #ggtitle("Resident species only") +
  theme_classic(base_size = 14) +
  theme(#axis.title.x = element_text(face = "bold"),
    #axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1) +
  guides(size = "none") # removing size legend 

ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 6)


# Rectrice map with points
MapData <- CollectiveData %>%
  filter(!is.na(Tail_Hg_ppm)) %>% 
  arrange(Tail_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

# setting breaks for the color ramp
my_breaks = c(0.001, 0.01, 0.1, 1, 10)
#breaks = c(0, 2, 6, 10)

tailmap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Banding_Station_Long, y = Banding_Station_Lat,
                                           color = Tail_Hg_ppm, size = Tail_Hg_ppm)) +
  scale_color_gradient(high = "red", low = "yellow", trans = "log", breaks = my_breaks,
                       labels = my_breaks) +
  scale_size(range = c(10,2), trans = "log") + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "Tail feather\nTHg (µg/g fw)") + 
  theme_classic(base_size = 14) +
  theme(#axis.title.x = element_text(face = "bold"),
    #axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1) +
  guides(size = "none") # removing size legend

ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 6)

ggarrange(bloodmap, bodymap, tailmap, labels = c("A)", "B)", "C)"), nrow = 3, ncol = 1)
ggview(device = "jpeg", units = "in", dpi = 1200, width = 8.5, height = 18)


# SAMPLE MAPS ------------------------------------------------------------------
library(devtools)
# VERY important to install the older version of the ddgridR package
#install_version("dggridR", version = "2.0.4", repos = "http://cran.us.r-project.org")
library(dggridR)
library(mapdata)

# filtering specific samples
MapData <- CollectiveData %>% 
  filter(!is.na(Blood_Hg_ppm), !is.na(Banding_Station_Long),
         !is.na(Banding_Station_Lat))

# constructing hex grid of a certain size
hexgrid <- dggridR::dgconstruct(res = 7)

MapData$cell <- dgGEO_to_SEQNUM(hexgrid, MapData$Banding_Station_Long, MapData$Banding_Station_Lat)$seqnum

cell_count <- MapData %>%
  group_by(cell) %>%
  summarise(count = n())

grid <- dgcellstogrid(hexgrid, cell_count$cell, frame = TRUE, wrapcells = TRUE)
grid <- merge(grid, cell_count, by = "cell")

bloodsamplemap <- ggmap(satellitemap) +
  #coord_cartesian() +
  geom_polygon(data = grid, aes(x = long, y = lat, group = group, fill = count)) +
  geom_path(data = grid, aes(x = long, y = lat, group = group), alpha = 0.8, color = "white") +
  #geom_hex(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #geom_bin2d(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(high = "red", low = "yellow", breaks = c(100, 200, 300, 400)) +
  labs(fill = "Whole blood\nsample size") +
  theme_classic(base_size = 14) +
  theme(#axis.title.x = element_text(face = "bold"),
    #axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1) +
  guides(size = "none") # removing size legend 

ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 6)


# Body feather sample map
MapData <- CollectiveData %>%
  filter(!is.na(Body_Hg_ppm), !is.na(Banding_Station_Long),
         !is.na(Banding_Station_Lat))

# constructing hex grid of a certain size
hexgrid <- dggridR::dgconstruct(res = 7)

MapData$cell <- dgGEO_to_SEQNUM(hexgrid, MapData$Banding_Station_Long, MapData$Banding_Station_Lat)$seqnum

cell_count <- MapData %>%
  group_by(cell) %>%
  summarise(count = n())

grid <- dgcellstogrid(hexgrid, cell_count$cell, frame = TRUE, wrapcells = TRUE)
grid <- merge(grid, cell_count, by = "cell")

bodysamplemap <- ggmap(satellitemap) +
  #coord_cartesian() +
  geom_polygon(data = grid, aes(x = long, y = lat, group = group, fill = count)) +
  geom_path(data = grid, aes(x = long, y = lat, group = group), alpha = 0.8, color = "white") +
  #geom_hex(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #geom_bin2d(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(high = "red", low = "yellow", breaks = c(50, 100, 150, 200, 250)) +
  labs(x = "Longitude", y = "Latitude", fill = "Body feather\nsample size") + 
  #ggtitle("Resident species only") +
  theme_classic(base_size = 14) +
  theme(#axis.title.x = element_text(face = "bold"),
    #axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1) +
  guides(size = "none") # removing size legend 

ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 6)


# Tail feather sample map
MapData <- CollectiveData %>%
  filter(!is.na(Tail_Hg_ppm), !is.na(Banding_Station_Long),
         !is.na(Banding_Station_Lat))

# constructing hex grid of a certain size
hexgrid <- dggridR::dgconstruct(res = 7)

MapData$cell <- dgGEO_to_SEQNUM(hexgrid, MapData$Banding_Station_Long, MapData$Banding_Station_Lat)$seqnum

cell_count <- MapData %>%
  group_by(cell) %>%
  summarise(count = n())

grid <- dgcellstogrid(hexgrid, cell_count$cell, frame = TRUE, wrapcells = TRUE)
grid <- merge(grid, cell_count, by = "cell")

tailsamplemap <- ggmap(satellitemap) +
  #coord_cartesian() +
  geom_polygon(data = grid, aes(x = long, y = lat, group = group, fill = count)) +
  geom_path(data = grid, aes(x = long, y = lat, group = group), alpha = 0.8, color = "white") +
  #geom_hex(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #geom_bin2d(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(high = "red", low = "yellow",breaks = c(50, 100, 150, 200, 250)) +
  labs(x = "Longitude", y = "Latitude", fill = "Tail feather\nsample size") +
  #ggtitle("Resident species only") +
  theme_classic(base_size = 14) +
  theme(#axis.title.x = element_text(face = "bold"),
    #axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1) +
  guides(size = "none") # removing size legend

ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 6)

ggarrange(bloodsamplemap, bodysamplemap, tailsamplemap, labels = c("a)", "b)", "C)"), nrow = 3, ncol = 1)
ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 18)


ggarrange(bloodmap, bloodsamplemap, bodymap, bodysamplemap, tailmap, tailsamplemap,
          labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
          nrow = 3, ncol = 2)
ggview(device = "jpeg", units = "in", dpi = 1200, width = 16.5, height = 18)

ggsave("Publication-Figures/Fig5_Facet_Map_AllTissueConcSample.jpg", dpi = 1200, width = 16.5, height = 18)
#ggsave("Publication-Figures/Fig5_Facet_Map_AllTissueConcSample.tiff", dpi = 600, width = 16.5, height = 18)


# Graphical abstract Hg
# All tissue map with points
MapData <- CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Concentration)) %>%
  arrange(Concentration) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

# setting breaks for the color ramp
my_breaks = c(0.001, 0.01, 0.1, 1, 10, 60)

GAbstractHg <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Banding_Station_Long, y = Banding_Station_Lat,
                                           color = Concentration, size = Concentration),
             position = position_jitter(width = 0.5, height = 0.5)) +
  scale_color_gradient(high = "red", low = "yellow", trans = "log", breaks = my_breaks,
                       labels = my_breaks) +
  scale_size(range = c(10, 2), trans = "log") + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "THg (µg/g)") +
  theme_classic(base_size = 14) +
  theme(#axis.title.x = element_text(face = "bold"),
    #axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1) +
  guides(size = "none") # removing size legend 

ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 6)

ggsave("Publication-Figures/GAbstractHg.jpg", dpi = 1200, width = 8, height = 6)


# Graphical abstract sample map
MapData <- CollectiveData %>%
  filter(!is.na(Banding_Station_Long),
         !is.na(Banding_Station_Lat))

# constructing hex grid of a certain size
hexgrid <- dggridR::dgconstruct(res = 7)

MapData$cell <- dgGEO_to_SEQNUM(hexgrid, MapData$Banding_Station_Long, MapData$Banding_Station_Lat)$seqnum

cell_count <- MapData %>%
  group_by(cell) %>%
  summarise(count = n())

grid <- dgcellstogrid(hexgrid, cell_count$cell, frame = TRUE, wrapcells = TRUE)
grid <- merge(grid, cell_count, by = "cell")

GAbstractN <- ggmap(satellitemap) +
  #coord_cartesian() +
  geom_polygon(data = grid, aes(x = long, y = lat, group = group, fill = count)) +
  geom_path(data = grid, aes(x = long, y = lat, group = group), alpha = 0.8, color = "white") +
  #geom_hex(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #geom_bin2d(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(high = "red", low = "yellow", breaks = c(200, 400, 600)) +
  labs(x = "Longitude", y = "Latitude", fill = "Sample size") +
  #ggtitle("Resident species only") +
  theme_classic(base_size = 14) +
  theme(#axis.title.x = element_text(face = "bold"),
    #axis.title.y = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1) +
  guides(size = "none") # removing size legend

ggview(device = "jpeg", units = "in", dpi = 1200, width = 8, height = 6)


ggarrange(GAbstractHg, GAbstractN, nrow = 1, ncol = 2)
ggview(device = "jpeg", units = "in", dpi = 1200, width = 16, height = 6)

ggsave("Publication-Figures/GAbstract_Satellite_Circle_AllTissue.jpg", dpi = 1200, width = 16, height = 6)
#ggsave("Publication-Figures/GAbstract_Satellite_Circle_AllTissue.tiff", dpi = 600, width = 16, height = 6)

