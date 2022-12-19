
# Script designed to produce graphs featured in publication and annex
# Requires objects from Sayers_Main_Script.R

library(tidyverse)
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(scales)

library(lmerTest)
library(ggeffects)
library(multcomp)

library(ggtext)
library(glue)

# resolve namespace conflicts and creating necessary functions
select <- dplyr::select
"%nin%" <- Negate("%in%")

# 1st place model from model selection
topmodel <- lmer(log(Hg_Concentration) ~ Trophic_Niche + Tissue_Type + Mining_Present_Yes_No +
                   (1 | SiteID) + (1 | Family/Species_Common_Name) + (1 | Date),
                 data = HgSamples, REML = F)

# PREDICTED TROPHIC NICHE -----------------------------------------------------------

# predicted mean points overlayed with raw data

# calculating tissue sample sizes for y axis
ss <- HgSamples %>%
  group_by(Trophic_Niche, Tissue_Type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Tissue_Type, values_from = n, values_fill = 0) %>%
  select(Trophic_Niche, n_Blood = Blood_Hg_ppm, n_Body = Body_Hg_ppm, n_Tail = Tail_Hg_ppm) %>% 
  mutate(Trophic_Niche.s = str_c(Trophic_Niche, " (n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"),
         n_Total = sum(n_Blood, n_Body, n_Tail))

# calculating predicted Hg values with top model structure
# type = "random" gives prediction intervals rather than confidence intervals
pr <- ggpredict(topmodel, terms = c("Trophic_Niche", "Tissue_Type", "Mining_Present_Yes_No"),
                type = "random", back.transform = T) %>% 
  rename(Trophic_Niche = x, Tissue_Type = group, Mining_Present_Yes_No = facet) %>%
  left_join(ss, by = "Trophic_Niche") %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("Yes", "No"),
                                           labels = c("ASGM present", "ASGM absent")))

# create final data frame with raw data and predicted means to plot
df <- HgSamples %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                           labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("Yes", "No"),
                                           labels = c("ASGM present", "ASGM absent"))) %>% 
  full_join(pr, by = c("Trophic_Niche", "Tissue_Type", "Mining_Present_Yes_No"))


ggplot() +
  # putting a dummy geom in front to get the correct order
  geom_point(data = pr, mapping = aes(x = predicted, y = reorder(Trophic_Niche.s, predicted)),
             color = "white", size = 0.01) +
  geom_point(data = df, mapping = aes(x = Hg_Concentration, y = reorder(Trophic_Niche.s, predicted),
                                      color = Tissue_Type, fill = Mining_Present_Yes_No),
             position = position_jitterdodge(jitter.height = 0.1, dodge.width = 0.6),
             alpha = 0.5, size = 2) +
  geom_errorbar(data = pr, mapping = aes(xmin = conf.low, xmax = conf.high,
                                         y = reorder(Trophic_Niche.s, predicted),
                                         shape = Mining_Present_Yes_No), color = "black",
                position = position_dodge(0.6), width = 0) +
  geom_point(data = pr, mapping = aes(x = predicted, y = reorder(Trophic_Niche.s, predicted),
                           shape = Mining_Present_Yes_No), color = "black",
             position = position_dodge(0.6), size = 3) +
  #geom_pointrange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(0.9)) +
  labs(x = "Predicted THg (µg/g)", y = "Trophic niche") +
  scale_x_continuous(expand = c(0.1, 0)
                     , trans = "log",
                     breaks = c(0, 0.01, 0.1, 1, 10, 60),
                     labels = c("0", "0.01", "0.1", "1", "10", "60")) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), guide = "none") +
  scale_shape_manual(limits = c("ASGM present", "ASGM absent"), values = c(17, 16), 
                     guide = "none") +
  scale_fill_discrete(guide = "none") +
  theme_minimal() +
  facet_grid(Mining_Present_Yes_No ~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"))

ggsave("FinalGraphs/Fig1_Predicted_Trophic_AllTissues_Hist.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Fig1_Predicted_Trophic_AllTissues_Hist.tiff", dpi = 1000, width = 10, height = 8)

# PREDICTED FAMILY -----------------------------------------------------------------

# predicted mean points overlayed with raw data

# calculating tissue sample sizes for y axis
ss <- HgSamples %>%
  group_by(Family, Tissue_Type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Tissue_Type, values_from = n, values_fill = 0) %>%
  select(Family, n_Blood = Blood_Hg_ppm, n_Body = Body_Hg_ppm, n_Tail = Tail_Hg_ppm) %>% 
  mutate(Family.s = str_c(Family, " (n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"),
         n_Total = sum(n_Blood, n_Body, n_Tail))

# calculating predicted Hg values with top model structure
# type = "random" gives prediction intervals rather than confidence intervals
pr <- ggpredict(topmodel, terms = c("Family", "Tissue_Type", "Mining_Present_Yes_No"),
                type = "random", back.transform = T) %>% 
  rename(Family = x, Tissue_Type = group, Mining_Present_Yes_No = facet) %>%
  left_join(ss, by = "Family") %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("Yes", "No"),
                                           labels = c("ASGM present", "ASGM absent"))) %>%
  filter(n_Total > 9)

# create final data frame with raw data and predicted means to plot
df <- HgSamples %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  transform(Mining_Present_Yes_No = factor(Mining_Present_Yes_No, levels = c("Yes", "No"),
                                           labels = c("ASGM present", "ASGM absent"))) %>% 
  full_join(pr, by = c("Family", "Tissue_Type", "Mining_Present_Yes_No")) %>%
  filter(n_Total > 9)


ggplot() +
  # putting a dummy geom in front to get the correct order
  geom_point(data = pr, mapping = aes(x = predicted, y = reorder(Family.s, predicted)),
             color = "white", size = 0.01) +
  geom_point(data = df, mapping = aes(x = Hg_Concentration, y = reorder(Family.s, predicted),
                                      color = Tissue_Type, fill = Mining_Present_Yes_No),
             position = position_jitterdodge(jitter.height = 0.1, dodge.width = 0.9),
             alpha = 0.5, size = 1) +
  geom_errorbar(data = pr, mapping = aes(xmin = conf.low, xmax = conf.high,
                                         y = reorder(Family.s, predicted),
                                        shape = Mining_Present_Yes_No), color = "black",
                position = position_dodge(0.9), width = 0) +
  geom_point(data = pr, mapping = aes(x = predicted, y = reorder(Family.s, predicted),
                                      shape = Mining_Present_Yes_No), color = "black",
             position = position_dodge(0.9), size = 2) +
  facet_grid(Mining_Present_Yes_No ~ Tissue_Type, scales = "free") +
  labs(x = "Predicted THg (µg/g)", y = "Family (n ≥ 10)") +
  scale_x_continuous(expand = c(0.1, 0), trans = "log",
                     breaks = c(0, 0.01, 0.1, 1, 10, 60),
                     labels = c("0", "0.01", "0.1", "1", "10", "60")) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), guide = "none") +
  scale_shape_manual(limits = c("ASGM present", "ASGM absent"), values = c(17, 16),
                     guide = "none") +
  scale_fill_discrete(guide = "none") +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"))

ggsave("FinalGraphs/Fig2_Predicted_Family_AllTissues_Hist.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Fig2_Predicted_Family_AllTissues_Hist.tiff", dpi = 1000, width = 10, height = 8)


# PREDICTED SITE --------------------------------------------------------------------

# predicted mean points overlayed with raw data

# calculating tissue sample sizes for y axis
ss <- HgSamples %>%
  mutate(Site_Name = if_else(Site_Name == "El Portal", str_c("El Yunque National Forest", "--", Site_Name), Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Río Icacos", str_c("El Yunque National Forest", "--", Site_Name), Site_Name)) %>% 
  mutate(Site_Name = if_else(Site_Name == "Torre Britton", str_c("El Yunque National Forest", "--", Site_Name), Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Azul Mine", "Mining Site A", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "La Torre", "Reserva Nacional Tambopata", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Paolita Mine", "Mining Site B", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Santa Rita Mine", "Mining Site C", Site_Name)) %>%
  mutate(Site_Name = str_c(Site_Name, ", ", Country)) %>%
  group_by(SiteID, Site_Name, Tissue_Type) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Tissue_Type, values_from = n, values_fill = 0) %>%
  select(SiteID, Site_Name, n_Blood = Blood_Hg_ppm, n_Body = Body_Hg_ppm, n_Tail = Tail_Hg_ppm) %>% 
  mutate(Site_Name.s = str_c(Site_Name, " (n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")"),
         n_Total = sum(n_Blood, n_Body, n_Tail))

# calculating predicted Hg values with top model structure
# type = "random" gives prediction intervals rather than confidence intervals
pr <- ggpredict(topmodel, terms = c("SiteID", "Tissue_Type"),
                type = "random", back.transform = T) %>% 
  rename(SiteID = x, Tissue_Type = group) %>%
  left_join(ss, by = "SiteID") %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>% 
  filter(n_Total > 24)

# create final data frame with raw data and predicted means to plot
df <- HgSamples %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type, levels = c("Blood_Hg_ppm", "Body_Hg_ppm", "Tail_Hg_ppm"),
                                 labels = c("Whole blood", "Body feather", "Tail feather"))) %>%
  left_join(pr, by = c("SiteID", "Tissue_Type")) %>%
  filter(n_Total > 24)

ggplot() +
  # putting a dummy geom in front to get the correct order
  geom_point(data = pr, mapping = aes(x = predicted, y = reorder(Site_Name.s, predicted)),
             color = "white", size = 0.01) +
  geom_point(data = df, mapping = aes(x = Hg_Concentration, y = reorder(Site_Name.s, predicted), 
             color = Tissue_Type),
             position = position_jitter(height = 0.1), alpha = 0.6, size = 2) +
  geom_errorbar(data = pr, mapping = aes(xmin = conf.low, xmax = conf.high,
                                         y = reorder(Site_Name.s, predicted)),
                color = "black", width = 0) +
  geom_point(data = pr, mapping = aes(x = predicted, y = reorder(Site_Name.s, predicted)),
             color = "black", size = 3) +
  labs(x = "Predicted THg (µg/g)", y = "Site (n ≥ 25)") +
  scale_x_continuous(expand = c(0.1, 0), trans = "log",
                     breaks = c(0, 0.01, 0.1, 1, 10, 60),
                     labels = c("0", "0.01", "0.1", "1", "10", "60")
                     ) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), guide = "none") +
  theme_minimal() +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"))

ggsave("FinalGraphs/Fig3_Predicted_Site_AllTissues_Hist.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Fig3_Predicted_Site_AllTissues_Hist.tiff", dpi = 1000, width = 10, height = 8)


# TROPHIC NICHE x PRIMARY HABITAT ---------------------------------------------

# columns of raw data arithmetic means ±95%CI

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Trophic_Niche), !is.na(HAB1), !is.na(Concentration)) %>%
  pivot_longer(c(Trophic_Niche, HAB1), names_to = "Category", values_to = "Type") %>%
  group_by(Category, Type, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration)) %>%
  data.table() %>% 
  data.table::dcast(Category + Type ~ Tissue_Type, value.var = c("n", "mean", "sd"),
                    fill = list(n = 0, mean = NA, sd = NA)) %>%
  select(Category, Type, n_Blood = n_Blood_Hg_ppm, n_Body = n_Body_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         mean_Blood = mean_Blood_Hg_ppm, mean_Body = mean_Body_Hg_ppm, mean_Tail = mean_Tail_Hg_ppm,
         sd_Blood = sd_Blood_Hg_ppm, sd_Body = sd_Body_Hg_ppm, sd_Tail = sd_Tail_Hg_ppm) %>% 
  mutate(Type = str_c(Type, " (n = ", n_Blood, ", ", n_Body, ", ", n_Tail, ")")) %>% 
  # renaming category values for later plotting
  mutate(Category = if_else(Category == "Trophic_Niche", "Trophic niche", "Primary habitat"))

p1 <- df %>%
  select(Category, Type, mean = mean_Blood, sd = sd_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Category, Type, mean = mean_Body, sd = sd_Body, n_Body = n_Body) %>%
  mutate(Tissue_Type = "Body feather")
p3 <- df %>%
  select(Category, Type, mean = mean_Tail, sd = sd_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(mean != 0) %>%
  # calculating percentile rank for each tissue
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>%
  group_by(Category, Type) %>% 
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather"))) %>% 
  transform(Category = factor(Category,
                              levels = c("Trophic niche", "Primary habitat")))

ggplot(final, mapping = aes(x = mean, y = reorder(Type, max_percentile),
                            fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0.5*mean, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.2) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_minimal() +
  facet_grid(Category ~ Tissue_Type, scales = "free", switch = "y") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold", size = 11),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"))

ggsave("FinalGraphs/Sayers_TrophicxHabitat_AllTissues_Hist.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_TrophicxHabitat_AllTissues_Hist.tiff", dpi = 1000, width = 10, height = 8)


# ORDER VARIATION ---------------------------------------------------------

# columns of raw data arithmetic means ±95%CI

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  # adding tissue type as a data field 
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Order), !is.na(Concentration)) %>%
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
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Order) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather")))
  
ggplot(final, mapping = aes(x = mean, y = reorder(Order, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Order") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_minimal() +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"))

ggsave("FinalGraphs/Sayers_Order_AllTissues_Hist.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Order_AllTissues_Hist.tiff", dpi = 1000, width = 10, height = 8)


# FAMILY VARIATION --------------------------------------------------------

#columns of raw data arithmetic means ±95%CI

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
  filter(sum > 9) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Family) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather")))

ggplot(final, mapping = aes(x = mean, y = reorder(Family, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0.5*mean, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Family (n ≥ 10)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_minimal() +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"))

ggsave("FinalGraphs/Sayers_Family_AllTissues_Hist.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Family_AllTissues_Hist.tiff", dpi = 1000, width = 10, height = 8)


# SPECIES VARIATION --------------------------------------------------------

# columns of raw data arithmetic means ±95%CI

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
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Species_Latin_Name) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather")))

ggplot(final, mapping = aes(x = mean, y = reorder(Species_Latin_Name, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0.5*mean, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Species (n ≥ 10)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_minimal() +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_markdown(), # this is the key step for italic species names
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"))

ggsave("FinalGraphs/Sayers_Species_AllTissues_Hist.jpg", dpi = 1000, width = 8, height = 10)
ggsave("FinalGraphs/Sayers_Species_AllTissues_Hist.tiff", dpi = 1000, width = 8, height = 10)




# SITE VARIATION ---------------------------------------------------------

# columns of raw data arithmetic means ±95%CI

# calculating tissue sample sizes for y axis
df <- CollectiveData %>% 
  mutate(Site_Name = if_else(Site_Name == "El Portal", str_c("El Yunque National Forest", "--", Site_Name), Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Río Icacos", str_c("El Yunque National Forest", "--", Site_Name), Site_Name)) %>% 
  mutate(Site_Name = if_else(Site_Name == "Torre Britton", str_c("El Yunque National Forest", "--", Site_Name), Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Azul Mine", "Mining Site A", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "La Torre", "Reserva Nacional Tambopata", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Paolita Mine", "Mining Site B", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Santa Rita Mine", "Mining Site C", Site_Name)) %>%
  mutate(Site_Name = str_c(Site_Name, ", ", Country)) %>%
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
  filter(n_Total > 24)

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
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = percent_rank(mean)) %>% 
  group_by(Site_Name) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Body feather", "Tail feather"))) 

ggplot(final, mapping = aes(x = mean, y = reorder(Site_Name, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0.5*mean, xmax = mean + sd),
                position = position_dodge(0.9), width = 0.15) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = "Site (n ≥ 25)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
  theme_minimal() +
  facet_grid(~ Tissue_Type, scales = "free") +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"))

ggsave("FinalGraphs/Sayers_Site_AllTissues_Hist.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Site_AllTissues_Hist.tiff", dpi = 1000, width = 10, height = 8)


# RISK ASSESSMENT ---------------------------------------------------------

# Risk graph for blood
bloodrisk <- CollectiveData %>%
  filter(!is.na(Blood_Hg_ppm), !is.na(Species_Latin_Name),
         Species_Code != "BIRD") %>%
  group_by(Species_Latin_Name) %>%
  mutate(Sample_Size = str_c("n = ", n()),
         Species_Latin_Name = glue("*{Species_Latin_Name}* ({Sample_Size})"),
         Risk_Levels = cut(Blood_Hg_ppm, breaks = c(0, 0.7, 1.2, 1.7, 2.2, Inf),
                           labels = c("<0.7 µg/g ww", "≥0.7 µg/g ww", "≥1.2 µg/g ww", "≥1.7 µg/g ww", "≥2.2 µg/g ww")),
         # reordering risk categories
         Risk_Levels = factor(Risk_Levels, levels = c("≥2.2 µg/g ww","≥1.7 µg/g ww", "≥1.2 µg/g ww", "≥0.7 µg/g ww", "<0.7 µg/g ww"))) %>%
  filter(n() > 4) %>%
  # Trying to create some sort of weighted average for the ordering of the y axis - this was sooo painstaking
  mutate(Ext = ifelse(Risk_Levels == "≥2.2 µg/g ww", 1, 0),
         Prop_Ext = sum(Ext)/n(),
         High = ifelse(Risk_Levels == "≥1.7 µg/g ww", 1, 0),
         Prop_High = sum(High)/n(),
         Med = ifelse(Risk_Levels == "≥1.2 µg/g ww", 1, 0),
         Prop_Med = sum(Med)/n(),
         Low = ifelse(Risk_Levels == "≥0.7 µg/g ww", 1, 0),
         Prop_Low = sum(Low)/n(),
         None = ifelse(Risk_Levels == "<0.7 µg/g ww", 1, 0),
         Prop_None = sum(None)/n(),
         Score = sum(Ext, High, Med, Low, None)) %>%
  ggplot(mapping = aes(y = reorder(reorder(reorder(reorder(
    reorder(reorder(reorder(Species_Latin_Name, desc(Species_Latin_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightgray")) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)",
       fill = "Whole blood\nrisk categories") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_markdown(), # this is the key step for italic species names
        legend.title = element_text(face = "bold"),
        #legend.position = "none",
        aspect.ratio = 1)

ggsave("FinalGraphs/Sayers_Risk_Blood.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Risk_Blood.tiff", dpi = 1000, width = 10, height = 8)

# Risk graph for body feathers
bodyrisk <- CollectiveData %>%
  filter(!is.na(Species_Latin_Name), !is.na(Body_Hg_ppm), Species_Code != "BIRD") %>%
  group_by(Species_Latin_Name) %>%
  mutate(Sample_Size = str_c("n = ", n()),
         Species_Latin_Name = glue("*{Species_Latin_Name}* ({Sample_Size})"),
         Risk_Levels = cut(Body_Hg_ppm, breaks = c(0, 2.4, 3.4, 4.5, 5.3, Inf),
                           labels = c("<2.4 µg/g fw", "≥2.4 µg/g fw", "≥3.4 µg/g fw", "≥4.5 µg/g fw", "≥5.3 µg/g fw")),
         # reordering risk categories
         Risk_Levels = factor(Risk_Levels, levels = c("≥5.3 µg/g fw", "≥4.5 µg/g fw",
                                                      "≥3.4 µg/g fw", "≥2.4 µg/g fw", "<2.4 µg/g fw"))) %>%
  filter(n() > 4) %>% 
  # Trying to create some sort of weighted average for the ordering of the y axis - this was sooo painstaking
  mutate(Ext = ifelse(Risk_Levels == "≥5.3 µg/g fw", 1, 0),
         Prop_Ext = sum(Ext)/n(),
         High = ifelse(Risk_Levels == "≥4.5 µg/g fw", 1, 0),
         Prop_High = sum(High)/n(),
         Med = ifelse(Risk_Levels == "≥3.4 µg/g fw", 1, 0),
         Prop_Med = sum(Med)/n(),
         Low = ifelse(Risk_Levels == "≥2.4 µg/g fw", 1, 0),
         Prop_Low = sum(Low)/n(),
         None = ifelse(Risk_Levels == "<2.4 µg/g fw", 1, 0),
         Prop_None = sum(None)/n(),
         Score = sum(Ext, High, Med, Low, None)) %>% 
  ggplot(mapping = aes(y = reorder(reorder(reorder(reorder(
    reorder(reorder(reorder(Species_Latin_Name, desc(Species_Latin_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightgray")) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)", fill = "Body feather\nrisk categories") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_markdown(), # this is the key step for italic species names
        legend.title = element_text(face = "bold"),
        #legend.position = "none",
        aspect.ratio = 1)

ggsave("FinalGraphs/Sayers_Risk_Body.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Risk_Body.tiff", dpi = 1000, width = 10, height = 8)

# Risk graph for rectrices
tailrisk <- CollectiveData %>%
  filter(!is.na(Species_Latin_Name), !is.na(Tail_Hg_ppm), Species_Code != "BIRD") %>%
  group_by(Species_Latin_Name) %>%
  mutate(Sample_Size = str_c("n = ", n()),
         Species_Latin_Name = glue("*{Species_Latin_Name}* ({Sample_Size})"),
         Risk_Levels = cut(Tail_Hg_ppm, breaks = c(0, 3, 4.7, 6.4, 7.7, Inf),
                           labels = c("<3 µg/g fw", "≥3 µg/g fw", "≥4.7 µg/g fw", "≥6.4 µg/g fw", "≥7.7 µg/g fw")),
         # reordering risk categories
         Risk_Levels = factor(Risk_Levels, levels = c("≥7.7 µg/g fw", "≥6.4 µg/g fw", "≥4.7 µg/g fw", "≥3 µg/g fw", "<3 µg/g fw"))) %>%
  filter(n() > 4) %>% 
  # Trying to create some sort of weighted average for the ordering of the y axis - this was sooo painstaking
  mutate(Ext = ifelse(Risk_Levels == "≥7.7 µg/g fw", 1, 0),
         Prop_Ext = sum(Ext)/n(),
         High = ifelse(Risk_Levels == "≥6.4 µg/g fw", 1, 0),
         Prop_High = sum(High)/n(),
         Med = ifelse(Risk_Levels == "≥4.7 µg/g fw", 1, 0),
         Prop_Med = sum(Med)/n(),
         Low = ifelse(Risk_Levels == "≥3 µg/g fw", 1, 0),
         Prop_Low = sum(Low)/n(),
         None = ifelse(Risk_Levels == "<3 µg/g fw", 1, 0),
         Prop_None = sum(None)/n(),
         Score = sum(Ext, High, Med, Low, None)) %>% 
  ggplot(mapping = aes(y = reorder(reorder(reorder(reorder(
    reorder(reorder(reorder(Species_Latin_Name, desc(Species_Latin_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightgray")) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)", fill = "Tail feather\nrisk categories") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_markdown(), # this is the key step for italic species names
        legend.title = element_text(face = "bold"),
        #legend.position = "none",
        aspect.ratio = 1)

ggsave("FinalGraphs/Sayers_Risk_Tail.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Risk_Tail.tiff", dpi = 1000, width = 10, height = 8)


ggarrange(bloodrisk, bodyrisk, tailrisk, labels = c("a)", "b)", "c)"), nrow = 3, ncol = 1)
ggsave("FinalGraphs/Fig4_Facet_Risk_AllTissues.jpg", dpi = 1000, width = 10, height = 15)
ggsave("FinalGraphs/Fig4_Facet_Risk_AllTissues.tiff", dpi = 800, width = 10, height = 15)


# MAPS ------------------------------------------------------------------
library(gcookbook)
library(sp)
library(raster)

# Google satellite imagery as a background
library(ggmap)
register_google(key = "AIzaSyAhVAGnjiPZgt0KtXSO_Od2j67CF6wEmD8", write = TRUE) # Chris' personal key

# Graphical abstract
# All tissue map with points
MapData <- CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Concentration)) %>%
  arrange(Concentration) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

satellitemap <- get_map(location = c(lon = -83, lat = 5), maptype = "satellite", source = "google", zoom = 4)

# setting breaks for the color ramp
my_breaks = c(0.01, 0.1, 1, 10, 60)

ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Banding_Station_Long, y = Banding_Station_Lat,
                                           color = Concentration, size = Concentration),
             position = position_jitter(width = 0.5, height = 0.5)) +
  scale_color_gradient(high = "red", low = "yellow", trans = "log", breaks = my_breaks,
                       labels = my_breaks) +
  scale_size(range = c(10, 2), trans = "log") + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "THg (µg/g)") +
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

ggsave("FinalGraphs/GAbstract_Satellite_Circle_AllTissue.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/GAbstract_Satellite_Circle_AllTissue.tiff", dpi = 1000, width = 10, height = 8)




# Blood map with points as circles
MapData <- CollectiveData %>%
  filter(!is.na(Blood_Hg_ppm)) %>% 
  arrange(Blood_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

# For google map, you have to give the center of the window you are looking at.
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid
# get the map info
satellitemap <- get_map(location = c(lon = -83, lat = 14), maptype = "satellite", source = "google", zoom = 4)

# setting breaks for the color ramp
my_breaks = c(0.01, 0.1, 1, 3)

bloodmap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat, color = Blood_Hg_ppm, size = Blood_Hg_ppm)) +
  #scale_color_viridis_c() +
  scale_color_gradient(high = "red", low = "yellow", trans = "log", breaks = my_breaks,
                       labels = my_breaks) +
  scale_size(range = c(10,2), trans = "log") + # setting the point size
  #labs(x = "Longitude", y = "Latitude", color = "Whole blood\nTHg (µg/g ww)") + 
  labs(color = "Whole blood\nTHg (µg/g ww)") +
  theme_classic() +
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

ggsave("FinalGraphs/Sayers_Satellite_Circle_Blood_Map.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Satellite_Circle_Blood_Map.tiff", dpi = 1000, width = 10, height = 8)

# Body feather map with points
MapData <- CollectiveData %>%
  filter(!is.na(Body_Hg_ppm)) %>% 
  arrange(Body_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

satellitemap <- get_map(location = c(lon = -83, lat = 5), maptype = "satellite", source = "google", zoom = 4)

# setting breaks for the color ramp
my_breaks = c(0.01, 0.1, 1, 10, 60)

bodymap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Banding_Station_Long, y = Banding_Station_Lat,
                                           color = Body_Hg_ppm, size = Body_Hg_ppm)) +
  scale_color_gradient(high = "red", low = "yellow", trans = "log", breaks = my_breaks,
                       labels = my_breaks) +
  scale_size(range = c(10,2), trans = "log") + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "Body feather\nTHg (µg/g fw)") + 
  #ggtitle("Resident species only") +
  theme_classic() +
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

ggsave("FinalGraphs/Sayers_Satellite_Circle_Body_Map.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Satellite_Circle_Body_Map.tiff", dpi = 1000, width = 10, height = 8)

# Rectrice map with points
MapData <- CollectiveData %>%
  filter(!is.na(Tail_Hg_ppm)) %>% 
  arrange(Tail_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

satellitemap <- get_map(location = c(lon = -83, lat = 5), maptype = "satellite", source = "google", zoom = 4)

# setting breaks for the color ramp
my_breaks = c(0.01, 0.1, 1, 10)
#breaks = c(0, 2, 6, 10)

tailmap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Banding_Station_Long, y = Banding_Station_Lat,
                                           color = Tail_Hg_ppm, size = Tail_Hg_ppm)) +
  scale_color_gradient(high = "red", low = "yellow", trans = "log", breaks = my_breaks,
                       labels = my_breaks) +
  scale_size(range = c(10,2), trans = "log") + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "Tail feather\nTHg (µg/g fw)") + 
  #ggtitle("Resident species only") +
  theme_classic() +
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

ggsave("FinalGraphs/Sayers_Satellite_Circle_Tail_Map.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Satellite_Circle_Tail_Map.tiff", dpi = 1000, width = 10, height = 8)

ggarrange(bloodmap, bodymap, tailmap, labels = c("a)", "b)", "c)"), nrow = 3, ncol = 1)
ggsave("FinalGraphs/Fig5_Facet_Map_AllTissues.jpg", dpi = 1000, width = 10, height = 15)
ggsave("FinalGraphs/Fig5_Facet_Map_AllTissues.tiff", dpi = 800, width = 10, height = 15)



# Sample size maps
MapData <- CollectiveData %>%
  filter(!is.na(Blood_Hg_ppm))

# For google map, you have to give the center of the window you are looking at.
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid
# get the map info
satellitemap <- get_map(location = c(lon = -83, lat = 14), maptype = "satellite", source = "google", zoom = 4)

bloodsamplemap <- ggmap(satellitemap) +
  coord_cartesian() +
  stat_binhex(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #geom_bin2d(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat), bins = 50) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(high = "red", low = "yellow", breaks = c(0, 50, 100, 150, 200, 250)) +
  labs(fill = "Whole blood\nsample size") +
  theme_classic() +
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

ggsave("FinalGraphs/Sayers_Satellite_Blood_Sample_Map.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Satellite_Blood_Sample_Map.tiff", dpi = 1000, width = 10, height = 8)


# Body feather sample map
MapData <- CollectiveData %>%
  filter(!is.na(Body_Hg_ppm))

satellitemap <- get_map(location = c(lon = -83, lat = 5), maptype = "satellite", source = "google", zoom = 4)

bodysamplemap <- ggmap(satellitemap) +
  coord_cartesian() +
  stat_binhex(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat), binwidth = c(1.75, 1.75)) +
  #geom_bin2d(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat), bins = 50) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(high = "red", low = "yellow", breaks = c(0, 50, 100, 150)) +
  labs(x = "Longitude", y = "Latitude", fill = "Body feather\nsample size") + 
  #ggtitle("Resident species only") +
  theme_classic() +
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

ggsave("FinalGraphs/Sayers_Satellite_Body_Sample_Map.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Satellite_Body_Sample_Map.tiff", dpi = 1000, width = 10, height = 8)

# Tail feather sample map
MapData <- CollectiveData %>%
  filter(!is.na(Tail_Hg_ppm))

satellitemap <- get_map(location = c(lon = -83, lat = 5), maptype = "satellite", source = "google", zoom = 4)

tailsamplemap <- ggmap(satellitemap) +
  coord_cartesian() +
  stat_binhex(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat)) +
  #geom_bin2d(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat), bins = 30) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(high = "red", low = "yellow", breaks = c(0, 50, 100, 150, 200)) +
  labs(x = "Longitude", y = "Latitude", fill = "Tail feather\nsample size") +
  #ggtitle("Resident species only") +
  theme_classic() +
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

ggsave("FinalGraphs/Sayers_Satellite_Tail_Sample_Map.jpg", dpi = 1000, width = 10, height = 8)
ggsave("FinalGraphs/Sayers_Satellite_Tail_Sample_Map.tiff", dpi = 1000, width = 10, height = 8)

ggarrange(bloodsamplemap, bodysamplemap, tailsamplemap, labels = c("a)", "b)", "c)"), nrow = 3, ncol = 1)
ggsave("FinalGraphs/Fig6_Facet_Sample_Map_AllTissues.jpg", dpi = 1000, width = 10, height = 15)
ggsave("FinalGraphs/Fig6_Facet_Sample_Map_AllTissues.tiff", dpi = 800, width = 10, height = 15)

