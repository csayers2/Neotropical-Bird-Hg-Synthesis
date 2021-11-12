# Script designed to produce graphs featured in publication and appendix
# Requires objects from Sayers_Main_Script.R

library(tidyverse)
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)

# defining statistical functions for graphs
gmean <- function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

se <- function(x, na.rm = TRUE){
  sd(x)/sqrt(length(x))
}

# resolve namespace conflicts and creating necessary functions
select <- dplyr::select
"%nin%" <- Negate("%in%")

# creating a filtered data set that excludes taxa that were poorly sampled
# and/or aren't associated with terrestrial or freshwater systems
# (ie. large marine piscivores)
# This only excludes 9 birds
GraphingData <- CollectiveData
#filter(Order %nin% c("Anseriformes", "Accipitriformes", "Falconiformes", "Strigiformes", "Suliformes"),
#       Family %nin% c("Laridae (Gulls, Terns, and Skimmers)"))

# my custom color-blind friendly palette that has a bit more pop than some of the stock options
# change order as necessary, but these color selections should have enough contrast
# red, blue, green, orange, pink, gray, yellow, purple, brown, black
mypal <- c("red3", "#0072B2", "forestgreen", "orange2", "orchid3","gray60",
           "gold1", "mediumpurple3", "chocolate4", "black")

# Brewer color-blind friendly palette
display.brewer.pal(n = 9, name = "Set1")

# Hexadecimal color specification 
# red, blue, green, purple, orange, yellow, brown, pink, gray, black
pal <- brewer.pal(n = 9, name = "Set1")
pal[10] <- "black"
pal

# ORDER VARIATION ---------------------------------------------------------

df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  filter(!is.na(Order), !is.na(Concentration)) %>%
  group_by(Order, Tissue_Type) %>%
  summarize(n = n(), gmean = gmean(Concentration), se = se(Concentration)) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(Order ~ Tissue_Type, value.var = c("n", "gmean", "se"), fill = NA) %>%
  select(Order, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         gmean_Blood = gmean_Blood_Hg_ppm, gmean_Contour = gmean_Contour_Hg_ppm, gmean_Tail = gmean_Tail_Hg_ppm,
         se_Blood = se_Blood_Hg_ppm, se_Contour = se_Contour_Hg_ppm, se_Tail = se_Tail_Hg_ppm) %>% 
  # this will display the samples sizes for blood, contour, then tail
  mutate(Order = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                         str_c(Order, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                         if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                 str_c(Order, " (n = ", n_Blood, ", ", n_Contour,")"),
                                 if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                         str_c(Order, " (n = ", n_Blood, ", ", n_Tail,")"),
                                         if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                 str_c(Order, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                 if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                         str_c(Order, " (n = ", n_Blood,")"),
                                                         if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                 str_c(Order, " (n = ", n_Contour,")"),
                                                                 if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                         str_c(Order, " (n = ", n_Tail,")"),
                                                                         Order))))))))

p1 <- df %>%
  select(Order, gmean = gmean_Blood, se = se_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Order, gmean = gmean_Contour, se = se_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(Order, gmean = gmean_Tail, se = se_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(gmean)) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(gmean))/length(gmean)) %>%
  group_by(Order) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                               levels = c("Whole blood", "Contour feather", "Tail feather")))

ggplot(final, mapping = aes(x = gmean, y = reorder(Order, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = gmean + se),
                position = position_dodge(0.9), width = 0.2) +
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

ggsave("Graphs/Sayers_Order_AllTissues_Hist.jpg", dpi = 800, width = 10, height = 8)


# FAMILY VARIATION --------------------------------------------------------

df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Family), !is.na(Concentration)) %>%
  group_by(Family, Tissue_Type) %>%
  summarize(n = n(), gmean = gmean(Concentration), se = se(Concentration)) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(Family ~ Tissue_Type, value.var = c("n", "gmean", "se"), fill = NA) %>%
  select(Family, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         gmean_Blood = gmean_Blood_Hg_ppm, gmean_Contour = gmean_Contour_Hg_ppm, gmean_Tail = gmean_Tail_Hg_ppm,
         se_Blood = se_Blood_Hg_ppm, se_Contour = se_Contour_Hg_ppm, se_Tail = se_Tail_Hg_ppm,) %>% 
  # this will display the samples sizes for blood, contour, then tail
  mutate(Family = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                          str_c(Family, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                          if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                  str_c(Family, " (n = ", n_Blood, ", ", n_Contour,")"),
                                  if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                          str_c(Family, " (n = ", n_Blood, ", ", n_Tail,")"),
                                          if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                  str_c(Family, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                  if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                          str_c(Family, " (n = ", n_Blood,")"),
                                                          if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                  str_c(Family, " (n = ", n_Contour,")"),
                                                                  if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                          str_c(Family, " (n = ", n_Tail,")"),
                                                                          Family))))))))

p1 <- df %>%
  select(Family, gmean = gmean_Blood, se = se_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Family, gmean = gmean_Contour, se = se_Contour, n_Contour = n_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(Family, gmean = gmean_Tail, se = se_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(gmean)) %>%
  # filtering by sample size so that the names fit on the axis
  group_by(Family) %>% 
  mutate(sum = sum(n_Blood, n_Contour, n_Tail, na.rm = T)) %>%
  filter(sum > 9) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(gmean))/length(gmean)) %>%
  group_by(Family) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Contour feather", "Tail feather")))

ggplot(final, mapping = aes(x = gmean, y = reorder(Family, max_percentile), fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = gmean + se),
                position = position_dodge(0.9), width = 0.2) +
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

ggsave("Graphs/Sayers_Family_AllTissues_Hist.jpg", dpi = 800, width = 10, height = 8)


# SPECIES VARIATION --------------------------------------------------------

df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Species_Common_Name), !is.na(Concentration)) %>%
  group_by(Species_Common_Name, Tissue_Type) %>%
  summarize(n = n(), gmean = gmean(Concentration), se = se(Concentration)) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(Species_Common_Name ~ Tissue_Type, value.var = c("n", "gmean", "se"), fill = NA) %>%
  select(Species_Common_Name, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         gmean_Blood = gmean_Blood_Hg_ppm, gmean_Contour = gmean_Contour_Hg_ppm, gmean_Tail = gmean_Tail_Hg_ppm,
         se_Blood = se_Blood_Hg_ppm, se_Contour = se_Contour_Hg_ppm, se_Tail = se_Tail_Hg_ppm,) %>% 
  mutate(Species_Common_Name = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                       str_c(Species_Common_Name, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                                       if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                               str_c(Species_Common_Name, " (n = ", n_Blood, ", ", n_Contour,")"),
                                               if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                       str_c(Species_Common_Name, " (n = ", n_Blood, ", ", n_Tail,")"),
                                                       if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                               str_c(Species_Common_Name, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                               if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                                       str_c(Species_Common_Name, " (n = ", n_Blood,")"),
                                                                       if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                               str_c(Species_Common_Name, " (n = ", n_Contour,")"),
                                                                               if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                                       str_c(Species_Common_Name, " (n = ", n_Tail,")"),
                                                                                       Species_Common_Name))))))))

p1 <- df %>%
  select(Species_Common_Name, gmean = gmean_Blood, se = se_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Species_Common_Name, gmean = gmean_Contour, se = se_Contour, n_Contour = n_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(Species_Common_Name, gmean = gmean_Tail, se = se_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(gmean)) %>%
  # filtering by sample size so that the names fit on the axis
  group_by(Species_Common_Name) %>% 
  mutate(sum = sum(n_Blood, n_Contour, n_Tail, na.rm = T)) %>%
  filter(sum > 9) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(gmean))/length(gmean)) %>%
  group_by(Species_Common_Name) %>%
  mutate(max_percentile = max(percentile)) %>%
  #filter(max_percentile >= 0.50) %>% 
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Contour feather", "Tail feather")))
  
ggplot(final, mapping = aes(x = gmean, y = reorder(Species_Common_Name, max_percentile),
                            fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = gmean + se),
                position = position_dodge(0.9), width = 0.2) +
  geom_col(position = "dodge") +
  labs(x = "THg (µg/g)", y = bquote("Species (n ≥ 10)")) +
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

ggsave("Graphs/Sayers_Species_AllTissues_Hist.jpg", dpi = 800, width = 8, height = 10)


# TROPHIC NICHE x PRIMARY HABITAT ---------------------------------------------

df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Trophic_Niche), !is.na(HAB1), !is.na(Concentration)) %>%
  pivot_longer(c(Trophic_Niche, HAB1),
               names_to = "Category", values_to = "Type") %>%
  group_by(Category, Type, Tissue_Type) %>%
  summarize(n = n(), gmean = gmean(Concentration), se = se(Concentration)) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(Category + Type ~ Tissue_Type, value.var = c("n", "gmean", "se"), fill = NA) %>%
  select(Category, Type, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         gmean_Blood = gmean_Blood_Hg_ppm, gmean_Contour = gmean_Contour_Hg_ppm, gmean_Tail = gmean_Tail_Hg_ppm,
         se_Blood = se_Blood_Hg_ppm, se_Contour = se_Contour_Hg_ppm, se_Tail = se_Tail_Hg_ppm,) %>% 
  mutate(Type = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                        str_c(Type, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                        if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                str_c(Type, " (n = ", n_Blood, ", ", n_Contour,")"),
                                if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                        str_c(Type, " (n = ", n_Blood, ", ", n_Tail,")"),
                                        if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                str_c(Type, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                        str_c(Type, " (n = ", n_Blood,")"),
                                                        if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                str_c(Type, " (n = ", n_Contour,")"),
                                                                if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                        str_c(Type, " (n = ", n_Tail,")"),
                                                                        Type)))))))) %>% 
  # renaming category values for later plotting
  mutate(Category = if_else(Category == "Trophic_Niche", "Trophic niche", "Primary habitat"))

p1 <- df %>%
  select(Category, Type, gmean = gmean_Blood, se = se_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Category, Type, gmean = gmean_Contour, se = se_Contour, n_Contour = n_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(Category, Type, gmean = gmean_Tail, se = se_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(gmean)) %>%
  # calculating percentile rank for each tissue
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(gmean))/length(gmean)) %>%
  group_by(Category, Type) %>% 
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Contour feather", "Tail feather"))) %>% 
  transform(Category = factor(Category,
                                 levels = c("Trophic niche", "Primary habitat")))

ggplot(final, mapping = aes(x = gmean, y = reorder(Type, max_percentile),
                            fill = Tissue_Type)) +
  geom_errorbar(aes(xmin = 0, xmax = gmean + se),
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

ggsave("Graphs/Sayers_TrophicxHabitat_AllTissues_Hist.jpg", dpi = 800, width = 10, height = 8)


# RISK ASSESSMENT ---------------------------------------------------------

# Risk graph for blood
bloodrisk <- GraphingData %>%
  filter(!is.na(Blood_Hg_ppm), !is.na(Species_Common_Name),
         Species_Code != "BIRD") %>%
  group_by(Species_Common_Name) %>%
  mutate(Species_Common_Name = str_c(Species_Common_Name, " (n = ", n(), ")"),
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
    reorder(reorder(reorder(Species_Common_Name, desc(Species_Common_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightskyblue")) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)",
       fill = "Whole blood\nrisk categories") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

ggsave("Graphs/Sayers_Risk_Blood.jpg", dpi = 800, width = 10, height = 8)

# Risk graph for contour feathers
contourrisk <- GraphingData %>%
  filter(!is.na(Species_Common_Name), !is.na(Contour_Hg_ppm), Species_Code != "BIRD") %>%
  group_by(Species_Common_Name) %>%
  mutate(Species_Common_Name = str_c(Species_Common_Name, " (n = ", n(), ")"),
         Risk_Levels = cut(Contour_Hg_ppm, breaks = c(0, 2.4, 3.4, 4.5, 5.3, Inf),
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
    reorder(reorder(reorder(Species_Common_Name, desc(Species_Common_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightskyblue")) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)", fill = "Contour feather\nrisk categories") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        aspect.ratio = 1)

ggsave("Graphs/Sayers_Risk_Contour.jpg", dpi = 800, width = 10, height = 8)

# Risk graph for rectrices
tailrisk <- GraphingData %>%
  filter(!is.na(Species_Common_Name), !is.na(Tail_Hg_ppm), Species_Code != "BIRD") %>%
  group_by(Species_Common_Name) %>%
  mutate(Species_Common_Name = str_c(Species_Common_Name, " (n = ", n(), ")"),
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
    reorder(reorder(reorder(Species_Common_Name, desc(Species_Common_Name)), Score), Prop_None),
    Prop_Low), Prop_Med), Prop_High), Prop_Ext), fill = Risk_Levels)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("black", "red3", "darkorange", "gold", "lightskyblue")) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) + 
  labs(x = "Proportion of individuals sampled", y = "Species (n ≥ 5)", fill = "Tail feather\nrisk categories") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        #axis.title.y = element_blank(),
        legend.title = element_text(face = "bold"),
        aspect.ratio = 1)

ggsave("Graphs/Sayers_Risk_Tail.jpg", dpi = 800, width = 10, height = 8)


ggarrange(bloodrisk, contourrisk, tailrisk, labels = c("A)", "B)", "C)"), nrow = 2, ncol = 2)
ggarrange(contourrisk, tailrisk, labels = c("A)", "B)"), nrow = 2, ncol = 1)
ggsave("Graphs/Sayers_Facet_Risk_Feathers.jpg", dpi = 800, width = 10, height = 8)




# MAPS ------------------------------------------------------------------
library(gcookbook)
library(sp)
library(raster)

# Google satellite imagery as a background
library(ggmap)
register_google(key = "Insert your Google API key here", write = TRUE)

# Blood map with points as circles
MapData <- GraphingData %>%
  filter(!is.na(Blood_Hg_ppm)) %>% 
  arrange(Blood_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

# For google map, you have to give the center of the window you are looking at.
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid
# get the map info
satellitemap <- get_map(location = c(lon = -83, lat = 14), maptype = "satellite", source = "google", zoom = 4)

bloodmap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Site_Long, y = Site_Lat, color = Blood_Hg_ppm, size = Blood_Hg_ppm)) +
  #scale_color_viridis_c() +
  scale_color_gradient(high = "red", low = "yellow") +
  scale_size(range = c(8,2)) + # setting the point size
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
  guides(size = F) # removing size legend 

ggsave("Graphs/Sayers_Satellite_Circle_Blood_Map.jpg", dpi = 800, width = 10, height = 8)

# Contour feather map with points
MapData <- GraphingData %>%
  filter(!is.na(Contour_Hg_ppm)) %>% 
  arrange(Contour_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

satellitemap <- get_map(location = c(lon = -83, lat = 5), maptype = "satellite", source = "google", zoom = 4)

contourmap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Banding_Station_Long, y = Banding_Station_Lat,
                                           color = Contour_Hg_ppm, size = Contour_Hg_ppm)) +
  scale_color_gradient(high = "red", low = "yellow") +
  scale_size(range = c(8,2)) + # setting the point size
  labs(x = "Longitude", y = "Latitude", color = "Contour feather\nTHg (µg/g fw)") + 
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
  guides(size = F) # removing size legend 

ggsave("Graphs/Sayers_Satellite_Circle_Contour_Map.jpg", dpi = 800, width = 10, height = 8)

# Rectrice map with points
MapData <- GraphingData %>%
  filter(!is.na(Tail_Hg_ppm)) %>% 
  arrange(Tail_Hg_ppm) %>% # this is so that the points will be plotted big -> small
  # making sure all individuals can be plotted
  mutate(Banding_Station_Lat = if_else(is.na(Banding_Station_Lat), Site_Lat, Banding_Station_Lat),
         Banding_Station_Long = if_else(is.na(Banding_Station_Long), Site_Long, Banding_Station_Long))

satellitemap <- get_map(location = c(lon = -83, lat = 5), maptype = "satellite", source = "google", zoom = 4)

tailmap <- ggmap(satellitemap) +
  geom_point(data = MapData, mapping = aes(x = Banding_Station_Long, y = Banding_Station_Lat,
                                           color = Tail_Hg_ppm, size = Tail_Hg_ppm)) +
  scale_color_gradient(high = "red", low = "yellow", breaks = c(0, 2, 6, 10)) +
  scale_size(range = c(8,2)) + # setting the point size
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
  guides(size = F) # removing size legend 

ggsave("Graphs/Sayers_Satellite_Circle_Tail_Map.jpg", dpi = 800, width = 10, height = 8)

ggarrange(bloodmap, contourmap, tailmap, labels = c("A)", "B)", "C)"), nrow = 2, ncol = 2)
ggsave("Graphs/Sayers_Facet_Map_AllTissues.jpg", dpi = 800, width = 15, height = 10)



# COEFFICIENT OF VARIATION ---------------------------------------------------------

# Order
df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Order), !is.na(Concentration)) %>%
  group_by(Order, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration), cv = sd/mean) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(Order ~ Tissue_Type, value.var = c("n", "cv"), fill = NA) %>%
  select(Order, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         cv_Blood = cv_Blood_Hg_ppm, cv_Contour = cv_Contour_Hg_ppm, cv_Tail = cv_Tail_Hg_ppm) %>% 
  # this will display the samples sizes for blood, contour, then tail
  mutate(Order = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                         str_c(Order, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                         if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                 str_c(Order, " (n = ", n_Blood, ", ", n_Contour,")"),
                                 if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                         str_c(Order, " (n = ", n_Blood, ", ", n_Tail,")"),
                                         if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                 str_c(Order, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                 if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                         str_c(Order, " (n = ", n_Blood,")"),
                                                         if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                 str_c(Order, " (n = ", n_Contour,")"),
                                                                 if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                         str_c(Order, " (n = ", n_Tail,")"),
                                                                         Order))))))))
p1 <- df %>%
  select(Order, cv = cv_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Order, cv = cv_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(Order, cv = cv_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(cv)) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(cv))/length(cv)) %>%
  group_by(Order) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Contour feather", "Tail feather")))

ggplot(final, mapping = aes(x = cv, y = reorder(Order, max_percentile), fill = Tissue_Type)) +
  geom_col() +
  labs(x = "Coefficient of THg variation", y = "Order") +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) +
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

ggsave("Graphs/Sayers_Order_CV_AllTissues_Hist.jpg", dpi = 800, width = 10, height = 8)


# Family
df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Family), !is.na(Concentration)) %>%
  group_by(Family, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration), cv = sd/mean) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(Family ~ Tissue_Type, value.var = c("n", "cv"), fill = NA) %>%
  select(Family, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         cv_Blood = cv_Blood_Hg_ppm, cv_Contour = cv_Contour_Hg_ppm, cv_Tail = cv_Tail_Hg_ppm) %>% 
  # this will display the samples sizes for blood, contour, then tail
  mutate(Family = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                         str_c(Family, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                         if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                 str_c(Family, " (n = ", n_Blood, ", ", n_Contour,")"),
                                 if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                         str_c(Family, " (n = ", n_Blood, ", ", n_Tail,")"),
                                         if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                 str_c(Family, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                 if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                         str_c(Family, " (n = ", n_Blood,")"),
                                                         if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                 str_c(Family, " (n = ", n_Contour,")"),
                                                                 if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                         str_c(Family, " (n = ", n_Tail,")"),
                                                                         Family))))))))
p1 <- df %>%
  select(Family, cv = cv_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Family, cv = cv_Contour, n_Contour = n_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(Family, cv = cv_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(cv)) %>%
  # filtering by sample size so that the names fit on the axis
  group_by(Family) %>% 
  mutate(sum = sum(n_Blood, n_Contour, n_Tail, na.rm = T)) %>%
  filter(sum > 9) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(cv))/length(cv)) %>%
  group_by(Family) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Contour feather", "Tail feather")))

ggplot(final, mapping = aes(x = cv, y = reorder(Family, max_percentile), fill = Tissue_Type)) +
  geom_col() +
  labs(x = "Coefficient of THg variation", y = "Family (n ≥ 10)") +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) +
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

ggsave("Graphs/Sayers_Family_CV_AllTissues_Hist.jpg", dpi = 800, width = 10, height = 8)

# Species
df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Species_Common_Name), !is.na(Concentration)) %>%
  group_by(Species_Common_Name, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration), cv = sd/mean) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(Species_Common_Name ~ Tissue_Type, value.var = c("n", "cv"), fill = NA) %>%
  select(Species_Common_Name, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         cv_Blood = cv_Blood_Hg_ppm, cv_Contour = cv_Contour_Hg_ppm, cv_Tail = cv_Tail_Hg_ppm) %>% 
  # this will display the samples sizes for blood, contour, then tail
  mutate(Species_Common_Name = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                          str_c(Species_Common_Name, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                          if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                  str_c(Species_Common_Name, " (n = ", n_Blood, ", ", n_Contour,")"),
                                  if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                          str_c(Species_Common_Name, " (n = ", n_Blood, ", ", n_Tail,")"),
                                          if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                  str_c(Species_Common_Name, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                  if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                          str_c(Species_Common_Name, " (n = ", n_Blood,")"),
                                                          if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                  str_c(Species_Common_Name, " (n = ", n_Contour,")"),
                                                                  if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                          str_c(Species_Common_Name, " (n = ", n_Tail,")"),
                                                                          Species_Common_Name))))))))
p1 <- df %>%
  select(Species_Common_Name, cv = cv_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Species_Common_Name, cv = cv_Contour, n_Contour = n_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(Species_Common_Name, cv = cv_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(cv)) %>%
  # filtering by sample size so that the names fit on the axis
  group_by(Species_Common_Name) %>% 
  mutate(sum = sum(n_Blood, n_Contour, n_Tail, na.rm = T)) %>%
  filter(sum > 9) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(cv))/length(cv)) %>%
  group_by(Species_Common_Name) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Contour feather", "Tail feather")))

ggplot(final, mapping = aes(x = cv, y = reorder(Species_Common_Name, max_percentile), fill = Tissue_Type)) +
  geom_col() +
  labs(x = "Coefficient of THg variation", y = "Species (n ≥ 10)") +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) +
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

ggsave("Graphs/Sayers_Species_CV_AllTissues_Hist.jpg", dpi = 800, width = 10, height = 8)

# Trophic Niche
df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(Trophic_Niche), !is.na(Concentration)) %>%
  group_by(Trophic_Niche, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration), cv = sd/mean) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(Trophic_Niche ~ Tissue_Type, value.var = c("n", "cv"), fill = NA) %>%
  select(Trophic_Niche, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         cv_Blood = cv_Blood_Hg_ppm, cv_Contour = cv_Contour_Hg_ppm, cv_Tail = cv_Tail_Hg_ppm) %>% 
  # this will display the samples sizes for blood, contour, then tail
  mutate(Trophic_Niche = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                       str_c(Trophic_Niche, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                                       if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                               str_c(Trophic_Niche, " (n = ", n_Blood, ", ", n_Contour,")"),
                                               if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                       str_c(Trophic_Niche, " (n = ", n_Blood, ", ", n_Tail,")"),
                                                       if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                               str_c(Trophic_Niche, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                               if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                                       str_c(Trophic_Niche, " (n = ", n_Blood,")"),
                                                                       if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                               str_c(Trophic_Niche, " (n = ", n_Contour,")"),
                                                                               if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                                       str_c(Trophic_Niche, " (n = ", n_Tail,")"),
                                                                                       Trophic_Niche))))))))
p1 <- df %>%
  select(Trophic_Niche, cv = cv_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(Trophic_Niche, cv = cv_Contour, n_Contour = n_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(Trophic_Niche, cv = cv_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(cv)) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(cv))/length(cv)) %>% 
  group_by(Trophic_Niche) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Contour feather", "Tail feather")))

ggplot(final, mapping = aes(x = cv, y = reorder(Trophic_Niche, max_percentile), fill = Tissue_Type)) +
  geom_col() +
  labs(x = "Coefficient of THg variation", y = "Trophic niche") +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) +
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

ggsave("Graphs/Sayers_Trophic_CV_AllTissues_Hist.jpg", dpi = 800, width = 10, height = 8)


# Primary habitat
df <- GraphingData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Contour_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  filter(!is.na(HAB1), !is.na(Concentration)) %>%
  group_by(HAB1, Tissue_Type) %>%
  summarize(n = n(), mean = mean(Concentration), sd = sd(Concentration), cv = sd/mean) %>%
  # turning the df into a data.table object so dcast can accept more than one value.var
  data.table() %>% 
  data.table::dcast(HAB1 ~ Tissue_Type, value.var = c("n", "cv"), fill = NA) %>%
  select(HAB1, n_Blood = n_Blood_Hg_ppm, n_Contour = n_Contour_Hg_ppm, n_Tail = n_Tail_Hg_ppm,
         cv_Blood = cv_Blood_Hg_ppm, cv_Contour = cv_Contour_Hg_ppm, cv_Tail = cv_Tail_Hg_ppm) %>% 
  # this will display the samples sizes for blood, contour, then tail
  mutate(HAB1 = if_else(!is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                 str_c(HAB1, " (n = ", n_Blood, ", ", n_Contour, ", ", n_Tail, ")"),
                                 if_else(!is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                         str_c(HAB1, " (n = ", n_Blood, ", ", n_Contour,")"),
                                         if_else(!is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                 str_c(HAB1, " (n = ", n_Blood, ", ", n_Tail,")"),
                                                 if_else(is.na(n_Blood) & !is.na(n_Contour) & !is.na(n_Tail),
                                                         str_c(HAB1, " (n = ", n_Contour, ", ", n_Tail,")"),
                                                         if_else(!is.na(n_Blood) & is.na(n_Contour) & is.na(n_Tail),
                                                                 str_c(HAB1, " (n = ", n_Blood,")"),
                                                                 if_else(is.na(n_Blood) & !is.na(n_Contour) & is.na(n_Tail),
                                                                         str_c(HAB1, " (n = ", n_Contour,")"),
                                                                         if_else(is.na(n_Blood) & is.na(n_Contour) & !is.na(n_Tail),
                                                                                 str_c(HAB1, " (n = ", n_Tail,")"),
                                                                                 HAB1))))))))
p1 <- df %>%
  select(HAB1, cv = cv_Blood, n_Blood = n_Blood) %>%
  mutate(Tissue_Type = "Whole blood")
p2 <- df %>%
  select(HAB1, cv = cv_Contour, n_Contour = n_Contour) %>%
  mutate(Tissue_Type = "Contour feather")
p3 <- df %>%
  select(HAB1, cv = cv_Tail, n_Tail = n_Tail) %>%
  mutate(Tissue_Type = "Tail feather")

final <- rbind(p1, p2, p3, fill = T) %>% 
  filter(!is.na(cv)) %>%
  # creating a system to better rank the y axis
  group_by(Tissue_Type) %>% # calculating percentile rank for each tissue
  mutate(percentile = trunc(rank(cv))/length(cv)) %>% 
  group_by(HAB1) %>%
  mutate(max_percentile = max(percentile)) %>%
  # this function is critical to order the facets properly
  transform(Tissue_Type = factor(Tissue_Type,
                                 levels = c("Whole blood", "Contour feather", "Tail feather")))

ggplot(final, mapping = aes(x = cv, y = reorder(HAB1, max_percentile), fill = Tissue_Type)) +
  geom_col() +
  labs(x = "Coefficient of THg variation", y = "Primary habitat") +
  scale_x_continuous(labels = scales::percent, expand = c(0,0)) +
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

ggsave("Graphs/Sayers_Trophic_CV_AllTissues_Hist.jpg", dpi = 800, width = 10, height = 8)
